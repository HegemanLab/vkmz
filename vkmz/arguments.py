#!/usr/bin/env python
"""vkmz.arguments

foo bar
"""


import argparse
import os


parser = argparse.ArgumentParser()
sub_parser = parser.add_subparsers(help="Select mode:", dest="mode")
sub_parser.required = True

# Tabular mode arguments
parse_tabular = sub_parser.add_parser("tabular", help="Use tabular data as input.")
parse_tabular.add_argument("--input", "-i", required=True, help="Path to tabular file.")

# XCMS-tabular mode arguments
parse_xcms = sub_parser.add_parser("xcms", help="Use XCMS data as input.")
parse_xcms.add_argument(
    "--data-matrix",
    "-xd",
    required=True,
    nargs="?",
    type=str,
    help="Path to XCMS data matrix file.",
)
parse_xcms.add_argument(
    "--sample-metadata",
    "-xs",
    required=True,
    nargs="?",
    type=str,
    help="Path to XCMS sample metadata file.",
)
parse_xcms.add_argument(
    "--variable-metadata",
    "-xv",
    required=True,
    nargs="?",
    type=str,
    help="Path to XCMS variable metadata file.",
)

# Arguments for all modes
for mode in [parse_tabular, parse_xcms]:
    mode.add_argument(
        "--output",
        "-o",
        required=True,
        nargs="?",
        type=str,
        help="Specify output file path.",
    )
    mode.add_argument(
        "--error",
        "-e",
        required=True,
        nargs="?",
        type=float,
        help="Mass error of MS data in parts-per-million.",
    )
    mode.add_argument(
        "--json", "-j", action="store_true", help="Set JSON flag to save JSON output."
    )
    mode.add_argument(
        "--sql", "-s", action="store_true", help="Set SQL flag to save SQL output."
    )
    mode.add_argument(
        "--metadata",
        "-m",
        action="store_true",
        help="Set metadata flag to save the tools metadata.",
    )
    mode.add_argument(
        "--database",
        "-db",
        nargs="?",
        default="databases/bmrb-light.tsv",
        help="Define database of known formula-mass pairs.",
    )
    mode.add_argument(
        "--prefix", nargs="?", type=str, help="Define tool directory path."
    )
    mode.add_argument(
        "--polarity",
        "-p",
        choices=["positive", "negative"],
        help="Set all polarities to positive or negative.",
    )
    mode.add_argument(
        "--neutral",
        "-n",
        action="store_true",
        help="Set masses in input data are neutral.",
    )
    mode.add_argument(
        "--alternate",
        "-a",
        action="store_true",
        help="Set to keep features with multiple predictions.",
    )
    mode.add_argument(
        "--charge", "-c", action="store_true", help="Set if input data contains charge."
    )

# create constants
args = parser.parse_args()
ALTERNATE = getattr(args, "alternate")
DATABASE = getattr(args, "database")
CHARGE = getattr(args, "charge")
JSON = getattr(args, "json")
MASS_ERROR = getattr(args, "error")
METADATA = getattr(args, "metadata")
MODE = getattr(args, "mode")
NEUTRAL = getattr(args, "neutral")
OUTPUT = getattr(args, "output")
POLARITY = getattr(args, "polarity")
PREFIX = getattr(args, "prefix")
# TODO: check if not against PEP8
if not PREFIX:
    PREFIX = os.path.abspath(os.path.dirname(__file__))
SQL = getattr(args, "sql")
# MASS and FORMULA are used as indexable dictionaries
MASS = []
FORMULA = []
try:
    with open(os.path.join(PREFIX, DATABASE), "r") as tabular:
        next(tabular)  # skip header
        for row in tabular:
            mass, formula = row.split()
            MASS.append(float(mass))
            FORMULA.append(formula)
except:
    print(f"An error occured while reading the {DATABASE} database.")
    raise
MAX_MASS_INDEX = len(MASS) - 1
#!/usr/bin/env python

import os
from vkmz.arguments import args, JSON, METADATA, MODE, parser, SQL
from vkmz.read import tabular as readTabular, xcmsTabular as readXcmsTabular
from vkmz.predict import predict
import vkmz.write as write
