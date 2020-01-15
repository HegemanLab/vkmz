#!/usr/bin/env python

import argparse
import os

parser = argparse.ArgumentParser()
sub_parser = parser.add_subparsers(help="Select mode:", dest="mode")
sub_parser.required = True

# Tabular mode arguments
parse_tabular = sub_parser.add_parser("tabular", help="Tabular data mode")
parse_tabular.add_argument("--input", "-i", required=True, help="Path to tabular file.")

# XCMS-tabular mode arguments
parse_xcms = sub_parser.add_parser("w4m-xcms", help="W4M-XCMS data mode")
parse_xcms.add_argument(
    "--data-matrix",
    "-xd",
    required=True,
    nargs="?",
    type=str,
    help="Path to XCMS data matrix file",
)
parse_xcms.add_argument(
    "--sample-metadata",
    "-xs",
    required=True,
    nargs="?",
    type=str,
    help="Path to XCMS sample metadata file",
)
parse_xcms.add_argument(
    "--variable-metadata",
    "-xv",
    required=True,
    nargs="?",
    type=str,
    help="Path to XCMS variable metadata file",
)

# Required modes for tabular and xcms
for mode in [parse_tabular, parse_xcms]:
    mode.add_argument(
        "--error",
        "-e",
        required=True,
        nargs="?",
        type=float,
        help="Mass error of MS data in parts-per-million",
    )

# Annotated molecular formula mode
parse_formula = sub_parser.add_parser(
    "formula", help="Annotated molecular formula mode"
)
parse_formula.add_argument(
    "--input", "-i", required=True, help="Path to tabular formula file."
)

# all modes
for mode in [parse_formula, parse_tabular, parse_xcms]:
    mode.add_argument(
        "--output",
        "-o",
        required=True,
        nargs="?",
        type=str,
        help="Specify output file path",
    )
    mode.add_argument(
        "--json", "-j", action="store_true", help="Set JSON flag to save JSON output"
    )
    mode.add_argument(
        "--sql", "-s", action="store_true", help="Set SQL flag to save SQL output"
    )
    mode.add_argument(
        "--metadata",
        "-m",
        action="store_true",
        help="Set metadata flag to save argument metadata",
    )
    mode.add_argument(
        "--database",
        "-db",
        nargs="?",
        default="databases/bmrb-light.tsv",
        help="Define path to custom database of known formula-mass pairs",
    )
    mode.add_argument(
        "--prefix",
        nargs="?",
        type=str,
        help='Define path prefix to support files ("d3.html" and database directory)',
    )
    mode.add_argument(
        "--polarity",
        "-p",
        choices=["positive", "negative"],
        help="Set flag to force polarity of all features to positive or negative",
    )
    mode.add_argument(
        "--neutral",
        "-n",
        action="store_true",
        help="Set flag if input data contains neutral feature mass instead of mz",
    )
    mode.add_argument(
        "--alternate",
        "-a",
        action="store_true",
        help="Set flag to keep features with multiple predictions",
    )
    mode.add_argument(
        "--impute-charge",
        "--impute",
        action="store_true",
        help='Set flag to impute "1" for missing charge information',
    )

# create constants
args = parser.parse_args()
JSON = getattr(args, "json")
METADATA = getattr(args, "metadata")
MODE = getattr(args, "mode")
SQL = getattr(args, "sql")
IMPUTE = getattr(args, "impute_charge")
POLARITY = getattr(args, "polarity")
ALTERNATE = getattr(args, "alternate")
DATABASE = getattr(args, "database")
if "error" in args:
    MASS_ERROR = getattr(args, "error")
else:
    MASS_ERROR = "NA"
NEUTRAL = getattr(args, "neutral")
OUTPUT = getattr(args, "output")
PREFIX = getattr(args, "prefix")
if not PREFIX:
    PREFIX = os.path.abspath(os.path.dirname(__file__))
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
    print(f"An error occurred while reading the {DATABASE} database.")
    raise
MAX_MASS_INDEX = len(MASS) - 1
