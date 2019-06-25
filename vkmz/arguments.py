#!/usr/bin/env python
"""read command line arguments and generate constants

vkmz runs in two modes: "tabular" and "xcms". xcms refers to Workflow4Metabolomics'
XCMS for Galaxy.

tabular mode arguments:
    usage: vkmz tabular [-h] --input INPUT --output [OUTPUT] --error [ERROR]
                        [--json] [--sql] [--metadata] [--database [DATABASE]]
                        [--prefix [PREFIX]] [--polarity {positive,negative}]
                        [--neutral] [--alternate] [--impute-charge]

    required arguments:
      --input [INPUT], -i [INPUT]
                            Path to tabular file

    optional arguments:
      --help, -h            tabular mode help message and exit

xcms mode arguments:
    usage: vkmz xcms [-h] --data-matrix [DATA_MATRIX] --sample-metadata
                     [SAMPLE_METADATA] --variable-metadata [VARIABLE_METADATA]
                     --output [OUTPUT] --error [ERROR] [--json] [--sql]
                     [--metadata] [--database [DATABASE]] [--prefix [PREFIX]]
                     [--polarity {positive,negative}] [--neutral] [--alternate]
                     [--impute-charge]

    required arguments:
      --data-matrix [DATA_MATRIX], -xd [DATA_MATRIX]
                            Path to XCMS data matrix file
      --sample-metadata [SAMPLE_METADATA], -xs [SAMPLE_METADATA]
                            Path to XCMS sample metadata file
      --variable-metadata [VARIABLE_METADATA], -xv [VARIABLE_METADATA]
                            Path to XCMS variable metadata file

    optional arguments:
      --help, -h            xcms mode help message and exit

mode shared arguments:
    required arguments:
      --output [OUTPUT], -o [OUTPUT]
                            Specify output file path.
      --error [ERROR], -e [ERROR]
                            Mass error of MS data in parts-per-million.

    optional arguments:
      --json, -j            Set JSON flag to save JSON output
      --sql, -s             Set SQL flag to save SQL output
      --metadata, -m        Set metadata flag to save argument metadata
      --database [DATABASE], -db [DATABASE]
                            Define path to custom database of known formula-mass
                            pairs
      --prefix [PREFIX]     Define path prefix to support files ("d3.html" and
                            database directory
      --polarity {positive,negative}, -p {positive,negative}
                            Set flag to force polarity of all features to positive
                            or negative
      --neutral, -n         Set flag if input data contains neutral feature mass
                            instead of mz
      --alternate, -a       Set flag to keep features with multiple predictions
      --impute-charge, --impute
                            Set flag to impute "1" for missing charge information
"""

import argparse
import os

parser = argparse.ArgumentParser()
sub_parser = parser.add_subparsers(help="Select mode:", dest="mode")
sub_parser.required = True

# Tabular mode arguments
parse_tabular = sub_parser.add_parser("tabular", help="Use tabular data as input")
parse_tabular.add_argument("--input", "-i", required=True, help="Path to tabular file.")

# XCMS-tabular mode arguments
parse_xcms = sub_parser.add_parser("xcms", help="Use XCMS data as input")
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

# Arguments for all modes
for mode in [parse_tabular, parse_xcms]:
    mode.add_argument(
        "--output",
        "-o",
        required=True,
        nargs="?",
        type=str,
        help="Specify output file path",
    )
    mode.add_argument(
        "--error",
        "-e",
        required=True,
        nargs="?",
        type=float,
        help="Mass error of MS data in parts-per-million",
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
ALTERNATE = getattr(args, "alternate")
DATABASE = getattr(args, "database")
IMPUTE = getattr(args, "impute_charge")
JSON = getattr(args, "json")
MASS_ERROR = getattr(args, "error")
METADATA = getattr(args, "metadata")
MODE = getattr(args, "mode")
NEUTRAL = getattr(args, "neutral")
OUTPUT = getattr(args, "output")
POLARITY = getattr(args, "polarity")
PREFIX = getattr(args, "prefix")
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
    print(f"An error occurred while reading the {DATABASE} database.")
    raise
MAX_MASS_INDEX = len(MASS) - 1
