# vkmz v1.4dev0

vkmz predicts molecular formulas by searching a known mass-formula dictionary
for a feature observed by a mass spectrometer. Elemental ratios forpredicted-features
are calculated to create the carbon-to-oxygen and carbon-to-hydrogen axis of a
van Krevelen Diagram (VKD). VKD's are a convenient visualization tool for
briefly conveying the constituence of a complex MS mixture (e.g., untargetted
plant metabolomics). As output predicted-feature are saved to a tabular file,
an interactive VKD web page, and other optional formats.

## Installation

vkmz requires Python version 3.6 or greater.

### setup.py

Clone or otherwise repo and go to the root directory. Install vkmz as any other
python package:
```
python3 setup.py install
```

### Conda

Version 1.4dev0 has not been wrapped for conda.

### Galaxy

Version 1.4dev0 has not been wrapped for Galaxy.

## Input Data

Can either parse a tabular file or Workflow4Metabolomics' XCMS tabular as input.

Input MS data can be given in two "modes", (1) tabular or (2) Workflow4Metabolomics'
XCMS for Galaxy (W4M-XCMS) files.

Tabular mode requires a single tabular file as input and  must include the columns
"sample_name", "polarity", "mz", "rt", and "intensity". Each row represents a 
feature. Optionally a "charge" column can exist.

W4M-XCMS mode requires the sample metadata, variable metadata, and data matrix
files generated with W4M-XCMS. Feature charge infomration can be read from the
variable metadata file if it has been annotated with CAMERA.

Polarity values should be either "positive" or "negative".

If feature charge information is present, features without charge information
will be removed. If CAMERA annotation is present, only monoisotopic features
will be kept. An argument flag (`--impute-charge`) can be set to disable removing
features without charge annotation. Users should be wary of false results when
using this non-default option.

## Output

vkmz always outputs tabular and html files. Optionally, vkmz can output JSON
and SQL as well.

## Command Line Interface

### Quick start

```
vkmz tabular --input test-data/tabular.tabular --output foo --error 10
vkmz xcms -xd test-data/datamatrix.tabular -xv test-data/variableMetadata.tabular -xs test-data/sampleMetadata.tabular -o foo -e 10 --impute
```

### All Arguments

vkmz runs in two modes: "tabular" and "xcms". xcms refers to Workflow4Metabolomics'
XCMS for Galaxy.

```
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
```
