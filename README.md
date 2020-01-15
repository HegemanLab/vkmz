# vkmz v1.4.5

vkmz is a visualization and annotation tool for high resolution liquid chromatography mass spectrometry and metabolomics.

## Overview

vkmz reads LC-MS data and creates a *van Krevelen Diagram* (VKD). A VKD plots moleclues on a graph by their elemental ratios. A tradional VKD plots a molecule's oxygen-to-carbon ratio (O:C) on the X-axis and its hydrogen-to-carbon ratio (H:C) on the Y-axis. Families of molecules (e.g., lipids, terpenoids) generally cluster together on a VKD ([Brockman et al. 2018](https://link.springer.com/article/10.1007/s11306-018-1343-y)). VKDs provide a brief visulaization of a samples chemical composition.

vkmz has three modes to read several types of LC-MS data. `tabular` and `w4m-xcms` mode annotate features by calculating their neutral mass and attempts to match the feature's neutral mass to a lookup table of known mass-formulas pairs within a user given error range. `tabular` mode reads a "generic" fromat of LC-MS data that can be adapted to most LC-MS data processing workflows. `w4m-xcms` mode works with data processed by [Workflow4Metabolomics' XCMS Galaxy tools](https://workflow4metabolomics.org/raw-data-pre-processing-with-xcms). `formula` mode is similar to `tabular` mode, but reads pre-annotated molecular formulas to generate the VKD.

The primary output of vkmz is an interactive VKD webpage and a tabular file with annotation and VKD information. Output can also be saved as JSON and SQL.

## Installation

vkmz can be installed directly as a Python package. Alternatively, vkmz can be installed with [BioConda](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/vkmz) or to a [Galaxy Server](#).

### Python package installation:

```
# download vkmz
#   manually from: https://github.com/HegemanLab/vkmz/releases
#   or with git clone:
git clone https://github.com/HegemanLab/vkmz

# enter the directory of vkmz
cd vkmz/

# install (may require sudo)
python3 setup.py install

# run
vkmz --help
```

## Input Data

### Data without Molecular Formulas

vkmz is built to parse LC-MS data from a single tabular file or from Workflow4Metabolomics' XCMS (W4M-XCMS) tabular files. W4M-XCMS, a Galaxy tool, is the primary input type which vkmz was built around.

*Tabular* mode requires a single tabular file as input and  must include the columns "sample_name", "polarity", "mz", "rt", and "intensity". Each row represents a feature. Polarity refers to voltage polarity. Optionally a "charge" column can exist if software, such as CAMERA, has annotated the charge of features.
  - See [test-data/tabular.tabular](test-data/tabular.tabular) for tabular input example

*W4M-XCMS* mode requires the sample metadata, variable metadata, and data matrix tabular files generated with W4M-XCMS. Alternatively, these files can be annotated by W4M-CAMERA first.
  - See [test-data/datamatrix.tabular](test-data/datamatrix.tabular), [test-data/sampleMetadata.tabular](test-data/sampleMetadata.tabular), and [test-data/variableMetadata.tabular](test-data/variableMetadata.tabular) for W4M-XCMS data input example

In either mode, polarity values should be either "positive" or "negative".

If feature charge annotation is present, features without charge information will be removed. If CAMERA annotation is present, only monoisotopic features will be kept. An argument flag (`--impute-charge`) can be set to disable removing features without charge annotation. Users should be wary of false results when using this non-default option.

### Data with Molecular Formulas

Alternatively, *Formula* mode allows vkmz to read molecular structures from the input data. This is useful for workflows involving annotation software.
  - See [test-data/annotation.tabular](test-data/annotation.tabular) for annotated tabular input example

## Useage

### Galaxy

If you are using vkmz in Galaxy, please see the wrapper for useage info.

### Command Line Interface

#### Quick Start

```
# Tabular mode
vkmz tabular --input test-data/tabular.tabular --output foo --error 10

# W4M-XCMS mode
vkmz w4m-xcms -xd test-data/datamatrix.tabular -xv test-data/variableMetadata.tabular -xs test-data/sampleMetadata.tabular -o foo -e 10 --impute

# Formula mode
vkmz formula -i test-data/annotation.tabular -o foo
```

#### Help Menu

Add `--help` to a command to learn argument options.
```
$ vkmz --help
usage: vkmz [-h] {tabular,w4m-xcms,formula} ...

positional arguments:
  {tabular,w4m-xcms,formula}
                        Select mode:
    tabular             Tabular data mode
    w4m-xcms            W4M-XCMS data mode
    formula             Annotated molecular formula mode

optional arguments:
  -h, --help            show this help message and exit
```

Specific modes also have --help info:
```
$ vkmz tabular --help
usage: vkmz tabular [-h] --input INPUT --error [ERROR] --output [OUTPUT]
                    [--json] [--sql] [--metadata] [--database [DATABASE]]
                    [--prefix [PREFIX]] [--polarity {positive,negative}]
                    [--neutral] [--alternate] [--impute-charge]

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        Path to tabular file.
  --error [ERROR], -e [ERROR]
                        Mass error of MS data in parts-per-million
  --output [OUTPUT], -o [OUTPUT]
                        Specify output file path
  --json, -j            Set JSON flag to save JSON output
  --sql, -s             Set SQL flag to save SQL output
  --metadata, -m        Set metadata flag to save argument metadata
  --database [DATABASE], -db [DATABASE]
                        Define path to custom database of known formula-mass
                        pairs
  --prefix [PREFIX]     Define path prefix to support files ("d3.html" and
                        database directory)
  --polarity {positive,negative}, -p {positive,negative}
                        Set flag to force polarity of all features to
                        positive or negative
  --neutral, -n         Set flag if input data contains neutral feature mass
                        instead of mz
  --alternate, -a       Set flag to keep features with multiple predictions
  --impute-charge, --impute
                        Set flag to impute "1" for missing charge information
```
