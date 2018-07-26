# VKMZ version 1.3.1

VKMZ is a metabolomics prediction and vizualization tool which creates van Krevelen diagrams from mass spectrometry data. A van Krevelen diagram (VKD) plots a molecule on a 2D scatterplot based on the molecule's oxygen to carbon ratio (O:C) against it's hydrogen to carbon ratio (H:C). Classes of metabolites cluster together on a VKD [0]. Plotting a complex mixture of metabolites on a VKD can be used to briefly convey untargeted metabolomics data.

VKMZ attempts to predict a molecular formula for each feature in LC-MS data. Each feature's mass is compared to a database of known formula masses. A prediction is made when a known mass is within the mass error range of an feature's uncharged (neutral) mass. A binary search algorithm is used to quickly make matches. Heristically generated databases for labeled and unlabeled metabolites are included [1]. VKMZ finds all predictions for an observed mass within the mass error. The prediction with the lowest delta (absolute difference between an feature's neutral mass and the predicted mass) is plotted. Features without predictions are discarded. Outputed is saved as a tabular and html file.

This software works best with, accurate, high resolution LC-MS data. A well calibrated LC-MS is essential for correct predictions. It is best to emperically derive mass error etiher from the data or from data using the same methods and spiked standards. Using low resolution data will result in false positive predictions, especially for large mass metabolites.

VKMZ can be used as a command line tool or on the Galaxy web platform [2]. A Galaxy wrapper for VKMZ is maintatined in this repository. VKMZ was developed on the Workflow4Metabolomics version of Galaxy [3].

## Using VKMZ command line

### Input modes

VKMZ has two input modes:
  1. `xcms` mode reads features from XCMS data
  2. `tsv` mode reads a specially formatted tabular file

Select a mode by declaring it as the first argument to `vkmz.py`.

> **Example:**
> ```
> python vkmz.py xcms [other parameters]
> ```

Different modes allow different parameters.

### Required parameters

#### xcms mode

xcms mode requires three tabular files generated by XCMS:
  * `--data-matrix [XCMS_DATA_MATRIX_FILE]`
  * `--sample-metadata [XCMS_SAMPLE_METADATAFILE]`
  * `--variable-metadata [XCMS_VARIABLE_METADATAFILE]`

##### xcms mode example:
```
python vkmz.py xcms --data-matrix test-data/datamatrix.tabular --sample-metadata test-data/sampleMetadata.tabular --variable-metadata test-data/variableMetadata.tabular [other parameters]
```

#### tsv mode

tsv mode requires a tabular file of a specific format as input:
  * `--input [TSV_FILE]`

The first five columns of the input tabular file must be:
>| sample_id | polarity | mz | rt | intensity |
>|-----------|----------|----|----|-----------|


#### All modes

Mass error of LC-MS in parts-per-million:
  * `--error [PPM_ERROR_NUMBER]`
    * It is critical to set the mass error correctly

Output name:
  * `--output [FILENAME]`
    * A `.tsv` and `.html` file will be generated by VKMZ with the given filename

### Optional parameters

Database:
  * `--database [DATABASE_FILE_PATH]`
    * Default is BMRB's monoisotopic heuristically generated database
    * Path is relative to `--directory`

Directory:
  * `--directory [TOOL_PATH]`
    * Explicitly define tool directory
    * Paths are relative if unset
    * Affects database and web page template paths

Forced Polarity:
  * `--polarity [positive|negative]`
    * Set all features to have either a positive or negative polarity
    * Overrides input files polarity information
    * Do not use this parameter on data containing both polarities

Neutral:
  * `--neutral`
    * Using this flag disables charged mass adjustment
    * Without this flag VKMZ adjusts a feature mass by adding or removing that mass of a proton based on the features charged polarity

Unique:
  * `--unique`
    * Remove features with multiple predictions from output

## Special thanks to

Adrian, Art, Eric, Jerry, Kevin, Renata, Stephen, Tim, and Yuan.

## Citations

0. Brockman et al. [doi:10.1007/s11306-018-1343-y](https://doi.org/10.1007/s11306-018-1343-y)
1. Hegeman et al. [doi:10.1021/ac070346t](https://doi.org/10.1021/ac070346t)
2. [Galaxy Project](https://galaxyproject.org/)
3. [Workflow4Metabolomics](http://workflow4metabolomics.org/)
4. Smith et al. [doi:10.1021/ac051437y](https://www.ncbi.nlm.nih.gov/pubmed/16448051)
