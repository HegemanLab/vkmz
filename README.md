## vkmz version 0.002

This code is in the process of being written. Original code from [VanKrevelenLocal](https://github.com/HegemanLab/VanKrevelenLocal). Code is based on the original ideas of Stephen Brockman and Eric Roden.

Here after releases will be used and development will be done in a branch.

### Basic usage

This program can take centroided mzXML files and create Van Krevelen diagrams using the plot.ly library.

Convert `.mzXML` files to `.csv` formatted for vkmzDriver:
```
python vkmzMZML.py -i test-data/example.mzXML -t 15
```
In this case an input (`-i`) mzXML file is chosen and a threshold is set. In this case every feature that has an intensity less than 15% of the maximum feature is removed.

Identify features and plot:
```
python vkmzDriver.py -i test-data/example.mzXML.csv
```

> Note that the lookup may not identify the true formula. The seen mass to charge ratio is compared to a list of chemicals with a known mass to charge to ratio. The first chemical with a known mass to charge ratio within a range from a seen mass is taken as the formula.

### Advanced Driver Options

##### PPM Error

PPM error can be set with the `--error` or `-e` argument.

The default error is 5 PPM. To set it to, say, 10 use:
```
python vkmzDriver.py -i test-data/example.mzXML.csv --error 10
```

##### Ratio Ouput

Identified compounds and their ratio an be saved as output. This information is saved as a csv file which can be read in a spreadsheet program.

The output file will be saved in the folder that vkmzDriver is run from with a filename similar to `ratios-20171222140540.csv`.

To add ouput declare the output argument:
```
python vkmzDriver.py -i test-data/example.mzXML.csv --output
```
 
