# vkmz version 0.005

vkmz can take centroided mzXML files and create Van Krevelen diagrams using the plot.ly library.

This code is in the process of being written. Original code comes from [VanKrevelenLocal](https://github.com/HegemanLab/VanKrevelenLocal). Code is based on the original ideas of Stephen Brockman and Eric Roden.

## Converting mzXML files for vkmz

Currently vkmz uses a program named `vkmzMZXML` to convert a `mzXML` into a `TSV` file.

To convert a file run something like:
```
python vkmzMZML.py -i test-data/example.mzXML -o test-data/example.mzXML.tsv -t 15
```
This will create a new TSV file in the directory of the mzXML file.

In this example two arguments are used:
  * The input argument `-i` specified a path to an mzXML file.
  * The output argument `-o` specifies the output file name.
  * Threshold (`-t`) was set to 15.
    * This means that features with an intensity less than 15% of the maximum intensity will be filtered out.

## Identify features in TSV file

The TSV file can be used to search for formula identification:
```
python vkmzIdentifier.py -i test-data/example.mzXML.tsv -o test-data/example.mzXML.identified.tsv
```

This will generate a TSV file that can be read in a spread sheet program or generate graphs with vkmzPlotter.

The output file will be saved in the folder that vkmzDriver is run from with a filename similar to `identified-20171222140540.tsv`.

> Note that the lookup may not identify the true formula. The seen mass to charge ratio is compared to a list of chemicals with a known mass to charge to ratio. The first chemical with a known mass to charge ratio within a range from a seen mass is taken as the formula.

### vkmzIdentifier has several useful arguments:

#### PPM error can be set with the `--error` or `-e` argument.

The default error is 5 PPM. To set it to, say, 10 use:
```
python vkmzIdentifier.py -i test-data/example.mzXML.tsv --error 10
```

#### Multiprocessing

vkmz' multiprocessing feauture divides the search time for masses in the lookup table by the number of CPU cores your comptuer has.

To turn on multiprocessig add the `-m` argument:
```
python vkmzIdentifier.py -i test-data/example.mzXML.tsv -m
```

#### Output TSV File Naming

A more descriptive name can be specified with the `-o` argument:
```
python vkmzIdentifier.py -i test-data/example.mzXML.tsv -o corylus
```
This would create a filename similar to `corylus.tsv`.

## Plotting Identified TSV File
```
python vkmzPlotter.py -i test-data/ratios-20171228162504.tsv
``` 
vkmzPlotter requires a TSV file input. By default it will make a 2D scatter plot with dots, for each identified mass, of a uniform size and with color based on it's retention time.

#### Plot Types
Select other plot types with the `-p` arugment:
```
python vkmzPlotter.py -i test-data/ratios-20171228162504.tsv -p 3d
``` 

This will create a 3D graph with N:C as the Z-axis.

#### Size and Size Algorithm
Base size can be adjusted with `-s [number]` and size algorithm can be selected with `-a [numnber`

For instance, to use 15 as the base size and to use algorithm 2 use:
```
python vkmzPlotter.py -i test-data/ratios.tsv -p 3d -s 15 -a 2
``` 

There are currently three size algorithms. All of them use base_size (`-s`) as a variable:
  0. `size = base_size`
  1. `size = base_size+4*base_size*feature+peak/(highest_peak-lowest_peak)`
  2. `size = base_size+2*log(base_size*feature+peak/(highest_peak-lowest_peak)`
