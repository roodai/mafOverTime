# Minor Allele Frequency Over Time
The script multiMafCalc.py takes .ped files and corresponding .map files as input and calculates minor allele frequency (MAF). MAFs for SNPs of interest are then plotted versus time.

## Dependencies
NumPy\
Pyplot from Matplotlib

## Outline
- Parse .ped files and compute allele counts
- Cross reference to .map files and extract SNP IDs
- Calculate allele frequencies
- Get time and SNP IDs of interest from user
- Plot MAF for SNP IDs of interest

## Procedure
##### Usage
The script can be ran as such:

```shell
 python multiMafCalc.py test1.map test1.ped test2.map test2.ped test3.map test3.ped test4.map test4.ped
 ```
##### Input
The input files have to loosely follow the .ped format where nucleotides are numerically represented in a tab-separated values file from the 7th field onwards. The SNP ID list also has to follow the .map file format where the SNP ID is the second field in a tab-separated file. The SNP IDs have to occur in corresponding order to columns in the .ped files, where one SNP represents two nucleotides. SNP IDs have to be consistent through all .map inputs. Input .ped files are to be included in order of

The script asks the user for input as follows:
```text
Time between datasets test1.ped and test2.ped: 1
Time between datasets test2.ped and test3.ped: 2
Time between datasets test3.ped and test4.ped: 3
Provide SNP(s) of interest: rs3094315 rs6696609
```
The time values have to numerical and SNP(s) of interest have to be found in all the .ped files included for the input to be read successfully  and for plotting to be possible.

##### Output
The script outputs a png image where MAF is plotted against time, with error bars that represent standard deviation.

#### Alternative Usage
As mentioned above different types of input files can be used, as long as they adhere to the aforementioned constraints. Likewise different numerical conditions between datasets can be used, as long as a cumulative sum is a correct representation of the over all condition. The input text and the x-axis of the plot will be labeled as time, regardless of condition between datasets.
