# LeisHLA

This pipeline enables the parsing, standardization and analysis of NetMHC(II)pan T cell epitope prediction data.  
More specifically, this pipeline was applied on leishmaniasis-associated HLA alleles. However, this pipeline can be applied to any possible pathogen for which HLA disease-associations have been identified.

## Requirements
This pipeline requires a Linux or OSX environment, however, it can also be called from Windows Subsystem for Linux.  
It was designed and tested for R version ≥ 4.0.3.  
Installation of R package dependencies are handled within the different R scripts, however, if this fails these can be manually installed. A list with the dependencies of each script are provided within these scripts.

## Usage
First call the Epitope_Standardisation_Pipeline.sh script on raw prediction data like this:  
./Epitope_Standardisation_Pipeline.sh [-c MHC CLASS 1 OR 2] [-i INPUT FILE(S)] [-o OUTPUT FILE] [-d DISTRIBUTION FILE] [-s] [-r REFERENCE FILE]  

with:  
-c  MHC CLASS  Specify the MHC class, provide either a 1 or a 2.  
-i  INPUT FILE(S)  Provide the names of the input files containing your NetMHC(II)pan results separated by a comma.  
-o  OUTPUT FILE  Provide the name of the file you want to save the parsed results to.  
-d  DISTRIBUTION  Specify an output name to generate a plot with the binding affinity distribution of your data. When this option is specified, it is useful to append an identifier such as the species name, ending with an underscore, to the start of the input file name.  
-s  STANDARDIZE  Provide this option if you want to standardize.  
-r  REFERENCE FILE  Provide the name of the reference data file with which to standardize the prediction data or with which to generate a distribution plot.  
  
Once both protective- and risk-associated HLA alleles have been standardized, the Epitope_Analysis_Pipeline.sh can be called from the commandline to generate some preliminary results, like:  
./Epitope_Analysis_Pipeline.sh [-t THRESHOLD] [-P PROTECTIVE FILE] [-R RISK FILE] [-o OUTPUT FILE] [-d DISTRIBUTION]  

with:  
-t  THRESHOLD  Provide a numerical threshold with which to identify the strong-binding epitopes  
-p  PROTECTIVE FILE  Provide the name of the inputfile containing your standardized NetMHC(II)pan results for the Protective-associated alleles  
-r  RISK FILE  Provide the name of the inputfile containing your standardized NetMHC(II)pan results for the Risk-associated alleles  
-o  OUTPUT FILE  Provide the name of the file you want to save the statistical output to.  
-d  DISTRIBUTION  Specify an output name to generate a plot with the relative binding affinity distribution of your data. When this option is specified, it is useful to append an identifier such as the species name, ending with an underscore, to the start of the input file name.  
