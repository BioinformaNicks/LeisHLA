#!/bin/bash

# This script parses the results from NetMHCpan/NetMHCIIpan to fit a tabular data format that R or Python can read in
# Future versions will include direct calls to R or Python scripts

# This function defines the usage statement of this script.
usage() {
  echo "Usage: ${0} [-t THRESHOLD] [-P PROTECTIVE FILE] [-R RISK FILE] [-o OUTPUT FILE] [-d DISTRIBUTION]"
  echo "Analyse the Protective-associated alleles and Risk-associated alleles and generate statistical and graphical output."
  echo " -t  THRESHOLD  Provide a numerical threshold with which to identify the strong-binding epitopes"
  echo " -p  PROTECTIVE FILE  Provide the name of the inputfile containing your standardized NetMHC(II)pan results for the Protective-associated alleles"
  echo " -r  RISK FILE  Provide the name of the inputfile containing your standardized NetMHC(II)pan results for the Risk-associated alleles"
  echo " -o  OUTPUT FILE  Provide the name of the file you want to save the statistical output to."
  echo " -d  DISTRIBUTION  Specify an output name to generate a plot with the relative binding affinity distribution of your data.
     When this option is specified, it is useful to append an identifier such as the species name, ending with an underscore, to the start of the input file name."
  exit 1
}

if [[ ! $@ =~ ^\-.+ ]];
  then
    usage
fi

# This block checks whether the necessary options have been provided, and reads in the parameters for each option.
while getopts 't:p:r:o:d:' opt; do
  case $opt in
    t) 
       THRESHOLD=${OPTARG}
       ;;
	p) 
       PROTFILE=${OPTARG}
       ;;
	r) 
	   RISKFILE=${OPTARG}
	   ;;
    o) 
       OUTPUTFILE=${OPTARG}
	   ;;
	d)
	   DISTRIBUTION='yes'
	   DISTOUTPUT=${OPTARG}
	   ;;
    ?) 
       usage
       ;;
  esac
done

# -t or -p or -r or -o can not be empty. If they are empty, usage + exit. 
if [[ ${THRESHOLD} == "" ]] || [[ ${PROTFILE} == "" ]] || [[ ${RISKFILE} == "" ]] || [[ ${OUTPUTFILE} == "" ]];
  then
    usage
fi

if [[ ${DISTRIBUTION} == 'yes' ]] && [[ ${DISTOUTPUT} == "" ]];
  then
    usage
fi

echo "The analysis will now begin. This procedure will take a moment. Go get some coffee."
Rscript --vanilla ./Scripts/Analysis.R ${THRESHOLD} ${PROTFILE} ${RISKFILE} ${OUTPUTFILE} ${DISTOUTPUT}
echo "The results of the statistical test have been saved in ${OUTPUTFILE}"

if [[ ${DISTRIBUTION} == "yes" ]];
  then
    echo "A plot of the relative binding affinity distributions has been saved in ${DISTOUTPUT}"
fi
