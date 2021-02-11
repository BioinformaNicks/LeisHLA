#!/bin/bash

# This script parses the results from NetMHCpan/NetMHCIIpan to fit a tabular data format that R or Python can read in
# Future versions will include direct calls to R or Python scripts

# This function defines the usage statement of this script.
usage() {
  echo "Usage: ${0} [-c MHC CLASS 1 OR 2] [-i INPUT FILE(S)] [-o OUTPUT FILE] [-d DISTRIBUTION FILE] [-s] [-r REFERENCE FILE]"
  echo "Parse your NetMHC(II)pan Epitope Prediction results and then stardardize these results according to a reference."
  echo " -c  MHC CLASS  Specify the MHC class, provide either a 1 or a 2."
  echo " -i  INPUT FILE(S)  Provide the names of the input files containing your NetMHC(II)pan results separated by a comma."
  echo " -o  OUTPUT FILE  Provide the name of the file you want to save the parsed results to."
  echo " -d  DISTRIBUTION  Specify an output name to generate a plot with the binding affinity distribution of your data.
     When this option is specified, it is useful to append an identifier such as the species name, ending with an underscore, to the start of the input file name."
  echo " -s  STANDARDIZE  Provide this option if you want to standardize."
  echo " -r  REFERENCE FILE  Provide the name of the reference data file with which to standardize the prediction data or with which to generate a distribution plot." 
  exit 1
}

if [[ ! $@ =~ ^\-.+ ]];
  then
    usage
fi

# This block checks whether the necessary options have been provided, and reads in the parameters for each option.
while getopts 'c:i:o:d:sr:' opt; do
  case $opt in
    c) 
       CLASS=${OPTARG}
       ;;
	i) 
       set -f
       IFS=, 
       INPUTFILES=($OPTARG) 
       ;;
    o) 
       OUTPUTFILE=($OPTARG)
	   ;;
	d)
	   DISTRIBUTION="yes"
	   DISTOUTPUT=${OPTARG}
	   ;;
    s)
	   STANDARDIZE="yes"
	   ;;
	r)
	   HUMAN=${OPTARG}
	   ;;
    ?) 
       usage
       ;;
  esac
done

# If -s is provided, a reference will always need to be provided as well by using -r and specifying a reference file. If this is not the case, usage + exit.
if [[ "${STANDARDIZE}" =~ ^(yes|y)$ ]] && [[ ${HUMAN} == "" ]];
  then
    usage
fi

# If -d is provided, a reference will always need to be provided as well by using -r and specifying a reference file. If this is not the case, usage + exit.
if [[ "${DISTRIBUTION}" =~ ^(yes|y)$ ]] && [[ ${HUMAN} == "" ]];
  then
    usage
fi

# -c -i or -o can not be empty. If they are empty, usage + exit. 
if [[ ${CLASS} == "" ]] || [[ ${INPUTFILES} == "" ]] || [[ ${OUTPUTFILE} == "" ]];
  then
    usage
fi

# the -c option (MHC CLASS) option needs to be set to either MHC class 1 or class 2. Otherwise exit.
if [ ${CLASS} -ne 1 ] && [ ${CLASS} -ne 2 ];
  then
    usage
fi

# This if block will parse NetMHCpan results when MHC Class is set to 1
if [ ${CLASS} -eq 1 ];
  then
    HEADERLINE=`cat ${INPUTFILES[0]} | head -n100 | grep 'Pos' | sed 's/MHC/HLA/g' | sed 's/Aff(nM)/Affinity/g' | sed 's/\s//' | sed 's/ \{1,\}/ /g' | awk '{print $2, $3, $4, $11, $13, $16}'`
    echo ${HEADERLINE} > ${OUTPUTFILE}

    for inputfile in "${INPUTFILES[@]}";
      do
        echo "Currently parsing: ${inputfile}"
	cat ${inputfile} | grep -v '#' | grep -v '^-' | grep -v 'nearest' | grep -v 'Number of' | grep -v 'Pos' | sed '/^[[:space:]]*$/d' | awk '{print $2, $3, $4, $11, $13, $16}' >> ${OUTPUTFILE}
      done
    echo "The epitope prediction data has been parsed and saved in ${OUTPUTFILE}"
fi

# This if block will parse NetMHCIIpan results when MHC Class is set to 2
if [ ${CLASS} -eq 2 ];
  then
    HEADERLINE=`cat ${INPUTFILES[0]} | head -n25 | grep 'Pos' | sed 's/MHC/HLA/g' | sed 's/Affinity(nM)/Affinity/g' | sed 's/\s//' | sed 's/ \{1,\}/ /g' | awk '{print $2, $3, $5, $7, $9, $12}'`
    echo ${HEADERLINE} > ${OUTPUTFILE}

    for inputfile in "${INPUTFILES[@]}";
      do
        echo "Currently parsing: ${inputfile}"
        cat ${inputfile} | grep -v '#' | grep -v '^-' | grep -v 'Number of' | grep -v 'Pos' | sed '/^[[:space:]]*$/d' | awk '{print $2, $3, $5, $7, $9, $12}' >> ${OUTPUTFILE}
      done
    echo "The epitope prediction data has been parsed and saved in ${OUTPUTFILE}"
fi

# This block will generate a distribution plot of the binding affinity by calling an R script
if [[ "${DISTRIBUTION}" =~ ^(yes|y)$ ]] && [[ ${HUMAN} != "" ]];
  then
    echo "The binding affinity distribution will now be plotted. This procedure will take a moment. Go get some coffee."
    Rscript --vanilla ./Scripts/Distributions.R ${OUTPUTFILE} ${HUMAN} ${DISTOUTPUT}
	echo "The binding affinity distribution plot has been saved to ${DISTOUTPUT}"
elif [[ "${DISTRIBUTION}" =~ ^(yes|y)$ ]] && [[ ${HUMAN} == "" ]];
  then
    usage
fi

# This block will perform the standardization of the epitope prediction data by calling an R script
if [[ "${STANDARDIZE}" =~ ^(yes|y)$ ]] && [[ ${HUMAN} != "" ]];
  then
    echo "Your data will now be standardized. This might take a moment. Go get some coffee."
    Rscript --vanilla ./Scripts/Standardisation.R ${HUMAN} ${OUTPUTFILE}
    echo "The epitope prediction data has been standardized according to the reference, and saved in ${OUTPUTFILE}"
elif [[ "${STANDARDIZE}" =~ ^(yes|y)$ ]] && [[ ${HUMAN} == "" ]];
  then
    usage
fi
