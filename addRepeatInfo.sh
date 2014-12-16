#!/bin/bash

if [[ $# -lt 1 ]]
then
    echo "usage:" 1>&2
    echo -en "./""\033[31m""addRepeatInfo.sh""\033[m" "\033[32m""in_file_name [ out_file_name ] [ repeat_table ] " "\033[m" 1>&2
    echo " "
    exit 1;
fi

# get the input file name:
FILE1=$1;

if [[ $# -lt 2 ]]
then
    # construct the output file name by replacing the extension:
    FILE2=$(echo $FILE1 | sed 's/.vcf/-repfilt.vcf/g');
else
    FILE2=$2
fi

if [[ $# -lt 3 ]]
then
    REPEAT_TABLE='./reference/TAIR10-chr.fas.out'
else
    REPEAT_TABLE=$3
fi

# the second argument is the table with repeat information:
awk -f awkFilterMasked.awk $FILE1 $REPEAT_TABLE > $FILE2;

# print the report to the command line window  :
echo "====================================================================="
echo "                        REPEAT MASKING";
echo -en $FILE1 "\n\t>  has been converted to >\n\t" $FILE2 "\n" ;

# count number of lines (SNPs)  :
nlo=$(cat < $FILE1| sed '/^\s*#/d'|wc -l );
echo -en "the old file contains\t"$nlo"\tSNPs\n";

nlo=$(cat < $FILE2| sed '/^\s*#/d'|wc -l );
echo -en "the new file contains\t"$nlo"\tSNPs\n";

echo "====================================================================="
