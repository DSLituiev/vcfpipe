#!/bin/sh
####################################
# TO DO : option to select all SNPs
####################################
# get the input file name:
FILE1=$1; 
FILE2="$2";
if [ $# -eq 1 ]  ||  [ -z "$2" ]
then
# construct the output file name by replacing the extension: ".vcf" > "-ems.vcf"  :
   FILE2=$(echo $1 | sed 's/.vcf/-ems.vcf/g'); 
fi

if ! [ -s $1 ]
then
    echo "The input file does not exist or has zero-size! Exiting" 1>&2
    echo `ls "$1"1` 1>&2
    exit 1;
fi

SCRPATH=$(dirname $(readlink -f $0));
# do the conversion with an 'awk' script 
$SCRPATH/select_ems_snps.awk $FILE1 > $FILE2;

# replace the delimiter from "\t" to ";"  :
# | sed  's/\t/;/g' > $FILE2;
# print the report to the command line window  :
echo "====================================================================="
echo -n $FILE1 "\t>  has been converted to >\t" $FILE2 "\n" ;

# count number of lines (SNPs)  :
numLinesOld=$(cat < $1| sed '/^\s*#/d'|wc -l );
numLinesNew=$(cat < $FILE2| sed '/^\s*#/d'|wc -l ); # |/^\s*/d
echo -n "the old file contains\t"$numLinesOld"\tSNPs\n";
echo -n "the new file contains\t"$numLinesNew"\tSNPs\n";
echo -n "proportion of EMS SNPs\t"  $(echo "scale=2; ((100 * $nln / $nlo ))" |bc) "\t%\n";
echo "====================================================================="
