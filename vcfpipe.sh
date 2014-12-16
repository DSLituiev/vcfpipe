#!/bin/sh
####################################
# get the input file name:

USAGE="Usage: \n\033[31mvcfpipe \033[m  \033[32m file_name.vcf\033[m [\033[32m -c\033[m|\033[32m-chr \033[m: to append 'Chr' flag]"

if [ $# -eq 0 ]
then
   echo "$USAGE"
   exit 1;
fi

SCRPATH=$(dirname $(readlink -f $0));

echo "scrpath: $SCRPATH"

FILENAME=$(echo $1 | sed 's/.vcf//g'); 

if [ -z "$1" ]  # If the file name is empty
then
  echo "The file name is missing! Aborting..";
  exit 1
else
     if [ ! -f  "$(eval echo $1 )" ];
     then
       echo "The file \"" $1 "\" has not been found! Aborting...";
       exit 1
     else
       echo "Processing " $FILENAME;
     fi
fi
 
 
# CHR_FLAG="$2";
# if [! -z "$CHR_FLAG" ] # If the flag is not empty
# then
#   ./appendChr.awk $FILENAME.vcf > $FILENAME-$CHR_FLAG.vcf
# fi
TEMP=`getopt -o c --long chr -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables.
while true ; do
    case "$1" in
        -c|--chr) CHR_FLAG="-chr" ; shift ;;
        --) shift ; break ;;
        *) echo $USAGE ; exit 1 ;;
    esac
done


if [ ! -z "$CHR_FLAG" ] # If the flag is not empty
then
      echo "=====================================================================";
      echo "appending 'chr':\t" $FILENAME.vcf " > " $FILENAME$CHR_FLAG.vcf;
      $SCRPATH/append_chr.awk $FILENAME.vcf > $FILENAME$CHR_FLAG.vcf;
fi

       # filter to keep only EMS SNPs:
sh $SCRPATH/select_ems_snps.sh  $FILENAME$CHR_FLAG.vcf $FILENAME-ems.vcf;
	# CALLS  ./onlysnpsems.awk
	# OUTPUT $FILENAME-ems.vcf

        # run variant effect predictor
sh $SCRPATH/run_vep.sh $FILENAME-ems.vcf
	# CALLS ./variant_effect_predictor/variant_effect_predictor.pl
        # in the QUAL column the VEP replaces "0" -> "."

	# repeat masker
# sh $SCRPATH/addrepeatinfo.sh $FILENAME-ems-annotation.vcf
	# CALLS awk -f filterMasked.awk  ./data/HL7-ems-annotation.vcf araTha5.fa.out > ./data/HL7-ems-annotation-filt5.vcf

        # produce CSV output
python3.4 $SCRPATH/vcftocsv_python/trnslvcf.py $FILENAME-ems-annotation.vcf $FILENAME-ems-annotation.csv -a 1 \
     -f 0.05 -u 2 -s "$SCRPATH/vcftocsv_python/SO_terms.csv"

# sh $SCRPATH/vcftocsv.sh $FILENAME-ems-annotation-repfilt.vcf
	# CALLS  awk -f awktrnslvcf.awk ./data/ABD159-ems-annotation-repfilt.vcf | less

	# copy files
# find . -name '*-ems-annotation-repfilt.csv' | cpio -updm /media/Data/Documents/MATLAB/SeqMapping_20140218/data

# cp ./data/*.csv /media/Data/Documents/MATLAB/SeqMapping_20140218/data
 
