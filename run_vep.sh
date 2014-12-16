#!/bin/bash
SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")

# set the species
SPECIES='arabidopsis_thaliana';


echo "====================================================================="
echo "VARIANT EFFECT PREDICTOR";
echo "species: " $SPECIES;
# get the input file name:
FILE1=$1; 
# construct the output file name by replacing "vcf" > "-annotation.vcf"  :
# FILEtab=$(echo $FILE1 | sed 's/.vcf/-tab.vcf/g'); 

FILE2="$2";
if [ -z "$FILE2" ]
then
   FILE2=$(echo $FILE1 | sed 's/.vcf/-annotation.vcf/g'); 
# else 
fi

# runs perl VEP script
# sed {'/^\s*$/d'} $FILE1 | sed {'s/;/\t/g'}  > $FILEtab;

# perl ./variant_effect_predictor/variant_effect_predictor.pl --cache  --vcf --species $SPECIES --host gramenedb.gramene.org --port 3306 --user anonymous --pass gramene -i $FILE1 -o $FILE2 --force_overwrite --sift=b

# perl ./variant_effect_predictor/variant_effect_predictor.pl --cache --offline -db_version 73 --vcf --species $SPECIES -i $FILE1 -o $FILE2 --force_overwrite --sift=b 

# old#73
# perl ./variant_effect_predictor/variant_effect_predictor.pl --vcf -db_version 73 --cache --offline --species  $SPECIES  -i $FILE1 -o $FILE2 --force_overwrite


perl "${SCRIPTPATH}"/variant_effect_predictor/variant_effect_predictor.pl --vcf -db_version 22 --cache --offline --species  $SPECIES  -i $FILE1 -o $FILE2 --force_overwrite


# -- fork 4 # run in parallel
# --regulatory 

# rm $FILEtab;

### WORKS! 
# perl variant_effect_predictor.pl -i attab_example.vcf --cache --species arabidopsis_thaliana --host gramenedb.gramene.org --port 3306 --user anonymous --pass gramene -o out0.vcf


# perl ./variant_effect_predictor/variant_effect_predictor.pl --cache -i $FILE1.vcf -o $FILE2.vcf  -fork 6  --regulatory --vcf --species arabidopsis_thaliana
