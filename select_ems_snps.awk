#!/usr/bin/awk -f
BEGIN{ OFS = "\t"; }
(/^($|[:space:]*#)/ || (length($4)==1 && length($5)==1) && ( ($4 ~ /G/ && $5 ~ /A/ )|| ($4 ~ /C/ && $5 ~ /T/ ) ) ) {
printf $0;
printf "\n";
}
