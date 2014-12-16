#!/usr/bin/awk -f

/^($|[:space:]*#)/ && NR==FNR {  # copy the VCF file
  print $0;
  next;
 } 
 
NR==FNR {  
printf( "Chr%s\n", $0 );
}

