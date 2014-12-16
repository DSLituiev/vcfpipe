#!/usr/local/bin/python3.4

import sys
from Locus import Locus
from readsodict import readsodict 
import re
# cmd line options: test500.vcf out.csv -a 1 -f 0.1

#####################################################
#####################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("inFile", 
                    help="input file (.vcf) with gene annotation")

parser.add_argument("outFile", nargs='?', type=str, default=r'',
                    help="output file (.csv); for stdout type '-' ")
                    
parser.add_argument("-s", "--soterms", default="SO_terms.csv",
                    help="a path to a .csv file with SO terms and prior definition")

parser.add_argument("-f", "--frequencyFilter", type=float, default=0.0,
                    help="lower threshold of frequency to output")

parser.add_argument("-a", "--altCountFilter", type=int, default=0,
                    help="output reads with alternative counts strictly more than the given value")

parser.add_argument("-t", "--totCountFilter", type=int, default=0,
                    help="output reads with total counts strictly more than the given value")            
   
parser.add_argument("-u", "--sumAltCountFilter", type=int, default=0,
                    help="output loci with alternative counts strictly more than the given value")
                    
parser.add_argument("-c", "--csvseparator", type=str, default= r';',
                    help="separator for the output .csv file [;]")     

args = parser.parse_args()
###############################################################################
if not args.outFile == r'-':
    if not args.outFile == r'':
        print('output file: %s' % args.outFile,  file=sys.stderr) 
        sys.stdout = open(args.outFile, 'w')
    else:
        print('output file: %s' %args.inFile.replace('.vcf','.csv'),  file=sys.stderr) 
        sys.stdout = open(args.inFile.replace('.vcf','.csv'), 'w')

###############################################################################
def applyFilter(args, loc):
    flag = False
    totAlt = 0
    for ii in range(0, len(loc.altCount)):
        totAlt += loc.altCount[ii]
        flag = flag or ((loc.altFrequency[ii] > args.frequencyFilter) \
            and (loc.altCount[ii] > args.altCountFilter) \
            and (loc.totCount[ii] > args.totCountFilter)) 
    flag = flag and totAlt > args.sumAltCountFilter
    # if flag:
    #     print('##########' , file=sys.stderr)
    return flag

#####################################################

SO_DICTIONARY = readsodict(args.soterms)

f = open(args.inFile)

# columnnames = f.readline();

commre = re.compile(r"^[ ]*#.*");

totLines = 0
skippedLines = 0
for line in f:
    if not (commre.match(line)):
        header = lastheaderline.split("\t");
        print("header: \n %s" % lastheaderline, file=sys.stderr, end='')
       
        numOfSamples = 0;
        repeatInfoFlag = False
        for sampleName in header[9:]:
            repeatInfoFlag = repeatInfoFlag or (sampleName.strip() == r'Repeat_Info')
            if not repeatInfoFlag:
                print("sample no %u " % numOfSamples, ":'"+sampleName.strip()+"'", sep = None, file=sys.stderr)
                numOfSamples += 1
            else:
                print('repeat info is available', sep = None, file=sys.stderr)
        # print('number of samples: %u' % numOfSamples , file=sys.stderr)
        #
        loc = Locus("", SO_DICTIONARY, numOfSamples, repeatInfoFlag);
        loc.printFieldNames(args.csvseparator)
        # process the line            
        prevLine = line
        totLines += 1
        loc = Locus(prevLine, SO_DICTIONARY, numOfSamples, repeatInfoFlag);
        if (applyFilter(args, loc)) :
            loc.printFields(args.csvseparator)
        else:
            skippedLines += 1
        #
        break
    else:
        lastheaderline = line
    
assert len(lastheaderline) > 0
print(lastheaderline, file=sys.stderr)
assert len(header) > 0
        
    # print(line, file=sys.stderr)
print('-----------------------------', file=sys.stderr)


for line in f:
        totLines += 1
        loc = Locus(line, SO_DICTIONARY, numOfSamples, repeatInfoFlag);
        if (applyFilter(args, loc)) :
            loc.printFields(args.csvseparator)
        else:
            skippedLines += 1
        
print( "out of %u lines," % totLines, "%u skipped" % skippedLines , "%u remained" % (totLines-skippedLines), file=sys.stderr)
     
f.close()

if not args.outFile == r'-':
    sys.stdout.close()
    
sys.stdout = sys.__stdout__
