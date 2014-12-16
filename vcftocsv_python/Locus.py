# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 14:18:38 2014

@author: Dmytro Lituiev
"""
import re
import sys

class Locus(object):
    outFields = ['chr', 'pos', 'rating', 'refAllele', 'altAllele',\
    'geneID', 'geneSO', \
    # , 'mutType',\
    'mutCDS', 'mutProt', 'mutPosCDS', 'mutPosProt', 'closestDist']    
    outCountFields = ['totCount', 'refCount', 'altCount']
    #===========
    # refCount = [0, 0];
    # altCount = [0, 0];
    # totCount = [0, 0];
    # refQuality = [0, 0];
    # altQuality = [0, 0];
    # altFrequency = [0.0, 0.0];
    #===========
    
    # mFlag = True; # True if the read comes from the mutant pool
    
    def __init__(self, line, SO_DICTIONARY, numOfSamples = 1, repeatInfoFlag = False):
        self.numOfSamples = numOfSamples
        self.repeatInfoFlag = repeatInfoFlag
        #===========
        self.refCount = [0]* self.numOfSamples;
        self.altCount = [0]* self.numOfSamples;
        self.totCount = [0]* self.numOfSamples;
        self.refQuality = [0]* self.numOfSamples;
        self.altQuality = [0]* self.numOfSamples;
        self.altFrequency = [0.0]* self.numOfSamples;
        #===========
    
        if not line:  
            self.chr = 0;    
            self.pos = 0;
            self.quality = 0;
            self.refAllele = "";
            self.altAllele = "";
            self.geneID = "";
            self.geneSO = "";
            self.mutType = "";
            self.mutCDS = "";
            self.mutProt = "";
            self.mutPosCDS = 0;
            self.mutPosProt = 0;
            self.repeatName = "";
            self.repeatType = "";
            self.rating = 0;
            self.closestDist = float("inf");
            return
        
        cols = line.split('\t')
        if (cols[0][0:3] == "Chr") or (cols[0][0:3] == "Chr"):
            self.chr = cols[0][3:]
        else:
            self.chr = cols[0]  # chromosome
        
        self.pos =   int(cols[1])    # Position
        self.refAllele = cols[3]     # ref allele    
        self.altAllele = cols[4]     # alt 
        
        if (cols[5] != "."):
            self.quality = cols[5]
        else:
            self.quality = 0
        
        ##
        def readSampleCounts(self, cols, ii):
            SUMMARY = cols[9+ii].split(":");
            # GT:DP:RO:QR:AO:QA:GL
            #  0  1  2  3  4  5  6
            # print('sample: %u'% ii, 'data: %s' % cols[9+ii], 'datalen %u' % len(cols[9+ii]), file=sys.stderr)
            if not (cols[9+ii]=='.') and len(cols[9+ii])>2:
                self.totCount[ii]   = int(SUMMARY[1])     # Depth
                self.refCount[ii]   = int(SUMMARY[2])     # Observations reference
                self.altCount[ii]   = int(SUMMARY[4])      # Observations alternative
                self.refQuality[ii] = SUMMARY[3]      # Quality reference
                self.altQuality[ii] = SUMMARY[5]      # Quality alternative
                
                self.altFrequency[ii] = float(self.altCount[ii]) / float(self.totCount[ii])
                
        def strtofloat(s):
            if s == '':
                return 0.0
            else:
                return float(s)
        def strtoint(s):
            if s == '':
                return 0
            else:
                return int(s)
        
        
        #ii = 0
        for ii in range(0,self.numOfSamples):
        # while ii<self.numOfSamples:            
            readSampleCounts(self, cols, ii)
        #    ii +=1
            
        ## gene annotation column    
        re1 = re.compile(";CSQ=")
        try:
                INFO = re1.split(cols[7])[1];
        except:
                print(" probably there is no annotation in this '.vcf' file, or its format is wrong."\
                "\n please use Ensemble VEP to add gene annotation\n______________________________ ", file=sys.stderr)
                raise
        genes = INFO.split(',')
        ratingHit = 0
        geneHit = ['0']; # genes[0];
        closestDist = 100000000;
        for g in genes:
            gd = g.split("|", 12)
            if (len(gd)>1) :
                descr = gd[4].split('&');                
                # print(gd, end = '\n', file=sys.stderr)
                currDist = strtoint(gd[-2])
                for dd in descr:
                    rating = SO_DICTIONARY.get(dd);
                    # print(rating, end = '\t')
                    if (rating!= None) and ( (rating[-1] > ratingHit) or ( (rating[-1] == ratingHit) and ( currDist < closestDist) )):
                        ratingHit = rating[-1];
                        gd[0] = SO_DICTIONARY.get(dd)[1]
                        geneHit = gd;
                        # print(closestDist, currDist, rating[-1] , ratingHit, end = '\n', file=sys.stderr)
                        closestDist = currDist
        if len(geneHit)<10:
            print( "geneHit entries are missing %s, %u" % (self.chr, self.pos) , file=sys.stderr)
            for i in range(len(geneHit)-1, 9): geneHit.append("")        


        self.geneSO     = geneHit[0][3:] # 'SO:' prefix is removed
        self.geneID     = geneHit[1]
        self.mutType    = geneHit[4]
        self.mutPosCDS  = geneHit[6]
        self.mutPosProt = geneHit[7]
        self.mutCDS     = geneHit[8]
        self.mutProt    = geneHit[9]
        self.closestDist   = geneHit[-2]
        self.rating     = ratingHit;
        
        ## repeat info
        
        if repeatInfoFlag or (len(cols)> 9+self.numOfSamples):
            REPEAT_INFO = cols[10+self.numOfSamples-1].split(";");
            self.repeatType = REPEAT_INFO[1];
            if  not (REPEAT_INFO[1] == "NO") :
                self.repeatName = REPEAT_INFO[0] # Repeat name
            else:
                self.repeatName = "" # Repeat name
        
    def printFields(self, SEP):
        for name in self.outFields:
            val = getattr(self, name)
            print(val, end = SEP)
        for ind in range(0, self.numOfSamples):
            for name in self.outCountFields:
                val = getattr(self, name)
                print(val[ind], end = SEP)
        if self.repeatInfoFlag:
            print(self.repeatName, self.repeatType, sep= SEP, end = SEP)                
        print('')
    
    def printFieldNames(self, SEP):
        print('num of samples: %u' % self.numOfSamples, file=sys.stderr)
        for name in self.outFields:
            print(name, end = SEP)
        for name in self.outCountFields:
            print('mt_'+name, end = SEP)
        if (self.numOfSamples > 1):
            for name in self.outCountFields:
               print('wt_'+name, end = SEP)        
        if self.repeatInfoFlag:
            print('repeatName', 'repeatType', sep= SEP, end = SEP)       
        print('')
        
    def getRowSqlite(self):
        ii = 0
        outList = [0]* ( len(self.outFields) + len(self.outCountFields)*self.numOfSamples -1)
        for name in self.outFields[1:]:
            outList[ii] = getattr(self, name)
            ii += 1
            
        for ind in range(0, self.numOfSamples):
            for name in self.outCountFields:
                outList[ii] = getattr(self, name)[ind]
                ii += 1
                
        return (self.chr, outList)
            
    def getRowNamesSqlite(self):
        ii = 0
        outList = [0]* ( len(self.outFields) + len(self.outCountFields)*self.numOfSamples -1)
        
        for name in self.outFields[1:]:
            outList[ii] = name
            ii += 1
        for name in self.outCountFields:
            outList[ii] = 'mt_'+name
            ii += 1
        if (self.numOfSamples > 1):
            for name in self.outCountFields:
               outList[ii] = 'wt_'+name
               ii += 1
               
        return outList
