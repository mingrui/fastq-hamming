#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  4 16:41:23 2018

@author: kevinmu

"""
import os
import shutil
import sys
import time
import gzip
import itertools
import distance
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# debugging constants
FILE_IO_TIME = 0
CALC_TIME = 0
WHOLE_LOOP_TIME = 0

infq1 = 'input/test1.1.fastq.gz'#sys.argv[1]
infq2 = 'input/test1.2.fastq.gz'#sys.argv[2]

prefix = 'test'#sys.argv[3]

if os.path.exists('output'):
    shutil.rmtree('output')
os.mkdir('output')

bcds = {}
with open("input/tang_barcode.txt",'r') as fi:
    fi.readline()
    for l in fi:
        l = l.strip().split('\t')
        bcds[l[1]] = l[0]
        
for bar in bcds:
    print "Processing " + bar + " " + bcds[bar] + "..."
    start_time = time.time()
    fh1 = gzip.open(infq1,'r')
    fh2 = gzip.open(infq2,'r')
    fq1iter = FastqGeneralIterator(fh1)
    fq2iter = FastqGeneralIterator(fh2)
    fo = prefix + '_' + bcds[bar]  + '_' + bar + '.fastq'
    fo = os.path.join('output', fo)
    wh = open(fo,"w")
    
    i,j=0,0
    st2= time.time()
    for (title1, seq1, qual1),(title2, seq2, qual2) in itertools.izip(fq1iter, fq2iter):
        i+=1
        if str(seq2[0:8])==bar or distance.hamming(str(seq2[0:8]), bar) < 3:#
            j+=1
            wh.write("@%s\n%s\n+\n%s\n" % (title1, seq1, qual1))
        if i%1000000==0:
            et2 = time.time() - st2
            print "Filtering " + bar + ", "+str(i) + " reads done, assigned rate " + str(j*1.0/i) +", taking " + str(et2) + " seconds"
            st2 = time.time()
    elapsed_time = time.time() - start_time
    print "Found " + str(j) + " reads in " + str(i) + " reads, consumed time {}".format(elapsed_time)
    WHOLE_LOOP_TIME += elapsed_time
    fh1.close()
    fh2.close()

print('WHOLE_LOOP_TIME, CALC_TIME, FILE_IO_TIME, IO_%_TIME: {}, {}, {}, {}'.format(WHOLE_LOOP_TIME, CALC_TIME, FILE_IO_TIME, FILE_IO_TIME/WHOLE_LOOP_TIME))

