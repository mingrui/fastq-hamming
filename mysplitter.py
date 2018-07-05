#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  4 13:36:27 2018

@author: kevinmu
"""
import os
import sys
import gzip
import itertools
import time
from Bio import SeqIO
import shutil
import psutil
from hamming_cython_solution import hamming_loop

# debugging constants
FILE_IO_TIME = 0
CALC_TIME = 0
WHOLE_LOOP_TIME = 0

start_time = time.time()

process = psutil.Process(os.getpid())
print('process memory resident set size: {}'.format(process.memory_info().rss))

# this file dict will remember opened files and keep the files open
FILE_DICT = {}

if os.path.exists('output'):
    shutil.rmtree('output')
os.mkdir('output')

bcds = {}
with open("input/tang_barcode.txt",'r') as fi:
    fi.readline()
    for l in fi:
        l = l.strip().split('\t')
        bcds[l[1]] = l[0]

infq1 = 'input/test1.1.fastq.gz'#sys.argv[1]
infq2 = 'input/test1.2.fastq.gz'#sys.argv[2]

prefix = 'test'#sys.argv[3]

'''
fh1 = gzip.open(infq1,'r')
fh2 = gzip.open(infq2,'r')
fq1iter = SeqIO.parse(fh1,'fastq')
fq2iter = SeqIO.parse(fh2,'fastq')

iter_start_time = time.time()
for rec1, rec2 in itertools.izip(fq1iter, fq2iter):
    pass
print('iter elapsed time: {}'.format(time.time() - iter_start_time))
fh1.close()
fh2.close()

WHOLE_LOOP_TIME += time.time() - start_time
start_time = time.time()
'''
fh1 = gzip.open(infq1,'r')
fh2 = gzip.open(infq2,'r')
fq1iter = SeqIO.parse(fh1,'fastq')
fq2iter = SeqIO.parse(fh2,'fastq')

i = 0
unassigned = 0
for rec1, rec2 in itertools.izip(fq1iter, fq2iter):
    cal_start_time = time.time()

    bcd = str(rec2.seq[0:8])
    #print(bcd)
    umi = str(rec2.seq[8:16])
    ci = None
    if bcd in bcds:
        ci = bcds[bcd]
    else:
        for bar in bcds:
            di = hamming_loop(bcd, bar)
            if di < 3:
                bcd = bar
                ci = bcds[bcd]
                break

    CALC_TIME += time.time()-cal_start_time

    if ci is None:
        unassigned += 1   
    else:
        io_start_time = time.time()

        fo = prefix + '_' + ci  + '_' + bcd + '.fastq'
        fo = os.path.join('output', fo)

        if fo in FILE_DICT:
            outfile = FILE_DICT[fo]
            SeqIO.write(rec1, outfile,'fastq')
        else:
            outfile = open(fo,"w+")
            FILE_DICT[fo] = outfile
            SeqIO.write(rec1, outfile,'fastq')

        FILE_IO_TIME += time.time() - io_start_time
        
    i = i + 1

    if i%10000==0:
        elapsed_time = time.time() - start_time
        print '\n' + str(i) + " has been processed, assigned rate: " + str(1-unassigned/100000.0) + ", consumed time: {}".format(elapsed_time)
        WHOLE_LOOP_TIME += elapsed_time
        start_time = time.time()
        print('memory info: {}'.format(process.memory_info()[0]))
        unassigned = 0

print('WHOLE_LOOP_TIME, CALC_TIME, FILE_IO_TIME, IO_%_TIME: {}, {}, {}, {}'.format(WHOLE_LOOP_TIME, CALC_TIME, FILE_IO_TIME, FILE_IO_TIME/WHOLE_LOOP_TIME))

# closing files
for k,v in FILE_DICT.iteritems():
    v.close()

fh1.close()
fh2.close()

