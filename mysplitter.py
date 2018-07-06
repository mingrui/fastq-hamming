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
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import shutil
import psutil
import multiprocessing

# cython
from util import Batch, split

# input
infq1 = sys.argv[1]
infq2 = sys.argv[2]
prefix = sys.argv[3]
process_num = sys.argv[4]

ALL_TIME = time.time()

# folder setup
if os.path.exists('output'):
    shutil.rmtree('output')
os.mkdir('output')

# bar code
print('')
print('Read Barcode...')
bcds = {}
with open("input/tang_barcode.txt",'r') as fi:
    fi.readline()
    for l in fi:
        l = l.strip().split('\t')
        bcds[l[1]] = l[0]

print('Convert iterator to list...')
# convert iter to list
fh1 = gzip.open(infq1,'r')
fh2 = gzip.open(infq2,'r')
fq1iter = FastqGeneralIterator(fh1)
fq2iter = FastqGeneralIterator(fh2)

batch_1= Batch(fq1iter, 1e7)
batch_2= Batch(fq2iter, 1e7)

for b1, b2 in itertools.izip(batch_1, batch_2):
    loop_time = time.time()

    print('\nStart split...')
    n_cpus = psutil.cpu_count()
    print('number of cpus: {}'.format(n_cpus))
    n = int(process_num)
    print('Numer of processes: {}'.format(n))

    to_list_t = time.time()
    fq1list = list(b1)
    fq2list = list(b2)
    to_list_t = time.time() - to_list_t
    print('to_list_time: {}'.format(to_list_t))

    len1 = len(fq1list)
    len2 = len(fq2list)
    jobs = []
    for i in range(n):
        fq1list_i = fq1list[len1*(i)/n:len1*(i+1)/n]
        fq2list_i = fq2list[len2*(i)/n:len2*(i+1)/n]
        process_i = multiprocessing.Process(target=split,
                                            args=(fq1list_i, fq2list_i, bcds, prefix, i))
        jobs.append(process_i)
    for j in jobs:
        j.start()
    for j in jobs:
        j.join()

    print('loop time: {}'.format(time.time() - loop_time))

fh1.close()
fh2.close()

ALL_TIME = time.time() - ALL_TIME
print('ALL_TIME: {}'.format(ALL_TIME))




