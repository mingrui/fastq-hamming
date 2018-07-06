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
import distance
import multiprocessing

# cython
from hamming_cython_solution import hamming_loop

# input
infq1 = sys.argv[1]
infq2 = sys.argv[2]
prefix = sys.argv[3]

# debugging constants
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
TO_LIST_TIME = time.time()
# convert iter to list
fh1 = gzip.open(infq1,'r')
fh2 = gzip.open(infq2,'r')
fq1iter = FastqGeneralIterator(fh1)
fq2iter = FastqGeneralIterator(fh2)

def iter_to_list(fq_iter, return_dict, return_key):
    return_dict[return_key] = list(fq_iter)

manager = multiprocessing.Manager()
return_dict = manager.dict()
process1 = multiprocessing.Process(target=iter_to_list,
                                  args=(fq1iter, return_dict, 1))
process2 = multiprocessing.Process(target=iter_to_list,
                                   args=(fq2iter, return_dict, 2))
process1.start()
process2.start()
process1.join()
process2.join()
fh1.close()
fh2.close()

TO_LIST_TIME = time.time() - TO_LIST_TIME

print('Start split...')
n_cpus = psutil.cpu_count()
print('number of cpus: {}'.format(n_cpus))

def split(sq1list, sq2list, proc_num):
    print "process initializing", multiprocessing.current_process()
    file_dict = {}
    i = 0
    unassigned = 0
    FILE_IO_TIME = 0
    CALC_TIME = 0
    for (title1, seq1, qual1),(title2, seq2, qual2) in zip(sq1list, sq2list):
        start_time = time.time()

        cal_start_time = time.time()

        bcd = seq2[0:8]
        #print(bcd)
        umi = seq2[8:16]
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

            rec1 = "@%s\n%s\n+\n%s\n" % (title1, seq1, qual1)
            fo = prefix + '_' + ci + '_' + bcd + '.fastq'
            fo = os.path.join('output', fo)

            if fo in file_dict:
                file_dict[fo].append(rec1)
            else:
                file_dict[fo] = [rec1]

            FILE_IO_TIME += time.time() - io_start_time

        i = i + 1

        if i%10000==0:
            elapsed_time = time.time() - start_time
            print str(i) + " has been processed, assigned rate: " + str(1-unassigned/100000.0) + ", consumed time: {}".format(elapsed_time)
            unassigned = 0

    for k,v in file_dict.iteritems():
        outfile = open(k, 'w+')
        outfile.write("\n".join(v))
        outfile.close()

# split list to 4 section, run each section on single cpu

pool = multiprocessing.Pool()

fq1list = return_dict[1]
fq2list = return_dict[2]
len1 = len(fq1list)
len2 = len(fq2list)
n = 8
jobs = []
for i in range(n):
    fq1list_i = fq1list[len1*(i)/n:len1*(i+1)/n]
    fq2list_i = fq2list[len2*(i)/n:len2*(i+1)/n]
    process_i = multiprocessing.Process(target=split,
                                        args=(fq1list_i, fq2list_i, i))
    jobs.append(process_i)
# for j in jobs:
#     j.start()
# for j in jobs:
#     j.join()
pool.map(jobs)

ALL_TIME = time.time() - ALL_TIME
LOOP_TIME = ALL_TIME - TO_LIST_TIME
print('ALL_TIME, TO_LIST_TIME, LOOP_TIME: {}, {}, {}'.format(ALL_TIME, TO_LIST_TIME, LOOP_TIME))




