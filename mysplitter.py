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

# debugging constants
WHOLE_LOOP_TIME = time.time()
ITER_TO_LIST_TIME = time.time()

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

infq1 = sys.argv[1]
infq2 = sys.argv[2]

prefix = sys.argv[3]


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

ITER_TO_LIST_TIME = time.time() - ITER_TO_LIST_TIME

n_cpus = psutil.cpu_count()
print 'number of cpus: {}'.format(n_cpus)

def split(sq1list,sq2list, proc_num):
    print "process initializing", multiprocessing.current_process()

    # used this to test cpu usage, 4 cpus are at 100%
    #x = 0
    #while x < 1e100:
    #    x +=1

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
            fo = str(proc_num) + '_' + ci + '_' + bcd + '.fastq'
            fo = os.path.join('output', fo)

            if fo in FILE_DICT:
                outfile = FILE_DICT[fo]
                outfile.write(rec1)
            else:
                outfile = open(fo, "w+")
                FILE_DICT[fo] = outfile
                outfile.write(rec1)

            FILE_IO_TIME += time.time() - io_start_time

        i = i + 1

        if i%10000==0:
            elapsed_time = time.time() - start_time
            print str(i) + " has been processed, assigned rate: " + str(1-unassigned/100000.0) + ", consumed time: {}".format(elapsed_time)
            # print('memory info: {}'.format(process.memory_info()[0]))
            unassigned = 0





fq1list = return_dict[1]
fq2list = return_dict[2]

len1 = len(fq1list)
fq1list_a = fq1list[:len1*1/4]
fq1list_b = fq1list[len1*1/4:len1*2/4]
fq1list_c = fq1list[len1*2/4:len1*3/4]
fq1list_d = fq1list[len1*3/4:]

len2 = len(fq2list)
fq2list_a = fq2list[:len2*1/4]
fq2list_b = fq2list[len2*1/4:len2*2/4]
fq2list_c = fq2list[len2*2/4:len2*3/4]
fq2list_d = fq2list[len2*3/4:]




process1 = multiprocessing.Process(target=split,
                                   args=(fq1list_a, fq2list_a, 1))
process2 = multiprocessing.Process(target=split,
                                   args=(fq1list_b, fq2list_b, 2))
process3 = multiprocessing.Process(target=split,
                                   args=(fq1list_c, fq2list_c, 3))
process4 = multiprocessing.Process(target=split,
                                   args=(fq1list_d, fq2list_d, 4))
process1.start()
process2.start()
process3.start()
process4.start()
process1.join()
process2.join()
process3.join()
process4.join()

WHOLE_LOOP_TIME = time.time() - WHOLE_LOOP_TIME
print('WHOLE_LOOP_TIME, ITER_TO_LIST_TIME: {}, {}'.format(WHOLE_LOOP_TIME, ITER_TO_LIST_TIME))

# closing files
for k,v in FILE_DICT.iteritems():
    v.close()



