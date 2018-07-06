import numpy as np
from itertools import islice
import multiprocessing
import time
from hamming_cython_solution import hamming_loop
import os


class Batch:
    def __init__(self, iterable, condition=(lambda x:False), limit=None):
        self.iterator = iter(iterable)
        self.condition = condition
        self.limit = limit
        try:
            self.current = next(self.iterator)
        except StopIteration:
            self.on_going = False
        else:
            self.on_going = True

    def group(self):
        yield self.current
        # start enumerate at 1 because we already yielded the last saved item
        for num, item in enumerate(self.iterator, 1):
            self.current = item
            if num == self.limit or self.condition(item):
                break
            yield item
        else:
            self.on_going = False

    def __iter__(self):
        while self.on_going:
            yield self.group()

def split(sq1list, sq2list, bcds, prefix, proc_num):
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
        outfile = open(k, 'a')
        outfile.write("\n".join(v))
        outfile.close()
