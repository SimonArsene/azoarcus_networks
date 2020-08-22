# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 17:57:27 2016

@author: Simon
"""

from Bio import SeqIO
import sys
import os
import azoarcus_seq_tools as sq
from ruffus import *
import argparse
import multiprocessing

parser = argparse.ArgumentParser(description="Extract Azoarcus meta data from sequences")
parser.add_argument("folder", type=str, help="absolute path to folder")
parser.add_argument("prefix", type=str, help="prefix common to all chunks")
parser.add_argument("nb_chunks", type=int, help="number of chunks")
parser.add_argument("nb_indexes", type=int, help="number of indexes")
parser.add_argument("start_read2", type=int, help="position of start of read 2")
parser.add_argument("-distance", dest="distance", type=str, help="either levenshtein or hamming (default = hamming)")
parser.add_argument("-N", type=int, default=1000000000, help="total number of reads")

args = parser.parse_args()

if args.distance == None:
    args.distance = "hamming"
if args.distance == "levenshtein":
    use_lev = 1
elif args.distance == "hamming":
    use_lev = 0
else:
    raise ValueError("%s is not a valid distance metric: use levenshtein or hamming instead."%args.distance)

# os.chdir(args.folder)

starting_files = ["%s%s.fastq"%(args.prefix, i) for i in range(args.nb_chunks)]
n = 0
bar = 0

@transform(starting_files, suffix(".fastq"), "_meta.csv")
def extract_meta_rcd(input_file, output_file):
    global n, bar
    with open(output_file, "w") as output:
        
        c = 0        
        
        for rcd in SeqIO.parse(input_file, "fastq"):
            c += 1
            n += 1
            meta = sq.get_meta_az(rcd, args.start_read2, args.nb_indexes, use_lev=use_lev)
            output.write(",".join(meta) + "\n")
            
            if n*args.nb_chunks/(args.N*0.01) > bar:
                sys.stdout.write("\rExtracting Azoarcus meta information |"+"#"*(bar)+" "*(100-bar)+"| %s per cent %s reads treated"%(bar, n))
                sys.stdout.flush()
                bar += 1
                    

if __name__ == '__main__':
    try:
        pipeline_run(multiprocess = multiprocessing.cpu_count()-2, verbose = 2)
    except Exception, e:
        print e.args
