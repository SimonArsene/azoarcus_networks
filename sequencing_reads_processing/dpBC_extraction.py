# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 11:52:40 2016

@author: Simon
"""

from Bio import SeqIO
import os
import sys
import dpBC_tools as az
from collections import Counter
from ruffus import *
import multiprocessing 
import argparse

parser = argparse.ArgumentParser(description="Extract indexes sequences from reads using linkers")
parser.add_argument("folder", type=str, help="absolute path to folder")
parser.add_argument("prefix", type=str, help="prefix common to all chunks")
parser.add_argument("nb_chunks", type=int, help="number of chunks")
parser.add_argument("position", type=int, help="position of the first linker")
parser.add_argument("-L", dest="linkers_list", type=str, nargs='+', help="linkers sequences")
parser.add_argument("-distance", dest="distance", type=str, help="either levenshtein or hamming (default = hamming)")
parser.add_argument("-N", type=int, default=1000000000, help="total number of reads")


args = parser.parse_args()

if args.linkers_list == None:
    args.linkers_list = ["TTCG", "TGAC", "ACCA", "CAAC"]

if args.distance == None:
    args.distance = "hamming"
if args.distance == "levenshtein":
    use_lev = 1
elif args.distance == "hamming":
    use_lev = 0
else:
    raise ValueError("%s is not a valid distance metric: use levenshtein or hamming instead."%args.distance)

os.chdir(args.folder)
starting_files = ["%s%s.fastq"%(args.prefix, i) for i in range(args.nb_chunks)]

n = 0
bar = 0

@transform(starting_files, suffix(".fastq"), ".txt")
def ruffus_extract_droplet_barcodes(input_file, output_file):
    global n, bar

    linkers_list = args.linkers_list

    # Number of barcodes to be extracted
    nb_BC = len(linkers_list) - 1
    
    # Output files
    output = []
    for i in range(nb_BC):
        output.append(open(input_file.split(".")[0] + "_dpBC_%s."%(i) + input_file.split(".")[1], "w"))

    # Counter which we call c
    c = Counter()

    # Position where the first linker is supposed to be found
    start = args.position
    
    for rcd in SeqIO.parse(input_file, "fastq"):
    
        c["total"] += 1
        n += 1

        # List of the linkers positions to be filled
        positions = []

        for i in range(len(linkers_list)):
            p = az.approx_pos(linkers_list[i], str(rcd.seq), start + i*20, 2, 1, use_lev=use_lev)
            positions.append(p)
            if p != -1:
                c[i] += 1

        for i in range(nb_BC):
            if positions[i] != -1 and positions[i+1] != -1:
                SeqIO.write(rcd[positions[i]+4:positions[i+1]], output[i], "fastq")
            elif positions[i] != -1:
                SeqIO.write(rcd[positions[i]+4:positions[i]+20], output[i], "fastq")
            else:
                SeqIO.write(rcd[:0], output[i], "fastq")
    
        if n*args.nb_chunks/(args.N*0.01) > bar:
            sys.stdout.write("\rExtracting droplet barcodes |"+"#"*(bar)+" "*(100-bar)+"| %s per cent %s reads treated"%(bar, n))
            sys.stdout.flush()
            bar += 1
             
    for f in output:
        f.close()
    
    with open(output_file, "w") as f:
        f.write(input_file+"\n")
        f.write(repr(c))
    
if __name__ == '__main__':
    try:
        pipeline_run(multiprocess = multiprocessing.cpu_count()-2, verbose = 2)
    except Exception, e:
        print e.args
