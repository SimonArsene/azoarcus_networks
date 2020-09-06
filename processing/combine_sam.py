# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 13:47:38 2016

@author: Simon
"""

from itertools import izip
from ruffus import *
import argparse
import multiprocessing 
import os
import sys


parser = argparse.ArgumentParser(description="Convert sam file to csv for aligned indexes with MQ filter")
parser.add_argument("folder", type=str, help="absolute path to folder")
parser.add_argument("prefix", type=str, help="prefix common to all chunks")
parser.add_argument("nb_chunks", type=int, help="number of chunks")
parser.add_argument("nb_indexes", type=int, help="number of indexes")
parser.add_argument("-N", type=int, default=1000000000, help="total number of reads")

args = parser.parse_args()

os.chdir(args.folder)
chunks = [args.prefix + str(i) for i in range(args.nb_chunks)]

starting_files = []

for i in range(args.nb_chunks):
    starting_files.append(["%s_dpBC_%s.sam"%(chunks[i], j) for j in range(args.nb_indexes)])
    
n = 0
bar = 0
    
@transform(starting_files, suffix(".sam"), ".csv")
def combine_dpBC(input_chunks, output_file):
    global n, bar
    with open(input_chunks[0].split(".")[0][:-2] + ".csv", "w") as output:
        c = 0
        if args.nb_indexes == 3:
            for l in izip(open(input_chunks[0]), open(input_chunks[1]), open(input_chunks[2])):
                
                seq_id = []
                dpBC = []
                ll = [""]*3
                
                for i in range(3):
                    ll[i] = l[i].split("\t")
                    seq_id.append(ll[i][0])
                    if int(ll[i][4]) >= 20:
                        dpBC.append(ll[i][2])
                    else:
                        dpBC.append("*")
    
                output.write(",".join(seq_id + ["#".join(dpBC)]) + "\n")
                
                c += 1
                n += 1
                if n*args.nb_chunks/(args.N*0.01) > bar:
                    sys.stdout.write("\rCombining SAM files |"+"#"*(bar)+" "*(100-bar)+"| %s per cent %s reads treated"%(bar, n))
                    sys.stdout.flush()
                    bar += 1 
                    
        elif args.nb_indexes == 4:
            for l in izip(open(input_chunks[0]), open(input_chunks[1]), open(input_chunks[2]), open(input_chunks[3])):
                
                seq_id = []
                dpBC = []
                ll = [""]*4
                
                for i in range(4):
                    ll[i] = l[i].split("\t")
                    seq_id.append(ll[i][0])
                    if int(ll[i][4]) >= 20:
                        dpBC.append(ll[i][2])
                    else:
                        dpBC.append("*")
    
                output.write(",".join(seq_id + ["#".join(dpBC)]) + "\n")
                
                c += 1
                n += 1
                if n*args.nb_chunks/(args.N*0.01) > bar:
                    sys.stdout.write("\rCombining SAM files |"+"#"*(bar)+" "*(100-bar)+"| %s per cent %s reads treated"%(bar, n))
                    sys.stdout.flush()
                    bar += 1
                
  
if __name__ == '__main__':
    try:
        pipeline_run(multiprocess = multiprocessing.cpu_count()-2)
    except Exception, e:
        print e.args
