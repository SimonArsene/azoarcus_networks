# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 12:11:36 2016

@author: Simon
"""

import subprocess as sub
import argparse
import multiprocessing

parser = argparse.ArgumentParser(description="Identify extracted droplet barcode indexes by aligning them to reference with Bowtie2")

parser.add_argument("folder", type=str, help="absolute path to folder")
parser.add_argument("prefix", type=str, help="prefix common to all chunks")
parser.add_argument("nb_chunks", type=int, help="number of chunks")
parser.add_argument("nb_indexes", type=int, help="number of indexes")
parser.add_argument("-I", dest="indexes_list", type=str, nargs='+', help="Bowtie2 index for indexes")


args = parser.parse_args()

base_name = args.folder

input_files = []

for i in range(args.nb_chunks):
    for k in range(args.nb_indexes):
        input_files.append(base_name + "%s%s_dpBC_%s.fastq"%(args.prefix, i, k))

if args.indexes_list == None:
    indexes = ["D:/EvoEvo/indexes/new_indexes/B_index", "D:/EvoEvo/indexes/new_indexes/C_index", "D:/EvoEvo/indexes/new_indexes/D_index"]
else:
    indexes = args.indexes_list
        
with open(base_name + "%s_mapping_log.txt"%(args.prefix), "w") as log:
    
    for f in input_files:    
        
        index = indexes[int(f.split(".")[0][-1])]
        
        log.write("File:%s\nIndex:%s\n"%(f.split("/")[-1], index.split("/")[-1]))       
        
        sub.call(["bowtie2-align-s", 
        "-x", index, 
        "-U", f, 
        "-S", f.split(".")[0] + ".sam", 
        "-p", "%s"%(multiprocessing.cpu_count()-2), 
        "--no-hd", "--very-sensitive", "--reorder", "-t"], 
        stdout=log, stderr=log)
        
        print("File %s done"%(f.split("/")[-1]))
    
with open(base_name + "%s_mapping_log.txt"%(args.prefix)) as f, open(base_name + "%s_mapping_log_cleared.txt"%(args.prefix), "w") as ff:
    for l in f:
        if not l.startswith("W"):
            ff.write(l)