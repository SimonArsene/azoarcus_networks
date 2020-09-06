# -*- coding: utf-8 -*-
"""
Created on Tue May 23 16:26:15 2017

@author: Simon
"""

import pandas as pd
import sys
sys.path.append("C:/Users/smnar/az_network_droplets/utils")
import dpBC_tools as az
import argparse
from itertools import izip

parser = argparse.ArgumentParser(description="Concatenate both reads in one file and rename reads according to sample barcode.")
parser.add_argument("R1", type=str, help="read 1 file")
parser.add_argument("R2", type=str, help="read 2 file")
parser.add_argument("O", type=str, help="output file")
parser.add_argument("--sample_list", dest="sp_list", type=str, help="file containing sample barcode sequences (must have spBC_seq and spBC fields)")
parser.add_argument("-distance", dest="distance", type=str, help="either levenshtein or hamming (default = hamming)")
parser.add_argument("-N", type=int, default=1000000000, help="total number of reads")

args = parser.parse_args()

if args.sp_list != None:
    df = pd.read_excel(args.sp_list)
    spBC_list = list(df["spBC_seq"])
    spBC_ix = list(df["spBC"])
else:
    spBC_list = ["XXXXXX"]
    spBC_ix = [0]
    
if args.distance == None:
    args.distance = "hamming"
if args.distance == "levenshtein":
    use_lev = 1
elif args.distance == "hamming":
    use_lev = 0
else:
    raise ValueError("%s is not a valid distance metric: use levenshtein or hamming instead."%args.distance)

n = 0
bar = 0

with open(args.O, "w") as f:
    c = 0
    i = 0
    spBC = -1
    rcd = ["@id","seq","+","quality"]
    for (r1, r2) in izip(open(args.R1), open(args.R2)):
        if i == 0:
            assert r1.split(" ")[0] == r2.split(" ")[0]
            rcd[0] = "@"+str(c)
            i += 1
        elif i == 1:
            r1 = r1.strip()
            r2 = r2.strip()
            spBC = az.approx_search(str(r1[:6]), spBC_list, 1, use_lev=use_lev)
            rcd[0] += "," + str(spBC)
            rcd[1] = r1 + r2
            i += 1
        elif i == 2:
            i += 1
        elif i == 3:
            r1 = r1.strip()
            r2 = r2.strip()
            rcd[3] = r1 + r2
            i = 0
            c += 1
            f.write("\n".join(rcd)+"\n")        
            n += 1
        if n/(args.N*0.01) > bar:
            sys.stdout.write("\rConcatenating reads|"+"#"*(bar)+" "*(100-bar)+"| %s per cent\t%s reads treated"%(bar, n))
            sys.stdout.flush()
            bar += 1
