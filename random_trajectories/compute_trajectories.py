from ruffus import *
import multiprocessing 
import argparse
import sys
import numpy as np
from random import sample
import datetime

parser = argparse.ArgumentParser(description="Compute network growth trajectoriesries")
parser.add_argument("N", type=int, help="Number of trajectories per set of values of the parameters")
parser.add_argument("nb_points", type=int, help="Number of values of link density to test")
#parser.add_argument("prefix", type=str, help="prefix common to all chunks")
#parser.add_argument("nb_chunks", type=int, help="number of chunks")
#parser.add_argument("position", type=int, help="position of the first linker")
parser.add_argument("-L", dest="list_spec", type=int, nargs='+', help="List of alphabet size")
#parser.add_argument("-distance", dest="distance", type=str, help="either levenshtein or hamming (default = hamming)")
#parser.add_argument("-N", type=int, default=1000000000, help="total number of reads")

args = parser.parse_args()
if args.list_spec == None:
    args.list_spec = [4]

def return_list_G(nb_spec):
    list_G = [(i, j) for i in range(nb_spec) for j in range(nb_spec)]
    return list_G

def transform_nt(x, list_G):
    if x == None:
        return x
    elif type(x) == list:
        return "".join(map(str, [int(g in x) for g in list_G]))
    else:
        return [list_G[i] for i in xrange(len(list_G)) if x[i] == '1']
    
def network_matrix(nt, rates, list_G):
    if type(nt) == str: nt = transform_nt(nt, list_G)
    n = len(nt)
    M = np.zeros((n, n))
    for i in xrange(n):
        for j in xrange(n):
            M[i,j] = rates[(nt[j][0], nt[i][1])]
    return M

def degree_all_nodes(x, norm, deg, rates, list_G):
    M = network_matrix(x, rates, list_G)
    if norm == 1:
        if M.sum() == 0:
            return np.ones(M.shape[0])*(1./x.count("1"))
        if deg == "in":
            return M.sum(1)/M.sum()
        else:
            return M.sum(0)/M.sum()
    else:
        if deg == "in":
            return M.sum(1)
        else:
            return M.sum(0)
        
def compute_perturbation(before, after, rates, list_G):
    y_before = degree_all_nodes(before, 1, 'in', rates, list_G)
    y_after = degree_all_nodes(after, 1, 'in', rates, list_G)
    x_before = transform_nt(before, list_G)
    x_after = transform_nt(after, list_G)
    y_before = {g:y_before[x_before.index(g)] for g in x_before}
    y_after = {g:y_after[x_after.index(g)] for g in x_after}
    common = list(set(y_before.keys()).intersection(set(y_after.keys())))
    common = transform_nt(transform_nt(common, list_G), list_G)    
    y = np.array([y_after[c] for c in common])
    if y.sum() > 0: 
        y = y/y.sum()
    else: 
        y = np.array([1./len(y)]*len(y))
    x = np.array([y_before[c] for c in common])
    if x.sum() > 0: 
        x = x/x.sum()
    else: 
        x = np.array([1./len(x)]*len(x))
    return np.abs(y-x).sum()

def random_trajectory(start, end, starting_nt, rates, list_G):
    X = [sample(starting_nt, 1)[0]]
    Y = [0]
    for i in range(end-start):
        ix = sample(np.arange(len(X[-1]))[(np.array(map(int, list(X[-1]))) == 0)], 1)[0]
        nt = list(X[-1])
        nt[ix] = "1"
        X.append("".join(nt))
        Y.append(Y[-1] + compute_perturbation(X[-2], X[-1], rates, list_G))
    return X, Y

def n_trajectories(start, end, N, rates, list_G):
    ix = [sorted(np.random.randint(0, len(list_G), size=start)) for _ in xrange(N*10)]
    starting_nt = ["".join(map(str, [1 if i in ix[j] else 0 for i in range(len(list_G))])) for j in xrange(N*10)]
    res = [[], []]
    for _ in xrange(N):
        X, Y = random_trajectory(start, end, starting_nt, rates, list_G)
        res[0].append(X)
        res[1].append(Y)
    return res

starting_files = ["network_trajectories_alph_%s.txt"%i for i in args.list_spec]

n = 0
bar = 0

@transform(starting_files, suffix(".txt"), ".csv")
def ruffus_compute_trajectories(input_file, output_file):
    global n, bar
    nb_spec = int(input_file.split(".")[0].split("_")[-1])
    with open(output_file, "w") as output:
        list_G = return_list_G(nb_spec)
        for k in [int(x*nb_spec**2) for x in np.linspace(0, 1, args.nb_points)]:
            for _ in xrange(args.N):
                ix = sample(range(len(list_G)), min(k, len(list_G)))
                rates = [1. if j in ix else 0. for j in xrange(len(list_G))]
                rates = {list_G[j]:rates[j] for j in xrange(len(list_G))}
                res = n_trajectories(2, min(nb_spec**2, 100), 1, rates, list_G)
                x = res[-1][0]
                nt_x = res[0][0]
                x = ["%.3f"%y for y in x] 
                if len(x) < 100: x += ["-"]*(100 - min(nb_spec**2, 100) + 1)
                output.write(",".join(map(str, [nb_spec, k] + x + ["_".join(nt_x), "".join([str(int(rates[g])) for g in list_G])]))+"\n")
                n += 1
                if n*len(args.list_spec)/(args.N*args.nb_points*len(args.list_spec)*0.01) > bar:
                    sys.stdout.write("\r|"+"#"*(bar)+" "*(100-bar)+"|")
                    sys.stdout.flush()
                    bar += 1
    
if __name__ == '__main__':
    try:
        pipeline_run(multiprocess = multiprocessing.cpu_count()-2, verbose = 2)
    except Exception, e:
        print e.args
