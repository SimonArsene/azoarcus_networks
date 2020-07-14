# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 15:17:08 2016

@author: Simon Arsène
"""

import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy import stats
import networkx as nx
import az_tools_ as az
from itertools import permutations, combinations

list_genotypes = [M+N for M in "ACUG" for N in "ACUG"]
G = [M+N for M in "ACUG" for N in "ACUG"]

#Remove junction letter in a list of genotypes
def rmv_junction(x):
    if type(x) == list:
        return [y[0] + y[-1] for y in x]
    return x[0] + x[-1]

#Transform list of genotypes into network string and vice-versa
def transform(x):
    list_genotypes = [M+N for M in "ACUG" for N in "ACUG"]
    if x == None:
        return x
    elif type(x) == list:
        return "".join([str(int(g in x)) for g in list_genotypes])
    else:
        return [list_genotypes[i] for i in range(16) if x[i] == '1']

#Transform an integer to binary representation with 16 digits
get_bin = lambda x: format(x, 'b').zfill(16)

#Given a binary representation, returns a boolean mask
bin_to_bool = lambda x: x == '1'

#Vectorized version of the bin to bool function
mask_bin_to_bool = np.vectorize(bin_to_bool)

def network_matrix(nt, rates):
    """
    Returns the network matrix of a network nt
    Takes a list of interaction weigths (in the same order as the list of nodes) as the "rates" argument
    """
    if type(nt) == str: nt = transform(nt)
    else: nt = transform(transform(nt))
    n = len(nt)
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            M[i,j] = rates[list_genotypes.index(nt[j][0] + nt[i][-1])]
    return M
        
def asymptotic_solution(x, rates, centrality=0, base_16=0):
    """
    Returns the first-order model asymptotic solution for a network x
    Takes a list of interaction weigths (in the same order as the list of nodes) as the "rates" argument
    """
    if type(x) == list:
        x = transform(x)
        M = network_matrix(x, rates=rates)
    elif type(x) == np.ndarray:
        M = x
        x = None
    else:
        M = network_matrix(x, rates=rates)
    L, V = np.linalg.eig(M)
    kmax = np.real(L).argmax()
    if centrality: 
        return (np.abs(np.real(V[:,kmax])), transform(x))
    if base_16: 
        temp = np.real(V[:,kmax])/np.real(V[:,kmax]).sum()
        res = []
        for i in range(16):
            if list_genotypes[i] in transform(x):
                res.append(temp[transform(x).index(list_genotypes[i])])
            else: 
                res.append(0)
        return np.array(res)
    return (np.real(V[:,kmax])/np.real(V[:,kmax]).sum(), transform(x))

def az_solve(nt, X0, tf, normed=True):
    """
    MAY NOT BE USED IN FINAL STUDY
    Take a network structure x (list of genotypes or network structure string), initial conditions x0 and a final time-point tf
    Returns integrated solution (normalized or not) based on 1st order model
    """
    def dx(x, t, M):
        return np.dot(M, x)
        
    if type(nt) == list:
        nt = transform(nt)
        M = network_matrix(nt)
    else:
        M = network_matrix(nt)
        
    X = odeint(lambda x, t: dx(x, t, M), X0, np.linspace(0, tf, 1000))
    if normed:
        return az.norm_array(X)
    else:
        return X

#Returns a Networkx graph object for a network
#Takes a list of interaction weigths (in the same order as the list of nodes) as the "rates" argument
def graph(x, rates):
    if type(x) == list:
        return graph(transform(x), rates)
    else:
        return nx.from_numpy_matrix(network_matrix(x, rates).T, create_using=nx.DiGraph())

#MAY NOT BE USED 
#Draw a network (only Watson-Crick interactions) with network string or list of genotypes using networkx drawing routine
def az_draw(x):
    if type(x) == str:
        g = graph(x, WC=True)
        nodes = dict(enumerate(transform(x)))
        g = nx.relabel_nodes(g, nodes)
        nx.draw(g, arrows=True, with_labels=True, node_color="w", node_size=600)
    if type(x) == list:
        az_draw(transform(x))
       
#Computes first eigenvalue of a network
def growth(x):
    if type(x) == list:
        x = transform(x)
        M = network_matrix(x)
    elif type(x) == np.ndarray:
        M = x
        x = None
    else:
        M = network_matrix(x)
    L, V = np.linalg.eig(M)
    return np.real(L).max()


def degree(g, x, norm, deg, rates):
    """
    Computes the in or out degree of a node "g" in a network "x"
    Takes a list of interaction weigths (in the same order as the list of nodes) as the "rates" argument
    The "norm" argument is used to turn on or off the normalisation
    If the network has no link at all, returns 1./(number of nodes)
    """
    if type(x) == list: 
        x = transform(x)
    if norm == 1:
        if network_matrix(x, rates).sum() == 0:
            return 1./x.count("1")
        if deg == "in":
            return network_matrix(x, rates).sum(1)[transform(x).index(g)]/network_matrix(x, rates).sum()
        else:
            return network_matrix(x, rates).sum(0)[transform(x).index(g)]/network_matrix(x, rates).sum()
    else:
        if deg == "in":
            return network_matrix(x, rates).sum(1)[transform(x).index(g)]
        else:
            return network_matrix(x, rates).sum(0)[transform(x).index(g)]
    
def compute_diff(A, B):
    #Auxiliary function to what_mut function
    return [-1*(a=='1' and b=='0') + 1*(a=='0' and b=='1') for a, b in zip(list(B), list(A))]

def what_mut(A, B):
    #Returns mutation is required to go from network B to network A
    D = compute_diff(A, B)
    if D.count(0) == 14 and sum(D) == 0:
        return (G[D.index(-1)], G[D.index(1)])
    else:
        if D.count(0) == 15 and sum(D) == 1:
            return ("-", G[D.index(1)])
        elif D.count(0) == 15 and sum(D) == -1:
            return (G[D.index(-1)], "-")
        else:
            return "-"

#NOT USED 
def WC_matrix_in_order(x):
    #Returns a WC matrix for list of genotypes x in same order as given
    WC_pairs = ["AU", "UA", "CG", "GC"]
    n = len(x)
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            A[i,j] = (x[j][0] + x[i][1]) in WC_pairs
    return A.astype("int")

#NOT USED 
def compute_all_topologies(size_nt):
    list_nt = map(list, list(combinations(G, size_nt)))    
    list_topo = []
    list_nt_per_topo = []
    for i, s in enumerate(map(lambda x: set(["".join(map(str, list(WC_matrix_in_order(y).flatten()))) for y in permutations(x)]), list_nt)):
        if s in list_topo:
            list_nt_per_topo[list_topo.index(s)].append(list_nt[i])
        else:
            list_topo.append(s)
            list_nt_per_topo.append([])
            list_nt_per_topo[list_topo.index(s)].append(list_nt[i])
    return (list_nt_per_topo)

#USED?
def where(element, list_of_lists):
    for i, l in enumerate(list_of_lists):
        if element in l: return i
    return -1

"""
Functions below are dedicated to pertubation analysis

SOME CLEAN UP IS REQUIRED HERE

"""

def compute_pscore(before, after, data_dict, method, th=0, rates=None, with_mut=0, pseudo_count=10**-2):
    """
    Computes shuffling score between networks before and after either with the data or with the model (th parameter)
    with_mut parameter allows to take or not into account the node which is added/deleted/substituted
    """
    if th:
        y_before = asymptotic_solution(before, rates=rates)[0]
        y_after = asymptotic_solution(after, rates=rates)[0]  
    else:
        y_before = data_dict[before]
        y_after = data_dict[after]
    y_before = {g:y_before[transform(before).index(g)] for g in transform(before)}
    y_after = {g:y_after[transform(after).index(g)] for g in transform(after)}
    common = list(set(y_before.keys()).intersection(set(y_after.keys())))
    common = transform(transform(common))
    
    if before.count("1") == after.count("1"):
        mut = (set(y_before.keys()).difference(set(y_after.keys())), set(y_after.keys()).difference(set(y_before.keys())))
        mut = [list(m)[0] for m in mut]
    
    # With mutated node (not up to date)
    if with_mut == 1:
        if before.count("1") == after.count("1"):
            pscore = [abs(y_after[c] - y_before[c]) for c in common] + [abs(y_after[mut[1]] - y_before[mut[0]])]
            return np.mean(pscore)
        else:
            return -1
    #Without mutated node
    else:
        y = np.array([y_after[c] for c in common])
        if y.sum() > 0: y = y/y.sum()
        else: y = np.array([1./len(y)]*len(y))
        x = np.array([y_before[c] for c in common])
        if x.sum() > 0: x = x/x.sum()
        else: x = np.array([1./len(x)]*len(x))
        if method == "mean_abs":
            return [common, x, y, np.mean(np.abs(y-x))]
        elif method == "euclidian":
            return [common, x, y, np.sqrt(((y-x)**2).sum())]
        elif method == "KL":
            return [common, x, y, stats.entropy(y + pseudo_count, x + pseudo_count)]
        
def compute_epi(data_, method, also_th=0, rates=None):
    """
    Returns perturbation data
    """
    epi = []
    temp = data_[(data_.nb_az.isin(range(2, 17))) & (data_.az_list != "0"*16)]
    temp_dict = {x.az_list:x.az_gd for name, x in temp.iterrows()}
    list_nt = list(temp.az_list.values)
    S, D = az.sample_graph_matrix_from_list(list_nt)
    S = nx.from_numpy_matrix(S + D)
    for e in S.edges_iter():
        for i, j in [(0, 1), (1, 0)]:
            b = list_nt[e[i]]
            a = list_nt[e[j]]
            mut = what_mut(a, b)
            with_mut = 0
            order, y_before, y_after, pscore = compute_pscore(b, a, temp_dict, method=method, th=0, with_mut=with_mut)
            if also_th:
                order_th, yth_before, yth_after, pscore_th = compute_pscore(b, a, temp_dict, method=method, th=1, rates=rates, with_mut=with_mut)
            else:
                order_th, yth_before, yth_after, pscore_th = order, y_before, y_after, pscore
            n = min(a.count("1"), b.count("1"))
            if method == "mean_abs":
                prandom = np.abs(az.norm_array(np.random.rand(n)) - az.norm_array(np.random.rand(n)))
            elif method == "euclidian":
                prandom = np.sqrt(((az.norm_array(np.random.rand(n)) - az.norm_array(np.random.rand(n)))**2).sum())
            elif method == "KL":
                prandom = stats.entropy(az.norm_array(np.random.rand(n)), az.norm_array(np.random.rand(n)))
            what = -1
            if "-" in mut: what = mut.index("-")
            epi.append([b, a, mut, pscore, b.count("1"), pscore_th, what, prandom.mean(), prandom.sum(), y_before, y_after, yth_before, yth_after, order, order_th])
    epi = pd.DataFrame(epi, columns=["before", "after", "mutation", "pscore", "nt_size_before", "pscore_th", "what", "prand", "prand_sum", "y_before", "y_after", "yth_before", "yth_after", "order_sp", "order_th_sp"])
    return epi

def transform_epi(size, epi, rates):
    temp = epi[(epi.what == 0) & (epi.nt_size_before == size)]
    temp["mutation"] = temp.mutation.apply(lambda x: x[1])
    df = []
    for name, x in temp[:].iterrows():
        v = transform(x.before)
        n = x.before.count('1')
        res = []
        d_a_out = sum([rates[G.index(x.mutation[0] + g[1])] for g in v])
        sigma_G = network_matrix(x.before, rates=rates).sum()
        d_v_in = [degree(g, x.before, norm=0, deg='in', rates=rates) for g in v]
        va = [g for g in transform(x.before) if np.isclose(rates[G.index(x.mutation[0] + g[1])], 0) == False]
        va_star = [g for g in transform(x.before) if g not in va]
        y = x.y_before
        y_ = x.y_after
        if sum(y_) == 0:
            y_ = np.array([1./n]*n)
        else:
            y_ = az.norm_array(y_)
        sv = y_ - y
        for i in range(n):
            res.append([x.before.count("1"), x.before, x.after, x.mutation, v[i], d_v_in[i], d_a_out, sigma_G, rates[G.index(x.mutation[0] + v[i][1])], sv[i], 
                        (rates[G.index(x.mutation[0] + v[i][1])] - d_a_out*d_v_in[i]/sigma_G)/(sigma_G + d_a_out),
                        sum([degree(g, x.before, norm=0, deg='in', rates=rates) for g in va]), 
                        sum([degree(g, x.before, norm=0, deg='in', rates=rates) for g in va_star])])
        df = df + res
    df = pd.DataFrame(df, columns=["nt_size","before", "after", "mutation", "node", "deg_in_v", "deg_out_a", "sigma_g", "edge_a_v", "p_v", "p_v_th", "deg_in_va", "deg_in_va_star"])
    return df