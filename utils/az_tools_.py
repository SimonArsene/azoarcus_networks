# from Bio import SeqIO
#import sys
#sys.path.append("D:/EvoEvo/dist/editdistance-0.2/")
#import editdistance as ed
import pandas as pd
import numpy as np
import re
#import swalign
from random import random
from itertools import chain
import os
import az_model_ as azm
from math import sqrt
import networkx as nx

WT_with_insert = ("GTGCCTTGCGCCGGGAAACCACGCAAGGAATGGTGTCAAATTCGGCGAAACCTAAGCGCCCGCNTGCAT"
                     + "CGGGCGTATGGCAACGCCGAGCCAAGCTTCGGCGCCNTGCAT"
                     + "GCGCCGATGAAGGTGTAGAGACTAGACGGCACCCACCTAAGGCACNTGCAT"
                     + "CGCTATGGTGAAGGCATAGTCCAGGGAGTGGCGAAAGTCACACAAACCGG")
                     
WT = ("GTGCCTTGCGCCGGGAAACCACGCAAGGGATGGTGTCAAATTCGGCGAAACCTAAGCGCCCGCCCGGGCGTATGGCAACGCCGAGCCAAGCTTCGGCGCCTGCGCCGATGAAGGT"
        + "GTAGAGACTAGACGGCACCCACCTAAGGCAAACGCTATGGTGAAGGCATAGTCCAGGGAGTGGCGAAAGTCACACAAACCGG")

G = [M+N for M in "ACUG" for N in "ACUG"]

def join_columns(df, sep, list_columns, new_column):
    
    for i in range(len(list_columns)):
        if i == 0:
            df[new_column] = df[list_columns[0]].map(str)
        else:
            df[new_column] = df[new_column].map(str) + sep + df[list_columns[i]].map(str)

def group(value, groups, dist_max):
    
    distances = []    
    for group in groups:
        distances.append(abs(value - group))
    m = min(distances)
    if m <= dist_max and distances.count(m) == 1:
        return groups[distances.index(m)]
    return -1
    
def norm_array_aux(x):
    #Rows of x = samples
    #Columns of x = observations
    return x.sum(1).repeat(x.shape[1]).reshape(x.shape)

def norm_array(x):
    if type(x) == list: x = np.array(x)
    if x.ndim > 1:
        y = x/norm_array_aux(x)
        return np.nan_to_num(y)
    else:
        if x.sum() > 0:
            return x/float(x.sum())
        else:
            return np.array([1./len(x)]*len(x))
    
def approx_pos(query, subject, position, window, dist_max, use_lev=0):
    
    assert (type(query) == str)
    assert (type(subject) == str)
    assert (type(position) == int)
    assert (type(window) == int)
    assert (position + len(query) <= len(subject))
    assert (len(subject) >= (len(query) + 2*window))
    
    lev = []
    pos_lev = []
    
    if position - window < 0: 
        min_i = 0
    else:
        min_i = position - window
    
    if position + window + len(query) > len(subject):
        max_i = len(subject) - len(query) + 1
    else:
        max_i = position + window + 1
    
    for i in range(min_i, max_i):
        pos_lev.append(i)
        if use_lev == 1:
            lev.append(levenshtein(query, subject[i:i+len(query)]))
        else:
            lev.append(hamming_distance(query, subject[i:i+len(query)]))
        
    if min(lev) <= dist_max and lev.count(min(lev)) == 1:
        return pos_lev[lev.index(min(lev))]
    else:
        return -1

def approx_search(query, search_list, distance, use_lev=0):
    
    assert (len(query) == len(search_list[0]))
    
    lev_dist = []
    
    for e in search_list:
        if use_lev == 1:
            lev_dist.append(levenshtein(query, e))
        else:
            lev_dist.append(hamming_distance(query, e))
        
    m = min(lev_dist)
    if m <= distance and lev_dist.count(m) == 1:
        return lev_dist.index(m)
    else:
        return -1

def frequency(data, bin, value):

    data_with_count = data.groupby([bin, value]).size().reset_index(name="count")
    freq = data_with_count.groupby(bin)["count"].apply(lambda x: x/float(x.sum())).reset_index(name=value + "_freq")   
    del freq["index"]
    freq = pd.concat([data_with_count, freq], axis=1)
    return freq

def heatmap(data, bin, value):
    """
    Takes a DataFrame, groups it according to value in "bin" and returns a "heatmap"
    """
    freq = frequency(data, bin, value)
    all_value = list(freq.groupby(value).size().reset_index(name="count")[value])
    nb_group = len(data.groupby(bin))

    tab = np.zeros((nb_group, len(all_value)), dtype='f')
    hm = pd.DataFrame(tab, columns=all_value)

    group_ID = pd.DataFrame(freq[bin].unique())
    group_ID.columns = ["group_ID"]
    
    size_group = data.groupby(bin).size().reset_index(name="group_size").get("group_size")

    id = 0
    for name, group in freq.groupby(bin):
        for f,s in zip(group[value + "_freq"], group[value]):
            if f > 0:
                hm[s][id] = f
        id += 1

    return pd.concat([group_ID, size_group, hm], axis=1)

def levenshtein(s, ss):
    #return ed.eval(s, ss)
    return old_levenshtein(s, ss)

def hamming_distance(a, b):
    assert len(a) == len(b)
    s, i = 0, 0
    n = len(a)
    while i < n:
        if a[i] != b[i]: 
            s += 1
        i += 1
    return s
"""
def align(a, b, dis=0, gep=-1, match=5, mismatch=-1, gap_penalty=-2):
    match = match
    mismatch = mismatch
    gap_penalty = gap_penalty
    gap_extension_penalty = gep
    
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring, gap_penalty=gap_penalty, gap_extension_penalty=gap_extension_penalty)
    
    al = sw.align(a, b)
    
    (aa, bb) = cigar_tr(a, b, al)
    ms = get_matching_string(aa, al)
    
    if dis != 0:
        print_aligned(aa, bb, ms, dis)
        
    return al
"""
def cigar_tr(a, b, al):
    if al.r_pos > al.q_pos:     
        ra = a[:al.r_pos]
        rb = "-"*(al.r_pos-al.q_pos) + b[:al.q_pos]
    elif al.r_pos < al.q_pos:     
        ra = "-"*(al.q_pos-al.r_pos) + a[:al.r_pos]
        rb = b[:al.q_pos]
    else:
        ra = ""
        rb = ""
    ca = al.r_pos
    cb = al.q_pos
    for i, sym in al.cigar:
        if sym == "M":
            ra += a[ca:ca+i]
            rb += b[cb:cb+i]
            ca += i
            cb += i
        elif sym == "I":
            ra += "-"*i
            rb += b[cb:cb+i]
            cb += i
        elif sym == "D":
            ra += a[ca:ca+i]
            rb += "-"*(i)
            ca += i
    ra += a[ca:] + "-"*(ca - al.r_end)
    rb += b[cb:] + "-"*(cb - al.q_end)
    return (ra,rb)

def print_aligned(a, b, ms, window=120):
    w = window
    n = len(a)/w

    for i in range(n):
        print(a[i*w:i*w+w])
        print(ms[i*w:i*w+w])
        print(b[i*w:i*w+w])
        #print("\n")
    
    print(a[n*w:])
    print(ms[n*w:])
    print(b[n*w:])
    
def get_cigar(ext_cigar):
    str_list = re.findall(r'\d+[A-Z]', ext_cigar)
    r = []
    for s in str_list:
        r.append((int(re.sub(r'[A-Z]', "", s)), re.sub(r'\d+', "", s)))
    return r

def get_matching_string(aa, al):
    cigar = get_cigar(al.extended_cigar_str)
    r = " "*max(al.r_pos, al.q_pos)
    for i, sym in cigar:
        if sym == "M":
            r += "|"*i
        elif sym in ["D", "I"]:
            r += " "*i
        elif sym == "X":
            r += "."*i
    r += " "*(len(aa) - len(r))
    return r
    
def find_seq_ordered(filename, target_list):
    seqID_list = []
    seq_list = []
    
    with open(filename, "r") as f:
        for i, l in enumerate(f):
            if i / 4 in target_list:
                if i % 4 == 0:
                    seqID_list.append(l.split(",")[0].strip("@"))
                elif i % 4 == 1:
                    seq_list.append(l.strip("\n"))

    return (seqID_list, seq_list)
    
    
def find_seq(handle, target_list):
    """
    Find the sequences of the reads whose ID is in the target_lists
    Take as an input a file handle, works on the normal handle rather than on BioPython parser (slower)
    """
    n = 0
    here = False
    seqID_list = []
    seq_list = []
    for line in handle:
        if here == True:
            seq_list.append(line.strip('\n'))
            n += 1
            here = False
            if n == len(target_list):
                break
        if here == False and line.startswith("@") and line.split(" ")[0].strip("@").strip('\n') in target_list:
            here = True
            seqID_list.append(line.split(" ")[0].strip("@"))
    r = []
    assert len(seqID_list) == len(seq_list)
    for i in range(len(seq_list)):
        r.append((seqID_list[i], seq_list[i]))
    return r


def old_levenshtein(s1, s2):
    """
    Compute a levenshtein distance between string s1 and string s2
    """
    l1 = len(s1)
    l2 = len(s2)

    matrix = [range(l1 + 1)] * (l2 + 1)
    for zz in range(l2 + 1):
      matrix[zz] = range(zz,zz + l1 + 1)
    for zz in range(0,l2):
      for sz in range(0,l1):
        if s1[sz] == s2[zz]:
          matrix[zz+1][sz+1] = min(matrix[zz+1][sz] + 1, matrix[zz][sz+1] + 1, matrix[zz][sz])
        else:
          matrix[zz+1][sz+1] = min(matrix[zz+1][sz] + 1, matrix[zz][sz+1] + 1, matrix[zz][sz] + 1)
    return matrix[l2][l1]

#Mimic library fusion
def one_fusion(params=None):
    """
    Mimic a single library fusion
    """
    if params == None:
        mix = [0, 0.05, 0.40, 0.50, 0.05]
    else:
        mix = params[:]
    for i in range(len(mix)-1):
        mix[i+1] += mix[i]
    x = random()
    p = -1
    for i in range(len(mix)):
        if x < mix[-i]:
            p = len(mix)-i
    return p

def many_fusions(wells, weights=None, nb_fusion=10000, params=None):
    """
    Mimic the whole library fusion
    """
    if weights == None:
        weights = np.array([1./len(wells)]*len(wells))
    
    fusion_dist = [one_fusion(params=params) for _ in xrange(nb_fusion)]
    return [azm.transform(list(chain.from_iterable([np.random.choice(wells, p=weights) for _ in range(p)]))) for p in fusion_dist]

#Here the mean is taken with the fractions computed with all 16 nodes
def mean_data_old(data, hp_sf, az_sf):
    temp = data[(data.az_total >= az_sf) & (data.hp_total >= hp_sf)]
    nt_list = []
    values = []
    az_prop = []
    for name, group in temp.groupby("az_list"):
        nt_list.append(name)
        values.append((norm_array(group.get(G).values)).mean(0))
        a = group.get("az_total").values.sum()
        b = group.get("hp_total").values.sum()
        az_prop.append(a/(a+b))
    temp = pd.DataFrame(np.vstack(values), columns=G)
    temp["az_list"] = nt_list
    temp["az_prop"] = az_prop
    temp["nb_az"] = temp.az_list.map(lambda x: x.count('1'))
    temp["az_gd"] = temp.apply(genotype_distribution, axis=1)
    return temp

def mean_data(data, hp_sf, az_sf):
    temp = data[(data.az_total >= az_sf) & (data.hp_total >= hp_sf)]
    nt_list = []
    values = []
    az_UMI = []
    conc = []
    for name, group in temp.groupby("az_list"):
        nt_list.append(name)
        nodes = azm.transform(name)
        y = (norm_array(group.get(nodes).values)).mean(0)
        values.append([y[nodes.index(g)] if g in nodes else 0. for g in G ])
        az_UMI.append(group.get("az_total").mean())
        c = np.mean(group["az_total"].values*group["nb_hp"].values*10.0/group["hp_total"].values)
        conc.append(c)
    temp = pd.DataFrame(np.vstack(values), columns=G)
    temp["az_list"] = nt_list
    temp["az_UMI"] = az_UMI
    temp["conc"] = conc
    temp["nb_az"] = temp.az_list.map(lambda x: x.count('1'))
    temp["az_gd"] = temp.apply(genotype_distribution, axis=1)
    return temp

def load_jan17():
    cwd = os.getcwd()
    os.chdir("D:/EvoEvo/miseq17/jan17/droplet/")
    data = pd.read_csv("data_jan17.csv", converters={'az_list':str}, index_col="dpBC")
    os.chdir(cwd)
    return data

def load_fitted_rates_jan17():
    df = pd.read_csv("D:/EvoEvo/miseq17/jan17/droplet/az_rates_fitted.csv", index_col=0)
    d = {}
    for col in df.columns:
        d[col] = (list(df[col].values))
    return d

def genotype_distribution(x):
    Y = list(x[azm.transform(x["az_list"])])
    if sum(Y) == 0:
        return [0.0]*(x["az_list"].count('1'))
    Y = [y/sum(Y) for y in Y]
    return Y

def list_hp(x):
    HP = ["hp_%s"%i for i in range(24)]
    t = x["hp_threshold"]
    return list(np.nonzero((x[HP]*(x[HP] >= t)[HP]).values)[0])

def list_az(x):
    E = [l.strip().split(",") for l in open("D:/EvoEvo/miseq17/jan17/droplet/emulsions_jan17.txt").readlines()]
    return azm.transform(list(set((",".join([",".join(E[i]) for i in x["hp_list"]])).split(","))))

def az_proportion_correct(x):
    if x["nb_az"] == 0 or x["az_total"] == 0:
        return 0
    return sum(x[azm.transform(x["az_list"])])/float(x["az_total"])

def hp_proportion_correct(x):
    if x["nb_hp"] == 0 or x["hp_total"] == 0:
        return 0
    ix = ["hp_%s"%i for i in x["hp_list"]]
    return sum(x[ix])/float(x["hp_total"])

def hp_threshold(x, t=0.1):
    return (x["hp_total"]*t)

def distance_to_predictions(x, norm=0):
    Y = list(x[azm.transform(x["az_list"])])
    if sum(Y) == 0:
        return 0
    Y = [y/sum(Y) for y in Y]
    Yth = list(azm.asymptotic_solution(x["az_list"])[0])
    if norm == 0:
        return sqrt(sum([(y-yth)**2 for y, yth in zip(Y,Yth)]))
    else:
        return sqrt(sum([(y-yth)**2 for y, yth in zip(Y,Yth)])/float(len(Yth)))
    
def compute_data(df, dpBC="dpBC", t=0.07):
    
    data = df.groupby([dpBC, "new_G"]).size().reset_index(name="nb_reads").pivot(index=dpBC, columns="new_G", values="nb_reads")

    data.fillna(0, inplace=True)
    if "hp_24" in data.columns:
        del data["hp_24"]

    data["hp_total"] = data.get(["hp_%s"%i for i in range(24)]).sum(1)
    data["az_total"] = data.get(G).sum(1)
    data["total"] = data["hp_total"] + data["az_total"]
    
    data["hp_threshold"] = data.apply(lambda x: hp_threshold(x, t=t), axis=1)

    data["hp_list"] = data.apply(list_hp, axis=1)
    data["az_list"] = data.apply(list_az, axis=1)

    data["nb_hp"] = data["hp_list"].apply(len)
    data["nb_az"] = data["az_list"].apply(lambda x: x.count('1'))

    data["az_correct"] = data.apply(az_proportion_correct, axis=1)
    data["hp_correct"] = data.apply(hp_proportion_correct, axis=1)

    #data["dth"] = data.apply(lambda x: distance_to_predictions(x, norm=0), axis=1)
    #data["dth_normed"] = data.apply(lambda x: distance_to_predictions(x, norm=1), axis=1)
    
    return data

def filter_and_mean(data, filters, only_correct=1, per_node=0):
    """
    Takes unfiltered data and apply filters (eg: {"A":[10, 20], "C":[10**10, 10**10], "E":[5, 5], "F":[10**10, 10**10]}) specific for each experiment.
    Either apply AZ UMI filter on "correct" UMI (only_correct = 1) or on total # AZ UMI.
    """
    def is_ok(x, only_correct):
        if per_node == 1:
            if only_correct == 1:
                return bool((x.hp_total >= filters[x.spBC_letter][0]) * (x.az_total*x.az_correct >= filters[x.spBC_letter][1]*x.nb_az))
            return bool((x.hp_total >= filters[x.spBC_letter][0]) * (x.az_total >= filters[x.spBC_letter][1]*x.nb_az))
        else:
            if only_correct == 1:
                return bool((x.hp_total >= filters[x.spBC_letter][0]) * (x.az_total*x.az_correct >= filters[x.spBC_letter][1]))
            return bool((x.hp_total >= filters[x.spBC_letter][0]) * (x.az_total >= filters[x.spBC_letter][1]))           

    M = data.apply(lambda x: is_ok(x, only_correct=only_correct), axis=1)
    data_temp = data[M]
    
    if len(data_temp > 0):
        #Remove networks with no UMI for all nodes (useless if only_correct=1)
        data_temp = data_temp[~data_temp.apply(genotype_distribution, axis=1).apply(lambda x: sum(x) == 0)]

        data_ = mean_data(data_temp, 0, 0)
        return (data_, data_temp)
    else:
        return "empty"

"""
Functions related to the "sample graph" = the graph in which vertices are networks and in which edges connect two networks that differ only by an addition or by a substitution
"""

def is_neighbor(a, b):
    diff = sum([x != y for x, y in zip(list(a), list(b))])
    if a.count("1") != b.count("1") and diff == 1:
        return "D"
    elif a.count("1") == b.count("1") and diff == 2:
        return "S"
    return "N"

def sample_graph_matrix_from_list(list_nt, sub=1.0, deletion=1.0):
    n = len(list_nt)
    S = np.zeros((n, n))
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            link = is_neighbor(list_nt[i], list_nt[j])
            if link == "S": 
                S[i, j] = sub
            elif link == "D":
                D[i, j] = deletion
    return (S, D)

def sample_graph_matrix_from_data(g, nt_size, data, sub=1.0, deletion=1.0):
    temp = data[data.nb_az == nt_size]
    list_nt = [x for x in list(temp.az_list.values) if x[G.index(g)] == "1"]
    return sample_graph_matrix_from_list(list_nt, sub=sub, deletion=deletion)

def inverse_presence(x, i):
    if x[i] == "1":
        y = "0"
    else:
        y = "1"
    return x[:i] + y + x[i+1:]

def neighborhood(x, g=None, sub=1, deletion=0):
    ix = [i for i in range(len(x)) if x[i] == '1']
    iy = [i for i in range(len(x)) if i not in ix]
    res = []
    for i in ix:
        for j in iy:
            res.append(inverse_presence(inverse_presence(x, i), j))
    if g != None:
        return [y for y in res if y[G.index(g)] == '1']
    return res

def sample_graph_for_series(epi, start=1, end=12):
    temp = epi[epi.what == 0]
    graph = {}
    for s in range(1, 12):
        list_nt = list(temp[temp.nt_size_before == s].before.unique())
        for nt in list_nt:
            graph[nt] = list(temp[temp.before == nt].after.values)
    gr = nx.DiGraph(graph)
    return gr

def epi_dict(epi, rates):
    epi_dict = {}
    temp = epi[epi.what == 0]
    def compute_aux(x):
        return azm.degree(x["mutation"][1], x["after"], deg="in", norm=0, rates=rates) + azm.degree(x["mutation"][1], x["after"], deg="out", norm=0, rates=rates)
    temp["total_degree"] = temp.apply(compute_aux, axis=1)
    for name, x in temp.groupby("before"):
        epi_dict[name] = {y.after:[y.pscore, y.pscore_th, y.total_degree] for _, y in x.iterrows()}
    return epi_dict
    

def compute_series(start, end, epi_dict, gr):
    sources = [node for node in gr.nodes() if node.count("1") == start]
    targets = [node for node in gr.nodes() if node.count("1") == end]

    tot = []
    for s in sources:
        for t in targets:
            for path in nx.all_simple_paths(gr, s, t):
                pscores, pscores_th, degree_path = [0], [0], [0]
                for i in range(len(path) - 1):
                    pscores.append(pscores[-1] + epi_dict[path[i]][path[i+1]][0])
                    pscores_th.append(pscores_th[-1] + epi_dict[path[i]][path[i+1]][1])
                    degree_path.append(epi_dict[path[i]][path[i+1]][-1])
                tot.append([path[0], path[-1], pscores[-1], pscores_th[-1], path, pscores, pscores_th, sum(degree_path), degree_path])
    paths = pd.DataFrame(tot, columns=["start", "end", "cumu_pscore", "cumu_pscore_th", "path", "pscores", "pscores_th", "cumu_degree", "degree_path"])
    
    return paths



def print_full(x):
    pd.set_option('display.max_rows', x.shape[0])
    pd.set_option('display.max_columns', x.shape[1])
    pd.set_option('display.width', None)
    print(x)
    pd.reset_option('display.max_rows')
    pd.reset_option('display.max_columns')
    pd.reset_option('display.width')

# Might not be useful as izip() does the same
def two_gen(handle_one, handle_two):
    entry_one = True
    entry_two = True
    while entry_one and entry_two:
        try:
            entry_one = handle_one.next()
        except StopIteration:
            entry_one = None
        try:
            entry_two = handle_two.next()
        except StopIteration:
            entry_two = None
        if entry_one == None or entry_two == None:
            break
        yield (entry_one, entry_two)

def three_gen(handle_one, handle_two, handle_three):
    entry_one = True
    entry_two = True
    entry_three = True
    while entry_one and entry_two and entry_three:
        try:
            entry_one = handle_one.next()
        except StopIteration:
            entry_one = None
        try:
            entry_two = handle_two.next()
        except StopIteration:
            entry_two = None
        try:
            entry_three = handle_three.next()
        except StopIteration:
            entry_three = None
        if entry_one == None or entry_two == None or entry_three == None:
            break
        yield (entry_one, entry_two, entry_three)

def merge_reads(rd_one, rd_two):
    counter = 0
    merged = open(rd_one[:-6] + rd_two[:-6] + ".fastq","w")
    for (g, d) in two_gen(SeqIO.parse(rd_one, "fastq"), SeqIO.parse(rd_two, "fastq")):
        counter += 1
        SeqIO.write(g+d, merged, "fastq")
        if counter % 10000 == 0:
            print("%s sequences treated"%counter)
    print(counter)
    merged.close()