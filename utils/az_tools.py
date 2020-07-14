import pandas as pd
import numpy as np
import networkx as nx
import az_model as azm

G = [M+N for M in "ACUG" for N in "ACUG"]
base_comp = {"A":"U", "U":"A", "C":"G", "G":"C"}

#Transform list of genotypes into network string and vice-versa
def transform(x):
    list_genotypes = [M+N for M in "ACUG" for N in "ACUG"]
    if x == None:
        return x
    elif type(x) == list:
        return "".join([str(int(g in x)) for g in list_genotypes])
    else:
        return [list_genotypes[i] for i in range(16) if x[i] == '1']

def norm_array_aux(x):
    #Rows of x = samples
    #Columns of x = observations
    return x.sum(1).repeat(x.shape[1]).reshape(x.shape)

def norm_array(x):
    if type(x) == list: 
        x = np.array(x)
    if x.ndim > 1:
        y = x/norm_array_aux(x)
        return np.nan_to_num(y)
    else:
        if x.sum() > 0:
            return x/float(x.sum())
        else:
            return np.array([1./len(x)]*len(x))
        
def genotype_distribution(x):
    # Returns as a list the node fractions for a given network
    Y = list(x[transform(x["az_list"])])
    if sum(Y) == 0:
        return [0.0]*(x["az_list"].count('1'))
    Y = [y/sum(Y) for y in Y]
    return Y

def filter_az_droplet_data(data, filters):
    """
    Process unfiltered azoarcus network droplets data by:
        1) Applying filters on number of azoarcus and hairpin UMIs (only UMIs which are associated with expected azoarcus genotypes are considered).
        2) Compute the mean of node fractions, total number of azoarcus UMIs and yield per azoarcus network when several droplets contain the same network
        
    Return a tuple containing the data before and after the mean node fractions per droplet have been computed.
    """
    def is_ok(x):
        return bool((x.hp_total >= filters[0]) * (x.az_total*x.az_correct >= filters[1]))

    M = data.apply(lambda x: is_ok(x), axis=1)
    data_temp = data[M]
    assert len(data_temp) > 0
    data_temp["yield"] = data_temp.apply(compute_yield, axis=1)
    data_temp["az_gd"] = data_temp.apply(genotype_distribution, axis=1)
    data_ = compute_mean_per_network(data_temp)
    
    return (data_, data_temp)

def compute_mean_per_network(data):
    """
    Compute the mean of node fractions, total number of azoarcus UMIs and yield per azoarcus network when several droplets contain the same network.
    """
    temp = data[:]
    nt_list = []
    values = []
    az_UMI = []
    conc = []
    for name, group in temp.groupby("az_list"):
        nt_list.append(name)
        nodes = transform(name)
        # Node fractions
        y = (norm_array(group.get(nodes).values)).mean(0)
        values.append([y[nodes.index(g)] if g in nodes else 0. for g in G])
        # Total number of azoarcus UMIs
        az_UMI.append(group.get("az_total").mean())
        # Yield
        conc.append(np.mean(group["yield"].values))
    temp = pd.DataFrame(np.vstack(values), columns=G)
    temp["az_list"] = nt_list
    temp["az_UMI"] = az_UMI
    temp["yield"] = conc
    temp["nb_az"] = temp.az_list.map(lambda x: x.count('1'))
    temp["az_gd"] = temp.apply(genotype_distribution, axis=1)
    return temp

def compute_yield(x):
    """
    Compute yield as measured with UMIs for a droplet x. 
    """
    return x["az_total"]*x["nb_hp"]*10.0/x["hp_total"]

# Perturbation analysis

def perturbation_analysis(node_fraction_data):
    """
    This function computes the perturbation analysis using node fraction data which can come from the measured node fractions in droplets or from a model (kinetic, eigenvector centrality or in-degree centrality).
    Results of the perturbation analysis is stored in "epi".
    """
    epi = []
    
    # Construction of the graph of azoarcus networks
    list_nt = node_fraction_data.keys()
    S = nx.from_numpy_matrix(az_networks_graph(list_nt))

    # We traverse the graph of azoarcus networks and compute perturbation metrics between each connected pair of azoarcus networks.
    for e in S.edges:
        for i, j in [(0, 1), (1, 0)]:
            b = list_nt[e[i]]
            a = list_nt[e[j]]
            mut = what_mut(a, b)
            if mut[0] == "-":
                mut = mut[1]
                order, y_before, y_after, pscore = compute_pscore(b, a, node_fraction_data[b], node_fraction_data[a])
                epi.append([b, a, mut, pscore, b.count("1"), y_before, y_after, order])
            
    epi = pd.DataFrame(epi, columns=["before", "after", "mutation", "pscore", "nt_size_before", "y_before", "y_after", "order_sp"])
    return epi

def az_networks_graph(list_nt):
    """
    Constructs the graph of azoarcus networks from a list of azoarcus networks. 
    In this graph, the vertices are azoarcus networks and edges connect two networks that differ only by the addition/deletion of a node.
    """
    n = len(list_nt)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            D[i, j] = is_neighbor(list_nt[i], list_nt[j])
    return D

def is_neighbor(a, b):
    """
    Returns whether or not two networks are neighbors. 
    Two networks are neighbors if and only if they differ by the deletion/addition of a node.
    """
    diff = sum([x != y for x, y in zip(list(a), list(b))])
    if a.count("1") != b.count("1") and diff == 1:
        return True
    return False

def compute_diff(A, B):
    #Auxiliary function to what_mut function
    return [-1*(a=='1' and b=='0') + 1*(a=='0' and b=='1') for a, b in zip(list(B), list(A))]

def what_mut(A, B):
    #Returns which "mutation" is required to go from network B to network A
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

def compute_pscore(before, after, nf_before, nf_after):
    # Computes perturbation metrics between networks before and after (which are neighbors hence they differ only by the addition/deletion of a node)
    y_before = {g:nf_before[transform(before).index(g)] for g in transform(before)}
    y_after = {g:nf_after[transform(after).index(g)] for g in transform(after)}
    common = list(set(y_before.keys()).intersection(set(y_after.keys())))
    common = transform(transform(common))
    y = norm_array(np.array([y_after[c] for c in common]))
    x = norm_array(np.array([y_before[c] for c in common]))
    return [common, x, y, np.mean(np.abs(y-x))]

def postprocess_perturbation_analysis(size, epi, rates):
    """
    Postprocess results from perturbation analysis:
        1. Selects only additions of a node
        2. Selects a given starting network size (4 nodes for example)
        3. Computes relevant quantities from network theory using interaction weights data given in the last argument "rates"
        4. Give the result in a per-node basis (function collapse_postprocess aggregates the result per pair of networks)
    """
    temp = epi[epi.nt_size_before == size]
    df = []
    for name, x in temp[:].iterrows():
        v = transform(x.before)
        n = x.before.count('1')
        res = []
        d_a_out = sum([rates[G.index(x.mutation[0] + g[1])] for g in v])
        sigma_G = azm.network_matrix(x.before, rates=rates).sum()
        d_v_in = [azm.degree(g, x.before, norm=0, deg='in', rates=rates) for g in v]
        va = [g for g in transform(x.before) if np.isclose(rates[G.index(x.mutation[0] + g[1])], 0) == False]
        va_star = [g for g in transform(x.before) if g not in va]
        y = x.y_before
        y_ = x.y_after
        if sum(y_) == 0:
            y_ = np.array([1./n]*n)
        else:
            y_ = norm_array(y_)
        sv = y_ - y
        for i in range(n):
            res.append([x.before.count("1"), x.before, x.after, x.mutation, v[i], d_v_in[i], d_a_out, sigma_G, rates[G.index(x.mutation[0] + v[i][1])], sv[i], 
                        (rates[G.index(x.mutation[0] + v[i][1])] - d_a_out*d_v_in[i]/sigma_G)/(sigma_G + d_a_out),
                        sum([azm.degree(g, x.before, norm=0, deg='in', rates=rates) for g in va]), 
                        sum([azm.degree(g, x.before, norm=0, deg='in', rates=rates) for g in va_star])])
        df = df + res
    df = pd.DataFrame(df, columns=["nt_size", "before", "after", "mutation", "node", "deg_in_v", "deg_out_a", "sigma_g", "edge_a_v", "p_v", "p_v_th", "deg_in_va", "deg_in_va_star"])
    return df

def collapse_postprocess(temp, rates):
    """
    Aggregates per pair of networks (before, after) the results from postprocessing of perturbation analysis.
    """
    temp["abs_p_v"] = np.abs(temp.p_v.values)
    temp = temp.groupby(["before", "after"]).agg({'abs_p_v':'sum', 
                                                  'sigma_g':'mean', 
                                                  'edge_a_v':lambda x: max(x), 
                                                  'mutation':lambda x: list(x)[0]}).reset_index()
    temp["edge_a_v"] = temp.mutation.apply(lambda x: rates[G.index(x[0]+base_comp[x[0]])])
    temp["sigma_e"] = temp.sigma_g.values/temp.edge_a_v.values
    temp["is_cat_innov"] = temp.apply(lambda x: is_catalytic_innovation(x.before, x.mutation), axis=1)
    temp["m"] = temp.apply(lambda x: sum([g[0] == x.mutation[0] for g in transform(x.before)]), axis=1)
    temp["n"] = temp.apply(lambda x: sum([g[1] == base_comp[x.mutation[0]] for g in transform(x.before)]), axis=1)
    return temp

def is_catalytic_innovation(nt, added_node):
    """
    The added node is a "catalytic innovation" for the network "nt" if no other node has the same IGS and if it has at least a target different from itself.
    """
    is_new_igs = sum([g[0] == added_node[0] for g in transform(nt)]) == 0
    has_target = sum([g[1] == base_comp[added_node[0]] for g in transform(nt)]) > 0
    return is_new_igs and has_target

def at_least_one_WC_link(nt):
    """
    Tests if the network nt has at least a Watson-Crick interaction.
    """
    WC_rates = [1 if g in ["AU", "UA", "GC", "CG"] else 0 for g in G]
    M = azm.network_matrix(nt, WC_rates)
    return M.sum() > 0

def dir_graph_cumulative_perturbation_analysis(epi, start=1, end=12):
    """
    Computes the directed graph of azoarcus networks with at least `start` nodes and at most `end` nodes for cumulative perturbation analysis where an interaction represents the addition of a node to a network.
    """
    graph = {}
    for s in range(start, end):
        list_nt = list(epi[epi.nt_size_before == s].before.unique())
        for nt in list_nt:
            graph[nt] = list(epi[epi.before == nt].after.values)
    gr = nx.DiGraph(graph)
    return gr

def cumulative_perturbation_analysis(start, end, epi_dict, gr):
    """
    Computes trajectories of cumulative perturbations starting from networks with number of nodes == `start` to networks with number of nodes == `end`.
    Takes two auxiliary data structures:
        - `gr`, the directed graph of azoarcus networks computed with `dir_graph_cumulative_perturbation_analysis`
        - `epi_dict`, a dictionnary storing perturbation score between a source network and all its targets for fast look-up
    """
    sources = [node for node in gr.nodes() if node.count("1") == start]
    targets = [node for node in gr.nodes() if node.count("1") == end]

    tot = []
    for s in sources:
        for t in targets:
            for path in nx.all_simple_paths(gr, s, t):
                pscores = [0]
                for i in range(len(path) - 1):
                    pscores.append(pscores[-1] + epi_dict[path[i]][path[i+1]])
                tot.append([path[0], path[-1], pscores[-1], path, pscores])
    paths = pd.DataFrame(tot, columns=["start", "end", "cumu_pscore", "path", "pscores"])
    
    return paths