# -*- coding: utf-8 -*-
"""
Created on Wed Feb 01 11:20:57 2017

@author: Simon
"""
import az_tools_ as az
import az_model_ as azm
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy import stats
#import graphviz as gv

#Color and other useful constants
golden = 1.618
colors=[u'#1b9e77', u'#d95f02', u'#7570b3', u'#e7298a', u'#66a61e', u'#e6ab02', u'#a6761d', u'#666666']
almost_black = "0.2"
# E = [l.strip().split(",") for l in open("initial_emulsions.txt").readlines()]
G = [M+N for M in "ACUG" for N in "ACUG"]
WC_rates = {'CG':1, 'GC':1, 'AU':1, 'UA':1} 
WC_rates = [WC_rates[g] if g in WC_rates.keys() else 0. for g in G]


#Draw a network (only Watson-Crick interactions) with network string or list of genotypes
def draw_nt(x, ax=None):
    if type(x) == str:
        g = azm.graph(x, WC=True)
        nodes = dict(enumerate(azm.transform(x)))
        g = nx.relabel_nodes(g, nodes)
        nx.draw(g, ax=ax, arrows=True, with_labels=True, node_color="w", node_size=600, edge_color="0.2", font_color="0.2")
    if type(x) == list:
        draw_nt(azm.transform(x))
   
#Barplot the asymptotic solution of a network     
def az_barplot(nt, ax=None, color=colors[0], rotation=90, rates=0):
    if type(nt) == str: nt = azm.transform(nt)
    else: nt = azm.transform(azm.transform(nt))
    ret = 0
    if ax == None:
        ret = 1
        f = plt.figure(figsize=(0.5*len(nt), 1.5))
        ax = plt.gca()
    ax.bar(range(len(nt)), azm.asymptotic_solution(nt, rates=rates)[0], color=color, lw=0, align='center')
    ha(ax)
    ax.set_xlim(-1, len(nt))
    ax.set_xticks(range(len(nt)))
    ax.set_xticklabels(nt, rotation=90)
    if ret == 1:
        return f

#Hide top and right axis and ticks
def ha(ax):
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    
def CI_linreg(X, Y, lr, alpha=0.05, x=None):
    """
    Compute the confidence and prediction intervals
    X, Y = data
    lr = linear regression (using scipy.stats.linregress)
    alpha = threshold
    x = range over which to compute the CI
    """
    if x == None:
        x = np.linspace(X.min(), X.max(), 100)
    n = len(X)
    Y_hat = lr.slope*X + lr.intercept
    y = lr.slope*x + lr.intercept
    s_err = np.sum((Y-Y_hat)**2)
    sd = np.sqrt(1.0/(n-2)*s_err)
    sxd = np.sum((X-X.mean())**2)
    sx = (x-X.mean())**2
    ttest = stats.t.ppf(1-alpha/2.0,n-2)
    dy = ttest*sd*np.sqrt(1 + 1.0/n + sx/sxd)
    dslope = ttest*sd*np.sqrt(1.0/n + sx/sxd)
    
    return (x, y - dy, y + dy, y - abs(dslope), y + abs(dslope))

def draw_CI_linreg(ax, X, Y, lr, color='r', alpha=0.05, x=None, CI_95=0, fade=1.):
    #Draw onto axes the confidence and prediction intervals
    CI = CI_linreg(X, Y, lr, alpha, x)
    ax.plot(CI[0], CI[0]*lr.slope+lr.intercept, '--', color=color, alpha=fade)
    if CI_95: ax.fill_between(CI[0], CI[1], CI[2], alpha=0.1, lw=0, facecolor=color)
    ax.fill_between(CI[0], CI[-2], CI[-1], alpha=fade/3, lw=0, facecolor=color)
    
    
#Usual plots collection

def fusion_dist(data, sf, var, color):
    temp = data[(data[var] >= sf)].get("nb_hp").values
    X, Y = np.unique(temp, return_counts=True)
    f = plt.figure()
    ax = plt.gca()
    ax.bar(X, Y, align='center', lw=0, color=color);
    ax.set_xticks(range(max(X)+1))
    ax.set_xlim(-1, 7)
    ax.set_xlabel("Number of hairpin-reporters per droplet-barcode", fontsize=12)
    ax.set_ylabel("Number of droplets", fontsize=12)
    return f
    
def correct_genotypes(data, sf, var, color):
    f = plt.figure()
    ax = plt.gca()
    ax.hist(data[(data[var] >= sf)].get("az_correct").values, bins=20, color=color)
    ax.set_xlim(-0.05, 1.05)
    ax.set_xlabel("Proportion of correct genotypes", fontsize=12)
    ax.set_ylabel("Number of droplets", fontsize=12)
    return f
    
def correct_vs_filter_value(data, var, color):
    Y = [data[(data[var] >= j)].get("az_correct").values.mean() for j in range(500)]
    X = range(500)
    f = plt.figure(figsize=(6, 6/golden))
    ax = plt.gca()
    ax.plot(X, Y, '.', ms=5, color=color)
    ax.set_xlabel("Size filter value")
    ax.set_ylabel("Proportion of correct genotypes")
    ax.set_xticks(np.arange(10)*50)
    ax.grid()
    return f
    
def nt_size_dist(data, sf, v, color, sim=0):
    temp = data[data[v] >= sf]
    
    #Unique
    X = [x.count('1') for x in temp.az_list.unique()]
    #All
    #X = [x.count('1') for x in temp.az_list.values]
    
    if sim == 1: f = plt.figure(figsize=(8, 4/golden))
    else: f = plt.figure(figsize=(4, 4/golden))
    
    if sim == 1: ax = f.add_subplot(121)
    else: ax = plt.gca()
    ax.hist(X, np.arange(18)-0.5, color=color)
    ax.set_xlim(-1,17)
    ax.set_xlabel("Network size")
    ax.set_ylabel("Counts of unique network")
    
    if sim == 1:
        ax = f.add_subplot(122)
        ax.hist([x.count("1") for x in az.many_fusions(E, nb_fusion=len(temp))], np.arange(18)-0.5, color=color, alpha=0.5)
        ax.set_xlim(-1,17)
        ax.set_xlabel("Network size")
        ax.set_ylabel("Counts of unique network")
    
    ax.set_title("%s, %s"%(len(temp), len(X)))
    return f

def nb_network_vs_filter_value(data, v, color, start=20):

    Y_u = [len(data[(data[v] >= i)].az_list.unique()) for i in range(start, 500)]
    Y = [len(data[(data[v] >= i)].az_list.values) for i in range(start, 500)]
    
    X = range(start, 500)

    f = plt.figure(figsize=(16, 8))
    
    ax = f.add_subplot(121)
    ax.plot(X, Y_u, '.', ms=5, color=color)
    ax.set_xlabel("Size filter value")
    ax.set_ylabel("Number of unique networks")
    ax.locator_params(nbins=12)
    ax.grid()
    
    ax = f.add_subplot(122)
    ax.plot(X, Y, '.', ms=5, color=color, alpha=0.7)
    ax.set_xlabel("Size filter value")
    ax.set_ylabel("Number of networks")
    ax.locator_params(nbins=12)
    ax.grid()
    return f

def error_bar_dist(data, sf, v, color):
    temp = data[(data[v] >= sf)]
    liste_nt = temp.groupby("az_list").size().reset_index(name="nb_replicates")
    liste_nt = list(liste_nt[liste_nt.nb_replicates > 1].az_list.values)
    H = []
    for i, nt in enumerate(liste_nt):
        X = temp[temp.az_list == nt]
        Y = np.vstack(X.apply(az.genotype_distribution, axis=1))
        H += list(Y.std(0))
    f = plt.figure()
    ax = plt.gca()
    ax.hist(H, np.linspace(0, 0.5), color=color)
    ax.set_title("Mean standard deviation = %.2f"%(np.mean(H)))
    ax.set_xlabel("Standard deviation")
    ax.set_ylabel("Counts")
    return f

def fitness_dist(data, sf, v, color, nt_size):
    temp = data[(data[v] >= sf) & (data.nb_az.isin(nt_size))]
    Y = np.hstack(temp.apply(az.genotype_distribution, axis=1).values)
    f = plt.figure()
    ax = plt.gca()
    ax.hist(Y, np.linspace(0, 1), color=color)
    ax.set_xlim(0,1)
    ax.set_xlabel("Fitness")
    ax.set_ylabel("Counts")
    return f

def correct_hp(data, sf, var, color):
    f = plt.figure()
    ax = plt.gca()
    ax.hist(data[(data[var] >= sf)].get("hp_correct").values, bins=20, color=color)
    ax.set_xlim(-0.05, 1.05)
    ax.set_xlabel("Proportion of correct hairpins", fontsize=12)
    ax.set_ylabel("Number of droplets", fontsize=12)
    return f

def compute_prop_exch(x):
    l = azm.transform(x["az_list"])
    r = []
    for M in [g[0] for g in l]:
        for N in [g[1] for g in l]:
            r.append(M+N)
    e = set(r).difference(set(l))
    ne_e = set(G).difference(set(l))
    if sum(x[list(ne_e)]) == 0:
        return 1.0
    return sum(x[list(e)])/float(sum(x[list(ne_e)]))

def prop_exch(data, sf, var, color):
    temp = data[(data[var] >= sf)]
    Y = temp.apply(compute_prop_exch, axis=1)
    f = plt.figure()
    ax = plt.gca()
    ax.hist(Y, bins=20, color=color)
    ax.set_xlim(-0.05, 1.05)
    ax.set_xlabel("Proportion of exchanged genotype", fontsize=12)
    ax.set_ylabel("Number of droplets", fontsize=12)
    return f

def umi_dist(data, sf, v, color, sim=0):
    temp = data[:]
    f = plt.figure(figsize=(4, 4/golden))
    ax = plt.gca()
    ax.hist(temp.total, 50, log=1, color=color)
    ax.set_xlabel("#UMI/dp")
    ax.set_ylabel("Counts")
    return f

def plot_first_order(data, sf=50, var="total", color=colors[1]):
    correct_vs_filter_value(data, var, color)
    nb_network_vs_filter_value(data, var, color)
    fusion_dist(data, sf, var, color)
    nt_size_dist(data, sf, var, color)
    correct_genotypes(data, sf, var, color)
    correct_hp(data, sf, var, color)
    error_bar_dist(data, sf, var, color)
    fitness_dist(data, sf, var, color, [6])
    prop_exch(data, sf, var, color)
    umi_dist(data, sf, var, color)
    

def dp_in_details(dpBC, data, t=0.06):
    
    def aux(x):
        order = "ACUG"
        a = [[], [], [], []]
        for g in x:
            a[order.index(g[0])].append(g)
        return a
    
    
    temp = data[data.index == dpBC]
    
    Y_hp = temp.get(["hp_%s"%i for i in range(24)]).values[0]
    
    f = plt.figure(figsize=(8, 2))
    ax = plt.gca()
    ax.bar(range(24), Y_hp, align='center', lw=0, color=colors[0])
    ax.plot(np.linspace(0, 24), [temp.hp_total.values[0]*t]*50, "--", alpha=0.5, color=colors[7])
    ax.set_xticks(range(24))
    ax.set_xlim(-1, 24)
    ax.set_ylabel("Counts")
    ax.set_xlabel("Hairpin")
    ax.set_title("%s,%s"%(temp.index.values[0], temp.nb_hp.values[0]))
    
    #hp = temp["hp_list"].values[0]
    az_exp = azm.transform(temp["az_list"].values[0])
    
    Y_az = temp.get(G).values[0]

    f = plt.figure(figsize=(6, 2))
    ax = plt.gca()
    for i in range(16):
        c = colors[3]
        if G[i] in az_exp:
            c = colors[1]
        ax.bar(i, Y_az[i], align='center', lw=0, color=c)
    
    ylim = ax.get_ylim()[-1]
    ax.text(17, ylim*0.90, "%.2f"%(temp.az_correct.values[0]))
    ax.text(17, ylim*0.80, "\n".join([", ".join(x) for x in aux(az_exp)]), va='top')
    ax.set_title(temp.index.values[0])
    ax.set_xticks(range(16))
    ax.set_xticklabels(G, rotation=90)
    ax.set_xlim(-1,16)
    ax.set_ylabel("Counts")
    ax.set_xlabel("Azoarcus")
    
def boxplot(ax, x, y, nbins=15, showmeans=True, log=False, sym='.'):
    if log:
        binning = np.logspace(np.log10(x).min(), np.log10(x).max(), nbins)
        midpoints = np.log10(0.5 * (binning[1:] + binning[:-1]))
        widths = [0.8 * (midpoints[1] - midpoints[0])]*len(midpoints)
        idx = np.digitize(x, binning)
        temp = [y[idx == i] for i in range(len(binning))]
        m = AutoLocator()
        ax.boxplot(temp[1:], positions = midpoints, widths=widths, sym=sym, showmeans=showmeans, meanline=True, flierprops=dict(marker='.', markerfacecolor="0.2", markersize=5, markeredgewidth=0.0))
        ax.set_xlim(np.log10(x).min(), np.log10(x).max())
        #ax.set_xticks(m.tick_values(np.log10(x).min(), np.log10(x).max()))
        ax.set_xticks(np.log10(m.tick_values(x.min(), x.max())[1:]))
        ax.set_xticklabels(m.tick_values(x.min(), x.max())[1:], rotation=90)
    else:
        h = (x.max() - x.min())*0.05
        if type(nbins) == int:
            binning = np.linspace(x.min()-h, x.max()+h, nbins)
        else:
            binning = np.array(nbins)
        midpoints = 0.5 * (binning[1:] + binning[:-1])
        widths = [0.8 * (midpoints[1] - midpoints[0])]*len(midpoints) 
        xmin, xmax = x.min()-h, x.max()+h
        idx = np.digitize(x, binning)
        temp = [y[idx == i] for i in range(len(binning))]
        ax.boxplot(temp[1:], positions = midpoints, widths=widths, sym=sym, showmeans=showmeans, meanline=True, flierprops=dict(marker='.', markerfacecolor="0.2", markersize=5, markeredgewidth=0.0))
        m = AutoLocator()
        ax.set_xlim(xmin, xmax)
        ax.set_xticks(m.tick_values(xmin, xmax))
        ax.set_xticklabels(m.tick_values(xmin, xmax))
        
def styled_boxplot(ax, x, y, nbins=15, plot=True, showmeans=True, log=False, sym='.', nb_points=1, showfliers=True, color=colors[1], alpha=0.3, ms=10, showcaps=True, alphafliers=0.1, whis=1.5, xlim=None, ticks=True):
    if type(nbins) == int:
        binning = np.linspace(x.min(), x.max(), nbins)
    else:
        binning = np.array(nbins)
        nbins = len(binning)
        
    midpoints = 0.5 * (binning[1:] + binning[:-1])
    widths = [0.8 * (midpoints[1] - midpoints[0])]*len(midpoints) 
    idx = np.digitize(x, binning)
    temp = [y[idx == i] for i in range(len(binning))]
    
    #Filter box with less than nb_points
    jdx = [i for i in range(len(temp)) if len(temp[i]) >= nb_points]
    temp = [temp[i] for i in range(nbins) if i in jdx]
    midpoints = [midpoints[i] for i in range(nbins-1) if i+1 in jdx]
    widths = [widths[i] for i in range(nbins-1) if i+1 in jdx]
    
    if plot:
        ax.boxplot(temp, positions = midpoints, widths=widths, sym=sym, showmeans=showmeans, meanline=False, showfliers=showfliers, showbox=True, showcaps=showcaps, whis=whis,
                   flierprops=dict(marker='.', markerfacecolor="0.2", markersize=5, markeredgewidth=0.0, alpha=alphafliers),
                   boxprops=dict(alpha=alpha), capprops=dict(alpha=alpha), whiskerprops=dict(alpha=alpha), 
                   medianprops=dict(alpha=0.), meanprops=dict(marker=".", markersize=ms, markerfacecolor=color, markeredgecolor=color))
        if xlim == None:
            xmin, xmax = midpoints[0] - widths[0], midpoints[-1] + widths[-1]
        else:
            xmin, xmax = xlim[0], xlim[1]
        if ticks:
            m = mpl.ticker.AutoLocator()
            ax.set_xticks(m.tick_values(xmin, xmax))
            ax.set_xticklabels(m.tick_values(xmin, xmax))
            ax.set_xlim(xmin, xmax)
    
    return (binning, midpoints, widths, jdx, temp)
"""        
def az_graphviz(nt, engine="neato", name="_", fmt="svg", dpi="70", colors=None, edge_colors=None, penwidth='1'):
    nt = azm.transform(azm.transform(nt))
    if type(nt) == str: nt = azm.transform(nt)
    g = azm.graph(nt, WC_rates)
    g = nx.relabel_nodes(g, {i:nt[i] for i in range(len(nt))})
    if name == "_": name = '_'.join(nt)
    gdot = gv.Digraph(engine=engine, name=name, format=fmt, edge_attr={"arrowhead":"open", "arrowsize":"0.7"}, 
                      node_attr={"fontname":"Calibri", "fontsize":"12", "width":"0.3", "height":"0.3", "fixedsize":"true", "style":"filled", 'penwidth':penwidth},
                      graph_attr={"nodesep":"0.2", "dpi":dpi})
    if colors == None:
        colors = {n:"white" for n in g.nodes()}
    if edge_colors == None:
        edge_colors = {n:"black" for n in g.nodes()}
    for nd in g.nodes():
        gdot.node(nd, color=edge_colors[nd], fillcolor=colors[nd])
    gdot.edges(g.edges())
    gdot.render()
    return gdot
    
def general_graphviz(g, directed=True, name="graph", engine="neato", fmt="svg", label=True, 
                     diameter="0.3", nodesep="0.2", arrowsize="0.7", size="2", dpi="70", ranksep="0.75",
                     shape="point", penwidth="0.2"):
    if directed == True:
        gdot = gv.Digraph(engine=engine, name=name, format=fmt, edge_attr={"arrowhead":"open", "arrowsize":arrowsize, "penwidth":penwidth}, 
                          node_attr={"fontname":"Calibri", "fontsize":"12", "width":diameter, "height":diameter, "fixedsize":"true", "shape":shape, "penwidth":penwidth},
                          graph_attr={"nodesep":nodesep, "center":"true", "size":size, "dpi":dpi, "ranksep":ranksep, "pad":"0.1"})
    else:
        gdot = gv.Graph(engine=engine, name=name, format=fmt, 
                          node_attr={"fontname":"Calibri", "fontsize":"12", "width":diameter, "height":diameter, "fixedsize":"true"},
                          graph_attr={"nodesep":nodesep, "center":"true", "size":size, "dpi":dpi})
    for nd in g.nodes():
        if label: gdot.node(str(nd))
        else: gdot.node(str(nd), label="")
    gdot.edges([map(str, e) for e in g.edges()])
    gdot.render()
"""    
def repr_lr(lr):
    return "slope = %.2g, p = %.2g, r = %.2g"%(lr.slope, lr.pvalue, lr.rvalue)

def custom_log_ticks(a, b, base_10=False):
    ticks = []
    for i in range(a, b):
        ticks += list(np.linspace(10**i, 10**(i+1), 10)[:-1])
    ticks += [10**b]
    aux = ["%.1g"%10**i for i in range(a, b+1)]
    if base_10:
        labels = [r"$10^{%.0f}$"%np.log10(t) if "%.1g"%t in aux else "" for t in ticks]
    else:
        labels = ["%.1g"%t if "%.1g"%t in aux else "" for t in ticks]
    return (np.log10(ticks)[::9], labels[::9], np.log10(ticks))

def what_pvalue(p):
    if p <= 0.001:
        return "***"
    if p <= 0.01:
        return "**"
    if p <= 0.05:
        return "*"
    return "NS"