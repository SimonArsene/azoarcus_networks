import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

#Color and other useful constants
golden = 1.618
colors = [u'#1b9e77', u'#d95f02', u'#7570b3', u'#e7298a', u'#66a61e', u'#e6ab02', u'#a6761d', u'#666666']
almost_black = "0.2"

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