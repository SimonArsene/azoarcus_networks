import az_tools as az
import numpy as np
from scipy.integrate import odeint

G = [M+N for M in "ACUG" for N in "ACUG"]

def network_matrix(nt, rates):
    """
    Returns the network matrix of a network nt
    Takes a list of interaction weigths (in the same order as the list of nodes) as the "rates" argument
    """
    if type(nt) == str: nt = az.transform(nt)
    else: nt = az.transform(az.transform(nt))
    n = len(nt)
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            M[i,j] = rates[G.index(nt[j][0] + nt[i][-1])]
    return M

def first_order_model(nt, rates):
    """
    Returns the first-order model asymptotic solution for a network nt.
    Takes a list of interaction weigths (in the same order as the list of nodes) as the "rates" argument
    """
    if type(nt) == list:
        nt = az.transform(nt)
        M = network_matrix(nt, rates=rates)
    elif type(nt) == np.ndarray:
        M = nt
        nt = None
    else:
        M = network_matrix(nt, rates=rates)
    L, V = np.linalg.eig(M)
    kmax = np.real(L).argmax()
    return (np.real(V[:,kmax])/np.real(V[:,kmax]).sum(), az.transform(nt))

def degree(g, x, norm, deg, rates):
    """
    Computes the in or out (controlled with the "deg" argument) degree of a node "g" in a network "x".
    Takes a list of interaction weigths (in the same order as the list of nodes) as the "rates" argument.
    The "norm" argument is used to turn on or off the normalisation.
    If the network has no link at all, returns 1./(number of nodes).
    """
    if type(x) == list: 
        x = az.transform(x)
    if norm == 1:
        if network_matrix(x, rates).sum() == 0:
            return 1./x.count("1")
        if deg == "in":
            return network_matrix(x, rates).sum(1)[az.transform(x).index(g)]/network_matrix(x, rates).sum()
        else:
            return network_matrix(x, rates).sum(0)[az.transform(x).index(g)]/network_matrix(x, rates).sum()
    else:
        if deg == "in":
            return network_matrix(x, rates).sum(1)[az.transform(x).index(g)]
        else:
            return network_matrix(x, rates).sum(0)[az.transform(x).index(g)]

def in_degree_centrality(nt, rates):
    if type(nt) == str: 
        nt = az.transform(nt)    
    return [degree(g, nt, 1, 'in', rates) for g in nt]
    
def network_matrix(nt, rates):
    """
    Returns the network matrix of a network nt
    Takes a list of interaction weigths (in the same order as the list of nodes) as the "rates" argument
    """
    if type(nt) == str: 
        nt = az.transform(nt)
    else: 
        nt = az.transform(az.transform(nt))
    n = len(nt)
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            M[i,j] = rates[G.index(nt[j][0] + nt[i][-1])]
    return M

def az_kinetic_model(nt, tf, ka_rates, kb_rates, normed=True, X0=None):
    """
    Takes a network structure "nt" (list of genotypes or network structure string), a final time-point "tf" and kinetic rates (ka and kb). 
    Returns integrated solution (normalized or not) of the resulting kinetic model.
    It is assumed that we start without any WXYZ catalysts and only WXY fragments but it can be specified with "X0" argument (initial conditions for the system).
    "tf" is in minutes.
    Concentration are in uM.
    """
    def dx(x, t, A, B):
        return np.dot(A, x) + B
        
    if type(nt) == list:
        nt = transform(nt)
        
    A = network_matrix(nt, ka_rates)
    B = np.dot(network_matrix(nt, kb_rates), np.ones(A.shape[0]))
    if X0 == None:
        X0 = np.zeros(A.shape[0])
    
    X = odeint(lambda x, t: dx(x, t, A, B), X0, np.linspace(0, tf, 1000))
    if normed:
        return az.norm_array(X[-1])
    else:
        return X[-1]
    
def eig_centrality(nt, rates):
    """
    Returns the eigenvector centrality for a network nt.
    Takes a list of interaction weigths (in the same order as the list of nodes) as the "rates" argument.
    """
    if type(nt) == list:
        nt = az.transform(nt)
        M = network_matrix(nt, rates=rates)
    elif type(nt) == np.ndarray:
        M = nt
        nt = None
    else:
        M = network_matrix(nt, rates=rates)
    L, V = np.linalg.eig(M)
    kmax = np.real(L).argmax()
    return (np.real(V[:,kmax])/np.real(V[:,kmax]).sum(), az.transform(nt))