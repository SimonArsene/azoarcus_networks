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

def levenshtein(s, ss):
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