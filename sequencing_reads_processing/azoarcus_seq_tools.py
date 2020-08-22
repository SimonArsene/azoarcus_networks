from Bio import SeqIO
import numpy as np
import dpBC_tools as az

"""
General strategy: 
1) junction_position > tag > genotype
2) genotype > insert_at_junctions > insert
"""

meta_az_header = "seq_id,G,I,x,y,z,UMI,is_hp,hp"
meta_dpBC_header = "seq_id_B,seq_id_C,seq_id_D,dpBC"

def get_meta_az(rcd, start_read2, nb_indexes, use_lev=1):
    """
    Returns meta data associated with a sequence: seq_id,G,I,x,y,z,UMI,is_hp,hp
    """
    j_positions = junction_position(rcd, use_lev=use_lev)
    tag_rcd = tag(rcd, j_positions)
    G = genotype(rcd, j_positions, tag_rcd, use_lev=use_lev)
    I = insert(rcd, insert_at_junctions(rcd, j_positions, tag_rcd), G)
    umi = umi_rd2(rcd, start_read2, nb_indexes, use_lev=use_lev)
    return [str(x) for x in ([rcd.id] + [G, I] + j_positions + [umi] + hp(rcd, use_lev=use_lev))]
    

def qlength(rcd):
    """
    Returns "quality length"
    """
    qscore = rcd.letter_annotations["phred_quality"]
    nq = np.array(qscore)
    q_mean = []
    for i in range(1,len(nq)):
        q_mean.append(abs(nq[:i].mean() - nq[i:].mean()))
    q_mean[0] = q_mean[1]
    return q_mean.index(max(q_mean))
    
def junction_position(rcd, use_lev):
    """
    Searchs for junctions positions
    Usual parameters: window of search = 20 and distance max = 2
    """
    x = "TAAGCGCCCG"
    y = "GCTTCGGCGC"
    z = "ACCTAAGGCA"
    x_pos = 59
    y_pos = 95
    z_pos = 141
    window = 20
    dist_max = 2
    
    return [az.approx_pos(x, str(rcd.seq), x_pos, window, dist_max, use_lev=use_lev), 
            az.approx_pos(y, str(rcd.seq), y_pos, window, dist_max, use_lev=use_lev), 
            az.approx_pos(z, str(rcd.seq), z_pos, window, dist_max, use_lev=use_lev)]
    
def tag(rcd, j_positions):
    """
    Looks at each junction if a TAG sequence can be found
    Returns [TAGx, TAGy, TAGz]
    """
    seq = str(rcd.seq)
    
    list_tag = ["CAT", "CCT", "CTT", "CGT"]
    tag = ["-","-","-"]
    
    for i in range(3):
        if j_positions[i] != -1:
            if seq[j_positions[i] + 10:j_positions[i] + 13] in list_tag:
                tag[i] = seq[j_positions[i] + 10:j_positions[i] + 13][-2]
    return tag
    
    
def genotype(rcd, j_positions, tag_rcd, use_lev):
    """
    Returns the full genotype formatted
    """
    junctions = ["x", "y", "z"]
    genotype = ["-","-","-"]
    
    if tag_rcd.count("-") == 2:
        for e in tag_rcd:
            if e != "-":
                genotype[2] = e
                genotype[1] = junctions[tag_rcd.index(e)]
    elif tag_rcd.count("-") == 1:
        if tag_rcd[0] == "C" and tag_rcd[1] != "-":
            genotype[2] = tag_rcd[1]
            genotype[1] = junctions[1]
        elif tag_rcd[0] == "C" and tag_rcd[2] != "-":
            genotype[2] = tag_rcd[2]
            genotype[1] = junctions[2]

    genotype[0] = igs(rcd, use_lev=use_lev)
    
    return "".join(genotype)

def insert_at_junctions(rcd, j_positions, tag):
    """
    Look at each junction if an insert can be found
    Returns [Ix, Iy, Iz]
    """
    seq = str(rcd.seq)
    
    insert = ["-","-","-"]
    for j in range(3):
        if tag[j] != "-":
            insert[j] = insert_aux(seq[j_positions[j] + 13:j_positions[j] + 18])
    return insert
        
   
def insert(rcd, insert_at_junctions, genotype):
    """
    Given the genotype (ie the junction that was chosen), returns the corresponding insert
    """
    junctions = ["x", "y", "z"]
    if genotype[1] != "-":
        j = junctions.index(genotype[1])
        return insert_at_junctions[j]
    else:
        return "-"

def insert_aux(seq):
    """
    Auxiliary function for insert_at_junctions
    Returns the first element encountered
    """
    list_insert = ["GGCAT", "GCAT", "CAT"]
    for insert in list_insert:
        if insert in seq:
            return insert
    return "-"
    
def igs(rcd, use_lev):
    """
    Returns the IGS
    Used by genotype function
    """
    if (type(rcd) != str):
        seq = str(rcd.seq)
    else:
        seq = rcd
    igs_seq = "GGGAAACCACG"
    igs_list = ["CAA", "CTA", "CCA", "CGA"]
    pos = az.approx_pos(igs_seq[4:], seq, 23, 5, 2, use_lev=use_lev)
    if pos != -1:
        igs = seq[pos+7:pos+10]
        if igs in igs_list:
            return igs[1]
    return "-"
    
def polyA(rcd):
    """
    Returns the distance of pA
    Usual parameters: pA = "AAAAA"
    """
    pA = "AAAAA"
    return rcd.seq[:301].find(pA)  
    
def polyG(rcd):
    """
    Returns the distance of pA
    Usual parameters: pA = "AAAAA"
    """
    pG = "GGGGG"
    return rcd.seq[:301].find(pG) 
    
def rBC(rcd, pA):
    rBC_list = ["ATAGCC", "CGAACC", "AGTGCC", "AAGGCC", "TCCGCC"]
    if pA > 6:
        return az.approx_search(str(rcd.seq[pA-6:pA]), rBC_list, 1)
    else:
        return -1
        
def umi_rd2(rcd, start_read2, nb_indexes, use_lev):
    uni_linker = "TACGCTACGGAACGA"
    p = az.approx_pos(uni_linker, str(rcd.seq[start_read2:]), 20*nb_indexes+4, 10, 2, use_lev=use_lev)
    if p != -1 and str(rcd.seq[start_read2+p+15-3:start_read2+p+15]) == "CGA":
        return str(rcd.seq[start_read2+p+15:start_read2+p+15+8])
    else:
        return "-"

#Out-dated because it uses pA distance        
def umi_rd1(rcd, pA):
    uni_linker_rev = "TCGTTCCGTAGCGTA"
    p = az.approx_pos(uni_linker_rev, str(rcd.seq[:301]), pA + 28, 10, 2)
    if p != -1 and str(rcd.seq[p-8-5:p-8]) == "AAAAA":
        return str(rcd.seq[p-8:p])
    else:
        return "-"

#Out-dated because bulk primer now has 3 indexes
def umi_bulk(rcd, start_read2):
    return str(rcd.seq[start_read2:start_read2+8])
    
def hp(rcd, use_lev):
    hp_list = ["TTTT","ACTG","AATT","GATA","GTGG","TATG","CCCA","AAGG","GGAA","ATAG",
               "TAAC","TGTA","TACA","CGTG","GCCG","GAAT","GGGT","CGCC","TTGG","GTCC",
               "AAAA","CTAA","ATGC","TGGT","GTAG"]
    before = "AAACCAGTTG"
    after = "CAACTGGTTT"
    p_before = az.approx_pos(before, str(rcd.seq), 21, 5, 2, use_lev=use_lev)
    p_after = az.approx_pos(after, str(rcd.seq), 35, 5, 2, use_lev=use_lev)
    is_hp = 0
    if p_before != -1 and p_after != -1:
        is_hp = 1
    hp_seq = "----"
    if str(rcd.seq[p_before+7:p_before+10]) == "TTG" and p_after == p_before + 14:
            hp_seq = str(rcd.seq[p_before+10:p_after])
    return [str(is_hp), str(az.approx_search(hp_seq, hp_list, 1))]