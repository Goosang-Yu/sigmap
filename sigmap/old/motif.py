import pandas as pd
import numpy as np
import re
from collections import Counter
from sigmap.const import dict_pstnp, list_best_features
from sigmap.utils import make_dimers_dict, make_trimers_dict



def calculate_dimer_counts(nuc, gap=0):
    """
    Calculate dimer counts for a given nucleotide sequence.
    Supports both contiguous dimers and gapped dimers.

    Parameters:
    - nuc (str): Nucleotide sequence.
    - gap (int): Gap between the two nucleotides in a dimer. Default is 0 (contiguous).

    Returns:
    - list: A list of counts for each dimer.
    """
    dim = [0] * 16
    dimers_dict = make_dimers_dict()

    for i in range(len(nuc) - (gap + 1)):
        dimer = nuc[i] + nuc[i + gap + 1]
        pos = dimers_dict[dimer]
        dim[pos] += 1

    return dim

def trimer_lgap(nuc,g):
    tri=[0 for i in range(64)]
    for i in range(len(nuc)-(g+2)):
        tri1=nuc[i]+""+nuc[i+g+1]+""+nuc[i+g+2]
        dict_pos_tri = make_trimers_dict()
        pos=dict_pos_tri[tri1]
        tri[pos]+=1
    return tri

def calculate_trimer_counts(nuc, gap=0, frame=1):
    """
    Calculate trimer counts for a given nucleotide sequence.
    Supports contiguous and gapped trimers.

    Parameters:
    - nuc (str): Nucleotide sequence.
    - gap (int): Gap between the nucleotides in a trimer. Default is 0 (contiguous).
    - frame (int): Frame of position of trimer for making gap. Default is 1. 1 or 2

    Returns:
    - list: A list of counts for each trimer.
    """
    tri = [0] * 64
    trimers_dict = make_trimers_dict()
    
    if frame == 1: 
        for i in range(len(nuc) - (gap + 2)):
            trimer = nuc[i] + nuc[i + 1] + nuc[i + gap + 2]
            pos = trimers_dict[trimer]
            tri[pos] += 1
            
    elif frame == 2: 
        for i in range(len(nuc) - (gap + 2)):
            trimer = nuc[i] + nuc[i + gap + 1] + nuc[i + gap + 2]
            pos = trimers_dict[trimer]
            tri[pos] += 1
    else:
        raise ValueError('Not available frame.')

    return tri

def pattern(nuc,pattern):
    return len(re.findall(pattern, nuc))

def GC_skew(nuc):
    g=nuc.count('G')
    c=nuc.count('C')
    return (c-g)/(c+g)

def AT_skew(nuc):
    a=nuc.count('A')
    t=nuc.count('T')
    return (a-t)/(a+t)

def pstnp(seq):
    """
    Calculate PSTNP values for the given sequence.

    Parameters:
    - seq (str): Input nucleotide sequence.

    Returns:
    - list: A list of PSTNP values.
    """
    return [dict_pstnp[seq[i:i+3]][i] for i in range(79)]

def eiip(seq):
    
    dict_val={'A': 0.1260,'C': 0.1335, 'G': 0.0806,'T': 0.1340}
    
    l1=[]
    dict2={}
    dict1={}
    dict_pos_tri = make_trimers_dict()
    
    for pat in dict_pos_tri.keys():
        dict1[pat]=0
        dict2[pat]=0
    for i in range(79):
        var=seq[i]+""+seq[i+1]+""+seq[i+2]
        dict1[var]+=1
    for i in range(79):
        etemp=dict_val[seq[i]]+dict_val[seq[i+1]]+dict_val[seq[i+2]]

        #Newly Added
        var=seq[i]+""+seq[i+1]+""+seq[i+2]

        res=etemp*(dict1[var])/79
        dict2[var]=res
        
    for k in dict2:
        l1.append(dict2[k])
        
    return l1


def motif_features(cdk:pd.DataFrame):
        
    #creating the first dataframe for initial columns (features)
    list_of_col_names=["f"+""+str(i) for i in range(0, 1598)]
    df_first=pd.DataFrame(columns=list_of_col_names)

    dimers_dict  = make_dimers_dict()
    dict_pos_tri = make_trimers_dict()

    #creating dataframe for independent dataset
    list_of_seq_all = cdk['Seq'].tolist()

    for i in range(len(list_of_seq_all)):
        list_to_add=[]
        nuc=str(list_of_seq_all[i])
        
        # Count all nucleotide occurrences
        nucleotide_counts = Counter(nuc)
        list_to_add.extend([nucleotide_counts.get(base, 0) for base in 'ATGC'])
        
        # Count dimer and trimer motif
        list_dim_wg=calculate_dimer_counts(nuc)
        list_trim_wg=calculate_trimer_counts(nuc)
        list_dimer_gaps = [calculate_dimer_counts(nuc, gap=g) for g in range(1, 6)]
        list_to_add += list_dim_wg + list_trim_wg + sum(list_dimer_gaps, [])
        
        # Calculate trimer counts for rgap and lgap in one loop each
        gap_values = [1, 2, 3, 7, 8, 9, 10, 15, 16, 17]
        list_trimer_rgap = [calculate_trimer_counts(nuc, gap=gap, frame=1) for gap in gap_values]
        list_trimer_lgap = [calculate_trimer_counts(nuc, gap=gap, frame=2) for gap in gap_values]
        list_to_add += sum(list_trimer_rgap + list_trimer_lgap, [])

        # Define patterns to count
        patterns = [
            "TTGAC", "TATAAT", "TTATAA", "TTGACA", "AACGAT", "ACAGTT", "AGGAGG", "TAAAAT", "TTGATT"
        ]

        # Calculate counts for patterns and skew
        pattern_counts = [pattern(nuc, pat) for pat in patterns]
        skew_counts = [GC_skew(nuc), AT_skew(nuc)]
        list_to_add.extend(pattern_counts + skew_counts)


        list_new=pstnp(nuc)
        list_new2=eiip(nuc)

        list_to_add+=list_new+list_new2

        final_list=np.array(list_to_add)
        df_first.loc[i]=final_list

    features_list2=["f"+str(i) for i in list_best_features if i< 1598]

    df_first=df_first[features_list2]
    
    return df_first