import pandas as pd

def make_dimers_dict():
    
    dimers_dict={}
    bases = "ACGT"
    pos=0
    for i in range(len(bases)):
        for j in range(len(bases)):
            dimers_dict[bases[i]+""+bases[j]]=pos
            pos+=1
            
    return dimers_dict


def make_trimers_dict():
    
    list_of_tri=[]
    dict_pos_tri={}
    bases = "ACGT"
    pos=0
    for i in range(len(bases)):
        for j in range(len(bases)):
            for k in range(len(bases)):
                tri=bases[i]+""+bases[j]+""+bases[k]
                list_of_tri.append(tri)
                dict_pos_tri[tri]=pos
                pos+=1
                
    return dict_pos_tri

def fasta2df(inFile:str) -> pd.DataFrame:
    # Extract file extension and initialize cdk DataFrame
    with open(inFile, "r") as f:
        lines = f.readlines()

    seq, s_id, s = [], [], ""
    for line in lines:
        if line.startswith('>'):
            if s:
                seq.append(s)
            s_id.append(line.strip())
            s = ""
        else:
            s += ''.join([char.capitalize() for char in line if char.upper() in {'A', 'C', 'G', 'T'}])

    if s:
        seq.append(s)
    
    df = pd.DataFrame()

    df['Sequence'] = seq
    df['Sequence_ID'] = s_id
    
    return df