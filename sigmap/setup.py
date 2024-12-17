import pandas as pd
import re


def readseq(input_data, scan=False) -> pd.DataFrame:
    
    if type(input_data) == str:
        with open(input_data) as f:
            records = f.read()
            
        records    = records.split('>')[1:]
        list_seq   = []
        list_seqid = []
        
        for fasta in records:
            array = fasta.split('\n')
            name, sequence = array[0].split()[0], re.sub('[^ATGCU-]', '', ''.join(array[1:]).upper())
            list_seqid.append('>'+name)
            list_seq.append(sequence)
            
        if len(list_seqid) == 0:
            f = open(input_data,"r")
            data1 = f.readlines()
            for each in data1:
                list_seq.append(each.replace('\n',''))
            for i in range (1,len(list_seq)+1):
                list_seqid.append(">Seq_"+str(i))
                
        final_df = pd.concat([pd.DataFrame(list_seqid),pd.DataFrame(list_seq)], axis=1)
    
    elif type(input_data) == pd.DataFrame:
        final_df = input_data.iloc[:, :2]
    
    else:
        raise TypeError('Not supported input type.')

    final_df.columns = ['Sequence_ID','Sequence']
    
    list_seq_processed = []
    
    if scan == True:
        for seq in final_df['Sequence']:
            if len(seq) < 81:
                buffer = 81 - len(seq)
                seq = seq + 'A'*buffer
            list_seq_processed.append(seq)
    else:
        for seq in final_df['Sequence']:
            if len(seq) < 81:
                buffer = 81 - len(seq)
                seq = seq + 'A'*buffer
            elif len(seq) > 81:
                seq = seq[:81]
            list_seq_processed.append(seq)
            
    final_df['Sequence'] = list_seq_processed
                      
    return final_df


def seq_pattern(data):
    win_len = 81
    list_seq_pattern = []
    
    for s in data['Sequence']:
        seq = s.upper()
        seq_pat=[]
        for pos in range(0, (len(seq) + 1 - win_len)):
            seq_pat += [seq[pos:pos+win_len]]
        list_seq_pattern.append(seq_pat)
    return list_seq_pattern
  
    
def seq_mutants(data):
    std = list("ATGC")
    list_mut = []
    for s in data['Sequence']:
        seq = s.upper()
        mut_seq=[]
        for pos in range(0,len(seq)):
            for nt in std:
                mut_seq += [seq[:pos] + nt + seq[pos+1:]]
        list_mut.append(mut_seq)
    return list_mut
  
def scanner(file1) -> pd.DataFrame:
    df_seq = readseq(file1, scan=True)
    bb = seq_pattern(df_seq)
    flat = [x for sublist in bb for x in sublist]
    df2 = pd.DataFrame(flat)
    cc = []
    dd = []
    ee = []
    for i in range(0,len(df_seq)):
        cc.append(len(df_seq['Sequence'][i]))
        if len(df_seq['Sequence'][i]) <= 81:
            dd.append(1)
        else:
            dd.append(abs(81-len(df_seq['Sequence'][i]))+1)
    df_seq['length'] = cc
    df_seq['Extra'] = dd
    for i in range(0,len(df_seq)):
        for j in range(0,df_seq.Extra[i]):
            ee.append(df_seq.ID[i])
    df3 = pd.concat([pd.DataFrame(ee),df2],axis=1)
    
    return df3


  
def designer(file1) -> pd.DataFrame:
    df_seq = readseq(file1)
    list_mut = seq_mutants(df_seq)
    flat = [x for sublist in list_mut for x in sublist]
    df2 = pd.DataFrame(flat)
    ee = []
    
    for i in range(0,len(df_seq)):
        for j in range(0,324):
            ee.append(df_seq['Sequence_ID'][i])
    df3 = pd.concat([pd.DataFrame(ee),df2],axis=1)
    
    return df3


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

    df['Sequence_ID'] = s_id
    df['Sequence'] = seq
    
    return df