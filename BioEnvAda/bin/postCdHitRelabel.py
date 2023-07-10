import pandas as pd
import sys


seqs = sys.argv[1]
clusters_file = sys.argv[2]
relabel = sys.argv[3]


def cd_hit_to_df(clusters_file):
    #read cdhit
    with open(clusters_file, 'r') as reps:
        rep = 'a'
        cluster= ['a']
        match_dic = {}
        rep_list=[]
        for line in reps:
            if line.startswith('>Cluster'):
                match_dic[rep] = cluster
                cluster=[]
                rep='a'
            else:
                id=line.split('>')[1].split('..')[0]
                cluster.append(id)
                if '*' in line:
                    #keys need to be unique -- genomeID_annotation
                    rep=id
                    rep_list.append(rep)
    match_dic[rep] = cluster
    match_dic.pop('a')

    print (match_dic.keys())
    #{'rep id':['id','id','id']}
    if len(match_dic) ==  len(rep_list):
        print('all cluster representatives collected:', len(match_dic))
    else:
        print('found reps: ', len(rep_list))
        print ('all keys: ', len(match_dic))
        print (rep_list)

    #gene_df=pd.DataFrame.from_dict(match_dic, orient='index')
    gene_df= pd.Series(match_dic).to_frame('seqIds') #df: index=rep, seqIds=[id,id,id]
    gene_df= gene_df.reset_index()

    return(gene_df)

def make_label(row, data_df):
#label MUST be shorter then 45 chars
    labels = data_df.loc[row,'seqIds']
    short_labels = []

    rep = data_df.loc[row,'index']

    max_label_len = 42-len(rep)
    no_of_ids = len(labels)
    indiv_id_len = int(max_label_len/no_of_ids)-1

    for elem in labels:
        short_label = elem
        if indiv_id_len < 2:
            short_label = '...' + str(len(labels))
        elif indiv_id_len < 6:
            if "Syn" in elem:
                short_label = elem.replace('Syn_', '')
            if 'Cya' in elem: 
                short_label = elem.replace('Cya_', '')
            short_label = short_label[:indiv_id_len]
        else:
            if "Syn" in elem:
                short_label = elem.replace('Syn', 'S')
            if 'Cya' in elem: 
                short_label = elem.replace('Cya', 'C')
            short_label = short_label[:indiv_id_len]
        
        if short_label[-1] == '_':
            short_label = short_label[:-1]
        
        if elem in data_df.loc[row,'index']:
            continue
        else:
            short_labels.append(short_label)

    if rep[-1] == '_':
        rep = rep[:-1]
    if rep[-1] == 'C' and rep[-2] == '_':
        rep = rep[:-2]
    if rep[-1] == 'K' and rep[-2] == 'C' and rep[-3] == '_':
        rep = rep[:-3]


    short_labels = list(set(short_labels))
    short_labels.insert(0,rep)
    label = '-'.join(short_labels)
    
    return(label)

def cluster_fasta(rep_df, seq_file):  
    print ('start cluster fasta')

    with open(seq_file, 'r') as seqs:
        label = 'a'
        seq= ['a']
        match_dic = {}
        for line in seqs:
            if line.startswith('>'):
                match_dic[label] = ''.join(seq)
                seq=[]
                label=line.replace('\n','').replace('>','')
            else:
                part=line.replace('\n','')
                seq.append(part)

        match_dic[label] =  ''.join(seq)
        match_dic.pop('a')

    #print (match_dic)
    seq_df=pd.DataFrame.from_dict(match_dic, orient='index', columns=['seqs'])
    seq_df= seq_df.reset_index()
    seq_df['index'] = seq_df['index'].str[:39]  
    print (seq_df)

    print ("blah")

    print (rep_df)
    df = pd.merge(rep_df, seq_df, on='index')
    print (df)

    if relabel == 'true':
        out_df = df[['cd_label','seqs']]
    else:
        out_df = df[['index','seqs']]



    
    #remove stopcodons for BuildTreeEvol
    stopcodons = ['TAG','TAA','TGA']
    SEQUENCES_COUNT =  out_df['seqs'].shape[0]

    for row in range(0, SEQUENCES_COUNT):
        fullseq= out_df['seqs'].iloc[row].replace('\n','').replace(' ','')
        codon = fullseq [-3:]
        if codon in stopcodons:
            out_df['seqs'].iloc[row] = fullseq[:-3]
            print("Stop codons were found in alignment in sequence %s. This causes problems for iqtree and therefore the last 3 nucleotides have been removed."%(row))
    
        if fullseq[-1:] == '*':
            out_df['seqs'].iloc[row] = fullseq[:-1]
            print("* were found in alignment in sequence %s. This causes problems for iqtree and therefore the last  element has been removed."%(row))
    


    #prep and write outfile

    file_extension = len(seq_file.split('.')[-1]) +1
    if relabel == 'true' :
        out_name = seq_file[: -file_extension]+'_clustered_relabeled.fasta'
    else:
        out_name = seq_file[: -file_extension]+'_clustered.fasta'


    out_df.iloc[:,0]= '>' + out_df.iloc[:,0] + '\n'

    rows = out_df.to_string(header=False,index=False,index_names=False).split('\n')

    with open (out_name, 'w') as m:
        for row in rows:
            row = row.replace('\\n','\n').replace(' ','')
            m.write( row + '\n')



df = cd_hit_to_df(clusters_file)
#print (df)


for row in range(len(df)):
    df.loc[row,'cd_label'] = make_label(row, df)


outname= seqs + '_representative_represented.csv'

df.to_csv(outname)


cluster_fasta(df, seqs)

