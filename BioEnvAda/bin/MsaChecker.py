import pandas as pd
import sys

msa = sys.argv[1] 
buildTreeEvo =  sys.argv[2] 
drop_empty = sys.argv[3]

made_new_file = False

with open(msa, 'r') as f:
    lines = f.readlines()

sequences_dict = {}
name = 'a'
sequence = 'a'
for line in lines:
    if line.startswith('>'):
        sequences_dict[name] = sequence
        name = line.strip('>\n')
        sequence = ''
    else:
        sequence += line.strip()


sequences_dict[name] = sequence
sequences_dict.pop('a')

msa_df = pd.DataFrame.from_dict(sequences_dict, orient='index',  columns=['seq'])

seq_df = pd.DataFrame(msa_df.seq.apply(list).tolist())
print(seq_df)

#check if seqs have all the same length
if seq_df.isnull().values.any() == True:
    raise ValueError("Sequences do not have all the same length, is this an MSA?")


#remove columns without minimal occupancy 
SEQUENCES_COUNT, RESIDUES_COUNT = seq_df.shape
OCCUPANCY = 1 - (seq_df == '-').sum() / SEQUENCES_COUNT

only_gap = []
cols = seq_df.columns

if drop_empty == 'false':
    print ("No occupancy check performed.")
else:
    for i in range (0,RESIDUES_COUNT):
        if drop_empty == 'true' :
            if OCCUPANCY[i]== 0:
                only_gap.append(cols[i])
        else:
            if OCCUPANCY[i] < float(drop_empty):
                only_gap.append(cols[i])
    seq_df.drop(labels = only_gap, axis =1, inplace =True)
    if only_gap!=[]:
        print ("Low occupancy column in MSA! Columns %s were removed and new file was created." %(only_gap))
        made_new_file =True 



#remove stopcodons for BuildTreeEvol
found_stopcodons=False
if buildTreeEvo == 'true':
    stopcodons = ['TAG','TAA','TGA']
    lastcodon_df = seq_df.iloc[:, -3:]

    for row in range(0, SEQUENCES_COUNT):
        if found_stopcodons == True:
            break

        codon= lastcodon_df.iloc[row].to_string(header=False, index=False).replace('\n','').replace(' ','')
        if codon in stopcodons:
            seq_df = seq_df.iloc[:, :-3]
            print("Stop codons were found in alignment in sequence %s. This causes problems for iqtree and therefore the last 3 nucleotides have been removed."%(row))
            made_new_file = True 
            found_stopcodons = True


#remove stars
seq_df.iloc[:,-1] = seq_df.iloc[:,-1].str.replace('*', '')


#add headers
#seq_ids = pd.DataFrame(alignment_file.iloc[::2].values, columns=['index'] )

msa_df=msa_df.reset_index()
seqid_df = msa_df['index']
seq_df = pd.concat([seqid_df,seq_df], axis=1)

print (seq_df)

#create_short_labels
#seq_df.ids = seq_df.ids.str.split('_CK', n=1, expand=True)[0]

#prep and write outfile
seq_df['index']= ">" + seq_df['index'] + '\n'
rows = seq_df.to_string(header=False,index=False,index_names=False).split('\n')

file_extension = len(msa.split('.')[-1]) +1

out_name = msa[: -file_extension]+'_checked.' + msa.split('.')[-1]

with open (out_name, 'w') as m:
    for row in rows:
        row = row.replace('\\n','\n').replace(' ','')
        m.write( row + '\n')
        print (row)