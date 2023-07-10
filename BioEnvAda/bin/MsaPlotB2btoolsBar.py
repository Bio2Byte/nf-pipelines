import os
import numpy as np
import re
import sys
import math
import matplotlib.pyplot as plt
import pandas as pd

os.environ['MPLCONFIGDIR'] = os.getcwd() + "/configs/"
import matplotlib as mpl
mpl.use("Agg")


msa = sys.argv[1]
predictions = sys.argv[2] 
stats = sys.argv[5]


split = sys.argv[6]
if split == 'true':
    split = 2000
else:
    split = int(split)


try:
    tools =sys.argv[3].replace(' ','').split(",")
except:
    tools =[]

selected_proteins = []
if sys.argv[4]!= 'false':
    try:
        selected_proteins = sys.argv[4].split(",")
    except:   
        raise ValueError("Selected proteins are not recognized")



#set plot titles and tool list for plotting
DynaMine = True
PREDICTION_TITLES = {}
if DynaMine == True:
    PREDICTION_TITLES = {
        'backbone': "DynaMine backbone dynamics",
        'sidechain': "DynaMine sidechain dynamics",
        'ppII': "DynaMine conformational propensities: ppII (polyproline II)",
        'coil': "DynaMine conformational propensities: Coil",
        'sheet': "DynaMine conformational propensities: Sheet",
        'helix': "DynaMine conformational propensities: Helix"}
if 'efoldmine' in tools:
    PREDICTION_TITLES['earlyFolding'] = "Early folding (EFoldMine)"
if 'disomine' in tools:
    PREDICTION_TITLES['disoMine'] = "Disorder (disoMine)"
if 'agmata' in tools:
    PREDICTION_TITLES['agmata'] = "Beta aggregation (AgMata)"

AXIS_TITLES = {
    "x": "Residue position in the MSA",
    "y": "Prediction values"}


#get occupancy


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

seq_df = pd.DataFrame.from_dict(sequences_dict, orient='index',  columns=['seq'])
seq_df = pd.DataFrame(seq_df.seq.apply(list).tolist())


SEQUENCES_COUNT, RESIDUES_COUNT = seq_df.shape
OCCUPANCY = 1 - (seq_df == '-').sum() / SEQUENCES_COUNT

#read json
data = pd.read_json(predictions)
data.drop(labels = 'sequence', axis =1, inplace =True)
pred_df= pd.DataFrame()

for pred in PREDICTION_TITLES.keys():
    dat = data.loc[[pred]]
    pred_df=pd.concat([pred_df, dat])

#get selected protein data 
full_label_selp = pd.DataFrame()
if selected_proteins != []:
    for i in selected_proteins:
        sel_df = pred_df.filter(like=i) 
        rows,cols = sel_df.shape
        if cols > 1:
            raise ValueError("Protein identifier is not unique:", i)
        if sel_df.empty:
            raise ValueError("Selected protein not found in MSA:", i)   
        else:
            if full_label_selp.empty:
                full_label_selp = sel_df
            else:
                full_label_selp=pd.concat([full_label_selp,sel_df],axis =1)

    
def plot_biophysical_msa(stats, suboccupancy, selected_prot , biophys_data, label, res_num):

    #Plot representation
    fig, (ax1,ax2) = plt.subplots(2,sharex=True, gridspec_kw={'height_ratios': [10, 1]})
    fig.set_figwidth(20)
    fig.set_figheight(10)
    
    figlabel = 'Predicted biophysical properties of the MSA: %s aligned residues from %s sequences \n Section %s of %s' %(RESIDUES_COUNT,SEQUENCES_COUNT, str(label[0]), str(label[1]))
    plt.suptitle(figlabel, fontsize=14)

    colors = plt.cm.tab10
    #remove nan values
    df = stats.dropna(axis=1)

    firstq = df['firstQuartile'].tolist()
    thirdq = df['thirdQuartile'].tolist()
    bottom = df['bottomOutlier'].tolist()
    top = df['topOutlier'].tolist()

    x=np.arange(0,len(bottom),1)
    ax1.fill_between(x, firstq, thirdq, alpha=0.3, color=colors(0), label='1st-3rd Quartiles')
    ax1.fill_between(x, bottom, top, alpha=0.1, color=colors(0), label='Outliers')
    ax1.plot(df['median'].tolist(), linewidth=1.25, color=colors(0), label='Median')
    
    ax2.bar(range(0,len(bottom)), suboccupancy,  width=1, color=colors(0))
    ax1.axis([0,len(bottom), min(bottom)-0.05, max(top)+0.05])

    #add residue numbering
    ticksevery = 20
    plt.xticks(np.arange(0,len(bottom),ticksevery),range(res_num[0],res_num[1],ticksevery))

    #Add cutoffs for predictions
    if biophys_data == 'backbone':
        ax1.axhline(y=1.0, linewidth= 1.5, linestyle='-.', label='Above: Membrane spanning', color=colors(0)) # color='green'
        ax1.axhline(y=0.8, linewidth= 1.5, linestyle='--', label='Above: Rigid', color=colors(0)) #color='orange'
        if min(bottom)-0.05 < 0.69:
            ax1.axhline(y=0.69, linewidth= 1.5, linestyle=':', label='Above: Context dependent Below: Flexible', color=colors(0)) #, color='red'
    if biophys_data == 'earlyFolding':
        ax1.axhline(y=0.169, linewidth= 1.5, linestyle='-.', label='Above: Likely to start folding', color=colors(0)) #, color='red'
    if biophys_data == 'disoMine':
        ax1.axhline(y=0.5, linewidth= 1.5, linestyle='-.', label='Above: Likely to be disordered', color=colors(0)) #, color='red'

    #Add the selected protein if there are some
    if selected_prot.empty == False:
        header = selected_prot.columns.tolist()
        for c in range (0,len(header)):
            row = selected_prot.loc[k, header[c]]
            ax1.plot(row, '-s', linewidth=1.5, color=colors(1+c), label=f'Prediction {selected_proteins[c]}', markersize=2)
        plot_name = re.sub(r'[^a-zA-Z0-9]', '_',selected_proteins[0] + '_' + PREDICTION_TITLES[biophys_data] +'_' + str(label[0])+'of'+str(label[1]) )
    else:
        plot_name = re.sub(r'[^a-zA-Z0-9]', '_', PREDICTION_TITLES[biophys_data] +'_' + str(label[0])+'of'+str(label[1])  )

    #add legend
    plt.figlegend(loc='lower center', ncol =4)
    fig.subplots_adjust(top=0.9, hspace = 0.0001)

    ax1.set_title(PREDICTION_TITLES[biophys_data])
    ax1.set_ylabel(AXIS_TITLES['y'])
    ax2.set_xlabel(AXIS_TITLES['x'])
    ax2.set_ylabel('Occupancy')

    msa_baseName = os.path.split(msa)[1].split('.',1)[0]
    plt.savefig( msa_baseName + '_msa_' + plot_name + '.png', dpi=300)


keys = list(PREDICTION_TITLES.keys())

df = pd.read_json(stats)
for k in keys:
    stats_df = pd.DataFrame.from_records(df['results'][k])

    if RESIDUES_COUNT > split:
        subplots = RESIDUES_COUNT/ split
        subplots = math.ceil(subplots)
        for elem in range(1,subplots+1):
            label = [elem, subplots]
            low = int(RESIDUES_COUNT * ((1/subplots)*(elem-1)))
            high = int(RESIDUES_COUNT * ((1/subplots)*elem))

            substats_df = stats_df.iloc[low:high]
            sub_occupancy = OCCUPANCY[low:high]
            res_num = [low, high]

            sub_label_selp = full_label_selp.iloc[low:high]

            print ('plotting:' , sub_label_selp )


            plot_biophysical_msa(substats_df, sub_occupancy, sub_label_selp, k, label , res_num)
    else:
        label = [1,1]
        res_num = [0, len(OCCUPANCY)]
        plot_biophysical_msa(stats_df,OCCUPANCY, full_label_selp, k, label, res_num)