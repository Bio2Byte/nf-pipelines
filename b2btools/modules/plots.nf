process plotBiophysicalFeatures {
    tag "${predictions.name}"

    input:
    path predictions
    val dynamine
    val efoldmine
    val disomine
    val agmata
    val psper

    output:
    path "*.png", emit: plots

    script:
    """
#!/usr/local/bin/python
import json
import os
os.environ['MPLCONFIGDIR'] = os.getcwd() + "/configs/"
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def plot_single_sequence(all_predictions, tools):
    keys = all_predictions.keys()
    nb_tools = 0
    if 'dynamine' in tools:
        nb_tools +=6
    if 'disomine' in tools:
        nb_tools +=1
    if 'efoldmine' in tools:
        nb_tools +=1
    if 'agmata' in tools:
        nb_tools +=1

    for seq_key in keys:
        predictions = all_predictions[seq_key]
        x_position = range(len(predictions['seq']))
        residues_count = len(predictions['seq'])

        fig, axs = plt.subplots(nb_tools)
        fig.suptitle(f'Single Sequence Predictions {seq_key}', fontsize=20)
        fig.set_figwidth(30)
        fig.set_figheight(8*nb_tools)

        number = 0

        if 'dynamine' in tools:

            ax1 = axs[number]
            ax2 = axs[number+1]
            ax3 = axs[number+2]
            ax4 = axs[number+3]
            ax5 = axs[number+4]
            ax6 = axs[number+5]
            number+=6

            backbone_pred = predictions['backbone']
            coil_pred = predictions['coil']
            sheet_pred = predictions['sheet']
            ppII_pred = predictions['ppII']
            helix_pred = predictions['helix']
            sidechain_pred = predictions['sidechain']

            ax1.plot(x_position, backbone_pred, label=seq_key)
            ax1.set_title('DynaMine backbone dynamics', fontsize=18)
            ax1.axis([0, residues_count, min(backbone_pred)-0.05, max(backbone_pred)+0.05])
            ax1.axhline(y=1.0, color='green', linewidth= 1.5, linestyle='-.', label='Above: Membrane spaning') #Membrane spaning
            ax1.axhline(y=0.8, color='orange', linewidth= 1.5, linestyle='-.', label='Above: Rigid') #Membrane spaning
            if min(backbone_pred)-0.05 < 0.69:
                ax1.axhline(y=0.69, color='red', linewidth= 1.5, linestyle='-.', label='Above: Context dependent \\\nBelow: Flexible') #context dependent (either rigide or flexible)
            ax1.legend(ncol=1, bbox_to_anchor =(1.1,0.5), loc='center left', fontsize=15)

            ax2.plot(x_position, sidechain_pred, label=seq_key)
            ax2.set_title('DynaMine sidechain dynamics', fontsize=18)
            ax2.axis([0, residues_count, min(sidechain_pred)-0.05, max(sidechain_pred)+0.05])
            ax2.legend(ncol=1, bbox_to_anchor =(1.1,0.5), loc='center left', fontsize=15)

            ax3.plot(x_position, coil_pred, label=seq_key)
            ax3.set_title('DynaMine conformational propensities: Coil', fontsize=18)
            ax3.axis([0, residues_count, min(coil_pred)-0.05, max(coil_pred)+0.05])
            ax3.legend(ncol=1, bbox_to_anchor =(1.1,0.5), loc='center left', fontsize=15)

            ax4.plot(x_position, sheet_pred, label=seq_key)
            ax4.set_title('DynaMine conformational propensities: Sheet', fontsize=18)
            ax4.axis([0, residues_count, min(sheet_pred)-0.05, max(sheet_pred)+0.05])
            ax4.legend(ncol=1, bbox_to_anchor =(1.1,0.5), loc='center left', fontsize=15)

            ax5.set_title('DynaMine conformational propensities: ppII (polyproline II)', fontsize=18)
            ax5.plot(x_position, ppII_pred, label=seq_key)
            ax5.axis([0, residues_count, min(ppII_pred)-0.05, max(ppII_pred)+0.05])
            ax5.legend(ncol=1, bbox_to_anchor =(1.1,0.5), loc='center left', fontsize=15)

            ax6.set_title('DynaMine conformational propensities: Helix', fontsize=18)
            ax6.plot(x_position, helix_pred, label=seq_key)
            ax6.axis([0, residues_count, min(helix_pred)-0.05, max(helix_pred)+0.05])
            ax6.legend(ncol=1, bbox_to_anchor =(1.1,0.5), loc='center left', fontsize=15)

        if 'efoldmine' in tools:
            ax7 = axs[number]
            number +=1
            earlyFolding_pred = predictions['earlyFolding']
            ax7.plot(x_position, earlyFolding_pred, label=seq_key)
            ax7.set_title('Early folding (EFoldMine)', fontsize=18)
            ax7.axhline(y=0.169, color='red', linewidth= 1.5, linestyle='-.', label='Above: Likely to start folding') #above: likely start protein folding process
            ax7.legend(ncol=1, bbox_to_anchor =(1.1,0.5), loc='center left', fontsize=15)

        if 'disomine' in tools:
            ax8 = axs[number]
            number +=1
            disomine_pred = predictions['disoMine']
            ax8.plot(x_position, disomine_pred, label=seq_key)
            ax8.set_title('Disorder (disoMine)', fontsize=18)
            ax8.axhline(y=0.5, color='red', linewidth= 1.5, linestyle='-.', label='Above: Likely to be disordered') #above: likely disordered
            ax8.legend(ncol=1, bbox_to_anchor =(1.1,0.5), loc='center left', fontsize=15)

        if 'agmata' in tools:
            ax9 = axs[number]
            agmata_pred = predictions['agmata']
            ax9.plot(range(len(agmata_pred)), agmata_pred, label=seq_key)
            ax9.set_title('Beta-aggregation (AgMata)', fontsize=18)
            ax9.legend(ncol=1, bbox_to_anchor =(1.1,0.5), loc='center left', fontsize=15)

        for ax in axs:
            ax.set_xlabel('Residue position', fontsize=15)
            ax.set_ylabel('Prediction values', fontsize=15)

        plt.tight_layout()

        fig.subplots_adjust(top=0.88+(nb_tools/100), hspace = 0.6)
        plt.savefig(f'${predictions.baseName}_{seq_key}.png')


def plot_psper(all_predictions):
    keys = all_predictions.keys()
    for seq_key in keys:
        predictions = all_predictions[seq_key]
        sequence = predictions['seq']
        x = range(len(sequence))

        tyr_x = []
        tyr_y = []
        arg_x = []
        arg_y = []

        for ind, residue in enumerate(sequence):
            if residue == 'Y':
                tyr_x.append(ind)
                tyr_y.append(0.5)
            if residue == 'R':
                arg_x.append(ind)
                arg_y.append(0.5)

        viterbi_path = predictions['viterbi']
        color_map = lambda step: '#ffb499' if step == ' RRM' else '#9999ff' if step == ' PLD' else '#bfbfbf' if step == ' SPACER' else '#99ff99' if step == ' OTHER' else '#fff'
        viterbi_colors = list(map(color_map, viterbi_path))

        fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, 1, gridspec_kw={'height_ratios': [1, 3, 3, 3, 3, 3]})
        fig.set_figwidth(30)
        fig.set_figheight(15)
        fig.suptitle(f"PSPer {seq_key} (PSPer score: {predictions['protein_score']:2.3}, length {len(sequence)} AA)", size=20, y=1.05)

        # AX1
        if predictions['protein_score'] >= 0.5611191832364677:
            ax1.set_title('Likely a phase-separation protein', loc='left')
        else:
            ax1.set_title('Likely not a phase-separation protein', loc='left')

        ax1.set_xlim([0, len(sequence) - 1])
        ax1.set_ylim([0, 1])
        ax1.scatter(tyr_x, tyr_y, marker='o', color='blue')
        ax1.scatter(arg_x, arg_y, marker='o', color='red')
        ax1.axes.yaxis.set_visible(False)
        ax1.set_xlabel('Residue position')

        # Complexity
        ax2.plot(x, predictions['complexity'], color='#808080')
        ax2.set_title('Complexity', loc='left')
        ax2.set_ylabel('Seq. complexity')
        ax2.set_xlabel('Residue position')
        ax2.set_xlim([0, len(sequence) - 1])
        ax2.set_ylim([0, max(predictions['complexity']) + 1])
        ax2.grid()

        # Tyr Enrichment
        ax3.plot(x, predictions['tyr'], color='#808080')
        ax3.set_title('Tyr', loc='left')
        ax3.set_ylabel('Tyr Enrichment')
        ax3.set_xlabel('Residue position')
        ax3.set_xlim([0, len(sequence) - 1])
        ax3.set_ylim([0, max(predictions['tyr']) + 1])
        ax3.grid()

        # Arg Enrichment
        ax4.plot(x, predictions['arg'], color='#808080')
        ax4.set_title('Arg', loc='left')
        ax4.set_ylabel('Arg Enrichment')
        ax4.set_xlabel('Residue position')
        ax4.set_xlim([0, len(sequence) - 1])
        ax4.set_ylim([0, max(predictions['arg']) + 1])
        ax4.grid()

        # RRM
        ax5.plot(x, predictions['RRM'], color='#808080')
        ax5.set_title('RRM', loc='left')
        ax5.set_ylabel('RRM')
        ax5.set_xlabel('Residue position')
        ax5.set_xlim([0, len(sequence) - 1])
        ax5.grid()

        # Disorder
        ax6.plot(x, predictions['disorder'], color='#808080')
        ax6.set_title('Disorder', loc='left')
        ax6.set_ylabel('Disorder')
        ax6.set_xlabel('Residue position')
        ax6.set_xlim([0, len(sequence) - 1])
        ax6.set_ylim([0, max(predictions['disorder']) + 1])
        ax6.grid()

        for i in x:
            ax1.axvspan(i, i+1, facecolor=viterbi_colors[i], alpha=0.5)
            ax2.axvspan(i, i+1, facecolor=viterbi_colors[i], alpha=0.5)
            ax3.axvspan(i, i+1, facecolor=viterbi_colors[i], alpha=0.5)
            ax4.axvspan(i, i+1, facecolor=viterbi_colors[i], alpha=0.5)
            ax5.axvspan(i, i+1, facecolor=viterbi_colors[i], alpha=0.5)
            ax6.axvspan(i, i+1, facecolor=viterbi_colors[i], alpha=0.5)

        axbox = ax1.get_position()

        custom_lines = [
            Line2D([0], [0], color='white', marker='o', markerfacecolor='blue'),
            Line2D([0], [0], color='white', marker='o', markerfacecolor='red'),
            Line2D([0], [0], color='#99ff99', lw=3),
            Line2D([0], [0], color='#bfbfbf', lw=3),
            Line2D([0], [0], color='#ffb499', lw=3),
            Line2D([0], [0], color='#9999ff', lw=3)
        ]

        fig.legend(
            custom_lines,
            ['Tyr', 'Arg', 'Other', 'Spacer', 'putative RRM', 'putative PLD'],
            loc='center',
            ncol=6,
            bbox_to_anchor=[axbox.x0 + 0.5*axbox.width, axbox.y1 + 0.12], bbox_transform=fig.transFigure
        )
        fig.tight_layout()
        plt.savefig(f'${predictions.baseName}_{seq_key}_psp.png')

with open('$predictions', 'r') as json_file:
    prediction_dict = json.loads(json_file.read())

tools = [${dynamine ? '"dynamine",' : ''} ${efoldmine ? '"efoldmine",' : ''} ${disomine ? '"disomine",' : ''} ${agmata ? '"agmata",' : ''} ${psper ? '"psp"' : ''}]

plot_single_sequence(prediction_dict, tools)

if 'psp' in tools:
    plot_psper(prediction_dict)
    """
}

process plotBiophysicalFeaturesOverview {
    tag "${msa.name}"

    input:
    path msa
    val efoldmine
    val disomine

    output:
    path "*.pdf", emit: documents
    path "*.png", emit: plots

    script:
    """
#!/usr/local/bin/python
import json
import os
import math
import numpy as np
import re

os.environ['MPLCONFIGDIR'] = os.getcwd() + "/configs/"
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import matplotlib.ticker as ticker

import Bio
from Bio import AlignIO
from b2bTools.multipleSeq.Predictor import MineSuiteMSA

tools = [${efoldmine ? '"efoldmine",' : ''} ${disomine ? '"disomine"' : ''}]

DynaMine = True
DisoMine = 'disomine' in tools
EFoldMine = 'efoldmine' in tools

msaSuite = MineSuiteMSA()

if DynaMine and not DisoMine and not EFoldMine:
    msaSuite.predictAndMapSeqsFromMSA(f'$msa', predTypes = ('dynamine'))
elif DisoMine and not DynaMine and not EFoldMine:
    msaSuite.predictAndMapSeqsFromMSA(f'$msa', predTypes = ('disoMine'))
elif EFoldMine and not DynaMine and not DisoMine:
    msaSuite.predictAndMapSeqsFromMSA(f'$msa', predTypes = ('eFoldMine'))
elif DynaMine and DisoMine and not EFoldMine:
    msaSuite.predictAndMapSeqsFromMSA(f'$msa', predTypes = ('disoMine', 'dynamine'))
elif DynaMine and EFoldMine and not DisoMine:
    msaSuite.predictAndMapSeqsFromMSA(f'$msa', predTypes = ('eFoldMine','dynamine'))
elif DisoMine and EFoldMine and not DynaMine:
    msaSuite.predictAndMapSeqsFromMSA(f'$msa', predTypes = ('eFoldMine', 'disoMine'))
elif DisoMine and EFoldMine and DynaMine:
    msaSuite.predictAndMapSeqsFromMSA(f'$msa', predTypes = ('eFoldMine', 'disoMine', 'dynamine'))

msaSuite.getDistributions()
jsondata_list = [msaSuite.alignedPredictionDistribs]

if DynaMine and not DisoMine and not EFoldMine:
    NB_SUBPLOTS = 6
    PREDICTION_TITLES = {
        'backbone': "DynaMine backbone dynamics",
        'sidechain': "DynaMine sidechain dynamics",
        'ppII': "DynaMine conformational propensities: ppII (polyproline II)",
        'coil': "DynaMine conformational propensities: Coil",
        'sheet': "DynaMine conformational propensities: Sheet",
        'helix': "DynaMine conformational propensities: Helix",
    }

    PREDICTION_POSITION = {
        'backbone':     0,
        'sidechain':    1,
        'ppII':         2,
        'coil':         3,
        'sheet':        4,
        'helix':        5,
    }

elif DisoMine and not DynaMine and not EFoldMine:
    NB_SUBPLOTS = 1
    PREDICTION_TITLES = {
        'disoMine': "Disorder (disoMine)"
    }

    PREDICTION_POSITION = {
        'disoMine': 0
    }

elif EFoldMine and not DisoMine and not DynaMine:
    NB_SUBPLOTS = 1

    PREDICTION_TITLES = {
        'earlyFolding': "Early folding (EFoldMine)"
    }

    PREDICTION_POSITION = {
        'earlyFolding': 0
    }

elif DynaMine and DisoMine and not EFoldMine:
    NB_SUBPLOTS = 7
    PREDICTION_TITLES = {
        'backbone': "DynaMine backbone dynamics",
        'sidechain': "DynaMine sidechain dynamics",
        'ppII': "DynaMine conformational propensities: ppII (polyproline II)",
        'coil': "DynaMine conformational propensities: Coil",
        'sheet': "DynaMine conformational propensities: Sheet",
        'helix': "DynaMine conformational propensities: Helix",
        'disoMine': "Disorder (disoMine)"
    }

    PREDICTION_POSITION = {
        'backbone':     0,
        'sidechain':    1,
        'ppII':         2,
        'coil':         3,
        'sheet':        4,
        'helix':        5,
        'disoMine':     6
    }

elif DynaMine and EFoldMine and not DisoMine:
    NB_SUBPLOTS = 7
    PREDICTION_TITLES = {
        'backbone': "DynaMine backbone dynamics",
        'sidechain': "DynaMine sidechain dynamics",
        'ppII': "DynaMine conformational propensities: ppII (polyproline II)",
        'coil': "DynaMine conformational propensities: Coil",
        'sheet': "DynaMine conformational propensities: Sheet",
        'helix': "DynaMine conformational propensities: Helix",
        'earlyFolding': "Early folding (EFoldMine)"
    }

    PREDICTION_POSITION = {
        'backbone':     0,
        'sidechain':    1,
        'ppII':         2,
        'coil':         3,
        'sheet':        4,
        'helix':        5,
        'earlyFolding': 6
    }
elif EFoldMine and DisoMine and not DynaMine:
    NB_SUBPLOTS = 2
    PREDICTION_TITLES = {
        'earlyFolding': "Early folding (EFoldMine)",
        'disoMine': "Disorder (disoMine)"
    }

    PREDICTION_POSITION = {
        'earlyFolding': 0,
        'disoMine':     1
    }
elif DynaMine and DisoMine and EFoldMine:
    NB_SUBPLOTS = 8
    PREDICTION_TITLES = {
        'backbone': "DynaMine backbone dynamics",
        'sidechain': "DynaMine sidechain dynamics",
        'ppII': "DynaMine conformational propensities: ppII (polyproline II)",
        'coil': "DynaMine conformational propensities: Coil",
        'sheet': "DynaMine conformational propensities: Sheet",
        'helix': "DynaMine conformational propensities: Helix",
        'earlyFolding': "Early folding (EFoldMine)",
        'disoMine': "Disorder (disoMine)"
    }

    PREDICTION_POSITION = {
        'backbone':     0,
        'sidechain':    1,
        'ppII':         2,
        'coil':         3,
        'sheet':        4,
        'helix':        5,
        'earlyFolding': 6,
        'disoMine':     7
    }

AXIS_TITLES = {
    "x": "Residue position in the MSA",
    "y": "Prediction values"
}

def plot_biophysical_msa(jsondata_list_interest, jsondata_list_selected, sequences, freq_gap, selected_prot):
    colors = ['blue', 'orange']
    residues_count = len(jsondata_list_interest[0]['backbone']['median'])
    sequences_count = len(sequences)

    #Plot representation
    fig, axs = plt.subplots(NB_SUBPLOTS)
    fig.set_figwidth(20)
    fig.set_figheight(50)

    plt.suptitle(f'Predicted biophysical properties of the MSA: {residues_count} aligned residues from {sequences_count} sequences', fontsize=14)

    # These for loops got too complicated, I have to think
    # something simpler to handle the None values in the data
    predictions = jsondata_list_interest[0].keys()
    for prediction_index, biophys_data in enumerate(predictions):
        if biophys_data == 'agmata':
            continue

        subplot_index_row = PREDICTION_POSITION[biophys_data]

        ax = axs[subplot_index_row]
        for data, col in zip(jsondata_list_interest, colors):
            none_idx = []

            for n in range(residues_count):
                if data[biophys_data]['median'][n] == None \
                        or data[biophys_data][
                    'firstQuartile'][n] == None \
                        or data[biophys_data][
                    'thirdQuartile'][n] == None:
                    none_idx.append(n)

            range_list = []
            for n in range(len(none_idx)):
                try:
                    if none_idx[n] + 1 != none_idx[n + 1]:
                        range_list.append(
                            (none_idx[n] + 1, none_idx[n + 1]))
                    else:
                        continue
                except:
                    if len(none_idx) == 1:
                        range_list.append((0, none_idx[0]))
                        range_list.append((none_idx[0] + 1, len(
                            data[biophys_data][
                                'median'])))

                    else:
                        range_list.append((0, none_idx[0]))
                        range_list.append((none_idx[-1] + 1, len(
                            data[biophys_data][
                                'median'])))

            # When there are None values in the data
            if range_list:
                for tuple in range_list:
                    x = np.arange(tuple[0], tuple[1], 1)
                    firstq = \
                        data[biophys_data][
                            'firstQuartile'][
                        tuple[0]:tuple[1]]
                    thirdq = \
                        data[biophys_data][
                            'thirdQuartile'][
                        tuple[0]:tuple[1]]
                    bottom = \
                        data[biophys_data][
                            'bottomOutlier'][
                        tuple[0]:tuple[1]]
                    top = \
                        data[biophys_data]['topOutlier'][
                        tuple[0]:tuple[1]]
                    ax.fill_between(
                        x, firstq, thirdq, alpha=0.3, color=col, label='1st-3rd Quartiles')
                    ax.fill_between(
                        x, bottom, top, alpha=0.1, color=col, label='Outliers')

            # When there aren't None values in the data
            else:
                x = np.arange(0, len(
                    data[biophys_data]['median']), 1)
                firstq = data[biophys_data][
                    'firstQuartile']
                thirdq = data[biophys_data][
                    'thirdQuartile']
                bottom = data[biophys_data][
                    'bottomOutlier']
                top = data[biophys_data]['topOutlier']
                ax.fill_between(
                    x, firstq, thirdq, alpha=0.3, color=col, label='1st-3rd Quartiles')
                ax.fill_between(
                    x, bottom, top, alpha=0.1, color=col, label='Outliers')

            ax.plot(data[biophys_data]['median'], linewidth=1.25, color=col, label='Median')

            #Add the selected protein if there are some
            colors_2 = ['magenta', 'blue', 'cyan'] #adapt in function of number of proteins of interest, now max 3 can be studied simultanously

            if len(selected_prot)>0:
                for count, prot in enumerate(selected_prot):
                    ax.plot(jsondata_list_selected[count][0][biophys_data], '-x', linewidth=1.5, color=colors_2[count], label=f'Prediction {prot}')

            #To represent the amount of sequences at every column of the MSA
            limits = [0.25,0.5,0.75,1]
            limit_previous = 0
            size = 20
            for limit in limits:
                plot_x = []
                plot_y = []
                for idx, freq in enumerate(freq_gap):
                    if limit_previous < freq <= limit:
                        plot_x.append(idx)
                        plot_y.append(data[biophys_data]['median'][idx])
                ax.scatter(plot_x,plot_y,marker='o',s=size, color=col)
                limit_previous = limit
                size += 15

        ax.set_title(PREDICTION_TITLES[biophys_data])

        ax.axis([0, residues_count-1, min(bottom)-0.05, max(top)+0.05])

        ax.set_ylabel(AXIS_TITLES['y'])
        ax.set_xlabel(AXIS_TITLES['x'])

        if biophys_data == 'backbone':
            ax.axhline(y=1.0, color='green', linewidth= 1.5, linestyle='-.', label='Above: Membrane spaning') #Membrane spaning
            ax.axhline(y=0.8, color='orange', linewidth= 1.5, linestyle='-.', label='Above: Rigid') #Membrane spaning
            if min(bottom)-0.05 < 0.69:
                ax.axhline(y=0.69, color='red', linewidth= 1.5, linestyle='-.', label='Above: Context dependent \\nBelow: Flexible') #context dependent (either rigide or flexible)
        if biophys_data == 'earlyFolding':
            ax.axhline(y=0.169, color='red', linewidth= 1.5, linestyle='-.', label='Above: Likely to start folding') #above: likely start protein folding process
        if biophys_data == 'disoMine':
            ax.axhline(y=0.5, color='red', linewidth= 1.5, linestyle='-.', label='Above: Likely to be disordered') #above: likely disordered
        ax.legend(ncol=1, bbox_to_anchor =(1.01,0.5), loc='center left')

    plt.tight_layout()
    fig.subplots_adjust(top=0.96, hspace = 0.2)

    plot_name = re.sub(r'[^\\w.-]+', '_', selected_prot[0].replace(' ', '_'))

    plt.savefig('${msa.baseName}' + '_' + plot_name + '_msa_biophysical_conservation.pdf')
    plt.savefig('${msa.baseName}' + '_' + plot_name + '_msa_biophysical_conservation.png')

    return fig, axs

alignment_file = AlignIO.read('$msa', 'fasta')
counter = [0] * alignment_file.get_alignment_length()
total_seq = len(alignment_file)
for prot in alignment_file:
    selected_prot_seq=list(prot.seq)
    for position,residue in enumerate(selected_prot_seq):
        if residue != "-":
            counter[position] +=1
freq_gap = [i/total_seq for i in counter]

for prot in alignment_file:
    jsondata_list = [msaSuite.alignedPredictionDistribs]
    jsondata_list_selected = []

    predictions_single_seq = msaSuite.allAlignedPredictions
    jsondata_list_selected.append([predictions_single_seq[prot.id]])

    sequences =  msaSuite.seqs
    plot_biophysical_msa(jsondata_list, jsondata_list_selected, sequences, freq_gap, [prot.id])
    """
}


process plotPhylogeneticTree {
    tag "${tree.name}"

    input:
    path tree

    output:
    path "${tree}.svg", emit: treePlot

    script:
    """
    xvfb-run ete3 view --image ${tree}.svg -t ${tree}
    """
}
