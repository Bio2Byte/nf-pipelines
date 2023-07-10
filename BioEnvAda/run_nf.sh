#!/bin/bash
now=`date +"%s"`
data=example 
~/nextflow run pipeline.nf \
    -resume \
    -profile standard,withdocker \
    --targetSequences sequences/$data.fasta  \
    --alignSequences \
    --buildTreeEvo \
    --efoldmine \
    --disomine \
    --plotBiophysicalFeatures \
    --fetchStructures \
    --csubst \
    --branchIds '1,5' \
    --eteEvol 'M8' \
    --outGroup 'Cya_NS01' \
    --selectedProteins  'AncNode14,Syn_BIOS_U3' \
    --buildLogo \
    --plotTree \
    >>  $data-$now.nflog
sessionName=$(sed -n '2s/.*\[\(.*\)\].*/\1/p' $data-$now.nflog)
~/nextflow log | grep $sessionName >> $data-$now.nflog


#-resume \
#-profile standard, withdocker, withsingularity, withconda \
#General Data Preparation
#    --type aa \
#    --clustering 1 \
#    --relabel \
#    --alignSequences \
#    --buildTreeEvo \
#B2Btools selection
#    --efoldmine \
#    --disomine \
#    --agmata \
#Evol predictions
#    --csubst \
#    --branchIds '1,2,3' \
#    --eteEvol\
#    --outGroup 'Cya_NS01_5_2B_1' \
#Generate Plots
#    --plotBiophysicalFeatures \
#    --selectedProteins  'Syn_WH5701 ,Syn_RS9917' \
#    --buildLogo \
#    --plotTree \
#    --fetchStructures \