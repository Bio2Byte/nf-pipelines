process predictBiophysicalFeatures {
    tag "${sequences.name}"

    input:
    path sequences

    output:
    path 'b2b_msa_results_*.json', emit: predictions
    path 'b2b_msa_stats_*.json', emit: stats

    // Dynamine runs always; psp has no msa function

    script:
    """
    #!/usr/local/bin/python

    import json
    from b2bTools.multipleSeq.Predictor import MineSuiteMSA

    msaSuite = MineSuiteMSA()
    msaSuite.predictAndMapSeqsFromMSA('$sequences', predTypes = ("dynamine", ${params.efoldmine ? '"eFoldMine",' : ''} ${params.disomine ? '"disoMine",' : ''} ${params.agmata ? '"agmata",' : ''}))
    
    predictions_single_seq = msaSuite.allAlignedPredictions
    json.dump(predictions_single_seq, open('b2b_msa_results_${sequences.baseName}.json', 'w'), indent=2)

    predictions=msaSuite.getDistributions()
    json.dump(predictions, open('b2b_msa_stats_${sequences.baseName}.json', 'w'), indent=2)
    """

    
}

process buildMultipleSequenceAlignmentNuc {

    label  'bigboy'
    tag "${sequences.name}"

    input:
    path sequences

    output:
    path '*.anuc', emit: msaNuc
    path "*.aaa", emit: multipleSequenceAlignment

    script:
    """
    macse -prog alignSequences -seq ${sequences} -out_NT ${sequences.baseName}.anuc -out_AA ${sequences.baseName}.aaa
    """
}

process buildMultipleSequenceAlignmentAA {

    label  'bigboy'
    tag "${sequences.name}"

    input:
    path sequences

    output:
    path "*.msa", emit: multipleSequenceAlignment

    script:
    """
    clustalo -i $sequences -o ${sequences}.msa
    """
}

process takeMultipleSequenceAlignment {
    tag "${sequences.name}"

    input:
    path sequences
    val buildTreeEvo 
    val qc

    output:
    path "${sequences.baseName}_checked*" , emit: multipleSequenceAlignment

    script:
    """
    python3 $projectDir/bin/MsaChecker.py $sequences $buildTreeEvo $qc
    """

}

process takeMultipleSequenceAlignmentNuc {
    tag "${sequences.name}"

    input:
    path sequences
    val buildTreeEvo 
    val qc

    output:
    path "${sequences.baseName}_checked*" , emit: multipleSequenceAlignmentNuc

    script:
    """
    python3 $projectDir/bin/MsaChecker.py $sequences $buildTreeEvo $qc
    """

}

process buildPhylogeneticTree {

    label  'bigboy'
    tag "${multipleSequenceAlignment.name}"

    input:
    path multipleSequenceAlignment

    output:
    path "*.treefile", emit: tree

    script:
    """
    iqtree -s $multipleSequenceAlignment -T AUTO -ntmax 4
    """
    //    FastTree $multipleSequenceAlignment > ${multipleSequenceAlignment}.tree

}

process buildPhylogeneticTreeEvol {

    label  'bigboy'
    tag "${multipleSequenceAlignmentNuc.name}"

    input:
    path multipleSequenceAlignmentNuc

    output:
    path "*.treefile", emit: tree
    path "${multipleSequenceAlignmentNuc}.*" , emit:iqtreefiles

    script:
    """
    iqtree -s $multipleSequenceAlignmentNuc -m ECMK07+F+R4 -B 1000 --seqtype CODON11 -nt AUTO --ancestral --rate 
    """

}

process buildLogo {
    tag "${multipleSequenceAlignment.name}"

    input:
    path multipleSequenceAlignment

    output:
    path "*_logo.png", emit: logo

    script:
    """
    weblogo --sequence-type protein --title "MSA logo" --first-index 0 --size large --format png_print < $multipleSequenceAlignment > ${multipleSequenceAlignment.simpleName}_logo.png
    """
}
