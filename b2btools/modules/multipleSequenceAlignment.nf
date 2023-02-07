process predictBiophysicalFeatures {
    tag "${sequences.name}"

    input:
    path sequences

    output:
    path '*.json', emit: predictions

    // Dynamine runs always; psp has no msa function

    script:
    """
    #!/usr/local/bin/python

    import json
    from b2bTools.multipleSeq.Predictor import MineSuiteMSA


    with open('$sequences', 'r') as file:
        print(file.read())

    msaSuite = MineSuiteMSA()
    msaSuite.predictAndMapSeqsFromMSA('$sequences', predTypes = (${params.dynamine ? '"dynamine",' : ''} ${params.efoldmine ? '"efoldmine",' : ''} ${params.disomine ? '"disomine",' : ''}))

    predictions=msaSuite.getDistributions()
    #jsondata_list = [msaSuite.alignedPredictionDistribs]  

    json.dump(predictions, open('b2b_msa_results_${sequences.baseName}.json', 'w'), indent=2)
    """

    
}

process buildMultipleSequenceAlignment {
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

    output:
    path "*.msa", emit: multipleSequenceAlignment

    script:
    """
    cp $sequences ${sequences}.msa
    """

}

process buildPhylogeneticTree {
    tag "${multipleSequenceAlignment.name}"

    input:
    path multipleSequenceAlignment

    output:
    path "*.tree", emit: tree

    script:
    """
    FastTree $multipleSequenceAlignment > ${multipleSequenceAlignment}.tree
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
    weblogo --sequence-type protein --title "MSA logo" --size large --format png_print < $multipleSequenceAlignment > ${multipleSequenceAlignment.simpleName}_logo.png
    """
}
