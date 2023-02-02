process predictBiophysicalFeatures {
    publishDir "${resultsDirectoryPath}", mode: 'copy'

    tag "${sequences.baseName}"

    input:
    path resultsDirectoryPath
    path sequences

    output:
    path '*.json', emit: predictions

    // Dynamine runs always; psp has no msa function

    script:
    """
    #!/usr/local/bin/python

    import matplotlib.pyplot as plt
    from b2bTools import MultipleSeq
    import json

    tool_list = [${params.efoldmine ? '"efoldmine,"' : ''} ${params.disomine ? '"disomine,"' : ''} ${params.agmata ? '"agmata"' : ''}]
    tool_list=[x for x in tool_list if x]

    msaSeq = MultipleSeq()
    msaSeq.from_aligned_file('$sequences',tools=tool_list)

    predictions = msaSeq.get_all_predictions_msa()
    json.dump(predictions, open('b2b_msa_results_${sequences.baseName}.json', 'w'), indent=2)
    """
}

process buildMultipleSequenceAlignment {
    publishDir "$resultsDirectoryPath", mode: 'copy'

    input:
    path resultsDirectoryPath
    path sequences

    output:
    path "*.msa", emit: multipleSequenceAlignment

    when:
    params.alignSingleSequences == true

    script:
    """
    clustalo -i $sequences -o ${sequences}.msa
    """
}

process takeMultipleSequenceAlignment {
    input:
    path resultsDirectoryPath
    path sequences

    output:
    path "*.msa", emit: multipleSequenceAlignment

    when:
    params.alignSingleSequences == false

    script:
    """
    cp $sequences ${sequences}.msa
    """

}

process buildPhylogeneticTree {
    publishDir "$resultsDirectoryPath", mode: 'copy'

    input:
    path resultsDirectoryPath
    path multipleSequenceAlignment

    output:
    path "*.tree", emit: tree

    script:
    """
    FastTree $multipleSequenceAlignment > ${multipleSequenceAlignment}.tree
    """
}

process buildLogo {
    publishDir "$resultsDirectoryPath", mode: 'copy'

    input:
    path resultsDirectoryPath
    path multipleSequenceAlignment

    output:
    path "*_logo.png", emit: logo

    script:
    """
    weblogo --sequence-type protein --title "MSA logo" --size large --format png_print < $multipleSequenceAlignment > ${multipleSequenceAlignment.simpleName}_logo.png
    """
}
