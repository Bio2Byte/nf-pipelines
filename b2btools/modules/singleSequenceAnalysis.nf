process predictBiophysicalFeatures {
    publishDir "$resultsDirectoryPath", mode: 'copy'
    tag "${sequences.baseName}"

    input:
    path resultsDirectoryPath
    path sequences
    val dynamine
    val efoldmine
    val disomine
    val agmata
    val psper

    output:
    path '*.json', emit: predictions
    path '*.index', emit: index

    script:
    """
    #!/usr/local/bin/python
    from b2bTools import SingleSeq
    import json

    def average(lst):
        return sum(lst) / len(lst)
    single_seq = SingleSeq("$sequences")

    tool_list = [${dynamine ? '"dynamine",' : ''} ${efoldmine ? '"efoldmine",' : ''} ${disomine ? '"disomine",' : ''} ${agmata ? '"agmata",' : ''} ${psper ? '"psp"' : ''}]
    tool_list=[x for x in tool_list if x]

    single_seq.predict(tools=tool_list)
    all_predictions = single_seq.get_all_predictions()
    json.dump(all_predictions, open('b2b_results_${sequences.baseName}.json', 'w'), indent=2)

    with open('b2b_results_${sequences.baseName}.index', 'w') as index_file:
        index_file.write("id,json_file,residues_count,avg_backbone,avg_coil,avg_helix,avg_ppII,avg_sheet,avg_sidechain,avg_earlyFolding \\n")#,avg_disoMine
        for sequence_key in all_predictions.keys():
            prediction = all_predictions[sequence_key]
            seq_len = len(prediction['seq'])
            avg_backbone = average(prediction['backbone'])
            avg_coil = average(prediction['coil'])
            #  avg_disoMine = average(prediction['disoMine'])
            avg_earlyFolding = average(prediction['earlyFolding'])
            avg_helix = average(prediction['helix'])
            avg_ppII = average(prediction['ppII'])
            avg_sheet = average(prediction['sheet'])
            avg_sidechain = average(prediction['sidechain'])
            index_line = "{0},b2b_results_${sequences.baseName}.json,{1},{2:.3f},{3:.3f},{4:.3f},{5:.3f},{6:.3f},{7:.3f},{8:.3f}\\n".format( #,{9:.3f}
                sequence_key,
                seq_len,
                avg_backbone,
                avg_coil,
                #avg_disoMine,
                avg_earlyFolding,
                avg_helix,
                avg_ppII,
                avg_sheet,
                avg_sidechain
            )
            index_file.write(index_line)
    """
}