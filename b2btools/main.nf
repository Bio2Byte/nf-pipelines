#!/usr/bin/env nextflow
// USAGE: nextflow run main.nf -resume -with-dag pipeline.png

params.targetSequences = "$launchDir/input_small.fasta"
params.groupBy = 100

targetSequencesFile = file(params.targetSequences)
allSequences = Channel.fromPath(params.targetSequences)

sequencesFiltered = allSequences
    .splitFasta( record: [id: true, seqString: true, sequence: true ])
    .filter { record -> record.seqString.size() >= 5 && record.seqString.size() < 2000 }

sequencesGrouped = sequencesFiltered
    .collectFile(name: "${targetSequencesFile.baseName}.fasta", newLine: true) {
        item -> '>' + item.id + '\n' + item.sequence + '\n'
    }
    .splitFasta( file: true, by: params.groupBy )

process predictBiophysicalFeatures {
    tag "${sequences.baseName}"

    input:
    path sequences

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
    single_seq.predict(tools=['dynamine', 'efoldmine', 'disomine', 'agmata'])
    all_predictions = single_seq.get_all_predictions()
    json.dump(all_predictions, open('b2b_results_${sequences.baseName}.json', 'w'), indent=2)
    with open('b2b_results_${sequences.baseName}.index', 'w') as index_file:
        index_file.write("id,json_file,residues_count,avg_backbone,avg_coil,avg_disoMine,avg_earlyFolding,avg_helix,avg_ppII,avg_sheet,avg_sidechain\\n")
        for sequence_key in all_predictions.keys():
            prediction = all_predictions[sequence_key]
            seq_len = len(prediction['seq'])
            avg_backbone = average(prediction['backbone'])
            avg_coil = average(prediction['coil'])
            avg_disoMine = average(prediction['disoMine'])
            avg_earlyFolding = average(prediction['earlyFolding'])
            avg_helix = average(prediction['helix'])
            avg_ppII = average(prediction['ppII'])
            avg_sheet = average(prediction['sheet'])
            avg_sidechain = average(prediction['sidechain'])
            index_line = "{0},b2b_results_${sequences.baseName}.json,{1},{2:.3f},{3:.3f},{4:.3f},{5:.3f},{6:.3f},{7:.3f},{8:.3f},{9:.3f}\\n".format(
                sequence_key,
                seq_len,
                avg_backbone,
                avg_coil,
                avg_disoMine,
                avg_earlyFolding,
                avg_helix,
                avg_ppII,
                avg_sheet,
                avg_sidechain
            )
            index_file.write(index_line)
    """
}

process compressPredictions {
    publishDir "results", mode: 'copy'

    input:
    path predictions
    path indexFile

    output:
    path "*.tar.gz"

    script:
    """
    tar -czvhf b2b_results.tar.gz $predictions $indexFile
    """
}

workflow b2bToolsAnalysis {
    take:
    sequencesGrouped

    main:
    predictBiophysicalFeatures(sequencesGrouped)

    emit:
    predictions = predictBiophysicalFeatures.out.predictions
    indexes = predictBiophysicalFeatures.out.index
}

workflow {
    b2bToolsAnalysis(sequencesGrouped)

    // Main workflow
    compressPredictions(
        b2bToolsAnalysis.out.predictions.collect(),
        b2bToolsAnalysis.out.indexes.collectFile(name: "${targetSequencesFile.baseName}.index", keepHeader: true)
    )
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
