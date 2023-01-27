#!/usr/bin/env nextflow
// USAGE: nextflow run main.nf -resume -with-dag pipeline.png

//USAGE: nextflow run main.nf  --targetSequences ../../example.fasta -profile standard,withdocker --dynamine --efoldmine --disomine --align_for_msa --msa


// watch index file! if disomine, efmine turned off it still wants to plot them, causes error

params.dynamine = false
params.efoldmine = false
params.disomine = false

params.agmata = false
params.fetchEsm = false

params.msa = false
params.align_for_msa = false

// sequence file

params.targetSequences = "$launchDir/example.fasta"
params.groupBy =10

targetSequencesFile = file(params.targetSequences)
allSequences = Channel.fromPath(params.targetSequences)

sequencesFiltered = allSequences
    .splitFasta( record: [id: true, seqString: true, sequence: true ])
    .filter { record -> record.seqString.size() >= 5 && record.seqString.size() < 2000 }
    .map {  record -> [id: record.id.replaceAll("[^a-zA-Z0-9]","_"), seqString: record.seqString, sequence: record.sequence] }

// TODO: Put the removed sequences into another channel and inform user about the excluded.

sequencesGrouped = sequencesFiltered
    .collectFile(name: "${targetSequencesFile.baseName}.fasta", newLine: true) {
        item -> '>' + item.id + '\n' + item.sequence + '\n'
    }
    .splitFasta( file: true, by: params.groupBy )


process createMultipleSequenceAlignment {
    input:
    path sequences

    output:
    path "*.msa", emit: multipleSequenceAlignment

    when:
    params.align_for_msa == true

    script:
    """
    clustalo -i $sequences -o ${sequences}.msa
    """
}

process takeMultipleSequenceAlignment {

    publishDir "results", mode: 'copy'
    input:
    path sequences

    output:
    path "*.msa", emit: multipleSequenceAlignment

    when:
    params.align_for_msa == false

    script:
    """
    cp $sequences ${sequences}.msa
    """

}

process buildPhylogeneticTree {
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

    publishDir "results", mode: 'copy'
    input:
    path multipleSequenceAlignment

    output:
    path "*_logo.png", emit: logo

    script:
    """
    weblogo --sequence-type protein --title "MSA logo" --size large --format png_print < $multipleSequenceAlignment > ${multipleSequenceAlignment.simpleName}_logo.png
    """
}

//render tree seems to take longest from all processes?
process renderTree {
    publishDir "results", mode: 'copy'
    input:
    path tree

    output:
    path "tree.png", emit: treePlot

    script:
    """
    #!/usr/bin/env Rscript
    tempLibPath <- tempdir()
    dir.create(tempLibPath) # create personal library
    .libPaths(tempLibPath)  # add to the path
    install.packages("ape")  # install like always
    library(ape)  # use library like always

    mytr <- read.tree("$tree")
    png("tree.png")
    plot(mytr)
    dev.off()
    """
}

process predictBiophysicalFeatures {
    tag "${sequences.baseName}"

    input:
    path sequences

    output:
    path '*.json', emit: predictions

    script:
    """
    python -m b2bTools  ${params.dynamine ? '-dynamine' : ''} ${params.efoldmine ? '-efoldmine' : ''} ${params.disomine ? '-disomine' : ''} ${params.agmata ? '-agmata' : ''} -file $sequences -output ${sequences}.json -identifier test
    """
}

process plotBiophysicalFeatures {
    publishDir "results", mode: 'copy'
    tag "${predictions.baseName}"

    input:
    path predictions

    output:
    path "*_predictions.png", emit: plots

    script:
    """
    #!/usr/bin/python3
    import matplotlib.pyplot as plt
    import json
    with open('$predictions', 'r') as json_file:
        prediction_dict = json.loads(json_file.read())
    for id, prediction in enumerate(prediction_dict['results']):
        fig, axs = plt.subplots(2, 4)
        ax1 = axs[0, 0]
        ax2 = axs[0, 1]
        ax3 = axs[0, 2]
        ax4 = axs[0, 3]
        ax5 = axs[1, 0]
        ax6 = axs[1, 1]
        ax7 = axs[1, 2]
        ax8 = axs[1, 3]
        fig.set_figwidth(30)
        fig.set_figheight(10)
        fig.suptitle(f"Single Sequence Predictions for: {prediction['proteinID']}")
        x_position = range(len(prediction['sequence']))
        backbone_pred = prediction['backbone']
        coil_pred = prediction['coil']
        sheet_pred = prediction['sheet']
        ppII_pred = prediction['ppII']
        helix_pred = prediction['helix']
        sidechain_pred = prediction['sidechain']
        #disomine_pred = prediction['disoMine']
        earlyFolding_pred = prediction['earlyFolding']
        ax1.plot(x_position, backbone_pred, label="Backbone")
        ax2.plot(x_position, sidechain_pred, label="Side chain")
        ax3.plot(x_position, coil_pred, label="Coil")
        ax4.plot(x_position, sheet_pred, label="Sheet")
        ax5.plot(x_position, ppII_pred, label="ppII")
        ax6.plot(x_position, helix_pred, label="Helix")
        #ax7.plot(x_position, disomine_pred, label="Disorder")
        ax8.plot(x_position, earlyFolding_pred, label="Early folding")
        ax1.set_title('DynaMine backbone dynamics')
        ax1.set_ylim([-0.2, 1.2])
        ax1.set_xlabel('residue index')
        ax1.set_ylabel('prediction values')
        ax1.axhspan(1, 1.2, alpha=0.3, color='red')
        ax1.axhspan(0.8, 1, alpha=0.5, color='pink')
        ax1.axhspan(0.69, 0.8, alpha=0.5, color='orange')
        ax1.axhspan(-0.2, 0.69, alpha=0.5, color='yellow')
        ax1.grid(axis='y')
        ax1.set_xlim([0, len(prediction['sequence']) - 1])
        ax2.set_title('DynaMine sidechain dynamics')
        ax2.set_ylim([-0.2, 1.2])
        ax2.set_xlabel('residue index')
        ax2.set_ylabel('prediction values')
        ax2.grid(axis='y')
        ax2.set_xlim([0, len(prediction['sequence']) - 1])
        ax3.set_title('DynaMine conformational propensities: Coil')
        ax3.set_ylim([-0.2, 1.2])
        ax3.set_xlabel('residue index')
        ax3.set_ylabel('prediction values')
        ax3.grid(axis='y')
        ax3.set_xlim([0, len(prediction['sequence']) - 1])
        ax4.set_title('DynaMine conformational propensities: Sheet')
        ax4.set_ylim([-0.2, 1.2])
        ax4.set_xlabel('residue index')
        ax4.set_ylabel('prediction values')
        ax4.grid(axis='y')
        ax4.set_xlim([0, len(prediction['sequence']) - 1])
        ax5.set_title('DynaMine conformational propensities: ppII (polyproline II)')
        ax5.set_ylim([-0.2, 1.2])
        ax5.set_xlabel('residue index')
        ax5.set_ylabel('prediction values')
        ax5.grid(axis='y')
        ax5.set_xlim([0, len(prediction['sequence']) - 1])
        ax6.set_title('DynaMine conformational propensities: Helix')
        ax6.set_ylim([-0.2, 1.2])
        ax6.set_xlabel('residue index')
        ax6.set_ylabel('prediction values')
        ax6.grid(axis='y')
        ax6.set_xlim([0, len(prediction['sequence']) - 1])
        ax7.set_title('Early folding (EFoldMine)')
        ax7.set_ylim([-0.2, 1.2])
        ax7.set_xlabel('residue index')
        ax7.set_ylabel('prediction values')
        ax7.axhspan(-0.2, 0.169, alpha=0.5, color='yellow')
        ax7.axhspan(0.169, 1.2, alpha=0.5, color='orange')
        ax7.grid(axis='y')
        ax7.set_xlim([0, len(prediction['sequence']) - 1])
        ax8.set_title('Disorder (disoMine)')
        ax8.set_ylim([-0.2, 1.2])
        ax8.set_xlabel('residue index')
        ax8.set_ylabel('prediction values')
        ax8.axhspan(0.5, 1.2, alpha=0.5, color='orange')
        ax8.axhspan(-0.2, 0.5, alpha=0.5, color='yellow')
        ax8.grid(axis='y')
        ax8.set_xlim([0, len(prediction['sequence']) - 1])
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3), fancybox=True, shadow=True, ncol=8)
        plt.subplots_adjust(hspace=0.4)
        plt.savefig(prediction['proteinID'] + '_predictions.png')
    """
}

process plotAgmata {
    tag "${predictions.baseName}"

    input:
    path predictions

    output:
    path "*_agmata_prediction.png", emit: plots

    when:
    params.agmata == true

    script:
    """
    #!/usr/bin/python3
    import matplotlib.pyplot as plt
    import json
    with open('$predictions', 'r') as json_file:
        prediction_dict = json.loads(json_file.read())
    for id, prediction in enumerate(prediction_dict['results']):
        fig, ax = plt.subplots(1, 1)
        fig.set_figwidth(30)
        fig.set_figheight(5)
        fig.suptitle('Agmata aggregation propensity')
        agmata_pred = prediction['agmata']
        ax.plot(range(len(agmata_pred)), agmata_pred, label="AgMata")
        ax.set_xlim([0, len(agmata_pred) - 1])
        ax.set_xlabel('residue index')
        ax.set_ylabel('prediction values')
        ax.grid(axis='y')
        plt.savefig(prediction['proteinID'] + '_agmata_prediction.png')
    """
}

process fetchStructure {
    publishDir "results", mode: 'copy'
    tag "$id"

    input:
    tuple val(id), val(sequence)

    output:
    path "*.pdb", emit: esmStructures

    when:
    params.fetchEsm == true

    script:
    """
    echo Folding sequence using ESM Atlas
    curl -X POST --data "$sequence" https://api.esmatlas.com/foldSequence/v1/pdb/ > ${id}.pdb
    echo Preview of the PDB content
    head ${id}.pdb
    tail ${id}.pdb
    """
}

process compressPredictions {
    publishDir "results", mode: 'copy'

    input:
    path predictions
    path plots

    path multipleSequenceAlignment
    path tree
    path treePlot
    path logo

    path esmStructures
    path agmata_plots

    output:
    path "*.tar.gz"

    script:
    """
    tar -czvhf ${multipleSequenceAlignment.simpleName}.tar.gz $tree $treePlot $multipleSequenceAlignment $logo $predictions $plots ${params.agmata ? agmata_plots : ''} ${params.fetchEsm ? esmStructures : ''}
    """
}


process sSeqPredictBiophysicalFeatures {
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


    tool_list = ['${params.dynamine ? '-dynamine' : ''}', '${params.efoldmine ? '-efoldmine' : ''}', '${params.disomine ? '-disomine' : ''}', '${params.agmata ? '-agmata' : ''}']
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

process sSeqCompressPredictions {
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

workflow sSeqB2bToolsAnalysis {
    take:
    sequencesGrouped

    main:
    sSeqPredictBiophysicalFeatures(sequencesGrouped)

    emit:
    predictions = sSeqPredictBiophysicalFeatures.out.predictions
    indexes = sSeqPredictBiophysicalFeatures.out.index
}

workflow multipleSequenceAlignmentAnalysis {
    take:
    allSequences

    main:
    if (params.align_for_msa == true){
        createMultipleSequenceAlignment(allSequences)
        multipleSequenceAlignment = createMultipleSequenceAlignment.out.multipleSequenceAlignment
    }
    else {
        takeMultipleSequenceAlignment(allSequences)
        multipleSequenceAlignment = takeMultipleSequenceAlignment.out.multipleSequenceAlignment
    }

    buildPhylogeneticTree(multipleSequenceAlignment)
    buildLogo(multipleSequenceAlignment)
    renderTree(buildPhylogeneticTree.out.tree)

    emit:

    multipleSequenceAlignment

    tree = buildPhylogeneticTree.out.tree
    treePlot = renderTree.out.treePlot
    logo = buildLogo.out.logo

}

workflow b2bToolsAnalysis {
    take:
    sequencesGrouped

    main:
    predictBiophysicalFeatures(sequencesGrouped)

    plotAgmata(predictBiophysicalFeatures.out.predictions)
    plotBiophysicalFeatures(predictBiophysicalFeatures.out.predictions)

    emit:
    predictions = predictBiophysicalFeatures.out.predictions
    plots = plotBiophysicalFeatures.out.plots
    agmata_plots = plotAgmata.out.plots
}

workflow {
    if (params.msa) {
        dummyAgmata = file('dummy')
        dummyEsm = file('dummy2')

        // First sub-workflow
        multipleSequenceAlignmentAnalysis(allSequences)
        // Second sub-workflow
        b2bToolsAnalysis(sequencesGrouped)
        // Third sub-workflow
        fetchStructure(sequencesFiltered.map { record -> [id: record.id, seqString: record.seqString.take(400)] })

        // Main workflow
        compressPredictions(
            b2bToolsAnalysis.out.predictions.collect(),
            b2bToolsAnalysis.out.plots.collect(),
            multipleSequenceAlignmentAnalysis.out.multipleSequenceAlignment,
            multipleSequenceAlignmentAnalysis.out.tree,
            multipleSequenceAlignmentAnalysis.out.treePlot,
            multipleSequenceAlignmentAnalysis.out.logo,

            // if agmata and fetchesm are not calculated, compressing is not happening , no error
            // set optional input files
            b2bToolsAnalysis.out.agmata_plots.collect().ifEmpty(dummyAgmata),
            fetchStructure.out.esmStructures.collect().ifEmpty(dummyEsm),
        )
    }
    else {
        sSeqB2bToolsAnalysis(sequencesGrouped)

        // sseq Main workflow
        sSeqCompressPredictions(
            sSeqB2bToolsAnalysis.out.predictions.collect(),
            sSeqB2bToolsAnalysis.out.indexes.collectFile(name: "${targetSequencesFile.baseName}.index", keepHeader: true)
        )
    }
}


workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
