#!/usr/bin/env nextflow

//USAGE: nextflow run main.nf  --targetSequences ../../example.fasta -profile standard,withdocker --dynamine --efoldmine --disomine --agmata --alignSingleSequences --msa

// watch index file! if disomine, efmine turned off it still wants to plot them, causes error

params.executionTimestamp = new java.text.SimpleDateFormat("yyyy_MM_dd_HH_mm_ss").format(new Date())

// Available predictors:
params.dynamine  = false
params.efoldmine = false
params.disomine  = false
params.agmata    = false
params.psper     = false

// MSA parameters:
// maybe flag for use all flags?
params.fetchEsm      = false
params.msa           = false
params.alignSingleSequences = false

// Input file
params.targetSequences      = "$launchDir/example.fasta"
targetSequencesFile         = file(params.targetSequences)
params.resultsDirectory     = "$projectDir/${params.executionTimestamp}"
params.compressedFile       = "$projectDir/${targetSequencesFile.simpleName}_${params.executionTimestamp}.tar.gz"
params.groupBy              = 10

log.info """\

================================
        LIST OF PARAMETERS
================================
            GENERAL

Launch dir       : $launchDir
Prject dir       : $projectDir
Execution time   : $params.executionTimestamp
Results dir      : $params.resultsDirectory
Compressed file  : $params.compressedFile
================================
            INPUT FILES

Input-file       : $params.targetSequences
JSON page size   : $params.groupBy
================================
            PREDICTORS

DynaMine         : $params.dynamine
DisoMine         : $params.disomine
EFoldMine        : $params.efoldmine
AgMata           : $params.agmata
PSP              : $params.psper
================================
        SINGLE SEQ or MSA

MSA mode         : $params.msa
Align for MSA    : $params.alignSingleSequences
Fetch structures : $params.fetchEsm
================================
"""

// Processing input
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

// Modules
include {
    predictBiophysicalFeatures as predictBiophysicalFeaturesForMSA;
    buildMultipleSequenceAlignment;
    takeMultipleSequenceAlignment;
    buildPhylogeneticTree;
    buildLogo;
} from "${launchDir}/modules/multipleSequenceAlignment"

include {
    fetchStructure;
} from "${launchDir}/modules/structures"

include {
    plotBiophysicalFeatures;
    plotAgmata;
    plotPhylogeneticTree;
} from "${launchDir}/modules/plots"

include {
    predictBiophysicalFeatures as predictBiophysicalFeaturesForSingleSeq;
} from "${launchDir}/modules/singleSequenceAnalysis"

// Processes
process compressPredictions {
    publishDir "results", mode: 'copy'

    input:
    // val outputs

    output:
    path "*.tar.gz"

    script:
    """
    echo Creating compressed file: ${params.compressedFile}
    echo Content of compressed: ${params.resultsDirectory}
    ls ${params.resultsDirectory}

    tar -czvf ${params.compressedFile} ${params.resultsDirectory}
    """
}

// Workflows
workflow workflowSingleSequences {
    take:
        sequencesGrouped

    main:
        predictBiophysicalFeaturesForSingleSeq(
            params.resultsDirectory,
            sequencesGrouped,
            params.dynamine,
            params.efoldmine,
            params.disomine,
            params.agmata,
            params.psper
        )

    emit:
        predictions = predictBiophysicalFeaturesForSingleSeq.out.predictions
        indexes     = predictBiophysicalFeaturesForSingleSeq.out.index
}

workflow workflowMSA {
    take:
        allSequences

    main:
        if (params.alignSingleSequences) {
            buildMultipleSequenceAlignment(params.resultsDirectory, allSequences)
            multipleSequenceAlignment = buildMultipleSequenceAlignment.out.multipleSequenceAlignment
        }
        else {
            takeMultipleSequenceAlignment(allSequences)
            multipleSequenceAlignment = takeMultipleSequenceAlignment.out.multipleSequenceAlignment
        }

        buildPhylogeneticTree(params.resultsDirectory, multipleSequenceAlignment)
        plotPhylogeneticTree(params.resultsDirectory, buildPhylogeneticTree.out.tree)
        buildLogo(params.resultsDirectory, multipleSequenceAlignment)

    emit:
        multipleSequenceAlignment
        tree = buildPhylogeneticTree.out.tree
        treePlot = plotPhylogeneticTree.out.treePlot
        logo = buildLogo.out.logo
}

workflow workflowMSABiophysicalFeatures {
    take:
        multipleSequenceAlignment

    main:
        predictBiophysicalFeaturesForMSA(params.resultsDirectory, multipleSequenceAlignment)
        plotBiophysicalFeatures(params.resultsDirectory, predictBiophysicalFeaturesForMSA.out.predictions)
        plotAgmata(params.resultsDirectory, predictBiophysicalFeaturesForMSA.out.predictions)

    emit:
        predictions = predictBiophysicalFeaturesForMSA.out.predictions
        plots = plotBiophysicalFeatures.out.plots
        agmata_plots = plotAgmata.out.plots
}



workflow workflowPredictions{
    take:
        allSequences
        sequencesFiltered
        sequencesGrouped

    main:
        if (params.msa) {
            log.info "Running specific workflow steps for Multiple Sequence Alignment"

            workflowMSA(allSequences)
            b2bPredictions = workflowMSABiophysicalFeatures(workflowMSA.out.multipleSequenceAlignment)
        }
        else {
            log.info "Running specific workflow steps for Single Sequences"

            b2bPredictions = workflowSingleSequences(sequencesGrouped)
        }

        log.info "Running generic workflow steps"

        fetchStructure(params.resultsDirectory, sequencesFiltered.map { record -> [id: record.id, seqString: record.seqString.take(400)] })
    emit:
        esmStructures = fetchStructure.out.esmStructures
        b2bBiophysicalPredictions = b2bPredictions.out.predictions
}

workflow {
    workflowPredictions(allSequences, sequencesFiltered, sequencesGrouped)
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Time to complete workflow execution: $workflow.duration"
    println "Execution status: ${workflow.success ? 'Success' : 'Failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
