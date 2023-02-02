#!/usr/bin/env nextflow

//USAGE: nextflow run fromMSA.nf  --targetSequences ../example.fasta -profile standard,withdocker --agmata--efoldmine --disomine --alignSingleSequences

params.executionTimestamp = new java.text.SimpleDateFormat("yyyy_MM_dd_HH_mm_ss").format(new Date())

// Available predictors:
params.efoldmine       = true
params.disomine        = true
params.agmata          = true
params.fetchStructures = true

// MSA parameters:
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

DynaMine         : ALWAYS
DisoMine         : $params.disomine
EFoldMine        : $params.efoldmine
AgMata           : $params.agmata
PSP              : NOT IMPLEMENTED
Fetch structures : $params.fetchStructures
================================
        SINGLE SEQ or MSA

MSA mode         : YES
Align for MSA    : $params.alignSingleSequences
================================
"""

// Processing input
allSequences = Channel.fromPath(params.targetSequences)

sequencesFiltered = allSequences
    .splitFasta( record: [id: true, seqString: true, sequence: true ])
    .filter { record -> record.seqString.size() >= 5 && record.seqString.size() < 2000 }
    .map {  record -> [id: record.id.replaceAll("[^a-zA-Z0-9]","_"), seqString: record.seqString, sequence: record.sequence] }

sequencesRemoved = allSequences
    .splitFasta( record: [id: true, seqString: true, sequence: true ])
    .filter { record -> record.seqString.size() < 5 || record.seqString.size() >= 2000 }

sequencesRemoved.view { it.id } // TODO: Put the removed sequences into another channel and inform user about the excluded.

sequencesGrouped = sequencesFiltered
    .collectFile(name: "${targetSequencesFile.baseName}.fasta", newLine: true) {
        item -> '>' + item.id + '\n' + item.sequence + '\n'
    }
    .splitFasta( file: true, by: params.groupBy )

// Modules
include {
    predictBiophysicalFeatures;
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
workflow {
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

    predictBiophysicalFeatures(params.resultsDirectory, multipleSequenceAlignment)
    // plotBiophysicalFeatures(params.resultsDirectory, predictBiophysicalFeatures.out.predictions)
    // plotAgmata(params.resultsDirectory, predictBiophysicalFeatures.out.predictions)

    fetchStructure(params.resultsDirectory, sequencesFiltered.map { record -> [id: record.id, seqString: record.seqString.take(400)] })

    // compressPredictions(
    //     workflowSingleSequences.out.predictions.view(),
    //     workflowSingleSequences.out.index.view(),
    //     fetchStructure.out.esmStructures.view()
    // )
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Time to complete workflow execution: $workflow.duration"
    println "Execution status: ${workflow.success ? 'Success' : 'Failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
