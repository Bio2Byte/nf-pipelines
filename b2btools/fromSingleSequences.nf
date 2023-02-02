#!/usr/bin/env nextflow

params.executionTimestamp = new java.text.SimpleDateFormat("yyyy_MM_dd_HH_mm_ss").format(new Date())

// Available predictors:
params.dynamine         = true
params.efoldmine        = true
params.disomine         = true
params.agmata           = true
params.psper            = true
params.fetchStructures  = true

// Input file
params.targetSequences      = "$launchDir/example.fasta"
targetSequencesFile         = file(params.targetSequences)
params.resultsDirectory     = "$projectDir/${params.executionTimestamp}"
params.compressedFile       = "$projectDir/${targetSequencesFile.simpleName}_${params.executionTimestamp}.tar.gz"
params.groupBy              = 10

log.info """\

Usage:

\$ nextflow run -resume fromSingleSequences.nf \

    -profile standard,withdocker \

    --targetSequences ../example.fasta \

    --dynamine \

    --efoldmine \

    --disomine \

    --agmata

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
Fetch structures : $params.fetchStructures
================================
        SINGLE SEQ or MSA

MSA mode         : No
Align for MSA    : No
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
    fetchStructure;
} from "${launchDir}/modules/structures"

include {
    plotBiophysicalFeatures;
} from "${launchDir}/modules/plots"

include {
    predictBiophysicalFeatures;
} from "${launchDir}/modules/singleSequenceAnalysis"

// Processes
process compressPredictions {
    input:
    val predictions
    val index
    val plots
    val esmStructures

    script:
    """
    echo Creating compressed file: ${params.compressedFile}
    echo Content of compressed: ${params.resultsDirectory}
    ls ${params.resultsDirectory}

    tar -czvf ${params.compressedFile} -C / ${params.resultsDirectory}
    """
}

// Workflows
workflow {
    predictBiophysicalFeatures(
        params.resultsDirectory,
        sequencesGrouped,
        params.dynamine,
        params.efoldmine,
        params.disomine,
        params.agmata,
        params.psper
    )

    plotBiophysicalFeatures(
        params.resultsDirectory,
        predictBiophysicalFeatures.out.predictions,
        params.dynamine,
        params.efoldmine,
        params.disomine,
        params.agmata,
        params.psper
    )
    fetchStructure(params.resultsDirectory, sequencesFiltered.map { record -> [id: record.id, seqString: record.seqString.take(400)] })

    compressPredictions(
        predictBiophysicalFeatures.out.predictions,
        predictBiophysicalFeatures.out.index,
        plotBiophysicalFeatures.out.plots,
        fetchStructure.out.esmStructures
    )
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Time to complete workflow execution: $workflow.duration"
    println "Execution status: ${workflow.success ? 'Success' : 'Failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
