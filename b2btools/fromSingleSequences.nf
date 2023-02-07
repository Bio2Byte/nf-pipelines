#!/usr/bin/env nextflow

// RUN: nextflow run -resume fromSingleSequences.nf -profile standard,withdocker --targetSequences ../example.fasta --dynamine --efoldmine --disomine --agmata
params.executionTimestamp = new java.text.SimpleDateFormat("yyyy_MM_dd_HH_mm_ss").format(new Date())

// Available predictors:
params.dynamine         = true
params.efoldmine        = false
params.disomine         = false
params.agmata           = false
params.psper            = false
params.fetchStructures  = false

// Input file
params.targetSequences      = "$launchDir/example.fasta"
targetSequencesFile         = file(params.targetSequences)
allSequences                = Channel.fromPath(params.targetSequences)

params.compressedFile       = "${targetSequencesFile.simpleName}_${params.executionTimestamp}"
params.groupBy              = 10

log.info """\

Usage:

\$ nextflow run -resume fromSingleSequences.nf \\
    -with-report report_fromSingleSequences.html \\
    -with-dag flowchart_fromSingleSequences.png \\
    -profile standard,withdocker \\
    --targetSequences ../example.fasta \\
    --dynamine \\
    --efoldmine \\
    --disomine \\
    --agmata \\
    --psper \\
    --fetchStructures

================================================================================
                                LIST OF PARAMETERS
================================================================================
                                GENERAL

Launch dir      : $launchDir
Prject dir      : $projectDir
Execution time  : $params.executionTimestamp
Compressed file : $params.compressedFile
================================================================================
                                INPUT FILES

Input-file (--targetSequences) : $params.targetSequences
JSON page size (--groupBy)     : $params.groupBy
================================================================================
                                OUTPUT FILES

Plots (--plotBiophysicalFeatures) : $params.plotBiophysicalFeatures
================================================================================
                                PREDICTORS

DynaMine (--dynamine)                : $params.dynamine
DisoMine (--disomine)                : $params.disomine
EFoldMine (--efoldmine)              : $params.efoldmine
AgMata (--agmata)                    : $params.agmata
PSP (--psper)                        : $params.psper
Fetch structures (--fetchStructures) : $params.fetchStructures
================================================================================
                                SINGLE SEQ or MSA

MSA mode      : No
Align for MSA : No
================================================================================
"""

if (!params.dynamine && !params.disomine && !params.efoldmine && !params.agmata && !params.psper) {
    error "Error: at least one biophysical predictor must be enabled."
}

// Processing input
allSequences
    .splitFasta( record: [header: true, seqString: true, sequence: true ])
    .branch {
        valid:  it.seqString.size() >= 5 && it.seqString.size() < 2000
        invalid: it.seqString.size() < 5 || it.seqString.size() >= 2000
    }.set { result }

result.valid.view { "VALID >${it.header}" }
result.invalid.view { "INVALID >${it.header}" }

sequencesSanitized = result.valid.map { record -> [
    header: record.header.replaceAll("[^a-zA-Z0-9]", "_"),
    seqString: record.seqString.toUpperCase(),
    sequence: record.sequence.toUpperCase()
] }

sequencesFiltered = sequencesSanitized.collectFile(name: "${targetSequencesFile.baseName}_filtered.fasta", newLine: true) {
    item -> '>' + item.header + '\n' + item.sequence + '\n'
}

sequencesRemoved = result.invalid.collectFile(name: "${targetSequencesFile.baseName}_sequences_ignored.fasta", newLine: true) {
    item -> '>' + item.header + '\n' + item.sequence + '\n'
}

sequencesGrouped = sequencesFiltered.splitFasta(record: [header: true, seqString: true, sequence: true ])
    .collectFile(name: "${targetSequencesFile.baseName}.fasta", newLine: true) {
        item -> '>' + item.header + '\n' + item.sequence + '\n'
    }
    .splitFasta( file: true, by: params.groupBy )

// Modules
include {
    fetchEsmAtlasStructure;
} from "${launchDir}/modules/structures"

include {
    plotBiophysicalFeatures;
} from "${launchDir}/modules/plots"

include {
    predictBiophysicalFeatures;
} from "${launchDir}/modules/singleSequenceAnalysis"

include { compressDirectory } from "${launchDir}/modules/utils"

// Workflows
workflow {
    predictBiophysicalFeatures(
        sequencesGrouped,
        params.dynamine,
        params.efoldmine,
        params.disomine,
        params.agmata,
        params.psper
    )

    if (params.plotBiophysicalFeatures) {
        plotBiophysicalFeatures(
            predictBiophysicalFeatures.out.predictions,
            params.dynamine,
            params.efoldmine,
            params.disomine,
            params.agmata,
            params.psper
        )

        plots = plotBiophysicalFeatures.out.plots
    } else {
        plots = Channel.empty()
    }

    if (params.fetchStructures) {
        fetchEsmAtlasStructure(sequencesSanitized.map { record -> [header: record.header, seqString: record.seqString.take(400)] })
        structures = fetchEsmAtlasStructure.out.esmStructures
    } else {
        structures = Channel.empty()
    }

    filesToCompress = Channel.empty().mix(
        predictBiophysicalFeatures.out.predictions,
        predictBiophysicalFeatures.out.index,
        plots,
        structures,
        sequencesRemoved
    ).collect()

    compressDirectory(params.compressedFile, filesToCompress)
}

workflow.onComplete {
    println "Pipeline completed at               : $workflow.complete"
    println "Time to complete workflow execution : $workflow.duration"
    println "Execution status                    : ${workflow.success ? 'Success' : 'Failed' }"
    println "Compressed file                     : $projectDir/${params.compressedFile}.tar.gz"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
    println "Details: \n ${workflow.errorReport}"
}
