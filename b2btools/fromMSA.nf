#!/usr/bin/env nextflow

// RUN: nextflow run fromMSA.nf --targetSequences ../input_example.fasta -profile standard,withdocker --efoldmine --disomine --alignSingleSequences
params.executionTimestamp = new java.text.SimpleDateFormat("yyyy_MM_dd_HH_mm_ss").format(new Date())

// Available predictors:
params.efoldmine       = false
params.disomine        = false
params.fetchStructures = false

// MSA parameters:
params.alignSingleSequences = false

// Input files
params.targetSequences  = "$launchDir/example.fasta"
allSequences            = Channel.fromPath(params.targetSequences)
targetSequencesFile     = file(params.targetSequences)
params.resultsDirectory = "$projectDir/${params.executionTimestamp}"
params.compressedFile   = "$projectDir/${targetSequencesFile.simpleName}_${params.executionTimestamp}.tar.gz"
params.groupBy          = 10

// Output files
params.plotBiophysicalFeatures = false
params.buildLogo               = false
params.buildTree               = false
params.plotTree                = false

log.info """\

Usage:

\$ nextflow run fromMSA.nf \\
    -with-report \\
    -with-dag flowchart.png \\
    -profile standard,withdocker \\
    --targetSequences ../input_example.fasta \\
    --plotBiophysicalFeatures \\
    --buildLogo \\
    --buildTree \\
    --plotTree \\
    --efoldmine \\
    --disomine \\
    --fetchStructures \\
    --alignSingleSequences

================================================================================
                                LIST OF PARAMETERS
================================================================================
                                GENERAL

Launch dir      : $launchDir
Project dir     : $projectDir
Execution time  : $params.executionTimestamp
Results dir     : $params.resultsDirectory
Compressed file : $params.compressedFile
================================================================================
                                INPUT FILES

Input-file (--targetSequences) : $params.targetSequences
JSON page size (--groupBy)     : $params.groupBy
================================================================================
                                OUTPUT FILES

Plots (--plotBiophysicalFeatures) : $params.plotBiophysicalFeatures
Logo (--buildLogo)                : $params.buildLogo
Phylo. Tree (--buildTree)         : $params.buildTree
Phylo. Tree plot (--plotTree)     : $params.plotTree
================================================================================
                                PREDICTORS

DynaMine                             : ALWAYS
DisoMine (--disomine)                : $params.disomine
EFoldMine (--efoldmine)              : $params.efoldmine
PSP                                  : NOT IMPLEMENTED
Fetch structures (--fetchStructures) : $params.fetchStructures
================================================================================
                                SINGLE SEQ or MSA

MSA mode                               : true
Align for MSA (--alignSingleSequences) : $params.alignSingleSequences
================================================================================
"""

// Processing input
allSequences
    .splitFasta( record: [id: true, seqString: true, sequence: true, desc: true ])
    .branch {
        valid:  it.seqString.size() >= 5 && it.seqString.size() < 2000
        invalid: it.seqString.size() < 5 || it.seqString.size() >= 2000
    }.set { result }

sequencesRemoved = result.invalid
    .collectFile(name: "${params.resultsDirectory}/${targetSequencesFile.baseName}_sequences_ignored.fasta", newLine: true) {
        item -> '>' + item.id + '_' + item.desc + '\n' + item.sequence + '\n'
    }.subscribe {
        println "Ignored entries are saved to file: $it"
    }

sequencesSanitized = result.valid.map { record -> [id: record.id.replaceAll("[^a-zA-Z0-9]", "_"), desc: record.desc.replaceAll("[^a-zA-Z0-9]", "_"), seqString: record.seqString, sequence: record.sequence] }
sequencesFiltered = sequencesSanitized.collectFile(name: "${targetSequencesFile.baseName}_filtered.fasta", newLine: true) {
    item -> '>' + item.id + '_' + item.desc + '\n' + item.sequence + '\n'
}

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
    plotBiophysicalFeaturesOverview as plotBiophysicalFeatures;
    plotPhylogeneticTree;
} from "${launchDir}/modules/plots"

// Processes
process compressPredictions {
    debug true

    input:
    val multipleSequenceAlignment
    val tree
    val treePlot
    val logo
    val predictions
    val plots
    val documents
    val esmStructures
    // val sequencesRemoved

    script:
    """
    echo Creating compressed file: ${params.compressedFile}
    echo Content of compressed: ${params.resultsDirectory}

    # tar -czvf ${params.compressedFile} -C / ${params.resultsDirectory}
    """
}

// Workflows
workflow {
    if (params.alignSingleSequences) {
        buildMultipleSequenceAlignment(params.resultsDirectory, sequencesFiltered)
        multipleSequenceAlignment = buildMultipleSequenceAlignment.out.multipleSequenceAlignment
    } else {
        takeMultipleSequenceAlignment(allSequences)
        multipleSequenceAlignment = takeMultipleSequenceAlignment.out.multipleSequenceAlignment
    }

    if (params.buildTree) {
        buildPhylogeneticTree(params.resultsDirectory, multipleSequenceAlignment)
        phylogeneticTree = buildPhylogeneticTree.out.tree

        if (params.plotTree) {
            plotPhylogeneticTree(params.resultsDirectory, phylogeneticTree)
            plottedPhylogeneticTree = plotPhylogeneticTree.out.treePlot
        } else {
            plottedPhylogeneticTree = Channel.empty()
        }
    } else {
        phylogeneticTree = Channel.empty()
        plottedPhylogeneticTree = Channel.empty()
    }

    if (params.buildLogo) {
        buildLogo(params.resultsDirectory, multipleSequenceAlignment)
        logo = buildLogo.out.logo
    } else {
        logo = Channel.empty()
    }

    predictBiophysicalFeatures(params.resultsDirectory, multipleSequenceAlignment)
    if (params.plotBiophysicalFeatures) {
        plotBiophysicalFeatures(
            params.resultsDirectory,
            multipleSequenceAlignment,
            params.efoldmine,
            params.disomine
        )
        plottedBiophysicalFeaturesInPNG = plotBiophysicalFeatures.out.plots
        plottedBiophysicalFeaturesInPDF = plotBiophysicalFeatures.out.documents
    } else {
        plottedBiophysicalFeaturesInPNG = Channel.empty()
        plottedBiophysicalFeaturesInPDF = Channel.empty()
    }

    if (params.fetchStructures) {
        fetchStructure(params.resultsDirectory, sequencesSanitized.map { record -> [id: record.id, desc: record.desc, seqString: record.seqString.take(400)] })
        structures = fetchStructure.out.esmStructures
    } else {
        structures = Channel.empty()
    }

    compressPredictions(
        multipleSequenceAlignment,
        phylogeneticTree,
        plottedPhylogeneticTree,
        logo,
        predictBiophysicalFeatures.out.predictions,
        plottedBiophysicalFeaturesInPNG,
        plottedBiophysicalFeaturesInPDF,
        structures,
        // sequencesRemoved
    )
}

workflow.onComplete {
    println "Pipeline completed at               : $workflow.complete"
    println "Time to complete workflow execution : $workflow.duration"
    println "Execution status                    : ${workflow.success ? 'Success' : 'Failed' }"
    println "Compressed file                     : $params.compressedFile"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
