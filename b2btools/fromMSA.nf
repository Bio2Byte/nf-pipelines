#!/usr/bin/env nextflow

// RUN: nextflow run fromMSA.nf --targetSequences ../input_example.fasta -profile standard,withdocker --efoldmine --disomine --alignSequences
params.executionTimestamp = new java.text.SimpleDateFormat("yyyy_MM_dd_HH_mm_ss").format(new Date())

// Available predictors:
params.dynamine        = true
params.efoldmine       = false
params.disomine        = false
params.fetchStructures = false

// MSA parameters:
params.alignSequences = false

// Input files
params.targetSequences  = "$launchDir/example.fasta"
targetSequencesFile     = file(params.targetSequences)
allSequences            = Channel.fromPath(params.targetSequences)

params.compressedFile   = "${targetSequencesFile.simpleName}_${params.executionTimestamp}"

// Output files
params.plotBiophysicalFeatures = false
params.buildLogo               = false
params.buildTree               = false
params.plotTree                = false

log.info """\

Usage:

\$ nextflow run fromMSA.nf \\
    -with-report report_fromMsa.html \\
    -with-dag flowchart_fromMSA.png \\
    -profile standard,withdocker \\
    --targetSequences ../input_example.fasta \\
    --plotBiophysicalFeatures \\
    --buildLogo \\
    --buildTree \\
    --plotTree \\
    --efoldmine \\
    --disomine \\
    --fetchStructures \\
    --alignSequences

================================================================================
                                LIST OF PARAMETERS
================================================================================
                                GENERAL

Launch dir      : $launchDir
Project dir     : $projectDir
Execution time  : $params.executionTimestamp
Compressed file : ${params.compressedFile}.tar.gz
================================================================================
                                INPUT FILES

Input-file (--targetSequences) : $params.targetSequences
================================================================================
                                OUTPUT FILES

Plots (--plotBiophysicalFeatures) : $params.plotBiophysicalFeatures
Logo (--buildLogo)                : $params.buildLogo
Phylo. Tree (--buildTree)         : $params.buildTree
Phylo. Tree plot (--plotTree)     : $params.plotTree
================================================================================
                                PREDICTORS

DynaMine                             : true
DisoMine (--disomine)                : $params.disomine
EFoldMine (--efoldmine)              : $params.efoldmine
PSP                                  : NOT IMPLEMENTED
Fetch structures (--fetchStructures) : $params.fetchStructures
================================================================================
                                SINGLE SEQ or MSA

MSA mode                         : true
Align for MSA (--alignSequences) : $params.alignSequences
================================================================================
"""

// Processing input
// TO BE DISCUSSED
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

// Modules
include {
    predictBiophysicalFeatures;
    buildMultipleSequenceAlignment;
    takeMultipleSequenceAlignment;
    buildPhylogeneticTree;
    buildLogo;
} from "${projectDir}/modules/multipleSequenceAlignment"

include {
    fetchEsmAtlasStructure
} from "${projectDir}/modules/structures"

include {
    plotBiophysicalFeaturesOverview as plotBiophysicalFeatures;
    plotPhylogeneticTree;
} from "${projectDir}/modules/plots"

include { compressDirectory } from "${projectDir}/modules/utils"

// Workflows
workflow {
    if (params.alignSequences) {
        buildMultipleSequenceAlignment(sequencesFiltered)
        multipleSequenceAlignment = buildMultipleSequenceAlignment.out.multipleSequenceAlignment
    } else {
        takeMultipleSequenceAlignment(sequencesFiltered)
        multipleSequenceAlignment = takeMultipleSequenceAlignment.out.multipleSequenceAlignment
    }

    if (params.buildTree) {
        buildPhylogeneticTree(multipleSequenceAlignment)
        phylogeneticTree = buildPhylogeneticTree.out.tree

        if (params.plotTree) {
            plotPhylogeneticTree(phylogeneticTree)
            plottedPhylogeneticTree = plotPhylogeneticTree.out.treePlot
        } else {
            println "Skipping Phylogenetic tree plot from MSA: ${targetSequencesFile}"

            plottedPhylogeneticTree = Channel.empty()
        }
    } else {
        println "Skipping Phylogenetic tree from MSA: ${targetSequencesFile}"
        println "Skipping Phylogenetic tree plot from MSA: ${targetSequencesFile}"

        phylogeneticTree = Channel.empty()
        plottedPhylogeneticTree = Channel.empty()
    }

    if (params.buildLogo) {
        buildLogo(multipleSequenceAlignment)
        logo = buildLogo.out.logo
    } else {
        println "Skipping Logo from MSA: ${targetSequencesFile}"

        logo = Channel.empty()
    }

    predictBiophysicalFeatures(multipleSequenceAlignment)
    if (params.plotBiophysicalFeatures) {
        plotBiophysicalFeatures(
            multipleSequenceAlignment,
            params.efoldmine,
            params.disomine
        )
        plottedBiophysicalFeaturesInPNG = plotBiophysicalFeatures.out.plots
        plottedBiophysicalFeaturesInPDF = plotBiophysicalFeatures.out.documents
    } else {
        println "Skipping Biophysical features plots from MSA: ${targetSequencesFile}"

        plottedBiophysicalFeaturesInPNG = Channel.empty()
        plottedBiophysicalFeaturesInPDF = Channel.empty()
    }

    if (params.fetchStructures) {
        // TO BE FIXED WHEN MSA input file (400 residues)
        fetchEsmAtlasStructure(sequencesSanitized.map { record -> [header: record.header, seqString: record.seqString.take(400)] })
        structures = fetchEsmAtlasStructure.out.esmStructures
    } else {
        println "Skipping fetching structures from EsmAtlas for sequences: ${targetSequencesFile}"

        structures = Channel.empty()
    }

    filesToCompress = Channel.empty().mix(
        multipleSequenceAlignment,
        phylogeneticTree,
        plottedPhylogeneticTree,
        logo,
        predictBiophysicalFeatures.out.predictions,
        plottedBiophysicalFeaturesInPNG,
        plottedBiophysicalFeaturesInPDF,
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
