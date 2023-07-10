#!/usr/bin/env nextflow

targetSequencesFile     = file(params.targetSequences)
allSequences            = Channel.fromPath(params.targetSequences)
params.compressedFile   = "${targetSequencesFile.simpleName}_${params.executionTimestamp}"

log.info """\

Usage:

\$ nextflow run pipeline.nf \\
    -profile standard,withdocker \\
    --targetSequences ../input_example.fasta \\
    --type 'aa' or 'nuc' \\
    --qc \\ 
    --clustering 0.85\\
    --relabel \\
    --alignSequences \\
    --efoldmine \\
    --disomine \\
    --agmata \\
    --fetchStructures \\
    --buildTreeEvo \\   
    --outGroup 'Species name to root your tree on' \\
    --csubst \\
    --branchIds '1,2,3'\\
    --eteEvol 'M7,M8' \\
    --selectedProteins 'your,proteins,as,str' \\
    --plotBiophysicalFeatures \\
    --buildLogo \\
    --plotTree

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

Input file (--targetSequences) : $params.targetSequences
Input file type (--type )      : $params.type  
================================================================================
                                FILTERING

Set minimal ooccupancy of position in MSA  (--qc 0.85)                            : $params.qc     
Clustering with CD-Hit (Siminarity percentage [0.85], word length 6)(--clustering): $params.clustering
Adapt labels to clustering (--relabel)                                            : $params.relabel
================================================================================
                                PREDICTORS

Align sequences for AA with Clustal, 
    for Nuc with MACSE(--alignSequences)        : $params.alignSequences

DynaMine                             : ALWAYS
DisoMine (--disomine)                : $params.disomine
EFoldMine (--efoldmine)              : $params.efoldmine
AgMata (--agmata)                    : $params.agmata
PSP                                  : NOT IMPLEMENTED

Fetch structures (--fetchStructures) : $params.fetchStructures

Phylo. Tree (--buildTreeEvo)         : $params.buildTreeEvo
Species name to root your tree on (--outGroup)   : $params.outGroup
Csubst (--csubst)                    : $params.csubst
CsubstSite (--branchIds)             : $params.branchIds
EteEvol (--eteEvol)                  : $params.eteEvol
================================================================================
                                OUTPUT FILES


Proteins to be highlighted in the plots (--selectedProteins)    : $params.selectedProteins
Plot B2btools (--plotBiophysicalFeatures)                       : $params.plotBiophysicalFeatures
Logo (--buildLogo)                : $params.buildLogo
Phylo. Tree plot (--plotTree)     : $params.plotTree
================================================================================
"""

// Processing input


allSequences.view()


seqsFiltered = allSequences
        .splitFasta( record: [header: true,  sequence: true ])
        .map { record -> [header: record.header.replaceAll("[^a-zA-Z0-9]", "_"),
            sequence: record.sequence.replaceAll("\n","").replaceAll("[^a-zA-Z-]", "X")] }

seqsQC = seqsFiltered.branch{
        valid:      it.sequence.count('X') / it.sequence.size() < 0.5
        invalid:    it.sequence.count('X') / it.sequence.size() >= 0.5
    }.set { result}

/* replace with other filters, eg:
        it.header.contains(params.filterValid) == false
        it.header.contains("_Pro_") == false
        it.seqString.size() >= 60 && it.seqString.size() < 2000
*/
result.invalid.view{ "INVALID >${it.header}" }

sequencesSanitized = result.valid

sequencesValid = sequencesSanitized.collectFile(name: "${targetSequencesFile.baseName}_filtered.fasta", newLine: true) {
    item -> '>' + item.header + '\n' + item.sequence + '\n'
}

sequencesRemoved = result.invalid.collectFile(name: "${targetSequencesFile.baseName}_sequences_ignored.fasta", newLine: true) {
    item -> '>' + item.header + '\n' + item.sequence + '\n'
}

// Modules
include {
    cdHitClustering;
    postClusteringLabels
    } from "${projectDir}/modules/cdhit"

include {
    predictBiophysicalFeatures;
    buildMultipleSequenceAlignmentAA;
    buildMultipleSequenceAlignmentNuc;
    takeMultipleSequenceAlignment;
    takeMultipleSequenceAlignmentNuc;
    buildPhylogeneticTree;
    buildPhylogeneticTreeEvol;
    buildLogo;
} from "${projectDir}/modules/multipleSequenceAlignment"

include {
    fetchEsmAtlasStructure
} from "${projectDir}/modules/structures"

include {
    plotBiophysicalFeatures;
    plotPhylogeneticTree;
} from "${projectDir}/modules/plots"

include {
    findRoot;
    runCsubst;
    runCsubstBranch;
    runEteEvol;
} from "${projectDir}/modules/dndsCsubst"

include { compressDirectory } from "${projectDir}/modules/utils"


workflow {

    if (params.clustering){
        cdHitClustering(sequencesValid, params.clustering)
        clusters =  cdHitClustering.out.clusters

        postClusteringLabels(sequencesValid,clusters, params.relabel)
        sequencesFiltered = postClusteringLabels.out.repSeqs
        representativeRepresented =  postClusteringLabels.out.representativeRepresented
    } else {
        sequencesFiltered = sequencesValid
        clusters = Channel.empty()
        representativeRepresented = Channel.empty()
    }

    if (params.alignSequences){
        if (params.type == 'nuc') {
            buildMultipleSequenceAlignmentNuc(sequencesFiltered)

            msaNuc = buildMultipleSequenceAlignmentNuc.out.msaNuc
            multipleSequenceAlignment = buildMultipleSequenceAlignmentNuc.out.multipleSequenceAlignment

            takeMultipleSequenceAlignmentNuc(msaNuc, params.buildTreeEvo, params.qc)
            multipleSequenceAlignmentNuc = takeMultipleSequenceAlignmentNuc.out.multipleSequenceAlignmentNuc
            
        } 
        if (params.type == 'aa') {
            multipleSequenceAlignmentNuc = Channel.empty()

            buildMultipleSequenceAlignmentAA(sequencesFiltered)
            msaAA = buildMultipleSequenceAlignmentAA.out.multipleSequenceAlignment
            takeMultipleSequenceAlignment(msaAA, params.buildTreeEvo, params.qc)
            multipleSequenceAlignment = takeMultipleSequenceAlignment.out.multipleSequenceAlignment

        } 
    } else {
        if (params.type == 'nuc') {
            multipleSequenceAlignment = Channel.empty()

            msaNuc =  sequencesFiltered
            takeMultipleSequenceAlignmentNuc(msaNuc, params.buildTreeEvo, params.qc)
            multipleSequenceAlignmentNuc = takeMultipleSequenceAlignmentNuc.out.multipleSequenceAlignmentNuc

        }
        if (params.type == 'aa') {
            multipleSequenceAlignmentNuc = Channel.empty()
            msaAA = sequencesFiltered
            takeMultipleSequenceAlignment(msaAA, params.buildTreeEvo, params.qc)
            multipleSequenceAlignment = takeMultipleSequenceAlignment.out.multipleSequenceAlignment

        }
        
    }


    if (params.buildTree) {
        buildPhylogeneticTree(multipleSequenceAlignment)
        phylogeneticTree = buildPhylogeneticTree.out.tree

        iqtreeFiles = Channel.empty()

        if (params.plotTree) {
            plotPhylogeneticTree(phylogeneticTree, params.plotTree )
            plottedPhylogeneticTree = plotPhylogeneticTree.out.treePlot
        } else {
            println "Skipping Phylogenetic tree plot from MSA: ${targetSequencesFile}"

            plottedPhylogeneticTree = Channel.empty()
        }
    } else if (params.buildTreeEvo) {
        buildPhylogeneticTreeEvol(multipleSequenceAlignmentNuc)
        phylogeneticTree = buildPhylogeneticTreeEvol.out.tree
        iqtreeFiles = buildPhylogeneticTreeEvol.out.iqtreefiles

        if (params.plotTree) {
            plotPhylogeneticTree(phylogeneticTree, params.plotTree)
            plottedPhylogeneticTree = plotPhylogeneticTree.out.treePlot
        } else {
            println "Skipping Phylogenetic tree plot from MSA: ${targetSequencesFile}"

            plottedPhylogeneticTree = Channel.empty()
        }
    } else {
        println "Skipping Phylogenetic tree from MSA: ${targetSequencesFile}"
        println "Skipping Phylogenetic tree plot from MSA: ${targetSequencesFile}"

        phylogeneticTree = Channel.empty()
        iqtreeFiles = Channel.empty()
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
            params.disomine,
            params.agmata,
            params.selectedProteins,
            predictBiophysicalFeatures.out.predictions,
            predictBiophysicalFeatures.out.stats,
            params.plotBiophysicalFeatures
        )
        plottedBiophysicalFeaturesInPNG = plotBiophysicalFeatures.out.plots
        
    } else {
        println "Skipping Biophysical features plots from MSA: ${targetSequencesFile}"
        plottedBiophysicalFeaturesInPNG = Channel.empty()
    }

    if (params.fetchStructures) {
        // Sequences need to be shorter then 400 residues!)

        fetchEsmAtlasStructure(multipleSequenceAlignment.splitFasta(record: [header: true, seqString: true]))
        structures = fetchEsmAtlasStructure.out.esmStructures
    } else {
        sequencesFiltered.view{ "Skipping fetching structures from EsmAtlas for sequences: "+ it}
        structures = Channel.empty()
    }

    if (params.csubst) {

        findRoot(phylogeneticTree, params.outGroup)
        rootedTree = findRoot.out.rootedTree

        runCsubst(multipleSequenceAlignmentNuc, rootedTree)
        csubstOutZip = runCsubst.out.csubstOut
    } else{
        rootedTree = Channel.empty()
        csubstOutZip = Channel.empty()
    }

    if (params.branchIds) {
        runCsubstBranch(multipleSequenceAlignmentNuc, rootedTree, csubstOutZip, params.branchIds)
        csubstBranchOutZip = runCsubstBranch.out.csubstBranchOut
    } else{
        csubstBranchOutZip = Channel.empty()
    }

    if (params.eteEvol) {
        if (params.csubst == false){
            findRoot(phylogeneticTree, params.outGroup )
            rootedTree = findRoot.out.rootedTree
        }


        modelList = params.eteEvol?.split(',') as List
        modelChannel = Channel.fromList(modelList)


        runEteEvol(multipleSequenceAlignmentNuc, rootedTree, modelChannel)
        eteOutZip = runEteEvol.out.eteOut.toList()

    } else{
        rootedTree = Channel.empty()
        eteOutZip = Channel.empty()
    }


    filesToCompress = Channel.empty().mix(
        sequencesFiltered,
        clusters,
        representativeRepresented,
        multipleSequenceAlignment,
        multipleSequenceAlignmentNuc,
        iqtreeFiles,
        plottedPhylogeneticTree,
        rootedTree,
        logo,
        predictBiophysicalFeatures.out.predictions,
        predictBiophysicalFeatures.out.stats,
        plottedBiophysicalFeaturesInPNG,
        csubstOutZip,
        csubstBranchOutZip,
        eteOutZip,
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
