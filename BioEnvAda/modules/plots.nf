
process plotBiophysicalFeatures{
    tag "${msa.name}"
    
    input:
    path msa
    val efoldmine
    val disomine
    val agmata
    val selected_proteins
    path predictions
    path stats
    val split

    output:
    path "*.png", emit: plots

    script:
    """
    python3 $projectDir/bin/MsaPlotB2btoolsBar.py $msa $predictions "${efoldmine ? 'efoldmine,' : ''}${disomine ? 'disomine,' : ''}${agmata ? 'agmata,' : ''}" $selected_proteins $stats $split
    """
}


process plotPhylogeneticTree {
    tag "${tree.name}"

    input:
    path tree
    val plotTree

    output:
    path "${tree}.png", emit: treePlot

    script:
    """
    python3 $projectDir/bin/treePlot.py $tree $plotTree
    """
}
