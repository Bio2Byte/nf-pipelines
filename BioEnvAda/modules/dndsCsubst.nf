process findRoot{
    tag "${outGroup}"
    input:
    path phylogeneticTree
    val outGroup

    output:
    path '*rooted.treefile', emit: rootedTree

    script:
    """
    python  $projectDir/bin/rootTree.py $phylogeneticTree $outGroup 
    """
}

process runCsubst{

    label  'tinysnail'
    tag "${multipleSequenceAlignmentNuc.name}"

    input:
    path multipleSequenceAlignmentNuc
    path rootedTree

    output:
    path "csubst_${multipleSequenceAlignmentNuc}_out.tar.gz" , emit: csubstOut

    script:
    """
    csubst analyze --alignment_file ${multipleSequenceAlignmentNuc}  --rooted_tree_file ${rootedTree} --iqtree_redo no 
    
    touch csubst_${multipleSequenceAlignmentNuc}_out.tar.gz
    tar -cvzf csubst_${multipleSequenceAlignmentNuc}_out.tar.gz --exclude={'./.*','csubst_${multipleSequenceAlignmentNuc}_out.tar.gz'} .
    """
} 

process runCsubstBranch{
    label  'tinysnail'
    tag "${multipleSequenceAlignmentNuc.name}"

    input:
    path multipleSequenceAlignmentNuc
    path rootedTree
    path csubstOutZip
    val branchIds

    output:
    path "csubst_branch_${multipleSequenceAlignmentNuc}_out.tar.gz" , emit: csubstBranchOut

    script:
    """
    tar -xvzf ${csubstOutZip}

    csubst site --alignment_file ${multipleSequenceAlignmentNuc}  --rooted_tree_file ${rootedTree} --branch_id ${branchIds} 
    
    touch csubst_branch_${multipleSequenceAlignmentNuc}_out.tar.gz
    tar -cvzhf csubst_branch_${multipleSequenceAlignmentNuc}_out.tar.gz --exclude={'./.*','*.tar.gz','csubst_branch_${multipleSequenceAlignmentNuc}_out.tar.gz'} .
    """
}


process runEteEvol{
    errorStrategy 'ignore'
    label  'tinysnail'
    tag "${multipleSequenceAlignmentNuc.name}"
    conda '/Users/sophie/miniconda3/envs/ete3'

    input:
    path multipleSequenceAlignmentNuc
    path rootedTree
    each eteEvol 

    output:
    path "ete_${multipleSequenceAlignmentNuc}_*_out.tar.gz" , emit: eteOut

    script:

    """
    python $projectDir/bin/EteEvol.py ${multipleSequenceAlignmentNuc.name} $rootedTree $multipleSequenceAlignmentNuc $eteEvol
    
    tar -cvzhf ete_${multipleSequenceAlignmentNuc}_${eteEvol}_out.tar.gz pamlwd/ plots/
    """



}
