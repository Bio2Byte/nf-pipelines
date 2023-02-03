process fetchStructure {
    publishDir "${resultsDirectory}", mode: 'copy'
    tag "${id}_${desc}.pdb"
    errorStrategy 'ignore'

    input:
    path resultsDirectory
    tuple val(id), val(desc), val(sequence)

    output:
    path "*.pdb", emit: esmStructures

    script:
    """
    echo Folding sequence using ESM Atlas

    curl -X POST --data "$sequence" https://api.esmatlas.com/foldSequence/v1/pdb/ > ${id}_${desc}.pdb
    """
}
