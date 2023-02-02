process fetchStructure {
    publishDir "${resultsDirectoryPath}", mode: 'copy'

    tag "$id"

    input:
    path resultsDirectoryPath
    tuple val(id), val(sequence)

    output:
    path "*.pdb", emit: esmStructures

    when:
    params.fetchStructures

    script:
    """
    echo Folding sequence using ESM Atlas
    curl -X POST --data "$sequence" https://api.esmatlas.com/foldSequence/v1/pdb/ > ${id}.pdb
    echo Preview of the PDB content
    head ${id}.pdb
    tail ${id}.pdb
    """
}
