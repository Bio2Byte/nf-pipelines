process fetchEsmAtlasStructure {
    tag "${header}"
    errorStrategy 'ignore'
    debug true

    input:
    tuple val(header), val(sequence)

    output:
    path "${header}.pdb", emit: esmStructures

    script:
    """
    echo Folding sequence using ESM Atlas
    echo To fetch: ${sequence.replaceAll('-', '').replaceAll('_', '')}

    curl -X POST --data "${sequence.replaceAll('-', '').replaceAll('_', '')}" https://api.esmatlas.com/foldSequence/v1/pdb/ > ${header}.pdb
    """
}
