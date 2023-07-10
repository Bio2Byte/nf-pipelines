process fetchEsmAtlasStructure {
    tag "${header}"
    errorStrategy 'ignore'

    input:
    tuple val(header), val(sequence)

    output:
    path "${header}.pdb", emit: esmStructures

    script:
    """
    echo Folding sequence using ESM Atlas
    sleep \$((RANDOM % 30))

    curl -X POST --data "${sequence.replaceAll('-', '').replaceAll('_', '')}" https://api.esmatlas.com/foldSequence/v1/pdb/ > ${header}.pdb

    find . -type f -name "*.pdb" -size -2k -exec bash -c 'mv "\$1" "\${1%.pdb}".fail' - '{}' +
    """
}
