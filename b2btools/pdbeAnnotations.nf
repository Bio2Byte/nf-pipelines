params.targetPdbIds = "$launchDir/pdb_list.txt"
targetPdbIdsFile = file(params.targetPdbIds)
allPdbIds = Channel.fromPath(params.targetPdbIds).splitCsv(header: false, strip: true).flatten().unique()

process fetchPdbSequences {
    publishDir "${projectDir}/results/${pdbId[1]}${pdbId[2]}/", mode: 'symlink'
    tag "${pdbId}.fasta"
    errorStrategy 'ignore'

    input:
        val pdbId
    output:
        path "${pdbId}.fasta", emit: entrySequences

    script:
    """
    curl https://www.rcsb.org/fasta/entry/${pdbId}/download/ --silent --output ${pdbId}.fasta
    """
}

process predictBiophysicalFeatures {
    publishDir "${projectDir}/results/${entrySequences.simpleName[1]}${entrySequences.simpleName[2]}/", mode: 'symlink'
    tag "${entrySequences.simpleName}.json"

    input:
        path entrySequences
    output:
        path "${entrySequences.simpleName}.json", emit: entryPredictions

    script:
    """
    #!/usr/local/bin/python
    import json
    from b2bTools import SingleSeq

    single_seq = SingleSeq("$entrySequences")
    single_seq.predict(tools=['dynamine', 'efoldmine', 'disomine'])
    all_predictions = single_seq.get_all_predictions()
    json.dump(all_predictions, open('${entrySequences.simpleName}.json', 'w'), indent=2)
    """
}

process formatFunPdbePrediction {
    publishDir "${projectDir}/results/${entryPredictions.simpleName[1]}${entryPredictions.simpleName[2]}/", mode: 'symlink'
    tag "${entryPredictions.simpleName}.funpdbe.json"
    cache false

    input:
        path entryPredictions
    output:
        path "${entryPredictions.simpleName}.funpdbe.json", emit: entryFunPdbePredictions

    script:
    """
    #!/usr/bin/python3
    import json
    import re

    predictions_dict = {}
    aa_dict = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
            'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
            'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
            'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}

    def mapResidueToFunPdbeFormat(b2b_predictions_dict, sequence_id):
        chain_predictions = b2b_predictions_dict[sequence_id]
        residues_list = []
        for position, residue in enumerate(chain_predictions['seq']):
            site_data_dict = []
            if 'backbone' in chain_predictions:
                backbone = chain_predictions['backbone'][position]

                site_data_dict.append({
                    'confidence_classification': 'null',
                    'confidence_score': 0.5,
                    'raw_score': round(backbone, 3) if isinstance(backbone, float) else None,
                    'site_id_ref': 1
                })

            if 'backbone' in chain_predictions:
                sidechain = chain_predictions['sidechain'][position]

                site_data_dict.append({
                    'confidence_classification': 'null',
                    'confidence_score': 0.5,
                    'raw_score': round(sidechain, 3) if isinstance(sidechain, float) else None,
                    'site_id_ref': 2
                })

            if 'earlyFolding' in chain_predictions:
                earlyFolding = chain_predictions['earlyFolding'][position]

                site_data_dict.append({
                    'confidence_classification': 'null',
                    'confidence_score': 0.5,
                    'raw_score': round(earlyFolding, 3) if isinstance(earlyFolding, float) else None,
                    'site_id_ref': 3
                })

            if 'disoMine' in chain_predictions:
                disoMine = chain_predictions['disoMine'][position]
                site_data_dict.append({
                    'confidence_classification': 'null',
                    'confidence_score': 0.5,
                    'raw_score': round(disoMine, 3) if isinstance(disoMine, float) else None,
                    'site_id_ref': 4
                })

            residue_dict = {
                'aa_type': aa_dict.get(residue.upper(), 'UNK'),
                'pdb_res_label': f'{position + 1}',
                'site_data': site_data_dict
            }

            residues_list.append(residue_dict)

        return residues_list


    with open('${entryPredictions}', 'r') as input_file:
        b2b_predictions_dict = json.loads(input_file.read())

    sequence_chains = []
    for sequence_id in b2b_predictions_dict.keys():

        pattern = r"_([A-Za-z])_"
        chains_list = re.findall(pattern, sequence_id)

        if chains_list:
            residues = mapResidueToFunPdbeFormat(b2b_predictions_dict, sequence_id)
            for chain in chains_list:
                chain_dict = {
                    'chain_label': chain,
                    'residues': residues
                }
                sequence_chains.append(chain_dict)


    predictions_dict['data_resource'] = 'dynamine'
    predictions_dict['resource_entry_url'] = 'http://dynamine.ibsquare.be/'
    predictions_dict['pdb_id'] = '${entryPredictions.simpleName}'
    predictions_dict['chains'] = sequence_chains
    predictions_dict['evidence_code_ontology'] = [
        { 'eco_code': 'ECO_0000203' },
        { 'eco_code': 'ECO_0000364' }
    ]
    predictions_dict['sites'] = [
        {
            'label': 'backbone',
            'site_id': 1,
            'source_accession': '${entryPredictions.simpleName}',
            'source_database': 'PDB',
            'source_release_date': '12/06/2018',
            'site_url': 'https://bio2byte.be/b2btools/dynamine'
        },
        {
            'label': 'sidechain',
            'site_id': 2,
            'source_accession': '${entryPredictions.simpleName}',
            'source_database': 'PDB',
            'source_release_date': '12/06/2018',
            'site_url': 'https://bio2byte.be/b2btools/dynamine'
        },
        {
            'label': 'efoldmine',
            'site_id': 3,
            'source_accession': '${entryPredictions.simpleName}',
            'source_database': 'PDB',
            'source_release_date': '12/06/2018',
            'site_url': 'https://bio2byte.be/b2btools/efoldmine'
        },
        {
            'label': 'disomine',
            'site_id': 4,
            'source_accession': '${entryPredictions.simpleName}',
            'source_database': 'PDB',
            'source_release_date': '12/06/2018',
            'site_url': 'https://bio2byte.be/b2btools/disomine'
        }
    ]
    predictions_dict['software_version'] = '3.0.5'

    with open('${entryPredictions.simpleName}.funpdbe.json', 'w') as output_file:
        json.dump(predictions_dict, output_file, indent=2)
    """
}

workflow {
    fetchPdbSequences(allPdbIds)
    predictBiophysicalFeatures(fetchPdbSequences.out.entrySequences)
    formatFunPdbePrediction(predictBiophysicalFeatures.out.entryPredictions)
}
