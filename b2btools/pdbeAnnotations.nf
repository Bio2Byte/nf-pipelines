params.targetPdbIds = "$launchDir/pdb_list.txt"
// targetPdbIdsFile = file(params.targetPdbIds)

println "Target: ${params.targetPdbIds}"

process extractPdb {
    publishDir "$launchDir", mode: 'symlink'
    tag "${uniprotPdbCsv.simpleName}"

    input:
        path uniprotPdbCsv

    output:
        path "${uniprotPdbCsv.simpleName}.flat.csv", emit: flatPdbCsv

    script:
    """
    #!/usr/local/bin/python
import csv

pdb_list = []
ignore_first = True
with open("${uniprotPdbCsv}", 'r') as csvfile:
    for row in csvfile.readlines():
        if ignore_first:
            ignore_first = False
            continue

        col1, *cols = row.split(";")
        pdb_list.append(col1.split(',')[1])

        for potential_pdbid in cols:
            if potential_pdbid and potential_pdbid != '' and potential_pdbid != '\\n':
                pdb_list.append(potential_pdbid)

with open("${uniprotPdbCsv.simpleName}.flat.csv", 'w') as csvfile:
    for pdb_id in pdb_list:
        csvfile.write(pdb_id + "\\n")
    csvfile.flush()
"""
}

process fetchPdb {
    // publishDir "${projectDir}/results/${pdbId[1]}${pdbId[2]}/", mode: 'symlink'
    errorStrategy 'ignore'

    tag "${pdbId}.pdb"
    input:
        val pdbId

    output:
        tuple val(pdbId), path ("${pdbId}.pdb"), emit: pdbTuple

    script:
    """
    curl -f -L --output ${pdbId}.pdb 'https://files.rcsb.org/download/${pdbId}.pdb'
    """
}

process buildFastaFromPdb {
    // publishDir "${projectDir}/results/${pdbId[1]}${pdbId[2]}/", mode: 'symlink'

    tag "${pdbId}.pdb"
    input:
        tuple val(pdbId), path(pdbFile)

    output:
        tuple val(pdbId), path(pdbFile), path("${pdbId}.fasta"), emit: entrySequencesTuple

    """
    #!/usr/local/bin/python
    from Bio.PDB import PDBParser, Polypeptide
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import generic_protein
    from Bio import SeqIO

    parser = PDBParser()
    structure = parser.get_structure("structure", "${pdbFile}")

    sequences = []
    for chain in structure.get_chains():
        seq = Polypeptide.Polypeptide(chain).get_sequence()
        seq = Seq(str(seq), generic_protein)
        seq_record = SeqRecord(seq, id='${pdbId}', description=chain.get_id())
        sequences.append(seq_record)

    SeqIO.write(sequences, "${pdbId}.fasta", "fasta")
    """
}

process fetchPdbSequences {
    // publishDir "${projectDir}/results/${pdbId[1]}${pdbId[2]}/", mode: 'symlink'
    tag "${pdbId}.fasta"

    input:
        tuple val(pdbId), path(pdbFile)

    output:
        tuple val(pdbId), path(pdbFile), path("${pdbId}.fasta"), emit: entrySequencesTuple

    script:
    """
    curl -L --silent --output ${pdbId}.fasta https://www.rcsb.org/fasta/entry/${pdbId}/download/
    """
}

process predictBiophysicalFeatures {
    // publishDir "${projectDir}/results/${pdbId[1]}${pdbId[2]}/", mode: 'symlink'
    tag "${pdbId}.json"
    errorStrategy 'ignore'

    input:
        tuple val(pdbId), path(pdbFile), path(entrySequences)

    output:
        tuple val(pdbId), path(pdbFile), path(entrySequences), path("${pdbId}.json"), emit: entryPredictionsTuple

    script:
    """
    #!/usr/local/bin/python
    import json
    from b2bTools import SingleSeq

    single_seq = SingleSeq("$entrySequences")
    single_seq.predict(tools=['dynamine', 'efoldmine'])
    all_predictions = single_seq.get_all_predictions()
    json.dump(all_predictions, open('${entrySequences.simpleName}.json', 'w'), indent=2)
    """
}

process formatFunPdbePrediction {
    publishDir "${projectDir}/results/${pdbId[1]}${pdbId[2]}/", mode: 'symlink'
    tag "${pdbId}.funpdbe.json"
    errorStrategy 'ignore'

    input:
        tuple val(pdbId), path(pdbFile), path(entrySequences), path(entryPredictions)

    output:
        path "${pdbId}.funpdbe.json", emit: entryFunPdbePredictions

    script:
    """
    #!/usr/local/bin/python
    import json
    import re
    from datetime import datetime
    from Bio.PDB import PDBParser

    parser = PDBParser()
    structure = parser.get_structure('${pdbId}', '${pdbFile}')
    deposition_date = datetime.strptime(structure.header['deposition_date'], '%Y-%m-%d')
    release_date = datetime.strptime(structure.header['release_date'], '%Y-%m-%d')

    predictions_dict = {}
    aa_dict = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
            'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
            'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
            'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}

    def mapResidueToFunPdbeFormat(b2b_predictions_dict, sequence_id, residue_numbers):
        chain_predictions = b2b_predictions_dict[sequence_id]
        residues_list = []

        for position, residue in enumerate(chain_predictions['seq']):
            aa_type = aa_dict.get(residue.upper(), 'UNK')
            if aa_type == 'UNK':
                continue

            site_data_dict = []
            if 'backbone' in chain_predictions:
                backbone = chain_predictions['backbone'][position]

                site_data_dict.append({
                    'confidence_classification': 'null',
                    'confidence_score': 0.5,
                    'raw_score': round(backbone, 3) if isinstance(backbone, float) else None,
                    'site_id_ref': 1
                })

            if 'sidechain' in chain_predictions:
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

            residue_dict = {
                'aa_type': aa_type,
                'pdb_res_label': f'{residue_numbers[position]}',
                'site_data': site_data_dict
            }

            residues_list.append(residue_dict)

        return residues_list


    with open('${entryPredictions}', 'r') as input_file:
        b2b_predictions_dict = json.loads(input_file.read())

    sequence_chains = []
    for chain in structure.get_chains():
        sequence_id = f'${pdbId}_{chain.get_id()}'

        residue_numbers = [resi.get_id()[1] for resi in chain.get_residues()]

        residues = mapResidueToFunPdbeFormat(b2b_predictions_dict, sequence_id, residue_numbers)
        chain_dict = {
            'chain_label': chain.get_id(),
            'residues': residues
        }
        sequence_chains.append(chain_dict)

    predictions_dict['data_resource'] = 'dynamine'
    predictions_dict['resource_entry_url'] = 'http://dynamine.ibsquare.be/'
    predictions_dict['resource_version'] = '3.0.5'
    predictions_dict['pdb_id'] = '${entryPredictions.simpleName}'
    predictions_dict['chains'] = sequence_chains
    predictions_dict['evidence_code_ontology'] = [
        { 'eco_code': 'ECO_0000203' },
        { 'eco_code': 'ECO_0000364' }
    ]
    predictions_dict['sites'] = [
        {
            'site_id': 1,
            'label': 'backbone',
            'source_accession': '${entryPredictions.simpleName}',
            'source_database': 'PDB',
            'source_release_date': release_date.strftime('%d/%m/%Y'),
            'site_url': 'https://bio2byte.be/b2btools/dynamine'
        },
        {
            'site_id': 2,
            'label': 'sidechain',
            'source_accession': '${entryPredictions.simpleName}',
            'source_database': 'PDB',
            'source_release_date': release_date.strftime('%d/%m/%Y'),
            'site_url': 'https://bio2byte.be/b2btools/dynamine'
        },
        {
            'site_id': 3,
            'label': 'efoldmine',
            'source_accession': '${entryPredictions.simpleName}',
            'source_database': 'PDB',
            'source_release_date': release_date.strftime('%d/%m/%Y'),
            'site_url': 'https://bio2byte.be/b2btools/efoldmine'
        }
    ]
    predictions_dict['software_version'] = '1.0'

    with open('${entryPredictions.simpleName}.funpdbe.json', 'w') as output_file:
        json.dump(predictions_dict, output_file, indent=2)
    """
}

process validateFunPdbePrediction {
    publishDir "${projectDir}/results/${entryFunPdbePredictions.simpleName[1]}${entryFunPdbePredictions.simpleName[2]}/", mode: 'copy'
    errorStrategy 'ignore'

    tag "${entryFunPdbePredictions}"

    input:
        path entryFunPdbePredictions

    output:
        path "${entryFunPdbePredictions.simpleName}.valid.funpdbe.json", optional: true, emit: validFunPdbe
        path "${entryFunPdbePredictions.simpleName}.invalid.funpdbe.txt", optional: true, emit: errorsFunPdbe

    script:
    """
    #!/usr/local/bin/python
    import sys
    import os

    from funpdbe_validator.validator import Validator
    from funpdbe_validator.residue_index import ResidueIndexes

    validator = Validator("dynamine")
    validator.load_schema()
    json_file_path = "${entryFunPdbePredictions}"
    validator.load_json(json_file_path)

    errors = []

    if validator.basic_checks() and validator.validate_against_schema():
        print("Passed data validations")
        residue_indexes = ResidueIndexes(validator.json_data)
        if residue_indexes.check_every_residue():
            print("Passed the index validation")
            os.symlink("${entryFunPdbePredictions}", "${entryFunPdbePredictions.simpleName}.valid.funpdbe.json")
            exit(0)
        else:
            error_message = "Failed index validation for %s: %s \\n" % (json_file_path, residue_indexes.mismatches)
            print(error_message)
            errors.append({ 'type': 'residues', 'message': error_message })

    error_message = "Failed data validations for %s: %s \\n" % (json_file_path, validator.error_log)
    print(error_message)
    errors.append({ 'type': 'basic_checks', 'message': error_message })

    with open("${entryFunPdbePredictions.simpleName}.invalid.funpdbe.txt", "w") as file_handler:
        file_handler.writelines(error_message)

    exit(1)
    """
}

process reportStatsSuccess {
    publishDir "${projectDir}/results", mode: 'copy'

    input:
        val validFunPdbeCount
        val failedFunPdbeCount

    output:
        path "reportSuccess.txt"

    script:
        """
        echo -e "Valid JSON in FunPDBe format: $validFunPdbeCount" > reportSuccess.txt
        echo -e "Failed JSON in FunPDBe format: $failedFunPdbeCount" >> reportSuccess.txt
        """
}

workflow {
    extractPdb(params.targetPdbIds)

    allPdbIds = extractPdb.out.flatPdbCsv.splitCsv(header: false, strip: true).flatten().unique()
    fetchPdb(allPdbIds)

    buildFastaFromPdb(fetchPdb.out.pdbTuple)
    predictBiophysicalFeatures(buildFastaFromPdb.out.entrySequencesTuple)
    formatFunPdbePrediction(predictBiophysicalFeatures.out.entryPredictionsTuple)
    validateFunPdbePrediction(formatFunPdbePrediction.out.entryFunPdbePredictions)

    // allValidFunPdbeFiles = validateFunPdbePrediction.out.validFunPdbe.collect()
    // allFailedFunPdbeFiles = validateFunPdbePrediction.out.errorsFunPdbe.collect()

    // reportStatsSuccess(
    //     allValidFunPdbeFiles.count(),
    //     allFailedFunPdbeFiles.count(),
    // )

    validateFunPdbePrediction.out.validFunPdbe
        .map { it.simpleName }
        .collectFile(name: "${projectDir}/results/allValidFunPdbeFiles.txt", newLine: true)
        .subscribe {
            println "Valid entries are saved to file: $it"
        }

    validateFunPdbePrediction.out.errorsFunPdbe
        .map { it.simpleName }
        .collectFile(name: "${projectDir}/results/allFailedFunPdbeFiles.txt", newLine: true)
        .subscribe {
            println "Failed entries are saved to file: $it"
        }
}
