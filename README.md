# NextFlow Pipelines (nf-pipelines)

<p align="center">
  Welcome to the repository of NextFlow pipelines created by the Bio2Byte group!
</p>
<p align="center" style="display: flex; justify-content: center; align-items: center">
  <a href="https://bio2byte.be"><img src="https://pbs.twimg.com/profile_images/1247824923546079232/B9b_Yg7n_400x400.jpg" width="105px"/></a>
  <a href="https://vub.be"><img src="https://upload.wikimedia.org/wikipedia/commons/thumb/3/31/Vrije_Universiteit_Brussel_logo.svg/1200px-Vrije_Universiteit_Brussel_logo.svg.png" alt="vub" width="345px"/></a>
</p>

## About bio2Byte

Proteins are the molecular machines that make cells work. They perform a wide variety of functions through interactions with each other and many additional molecules. Traditionally, proteins are described in a single static state (a picture). It is now increasingly recognized that many proteins can adopt multiple states and move between these conformational states dynamically (a movie).

We investigate how the dynamics, conformational states and available experimental data of proteins relates to their amino acid sequence. Underlying physical and chemical principles are computationally unravelled through data integration, analysis and machine learning, so connecting them to biological events and improving our understanding of the way proteins work.

The Bio2Byte group is primarily situated at the Interuniversity Institute of Bioinformatics in Brussels, a collaborative inter-faculty institute between the "Vrije Universiteit Brussel" (VUB) and the "Université Libre de Bruxelles" (ULB). It is located at the the ULB side of the "Pleinlaan/La Plaine" campus on the 6th floor of the C building.

At the VUB, the group is linked to Structural Biology Brussels at the Bioengineering sciences department, as well as the departments of Computer Science and Chemistry at the Faculty of Sciences, plus to the Biomedical sciences at the Faculty of Medicine and Pharmacy.

### About the biophysical predictors (b2btools)

The `b2btools` are a set of in-house developed biophysical predictors that enable the exploration of the 'biophysical variation' of proteins along its sequence. The predictions reflect 'emerging' properties, so what the sequence is capable of, not necessarily what it will do in a particular context, for example when it adopts a specific fold. Studying the biophysical properties of a protein is relevant as these properties, like the dynamics of a protein, are conserved by evolution in order to preserve the protein's function.

- Backbone and sidechain dynamics (**DynaMine**)
- Conformational propensities (sheet, helix, coil, polyproline II) (**DynaMine**)
- Early folding propensities (**EFoldMine**)
- Disorder propensities (**DisoMine**)
- Beta aggregation propensity (**AgMata**)
- Phase separation propensity (**PPser**)

<hr>

## Get started

> NextFlow is a reactive workflow framework and a programming DSL that eases writing computational pipelines with complex data.

NextFlow must be available on the system where the pipeline is going to be executed. For local environments, please follow the instructions provided by the official docs: [NextFlow - Get started](https://www.nextflow.io/docs/latest/getstarted.html).

If you are working on a HPC environment with software dependencies handled as modules, for instance VSC clusters, NextFlow is available by executing:

```shell
$ module load Nextflow/22.04.0
```

To check that everything is ready:

```shell
$ nextflow -v
nextflow version 22.04.0.5697
```

### Basic concepts

The basic command to run our [NextFlow pipelines](https://www.nextflow.io/docs/latest/index.html) is:

```shell
$ nextflow run Bio2Byte/nf-pipelines -r main -main-script /path/to/script.nf
```

Given NextFlow able to run pipelines from GitHub, the published ones here are available by executing "`nextflow run Bio2Byte/nf-pipelines`" with the flag "`-r main`" to indicate NextFlow to fetch the `main` branch. By adding the last flag "`-main-script`" you choose the pipeline.

#### Executors

NextFlow is a powerful pipeline framework that runs on different systems such as your local workspace or HPC clusters. To aim these different [NextFlow executors](https://www.nextflow.io/docs/latest/executor.html), we defined profiles inside the configuration files:

- For local development: "`-profile standard,...`"
- For [Slurm](https://slurm.schedmd.com/documentation.html)-based HPC clusters (such as [VUB's Hydra clusters](https://hpc.vub.be/docs/infrastructure/#hydra)): "`-profile hydra,...`"

#### Containerization

In addition, NextFlow can run processes using either [Docker](https://docs.docker.com/) or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html) images. For instance, VUB-HPC prefers Singularity to Docker when running containerized software ([read more](https://docs.vscentrum.be/en/latest/software/singularity.html?highlight=docker#can-i-run-containers-on-the-hpc-systems)). This feature is a game changer because you do not have to install any dependency to run our code. Each step of the pipeline (the "`process`" blocks) has a Docker image defined inside the configuration file. In case that using Singularity is a requirement, NextFlow will convert automatically Docker images to Singularity ones.

- Using Docker: "`-profile withdocker,...`"
- Using Singularity: "`-profile withsingularity,...`"

### In summary

If you are working on your local environment combine both `standard` and `withdocker` profiles:

```shell
$ nextflow run Bio2Byte/nf-pipelines \
  -r main \
  -main-script /path/to/script.nf \
  -profile standard,withdocker \
  [pipeline parameters]
```

If you are working on HPC clusters combine both `hydra` and `withsingularity` profiles:

```shell
$ nextflow run Bio2Byte/nf-pipelines \
  -r main \
  -main-script /path/to/script.nf \
  -profile standard,withsingularity \
  [pipeline parameters]
```

## Published Pipelines

### Biophysical features predicted from sequence only

This pipeline provides you structural predictions for protein sequences in FASTA format.

- Input:
  - At least one valid protein sequence file in FASTA format.
- Output:
  - Compressed file in format "`.tar.gz`" that includes:
    - Biophysical predictions in paginated JSON files (depending on the "`--sequencesPerJson`" parameter)
    - Ignored sequences in FASTA format of proteins with less than 5 or more than 2000 residues.
    - (optional) Structures in PDB format from ESM Atlas API (only for first 400 residues of each sequence)
    - (optional) MatPlotLib plots in PNG format:
      - "`<Sequence's name>.png`" contains plots of:
        - Backbone and sidechain dynamics (**DynaMine**)
        - Conformational propensities (sheet, helix, coil, polyproline II) (**DynaMine**)
        - Early folding propensities (**EFoldMine**)
        - Disorder propensities (**DisoMine**)
        - Beta aggregation propensity (**AgMata**)
      - "`<Sequence's name>_psp.png`" contains plots of:
        - Phase separation propensity (**PPser**)

#### Quick example

For these sequences in FASTA format:

```fasta
>SEQ_BAD
MA
>SEQ_ADR
MAAAALVLVLLAVLAMAAAALVLVLLAVLAMAAAALVLVLLAVLAMAAAALVLVLLAVLA
>SEQ_SOP
MALALSALALLALALAMAAAALVLVLLAVLAMAAAALVLVLLAVLAMAAAALVLVLLAVLAMAAAALVLVLLAVLAMAAAALVLVLLAVLA
```

The NextFlow command line example is:

```shell
$ nextflow run fromSingleSequences.nf \
    -profile standard,withdocker \
    --targetSequences ../example.fasta \
    --dynamine \
    --efoldmine \
    --disomine \
    --agmata \
    --psper \
    --fetchStructures \
    --plotBiophysicalFeatures
```

This pipeline will generate a compressed file containing:

```shell
example_2023_02_08_11_13_22.tar.gz
├── SEQ_ADR.pdb
├── SEQ_SOP.pdb
├── b2b_results_example.1.index
├── b2b_results_example.1.json
├── b2b_results_example.1_SEQ_ADR.png
├── b2b_results_example.1_SEQ_ADR_psp.png
├── b2b_results_example.1_SEQ_SOP.png
├── b2b_results_example.1_SEQ_SOP_psp.png
└── example_sequences_ignored.fasta
```

#### Parameters

- Input/Output:
  - Sequences file in FASTA format: `--targetSequences`
    - Description: Sequence must contain at least one protein sequence.
  - Sequences per results file in JSON format: `--sequencesPerJson`
    - Description: Biophysical predictions are saved in JSON format. By default each file contains up to 10 sequences.
- Predictions:
  - Dynamine - Backbone dynamics: `--dynamine`
  - EFoldMine - Early folding regions: `--efoldmine`
  - Disomine - Disordered regions: `--disomine`
  - AgMata - Regions prone to Beta-aggregation: `--agmata`
  - PSPer - Phase Separating Protein: `--psper`
  - ESM Atlas structures: `--fetchStructures`
    - Description: Resource provided by the "ESM Metagenomic Atlas". Fetching the folding of a sequence runs ESMFold (`esm.pretrained.esmfold_v1`) which provides fast and accurate atomic level structure prediction directly from the individual sequence of a protein.
- Visualizations:
  - Biophysical features charts: `--plotBiophysicalFeatures`
    - Description: Include `MatPlotLib` plots in PNG format.

#### Diagram (Mermaid format)

[![](https://mermaid.ink/img/pako:eNp9VEtzmzAQ_isenciM40HCgMShM20evXU6dU7BOSggB80gRJGcxInz3yvAEpadhhO737fvXb2DQpYMZGBTy5eiop2e3V2vm5n52jAIriraNKxebDopflNdXVwcMBjkqq25vqVK0werRUH-2NGmqJwmCvJnzl6cvAyCmXMSn4DJMZgGuaCtw3CQF7KuWaFvec2cmnyuhuGn-UH4Hzr6nB4dJwS93GHsSV7qMPUknLcdK3mhf3DZVjvFC1rfMqq3HVMPlkSOTVDoSdCTkCd5KaJl3tbyqzgo9tuKknzDdFHdKPFd11StdLctegtnkE5bwESrd1MwMxLBXydX0yym8XuFRNAQhGmGUte8MzzZ7WycyCsrnF1eftvTul6xv1vWFEztTY9sr3rQxLNxRzH2xcj6HcWlXblRTHxyanduCKtszBVtuOZvrDTBsd1Dn2LWSLNuYMDQd0rsivoWf5iQz4MBwm5bRwvoKjyUCJFbUN_Jz05u2zGqczJUui93DRW8YR60HCC2kXV5hsWjGVfyDErGKTwJqqkHjH1qVcs6T48Pk1meKhyDnCbpuCg8S3LC4GmSE4T8JCcg8pKc9GMz-jtRx0P4YvbI7dZhfZDbn7FFyl6N5xGlp9WftgNZBXE3NPjbmMdJ3cmrw7Hs-9NxFzUw7B2xsn_IflHBPNJhfyIE5kCwTlBemjf-vYfXQFdMsDXIzO8jVeZv3XwY3rYtqWY3JTdnCTJTEZsDutVytWsKK4-ca06fOipAtqG1MtqWNiB7B68ggylekDSBOMRpgkiUzsEOZJcQLiAmOI4TvMRREsGPOXiT0ngIFyQiSUhSgkKCCYGDt_sB60N-_ANEmOkf?type=png)](https://mermaid.live/edit#pako:eNp9VEtzmzAQ_isenciM40HCgMShM20evXU6dU7BOSggB80gRJGcxInz3yvAEpadhhO737fvXb2DQpYMZGBTy5eiop2e3V2vm5n52jAIriraNKxebDopflNdXVwcMBjkqq25vqVK0werRUH-2NGmqJwmCvJnzl6cvAyCmXMSn4DJMZgGuaCtw3CQF7KuWaFvec2cmnyuhuGn-UH4Hzr6nB4dJwS93GHsSV7qMPUknLcdK3mhf3DZVjvFC1rfMqq3HVMPlkSOTVDoSdCTkCd5KaJl3tbyqzgo9tuKknzDdFHdKPFd11StdLctegtnkE5bwESrd1MwMxLBXydX0yym8XuFRNAQhGmGUte8MzzZ7WycyCsrnF1eftvTul6xv1vWFEztTY9sr3rQxLNxRzH2xcj6HcWlXblRTHxyanduCKtszBVtuOZvrDTBsd1Dn2LWSLNuYMDQd0rsivoWf5iQz4MBwm5bRwvoKjyUCJFbUN_Jz05u2zGqczJUui93DRW8YR60HCC2kXV5hsWjGVfyDErGKTwJqqkHjH1qVcs6T48Pk1meKhyDnCbpuCg8S3LC4GmSE4T8JCcg8pKc9GMz-jtRx0P4YvbI7dZhfZDbn7FFyl6N5xGlp9WftgNZBXE3NPjbmMdJ3cmrw7Hs-9NxFzUw7B2xsn_IflHBPNJhfyIE5kCwTlBemjf-vYfXQFdMsDXIzO8jVeZv3XwY3rYtqWY3JTdnCTJTEZsDutVytWsKK4-ca06fOipAtqG1MtqWNiB7B68ggylekDSBOMRpgkiUzsEOZJcQLiAmOI4TvMRREsGPOXiT0ngIFyQiSUhSgkKCCYGDt_sB60N-_ANEmOkf)

### Biophysical features predicted from Multiple Sequence Alignment (MSA)

This pipeline provides you predictions of the biophysical features of a MSA, computation of the MSA biophysical and sequence conservation, and 2D visualization of these values.

- Input:
  - Either:
    - At Multiple Sequence Alignment file in FASTA format
    - At least **3** valid protein sequences in a single file in FASTA format
- Output:
  - Compressed file in format "`.tar.gz`" that includes:
    - MSA file in FASTA format.
    - Biophysical predictions in JSON format
    - Ignored sequences in FASTA format of proteins with less than 5 or more than 2000 residues.
    - (optional) Structures in PDB format from ESM Atlas API (only for first 400 residues of each sequence)
    - (optional) MatPlotLib plots in both PNG and PDF formats
    - (optional) WebLogo representation of the MSA
    - (optional) Phylogenetic tree in text format
    - (optional) Phylogenetic tree rendered in SVG image format

#### Quick examples

Given this pipelines supports two types of input files, the next section will describe both scenarios:

- Running the pipeline for a MSA input file in FASTA format
- Running the pipeline for a set of sequences in FASTA format

##### Running the pipeline for a MSA input file

For this MSA in FASTA format:

```fasta
>1ymg_A  THE CHANNEL ARCHITE
SASFWRAICAEFFASLFYVFFGLGASLRW-----AG------P---------lHVLQVAL
AFGLALATLVQAVGHISGAHVNPAVTFAFLVGSQMSLLRAICYMVAQLLGAVAGAAVLYS
VT--PPAvRGNlALNTLHPGVSVGQATIVEIFLTLQFVLCIFATYDERRNGRLGSVALAV
GFSLTLGHLFGMYYTGAGMNPARSFAPAILTR------NFTNHWVYWVGPVIGAGLGSLL
YDFLLFPRLKSVSERLSILKG

>2d57_A  DOUBLE LAYERED 2D C
TQAFWKAVTAEFLAMLIFVLLSVGSTINW-----GG-SENPLP---------VDMVLISL
CFGLSIATMVQCFGHISGGHINPAVTVAMVCTRKISIAKSVFYITAQCLGAIIGAGILYL
VT--PPSVVGGLGVTTVHGNLTAGHGLLVELIITFQLVFTIFASCDSKRTDVTGSVALAI
GFSVAIGHLFAINYTGASMNPARSFGPAVIMG------NWENHWIYwVGPIIGAVLAGAL
YEYVF--------------CP

>2f2b_A  CRYSTAL STRUCTURE O
MVSLTKRCIAEFIGTFILVFFGAGSAAVTLMIASGGTSPNPFNIGIGLLGGLGDWVAIGL
AFGFAIAASIYALGNISGCHINPAVTIGLWSVKKFPGREVVPYIIAQLLGAAFGSFIFLQ
CAGIGAATVGGLGATAPFPGISYWQAMLAEVVGTFLLMITIMGIAvDERAP-KGFAGIII
GLTVAGIITTLGNISGSSLNPARTFGPYLNDMifagtDlWNYYSIYvIGPIVGAVLAALT
YQYL---------------TS
```

The NextFlow command line example is:

```shell
$ nextflow run fromMSA.nf \
    -profile standard,withdocker \
    --targetSequences ./simple_alignment.fasta \
    --plotBiophysicalFeatures \
    --buildLogo \
    --buildTree \
    --plotTree \
    --efoldmine \
    --disomine \
    --fetchStructures
```

This pipeline will generate a compressed file containing:

```shell
simple_alignment_2023_02_08_12_57_28.tar.gz
├── 1ymg_A__THE_CHANNEL_ARCHITE.pdb
├── 2d57_A__DOUBLE_LAYERED_2D_C.pdb
├── 2f2b_A__CRYSTAL_STRUCTURE_O.pdb
├── b2b_msa_results_simple_alignment_filtered.fasta.json
├── simple_alignment_filtered.fasta.msa
├── simple_alignment_filtered.fasta.msa.tree
├── simple_alignment_filtered.fasta.msa.tree.svg
├── simple_alignment_filtered.fasta_1ymg_A__THE_CHANNEL_ARCHITE_msa_biophysical_conservation.pdf
├── simple_alignment_filtered.fasta_1ymg_A__THE_CHANNEL_ARCHITE_msa_biophysical_conservation.png
├── simple_alignment_filtered.fasta_2d57_A__DOUBLE_LAYERED_2D_C_msa_biophysical_conservation.pdf
├── simple_alignment_filtered.fasta_2d57_A__DOUBLE_LAYERED_2D_C_msa_biophysical_conservation.png
├── simple_alignment_filtered.fasta_2f2b_A__CRYSTAL_STRUCTURE_O_msa_biophysical_conservation.pdf
├── simple_alignment_filtered.fasta_2f2b_A__CRYSTAL_STRUCTURE_O_msa_biophysical_conservation.png
└── simple_alignment_filtered_logo.png
```

##### Running the pipeline for a set of sequences in FASTA format

For these sequences in FASTA format:

```fasta
>random_sequence_1 consisting of 25 residues.
RGGMSIQGTFVR

>random_sequence_2 consisting of 25 residues.
CE

>random_sequence_3 consisting of 25 residues.
LKLPSFD

>random_sequence_4 consisting of 25 residues.
PKMCQMTDHKEYQGSALGSGS

>random_sequence_5 consisting of 25 residues.
REDTWATASAACLITFNVSPDCMQV

>random_sequence_6 consisting of 25 residues.
IHYTTTPDVICLW

>random_sequence_7 consisting of 25 residues.
AAMQLCPQAYMKQPWTNVMIQE

>random_sequence_8 consisting of 25 residues.
FSSFHMHTMHLLSPLNTD

>random_sequence_9 consisting of 25 residues.
IIKYGAV

>random_sequence_10 consisting of 25 residues.
FPDHDDCGCLWFYQATWATKCLKEL
```

The NextFlow command line example is:

```shell
$ nextflow run fromMSA.nf \
    -profile standard,withdocker \
    --targetSequences ./10x25.fasta \
    --plotBiophysicalFeatures \
    --alignSequences \
    --buildLogo \
    --buildTree \
    --plotTree \
    --efoldmine \
    --disomine \
    --fetchStructures
```

This pipeline will generate a compressed file containing:

```shell
10x25_2023_02_08_13_20_44.tar.gz
├── 10x25_filtered.fasta.msa
├── 10x25_filtered.fasta.msa.tree
├── 10x25_filtered.fasta.msa.tree.svg
├── 10x25_filtered.fasta_random_sequence_10_consisting_of_25_residues__msa_biophysical_conservation.pdf
├── 10x25_filtered.fasta_random_sequence_10_consisting_of_25_residues__msa_biophysical_conservation.png
├── 10x25_filtered.fasta_random_sequence_1_consisting_of_25_residues__msa_biophysical_conservation.pdf
├── 10x25_filtered.fasta_random_sequence_1_consisting_of_25_residues__msa_biophysical_conservation.png
├── 10x25_filtered.fasta_random_sequence_3_consisting_of_25_residues__msa_biophysical_conservation.pdf
├── 10x25_filtered.fasta_random_sequence_3_consisting_of_25_residues__msa_biophysical_conservation.png
├── 10x25_filtered.fasta_random_sequence_4_consisting_of_25_residues__msa_biophysical_conservation.pdf
├── 10x25_filtered.fasta_random_sequence_4_consisting_of_25_residues__msa_biophysical_conservation.png
├── 10x25_filtered.fasta_random_sequence_5_consisting_of_25_residues__msa_biophysical_conservation.pdf
├── 10x25_filtered.fasta_random_sequence_5_consisting_of_25_residues__msa_biophysical_conservation.png
├── 10x25_filtered.fasta_random_sequence_6_consisting_of_25_residues__msa_biophysical_conservation.pdf
├── 10x25_filtered.fasta_random_sequence_6_consisting_of_25_residues__msa_biophysical_conservation.png
├── 10x25_filtered.fasta_random_sequence_7_consisting_of_25_residues__msa_biophysical_conservation.pdf
├── 10x25_filtered.fasta_random_sequence_7_consisting_of_25_residues__msa_biophysical_conservation.png
├── 10x25_filtered.fasta_random_sequence_8_consisting_of_25_residues__msa_biophysical_conservation.pdf
├── 10x25_filtered.fasta_random_sequence_8_consisting_of_25_residues__msa_biophysical_conservation.png
├── 10x25_filtered.fasta_random_sequence_9_consisting_of_25_residues__msa_biophysical_conservation.pdf
├── 10x25_filtered.fasta_random_sequence_9_consisting_of_25_residues__msa_biophysical_conservation.png
├── 10x25_filtered_logo.png
├── 10x25_sequences_ignored.fasta
├── b2b_msa_results_10x25_filtered.fasta.json
├── random_sequence_10_consisting_of_25_residues_.pdb
├── random_sequence_1_consisting_of_25_residues_.pdb
├── random_sequence_3_consisting_of_25_residues_.pdb
├── random_sequence_4_consisting_of_25_residues_.pdb
├── random_sequence_5_consisting_of_25_residues_.pdb
├── random_sequence_6_consisting_of_25_residues_.pdb
├── random_sequence_7_consisting_of_25_residues_.pdb
├── random_sequence_8_consisting_of_25_residues_.pdb
└── random_sequence_9_consisting_of_25_residues_.pdb
```

#### Parameters

- Input/Output:
  - Sequences file in FASTA format: `--targetSequences`
    - Description: Sequence must contain at least one protein sequence.
  - Should pipeline align sequences? `--alignSequences`
    - Description: If your input file is a set of unaligned sequences, with this flag the pipeline will build a MSA file using [Clustal Omaga](https://www.ebi.ac.uk/Tools/msa/clustalo/) before starting the regular workflow. Clustal Omega is a new multiple sequence alignment program that uses seeded guide trees and HMM profile-profile techniques to generate alignments between three or more sequences.
- Predictions:
  - EFoldMine - Early folding regions: `--efoldmine`
  - Disomine - Disordered regions: `--disomine`
  - ESM Atlas structures: `--fetchStructures`
    - Description: Resource provided by the "ESM Metagenomic Atlas". Fetching the folding of a sequence runs ESMFold (`esm.pretrained.esmfold_v1`) which provides fast and accurate atomic level structure prediction directly from the individual sequence of a protein.
  - MSA Logo: `--buildLogo`
    - Description: Create a Logo representation using [WebLogo](https://weblogo.berkeley.edu/). Sequence logos are a graphical representation of an amino acid or nucleic acid multiple sequence alignment developed by Tom Schneider and Mike Stephens.
  - Phylogenetic tree: `--buildTree`
    - Description: [FastTree](http://www.microbesonline.org/fasttree/) infers approximately-maximum-likelihood phylogenetic trees from alignments of protein sequences.
  - Phylogenetic tree image in SVG format: `--plotTree`
    - Description: Build a SVG image of the phylogenetic tree. It requires the flag `--buildTree`.
- Visualizations:
  - Biophysical features charts: `--plotBiophysicalFeatures`
    - Description: Include `MatPlotLib` plots in PNG format.

#### Diagram (Mermaid format)
[![](https://mermaid.ink/img/pako:eNqVVU1v2zAM_SuBTi6QBpbt-OswoGuaYcBaFEtPc3pQbSYWJlmeJbdN0_73KXalxv3a4pNFvkeK1KO0RbkoAKVoxcRdXpJGja5my2qkv9p1nNOSVBWwyaoR_JKo8ujo2YedTNaMqjmRilwbq-dkNw2p8tJafCe7pXBn14HjjGyQ6StnuO-MnIyT2vpiJ8sFY5CrOWVgzcn7ZuxmivyG85YpWjNYwJ8WqhxOGF1XHCp1bXA4u2kpKy7LDRNrqEDR_KoBsH4vq5lQH7r9nv5DrIW1BVndQEFz9ZWKutxImhM2B6LaBqQFTfdLxYPCcdTl_IwdD3uDk2wFKi_PJD9RjMiFatp8xzAEb-8ogddqY5N5-iA5vX85Qc829MXm72_PCzSA6xKlnNFG40SzsXkGZbmj4-Mvj4Qx0375qPdq9rxzaoZh9svpcOkbFfXLwOimX4ZDcGSE06WVJueCVFTRByh08tiIaQjR0lHQdAjsDoMmRmdDxk_g4rYjeLagvlz-keR2wS22q_6xfiWrHcSzyushWgkK3ujz0MS-FWyH1bHEoSECq-_nkzuEHFndd1hYCVZwWsHAF3a-gkrxxhXtN-OdyfheXc7mg4r-g3HxbZ_xiWxwbAevrx0ndvJ6lhk4uR_Rcw_vlPdvibzENzNkh6jXzErfhPJKnD5P6Y4S2FHuEGaAodjdmheEwwBkTniKxohDwwkt9Aux3bmXSJXAYYlS_XtDpP5bVk8a19YFUXBWUH0foFT3A8aItEosNlWO0hVhEgxoRsm6Idxaa1KhdIvuUYqjeJJEIY7dOAq9xI_GaIPSY4wnOE7i6TSMg9gPffw0Rg9C6AjuJPGT0NWc0HXDIEziLtyvztlvArotnfevXPfYPf0FFM4oUg?type=png)](https://mermaid.live/edit#pako:eNqVVU1v2zAM_SuBTi6QBpbt-OswoGuaYcBaFEtPc3pQbSYWJlmeJbdN0_73KXalxv3a4pNFvkeK1KO0RbkoAKVoxcRdXpJGja5my2qkv9p1nNOSVBWwyaoR_JKo8ujo2YedTNaMqjmRilwbq-dkNw2p8tJafCe7pXBn14HjjGyQ6StnuO-MnIyT2vpiJ8sFY5CrOWVgzcn7ZuxmivyG85YpWjNYwJ8WqhxOGF1XHCp1bXA4u2kpKy7LDRNrqEDR_KoBsH4vq5lQH7r9nv5DrIW1BVndQEFz9ZWKutxImhM2B6LaBqQFTfdLxYPCcdTl_IwdD3uDk2wFKi_PJD9RjMiFatp8xzAEb-8ogddqY5N5-iA5vX85Qc829MXm72_PCzSA6xKlnNFG40SzsXkGZbmj4-Mvj4Qx0375qPdq9rxzaoZh9svpcOkbFfXLwOimX4ZDcGSE06WVJueCVFTRByh08tiIaQjR0lHQdAjsDoMmRmdDxk_g4rYjeLagvlz-keR2wS22q_6xfiWrHcSzyushWgkK3ujz0MS-FWyH1bHEoSECq-_nkzuEHFndd1hYCVZwWsHAF3a-gkrxxhXtN-OdyfheXc7mg4r-g3HxbZ_xiWxwbAevrx0ndvJ6lhk4uR_Rcw_vlPdvibzENzNkh6jXzErfhPJKnD5P6Y4S2FHuEGaAodjdmheEwwBkTniKxohDwwkt9Aux3bmXSJXAYYlS_XtDpP5bVk8a19YFUXBWUH0foFT3A8aItEosNlWO0hVhEgxoRsm6Idxaa1KhdIvuUYqjeJJEIY7dOAq9xI_GaIPSY4wnOE7i6TSMg9gPffw0Rg9C6AjuJPGT0NWc0HXDIEziLtyvztlvArotnfevXPfYPf0FFM4oUg)

<hr>

## Links

- [Bio2Byte](https://bio2byte.be)
- [Bio2Byte online predictors](https://bio2byte.be/b2btools)
- [Feedback or Questions](https://www.bio2byte.be/b2btools/feedback)
- [Bio2Byte tools package](https://pypi.org/project/b2bTools/)
- [Bio2Byte Online Notebooks](https://github.com/Bio2Byte/public_notebooks)

## Citations

**Implementation of the b2btools to study the protein biophysical features and their conservation**

> Kagami, L. P., Orlando, G., Raimondi, D., Ancien, F., Dixit, B., Gavaldá-García, J., Ramasamy, P., Roca-Martínez, J., Tzavella, K., & Vranken, W. (2021). b2bTools: Online predictions for protein biophysical features and their conservation. Nucleic Acids Research, 49(W1), W52–W59. https://doi.org/10.1093/nar/gkab425

**DynaMine**

> Cilia, E., Pancsa, R., Tompa, P., Lenaerts, T., & Vranken, W. F. (2013). From protein sequence to dynamics and disorder with DynaMine. Nature Communications, 4(1), 2741–2741. https://doi.org/10.1038/ncomms3741

**EFoldMine**

> Raimondi, D., Orlando, G., Pancsa, R., Khan, T., & Vranken, W. F. (2017). Exploring the Sequence-based Prediction of Folding Initiation Sites in Proteins. Scientific Reports, 7(1), 8826–8826. https://doi.org/10.1038/s41598-017-08366-3

**Disomine**

> Orlando, G., Raimondi, D., Codicè, F., Tabaro, F., & Vranken, W. (2022). Prediction of Disordered Regions in Proteins with Recurrent Neural Networks and Protein Dynamics. Journal of Molecular Biology, 434(12), 167579. https://doi.org/10.1016/j.jmb.2022.167579

**AgMata**
> Orlando, G., Silva, A., Macedo-Ribeiro, S., Raimondi, D., & Vranken, W. (2020). Accurate prediction of protein beta-aggregation with generalized statistical potentials. Bioinformatics, 36(7), 2076–2081. https://doi.org/10.1093/bioinformatics/btz912

**PSPer**
> Orlando, G., Raimondi, D., Tabaro, F., Codicè, F., Moreau, Y., & Vranken, W. F. (2019). Computational identification of prion-like RNA-binding proteins that form liquid phase-separated condensates. Bioinformatics, 35(22), 4617–4623. https://doi.org/10.1093/bioinformatics/btz274

<hr>

## Copyright

© Wim Vranken, Bio2Byte group, Vrije Universiteit Brussel (VUB), Brussels, Belgium.