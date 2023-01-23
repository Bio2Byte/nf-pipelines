# nf-pipelines
In this repository we store our shared Nextflow pipelines.

<p align="center">
  <a href="https://bio2byte.be"><img src="https://pbs.twimg.com/profile_images/1247824923546079232/B9b_Yg7n_400x400.jpg" width="48px"/></a>
  <a href="https://vub.be"><img src="https://upload.wikimedia.org/wikipedia/commons/thumb/3/31/Vrije_Universiteit_Brussel_logo.svg/1200px-Vrije_Universiteit_Brussel_logo.svg.png" alt="vub" width="120px"/></a>
</p>

## About Bio2Byte

Proteins are the molecular machines that make cells work. They perform a wide variety of functions through interactions with each other and many additional molecules. Traditionally, proteins are described in a single static state (a picture). It is now increasingly recognised that many proteins can adopt multiple states and move between these conformational states dynamically (a movie).

We investigate how the dynamics, conformational states and available experimental data of proteins relates to their amino acid sequence. Underlying physical and chemical principles are computationally unravelled through data integration, analysis and machine learning, so connecting them to biological events and improving our understanding of the way proteins work.

The Bio2Byte group is primarily situated at the Interuniversity Institute of Bioinformatics in Brussels, a collaborative interfaculty institute between the Vrije Universiteit Brussel and the Université Libre de Bruxelles. It is located at the the ULB side of the Pleinlaan/La Plaine campus on the 6th floor of the C building.

At the VUB, the group is linked to Structural Biology Brussels at the Bioengineering sciences department, as well as the departments of Computer Science and Chemistry at the Faculty of Sciences, plus to the Biomedical sciences at the Faculty of Medicine and Pharmacy.

## Available pipelines

### bio2Byte - Biophysical features predicted from sequence only

Options:

- single sequence: `--sseq`
- multiple sequence alignment: `--msa`
- predictors:
  - Backbone dynamics: `--dynamine`
  - Early folding regions: `--efoldmine`
  - Disordered regions: `--disomine`
  - Regions prone to Beta-aggregation: `--agamata`

Command example:

```console
nextflow run Bio2Byte/nf-pipelines -r main -main-script b2btools/main.nf
```
