process cdHitClustering{
    tag "${sequencesValid.name}"

    input:
    path sequencesValid
    val clustering

    output:
    path "*_clustered.clstr", emit: clusters

    script:

    //thresholds = -n 8,9,10 for thresholds 0.90 ~ 1.0
    //-n 7      for thresholds 0.88 ~ 0.9
    //-n 6      for thresholds 0.85 ~ 0.88
    //-n 5      for thresholds 0.80 ~ 0.85
    //-n 4      for thresholds 0.75 ~ 0.8
    """
    if (($clustering >=0.9)); then
        n=10
    elif (($clustering >=0.88)); then
        n=7
    elif (($clustering >=0.85)); then
        n=6
    elif (($clustering >=0.8)); then
        n=5
    else
        n=4
    fi
    cd-hit-est -i ${sequencesValid} -o ${sequencesValid.baseName}_${clustering}_clustered  -c $clustering -n \$n -d 40  
    """

}

process postClusteringLabels{
    tag "${sequencesValid.name}"

    input:
    path sequencesValid
    path clusters
    val relabel

    output:
    path "*.fasta" , emit: repSeqs
    path "*_representative_represented.csv" , emit: representativeRepresented

    script:

    """
    python $projectDir/bin/postCdHitRelabel.py $sequencesValid $clusters $relabel
    """


}