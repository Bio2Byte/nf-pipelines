process compressDirectory {
    publishDir "$projectDir", mode: "copy"
    debug true
    input:
    val compressedFileName
    path outputFiles

    output:
    path "${compressedFileName}.tar.gz"

    script:
    """
    echo ${outputFiles}
    tar -czvhf ${compressedFileName}.tar.gz ${outputFiles}
    """
}