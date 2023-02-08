process compressDirectory {
    publishDir "$launchDir", mode: "copy"

    input:
    val compressedFileName
    path outputFiles

    output:
    path "${compressedFileName}.tar.gz"

    script:
    """
    tar -czvhf ${compressedFileName}.tar.gz ${outputFiles}
    """
}