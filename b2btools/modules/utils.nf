process compressDirectory {
    publishDir "$projectDir", mode: "copy"

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