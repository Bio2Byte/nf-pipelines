process compressDirectory {
    publishDir "$projectDir/results", mode: "copy"

    input:
    val compressedFileName
    path outputFiles

    output:
    path "${compressedFileName}.tar.gz"

    script:
    """
    tar -czh -f ${compressedFileName}.tar.gz ${outputFiles}
    """
}