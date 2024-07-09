#' @importFrom basilisk BasiliskEnvironment
env <- BasiliskEnvironment(
    envname = "igblast_wrap_basilisk",
    pkgname = "scifer",
    packages = c(
        "python==3.9.19",
        "dnaio==1.2.1",
        "igblast==1.22.0",
        "cffi==1.16.0",
        "isal==1.6.1"
    ),
    channels = c("bioconda", "conda-forge")
)
