#' @importFrom basilisk BasiliskEnvironment
env <- BasiliskEnvironment(
    envname = "igblast_wrap_basilisk",
    pkgname = "scifer",
    packages = c(
        "python==3.9.19",
        "dnaio==1.2.1",
        "igblast==1.22.0",
        "cffi==1.16.0",
        "python-isal==1.6.1",
        "pip==24.0",
        "pycparser==2.22",
        "setuptools==70.1.1"
    ),
    channels = c("bioconda", "conda-forge")
)
