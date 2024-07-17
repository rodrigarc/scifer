#' @importFrom basilisk BasiliskEnvironment
env <- BasiliskEnvironment(
    envname = "igblast_wrap_basilisk",
    pkgname = "scifer",
    packages = c(
        "python==3.9.19",
        "igblast==1.22.0",
        "cffi==1.16.0",
        "python-isal==1.6.1",
        "pip==24.0",
        "pycparser==2.22",
        "setuptools==70.1.1",
        "wheel==0.43.0",
        "xopen==2.0.1",
        "python-zlib-ng==0.4.3",
        "zstandard==0.22.0"
    ),
    channels = c("bioconda", "conda-forge"),
    pip = c("dnaio==1.2.1")
)
