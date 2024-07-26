#' @importFrom basilisk BasiliskEnvironment
#' @importFrom basilisk.utils isMacOSX isMacOSXArm
env_unix <- list(
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
        "zstandard==0.22.0"),
    channels = c("bioconda", "conda-forge"),
    pip = c("dnaio==1.2.1")
)

env_osxArm <- list(
  packages = c(
    "python==3.9.19",
    "cffi==1.16.0",
    "python-isal==1.6.1",
    "pip==24.0",
    "pycparser==2.22",
    "setuptools==70.1.1",
    "wheel==0.43.0",
    "xopen==2.0.1",
    "python-zlib-ng==0.4.3",
    "zstandard==0.22.0"),
  channels = c("bioconda", "conda-forge"),
  pip = c("dnaio==1.2.1")
)

env_windows <- list(
  packages = c(
    "python==3.9.19",
    "cffi==1.16.0",
    "python-isal==1.6.1",
    "pip==24.0",
    "pycparser==2.22",
    "setuptools==70.1.1",
    "wheel==0.43.0",
    "xopen==2.0.1",
    "python-zlib-ng==0.4.3",
    "zstandard==0.22.0"),
  channels = c("bioconda", "conda-forge"),
  pip = c("dnaio==1.2.1")
)

# Switch environment

if (basilisk.utils::isWindows()) {
  env <- env_windows
} else if (basilisk.utils::isMacOSXArm()) {
  env <- env_osxArm 
  } else if (Sys.info()[["sysname"]] == "Linux"|
           basilisk.utils::isMacOSX()) {
  env <- env_unix
} else {
  stop("Unsupported operating system or architecture.\n")
}

env <- BasiliskEnvironment(
  envname = "igblast_wrap_basilisk",
  pkgname = "scifer",
  packages = env$packages,
  channels = env$channels,
  pip = env$pip
)
