#' @import basilisk
#' @importFrom here here
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

env_unix_arch64 <- list(
  packages = c(
    "python==3.9.19",
    "igblast==1.22.0",
    "cffi==1.16.0",
    "python-isal==1.6.1",
    "pip==24.0",
    "pycparser==2.22",
    "xopen==2.0.1",
    "setuptools==70.1.1",
    "wheel==0.43.0",
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

if (basilisk::isWindows()) {
  # Not implemented yet, let the user download igblast rather than downloding it
  # Define the URL and destination path for Windows
  #url_igblast <- "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-1.22.0-x64-win64.tar.gz"
  #dest_dir <- system.file("script", package = "scifer")
  #dest_file <- here::here(dest_dir, "ncbi-igblast-1.22.0-x64-win64.tar.gz")
  #print("Running on a Windows machine.")
  # Download the file
  # Create the destination directory if it doesn't exist
  #if (!dir.exists(dest_dir)) {
  #  stop("Destination directory does not exist.")
  #}
  #print(paste0("Downloading the igblast executable to ", dest_file, "..."))
  #download.file(url = url_igblast, destfile = dest_file, mode = "wb")
  # Check if the file has been downloaded correctly
  #if (file.exists(dest_file)) {
  #  print("Igblast downloaded successfully. Unzipping...")
  #  R.utils::gunzip(dest_file, overwrite=TRUE, remove=FALSE)
  #} else {
  #  stop("Failed to download the igblast executable.\n")
  #}
  # set the windows environment
  env <- env_windows
} else if (isMacOSXArm()) {
  # set the OSXArm environment to pass the tests
  # at the momment igblast is not runnable from arm64 architecture
  env <- env_osxArm
} else if(isLinuxAarch64()){
  # set env for linux with arm64 architecture
  env <- env_unix_arch64 
} else if (isLinux() | isMacOSX()) {
  # set env for linux and intel MacOS OS
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