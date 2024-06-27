#' @import basilisk


env <- BasiliskEnvironment(envname="igblast_wrap_basilisk",
                              pkgname = "scifer",
                              packages=c("python=3.9", "dnaio==1.2.1", "igblast==1.22.0")
)
