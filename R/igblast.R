#' Run IgDiscover for IgBlast using reticulate and conda environment
#'
#' @param database Vector containing the database for VDJ sequences
#' @param fasta Vector containing the sequences, usually a column from a data.frame. eg. df$sequences
#' 
#' @return Creates a data frame with the Igblast analysis where each row is the tested sequence with columns containing the results for each sequence
#' @import reticulate here basilisk
#' @examples
#' ## Example with test sequences
#' igblast(
#'     database = system.file("inst/extdata/test_fasta/KIMDB_rm"),
#'                                 ,
#'     fasta = system.file("inst/extdata/test_fasta/test_igblast.txt")
#' )
#' @export


# library(reticulate)
# library(here)
# library(basilisk)

#check which environments are installed
#conda_list()
# conda env created through the terminal
# conda create -n igblast-wrap python=3.9 dnaio igblast
# cond env created with reticulate
#conda_create(envname = "igblast-wrap", packages = c("dnaio", "igblast"), channel = "bioconda", python_version = 3.9)

igblast <- function(database = "path/to/folder", fasta = "path/to/file", threads = 1) {
    db_dir <- here(database)
    fasta_dir <- here(fasta) 
  
     #db_dir <- here("inst", "extdata", "test_fasta", "KIMDB_rm")
     #fasta_dir <- here("inst", "extdata", "test_fasta")
    
    if (dir.exists(db_dir)) {
      print(paste("The database directory", db_dir, "exists."))
      } else {
      stop(paste("The database directory", db_dir, "does not exist."))
    }
  
    if (file.exists(fasta_dir)) {
      print(paste("The fasta file directory", fasta_dir, "exists."))
      } else {
      stop(paste("The fasta file directory", fasta_dir, "does not exist."))
    }
    
    # #my_example_function <- function(ARG_VALUE_1, ARG_VALUE_2) { 
    # igblast_wrap <- BasiliskEnvironment(envname="igblast_wrap_basilisk",
    #                                     pkgname = "scifer",
    #                                     packages=c("python=3.9", "dnaio==1.2.1", "igblast==1.22.0")
    # )

    # env <- BasiliskEnvironment(envname="igblast_wrap_basilisk",
    #                            pkgname = "scifer",
    #                            packages=c("python=3.9", "dnaio==1.2.1", "igblast==1.22.0")
    # )
    
    
      proc <- basiliskStart(env)
      on.exit(basiliskStop(proc))
      
      some_useful_thing <- basiliskRun(proc, fun=function(arg1, arg2, arg3) {
        # csv <- reticulate::import("csv") 
        # errno <- reticulate::import("errno")
        # logging <- reticulate::import( "logging")
        # multiprocessing <- reticulate::import("multiprocessing")
        # os <- reticulate::import("os")
        # re <- reticulate::import("re")
        # shlex <- reticulate::import("shlex")
        # subprocess <- reticulate::import("subprocess")
        # sys <- reticulate::import("sys")
        # tempfile <- reticulate::import("tempfile")
        # time <- reticulate::import("time")
        # dnaio <- reticulate::import("dnaio")
        
        try(source_python("inst/script/igblastwrap.py"), silent = TRUE)
        #source_python("inst/script/igblastwrap.py")
        #df <- system(run_igblastwrap(database, fasta) 
        df <- system(paste("python inst/script/igblastwrap.py --threads", threads,  database, fasta , sep = " "),
                     intern = TRUE)  
        
        results_airr = read.table(text = df, header = TRUE, sep = "\t")
        
        return(final_output = results_airr)
        
         
         
      }, arg1=database, arg2=fasta, arg3=threads)
      
      some_useful_thing
    }
    
    
    
    
    
#     
#     
#     if (py_available()==T) {print("Conda environment already loaded")} else {
#     conda_envs <- conda_list()
#     print("Activating conda environment")
#     
#       if ("igblast-wrap3" %in% conda_envs$name) {
#             use_condaenv("igblast-wrap3")
#           print("Conda environment active")
# 
#           } else { print("Creating conda environment in python version 3.9") 
#                   conda_create(envname = "igblast-wrap3", packages = c("dnaio", "igblast"),
#                                 channel = "bioconda", python_version = 3.9)
#                     use_condaenv("igblast-wrap3")
#                     print("Conda environment created and active")
#           }
#     } #basilisk end
# 
#     csv <- import("csv") 
#     errno <- import("errno")
#     logging <- import( "logging")
#     multiprocessing <- import("multiprocessing")
#      os <- import("os")
#      re <- import("re")
#      shlex <- import("shlex")
#      subprocess <- import("subprocess")
#      sys <- import("sys")
#      tempfile <- import("tempfile")
#     time <- import("time")
#     dnaio <- import("dnaio")
#     # from argparse import ArgumentParser
#     # from contextlib import ExitStack
#     # from dataclasses import dataclass
#     # from io import StringIO
#     # from itertools import islice
#     # from pathlib import Path
#     # from typing import Dict, Optional
#     
#   
#   #python_script<- 
#   try(source_python("inst/script/igblastwrap.py"), silent = TRUE)
#   #source_python("inst/script/igblastwrap.py")
#   #df <- system(run_igblastwrap(database, fasta) 
#   df <- system(paste("python inst/script/igblastwrap.py --threads", threads,  database, fasta , sep = " "),
#                         intern = TRUE)  
#                
#   results_airr = read.table(text = df, header = TRUE, sep = "\t")
#   
#   return(final_output = results_airr)
#   
# }

# db1 <- 'inst/extdata/test_fasta/KIMDB_rm'
# fasta1 <- 'inst/extdata/test_fasta/test_igblast.txt'
# 
# igblast(db1, fasta1,)









