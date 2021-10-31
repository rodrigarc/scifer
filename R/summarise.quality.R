#' B cell receptor immune repertoire quality control and analysis
#'
#' RepertoiR was created to do quality control of VDJ sequencing, combining single-cell sorted data and specificity with sanguer sequences.
#'
#' @param x Numeric vector.
#'
#' @return Factor variable.
#'
#' @examples
#'sf <- summarise.quality("abi/folder/path", secondary.peak.ratio = 0.33,
#'                     trim.cutoff = 0.01, processors = 1)
#'
#' @export summarise.quality
get.processors <- function(processors){

  if(Sys.info()["sysname"] == 'Windows'){
    # mclapply is not supported on windows
    # so we give a single processor,
    # in which case mclapply calls fall back
    # on lapply
    return(1)
  }

  if(is.null(processors)){
    processors = detectCores(all.tests = FALSE, logical = FALSE)
  }

  return(processors)

}

summarise.quality <- function(input.folder, trim.cutoff = 0.01, secondary.peak.ratio = 0.33, write.secondary.peak.files = FALSE, processors = NULL){

  processors = get.processors(processors)

  print("Looking for .ab1 files...")
  abi.fnames = list.files(input.folder, pattern = "\\.ab1$", full.names = T, recursive = T)

  print(sprintf(("Found %d .ab1 files..."), length(abi.fnames)))


  print("Loading reads...")
  abi.seqs = mclapply(abi.fnames, read.abif, mc.cores = processors)

  print("Calculating read summaries...")
  # now make a data.frame of summaries of all the files
  summaries.dat = mclapply(abi.seqs,
                           summarise.abi.file,
                           trim.cutoff = trim.cutoff,
                           secondary.peak.ratio = secondary.peak.ratio,
                           processors = 1,
                           mc.cores = processors
  )

  print("Cleaning up")
  summaries = mclapply(summaries.dat, function(x) x[["summary"]], mc.cores = processors)
  summaries = do.call(rbind, summaries)

  folder.names = basename(dirname(abi.fnames))
  file.names = basename(abi.fnames)

  summaries = cbind.data.frame("file.path" = as.character(abi.fnames), "folder.name" = as.character(folder.names), "file.name" = file.names, summaries, stringsAsFactors = FALSE)
  #edited by rodrigo
  qual_scores = mclapply(summaries.dat, function(x) x[["quality_score"]], mc.cores = processors)
  names(qual_scores) = as.character(abi.fnames)
  #edited by rodrigo
  return(list("summaries" = summaries, "quality_scores" = qual_scores))

}

loadread <- function(fname, trim, trim.cutoff, revcomp, max.secondary.peaks, secondary.peak.ratio, min.length, processors){

  read.abi = read.abif(fname)

  s = summarise.abi.file(read.abi, trim.cutoff, secondary.peak.ratio, processors = processors)

  summary = s$summary

  # here we store the secondary peaks by storing them in a single sequence
  # as ambiguity codes. Note, this is version of a read that we use.
  # So in this package, a read has an ambiguit whereever there's a
  # secondary peak
  d = c(DNAStringSet(s$read@primarySeq), DNAStringSet(s$read@secondarySeq))
  read = ConsensusSequence(d)[[1]]

  if(trim == TRUE){
    trim.start = summary["trim.start"]
    trim.finish = summary["trim.finish"]
    sp = summary["trimmed.secondary.peaks"]

  }else if(trim == FALSE){
    trim.start = 1
    trim.finish = length(read)
    sp = summary["raw.secondary.peaks"]
  }

  # FILTER read based on user specified limits
  read = read[trim.start:trim.finish]

  if(!is.null(max.secondary.peaks)){
    if(sp > max.secondary.peaks){
      read = NULL
    }
  }

  if(length(read) < min.length){
    read = NULL
  }

  if(!is.null(read)) {
    if(revcomp == TRUE){
      read = reverseComplement(read)
    }
  }
  return(list('read' = read, summary = summary))

}
