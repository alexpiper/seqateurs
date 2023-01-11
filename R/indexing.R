# Summarise_index ---------------------------------------------------------

#' Summarise index
#'
#' @param fq A path to a fastq file
#' @param qualityType The quality encoding of the fastq file
#'
#' @return
#' @export
#' @import ShortRead
#' @import tidyr
#' @import dplyr
#' @import stringr
#'
#' @examples
summarise_index <- function(fq, qualityType = "Auto") {
  fF <- ShortRead::FastqStreamer(fq)
  on.exit(close(fF))

  #Stream through fasta
  while( length(suppressWarnings(fqF <- yield(fF, qualityType = qualityType)))) {
    idF <- ShortRead::id(fqF)
  }

  # Check if index is present in fastq header
  index_check <- stringr::str_extract(as.character(idF[[1]]), pattern="(?!:)(?:.(?!:))+$")
  if(stringr::str_detect(index_check, "[A-Z]+[A-Z]")){
    index <- stringr::str_extract(as.character(idF), pattern="(?!:)(?:.(?!:))+$") %>%
      table() %>%
      as.data.frame() %>%
      tidyr::separate(col=1, into=c("index", "index2")) %>%
      dplyr::arrange(desc(Freq))
  } else{
    warning(paste0("Index sequences were not present in fastq header for sample: ", fq))
    index <- data.frame(index = NA_character_, index2=NA_character_, Freq = NA_integer_)
  }
  return(index)
}


#' Create_mismatch
#'
#' @param dna
#' @param dist
#' @param ...
#'
#' @return
#' @importFrom data.table tstrsplit
#' @importFrom data.table CJ
#' @import stringdist
#' @export
#'
#' @examples
create_mismatch <- function(dna, dist=1, ...) {
  all_bases <- c("A", "T", "C", "G")
  l <- data.table::tstrsplit(dna, "", fixed = TRUE)
  l <- lapply(l, function(x) all_bases)
  r <- Reduce(paste0, do.call(data.table::CJ, l))
  return(r[which(stringdist::stringdist(dna, r, method = "hamming") <= dist)])
}

