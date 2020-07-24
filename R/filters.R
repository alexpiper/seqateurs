# Get Reading frame ---------------------------------------------------------------

#' Get Reading frame of sequences
#'
#' @param x Sequences in DNAStringset or DNAbin format
#' @param genetic.code A genetic code for the Amino acid translation. See all known codes at GENETIC_CODE_TABLE
#' Default is the invertebrate mitochondrial code 'SGC4'
#' @param forward whether the forward complement should be returned
#' @param reverse Whether the reverse complemement should be returned
#' @param resolve_draws How draws should be resolved when multiple possible frames produce sequences with no stop codons.
#'
#' @return
#' @export
#'
#' @import ape
#' @import Biostrings
#'
#' @examples
get_reading_frame <- function(x, genetic.code = "SGC4", forward=TRUE, reverse=FALSE, resolve_draws="majority") {
  # Convert to DNAStringSet
  if (is(x, "DNAbin")) {
    x <- x %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
  }
  #Check for N's
  if(hasOnlyBaseLetters(x) == FALSE) {stop("Error: Sequences contain ambiguities")}
  #discards <- sapply(x, hasOnlyBaseLetters)
  #x <- insect::subset.DNAbin(x, subset = !discards)

  if(forward==TRUE) {
    F_frames <- lapply(1:3, function(pos) subseq(x, start=pos))
  }
  if(reverse==TRUE) {
    R_frames <- lapply(1:3, function(pos) subseq(reverseComplement(x), start=pos))
  }
  #Translate all reading frames
  suppressWarnings(translated <- lapply(F_frames, translate, genetic.code = getGeneticCode(genetic.code)))
  #select the reading frames that contain 0 stop codons, or return NA
  reading_frame <- vector("integer", length=length(x))
  for (i in 1:length(x)){
    fvec = c(str_count(as.character(translated[[1]][i]), "\\*"),
             str_count(as.character(translated[[2]][i]), "\\*"),
             str_count(as.character(translated[[3]][i]), "\\*"))
    if(sum(fvec==0)==1){
      reading_frame[i] <- which(fvec==0)
    } else if(sum(fvec==0)>1) {
      reading_frame[i] <- 0
    }else if(sum(fvec==0)==0) {
      reading_frame[i] <- NA
    }
  }
  if (resolve_draws == "majority") {
    reading_frame[reading_frame==0] <- reading_frame[which.max(tabulate(reading_frame))]
  } else if (resolve_draws == "remove") {
    reading_frame[reading_frame==0] <- NA
  }
  return(reading_frame)
}

# Codon_filter ------------------------------------------------------------

#' Filter sequences containing stop codons
#'
#' @param x Sequences in DNAStringset or DNAbin format
#' @param genetic.code A genetic code for the Amino acid translation. See all known codes at GENETIC_CODE_TABLE
#' Default is the invertebrate mitochondrial code 'SGC4'
#' @param forward Whether the forward complement should be used
#' @param reverse Whether the reverse complement should be used
#' @param remove.ambiguities Whether ambiguous bases (non ATGC) bases should be removed
#'
#' @return
#' @export
#' @import ape
#' @import Biostrings
#'
#' @examples
codon_filter <- function(x, genetic.code = "SGC4", forward=TRUE, reverse=FALSE, remove.ambiguities=TRUE){
  # Convert to DNAStringSet
  if (is(x, "DNAbin")) {
    x <- x %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
  }
  #Check for N's
  if(hasOnlyBaseLetters(x) == FALSE & remove.ambiguities == FALSE) {
    stop("Error: Sequences containing ambiguities, set remove.ambiguities to TRUE")
  } else if(hasOnlyBaseLetters(x) == FALSE & remove.ambiguities==TRUE) {
    discards <- sapply(x, hasOnlyBaseLetters)
    discards <- discards[discards == FALSE]
    x <- x[which(!names(x) %in% names(discards))]
    message(paste0(length(discards), " Sequences containing ambiguities were removed"))
  }

  #Get reading frames
  frames <- get_reading_frame(x, genetic.code = genetic.code, forward = forward, reverse = reverse)

  out <- x[!is.na(frames)]
  message(paste0(length(x) - length(out), " Sequences containing stop codons removed"))
  return(out)
}
