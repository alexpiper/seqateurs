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

#' Filter sequences containing stop codons
#'
#' @param x Sequences in DNAStringset or DNAbin format
#' @param genetic_code A genetic code for the Amino acid translation. set to 'SGC4' for Invertebrate mitochondrial or see all known codes at Biostrings::GENETIC_CODE_TABLE
#' @param tryrc Whether the reverse complemement should be evaluated if no frame without stop codons was found in the forward orientation.
#' @param resolve_draws How draws should be resolved when multiple possible frames produce sequences with no stop codons.
#' Options are "remove" to completely remove the sequence, or "majority" to pick the most common frame from the entire alignment.
#'
#' @return
#' @export
#' @importFrom ape as.DNAbin
#' @importFrom Biostrings reverseComplement
#' @importFrom DECIPHER RemoveGaps
#' @importFrom methods is
#'
codon_filter <- function(x, genetic_code = NULL, tryrc=TRUE, resolve_draws="majority"){
  if(is.null(genetic_code)){
    stop("genetic_code must not be NULL, set to 'SGC4' for Invertebrate mitochondrial or see Biostrings::GENETIC_CODE_TABLE for other options")
  }
  # Convert to DNAStringSet
  if (methods::is(x, "DNAbin")) {
    x <- DNAbin2DNAstringset(x, remove_gaps=FALSE)
    format <- "DNAbin"
  } else if(methods::is(x, "DNAStringSet")){
    format <- "DNAStringSet"
  } else {
    stop("x must be a DNAbin or DNAStringSet")
  }

  #Get reading frames
  frames <- get_reading_frame(DECIPHER::RemoveGaps(x), genetic_code = genetic_code, tryrc = tryrc, resolve_draws = resolve_draws)

  # Check if any sequences need RC (negative reading frame)
  to_rc <- sign(frames)==-1
  to_rc[is.na(to_rc)] <- FALSE
  if (any(to_rc)){
    message(sum(to_rc), " Sequences reverse complemented as no forward match was found")
    x[to_rc] <- Biostrings::reverseComplement(x[to_rc])
  }

  if(format=="DNAbin"){
    out <- ape::as.DNAbin(x[!is.na(frames)])
  } else if(format=="DNAStringSet"){
    out <- x[!is.na(frames)]
  }
  message(paste0(length(x) - length(out), " Sequences containing stop codons removed"))
  return(out)
}
