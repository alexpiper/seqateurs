# ps_to_fasta --------------------------------------------------------------
#' output a FASTA file from a phyloseq object
#'
#' @description This function outputs a FASTA-formatted text file from a \code{phyloseq} object.
#' This code was modified from \code{reltools} package https://github.com/DanielSprockett/reltools by Daniel Sprocket
#'
#' @param ps A \code{phyloseq} object that contains \code{\link[phyloseq]{refseq}}.
#' If there the \code{refseq} slot is not filled, this function will try pull the
#' sequences from \code{\link[phyloseq]{get_taxa}}
#'
#' @param file (optional) A file name that ends in ".fasta" or ".fa".
#' If a file name is not supplied, the file will be named after the phyloseq object.
#'
#' @param seqnames (optional) A taxonomic rank from the \code{\link[phyloseq]{tax_table}} which will be used to name the sequences.
#' alternatively "unique" will uniquely name the sequences as \code{ASV_#},
#' "sequence" will name the sequences in the fasta exactly the same same as the sequence, usfull for LULU curation.
#'
#' @param width (Default 1000) The number of characters in each fasta line before wrapping occurs
#'
#' @param ... (Optional) Any further paramaters to be passed to \code{writeXStringSet}
#'
#' @return This function saves a FASTA-formatted text file from the input \code{phyloseq} object.
#' @export
#' @import phyloseq
#' @import Biostrings
#'
#' @examples
#' save_fasta(ps)
#' save_fasta(ps = ps, file = "sequences.fasta", rank = "Genus")
ps_to_fasta <- function(ps, out.file = NULL, seqnames = "unique", width = 1000, ...){
  if (is.null(ps)){
    message("Phyloseq object not found.")
  }

  if (is.null(out.file)){
    out.file <- paste0(deparse(substitute(ps)), ".fa")
  }

  if (!is.null(phyloseq::refseq(ps, errorIfNULL = FALSE))){
    seqs <- Biostrings::DNAStringSet(as.vector(phyloseq::refseq(ps)))
  } else{
    message("refseq() not found. Using tax_table rownames for sequences.")
    if (sum(grepl("[^ACTG]", rownames(phyloseq::tax_table(ps)))) > 0){
      stop("Error: Taxa do not appear to be DNA sequences.")
    }
    seqs <- Biostrings::DNAStringSet(colnames(phyloseq::get_taxa(ps)))
  }

  if(seqnames %in% phyloseq::rank_names(ps)){
    names(seqs) <- make.unique(unname(phyloseq::tax_table(ps)[,seqnames]), sep = "_")
  } else if(seqnames == "unique" || is.null(seqnames) ){
    message("Rank not found. Naming sequences sequentially (i.e. ASV_#).")
    names(seqs) <- paste0("ASV_", 1:phyloseq::ntaxa(ps))
  } else if(seqnames == "sequence"){
    names(seqs) <- as.character(seqs)
  }

  Biostrings::writeXStringSet(seqs, filepath = out.file, width=width, ... = ...)
  message(paste0(phyloseq::ntaxa(ps), " sequences written to <", out.file, ">."))
}


#' ps_to_dnabin
#'
#'
#' @param ps A \code{phyloseq} object that contains \code{\link[phyloseq]{refseq}}.
#' If there the \code{refseq} slot is not filled, this function will try pull the
#' sequences from \code{\link[phyloseq]{get_taxa}}
#'
#'
#' @param seqnames (optional) A taxonomic rank from the \code{\link[phyloseq]{tax_table}} which will be used to name the sequences.
#' alternatively "unique" will uniquely name the sequences as \code{ASV_#},
#' "sequence" will name the sequences in the fasta exactly the same same as the sequence, usfull for LULU curation.
#'
#' @export
#' @import phyloseq
#' @import Biostrings
#' @importFrom ape as.DNAbin
#' @importFrom stringr str_to_sentence
#'
#' @examples
#' save_fasta(ps)
#' save_fasta(ps = ps, file = "sequences.fasta", rank = "Genus")
ps_to_dnabin <- function (ps, seqnames = "unique") {
  if (is.null(ps)) {
    message("Phyloseq object not found.")
  }
  if (!is.null(phyloseq::refseq(ps, errorIfNULL = FALSE))) {
    seqs <- ape::as.DNAbin(Biostrings::DNAStringSet(as.vector(phyloseq::refseq(ps))))
  } else {
    message("refseq() not found. Using tax_table rownames for sequences.")
    if (sum(grepl("[^ACTG]", rownames(phyloseq::tax_table(ps)))) >
        0) {
      stop("Error: Taxa do not appear to be DNA sequences.")
    }
    seqs <- ape::as.DNAbin(Biostrings::DNAStringSet(colnames(phyloseq::get_taxa(ps))))
  }

  seqnames <- stringr::str_to_sentence(seqnames)
  if (seqnames %in% phyloseq::rank_names(ps)) {
    names(seqs) <- make.unique(unname(phyloseq::tax_table(ps)[, seqnames]), sep = "_")
  } else if (seqnames == "unique" || is.null(seqnames)) {
    message("Rank not found. Naming sequences sequentially (i.e. ASV_#).")
    names(seqs) <- paste0("ASV_", 1:phyloseq::ntaxa(ps))
  } else if (seqnames == "sequence") {
    names(seqs) <- as.character(seqs)
  }
  return(seqs)
}


# Fast melt ---------------------------------------------------------------
#' Fast melt
#'
#' @param physeq
#'
#' @return
#' @import data.table
#' @import phyloseq
#' @export
#'
#' @examples
fast_melt  <- function(physeq){
  # supports "naked" otu_table as `physeq` input.
  otutab = as(phyloseq::otu_table(physeq), "matrix")
  if(!phyloseq::taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "taxaID")
  # Enforce character taxaID key
  otudt[, taxaIDchar := as.character(taxaID)]
  otudt[, taxaID := NULL]
  setnames(otudt, "taxaIDchar", "taxaID")
  # Melt count table
  mdt = data.table::melt.data.table(otudt,
                        id.vars = "taxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  # Remove zeroes, NAs
  mdt <- mdt[count > 0][!is.na(count)]
  # Calculate relative abundance
  mdt[, RelativeAbundance := count, by = SampleID]
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(phyloseq::tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "taxaID")
    # Enforce character taxaID key
    taxdt[, taxaIDchar := as.character(taxaID)]
    taxdt[, taxaID := NULL]
    setnames(taxdt, "taxaIDchar", "taxaID")
    # Join with tax table
    setkey(taxdt, "taxaID")
    setkey(mdt, "taxaID")
    mdt <- taxdt[mdt]
  }
  return(mdt)
}


# summarise_taxa ----------------------------------------------------------

#' Summarise_taxa
#'
#' @param physeq
#' @param rank
#' @param group_by
#'
#' @return
#' @import data.table
#' @import phyloseq
#' @export
#'
#' @examples
summarise_taxa <-  function(physeq, rank, group_by = NULL){
  rank <- rank[1]
  if(!rank %in% phyloseq::rank_names(physeq)){
    message("The argument to `rank` was:\n", rank,
            "\nBut it was not found among taxonomic ranks:\n",
            paste0(phyloseq::rank_names(physeq), collapse = ", "), "\n",
            "Please check the list shown above and try again.")
  }
  if(!is.null(group_by)){
    group_by <- group_by[1]
    if(!group_by %in% phyloseq::sample_variables(physeq)){
      message("The argument to `group_by` was:\n", group_by,
              "\nBut it was not found among sample variables:\n",
              paste0(phyloseq::sample_variables(physeq), collapse = ", "), "\n",
              "Please check the list shown above and try again.")
    }
  }
  # Start with fast melt
  mdt = fast_melt(physeq)
  if(!is.null(group_by)){
    # Add the variable indicated in `group_by`, if provided.
    sdt = data.table(SampleID = phyloseq::sample_names(physeq),
                     var1 = phyloseq::get_variable(physeq, group_by))
    setnames(sdt, "var1", group_by)
    # Join
    setkey(sdt, SampleID)
    setkey(mdt, SampleID)
    mdt <- sdt[mdt]
  }
  # Summarize
  summarydt = mdt[, list(totalRA = sum(RelativeAbundance)),
                  by = c(rank, group_by)]
  return(summarydt)
}

# Cluster OTUS ------------------------------------------------------------

#' OTU clustering of DADA2 seqtab
#'
#' @param x The input object can be a DADA2 generated ASV table with columns as sequences and rows as samples,
#' a phyloseq object with or without refseqs, or a DNAbin or DNAStringSet.
#' @param method An agglomeration method to parse to DECIPHER::IdClusters.
#' This should be (an abbreviation of) one of "complete", "single", "UPGMA", "WPGMA", "NJ", "ML", or "inexact".
#' (See help page on DECIPHER::IdClusters for more information.)
#' @param similarity A similarity threshold to cluster at. Must be a number between 0 and 1
#' @param cores The number of processsor cores to use
#'
#' @return
#' @export
#' @import Biostrings
#' @import DECIPHER
#' @import dplyr
#'
#' @examples
cluster_otus <- function(x, method="complete", similarity=0.97, cores=1) {

  if(is(x, "matrix")| is(x, "data.frame")){
    asv_sequences <- colnames(x)
  } else if(is(x, "phyloseq") & !is.null(phyloseq::refseq(x, errorIfNULL = FALSE))){
    asv_sequences <- as.vector(phyloseq::refseq(x))
  } else if(is(x, "phyloseq") & is.null(phyloseq::refseq(x, errorIfNULL = FALSE))){
    message("refseq() not found. Using tax_table rownames for sequences.")
    if (sum(grepl("[^ACTG]", rownames(phyloseq::tax_table(x)))) > 0) {
      stop("Error: Taxa do not appear to be DNA sequences.")
    }
    asv_sequences <- colnames(phyloseq::get_taxa(x))
  } else if(is(x, "DNAStringSet") ){
    asv_sequences <- as.character(x)
  } else if(is(x, "DNAbin")){
    asv_sequences <- taxreturn::DNAbin2char(x)
  } else{
    stop("Error: Taxa do not appear to be DNA sequences.")
  }
  seqs <- Biostrings::DNAStringSet(asv_sequences)

  # define cutoffs for clustering
  if(!dplyr::between(similarity, 0, 1)){
    stop("similarity must be a number between 0 and 1")
  }
  cutoff <- 1 - similarity

  ## Find clusters of ASVs to form the new OTUs
  aln <- DECIPHER::AlignSeqs(seqs, processors = cores)
  d <- DECIPHER::DistanceMatrix(aln, processors = cores)
  otus <- DECIPHER::IdClusters(
    d,
    method = method,
    cutoff = cutoff, # use `cutoff = 0.03` for a 97% OTU
    processors = cores) %>%
    dplyr::mutate(sequence = asv_sequences)  %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(cluster_size = n_distinct(sequence)) %>%
    dplyr::ungroup()

  # Return cluster memberships
  return(otus)
}
# Propagate taxonomic assignments to species level ------------------------

#' Propagate taxonomy
#'
#' @param tax
#' @param from
#'
#' @return
#' @export
#'
#' @examples
propagate_tax <- function(tax, from = "Family") {
  .Deprecated(new="seqateurs::na_to_unclassified", package="seqateurs", old="seqateurs::propagate_tax")
  col.prefix <- substr(colnames(tax), 1, 1) # Assumes named Kingdom, ...

  # Highest level to propagate from
  if (from == "Phylum") (start <- 2)
  if (from == "Class") (start <- 3)
  if (from == "Order") (start <- 4)
  if (from == "Family") (start <- 5)
  if (from == "Genus") (start <- 6)
  if (from == "Species") (start <- 7)

  # Propagate
  for (col in seq(start, ncol(tax))) {
    prop <- is.na(tax[, col]) & !is.na(tax[, col - 1])
    newtax <- tax[prop, col - 1]
    needs.prefix <- !grepl("^[A-z]__", newtax)
    newtax[needs.prefix] <- paste(col.prefix[col - 1], newtax[needs.prefix], sep = "__")
    tax[prop, col] <- newtax
  }
  tax
}


# NA_to_unclassified ------------------------------------------------------


#' Convert NA classifications to the their highest successful classification
#'
#' @param x A matrix, data frame or phyloseq object
#' @param sep The seperator between the rank, and the name
#' i.e. G__Drosophila would be the default for an OTU that could only be assigned to the genus Drosophila
#' @param rownames Whether rownames should be returned, or instead an OTU column
#' @param quiet Whether progress should be printed to console
#'
#' @seealso unclassified_to_na
#' @return
#' @export
#' @import tibble
#' @import dplyr
#' @import phyloseq
#' @import stringr
#'
#' @examples
na_to_unclassified <- function(x, sep="__", rownames=TRUE, quiet=FALSE){
  if(any(class(x) == "phyloseq")){
    tax <- phyloseq::tax_table(x) %>%
      as("matrix")%>%
      as.data.frame()
    ps <- TRUE
    ranks <- phyloseq::rank_names(x)
    if(!quiet){ message("Input is a phyloseq object")}
  } else if (any(class(x) %in% c("matrix", "data.frame"))){
    tax <- as.data.frame(x)
    ranks <- colnames(x)
    ps <- FALSE
  } else(stop("x must be a matrix, data frame or phyloseq object"))

  replacements <- tax %>%
    dplyr::bind_cols(lowest_classified(tax, return="both", rownames=TRUE) %>%
                       dplyr::transmute(lowest = paste0(substr(rank, start = 1, stop = 1), sep, name))) %>%
    tibble::rownames_to_column("OTU") %>%
    dplyr::mutate(across(ranks, function(x){coalesce(x, lowest)})) %>%
    dplyr::select(-lowest)

  #return modified table
  if(ps){
    tax_table(x) <- replacements %>%
      tibble::column_to_rownames("OTU")%>%
      as("matrix") %>%
      phyloseq::tax_table()
    out <- x
  } else if (!ps && rownames==TRUE){
    out <- replacements %>%
      tibble::column_to_rownames("OTU")
  } else {
    out <- replacements
  }
  return(out)
}

# Unclassified_to_na ------------------------------------------------------


#' Convert classifications that are annotated with their highest successful classification back ot NA
#'
#' @param x A matrix, data frame or phyloseq object
#' @param sep The seperator between the rank, and the name
#' i.e. G__Drosophila would be the default for an OTU that could only be assigned to the genus Drosophila
#' @param rownames Whether rownames should be returned, or instead an OTU column
#' @param quiet Whether progress should be printed to console
#'
#' @seealso na_to_unclassified
#' @return
#' @export
#'
#' @import tibble
#' @import dplyr
#' @import phyloseq
#' @import stringr
#' @import tidyr
#'
#' @examples
unclassified_to_na <- function(x, sep="__", rownames=TRUE, quiet=FALSE){
  if(any(class(x) == "phyloseq")){
    tax <- phyloseq::tax_table(x)%>%
      as("matrix") %>%
      as.data.frame()
    ps <- TRUE
    ranks <- phyloseq::rank_names(x)
    if(!quiet){ message("Input is a phyloseq object")}
  } else if (any(class(x) %in% c("matrix", "data.frame"))){
    tax <- as.data.frame(x)
    if(!"OTU" %in% colnames(tax)){
      tax <- tax %>%
        tibble::rownames_to_column("OTU")
    }

    ranks <- colnames(x)[!colnames(x)=="OTU"]
    ps <- FALSE
  } else(stop("x must be a matrix, data frame or phyloseq object"))
  replacements <- tax %>%
    tidyr::pivot_longer(cols=all_of(ranks),
                        names_to = "rank",
                        values_to = "name") %>%
    dplyr::filter(!stringr::str_detect(name, "__")) %>%
    tidyr::pivot_wider(names_from="rank",
                       values_from= "name")

  #return modified table
  if(ps){
    tax_table(x) <- replacements %>%
      tibble::column_to_rownames("OTU")%>%
      as("matrix") %>%
      phyloseq::tax_table()
    out <- x
  } else if (!ps && rownames==TRUE){
    out <- replacements %>%
      tibble::column_to_rownames("OTU")
  } else {
    out <- replacements
  }
  return(out)
}


# Lowest_classified_rank --------------------------------------------------


#' Get the lowest classified ranks for the OTUs
#'
#' @param x A matrix, data frame or phyloseq object
#' @param sep The seperator between the rank, and the name
#' i.e. G__Drosophila would be the default for an OTU that could only be assigned to the genus Drosophila
#' @param return The type of values to return.
#' Options include 'rank' which returns a character vector listing the lowest taxonomic each otu was classified to, i.e. Species, Genus etc
#' 'taxon' returns a character vector listing the lowest taxonomic name each otu was classified to i.e. Drosophila_suzukii, Diptera etc
#' 'both' returns a data frame containing both of the above
#' @param rownames Whether rownames should be returned, or instead an OTU column
#'
#' @return
#' @export
#' @seealso na_to_unclassified unclassified_to_na
#' @import phyloseq
#' @import stringr
#'
#' @examples
lowest_classified <- function(x, sep="__", return="rank", rownames=TRUE){
  if(any(class(x) == "phyloseq")){
    tax <- phyloseq::tax_table(x)%>%
      as("matrix") %>%
      as.data.frame()
    ps <- TRUE
    ranks <- phyloseq::rank_names(x)
  } else if (any(class(x) %in% c("matrix", "data.frame"))){
    tax <- as.data.frame(x)
    ranks <- colnames(x)
  } else(stop("x must be a matrix, data frame or phyloseq object"))
  #Check that NA's are present
  if(any(stringr::str_detect(tax, sep))){
    tax <- unclassified_to_na(tax, sep=sep, rownames=rownames)
  }
  # get lowest classified
  revtax <- rev(tax)
  logidf <- !is.na(revtax)
  keepvec <- unname(apply(logidf, 1, which.max))
  if(return=="rank"){
    lowest <- colnames(logidf)[keepvec]
  } else if(return=="taxon"){
    lowest <- mapply(function(x,y) revtax[x,y], x=seq_len(nrow(revtax)), y=keepvec)
  } else if(return=="both"){
    lowest <- data.frame(
      OTU = ifelse(has_rownames(tax), rownames(tax), tax$OTU),
      rank = colnames(logidf)[keepvec],
      name = mapply(function(x,y) revtax[x,y], x=seq_len(nrow(revtax)), y=keepvec)
    )
  }
  return(lowest)
}


# Coalescing join ---------------------------------------------------------

#' Coalescing join
#'
#' @description A join to combine datasets containing identical non-key columns in varying states of completeness
#' @param x The primary data frame or tibble from which values will be prefered
#' @param y The secondary data frame or tibble from which values will only be used when values in x are NA
#' @param by A character vector of columns to join by.
#' @param suffix The suffixes from the dplyr join to remove.
#' @param join The type of dplyr join to perform
#' @param ...
#'
#' @import purrr
#' @import dplyr
#'
#' @return
#' @export
#'
#' @examples
coalesce_join <- function(x, y,
                          by = NULL, suffix = c(".x", ".y"),
                          join = dplyr::full_join, ...) {
  joined <- join(x, y, by = by, suffix = suffix, ...)
  # names of desired output
  cols <- union(names(x), names(y))

  to_coalesce <- names(joined)[!names(joined) %in% cols]
  suffix_used <- suffix[ifelse(endsWith(to_coalesce, suffix[1]), 1, 2)]
  # remove suffixes and deduplicate
  to_coalesce <- unique(substr(
    to_coalesce,
    1,
    nchar(to_coalesce) - nchar(suffix_used)
  ))

  coalesced <- purrr::map_dfc(to_coalesce, ~dplyr::coalesce(
    joined[[paste0(.x, suffix[1])]],
    joined[[paste0(.x, suffix[2])]]
  ))
  names(coalesced) <- to_coalesce

  dplyr::bind_cols(joined, coalesced)[cols]
}


