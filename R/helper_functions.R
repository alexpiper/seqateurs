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
ps_to_fasta <- function(ps = ps, out.file = NULL, seqnames = "unique", width = 1000, ...){
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
#' Thanks to Michael mclaren and Erik Wright for code
#' @param seqtab A DADA2 generated ASV table with columns as sequences and rows as samples
#' @param method An agglomeration method to parse to DECIPHER::IdClusters.
#' This should be (an abbreviation of) one of "complete", "single", "UPGMA", "WPGMA", "NJ", "ML", or "inexact".
#' (See help page on DECIPHER::IdClusters for more information.)
#' @param similarity A similarity to cluster at
#' @param rename Whether clusters should be renamed to OTU<cluster #>
#' @param cores The number of processsor cores to use
#'
#' @return
#' @export
#' @import Biostrings
#' @import DECIPHER
#' @import tibble
#'
#' @examples
cluster_otus <- function(seqtab, method="complete", similarity=0.97, rename=FALSE, cores=1) {

  asv_sequences <- colnames(seqtab)
  sample_names <- rownames(seqtab)
  dna <- Biostrings::DNAStringSet(asv_sequences)

  cutoff <- 100 - similarity

  ## Find clusters of ASVs to form the new OTUs
  aln <- DECIPHER::AlignSeqs(dna, processors = cores)
  d <- DECIPHER::DistanceMatrix(aln, processors = cores)
  clusters <- DECIPHER::IdClusters(
    d,
    method = method,
    cutoff = cutoff, # use `cutoff = 0.03` for a 97% OTU
    processors = cores)

  # add sequences to the `clusters` data frame
  clusters <- clusters %>%
    tibble::add_column(sequence = asv_sequences)

  # Add rowsums
  merged_seqtab <- seqtab %>%
    t %>%
    rowsum(clusters$cluster) %>%
    t
  #Rename
  if(rename){
    colnames(merged_seqtab) <- paste0("OTU", colnames(merged_seqtab))
  }
  return(merged_seqtab)
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


