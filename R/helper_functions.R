
# Create mismatch ---------------------------------------------------------


#' Create_mismatch
#'
#' @param dna
#' @param dist
#' @param ...
#'
#' @return
#' @import data.table
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



# Convert to proportions --------------------------------------------------

#' Convert phyloseq table to proportions
#'
#' @param x
#' @param thresh
#' @param na_rm
#' @param ...
#'
#' @return
#'
#' @examples
proportions <- function(x, thresh = NA, na_rm = FALSE, ...) {
  xprop <- (x / sum(x)) # Convert to proportions
  xprop[xprop <= thresh] <- NA ## remove taxa under this level
  xprop2 <- (xprop / sum(xprop, na.rm = na_rm))
  return(xprop2)
}


# ps_to_fasta --------------------------------------------------------------


#' outputs a FASTA file from a phyloseq object
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
#' @param rank (optional) A taxonomic rank from the \code{\link[phyloseq]{tax_table}} which will be used to name the sequences.
#' If no rank is supplied, samples will be named \code{ASV_#}
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

ps_to_fasta <- function(ps = ps, out.file = NULL, rank = NULL, width = 1000, ...){

  if (is.null(ps)){
    message("Phyloseq object not found.")
  }

  if (is.null(out.file)){
    out.file <- paste0(deparse(substitute(ps)), ".fasta")
  }

  if (!is.null(refseq(ps, errorIfNULL = FALSE))){
    seqs <- DNAStringSet(as.vector(refseq(ps)))
  } else{
    message("refseq() not found. Using taxa names for sequences.")
    if (sum(grepl("[^ACTG]", rownames(tax_table(ps)))) > 0){
      stop("Error: Taxa do not appear to be DNA sequences.")
    }
    seqs <- Biostrings::DNAStringSet(colnames(get_taxa(ps)))
  }

  if (is.null(rank) || !rank %in% rank_names(ps)){
    message("Rank not found. Naming sequences sequentially (i.e. ASV_#).")
    names(seqs) <- paste0("ASV_", 1:ntaxa(ps))
  } else {
    names(seqs) <- make.unique(unname(tax_table(ps)[,rank]), sep = "_")
  }

  Biostrings::writeXStringSet(seqs, filepath = out.file, width=width, ... = ...)
  message(paste0(ntaxa(ps), " sequences written to <", out.file, ">."))
}



# Fast melt ---------------------------------------------------------------


#' Fast melt
#'
#' @param physeq
#'
#' @return
#' @import data.table
#' @export
#'
#' @examples
fast_melt  <- function(physeq){
  # supports "naked" otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "taxaID")
  # Enforce character taxaID key
  otudt[, taxaIDchar := as.character(taxaID)]
  otudt[, taxaID := NULL]
  setnames(otudt, "taxaIDchar", "taxaID")
  # Melt count table
  mdt = melt.data.table(otudt,
                        id.vars = "taxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  # Remove zeroes, NAs
  mdt <- mdt[count > 0][!is.na(count)]
  # Calculate relative abundance
  mdt[, RelativeAbundance := count, by = SampleID]
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
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
#' @param Rank
#' @param GroupBy
#'
#' @return
#' @import data.table
#' @export
#'
#' @examples
summarise_taxa <-  function(physeq, Rank, GroupBy = NULL){
  Rank <- Rank[1]
  if(!Rank %in% rank_names(physeq)){
    message("The argument to `Rank` was:\n", Rank,
            "\nBut it was not found among taxonomic ranks:\n",
            paste0(rank_names(physeq), collapse = ", "), "\n",
            "Please check the list shown above and try again.")
  }
  if(!is.null(GroupBy)){
    GroupBy <- GroupBy[1]
    if(!GroupBy %in% sample_variables(physeq)){
      message("The argument to `GroupBy` was:\n", GroupBy,
              "\nBut it was not found among sample variables:\n",
              paste0(sample_variables(physeq), collapse = ", "), "\n",
              "Please check the list shown above and try again.")
    }
  }
  # Start with fast melt
  mdt = fast_melt(physeq)
  if(!is.null(GroupBy)){
    # Add the variable indicated in `GroupBy`, if provided.
    sdt = data.table(SampleID = sample_names(physeq),
                     var1 = get_variable(physeq, GroupBy))
    setnames(sdt, "var1", GroupBy)
    # Join
    setkey(sdt, SampleID)
    setkey(mdt, SampleID)
    mdt <- sdt[mdt]
  }
  # Summarize
  summarydt = mdt[, list(totalRA = sum(RelativeAbundance)),
                  by = c(Rank, GroupBy)]
  return(summarydt)
}


# Compositional functions -------------------------------------------------

#' Geometric mean
#'
#' @param x
#' @param na.rm
#'
#' @return
#' @export
#'
#' @examples
gm_mean <- function(x, na.rm=TRUE){
  # The geometric mean, with some error-protection bits.
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm=na.rm) / length(x))
}



#' Centred log-ratio transformation
#'
#' @param x
#' @param base
#'
#' @return
#' @export
#'
#' @examples
clr <- function(x, base=2){
  x <- log((x / gm_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
}


# Create Samplesheet ------------------------------------------------------

#' Create Samplesheet
#'
#' @param SampleSheet
#' @param runParameters
#' @param format
#' @param Fprimer
#' @param Rprimer
#' @param Ftwintag
#' @param Rtwintag
#'
#' @return
#' @export
#' @import XML
#' @import dplyr
#' @import readr
#' @import magrittr
#' @import tidyr
#' @import tibble
#'
#' @examples
create_samplesheet <- function(SampleSheet, runParameters, format = "miseq"){

  if (format=="miseq"){
    sample_sheet <- readr::read_csv(SampleSheet, skip=20, col_types = cols(
      Sample_ID = col_character(),
      Sample_Name = col_character(),
      Sample_Plate = col_double(),
      Sample_Well = col_character(),
      I7_Index_ID = col_character(),
      index = col_character(),
      I5_Index_ID = col_character(),
      index2 = col_character(),
      Sample_Project = col_character(),
      Description = col_logical()
    ))

    sample_header <- readr::read_csv(SampleSheet, n_max=19, col_types = cols(
      `[Header]` = col_character()
    )) %>%
      dplyr::select(1:2) %>%
      magrittr::set_colnames(c("var", "value")) %>%
      tidyr::drop_na(var) %>%
      dplyr::mutate(var = str_replace(make.unique(var), ".1", "_R")) %>%
      dplyr::mutate(var = str_replace(var, " ", "_")) %>%
      tibble::column_to_rownames("var") %>%
      t() %>%
      tibble::as_tibble() %>%
      dplyr::select(Investigator_Name,
             Project_Name,
             Experiment_Name,
             Assay,
             Adapter)

    reads <- readr::read_csv(SampleSheet,skip=12, n_max=2, col_types = cols(
      `[Reads]` = col_number() ))%>%
      pull(`[Reads]`)
    reads <- tibble(Fread = reads[1], Rread = reads[2])

    xmlFromRunParameters <- XML::xmlParse(runParameters)
    run_params <- XML::xmlToDataFrame(nodes = XML::getNodeSet(xmlFromRunParameters, "/RunParameters")) %>%
      as.data.frame() %>%
      dplyr::mutate(FlowCellExpiry = FlowcellRFIDTag %>%
               str_replace("^.{0,23}", "") %>%
               str_replace(".{0,9}$", "") %>%
               as.Date(),
             ReagentKitExpiry = ReagentKitRFIDTag %>%
               str_replace("^.{0,23}", "") %>%
               str_replace(".{0,9}$", "") %>%
               as.Date(),
             PR2Expiry = PR2BottleRFIDTag %>%
               str_replace("^.{0,23}", "") %>%
               str_replace(".{0,9}$", "") %>%
               as.Date(),
             FCID = Barcode %>%
               str_replace("^.{0,10}", ""),
             RunStartDate = lubridate::ymd(RunStartDate)
      ) %>%
      dplyr::select(
        RunID,
        ScannerID,
        RunNumber,
        FCID,
        RunStartDate,
        PR2BottleBarcode,
        ReagentKitBarcode,
        FlowCellExpiry,
        ReagentKitExpiry,
        PR2Expiry,
        MostRecentWashType) %>%
      dplyr::mutate_if(is.factor, as.character)
  } else {
    message("Warning: Only miseq currently implemented")
    return(NULL)
  }

  combined <- sample_sheet %>%
    cbind(sample_header) %>%
    cbind(reads) %>%
    cbind(run_params)

    return(combined)
}


