

#### ALL FUNCTIONS TO BE MOVED

# Get Qual Stats ----------------------------------------------------------


#' Get qual stats
#'
#' @param input
#' @param n
#'
#' @return
#' @export
#'
#' @examples
get_qual_stats <- function (input, n = 5e+05) {

  statdf <- data.frame(Cycle = integer(0), Mean = numeric(0),
                       Q25 = numeric(0), Q50 = numeric(0), Q75 = numeric(0),
                       Cum = numeric(0), file = character(0))
  anndf <- data.frame(minScore = numeric(0), label = character(0),
                      rclabel = character(0), rc = numeric(0), file = character(0))
  #Check files arent empty
  fl <- input[file.size(input) > 28]
  if (length(input) > length(fl)) {
    message("Warning, the following files were empty: \n")
    message(paste0(input[!input %in% fl], " \n"))
  }
  # Start loop
  FIRST <- TRUE
  for (f in fl[!is.na(fl)]) {

    # f <- fl[!is.na(fl)][[1]]

    srqa <- qa(f, n = n)
    df <- srqa[["perCycle"]]$quality
    rc <- sum(srqa[["readCounts"]]$read)
    if (rc >= n) {
      rclabel <- paste("Reads >= ", n)
    }  else {      rclabel <- paste("Reads: ", rc)
    }
    means <- rowsum(df$Score * df$Count, df$Cycle)/rowsum(df$Count,  df$Cycle)

    get_quant <- function(xx, yy, q) {
      xx[which(cumsum(yy)/sum(yy) >= q)][[1]]
    }
    q25s <- by(df, df$Cycle, function(x) get_quant(x$Score, x$Count, 0.25), simplify = TRUE)
    q50s <- by(df, df$Cycle, function(x) get_quant(x$Score, x$Count, 0.5), simplify = TRUE)
    q75s <- by(df, df$Cycle, function(x) get_quant(x$Score, x$Count, 0.75), simplify = TRUE)
    cums <- by(df, df$Cycle, function(x) sum(x$Count), simplify = TRUE)
    if (!all(sapply(list(names(q25s), names(q50s), names(q75s), names(cums)), identical, rownames(means)))) {
      stop("Calculated quantiles/means weren't compatible.")
    }
    if (FIRST) {
      plotdf <- cbind(df, file = basename(f))
      FIRST <- FALSE
    }
    else {
      plotdf <- rbind(plotdf, cbind(df, file = basename(f)))
    }

    statdf <- rbind(statdf, data.frame(Cycle = as.integer(rownames(means)),
                                       Mean = means, Q25 = as.vector(q25s), Q50 = as.vector(q50s),
                                       Q75 = as.vector(q75s), Cum = 10 * as.vector(cums)/rc,
                                       file = basename(f)))
    anndf <- rbind(anndf, data.frame(minScore = min(df$Score),
                                     label = basename(f), rclabel = rclabel, rc = rc,
                                     file = basename(f)))
  }

  anndf$minScore <- min(anndf$minScore)

  get_ee <- function(Q){
    ee <- 10^(-Q/ 10)
    return(ee)
  }


  plotdf.summary <- aggregate(Count ~ Cycle + Score, plotdf, sum)
  plotdf.summary$label <- paste(nrow(anndf), "files (aggregated)")
  means <- rowsum(plotdf.summary$Score * plotdf.summary$Count, plotdf.summary$Cycle)/rowsum(plotdf.summary$Count, plotdf.summary$Cycle)
  EEmeans <- get_ee(rowsum(plotdf.summary$Score * plotdf.summary$Count, plotdf.summary$Cycle)/rowsum(plotdf.summary$Count, plotdf.summary$Cycle))
  q10s <- by(plotdf.summary, plotdf.summary$Cycle, function(x) get_quant(x$Score, x$Count, 0.1), simplify = TRUE)
  EE10s <- by(plotdf.summary, plotdf.summary$Cycle, function(x) get_ee(get_quant(x$Score, x$Count, 0.1)), simplify = TRUE)
  q25s <- by(plotdf.summary, plotdf.summary$Cycle, function(x) get_quant(x$Score, x$Count, 0.25), simplify = TRUE)
  EE25s <- by(plotdf.summary, plotdf.summary$Cycle, function(x) get_ee(get_quant(x$Score, x$Count, 0.25)), simplify = TRUE)
  q50s <- by(plotdf.summary, plotdf.summary$Cycle, function(x) get_quant(x$Score, x$Count, 0.5), simplify = TRUE)
  EE50s <- by(plotdf.summary, plotdf.summary$Cycle, function(x) get_ee(get_quant(x$Score, x$Count, 0.50)), simplify = TRUE)
  q75s <- by(plotdf.summary, plotdf.summary$Cycle, function(x) get_quant(x$Score, x$Count, 0.75), simplify = TRUE)
  EE75s <- by(plotdf.summary, plotdf.summary$Cycle, function(x) get_ee(get_quant(x$Score, x$Count, 0.75)), simplify = TRUE)
  q90s <- by(plotdf.summary, plotdf.summary$Cycle, function(x) get_quant(x$Score, x$Count, 0.9), simplify = TRUE)
  EE90s <- by(plotdf.summary, plotdf.summary$Cycle, function(x) get_ee(get_quant(x$Score, x$Count, 0.9)), simplify = TRUE)
  cums <- by(plotdf.summary, plotdf.summary$Cycle, function(x) sum(x$Count), simplify = TRUE)
  statdf.summary <- data.frame(Cycle = as.integer(rownames(means)),
                               QMean = means,
                               EEmean = EEmeans,
                               Q10 = as.vector(q10s),
                               EE10 = as.vector(EE10s),
                               Q25 = as.vector(q25s),
                               EE25 = as.vector(EE25s),
                               Q50 = as.vector(q50s),
                               EE50 = as.vector(EE50s),
                               Q75 = as.vector(q75s),
                               EE75 = as.vector(EE75s),
                               Q90 = as.vector(q90s),
                               EE90 = as.vector(EE90s),
                               reads = 10 * as.vector(cums)/rc)
  attr(statdf.summary, "qsummary") <- "qsummary"
  return(statdf.summary)
}


# Qual plot ---------------------------------------------------------------

#' Quality plot
#'
#' @param input
#' @param n
#'
#' @return
#' @export
#'
#' @examples
qual_plot <- function(input, n = 5e+05) {
    summary <- get_qual_stats(input, n = n )
    plot <- summary %>%
      dplyr::select(Cycle, reads, starts_with("Q")) %>%
      tidyr::pivot_longer(cols = starts_with("Q")) %>%
      ggplot(aes(x=Cycle, y=value, colour=name)) +
        #geom_point(size=1)  +
        geom_line(data = summary, aes(y = QMean), color = "#66C2A5") +
        geom_line(data = summary, aes(y = Q25), color = "#FC8D62", size = 0.25, linetype = "dashed") +
        geom_line(data = summary, aes(y = Q50), color = "#FC8D62", size = 0.25) +
        geom_line(data = summary, aes(y = Q75), color = "#FC8D62",  size = 0.25, linetype = "dashed") +
        #geom_tile(aes(fill = reads))
        labs(x = "Reads position", y = "Quality Score") +
        ylim(c(0, NA))
  return(plot)
}


#' MaxEE plot
#'
#' @param input
#' @param n
#'
#' @return
#' @export
#'
#' @examples
maxEEplot <- function(input, n = 5e+05) {
  summary <- get_qual_stats(input, n = n )
  plot <- summary %>%
    dplyr::select(Cycle, starts_with("EE")) %>%
    tidyr::pivot_longer(cols = starts_with("EE")) %>%
    dplyr::group_by(name) %>%
    dplyr::mutate(cumsumEE = cumsum(value)) %>%
    ggplot(aes(x=Cycle, y=log10(cumsumEE), colour=name)) +
      geom_point(size=1) +
      geom_hline(yintercept = log10(1), color = "red") +
      geom_hline(yintercept = log10(2), color = "red") +
      geom_hline(yintercept = log10(3), color = "red") +
      geom_hline(yintercept = log10(5), color = "red") +
      geom_hline(yintercept = log10(7), color = "red") +
      geom_text(label = "MaxEE=1", aes(x = 0, y = log10(1), hjust = 0, vjust = 0), color = "red") +
      geom_text(label = "MaxEE=2", aes(x = 0, y = log10(2), hjust = 0, vjust = 0), color = "red") +
      geom_text(label = "MaxEE=3", aes(x = 0, y = log10(3), hjust = 0, vjust = 0), color = "red") +
      geom_text(label = "MaxEE=5", aes(x = 0, y = log10(5), hjust = 0, vjust = 0), color = "red") +
      geom_text(label = "MaxEE=7", aes(x = 0, y = log10(7), hjust = 0, vjust = 0), color = "red") +
      labs(x = "Reads position", y = "Log10 Cumulative expected errors")
  return(plot)
}



# Summarise phyloseq ------------------------------------------------------


#' Summarise phyloseq table
#'
#' @param physeq
#' @param Rank
#' @param GroupBy
#'
#' @return
#'
#' @examples
summarize_taxa <- function(physeq, Rank, GroupBy = NULL) {
  Rank <- Rank[1]
  if (!Rank %in% rank_names(physeq)) {
    message(
      "The argument to `Rank` was:\n", Rank,
      "\nBut it was not found among taxonomic ranks:\n",
      paste0(rank_names(physeq), collapse = ", "), "\n",
      "Please check the list shown above and try again."
    )
  }
  if (!is.null(GroupBy)) {
    GroupBy <- GroupBy[1]
    if (!GroupBy %in% sample_variables(physeq)) {
      message(
        "The argument to `GroupBy` was:\n", GroupBy,
        "\nBut it was not found among sample variables:\n",
        paste0(sample_variables(physeq), collapse = ", "), "\n",
        "Please check the list shown above and try again."
      )
    }
  }
  # Start with fast melt
  mdt <- fast_melt(physeq)
  if (!is.null(GroupBy)) {
    # Add the variable indicated in `GroupBy`, if provided.
    sdt <- data.table(
      SampleID = sample_names(physeq),
      var1 = get_variable(physeq, GroupBy)
    )
    setnames(sdt, "var1", GroupBy)
    # Join
    setkey(sdt, SampleID)
    setkey(mdt, SampleID)
    mdt <- sdt[mdt]
  }
  # Summarize
  summarydt <- mdt[, list(totalRA = sum(RelativeAbundance)),
    by = c(Rank, GroupBy)
  ]
  return(summarydt)
}


# Create mismatch ---------------------------------------------------------


#' Create_mismatch
#'
#' @param dna
#' @param dist
#' @param ...
#'
#' @return
#'
#' @examples
create_mismatch <- function(dna, dist, ...) {
  all_bases <- c("A", "T", "C", "G")
  l <- tstrsplit(dna, "", fixed = TRUE)
  l <- lapply(l, function(x) all_bases)
  r <- Reduce(paste0, do.call(CJ, l))
  return(r[which(stringdist(dna, r, method = "hamming") <= dist)])
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
    seqs <- DNAStringSet(colnames(get_taxa(ps)))
  }

  if (is.null(rank) || !rank %in% rank_names(ps)){
    message("Rank not found. Naming sequences sequentially (i.e. ASV_#).")
    names(seqs) <- paste0("ASV_", 1:ntaxa(ps))
  } else {
    names(seqs) <- make.unique(unname(tax_table(ps)[,rank]), sep = "_")
  }

  writeXStringSet(seqs, filepath = out.file, width=width, ... = ...)
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


