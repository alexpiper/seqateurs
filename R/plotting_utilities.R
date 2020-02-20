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
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#'
#' @examples
plot_quality <- function(input, n = 5e+05) {
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
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#'
#' @examples
plot_maxEE <- function(input, n = 5e+05) {
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


#' Plot readlengths of fastq
#'
#' @param dir A directory of fastq files
#' @param aggregate Whether plot should be faceted by sample or aggregated together
#' @param sample How many reads to sample from each fastq, default 1e6
#'
#' @return
#' @export
#'
#' @examples
plot_lengths <- function(dir, aggregate = FALSE, sample = 1e6){
  Forward <-  sort(list.files(dir, pattern="_R1_", full.names = TRUE))
  Reverse <-  sort(list.files(dir, pattern="_R2_", full.names = TRUE))

  lengths <- vector("list", length=length(Forward))
  for(i in 1:length(Forward)){
    #Forwrard reads
    sampler <- FastqSampler(Forward[[i]], n=sample )
    Ffq <- as.data.frame(table(width(yield(sampler))), stringsAsFactors=FALSE) %>%
      magrittr::set_colnames(c("Readlength", "Forward")) %>%
      mutate_if(is.character, as.integer)
    #Reverse reads
    sampler <- FastqSampler(Reverse[[i]], n=sample )
    Rfq <- as.data.frame(table(width(yield(sampler))), stringsAsFactors=FALSE) %>%
      magrittr::set_colnames(c("Readlength", "Reverse"))%>%
      mutate_if(is.character, as.integer)

    lengths[[i]] <- tibble(Readlength=seq(1, max(max(Ffq$Readlength), max(Rfq$Readlength)),1)) %>%
      left_join(Ffq, by="Readlength") %>%
      left_join(Rfq, by="Readlength")
    lengths[[i]]$source <- Forward[i]
  }

  p <- do.call(rbind, lengths) %>%
    group_by(Readlength, source) %>%
    summarise(Forward = sum(Forward), Reverse = sum(Reverse)) %>%
    pivot_longer(cols=3:4,
                 names_to = "Direction",
                 values_to = "Count") %>%
    ggplot(aes(x=Readlength, y=Count, group=Direction, fill=Direction)) +
    geom_bar(stat="identity", colour="black")

  if(aggregate==TRUE){
    p <- p  + facet_grid(~Direction)
  } else {
    p <- p  + facet_grid(source~Direction)
  }
  return(p)
}



# estimate_trunclen ------------------------------------------------------------

#' Automatically estimate trunclen for filtering
#'
#' @param fwd
#' @param rev
#' @param threshold
#' @param maxlength
#' @param minlength
#' @param qa_sample
#' @param overlap_sample
#' @param quiet
#'
#' @return
#' @export
#'
#' @examples
estimate_trunclen <- function (fwd, rev = NULL, threshold = 25, maxlength, minlength, minoverlap=20, qa_sample = 5e+05, overlap_sample = 100, quiet = FALSE) {
  fqa <- get_qual_stats(fwd, n = qa_sample)
  if (!is.null(rev)) {
    rqa <- get_qual_stats(rev, n = qa_sample)
  }
  if (missing(maxlength)) {
    maxlength <- max(max(fqa$Cycle), max(rqa$Cycle))
    message(paste0("Maxlength is emply, setting to ",
                   maxlength, " bp"))
  }
  if (missing(minlength) & !is.null(rev)) {
    overlap <- unlist(mapply(n_overlap, fwd, rev, sample = overlap_sample))
    ux <- unique(overlap)
    tab <- tabulate(match(overlap, ux))
    overlap_mode <- ux[tab == max(tab)]
    readlength <- min(max(fqa$Cycle), max(rqa$Cycle))

    # Min-trunclength
    minlength <- readlength - (overlap_mode/2 - minoverlap)
    message(paste0("Minlength is empty, estimated minimum trunlen of ", minlength, "bp to maintain ", minoverlap, "bp overlap between reads"))
  }  else if (missing(minlength) & is.null(rev)) {
    minlength <- 0
    warning(paste0("Minlength is empty and has been set to ",
                   minlength, " bp"))
  }
  fwd_trunc <- suppressWarnings(fqa %>% filter(QMean < threshold) %>% pull(Cycle) %>%
                                  min())
  rev_trunc <- suppressWarnings(rqa %>% filter(QMean < threshold) %>% pull(Cycle) %>%
                                  min())
  if (fwd_trunc > maxlength) {
    fwd_trunc <- maxlength
  }
  if (rev_trunc > maxlength) {
    rev_trunc <- maxlength
  }
  if (fwd_trunc < minlength) {
    fwd_trunc <- minlength
  }
  if (rev_trunc < minlength) {
    rev_trunc <- minlength
  }
  out <- c(fwd_trunc, rev_trunc)
  return(out)
}


#' Minimum overlap between a query and reference
#'
#' @param query
#' @param ref
#'
#' @return
#' @export
#'
#' @examples
align_overlap <- function(query, ref) {
  al <- dada2::nwalign(query, ref)
  fwd_gaps <- length(gregexpr("-", al[[1]])[[1]])
  rev_gaps <- length(gregexpr("-", al[[2]])[[1]])

  overlap <- nchar(al[[1]]) - (fwd_gaps + rev_gaps)
  return(overlap)
}


#' Number of bases overlapping between forward and reverse reads
#'
#' @param fwd
#' @param rev
#' @param sample
#'
#' @return
#' @export
#'
#' @examples
n_overlap <- function(fwd, rev, sample = 100) {

  #Sample reads
  fF <- FastqStreamer(fwd, n=sample)
  on.exit(close(fF))
  set.seed(123L);fqF <- yield(fF)

  fR <- FastqStreamer(rev, n=sample)
  on.exit(close(fR), add=TRUE)
  set.seed(123L);fqR <- yield(fR)

  #Nfilter
  keep <- nFilter()(fqF) & nFilter()(fqR)
  fqF <- fqF[keep]
  fqR <- fqR[keep]

  #Check ids match
  idF <- trimTails(ShortRead::id(fqF), 1, " ")
  idR <- trimTails(ShortRead::id(fqR), 1, " ")

  if(!all(idF == idR)){stop("Error: paired end reads dont match")}

  #Get overlap
  overlap <- mapply(align_overlap, as.character(sread(fqF)),  rc(as.character(sread(fqR))))
  overlap <- as.numeric(overlap)
  return(overlap)
}

