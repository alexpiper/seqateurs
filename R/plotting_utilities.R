# Get Qual Stats ----------------------------------------------------------

#' Get qual stats
#'
#' @param input
#' @param n The number of reads to sample
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
    fq <- FastqFile(f)
    rc <- countFastq(f)$records
    close(fq)
    coll <- QACollate(QAFastqSource(f, n=n), QAQualityByCycle())
    srqa <- qa2(coll, BPPARAM=SerialParam(), verbose=FALSE)
    df <- data.frame(srqa$QAQualityByCycle@values)
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
    } else {
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
#' @param n The number of reads to sample
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
#' @import ShortRead
#' @import magrittr
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import ggplot2
#'
#' @examples
plot_lengths <- function(dir, aggregate = FALSE, sample = 1e6){
  Forward <-  sort(list.files(dir, pattern="_R1_", full.names = TRUE))
  Reverse <-  sort(list.files(dir, pattern="_R2_", full.names = TRUE))

  lengths <- vector("list", length=length(Forward))
  for(i in 1:length(Forward)){
    #Forwrard reads
    sampler <- ShortRead::FastqSampler(Forward[[i]], n=sample )
    Ffq <- as.data.frame(table(width(ShortRead::yield(sampler))), stringsAsFactors=FALSE) %>%
      magrittr::set_colnames(c("Readlength", "Forward")) %>%
      dplyr::mutate_if(is.character, as.integer)
    #Reverse reads
    sampler <- ShortRead::FastqSampler(Reverse[[i]], n=sample )
    Rfq <- as.data.frame(table(width(ShortRead::yield(sampler))), stringsAsFactors=FALSE) %>%
      magrittr::set_colnames(c("Readlength", "Reverse"))%>%
      dplyr::mutate_if(is.character, as.integer)

    lengths[[i]] <- tibble::tibble(Readlength=seq(1, max(max(Ffq$Readlength), max(Rfq$Readlength)),1)) %>%
      dplyr::left_join(Ffq, by="Readlength") %>%
      dplyr::left_join(Rfq, by="Readlength")
    lengths[[i]]$source <- Forward[i]
  }

  p <- do.call(rbind, lengths) %>%
    dplyr::group_by(Readlength, source) %>%
    dplyr::summarise(Forward = sum(Forward), Reverse = sum(Reverse)) %>%
    tidyr::pivot_longer(cols=3:4,
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
#' @param fwd A path to a fastq file containing forward reads
#' @param rev A path to a fastq file containing reverse reads
#' @param threshold A Quality threshold used to signal the threshold to truncate the reads at
#' @param maxlength The maximum length for the reads to be trimmed at. Suggested length of amplicon if known.
#' @param minlength The minimum length for reads to be trimmed at to allow merging.
#' @param qa_sample The number of reads to sample for calculating qualities.
#' @param overlap_sample The number of reads to sample for aligning and estimating minlength
#' @param quiet Whether progress should be printed to console
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
#' @importFrom dada2 nwalign
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
#' @param fwd A path to a fastq file containing forward reads
#' @param rev A path to a fastq file containing reverse reads
#' @param sample The number of reads to sample from the fastq files
#'
#' @return
#' @export
#' @import ShortRead
#' @import Biostrings
#'
#' @examples
n_overlap <- function(fwd, rev, sample = 100) {

  #Sample reads
  fF <- ShortRead::FastqStreamer(fwd, n=sample)
  on.exit(close(fF))
  set.seed(123L);fqF <- ShortRead::yield(fF)

  fR <- ShortRead::FastqStreamer(rev, n=sample)
  on.exit(close(fR), add=TRUE)
  set.seed(123L);fqR <- ShortRead::yield(fR)

  #Get reads where both contain no N's
  keep <- (Biostrings::alphabetFrequency(sread(fqF), baseOnly=TRUE)[,"other"] <= 0
           ) & (Biostrings::alphabetFrequency(sread(fqR), baseOnly=TRUE)[,"other"] <= 0)

  fqF <- fqF[keep]
  fqR <- fqR[keep]

  #Check ids match
  idF <- ShortRead::trimTails(ShortRead::id(fqF), 1, " ")
  idR <- ShortRead::trimTails(ShortRead::id(fqR), 1, " ")

  if(!all(idF == idR)){stop("Error: paired end reads dont match")}

  #Get overlap
  overlap <- mapply(align_overlap, as.character(sread(fqF)),  rc(as.character(sread(fqR))))
  overlap <- as.numeric(overlap)
  return(overlap)
}

#' Function to make ordination plots with margin boxplots
#'
#' @param phyloseq a phyloseq object
#' @param Dissimilarity index, partial match to "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger", "aitchison", or "robust.aitchison". Parsed to vegan::vegdist()
#' @param factor1 factor used for fill aesthetic
#' @param factor2 factor used for fill aesthetic
#' @param ellipse Which factor to use for ellipses. Value must be 'factor1', 'factor2', 'none' or 'both'
#' @param legend Which factor to plot legends for. Value must be 'factor1', 'factor2', 'none' or 'both'
#' @param marginplot Which factor to plot margin plots for. Value must be 'factor1', 'factor2', 'none' or 'both'
#' @param margintext Whether text should be plotted for margin plots
#' @param palette A colour palette compatible with RColorBrewer
#' @param title Add a title to the plot
#' @param formula Formula parsed to vegan::adonis2 for adding stats annotations.
#'
#' @return
#' @export
#' @import vegan
#' @import phyloseq
#' @import ggplot2
#' @import patchwork
#' @import RColorBrewer
#'
#' @examples
#
plot_ordination_margins <- function(ps, method="jac", factor1, factor2,
                                    ellipse="factor1", legend="both",marginplot="factor1",
                                    margintext=TRUE, palette="Set2", title=NULL,
                                    formula=NULL){
  # Get otutab and filter for empty rows and columns
  otutab <- phyloseq::otu_table(ps) %>%
    as("matrix")
  otutab <- otutab[rowSums(otutab) >0,colSums(otutab) > 0]

  # Create distance matrix
  if (method == "aitchison"){
    distmat <- as.matrix(vegan::vegdist(otutab, method="aitchison", pseudocount=1)) # replace zeros with a pseudocount of 1

  } else {
    distmat <- as.matrix(vegan::vegdist(otutab, method=method,binary = T))
  }

  # Ordinate distance matrix using pcoa
  ord <- phyloseq::ordinate(ps, 'PCoA', distance=as.dist(distmat))

  # Get ordtable
  ordtable <- phyloseq::plot_ordination(ps, ord, justDF = TRUE)

  # Plot main ordination
  p1 <- ordtable %>%
    ggplot(aes(x = Axis.1, y = Axis.2, fill=!!enquo(factor1), shape=!!enquo(factor2), colour=!!enquo(factor1)))

  if(ellipse =="factor1"){
    p1 <- p1 +
      stat_ellipse(type="norm", alpha=1, aes(group=!!enquo(factor1)))
  } else if(ellipse =="factor2"){
    p1 <- p1 +
      stat_ellipse(type="norm", alpha=1, aes(group=!!enquo(factor2)))
  } else if(ellipse =="both"){
    p1 <- p1 +
      stat_ellipse(type="norm", alpha=1, aes(group=!!enquo(factor1)))+
      stat_ellipse(type="norm", alpha=1, aes(group=!!enquo(factor2)))
  } else if(ellipse == "none"){
    p1 <- p1
  } else {
    stop("value of 'ellipse' must be 'factor1', 'factor2', 'none' or 'both'")
  }
  p1 <- p1 +
    geom_hline(yintercept = 0, linetype=2, alpha=0.5) +
    geom_vline(xintercept = 0, linetype=2, alpha=0.5) +
    geom_point(size=3, alpha=1, colour="black") +
    scale_shape_manual(values=c(21,22,24))+
    theme_bw() +
    theme(legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
    )

  if(legend =="factor1"){
    p1 <- p1 +
      guides(shape=FALSE)
  } else if(legend =="factor2"){
    p1 <- p1 +
      guides(fill = FALSE)
  } else if(legend =="both"){
    p1 <- p1
  } else if(legend == "none"){
    p1 <- p1+
      guides(fill = FALSE, colour = FALSE,  shape=FALSE)
  } else {
    stop("value of 'legend' must be 'factor1', 'factor2', 'none' or 'both'")
  }

  # Add percentage variation explained
  eigvec <- ord$values$Eigenvalues
  fracvar <- eigvec[1:2]/sum(eigvec)
  percvar <- round(100 * fracvar, 1)
  strivar <- as(c(p1$label$x, p1$label$y), "character")
  strivar <- paste0(strivar, "   [", percvar, "%]")
  p1 <- p1 + xlab(strivar[1]) + ylab(strivar[2])

  # Get limits of ordination to define ranges of margin plots
  xlim <- layer_scales(p1)$x$range$range
  ylim <- layer_scales(p1)$y$range$range

  # Create margin plots
  p1_axis1 <- ordtable %>%
    ggplot(aes(y = !!enquo(factor1), x = Axis.1, fill=!!enquo(factor1), group=!!enquo(factor1))) +
    geom_boxplot(outlier.shape = 21)+
    scale_x_continuous(limits = xlim)+
    theme_bw()+
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
    )

  p1_axis2 <- ordtable %>%
    ggplot(aes(x = !!enquo(factor1), y = Axis.2, fill=!!enquo(factor1))) +
    geom_boxplot(outlier.shape = 21)+
    scale_y_continuous(limits = ylim)+
    theme_bw()+
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
    )

  p2_axis1 <- ordtable %>%
    ggplot(aes(y = !!enquo(factor2), x = Axis.1, group=!!enquo(factor2))) +
    geom_boxplot(outlier.shape = 21)+
    scale_x_continuous(limits = xlim)+
    theme_bw()+
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
    )

  p2_axis2 <- ordtable %>%
    ggplot(aes(x = !!enquo(factor2), y = Axis.2)) +
    geom_boxplot(outlier.shape = 21)+
    scale_y_continuous(limits = ylim)+
    theme_bw()+
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
    )

  # remove margin text
  if(margintext){
    p1_axis1 <- p1_axis1 +
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank()
      )
    p1_axis2 <- p1_axis2 +
      theme(
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank()
      )
    p2_axis1 <- p2_axis1 +
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank()
      )
    p2_axis2 <- p2_axis2 +
      theme(
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank()
      )
  } else if(!margintext){
    p1_axis1 <- p1_axis1 +
      theme(
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()
      )
    p1_axis2 <- p1_axis2 +
      theme(
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()
      )
    p2_axis1 <- p2_axis1 +
      theme(
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()
      )
    p2_axis2 <- p2_axis2 +
      theme(
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()
      )
  }

  # Add colour palettes
  if(palette == "viridis"){
    p1 <- p1 +
      scale_fill_viridis_d()+
      scale_colour_viridis_d()
    p1_axis1 <- p1_axis1 +
      scale_fill_viridis_d()
    p1_axis2 <- p1_axis2 +
      scale_fill_viridis_d()
  }else {
    p1 <- p1 +
      scale_fill_brewer(palette=palette) +
      scale_colour_brewer(palette=palette)
    p1_axis1 <- p1_axis1 +
      scale_fill_brewer(palette=palette)
    p1_axis2 <- p1_axis2 +
      scale_fill_brewer(palette=palette)
  }

  # Add main plot title
  if(is.character(title) & marginplot =="factor1"){
    p1_axis1 <- p1_axis1 +
      ggtitle(label = title)
  } else if(is.character(title) & marginplot %in% c("factor2", "both")){
    p2_axis1 <- p2_axis1 +
      ggtitle(label = title)
  } else if(is.character(title) & marginplot == "none"){
    p1 <- p1 +
      ggtitle(label = title)
  }

  if(!is.null(formula)){
    # Run adonis
    metadata <- sample_data(phyloseq) %>%
      as("data.frame")

    #### permanova tests for differences in beta diversity
    y <- as.dist(distmat[metadata$sample_name, metadata$sample_name])

    # test for differences driven by endophyte, taking into account flow cell and month as variables
    res <- vegan::adonis2(formula, method="euclidean", data=metadata, permutations = 999)
    blurb <- broom::tidy(res) %>%
      filter(!is.na(p.value))%>%
      mutate(disp = paste0(term,"(",df,") R2: ", round(R2,3), " p: ", p.value)) %>%
      pull(disp) %>%
      paste(collapse="\n")
    blurb <- paste0("Adonis\n",blurb)
    text_stat <- ggplot() +
      annotate("text", x = 0, y = 25, label = blurb) +
      theme_void()
  } else {
    text_stat <- ggplot() +
      annotate("text", x = 0, y = 25, label = "") +
      theme_void()
  }

  if(marginplot =="factor1"){
    # Set up custom layout for patchwork
    design <- "
      12
      34
    "
    out <- p1_axis1 + text_stat + p1 + p1_axis2 +plot_layout(design = design, widths = c(6,1), heights=c(1,6))
  } else if(marginplot =="factor2"){
    # Set up custom layout for patchwork
    design <- "
      12
      34
    "
    out <- p2_axis1 + text_stat + p1 + p2_axis2 + plot_layout(design = design, widths = c(6,1), heights=c(1,6))
  } else if(marginplot =="both"){
    # Set up custom layout for patchwork
    design <- "
      12#
      3##
      456
    "
    out <- p2_axis1 + text_stat + p1_axis1 + p1 + p1_axis2 + p2_axis2 + plot_layout(design = design, widths = c(6,1,1), heights=c(1,1,6))  } else if(marginplot == "none"){
      out <- p1
    } else {
      stop("value of 'legend' must be 'factor1', 'factor2', 'none' or 'both'")
    }
  return(out)
}

