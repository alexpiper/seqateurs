# Transformations ---------------------------------------------------------

#' Convert phyloseq table to proportions
#' #'
#' @param x An OTU Table from Phyloseq
#' @param thresh A minimum abundance threshold to filter taxa
#' @param na_rm Whether NA's should be removed
#' @param ...
#'
#' @return
#'
#' @examples
proportions <- function(x, thresh = NA, na_rm = FALSE, ...) {
  xprop <- (x / sum(x))
  xprop[xprop <= thresh] <- NA
  xprop2 <- (xprop / sum(xprop, na.rm = na_rm))
  return(xprop2)
}

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
