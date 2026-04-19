#' truncProxy: Proximal Weighting Estimation for Dependent Left Truncation
#'
#' Implements proximal weighting estimators for the expectation of an
#' arbitrarily transformed event time under dependent left truncation, with
#' optional inverse probability of censoring weighting to handle right
#' censoring. The methods leverage proxy variables to handle dependent left
#' truncation in settings where dependence-inducing factors are not fully
#' observed.
#'
#' @keywords internal
#' @useDynLib truncProxy, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
