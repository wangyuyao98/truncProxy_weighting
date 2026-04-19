#' Proximal Weighting Estimation Under Dependent Left Truncation
#'
#' Computes the proximal weighting estimator for the expectation of an
#' arbitrarily transformed event time under dependent left truncation.
#'
#' @param nu A user-supplied transformation of the event time. For example,
#'   `function(t) as.numeric(t > tau)` results in estimating a survival
#'   probability beyond time `tau`.
#' @param dat A data frame containing the observed data.
#' @param time.name Name of the event time variable.
#' @param Q.name Name of the left truncation time variable.
#' @param W1.name Name of the truncation proxies.
#' @param W2.name Name of the event time proxies.
#' @param Z.name Name of the measured covariates that are directly associated with both the event time and the truncation time.
#' @param weights Optional nonnegative case weights. Defaults to equal weights.
#' @param trim.min Lower bound used to stabilize the estimated bridge weights.
#' @param trim.max Upper bound used to stabilize the estimated bridge weights.
#'
#' @return A numeric scalar containing the proximal weighting estimator.
#' @export
PQB_estimator <- function(nu, dat, time.name, Q.name, W1.name, W2.name, Z.name,
                          weights = rep(1, nrow(dat)),
                          trim.min = 1e-7, trim.max = 1e7) {
  TT <- dat[, time.name]
  Q <- dat[, Q.name]
  q_Z <- as.matrix(cbind(1, dat[, W1.name], dat[, Z.name]))
  colnames(q_Z) <- c("1", W1.name, Z.name)
  q_W <- as.matrix(cbind(1, dat[, W2.name], dat[, Z.name]))
  colnames(q_W) <- c("1", W2.name, Z.name)

  v <- sort(unique(Q))
  qtime.id <- (min(TT) <= v & v <= max(TT))
  qtime <- v[qtime.id]

  qfit <- qbridge(TT = TT, Q = Q, qtime = qtime, Z = q_Z, W = q_W, weights = weights)
  B_int <- qfit$B_int
  qtime <- qfit$qtime

  bT_hat <- exp(B_int_eval(TT, B_int, qtime))
  bT_hat <- pmin(pmax(bT_hat, trim.min), trim.max)

  numer <- mean(weights * bT_hat * nu(TT))
  denom <- mean(weights * bT_hat)

  numer / denom
}
