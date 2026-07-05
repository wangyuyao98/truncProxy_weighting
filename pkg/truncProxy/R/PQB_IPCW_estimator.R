#' IPCW-Adjusted Proximal Weighting Estimation under Dependent Left Truncation and Random Right Censoring
#'
#' Computes the IPCW-adjusted proximal weighting estimator for the expectation
#' of an arbitrarily transformed event time under dependent left truncation
#' and random right censoring.
#'
#' @param nu A user-supplied transformation of the event time. For example,
#'   `function(t) as.numeric(t > tau)` results in estimating a survival
#'   probability beyond time `tau`.
#' @param t0 A cutoff such that `nu(t) = nu(min(t, t0))`. For example, for
#'   survival probability at time `tau`, one can set `t0 = tau`.
#' @param dat A data frame containing the observed data.
#' @param time.name Name of the observed event or censoring time variable.
#' @param Q.name Name of the left truncation time variable.
#' @param event.name Name of the event indicator variable. Set to `NULL` when
#'   there is no right censoring.
#' @param W1.name Name of the truncation proxies.
#' @param W2.name Name of the event time proxies.
#' @param Z.name Name of the measured covariates that are directly associated with both the event time and the truncation time.
#' @param weights Optional nonnegative case weights. Defaults to equal weights.
#' @param trim.min Lower bound used to stabilize the estimated bridge weights.
#' @param trim.max Upper bound used to stabilize the estimated bridge weights.
#' @param IPCW_time_varying Logical; if `TRUE`, use time-varying IPCW weights in
#'   the bridge estimation step.
#' @param trim.IPCW Lower bound used to stabilize the denominators of IPCW weights.
#'
#' @return A numeric scalar containing the IPCW-adjusted proximal weighting
#'   estimator.
#' @export
PQB_IPCW_estimator <- function(nu, t0, dat, time.name, Q.name, event.name = NULL,
                               W1.name, W2.name, Z.name,
                               weights = rep(1, nrow(dat)),
                               trim.min = 1e-7, trim.max = 1e7,
                               IPCW_time_varying = FALSE, trim.IPCW = 1e-7) {
  names_used <- c(time.name, Q.name, W1.name, W2.name, Z.name)
  if (sum(names_used == "weights") > 0) {
    stop("The names of the variables cannot be 'weights', which is used in the middle of computation")
  } else if (sum(names_used == "IPCW_weights") > 0) {
    stop("The names of the variables cannot be 'IPCW_weights', which is used in the middle of computation")
  }

  if (is.null(event.name)) {
    if (sum(names_used == "delta1") > 0) {
      stop("The names of the variables cannot be 'delta1', which is used in the middle of computation")
    }
    event.name <- "delta1"
    dat$delta1 <- rep(1, nrow(dat))
  }

  Y <- dat[, time.name] - dat[, Q.name]
  delta_C <- 1 - dat[, event.name]
  fitD <- survival::survfit(
    survival::Surv(Y, delta_C) ~ 1,
    type = "kaplan-meier",
    weights = weights
  )
  stepfunD <- stats::stepfun(fitD$time, c(1, fitD$surv))

  if (IPCW_time_varying) {
    X <- dat[, time.name]
    Q <- dat[, Q.name]
    q_Z <- as.matrix(cbind(1, dat[, W1.name], dat[, Z.name]))
    colnames(q_Z) <- c("1", W1.name, Z.name)
    q_W <- as.matrix(cbind(1, dat[, W2.name], dat[, Z.name]))
    colnames(q_W) <- c("1", W2.name, Z.name)

    v <- sort(unique(Q))
    qtime.id <- (min(X) <= v & v <= max(X))
    qtime <- v[qtime.id]

    IPCW_weights_mx <- matrix(nrow = nrow(dat), ncol = length(qtime))
    delta_t_mx <- matrix(nrow = nrow(dat), ncol = length(qtime))
    for (j in seq_along(qtime)) {
      temp_time <- qtime[j]
      delta_t_mx[, j] <- as.numeric(dat[, event.name] == 1 | X > temp_time)
      IPCW_weights_mx[, j] <- delta_t_mx[, j] / pmax(stepfunD(pmin(X, temp_time) - Q), trim.IPCW)
    }
    if (anyNA(IPCW_weights_mx)) {
      stop("Time-varying IPCW weights contain NA values; check the censoring KM fit.")
    }

    qfit <- qbridge_IPCW(
      TT = X,
      Q = Q,
      qtime = qtime,
      Z = q_Z,
      W = q_W,
      IPCW_weights = IPCW_weights_mx,
      weights = weights
    )
    if (qfit$n_pinv_fail > 0) {
      warning(sprintf(
        "qbridge_IPCW: pseudo-inverse failed at %d of %d time points; the corresponding increments were skipped.",
        qfit$n_pinv_fail, length(qtime)
      ))
    }
    B_int <- qfit$B_int
    qtime <- qfit$qtime

    bT_hat <- exp(B_int_eval(X, B_int, qtime))
    bT_hat <- pmin(pmax(bT_hat, trim.min), trim.max)

    Q_max <- max(Q)
    tau_delta <- max(Q_max, t0)
    delta_t <- as.numeric(dat[, event.name] == 1 | X > tau_delta)
    IPCW_weights <- delta_t / pmax(stepfunD(pmin(X, tau_delta) - Q), trim.IPCW)
    if (anyNA(IPCW_weights)) {
      stop("IPCW weights contain NA values; check the censoring KM fit.")
    }
    numer <- mean(weights * IPCW_weights * bT_hat * nu(X))
    denom <- mean(weights * IPCW_weights * bT_hat)

    est <- numer / denom
  } else {
    X <- dat[, time.name]
    Q <- dat[, Q.name]

    Q_max <- max(Q)
    tau_delta <- max(Q_max, t0)
    delta_t <- as.numeric(dat[, event.name] == 1 | X > tau_delta)
    IPCW_weights <- delta_t / pmax(stepfunD(pmin(X, tau_delta) - Q), trim.IPCW)
    if (anyNA(IPCW_weights)) {
      stop("IPCW weights contain NA values; check the censoring KM fit.")
    }

    est <- PQB_estimator(
      nu = nu,
      dat = dat,
      time.name = time.name,
      Q.name = Q.name,
      W1.name = W1.name,
      W2.name = W2.name,
      Z.name = Z.name,
      weights = weights * IPCW_weights,
      trim.min = trim.min,
      trim.max = trim.max
    )
  }

  est
}
