# Evaluate the integrated bridge process at subject-specific time points.
B_int_eval <- function(time.eval, B_int, qtime) {
  if (length(time.eval) != nrow(B_int)) {
    stop("The number of time points does not equal the number of subjects!")
  }

  if (ncol(B_int) != (length(qtime) + 1)) {
    stop(
      paste(
        "ncol(B_int) should equal length(qtime) + 1, where the last column",
        "of B_int corresponds to the time in the initial condition."
      )
    )
  }

  id <- findInterval(time.eval, c(0, qtime, Inf))
  id <- pmin(pmax(id, 1), ncol(B_int))
  B_int[cbind(seq_along(time.eval), id)]
}
