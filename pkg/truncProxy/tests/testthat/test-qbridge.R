make_bridge_inputs <- function(n = 80L) {
  dat <- make_example_data(n = n)
  q_Z <- cbind(1, dat$W1, dat$Z)
  q_W <- cbind(1, dat$W2, dat$Z)
  v <- sort(unique(dat$Q))
  qtime <- v[min(dat$TT) <= v & v <= max(dat$TT)]
  list(dat = dat, q_Z = q_Z, q_W = q_W, qtime = qtime)
}

test_that("qbridge validates input dimensions", {
  inp <- make_bridge_inputs()
  dat <- inp$dat

  expect_error(
    truncProxy:::qbridge(TT = dat$TT, Q = dat$Q[-1], qtime = inp$qtime,
            Z = inp$q_Z, W = inp$q_W, weights = rep(1, nrow(dat))),
    "Length of Q"
  )
  expect_error(
    truncProxy:::qbridge(TT = dat$TT, Q = dat$Q, qtime = inp$qtime,
            Z = inp$q_Z[-1, ], W = inp$q_W, weights = rep(1, nrow(dat))),
    "rows in Z"
  )
  expect_error(
    truncProxy:::qbridge(TT = dat$TT, Q = dat$Q, qtime = inp$qtime,
            Z = inp$q_Z, W = inp$q_W, weights = rep(1, nrow(dat) - 1)),
    "Length of weights"
  )
})

test_that("qbridge_IPCW validates the IPCW weight matrix dimensions", {
  inp <- make_bridge_inputs()
  dat <- inp$dat
  bad_ipcw <- matrix(1, nrow = nrow(dat), ncol = length(inp$qtime) - 1)

  expect_error(
    truncProxy:::qbridge_IPCW(TT = dat$TT, Q = dat$Q, qtime = inp$qtime,
                 Z = inp$q_Z, W = inp$q_W, IPCW_weights = bad_ipcw,
                 weights = rep(1, nrow(dat))),
    "IPCW_weights"
  )
})

test_that("qbridge reports zero pseudo-inverse failures on well-behaved data", {
  inp <- make_bridge_inputs()
  dat <- inp$dat

  qfit <- truncProxy:::qbridge(TT = dat$TT, Q = dat$Q, qtime = inp$qtime,
                  Z = inp$q_Z, W = inp$q_W, weights = rep(1, nrow(dat)))
  expect_identical(qfit$n_pinv_fail, 0L)
})

test_that("qbridge_IPCW with unit IPCW weights matches qbridge", {
  inp <- make_bridge_inputs()
  dat <- inp$dat
  w <- rep(1, nrow(dat))
  unit_ipcw <- matrix(1, nrow = nrow(dat), ncol = length(inp$qtime))

  fit_plain <- truncProxy:::qbridge(TT = dat$TT, Q = dat$Q, qtime = inp$qtime,
                       Z = inp$q_Z, W = inp$q_W, weights = w)
  fit_ipcw <- truncProxy:::qbridge_IPCW(TT = dat$TT, Q = dat$Q, qtime = inp$qtime,
                           Z = inp$q_Z, W = inp$q_W, IPCW_weights = unit_ipcw,
                           weights = w)

  expect_equal(fit_ipcw$B_int, fit_plain$B_int, tolerance = 1e-12)
  expect_equal(fit_ipcw$dB, fit_plain$dB, tolerance = 1e-12)
})

test_that("B_int_eval clamps evaluation times outside the qtime range", {
  inp <- make_bridge_inputs()
  dat <- inp$dat

  qfit <- truncProxy:::qbridge(TT = dat$TT, Q = dat$Q, qtime = inp$qtime,
                  Z = inp$q_Z, W = inp$q_W, weights = rep(1, nrow(dat)))
  B_int <- qfit$B_int
  qtime <- qfit$qtime

  times <- rep(max(qtime), nrow(dat))
  times[1] <- max(qtime) * 2      # beyond the last jump
  times[2] <- 0                   # before the first jump

  vals <- truncProxy:::B_int_eval(times, B_int, qtime)
  expect_length(vals, nrow(dat))
  expect_false(anyNA(vals))
  expect_identical(vals[1], B_int[1, ncol(B_int)])
  expect_identical(vals[2], B_int[2, 1])
})
