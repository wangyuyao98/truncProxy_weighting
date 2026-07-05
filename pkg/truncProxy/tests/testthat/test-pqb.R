expected_example_results <- c(
  pqb = 0.46829407648884291,
  pqb_ipcw_case = 0.41083109830851139,
  pqb_ipcw_time_varying = 0.46108075911803298
)

test_that("PQB_estimator matches a fixed regression target", {
  dat <- make_example_data()
  nu <- function(t) as.numeric(t > 1)

  est <- PQB_estimator(
    nu = nu,
    dat = dat,
    time.name = "TT",
    Q.name = "Q",
    W1.name = "W1",
    W2.name = "W2",
    Z.name = "Z"
  )

  expect_equal(est, unname(expected_example_results["pqb"]), tolerance = 1e-8)
})

test_that("PQB_IPCW_estimator with case weights matches a fixed regression target", {
  dat <- make_example_data()
  nu <- function(t) as.numeric(t > 1)

  est <- PQB_IPCW_estimator(
    nu = nu,
    t0 = 1,
    dat = dat,
    time.name = "X",
    Q.name = "Q",
    event.name = "delta",
    W1.name = "W1",
    W2.name = "W2",
    Z.name = "Z",
    IPCW_time_varying = FALSE
  )

  expect_equal(est, unname(expected_example_results["pqb_ipcw_case"]), tolerance = 1e-8)
})

test_that("PQB_IPCW_estimator with time-varying weights matches a fixed regression target", {
  dat <- make_example_data()
  nu <- function(t) as.numeric(t > 1)

  est <- PQB_IPCW_estimator(
    nu = nu,
    t0 = 1,
    dat = dat,
    time.name = "X",
    Q.name = "Q",
    event.name = "delta",
    W1.name = "W1",
    W2.name = "W2",
    Z.name = "Z",
    IPCW_time_varying = TRUE
  )

  expect_equal(est, unname(expected_example_results["pqb_ipcw_time_varying"]), tolerance = 1e-8)
})
