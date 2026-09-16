# Assemble Table 1 of the paper (LaTeX) from the saved main_simu3.R summaries.
#
# Expects, in results/:
#   simu3_table1_n500_tvTRUE.csv    simu3_table1_n1000_tvTRUE.csv
#   simu3_table1_n500_tvFALSE.csv   simu3_table1_n1000_tvFALSE.csv
# (produce them with SIMU_N / SIMU_TV / SIMU_OUT, see README). Missing files
# yield "--" placeholders so a partial table can be generated at any point.
#
#   Rscript dev/make_table1.R            # prints LaTeX to stdout
#   Rscript dev/make_table1.R > results/table1.tex

read_run <- function(n, tv) {
  f <- sprintf("results/simu3_table1_n%d_tv%s.csv", n, if (tv) "TRUE" else "FALSE")
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f); rownames(d) <- d$method; d
}
runs <- list(`500` = list(tv = read_run(500, TRUE), cw = read_run(500, FALSE)),
             `1000` = list(tv = read_run(1000, TRUE), cw = read_run(1000, FALSE)))

# The estimators that do not use IPCW weights must agree between the tv and cw
# runs (same seeds); check that before trusting the table.
for (n in names(runs)) {
  a <- runs[[n]]$tv; b <- runs[[n]]$cw
  if (!is.null(a) && !is.null(b)) {
    common <- c("PQB_noC", "IPQWo", "naive", "PL", "KM")
    d <- max(abs(a[common, "bias"] - b[common, "bias"]))
    if (d > 1e-10) warning(sprintf("n=%s: tv and cw runs disagree on IPCW-free estimators (max %.2e)", n, d))
  }
}

# Table rows: paper label -> (which run, method name in summary_df, commented?)
rows <- list(
  list(lab = "PQB",       run = "tv", m = "PQB"),
  list(lab = "IPQW",      run = "tv", m = "IPQW"),
  list(lab = "IPQW2",     run = "tv", m = "IPQW2",   comment = TRUE),
  list(lab = "PQB-cw",    run = "cw", m = "PQB"),
  list(lab = "IPQW-cw",   run = "cw", m = "IPQW"),
  list(lab = "IPQW2-cw",  run = "cw", m = "IPQW2",   comment = TRUE),
  list(lab = "PL",        run = "tv", m = "PL"),
  list(lab = "KM",        run = "tv", m = "KM"),
  list(lab = "naive",     run = "tv", m = "naive"),
  list(midrule = TRUE),
  list(lab = "IPQW-U",    run = "tv", m = "IPQWu"),
  list(lab = "IPQW-U-cw", run = "cw", m = "IPQWu"),
  list(lab = "PQB-noC",   run = "tv", m = "PQB_noC", comment = TRUE),
  list(lab = "IPQW-o",    run = "tv", m = "IPQWo")
)

cell <- function(run, m) {
  if (is.null(run) || !(m %in% rownames(run))) return(rep("--", 4))
  r <- run[m, ]
  c(sprintf("%.4f", r$bias), sprintf("%.4f", r$empSE), sprintf("%.4f", r$bootSE), sprintf("%.3f", r$coverage))
}

cat("\\begin{table}[ht]\n\\centering\n")
cat("\\caption{Simulation results for different estimators under the case with right censoring. Each observed data set has sample size 500 or 1000, and 500 data sets are simulated for each sample size.}\n")
cat("\\label{tab:simu2_cens}\n")
cat("\\begin{tabular}{lrcccrccc}\n  \\toprule\n")
cat("  & \\multicolumn{4}{c}{$n = 500$} & \\multicolumn{4}{c}{$n = 1000$} \\\\\n")
cat("  \\cmidrule(lr){2-5} \\cmidrule(lr){6-9}\n")
cat("  Method & Bias & SD & bootSE & CP & Bias & SD & bootSE & CP \\\\\n  \\midrule\n")
for (r in rows) {
  if (isTRUE(r$midrule)) { cat("  \\midrule\n"); next }
  cells <- c(cell(runs[["500"]][[r$run]], r$m), cell(runs[["1000"]][[r$run]], r$m))
  cat(sprintf("  %s%s & %s \\\\\n", if (isTRUE(r$comment)) "% " else "", r$lab, paste(cells, collapse = " & ")))
}
cat("  \\bottomrule\n\\end{tabular}\n\\end{table}\n")
