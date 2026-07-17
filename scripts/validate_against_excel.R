## Compares this R model's deterministic lifetime totals against the *cached* values already
## computed inside the source workbook (Analysis sheet, "Aggregated Result" table, rows 6-13 /
## 26-34: lifetime discounted-and-half-cycle-corrected cost/QALY/LY by age band and arm).
##
## This is a structural rebuild, not a cell-by-cell transliteration (see README.md), so exact
## agreement isn't expected -- but the model should land within a few percent of the workbook's
## own numbers if the state-transition, cost and QALY logic have been carried over correctly.
##
## Usage: Rscript scripts/validate_against_excel.R

here <- if (interactive()) getwd() else normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
if (!file.exists(file.path(here, "R", "run_model.R"))) here <- getwd()
for (f in list.files(file.path(here, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
params <- load_parameters(file.path(here, "data", "parameters.csv"))

## benchmarks: workbook "Analysis" sheet, rows 9-13 (SOC/NRT cost) and columns H/J (SOC/NRT QALY),
## L/N (SOC/NRT LY) -- read directly from the cached (data_only) workbook, 06.07.26_Final_PSA_and_USA.xlsm
benchmarks <- data.frame(
  band = c("20_29", "30_39", "40_49", "50_59", "Over60"),
  start_age = c(25, 35, 45, 55, 65),
  cost_soc = c(885436478.69, 3306577954.89, 4544652189.38, 3581068122.28, 1135600366.35),
  cost_nrt = c(891276302.17, 3288844687.35, NA, NA, NA),
  qaly_soc = c(4896.821108, 16539.046006, 19618.850603, 12578.114625, 2868.548322),
  qaly_nrt = c(4931.539838, NA, NA, NA, NA)
)

cat(sprintf("%-8s %14s %14s %10s | %10s %10s %8s\n",
            "band", "cost (R)", "cost (xlsx)", "diff%", "qaly (R)", "qaly (xlsx)", "diff%"))
for (k in seq_len(nrow(benchmarks))) {
  b <- benchmarks[k, ]
  soc <- run_markov(params, "SOC", b$start_age)
  cat(sprintf(
    "%-8s %14.0f %14.0f %9.2f%% | %10.2f %10.2f %7.2f%%\n",
    b$band, soc$total_cost, b$cost_soc, 100 * (soc$total_cost / b$cost_soc - 1),
    soc$total_qaly, b$qaly_soc, 100 * (soc$total_qaly / b$qaly_soc - 1)
  ))
}

cat("\n20-29 NRT arm:\n")
nrt <- run_markov(params, "NRT", 25)
cat(sprintf(
  "  cost (R) %.0f vs (xlsx) %.0f  [%.2f%%]\n", nrt$total_cost, benchmarks$cost_nrt[1],
  100 * (nrt$total_cost / benchmarks$cost_nrt[1] - 1)
))
cat(sprintf(
  "  qaly (R) %.2f vs (xlsx) %.2f  [%.2f%%]\n", nrt$total_qaly, benchmarks$qaly_nrt[1],
  100 * (nrt$total_qaly / benchmarks$qaly_nrt[1] - 1)
))
