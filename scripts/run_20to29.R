## Deterministic run: NRT vs SOC, cohort starting age 20-29 (start_age = 25), reproducing
## "Markov model (NRT 20to29)" and "Markov model (SOC 20to29)" from the source workbook.
##
## Usage:  Rscript scripts/run_20to29.R   (run from the repo root)

here <- if (interactive()) getwd() else normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
if (!file.exists(file.path(here, "R", "run_model.R"))) here <- getwd() # fallback when sourced directly

for (f in list.files(file.path(here, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)

params <- load_parameters(file.path(here, "data", "parameters.csv"))

start_age <- 25 # start_page!L10 in the workbook (the 20-29 cohort's representative starting age)
result <- run_comparison(params, start_age)

fmt <- function(x) formatC(x, format = "f", digits = 2, big.mark = ",")

cat(sprintf("Cohort starting age: %.0f (age band 20-29)\n", start_age))
cat(sprintf("Cohort size (this age band): %s\n\n", fmt(get_param(params, "Cohort") * get_param(params, "Cohort_20_30"))))

cat("                         NRT              SOC\n")
cat(sprintf("Total discounted cost   %14s   %14s\n", fmt(result$NRT$total_cost), fmt(result$SOC$total_cost)))
cat(sprintf("Total discounted QALYs  %14s   %14s\n", fmt(result$NRT$total_qaly), fmt(result$SOC$total_qaly)))
cat(sprintf("Total discounted LYs    %14s   %14s\n", fmt(result$NRT$total_ly), fmt(result$SOC$total_ly)))
cat(sprintf("Total deaths (cohort)   %14s   %14s\n\n", fmt(result$NRT$total_deaths), fmt(result$SOC$total_deaths)))

cat(sprintf("Incremental cost:  %s\n", fmt(result$incremental_cost)))
cat(sprintf("Incremental QALYs: %s\n", fmt(result$incremental_qaly)))
cat(sprintf("ICER (cost per QALY gained): %s\n", fmt(result$icer)))
