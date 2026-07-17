## Parameter loading and access layer.
##
## Reads data/parameters.csv (extracted from the "Parameters" sheet of the source
## workbook: 06.07.26_Final_PSA_and_USA.xlsm). Column meanings mirror the Excel sheet:
##   mean        -> Parameters!D (deterministic value used as the point estimate)
##   se          -> Parameters!E (standard error, used for PSA draws)
##   distribution-> Parameters!F ("Log-normal", "Direchlet"/Dirichlet via Gamma, or blank/fixed)
##   lower/upper -> Parameters!H:I (alpha/lower or beta/upper bound, distribution-dependent)
##   ln_mean/ln_sd -> Parameters!J:K (precomputed log-normal mean/sd = log(mean), used directly
##                     for Log-normal draws so PSA reproduces the workbook's own parameterisation)
##
## Deterministic mode simply returns `mean`. PSA mode (not run by default - see below) draws
## from the distribution named in `distribution` using `se`/`ln_sd`/bounds, exactly mirroring
## the Excel "Probabilistic" column formulas (NORM.INV/EXP for Log-normal, GAMMA.INV for the
## Dirichlet blocks). The parameter table carries everything a PSA run needs; only the sampling
## step (draw_parameters()) needs to be turned on.

load_parameters <- function(path = file.path("data", "parameters.csv")) {
  params <- read.csv(path, stringsAsFactors = FALSE)
  params$distribution[is.na(params$distribution)] <- ""
  rownames(params) <- params$name
  params
}

## Deterministic accessor: mirrors Parameters!B (IF(probabilistic, C, D)) with probabilistic = 0.
get_param <- function(params, name) {
  if (!name %in% rownames(params)) {
    stop(sprintf("Unknown parameter: '%s'", name))
  }
  params[name, "mean"]
}

get_params <- function(params, names) {
  vapply(names, function(n) get_param(params, n), numeric(1))
}

## PSA draw for a single parameter, following the same distributional forms used in the
## Parameters sheet's "Probabilistic" column (col C). Not called by the deterministic run;
## kept here so the parameter table -> distribution wiring only needs to be written once.
draw_parameter <- function(params, name) {
  row <- params[name, ]
  dist <- tolower(trimws(row$distribution))
  if (dist == "log-normal") {
    exp(rnorm(1, mean = row$ln_mean, sd = row$ln_sd))
  } else if (dist %in% c("direchlet", "dirichlet")) {
    ## Excel draws independent Gamma(alpha, 1) per component and renormalises across the
    ## Dirichlet group afterwards (handled at the group level, not per-parameter).
    rgamma(1, shape = row$lower, rate = 1) # `lower` holds alpha (mean*1000) for these rows
  } else if (dist == "beta") {
    a <- row$lower
    b <- row$upper
    rbeta(1, a, b)
  } else {
    row$mean # fixed / structural parameter, no uncertainty specified
  }
}
