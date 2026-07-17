## Age-banding helpers.
##
## The source workbook looks up several families of age/decade-specific parameters via
## VLOOKUP against "time varying ..." tables. Those tables were confirmed (by inspecting
## their formulas) to be piecewise-constant within each band -- i.e. they are just a
## decade-bucketed re-expression of Parameters-sheet values, not independent data. Rather
## than importing the lookup tables verbatim, we replicate the bucketing rule and read the
## underlying named parameter directly, which keeps every number traceable to one row in
## data/parameters.csv.
##
## Four distinct banding schemes are used in the workbook (verified against the relevant
## Parameters rows / defined names):
##   - all-cause mortality & disease-specific mortality probabilities: 8 decade bands
##   - disease incidence, prevalence, and secondary-event transition rates: 5 bands
##     (20-29, 30-39, 40-49, 50-59, Over60 -- the model does not refine risk beyond 60)
##   - IHD mortality relative-risk: 6 bands (20-54, 55-64, 65-74, 75-84, 85-89, 90+)
##   - Stroke(AS) mortality relative-risk: 3 bands (20-39, 40-79, 80+)
##   - Lung-cancer mortality relative-risk: 2 bands (20-70, 70+)
##   - OAT baseline utility: 4 bands (<26, 26-40, 41-60, 60+)

band_mortality8 <- function(age) {
  cut(age,
      breaks = c(-Inf, 29, 39, 49, 59, 69, 79, 89, Inf),
      labels = c("20_29", "30_39", "40_49", "50_59", "60_69", "70_79", "80_89", "90"),
      right = TRUE) |> as.character()
}

band_disease5 <- function(age) {
  cut(age,
      breaks = c(-Inf, 29, 39, 49, 59, Inf),
      labels = c("20_29", "30_39", "40_49", "50_59", "Over60"),
      right = TRUE) |> as.character()
}

band_ihd_mortality_rr <- function(age) {
  cut(age,
      breaks = c(-Inf, 54, 64, 74, 84, 89, Inf),
      labels = c("20_54", "55_64", "65_74", "75_84", "85_89", "90above"),
      right = TRUE) |> as.character()
}

band_strokeAS_mortality_rr <- function(age) {
  cut(age,
      breaks = c(-Inf, 39, 79, Inf),
      labels = c("20_39", "40_79", "80above"),
      right = TRUE) |> as.character()
}

band_lc_mortality_rr <- function(age) {
  ifelse(age <= 70, "20_70", "70above")
}

band_oat_utility <- function(age) {
  cut(age,
      breaks = c(-Inf, 25, 40, 60, Inf),
      labels = c("lessthan25", "26_40", "41_60", "above_60"),
      right = TRUE) |> as.character()
}
