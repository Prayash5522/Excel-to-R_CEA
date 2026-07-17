## Age-varying rate / probability functions.
##
## These rebuild, from primitive Parameters-sheet values, the same quantities the workbook's
## "time varying mortality/incidence/prevalence" helper sheets compute -- confirmed cell by
## cell against the Markov model (NRT 20to29) sheet. Two conversions recur throughout the
## workbook and are kept explicit here so PSA can later resample the *primitive* rates and
## have everything downstream recompute automatically:
##   rate_to_prob(rate)          : annual hazard rate -> annual probability, 1 - exp(-rate)
##   rate_to_prob(rate * RR)     : same, with a relative-risk multiplier applied to the rate
##
## Source cross-checks (Markov sheet cell / Parameters row):
##   ET  = 1-EXP(-(mortality_rate lookup))                    -> rOAT_allcausemortality_<band8>
##   EU  = Parameters!p_IHD_mortality_<band8>                 -> 1-exp(-allcause*RR_IHD_mortality_<rr6>)
##   EW/EX/EY = Stroke AS/MS/SS mortality (same pattern, RR bands off Stroke(AS) breakpoints;
##              MS and SS carry an extra "compared to AS" RR multiplier)
##   EV  = COPD mortality -> RR_COPD_mortality has no age bands (applies at every age)
##   EZ  = Lung-cancer mortality -> RR_lungcancer_mortality_<rr2>
##   FS.. = primary disease incidence by smoking intensity -> 1-EXP(-(inc_<band5>_<disease>_<intensity>))
##   GH..HP = secondary/cross-disease transition rates -> inc_<band5>_<pathway>_<intensity>,
##            used directly as annual probabilities (workbook applies no exp() transform here)
##   FA  = age-adjusted OAT utility -> u_OAT_<utility band>

rate_to_prob <- function(rate) 1 - exp(-rate)

## intensity name differs by parameter family: primary incidence params abbreviate
## "moderate" to "mod"; every other family spells it out in full.
intensity_suffix <- function(intensity, family = c("incidence", "other")) {
  family <- match.arg(family)
  if (family == "incidence") {
    c(light = "light", moderate = "mod", heavy = "heavy")[intensity]
  } else {
    c(light = "light", moderate = "moderate", heavy = "heavy")[intensity]
  }
}

mortality_allcause_prob <- function(params, age) {
  band <- band_mortality8(age)
  rate_to_prob(get_param(params, paste0("rOAT_allcausemortality_", band)))
}

mortality_disease_prob <- function(params, age, disease = c("IHD", "COPD", "StrokeAS", "StrokeMS", "StrokeSS", "LC")) {
  disease <- match.arg(disease)
  band8 <- band_mortality8(age)
  allcause <- get_param(params, paste0("rOAT_allcausemortality_", band8))
  rr <- switch(disease,
    IHD = get_param(params, paste0("RR_IHD_mortality_", band_ihd_mortality_rr(age))),
    COPD = get_param(params, "RR_COPD_mortality"),
    StrokeAS = get_param(params, paste0("RR_Stroke__AS__mortality_", band_strokeAS_mortality_rr(age))),
    StrokeMS = get_param(params, paste0("RR_Stroke__AS__mortality_", band_strokeAS_mortality_rr(age))) *
      get_param(params, "RR_Stroke__MS__mortality__compared_to_AS"),
    StrokeSS = get_param(params, paste0("RR_Stroke__AS__mortality_", band_strokeAS_mortality_rr(age))) *
      get_param(params, "RR_Stroke__SS__mortality__compared_to_AS"),
    LC = get_param(params, paste0("RR_lungcancer_mortality_", band_lc_mortality_rr(age)))
  )
  rate_to_prob(allcause * rr)
}

## Primary incidence (healthy -> disease), by smoking intensity, for a single disease.
incidence_primary_prob <- function(params, age, disease = c("IHD", "stroke", "COPD", "lungcancer"),
                                    intensity = c("light", "moderate", "heavy")) {
  disease <- match.arg(disease)
  intensity <- match.arg(intensity)
  band <- band_disease5(age)
  suf <- intensity_suffix(intensity, "incidence")
  name <- paste0("inc_", band, "_", disease, "_", suf)
  rate_to_prob(get_param(params, name))
}

## Secondary / cross-disease transition probabilities (already probabilities in the workbook;
## no rate_to_prob() transform applied at the point of use).
transition_prob <- function(params, age, pathway, intensity = c("light", "moderate", "heavy")) {
  intensity <- match.arg(intensity)
  band <- band_disease5(age)
  name <- paste0("inc_", band, "_", pathway, "_", intensity)
  get_param(params, name)
}

utility_oat <- function(params, age) {
  get_param(params, paste0("u_OAT_", band_oat_utility(age)))
}

utility_oat_quit <- function(params, age) {
  utility_oat(params, age) + get_param(params, "u_quitting")
}
