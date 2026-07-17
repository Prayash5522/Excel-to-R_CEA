## Core cohort Markov engine for the OAT + smoking-cessation model.
##
## Rebuilds the state-transition logic of "Markov model (NRT 20to29)" / "Markov model
## (SOC 20to29)" as a generic, age-band-parameterised R model (works for any start_age).
## Verified against the workbook by reading the row-10 (steady-cycle) formulas for every
## column and confirming the formula pattern is identical for every cycle >= 1 (row 10 vs
## row 79 formulas are structurally identical once row-relative cell references are
## normalised - only the literal cycle/age values differ). Cycle 0 (row 4) is the one
## structurally different row and is handled separately by init_cohort().
##
## STATE SPACE (mirrors Excel columns D:BH)
##   Smoking/quit "healthy" states (no diagnosed disease):
##     smk[light|moderate|heavy]                        -- Excel D,E,F
##     red[light|moderate|heavy][y1..y5,y6p]             -- Excel G:L, M:R, S:X (reduced-smoking tunnel)
##     qsmk[y1..y5,y6p]                                   -- Excel Y:AD  (quit-from-full-smoking tunnel)
##     qred[y1..y5,y6p]                                   -- Excel AE:AJ (quit-from-reduced tunnel)
##   Disease states, 4 smoking-status variants each, for IHD / StrokeAS / StrokeMS / StrokeSS / COPD / LC:
##     dis[[disease]][smk|red|qsmk|qred]                  -- Excel AK:BH
##
## TOPOLOGY (confirmed from the G10/Y10/AE10/D10/AK10/AL10 formulas)
##   smk <-> qsmk   (p_quit_rate_<arm> forward; p_relapse_rate_quit_<arm>_<yearband> backward)
##   smk  -> red    (p_primary_outcome_<arm>, one-way entry into the reduced track)
##   red <-> smk    (p_relapse_rate_primary_outcome_<arm>_<yearband>, reduced relapsing to full smoking)
##   red <-> qred   (p_reduced_to_quit_rate_<arm> forward; p_relapse_rate_quit_to_reduced_<arm>_<yearband> backward)
## The red tunnel is the only one with *two* simultaneous exits each cycle (to smk via relapse,
## to qred via reduced_to_quit_rate) -- both are subtracted from the same survivor base, not
## applied sequentially (confirmed against the H10 formula: "(surv) - (surv)*reduced_to_quit -
## (surv)*relapse", i.e. both rates net off the same base).
##
## Disease states reuse this exact topology and exit-rate structure (confirmed against the
## fully-decoded IHD block: AK/AL/AM/AN formulas), with disease-specific mortality substituted
## for all-cause mortality, and cross-disease "pathway" transitions substituted for primary
## disease incidence.
##
## FIRST-YEAR EXCESS MORTALITY: newly-incident or newly-recurrent IHD and stroke cases use
## p_mortality_IHD_first_year / p_mortality_stroke_first_year in the cycle they occur, instead
## of the disease's chronic mortality probability; COPD and lung cancer have no such workbook
## parameter (none exists in the Parameters sheet) and use chronic mortality from cycle one.
##
## INTENSITY TRACKING: states that track smoking intensity (smk, red[i]) use intensity-specific
## rates directly. States that don't (qsmk, qred, every disease-state bucket) fall back to the
## population-average light/moderate/heavy mix (prv_light_smoker etc.) -- confirmed against the
## AK10 formula, which weights secondary-transition risk by prv_light_smoker/... even though AK
## is the "still fully smoking" IHD bucket.

diseases <- c("IHD", "StrokeAS", "StrokeMS", "StrokeSS", "COPD", "LC")
intensities <- c("light", "moderate", "heavy")
tunnel_years <- c("y1", "y2", "y3", "y4", "y5", "y6p")
statuses <- c("smk", "red", "qsmk", "qred")

## disease reached by primary incidence from a healthy state, and the `disease` argument
## incidence_primary_prob() expects for it (StrokeMS/StrokeSS are only reachable via
## progression, never directly diagnosed)
primary_incidence_disease <- c(IHD = "IHD", StrokeAS = "stroke", COPD = "COPD", LC = "lungcancer")
## reverse map: RR family name (used in RR_OAT_to_<x>_...) for each disease state
rr_family <- c(IHD = "IHD", StrokeAS = "stroke", StrokeMS = "stroke", StrokeSS = "stroke", COPD = "COPD", LC = "lungcancer")

## outgoing secondary/cross-disease pathways per disease state: target disease, the workbook's
## pathway name (feeds data/parameters.csv names inc_<ageband>_<pathway>_<intensity>), and which
## RR_OAT_to_<rr>_... family protects against arriving at the target via this pathway.
disease_pathways <- list(
  IHD = list(
    list(target = "StrokeAS", path = "IHD_to_stroke", rr = "stroke"),
    list(target = "COPD", path = "IHD_to_COPD", rr = "COPD"),
    list(target = "LC", path = "IHD_to_lungcancer", rr = "lungcancer"),
    list(target = "IHD", path = "IHD_to_secondaryIHD", rr = "IHD")
  ),
  StrokeAS = list(
    list(target = "StrokeMS", path = "Stroke_AS_to_secondary_stroke", rr = "stroke"),
    list(target = "LC", path = "stroke_to_lungcancer", rr = "lungcancer")
  ),
  StrokeMS = list(
    list(target = "StrokeSS", path = "Stroke_MS_to_stroke_SS", rr = "stroke"),
    list(target = "LC", path = "stroke_to_lungcancer", rr = "lungcancer")
  ),
  StrokeSS = list(), # terminal: no further progression modelled
  COPD = list(
    list(target = "StrokeAS", path = "COPD_to_stroke", rr = "stroke"),
    list(target = "LC", path = "COPD_to_lungcancer", rr = "lungcancer")
  ),
  LC = list() # terminal
)

first_year_mortality_param <- c(
  IHD = "p_mortality_IHD_first_year", StrokeAS = "p_mortality_stroke_first_year",
  StrokeMS = "p_mortality_stroke_first_year", StrokeSS = "p_mortality_stroke_first_year",
  COPD = NA_character_, LC = NA_character_
)

quit_mortality_rr_param <- c(
  IHD = "RR_mortality_quit_CVD", StrokeAS = "RR_mortality_quit_CVD",
  StrokeMS = "RR_mortality_quit_CVD", StrokeSS = "RR_mortality_quit_CVD",
  COPD = "RR_mortality_quit_COPD", LC = "RR_mortality_quit_LC"
)

## ---- arm-specific parameter name resolution -----------------------------------------------
## Both arms share the model above; only the parameter *names* differ (NRT vs soc), plus the
## NRT arm carries an NRT drug-cost term the SOC arm does not (see cost_qaly.R).
arm_param <- function(arm, key) {
  suffix <- if (arm == "NRT") "NRT" else "soc"
  switch(key,
    primary_outcome        = paste0("p_primary_outcome_", suffix),
    quit_rate               = paste0("p_quit_rate_", suffix),
    reduced_to_quit_rate    = paste0("p_reduced_to_quit_rate_", suffix),
    relapse_primary_y1      = paste0("p_relapse_rate_primary_outcome_", suffix, "_year_1"),
    relapse_primary_y2_5    = paste0("p_relapse_rate_primary_outcome_", suffix, "_year_2n5"),
    relapse_primary_y6p     = paste0("p_relapse_rate_primary_outcome_", suffix, "_year_6nbeyond"),
    relapse_quit_y1         = if (arm == "NRT") "p_relapse_rate_quit_NRT__year_1" else "p_relapse_rate_quit_soc__year_1",
    relapse_quit_y2_5       = paste0("p_relapse_rate_quit_", suffix, "_year2_5"),
    relapse_quit_y6p        = paste0("p_relapse_rate_quit_", suffix, "_year5nabove"),
    relapse_primary_disease = paste0("p_relapse_rate_primary_outcome_", suffix), # bare, no year suffix - used by disease states (no tunnel-year tracking)
    relapse_quit_disease    = if (arm == "NRT") "p_relapse_rate_quit_NRT__disease_state" else "p_relapse_rate_quit_soc__disease_state",
    relapse_qred_y1         = paste0("p_relapse_rate_quit_to_reduced_", suffix, "_year1"),
    relapse_qred_y2_5       = paste0("p_relapse_rate_quit_to_reduced_", suffix, "_year2_5"),
    relapse_qred_y6p        = paste0("p_relapse_rate_quit_to_reduced_", suffix, "_year5nabove"),
    relapse_qred_disease    = paste0("p_relapse_rate_quit_to_reduced_", suffix, "_disease_state"),
    stop(sprintf("unknown arm_param key: %s", key))
  )
}

## ---- small numeric helpers ------------------------------------------------------------------

## population-mix-weighted rate, used wherever the workbook can't track smoking intensity
## within a pool (quit tunnels, all disease states) and falls back to the population-average
## light/moderate/heavy split (prv_light_smoker etc.)
mix_weighted <- function(params, rate_by_intensity) {
  w <- c(
    light = get_param(params, "prv_light_smoker"),
    moderate = get_param(params, "prv_moderate_smoker"),
    heavy = get_param(params, "prv_heavy_smoker")
  )
  sum(w * rate_by_intensity[names(w)])
}

## 4-disease primary incidence probability vector for one smoking intensity, at this age
primary_risk_vec <- function(params, age, intensity) {
  stats::setNames(
    vapply(names(primary_incidence_disease), function(dname) {
      incidence_primary_prob(params, age, primary_incidence_disease[[dname]], intensity)
    }, numeric(1)),
    names(primary_incidence_disease)
  )
}

## reduction factor (RR) applied to primary incidence risk for a reduced/quit smoker
oat_rr <- function(params, disease_key, status, intensity, first_five) {
  nm <- paste0(
    "RR_OAT_to_", disease_key, "_", status, "_", intensity,
    if (first_five) "_first_five" else ""
  )
  get_param(params, nm)
}

## intensity-specific 4-disease risk vector, RR-adjusted for reduced/quit status
rr_adjusted_risk <- function(params, age, intensity, status, first_five) {
  raw <- primary_risk_vec(params, age, intensity)
  vapply(names(raw), function(dname) {
    raw[[dname]] * oat_rr(params, rr_family[[dname]], status, intensity, first_five)
  }, numeric(1))
}

## population-mix-weighted, RR-adjusted 4-disease risk vector (for qsmk/qred, which don't
## track intensity)
rr_adjusted_risk_mixed <- function(params, age, status, first_five) {
  m <- vapply(intensities, function(i) rr_adjusted_risk(params, age, i, status, first_five), numeric(length(primary_incidence_disease)))
  # m: rows = disease, cols = intensity
  apply(m, 1, function(row_by_intensity) mix_weighted(params, stats::setNames(row_by_intensity, intensities)))
}

survive_healthy <- function(n, mort_prob, risk_vec) n * (1 - mort_prob - sum(risk_vec))

## Tunnel bucket step for a y1..y6p vector. `mort5`/`mortL` are the mortality probability
## applied to the first-five-year buckets vs the year-6+ bucket (these differ only for the
## quit tunnels, which pick up the RR_mortality_quit long-run benefit once in the y6p bucket);
## `risk5`/`riskL` are 4-disease risk vectors (first-five vs year-6+ tier). `rate1/rate25/rate6p`
## is the primary exit rate (e.g. relapse back to smk/red) by bucket tier; `secondary_rate`
## (constant across buckets, may be 0) is an optional simultaneous second exit route (only the
## reduced-smoking tunnel uses this, exiting to the quit tunnel). Returns bucket survivors
## (pre-exit, needed by callers that route a *different* secondary exit off the same base),
## remaining y2..y6p carryover, and the two total outflows.
tunnel_step <- function(n_prev, mort5, mortL, risk5, riskL, rate1, rate25, rate6p, secondary_rate = 0) {
  surv <- c(
    y1 = survive_healthy(n_prev[["y1"]], mort5, risk5),
    y2 = survive_healthy(n_prev[["y2"]], mort5, risk5),
    y3 = survive_healthy(n_prev[["y3"]], mort5, risk5),
    y4 = survive_healthy(n_prev[["y4"]], mort5, risk5),
    y5 = survive_healthy(n_prev[["y5"]], mort5, risk5),
    y6p = survive_healthy(n_prev[["y6p"]], mortL, riskL)
  )
  exit_rate <- c(y1 = rate1, y2 = rate25, y3 = rate25, y4 = rate25, y5 = rate25, y6p = rate6p) + secondary_rate
  primary_share <- c(y1 = rate1, y2 = rate25, y3 = rate25, y4 = rate25, y5 = rate25, y6p = rate6p) / exit_rate
  primary_share[exit_rate == 0] <- 0
  total_exit <- surv * exit_rate
  out_primary <- sum(total_exit * primary_share)
  out_secondary <- sum(total_exit * (1 - primary_share))
  remain <- surv * (1 - exit_rate)
  list(
    deaths = sum(n_prev[c("y1", "y2", "y3", "y4", "y5")]) * mort5 + n_prev[["y6p"]] * mortL,
    surv_buckets = surv,
    out_primary = out_primary,
    out_secondary = out_secondary,
    carry_y2_y6p = c(
      y2 = remain[["y1"]], y3 = remain[["y2"]], y4 = remain[["y3"]],
      y5 = remain[["y4"]], y6p = remain[["y5"]] + remain[["y6p"]]
    )
  )
}

## sum of a 4-disease risk vector over the y1-y5 (first-five) and y6p (long-term) buckets,
## producing per-disease new-case totals for a tunnel
tunnel_new_cases <- function(n_prev, risk5, riskL) {
  early <- sum(n_prev[c("y1", "y2", "y3", "y4", "y5")])
  risk5 * early + riskL * n_prev[["y6p"]]
}

## ---- cohort initialisation (Excel row 4) -----------------------------------------------------

new_state <- function() {
  z6 <- function() stats::setNames(numeric(6), tunnel_years)
  list(
    smk = stats::setNames(numeric(3), intensities),
    red = stats::setNames(lapply(intensities, function(i) z6()), intensities),
    qsmk = z6(),
    qred = z6(),
    dis = stats::setNames(
      lapply(diseases, function(d) stats::setNames(numeric(4), statuses)),
      diseases
    )
  )
}

init_cohort <- function(params, start_age) {
  band <- band_disease5(start_age)
  cohort_suffix <- switch(band,
    "20_29" = "20_30", "30_39" = "31_40", "40_49" = "41_50",
    "50_59" = "51_60", "Over60" = "over60"
  )
  total <- get_param(params, "Cohort") * get_param(params, paste0("Cohort_", cohort_suffix))

  w <- c(
    light = get_param(params, "prv_light_smoker"),
    moderate = get_param(params, "prv_moderate_smoker"),
    heavy = get_param(params, "prv_heavy_smoker")
  )

  prv_disease <- function(disease_incprefix) {
    vapply(intensities, function(i) {
      get_param(params, paste0("prv_", band, "_", disease_incprefix, "_", intensity_suffix(i, "incidence")))
    }, numeric(1))
  }
  ihd0 <- total * sum(prv_disease("IHD") * w)
  strokeAS0 <- total * sum(prv_disease("stroke") * w)
  copd0 <- total * sum(prv_disease("COPD") * w)
  lc0 <- total * sum(prv_disease("lungcancer") * w)
  diseased_total <- ihd0 + strokeAS0 + copd0 + lc0

  s <- new_state()
  s$smk["light"] <- total * w[["light"]] - diseased_total * w[["light"]]
  s$smk["moderate"] <- (total - diseased_total) * w[["moderate"]]
  s$smk["heavy"] <- (total - diseased_total) * w[["heavy"]]
  s$dis$IHD["smk"] <- ihd0
  s$dis$StrokeAS["smk"] <- strokeAS0
  s$dis$COPD["smk"] <- copd0
  s$dis$LC["smk"] <- lc0
  s
}

## ---- one annual cycle --------------------------------------------------------------------------

step_cycle <- function(s, age, params, arm) {
  ET <- mortality_allcause_prob(params, age)
  RRmq <- get_param(params, "RR_mortality_quit")

  primary_outcome <- get_param(params, arm_param(arm, "primary_outcome"))
  quit_rate <- get_param(params, arm_param(arm, "quit_rate"))
  reduced_to_quit_rate <- get_param(params, arm_param(arm, "reduced_to_quit_rate"))
  rel_primary_y1 <- get_param(params, arm_param(arm, "relapse_primary_y1"))
  rel_primary_y2_5 <- get_param(params, arm_param(arm, "relapse_primary_y2_5"))
  rel_primary_y6p <- get_param(params, arm_param(arm, "relapse_primary_y6p"))
  rel_quit_y1 <- get_param(params, arm_param(arm, "relapse_quit_y1"))
  rel_quit_y2_5 <- get_param(params, arm_param(arm, "relapse_quit_y2_5"))
  rel_quit_y6p <- get_param(params, arm_param(arm, "relapse_quit_y6p"))
  rel_primary_disease <- get_param(params, arm_param(arm, "relapse_primary_disease"))
  rel_quit_disease <- get_param(params, arm_param(arm, "relapse_quit_disease"))
  rel_qred_y1 <- get_param(params, arm_param(arm, "relapse_qred_y1"))
  rel_qred_y2_5 <- get_param(params, arm_param(arm, "relapse_qred_y2_5"))
  rel_qred_y6p <- get_param(params, arm_param(arm, "relapse_qred_y6p"))
  rel_qred_disease <- get_param(params, arm_param(arm, "relapse_qred_disease"))

  p_relapse_dest <- c(
    light = get_param(params, "p_light_relapse"),
    moderate = get_param(params, "p_moderate_relapse"),
    heavy = get_param(params, "p_heavy_relapse")
  )

  new_s <- new_state()
  new_cases <- stats::setNames(numeric(length(diseases)), diseases)
  deaths <- 0
  first_year_excess_deaths <- 0

  route_new_case <- function(disease, status, amount) {
    fy_param <- first_year_mortality_param[[disease]]
    if (!is.na(fy_param)) {
      p_fy <- get_param(params, fy_param)
      deaths <<- deaths + amount * p_fy
      first_year_excess_deaths <<- first_year_excess_deaths + amount * p_fy
      amount <- amount * (1 - p_fy)
    }
    new_cases[disease] <<- new_cases[disease] + amount
    new_s$dis[[disease]][status] <<- new_s$dis[[disease]][status] + amount
  }

  ## ---------------- smk / red (intensity-tracked healthy states) ------------------------
  smk_survivors <- stats::setNames(numeric(3), intensities)
  red_out_primary <- stats::setNames(numeric(3), intensities) # relapse -> smk
  red_out_secondary <- stats::setNames(numeric(3), intensities) # reduced_to_quit -> qred

  for (i in intensities) {
    risk_i <- primary_risk_vec(params, age, i)
    surv <- survive_healthy(s$smk[[i]], ET, risk_i)
    smk_survivors[i] <- surv
    deaths <- deaths + s$smk[[i]] * ET
    for (dname in names(primary_incidence_disease)) route_new_case(dname, "smk", s$smk[[i]] * risk_i[[dname]])

    risk5 <- rr_adjusted_risk(params, age, i, "reduced", TRUE)
    riskL <- rr_adjusted_risk(params, age, i, "reduced", FALSE)
    ## the reduced-smoking tunnel gets no all-cause mortality benefit at any tunnel year
    ## (that benefit is modelled only for full cessation) - mort5 == mortL == ET throughout.
    tun <- tunnel_step(s$red[[i]], ET, ET, risk5, riskL, rel_primary_y1, rel_primary_y2_5, rel_primary_y6p,
                        secondary_rate = reduced_to_quit_rate)
    deaths <- deaths + tun$deaths
    nc <- tunnel_new_cases(s$red[[i]], risk5, riskL)
    for (dname in names(nc)) route_new_case(dname, "red", nc[[dname]])

    red_out_primary[i] <- tun$out_primary
    red_out_secondary[i] <- tun$out_secondary
    new_s$red[[i]][c("y2", "y3", "y4", "y5", "y6p")] <- tun$carry_y2_y6p
  }

  ## ---------------- qsmk / qred (pooled, mix-weighted healthy states) --------------------
  ## y6p (year 6+) picks up the long-run all-cause mortality benefit of quitting (RR_mortality_quit)
  qsmk_risk5 <- rr_adjusted_risk_mixed(params, age, "quit", TRUE)
  qsmk_riskL <- rr_adjusted_risk_mixed(params, age, "quit", FALSE)
  tun_qsmk <- tunnel_step(s$qsmk, ET, ET * RRmq, qsmk_risk5, qsmk_riskL, rel_quit_y1, rel_quit_y2_5, rel_quit_y6p)
  deaths <- deaths + tun_qsmk$deaths
  nc_qsmk <- tunnel_new_cases(s$qsmk, qsmk_risk5, qsmk_riskL)
  for (dname in names(nc_qsmk)) route_new_case(dname, "qsmk", nc_qsmk[[dname]])
  new_s$qsmk[c("y2", "y3", "y4", "y5", "y6p")] <- tun_qsmk$carry_y2_y6p

  qred_risk5 <- rr_adjusted_risk_mixed(params, age, "quit", TRUE)
  qred_riskL <- rr_adjusted_risk_mixed(params, age, "quit", FALSE)
  tun_qred <- tunnel_step(s$qred, ET, ET * RRmq, qred_risk5, qred_riskL, rel_qred_y1, rel_qred_y2_5, rel_qred_y6p)
  deaths <- deaths + tun_qred$deaths
  nc_qred <- tunnel_new_cases(s$qred, qred_risk5, qred_riskL)
  for (dname in names(nc_qred)) route_new_case(dname, "qred", nc_qred[[dname]])
  new_s$qred[c("y2", "y3", "y4", "y5", "y6p")] <- tun_qred$carry_y2_y6p

  ## ---------------- redistribute healthy populations for next cycle ----------------------
  for (i in intensities) {
    smk_out_primary <- smk_survivors[[i]] * primary_outcome
    smk_out_quit <- smk_survivors[[i]] * quit_rate
    new_s$smk[i] <- smk_survivors[[i]] - smk_out_primary - smk_out_quit +
      red_out_primary[[i]] + tun_qsmk$out_primary * p_relapse_dest[[i]]
    new_s$red[[i]]["y1"] <- smk_out_primary + tun_qred$out_primary * p_relapse_dest[[i]]
  }
  new_s$qsmk["y1"] <- sum(smk_survivors) * quit_rate
  new_s$qred["y1"] <- sum(red_out_secondary)

  ## ---------------- disease states --------------------------------------------------------
  for (d in diseases) {
    qmr <- get_param(params, quit_mortality_rr_param[[d]])
    chronic_mort <- mortality_disease_prob(params, age, d)
    out_pw <- disease_pathways[[d]]

    pathway_rate <- function(status) {
      vapply(out_pw, function(pw) {
        base <- mix_weighted(params, vapply(intensities, function(i) transition_prob(params, age, pw$path, i), numeric(1)))
        rr <- if (status == "smk") {
          1
        } else {
          mix_weighted(params, vapply(intensities, function(i) oat_rr(params, pw$rr, if (status == "red") "reduced" else "quit", i, FALSE), numeric(1)))
        }
        base * rr
      }, numeric(1))
    }

    for (status in statuses) {
      n <- s$dis[[d]][[status]]
      mort_prob <- if (status %in% c("qsmk", "qred")) chronic_mort * qmr else chronic_mort
      rates <- if (length(out_pw) > 0) pathway_rate(status) else numeric(0)
      surv <- n * (1 - mort_prob - sum(rates))
      deaths <- deaths + n * mort_prob

      out_primary <- switch(status,
        smk = surv * primary_outcome,
        red = surv * rel_primary_disease, # disease-state "red" relapses to "smk" (no year tunnel)
        0
      )
      out_secondary <- switch(status,
        smk = surv * quit_rate,
        red = surv * reduced_to_quit_rate,
        qsmk = surv * rel_quit_disease, # relapses back to smk
        qred = surv * rel_qred_disease, # relapses back to red
        0
      )
      remain <- surv - out_primary - out_secondary

      new_s$dis[[d]][[status]] <- new_s$dis[[d]][[status]] + remain
      if (status == "smk") {
        new_s$dis[[d]]["red"] <- new_s$dis[[d]]["red"] + out_primary
        new_s$dis[[d]]["qsmk"] <- new_s$dis[[d]]["qsmk"] + out_secondary
      } else if (status == "red") {
        new_s$dis[[d]]["smk"] <- new_s$dis[[d]]["smk"] + out_primary
        new_s$dis[[d]]["qred"] <- new_s$dis[[d]]["qred"] + out_secondary
      } else if (status == "qsmk") {
        new_s$dis[[d]]["smk"] <- new_s$dis[[d]]["smk"] + out_secondary
      } else if (status == "qred") {
        new_s$dis[[d]]["red"] <- new_s$dis[[d]]["red"] + out_secondary
      }

      if (length(out_pw) > 0) {
        for (k in seq_along(out_pw)) {
          route_new_case(out_pw[[k]]$target, status, n * rates[[k]])
        }
      }
    }
  }

  list(state = new_s, deaths = deaths, new_cases = new_cases, first_year_excess_deaths = first_year_excess_deaths)
}
