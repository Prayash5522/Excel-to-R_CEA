## Per-cycle cost and QALY accounting, plus half-cycle-corrected discounting.
##
## Rebuilt from the CS (cost), CT/CU (QALY), CV (life years), DI:DL (discounted), and DY:EB
## (half-cycle-corrected discounted) columns of the Markov sheet (formulas read directly from
## the workbook, see cost_qaly_cells.txt in the extraction notes). Confirmed pattern:
##   discounted_X[t]      = X[t] / (1 + discount_rate)^cycle[t]
##   half_cycle_X[t]       = (discounted_X[t] + discounted_X[t-1]) / 2         for t = 1..n_cycles
##   total_X                = sum_{t=1}^{n_cycles} half_cycle_X[t]             (cycle 0 excluded,
##                              exactly as the workbook's DY column is blank at cycle 0)
## Costs use the cost discount rate (cDR); QALYs and life-years use the outcome discount rate
## (oDR) -- both read from the Parameters sheet (both are 4% in the source workbook).

total_alive <- function(s) {
  sum(s$smk) + sum(vapply(s$red, sum, numeric(1))) + sum(s$qsmk) + sum(s$qred) +
    sum(vapply(s$dis, sum, numeric(1)))
}

## Per-cycle cost, following the CS formula: (population in each pool) * (drug/OAT cost for
## that pool), plus a flat per-event cost on newly-incident IHD and stroke events.
cycle_cost <- function(s, new_cases, params, arm) {
  suffix <- if (arm == "NRT") "recommended_dose" else NULL
  nrt_cost <- function(pool) if (arm == "NRT") get_param(params, paste0("c_NRT_", pool, "_", suffix)) else 0
  c_oat <- get_param(params, "c_OAT")

  c_light <- nrt_cost("light_smoker") + c_oat
  c_moderate <- nrt_cost("moderate_smoker") + c_oat
  c_heavy <- nrt_cost("heavy_smoker") + c_oat
  c_quit <- nrt_cost("quiter") + c_oat
  c_disease_group <- nrt_cost("disease_group") + c_oat

  c_disease <- c(
    IHD = get_param(params, "c_IHD"), StrokeAS = get_param(params, "c_stroke__AS"),
    StrokeMS = get_param(params, "c_stroke__MS"), StrokeSS = get_param(params, "c_stroke__SS"),
    COPD = get_param(params, "c_COPD"), LC = get_param(params, "c_lung_cancer")
  )

  healthy_cost <- s$smk[["light"]] * c_light + s$smk[["moderate"]] * c_moderate + s$smk[["heavy"]] * c_heavy +
    sum(s$red$light) * c_light + sum(s$red$moderate) * c_moderate + sum(s$red$heavy) * c_heavy +
    (sum(s$qsmk) + sum(s$qred)) * c_quit

  disease_cost <- sum(vapply(diseases, function(d) {
    (s$dis[[d]][["smk"]] + s$dis[[d]][["red"]]) * (c_disease_group + c_disease[[d]]) +
      (s$dis[[d]][["qsmk"]] + s$dis[[d]][["qred"]]) * (c_quit + c_disease[[d]])
  }, numeric(1)))

  event_cost <- new_cases[["IHD"]] * get_param(params, "c_eventIHD") +
    (new_cases[["StrokeAS"]] + new_cases[["StrokeMS"]] + new_cases[["StrokeSS"]]) * get_param(params, "c_eventstroke")

  healthy_cost + disease_cost + event_cost
}

## Per-cycle QALYs, following the CT formula: baseline age-adjusted utility for smk/red
## (FA) vs qsmk/qred (FB = FA + u_quitting), capped downward by disease-specific utility for
## anyone with a diagnosed disease. Lung cancer utility depends on whether the case is newly
## incident this cycle (u_lung_cancer_1styear) or established (u_lung_cancer_2ndyrnbeyond) --
## the workbook applies the age-adjusted utility (FA) uniformly for lung cancer regardless of
## smoking/quit status, which this replicates (see cost_qaly_cells.txt: the CT formula's lung
## cancer terms use `FA`, not the smk/quit-specific FA/FB split used for IHD/stroke/COPD).
cycle_qaly <- function(s, new_cases, age, params) {
  fa <- utility_oat(params, age)
  fb <- utility_oat_quit(params, age)

  u_disease <- c(
    IHD = get_param(params, "u_IHD"), StrokeAS = get_param(params, "u_stroke_AS"),
    StrokeMS = get_param(params, "u_stroke_MS"), StrokeSS = get_param(params, "u_stroke_SS"),
    COPD = get_param(params, "u_COPD")
  )
  u_lc_new <- get_param(params, "u_lung_cancer_1styear")
  u_lc_established <- get_param(params, "u_lung_cancer_2ndyrnbeyond")

  healthy_qaly <- (sum(s$smk) + sum(vapply(s$red, sum, numeric(1)))) * fa +
    (sum(s$qsmk) + sum(s$qred)) * fb

  disease_qaly <- sum(vapply(c("IHD", "StrokeAS", "StrokeMS", "StrokeSS", "COPD"), function(d) {
    u <- u_disease[[d]]
    (s$dis[[d]][["smk"]] + s$dis[[d]][["red"]]) * min(fa, u) +
      (s$dis[[d]][["qsmk"]] + s$dis[[d]][["qred"]]) * min(fb, u)
  }, numeric(1)))

  lc_total <- sum(s$dis$LC)
  lc_new <- min(new_cases[["LC"]], lc_total)
  lc_established <- lc_total - lc_new
  lc_qaly <- lc_new * min(fa, u_lc_new) + lc_established * min(fa, u_lc_established)

  healthy_qaly + disease_qaly + lc_qaly
}

## ---- discounting & aggregation --------------------------------------------------------------

## Half-cycle-corrected, discounted totals across the whole run. `raw` is a numeric vector of
## per-cycle values (cost or QALY or life-years), one per cycle 0..n_cycles, `cycles` the
## matching cycle-number vector, `rate` the applicable discount rate.
half_cycle_discounted_total <- function(raw, cycles, rate) {
  discounted <- raw / (1 + rate)^cycles
  n <- length(discounted)
  half_cycle <- (discounted[-1] + discounted[-n]) / 2 # pairs (t-1, t) for t = 1..n_cycles
  sum(half_cycle)
}
