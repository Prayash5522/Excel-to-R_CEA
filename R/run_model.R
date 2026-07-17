## Top-level driver: runs one arm (NRT or SOC) for one starting age band from cycle 0 to
## `n_cycles`, cycle-by-cycle, and returns per-cycle totals ready for cost-effectiveness
## aggregation. 76 cycles (ages start_age .. start_age+75) matches the workbook, which runs
## every age-band tab out to age 100.

run_markov <- function(params, arm = c("NRT", "SOC"), start_age, n_cycles = 75) {
  arm <- match.arg(arm)

  s <- init_cohort(params, start_age)
  ages <- start_age + 0:n_cycles

  alive <- numeric(n_cycles + 1)
  cost <- numeric(n_cycles + 1)
  qaly <- numeric(n_cycles + 1)
  deaths <- numeric(n_cycles + 1)
  new_cases_mat <- matrix(0, nrow = n_cycles + 1, ncol = length(diseases), dimnames = list(NULL, diseases))

  alive[1] <- total_alive(s)
  cost[1] <- cycle_cost(s, stats::setNames(numeric(length(diseases)), diseases), params, arm)
  qaly[1] <- cycle_qaly(s, stats::setNames(numeric(length(diseases)), diseases), ages[1], params)

  for (t in seq_len(n_cycles)) {
    step <- step_cycle(s, ages[t + 1], params, arm)
    s <- step$state
    alive[t + 1] <- total_alive(s)
    deaths[t + 1] <- step$deaths
    new_cases_mat[t + 1, ] <- step$new_cases[diseases]
    cost[t + 1] <- cycle_cost(s, step$new_cases, params, arm)
    qaly[t + 1] <- cycle_qaly(s, step$new_cases, ages[t + 1], params)
  }

  cDR <- get_param(params, "cDR")
  oDR <- get_param(params, "oDR")
  cycles <- 0:n_cycles

  list(
    arm = arm,
    start_age = start_age,
    ages = ages,
    alive = alive,
    cost = cost,
    qaly = qaly,
    deaths = deaths,
    new_cases = new_cases_mat,
    final_state = s,
    total_cost = half_cycle_discounted_total(cost, cycles, cDR),
    total_qaly = half_cycle_discounted_total(qaly, cycles, oDR),
    total_ly = half_cycle_discounted_total(alive, cycles, oDR),
    total_deaths = sum(deaths)
  )
}

## Deterministic, two-arm comparison for a single starting age band (mirrors the "Markov
## model (NRT 20to29)" vs "Markov model (SOC 20to29)" tab pair; pass a different `start_age`
## to reproduce the 30-39 / 40-49 / 50-59 / Over60 tabs with the same code).
run_comparison <- function(params, start_age) {
  nrt <- run_markov(params, "NRT", start_age)
  soc <- run_markov(params, "SOC", start_age)

  inc_cost <- nrt$total_cost - soc$total_cost
  inc_qaly <- nrt$total_qaly - soc$total_qaly

  list(
    NRT = nrt, SOC = soc,
    incremental_cost = inc_cost,
    incremental_qaly = inc_qaly,
    icer = if (inc_qaly != 0) inc_cost / inc_qaly else NA_real_
  )
}
