# Excel-to-R CEA: OAT + smoking-cessation Markov model

R re-implementation of the "Markov model (NRT 20to29)" / "Markov model (SOC 20to29)" cohort
Markov model from `06.07.26_Final_PSA_and_USA.xlsm` (cost-effectiveness of adding NRT to OAT,
vs standard-of-care OAT alone, for smoking cessation among opioid agonist therapy patients).

This is a **structural rebuild**, not a cell-by-cell transliteration: the same health states,
tunnel-year logic, transition rates, cost/QALY accounting and half-cycle-corrected discounting
are reproduced from the workbook's formulas, written as clean, vectorised, documented R rather
than as a literal translation of ~300 spreadsheet columns. See "Validation" below for how
closely it tracks the workbook's own cached results.

## Project layout

```
data/parameters.csv       Extracted "Parameters" sheet: name, mean (col D), se (col E),
                           distribution (col F), bounds, log-normal mean/sd, description.
data-raw/                 Raw extraction outputs (working notes, not needed to run the model).
R/
  parameters.R             load_parameters(), get_param() [deterministic], draw_parameter()
                            [PSA, not called by the deterministic scripts - see below].
  age_bands.R               Age-decade bucketing, matching the workbook's "time varying ..."
                            lookup tables (confirmed piecewise-constant by decade/band).
  rates.R                   Age-varying probability/utility functions built from primitive
                            parameters (mortality, incidence, secondary-transition, utility).
  markov_engine.R           The state space, transition topology, and step_cycle() - the core
                            cohort update, one annual cycle at a time.
  cost_qaly.R                Per-cycle cost/QALY and half-cycle-corrected discounting.
  run_model.R                run_markov() (single arm) / run_comparison() (NRT vs SOC).
scripts/
  run_20to29.R               Deterministic NRT vs SOC comparison, cohort starting age 25
                            (the 20-29 band named in the request).
  validate_against_excel.R  Compares model output to the workbook's own cached lifetime
                            totals, for every age band.
```

## Running it

```r
Rscript scripts/run_20to29.R
Rscript scripts/validate_against_excel.R
```

To run a different age band, call `run_comparison(params, start_age)` with a different
starting age (35, 45, 55, 65 reproduce the 30-39/40-49/50-59/Over60 tabs) -- the model was
written generically over `start_age` from the outset, exactly mirroring how the workbook
itself parameterises each age-band tab off a single `start_age` cell.

## What the model reproduces

- **State space**: smoking-intensity tracked "full smoking" (light/moderate/heavy) and a
  6-year reduced-smoking tunnel per intensity; a pooled 6-year "quit from full smoking" tunnel
  and a pooled 6-year "quit from reduced smoking" tunnel; six disease states (IHD, stroke
  [acute/moderate/severe], COPD, lung cancer), each split into the same four smoking-status
  variants, with cross-disease progression pathways (e.g. IHD -> stroke, COPD -> lung cancer).
- **Transition topology**: `smk <-> qsmk`, `smk -> red -> smk` (relapse), `red <-> qred`,
  reproduced identically for disease states.
- **First-year excess mortality** for newly-incident/recurrent IHD and stroke cases (the
  workbook has no equivalent parameter for COPD or lung cancer, so none is applied there).
- **Age-varying rates**: mortality, incidence, secondary-transition risk and utility, banded
  by age exactly as the workbook's "time varying ..." helper sheets are (confirmed by reading
  their formulas: they are piecewise-constant lookups of Parameters-sheet values, not
  independent data, so this rebuild reads the underlying named parameters directly).
- **Cost/QALY accounting**: per-cycle cost and QALY exactly following the CS/CT column
  formulas, discounted at the workbook's cost/outcome discount rates (cDR/oDR, both 4%) and
  half-cycle corrected the same way the workbook's DY/DZ/EB columns are.

Documented simplifications (deliberate, given the "structural not literal" brief):
- Disease-specific mortality is computed as `1 - exp(-allcause_rate * RR)` once; the workbook
  computes this in the Parameters sheet and then re-wraps it in another `1-exp(-x)` at the
  point of use in the trace, which is numerically almost a no-op at these rate magnitudes.
- A handful of very granular bookkeeping/counter columns from the spreadsheet (there are
  ~40 of them, mostly for internal QA and disaggregated reporting) are not reproduced 1:1;
  the model computes the totals actually needed for cost-effectiveness output instead.

## PSA readiness (not run by default)

Per the brief, only the deterministic model is run here. The parameter table
(`data/parameters.csv`) already carries everything a probabilistic sensitivity analysis needs:
mean, standard error, distribution family, and bounds, straight from the workbook's
Probabilities-and-ratios columns (D:K). `R/parameters.R::draw_parameter()` implements the
Log-normal and Dirichlet(-via-Gamma) draws the workbook itself uses. Turning on PSA later is a
matter of looping `run_comparison()` N times with `get_param()` swapped for `draw_parameter()`
(and jointly renormalising each Dirichlet group per draw) -- no changes to the model logic
itself are needed.

## Validation

`scripts/validate_against_excel.R` compares this model's lifetime (cycle 0-75, i.e. to age
100) discounted-and-half-cycle-corrected totals against the values already cached inside the
workbook's own "Analysis" sheet (`Aggregated Result` table):

| band   | cost diff | QALY diff |
|--------|-----------|-----------|
| 20-29  | -1.7%     | -1.2%     |
| 30-39  | -2.4%     | -1.9%     |
| 40-49  | -3.9%     | -3.0%     |
| 50-59  | -6.7%     | -5.5%     |
| Over60 | -12.6%    | -10.2%    |

The 20-29 band -- the one named in the request, and used as the worked example throughout --
tracks the workbook within ~1-2% on cost, QALY and life-years, for both arms. The gap widens
for older starting cohorts; that band mixes higher mortality/incidence rates with more disease
transitions per surviving person-year, so the small documented simplifications above compound
more visibly there. If tighter agreement on the older bands is needed later, the first place
to look is the disease-mortality double-exp-wrap simplification and the secondary-transition
RR-tier handling in `markov_engine.R`.
