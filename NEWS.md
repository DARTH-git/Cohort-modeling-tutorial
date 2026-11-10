# Cohort-modeling-tutorial 0.2.0 (2026-11-10)

This repository is no longer maintained. It has been superseded by
[cohort-modeling-tutorial-intro](https://github.com/DARTH-git/cohort-modeling-tutorial-intro)
and
[cohort-modeling-tutorial-timedep](https://github.com/DARTH-git/cohort-modeling-tutorial-timedep),
which accompany the published articles. This release only fixes known errors.

## Bug fixes

* The analysis scripts use `darthtools::gen_wcc()` for Simpson's 1/3rd
  within-cycle correction. `darthtools` up to v0.4.0 assigned the 4/3 and 2/3
  weights by position in the 1-based R vector instead of by time, so the
  weights were `1/3, 2/3, 4/3, ..., 1/3` instead of `1/3, 4/3, 2/3, ..., 1/3`.
  This is fixed in darthtools v0.4.1. The scripts now stop with an
  explanatory message if an older `darthtools` is installed. The manuscript
  files in `Manuscripts/` are kept unchanged and contain the same error in
  their `within-cycle-vector` code chunks. Reported by Lauren N. Steimle and
  @ambtnhlr
  ([DARTH-git/cohort-modeling-tutorial-intro#2](https://github.com/DARTH-git/cohort-modeling-tutorial-intro/issues/2)).
* `calculate_ce_out()` of the age- and history-dependent model
  (`R/Functions STM_03.R`) used the old names of the Weibull parameters
  (`r_S1S2_lambda`, `r_S1S2_gamma`), so `analysis/STM_03.R` stopped with
  "object 'r_S1S2_lambda' not found". It now uses `r_S1S2_scale` and
  `r_S1S2_shape`.
* Fixed comments that described `u_D` as the "utility when Healthy" (closes #2).

## Documentation

* README: notice that the repository has been superseded, the correction
  above, and fixed links to the analysis scripts.

# Cohort-modeling-tutorial 0.1 (2020-11-13)

* Code used for the first revision of the cSTM tutorial
  (arXiv:2001.07824v1).
