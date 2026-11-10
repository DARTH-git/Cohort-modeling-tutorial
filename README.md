# Cohort state-transition models in R: A Tutorial

> **This repository is no longer maintained.** It holds the code for the 2020 preprint (arXiv:2001.07824v1). That work was revised and published as two articles, each with its own maintained, tested repository:
>
> - Alarid-Escudero F, et al. An Introductory Tutorial on Cohort State-Transition Models in R Using a Cost-Effectiveness Analysis Example. *Medical Decision Making*. 2023;43(1):3-20. https://doi.org/10.1177/0272989X221103163. Code: [cohort-modeling-tutorial-intro](https://github.com/DARTH-git/cohort-modeling-tutorial-intro)
> - Alarid-Escudero F, et al. A Tutorial on Time-Dependent Cohort State-Transition Models in R using a Cost-Effectiveness Analysis Example. *Medical Decision Making*. 2023;43(1):21-41. https://doi.org/10.1177/0272989X221121747. Code: [cohort-modeling-tutorial-timedep](https://github.com/DARTH-git/cohort-modeling-tutorial-timedep)
>
> Please use those repositories instead.

> **Correction: within-cycle correction weights.** The analysis scripts use `darthtools::gen_wcc()` for Simpson's 1/3rd within-cycle correction. `darthtools` up to v0.4.0 assigned the weights by vector position instead of by time (`1/3, 2/3, 4/3, …` instead of `1/3, 4/3, 2/3, …`). This is fixed in [`darthtools` v0.4.1](https://github.com/DARTH-git/darthtools/releases/tag/v0.4.1), and the scripts now stop if an older `darthtools` is installed. Update it with `remotes::install_github("DARTH-git/darthtools")`. The manuscript files in [`Manuscripts`](Manuscripts) are kept unchanged and contain the same error in their within-cycle correction code chunks. We thank Lauren N. Steimle and [ambtnhlr](https://github.com/ambtnhlr) for reporting this error ([DARTH-git/cohort-modeling-tutorial-intro#2](https://github.com/DARTH-git/cohort-modeling-tutorial-intro/issues/2)).

This GitHub repository provides the code of the three different implementations of cohort state-transition model (cSTM), explained in the following manuscript: 

- Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM, Pechlivanoglou P, Jalal H. [Cohort state-transition models in R: A Tutorial](http://arxiv.org/abs/2001.07824). arXiv:200107824v1. 2020:1-31.

We recommend to first read the manuscript before using the code. Understanding of the use of the multidimensional dynamics array described in 

- Krijkamp EM, Alarid-Escudero F, Enns EA, Pechlivanoglou P, Hunink MGM, Yang A, Jalal HJ. A multidimensional array representation of state-transition model dynamics. Medical Decision Making, 2020;40(2):242-248. https://doi.org/10.1177/0272989X19893973

is recommended.

The [`analysis`](analysis) folder includes three different scripts corresponding to the different implementations of the cSTM (the [`R`](R) folder has the functions they use):
   - [`STM_01.R`](analysis/STM_01.R): time-independent cSTM 
   - [`STM_02.R`](analysis/STM_02.R): age-dependent cSTM with calculations of epidemiological outcomes
   - [`STM_03.R`](analysis/STM_03.R): age- and history-dependent cSTM

In addition we include the [`data`](data) folder with age-specific mortality for the age-dependent and age- and history-dependent cSTM.

To run the cost-effectiveness analysis, you require The Decision Analytic Modeling Package ([`dampack`](https://github.com/DARTH-git/dampack)), an R package for analyzing and visualizing the health economic outputs of decision models.

## Full list of Contributors:

  * [Fernando Alarid-Escudero](https://github.com/feralaes)
  
  * [Eline Krijkamp](https://github.com/krijkamp) 

  * [Eva Enns](https://github.com/evaenns)
  
  * [Alan Yang](https://github.com/alanyang0924)  
  
  * [Myriam Hunink](http://www.erasmus-epidemiology.nl/people/profile.php?id=45)
  
  * [Petros Pechlivanoglou](https://github.com/ppehli)
  
  * [Hawre Jalal](https://github.com/hjalal)
