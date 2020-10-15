# Cohort state-transition models in R: A Tutorial
This GitHub repository provides the code of the three different implementations of cohort state-transition model (cSTM), explained in the following manuscript: 

- Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM, Pechlivanoglou P, Jalal H. [Cohort state-transition models in R: A Tutorial](http://arxiv.org/abs/2001.07824). arXiv:200107824v1. 2020:1-31.

We recommend to first read the manuscript before using the code. Understanding of the use of the multidimensional dynamics array described in 

- Krijkamp EM, Alarid-Escudero F, Enns EA, Pechlivanoglou P, Hunink MGM, Yang A, Jalal HJ. A multidimensional array representation of state-transition model dynamics. Medical Decision Making, 2020;40(2):242-248. https://doi.org/10.1177/0272989X19893973

is recommended.

The [`R`](https://github.com/DARTH-git/Cohort-modeling-tutorial/tree/master/R) folder includes three different scripts corresponding to the different implementations of the cSTM:
   - [`STM_01.R`](https://github.com/DARTH-git/Cohort-modeling-tutorial/blob/master/R/STM_01.R): time-independent cSTM 
   - [`STM_02.R`](https://github.com/DARTH-git/Cohort-modeling-tutorial/blob/master/R/STM_02.R): age-dependent cSTM with calculations of epidemiological outcomes
   - [`STM_03.R`](https://github.com/DARTH-git/Cohort-modeling-tutorial/blob/master/R/STM_03.R): age- and history-dependent cSTM

In addition we include the [`data`](https://github.com/DARTH-git/Cohort-modeling-tutorial/tree/master/data) folder with age-specific mortality for the age-dependent and age- and history-dependent cSTM.

To run the cost-effectiveness analysis, you require The Decision Analytic Modeling Package ([`dampack`](https://github.com/DARTH-git/dampack)), an R package for analyzing and visualizing the health economic outputs of decision models.

## Full list of Contributors:

  * [Fernando Alarid-Escudero](https://github.com/feralaes)
  
  * [Eline Krijkamp](https://github.com/krijkamp) 

  * [Eva Enns](https://github.com/evaenns)
  
  * [Alan Yang](https://github.com/alanyang0924)  
  
  * [Myriam Hunink](http://www.erasmus-epidemiology.nl/people/profile.php?id=45)
  
  * [Petros Pechlivanoglou](https://github.com/ppehli)
  
  * [Hawre Jalal](https://github.com/hjalal)
