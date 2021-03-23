##############################################################################
### Cohort State-Transition Models in R                                    ###
##############################################################################
# This code forms the basis for the state-transition model of the article: 
# 'cohort State-transition Models in R' 
# Authors: 
# - Fernando Alarid-Escudero <fernando.alarid@cide.edu>
# - Eline Krijkamp
# - Eva A. Enns
# - Alan Yang
# - M.G. Myriam Hunink
# - Petros Pechlivanoglou
# - Hawre Jalal
# Please cite the article when using this code
#
# To program this tutorial we made use of 
# R version 4.0.2 (2020-06-22)
# Platform: 64-bit operating system, x64-based processor
# Running under: Windows 10
# RStudio: Version 1.3.1073 2009-2020 RStudio, Inc

##############################################################################
############################# Code of Appendix ############################### 
##############################################################################
# Implements a time-independent Sick-Sicker cSTM model                       #
# Standard of Care (SoC): best available care for the patients with the disease. This scenario reflects the natural history of the disease progressions
# Strategy A: treatment A is given to all sick patients, patients in sick and sicker, but does only improves the utility of those being sick.
# Strategy B: treatment B reduces disease progression from sick to sicker. However, it is not possible to distinguish those sick from sicker and therefore all individuals in one of the two sick states get the treatment.  
# Strategy AB: This strategy combines treatment A and treatment B. The disease progression is reduced and Sick individuals has an improved utility. 

################################ Initial setup ############################### 
rm(list = ls())    # remove any variables in R's memory 

### Install packages
# install.packages("dplyr")      # to manipulate data
# install.packages("data.table") # to manipulate data
# install.packages("tidyr")      # to manipulate data
# install.packages("reshape2")   # to manipulate data
# install.packages("ggplot2")    # to visualize data
# install.packages("gridExtra")  # to visualize data
# install.packages("scales")     # for dollar signs and commas
# install.packages("boot")       # to handle log odds and log odds ratios
# install.packages("devtools")   # to ensure compatibility among packages
# install.packages("dampack")    # for CEA and calculate ICERs
# devtools::install_github("DARTH-git/darthtools") # to install darthtools from GitHub
# install.packages("doParallel") # to handle parallel processing

### Load packages
library(dplyr)    
library(data.table)
library(tidyr)
library(reshape2) 
library(ggplot2) 
library(gridExtra)
library(scales)    
library(boot)
library(dampack) 
library(darthtools)
library(doParallel)

### Load supplementary functions
source("R/Functions.R")

################################ Model input ################################# 
## General setup
n_age_init  <- 25                         # age at baseline
n_age_max   <- 100                        # maximum age of follow up
n_cycles    <- n_age_max - n_age_init     # time horizon, number of cycles
v_names_states <- c("H", "S1", "S2", "D") # the 4 health states of the model:
                                          # Healthy (H), Sick (S1), Sicker (S2), Dead (D)
n_states    <- length(v_names_states)     # number of health states 

# Discounting factors
d_c         <- 0.03                       # discount rate for costs 
d_e         <- 0.03                       # discount rate for QALYs

# Discount weight (equal discounting is assumed for costs and effects)
v_dwc  <- 1 / ((1 + d_e) ^ (0:n_cycles))
v_dwe  <- 1 / ((1 + d_c) ^ (0:n_cycles))

# Strategies
v_names_str <- c("Standard of care",      # store the strategy names
                 "Strategy A", 
                 "Strategy B",
                 "Strategy AB") 
n_str       <- length(v_names_str)        # number of strategies

# Within-cycle correction (WCC) using Simpson's 1/3 rule
v_wcc <- darthtools::gen_wcc(n_cycles = n_cycles, 
                             method = "Simpson1/3") # vector of wcc

###################### Probabilistic Sensitivity Analysis (PSA) ##################### 
# List of input parameters
l_params_all <- as.list(data.frame(
  # Transition probabilities (per cycle), hazard ratios
  r_HD        = 0.002, # constant rate of dying when Healthy (all-cause mortality)
  p_HS1       = 0.15,  # probability to become Sick when Healthy conditional on surviving
  p_S1H       = 0.5,   # probability to become Healthy when Sick conditional on surviving
  p_S1S2      = 0.105, # probability to become Sicker when Sick conditional on surviving
  hr_S1       = 3,     # hazard ratio of death in Sick vs Healthy 
  hr_S2       = 10,    # hazard ratio of death in Sicker vs Healthy 
  # Effectiveness of treatment B 
  hr_S1S2_trtB = 0.6,  # hazard ratio of becoming Sicker when Sick under B under treatment B
  ## State rewards
  # Costs
  c_H    = 2000,  # cost of remaining one cycle in Healthy 
  c_S1   = 4000,  # cost of remaining one cycle in Sick 
  c_S2   = 15000, # cost of remaining one cycle in Sicker 
  c_D    = 0,     # cost of being dead (per cycle)
  c_trtA = 12000, # cost of treatment A
  c_trtB = 13000, # cost of treatment B
  # Utilities
  u_H    = 1,     # utility when Healthy 
  u_S1   = 0.75,  # utility when Sick 
  u_S2   = 0.5,   # utility when Sicker
  u_D    = 0,     # utility when Dead 
  u_trtA = 0.95   # utility when being treated with A
))

# store the parameter names into a vector
v_names_params <- names(l_params_all)

# Load model function
source("R/Functions STM_01.R")
# Test function
calculate_ce_out(l_params_all)

# Function to generate PSA input dataset
gen_psa <- function(n_sim = 1000, seed = 071818){
  set.seed(seed) # set a seed to be able to reproduce the same results
  df_psa <- data.frame(
    # Transition probabilities (per cycle), hazard ratios
    r_HD    = rlnorm(n_sim, meanlog = log(0.002), sdlog = 0.01), # constant rate of dying when Healthy (all-cause mortality)
    p_HS1   = rbeta(n_sim, shape1 = 30, shape2 = 170),        # probability to become Sick when Healthy conditional on surviving
    p_S1H   = rbeta(n_sim, shape1 = 60, shape2 = 60) ,        # probability to become Healthy when Sick conditional on surviving
    p_S1S2  = rbeta(n_sim, shape1 = 84, shape2 = 716),        # probability to become Sicker when Sick conditional on surviving
    hr_S1   = rlnorm(n_sim, meanlog = log(3),  sdlog = 0.01), # hazard ratio of death in Sick vs Healthy 
    hr_S2   = rlnorm(n_sim, meanlog = log(10), sdlog = 0.02), # hazard ratio of death in Sicker vs Healthy 
    
    # Effectiveness of treatment B 
    hr_S1S2_trtB = rlnorm(n_sim, meanlog = log(0.6), sdlog = 0.02), # hazard ratio of becoming Sicker when Sick under B under treatment B
    
    # State rewards
    # Costs
    c_H     = rgamma(n_sim, shape = 100, scale = 20),     # cost of remaining one cycle in Healthy 
    c_S1    = rgamma(n_sim, shape = 177.8, scale = 22.5), # cost of remaining one cycle in Sick 
    c_S2    = rgamma(n_sim, shape = 225, scale = 66.7),   # cost of remaining one cycle in Sicker 
    c_D     = 0,                                          # cost of being dead (per cycle)
    c_trtA  = rgamma(n_sim, shape = 73.5, scale = 163.3), # cost of treatment A
    c_trtB  = rgamma(n_sim, shape = 86.2, scale = 150.8), # cost of treatment B
    
    # Utilities
    u_H     = rbeta(n_sim, shape1 = 200, shape2 = 3),     # utility when Healthy 
    u_S1    = rbeta(n_sim, shape1 = 130, shape2 = 45),    # utility when Sick 
    u_S2    = rbeta(n_sim, shape1 = 230, shape2 = 230),   # utility when Sicker
    u_D     = 0,                                          # utility when Dead 
    u_trtA  = rbeta(n_sim, shape1 = 300, shape2 = 15)     # utility when being treated with A
  )
  return(df_psa)
}
# Try it
gen_psa(10) 

# Number of simulations
n_sim <- 1000

# Generate PSA input dataset
df_psa_input <- gen_psa(n_sim = n_sim)
# First six observations
head(df_psa_input)

# Histogram of parameters
ggplot(melt(df_psa_input, variable.name = "Parameter"), aes(x = value)) +
  facet_wrap(~Parameter, scales = "free") +
  geom_histogram(aes(y = ..density..)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) + 
  theme_bw(base_size = 16) + 
  theme(axis.text = element_text(size=6)) 

# Initialize dataframes with PSA output 
# Dataframe of costs
df_c <- as.data.frame(matrix(0, 
                             nrow = n_sim,
                             ncol = n_str))
colnames(df_c) <- v_names_str
# Dataframe of effectiveness
df_e <- as.data.frame(matrix(0, 
                             nrow = n_sim,
                             ncol = n_str))
colnames(df_e) <- v_names_str

## Conduct probabilistic sensitivity analysis
# Run Markov model on each parameter set of PSA input dataset
for(i in 1:n_sim){
  l_out_temp <- calculate_ce_out(df_psa_input[i, ])
  df_c[i, ] <- l_out_temp$Cost
  df_e[i, ] <- l_out_temp$Effect
  # Display simulation progress
  if(i/(n_sim/10) == round(i/(n_sim/10), 0)) { # display progress every 10%
    cat('\r', paste(i/n_sim * 100, "% done", sep = " "))
  }
}

## Visualize PSA results and CEA
# Create PSA object for dampack
l_psa <- make_psa_obj(cost          = df_c, 
                      effectiveness = df_e, 
                      parameters    = df_psa_input, 
                      strategies    = v_names_str)

# Create PSA graphs
# Vector with willingness-to-pay (WTP) thresholds.
v_wtp <- seq(0, 300000, by = 10000)

# Cost-Effectiveness Scatter plot
plot(l_psa)

# Conduct CEA with probabilistic output
# Compute expected costs and effects for each strategy from the PSA
df_out_ce_psa <- summary(l_psa)

# Calculate incremental cost-effectiveness ratios (ICERs)
df_cea_psa <- calculate_icers(cost       = df_out_ce_psa$meanCost, 
                              effect     = df_out_ce_psa$meanEffect,
                              strategies = df_out_ce_psa$Strategy)
df_cea_psa

# Plot cost-effectiveness frontier
plot(df_cea_psa)

# Cost-effectiveness acceptability curves (CEACs) and frontier (CEAF)
ceac_obj <- ceac(wtp = v_wtp, psa = l_psa)
# Regions of highest probability of cost-effectiveness for each strategy
summary(ceac_obj)
# CEAC & CEAF plot
plot(ceac_obj)

# Expected Loss Curves (ELCs)
elc_obj <- calc_exp_loss(wtp = v_wtp, psa = l_psa)
elc_obj
# ELC plot
plot(elc_obj, log_y = FALSE)

# Expected value of perfect information (EVPI)
evpi <- calc_evpi(wtp = v_wtp, psa = l_psa)
# EVPI plot
plot(evpi, effect_units = "QALY")

