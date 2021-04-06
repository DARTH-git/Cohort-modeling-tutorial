##############################################################################
### Cohort State-Transition Models in R                                    ###
##############################################################################
# This code forms the basis for the state-transition model of the article: 
# 'An Introductory Tutorial to Cohort State-Transition Models in R' 
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

# Strategies
v_names_str <- c("Standard of care",      # store the strategy names
                 "Strategy A", 
                 "Strategy B",
                 "Strategy AB") 
n_str       <- length(v_names_str)        # number of strategies

# Within-cycle correction (WCC) using Simpson's 1/3 rule
v_wcc <- darthtools::gen_wcc(n_cycles = n_cycles, 
                             method = "Simpson1/3") # vector of wcc

## Transition probabilities (per cycle), hazard ratios 
r_HD   <- 0.002 # constant rate of dying when Healthy (all-cause mortality)
p_HS1  <- 0.15  # probability to become Sick when Healthy conditional on surviving
p_S1H  <- 0.5   # probability to become Healthy when Sick conditional on surviving
p_S1S2 <- 0.105 # probability to become Sicker when Sick conditional on surviving
hr_S1  <- 3     # hazard ratio of death in Sick vs Healthy 
hr_S2  <- 10    # hazard ratio of death in Sicker vs Healthy 

# Effectiveness of treatment B 
hr_S1S2_trtB <- 0.6  # hazard ratio of becoming Sicker when Sick under treatment B

## State rewards
# Costs
c_H    <- 2000  # cost of remaining one cycle in Healthy 
c_S1   <- 4000  # cost of remaining one cycle in Sick 
c_S2   <- 15000 # cost of remaining one cycle in Sicker 
c_D    <- 0     # cost of being dead (per cycle)
c_trtA <- 12000 # cost of treatment A
c_trtB <- 13000 # cost of treatment B
# Utilities
u_H    <- 1     # utility when Healthy 
u_S1   <- 0.75  # utility when Sick 
u_S2   <- 0.5   # utility when Sicker
u_D    <- 0     # utility when Dead 
u_trtA <- 0.95  # utility when being treated with A

# Discount weight (equal discounting is assumed for costs and effects)
v_dwc  <- 1 / ((1 + d_e) ^ (0:n_cycles))
v_dwe  <- 1 / ((1 + d_c) ^ (0:n_cycles))

### Process model inputs
## Transition probabilities to the Dead state
# compute mortality rates
r_S1D <- r_HD * hr_S1        # Mortality in the Sick state
r_S2D <- r_HD * hr_S2        # Mortality in the Sick state
# transform rates to probabilities
p_HD  <- rate_to_prob(r_HD)  # Mortality risk in the Healthy state
p_S1D <- rate_to_prob(r_S1D) # Mortality risk in the Sick state
p_S2D <- rate_to_prob(r_S2D) # Mortality risk in the Sicker state

## Transition probability of becoming Sicker when Sick for treatment B
# transform probability to rate
r_S1S2      <- prob_to_rate(p = p_S1S2)
# apply hazard ratio to rate to obtain transition rate of becoming Sicker when Sick for treatment B
r_S1S2_trtB <- r_S1S2 * hr_S1S2_trtB
# transform rate to probability
p_S1S2_trtB <- rate_to_prob(r = r_S1S2_trtB) # probability to become Sicker when Sick 
                                             # under treatment B conditional on surviving

####################### Construct state-transition models ######################
## Initial state vector
# All starting healthy
v_s_init <- c(H = 1, S1 = 0, S2 = 0, D = 0) # initial state vector
v_s_init

## Initialize cohort trace for cSTM for strategies SoC and A
m_M <- matrix(0, 
              nrow = (n_cycles + 1), ncol = n_states,
              dimnames = list(0:n_cycles, v_names_states))
# Store the initial state vector in the first row of the cohort trace
m_M[1, ] <- v_s_init
## Initialize cohort trace for strategies B and AB
m_M_strB <- m_M # structure and initial states remain the same.

## Initialize transition probability matrix 
# all transitions to a non-death state are assumed to be conditional on survival 
m_P <- matrix(0, 
              nrow = n_states, ncol = n_states, 
              dimnames = list(v_names_states, 
                              v_names_states)) # define row and column names
## Fill in matrix
# From H
m_P["H", "H"]   <- (1 - p_HD) * (1 - p_HS1)
m_P["H", "S1"]  <- (1 - p_HD) * p_HS1 
m_P["H", "D"]   <- p_HD
# From S1
m_P["S1", "H"]  <- (1 - p_S1D) * p_S1H
m_P["S1", "S1"] <- (1 - p_S1D) * (1 - (p_S1H + p_S1S2))
m_P["S1", "S2"] <- (1 - p_S1D) * p_S1S2
m_P["S1", "D"]  <- p_S1D
# From S2
m_P["S2", "S2"] <- 1 - p_S2D
m_P["S2", "D"]  <- p_S2D
# From D
m_P["D", "D"]   <- 1

### For strategies B and AB
## Initialize transition probability array for strategies B and AB
m_P_strB <- m_P
## Only need to update the probabilities involving the transition from Sick to Sicker, p_S1S2
# From S1
m_P_strB["S1", "S1"] <- (1 - p_S1D) * (1 - (p_S1H + p_S1S2_trtB))
m_P_strB["S1", "S2"] <- (1 - p_S1D) * p_S1S2_trtB

### Check if transition probability matrices are valid
## Check that transition probabilities are [0, 1]
check_transition_probability(m_P,      verbose = TRUE)
check_transition_probability(m_P_strB, verbose = TRUE)
## Check that all rows sum to 1
check_sum_of_transition_array(m_P,      n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
check_sum_of_transition_array(m_P_strB, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)

#### Run Markov model ####
# Iterative solution of time-independent cSTM
for(t in 1:n_cycles){
  ## Fill in cohort trace
  # For strategies SoC and A
  m_M[t + 1, ]      <- m_M[t, ]      %*% m_P
  # For strategies B and AB
  m_M_strB[t + 1, ] <- m_M_strB[t, ] %*% m_P_strB
}

## Store the cohort traces in a list
l_m_M <- list(m_M,
              m_M,
              m_M_strB,
              m_M_strB)
names(l_m_M) <- v_names_str

#### Plot Outputs ####
### Plot the cohort trace for strategies SoC and A
plot_trace(m_M)

#### State Rewards ####
## Vector of state utilities under strategy SoC
v_u_SoC    <- c(H  = u_H, 
                S1 = u_S1, 
                S2 = u_S2, 
                D  = u_D)
## Vector of state costs under strategy SoC
v_c_SoC    <- c(H  = c_H, 
                S1 = c_S1,
                S2 = c_S2, 
                D  = c_D)
## Vector of state utilities under strategy A
v_u_strA   <- c(H  = u_H, 
                S1 = u_trtA, 
                S2 = u_S2, 
                D  = u_D)
## Vector of state costs under strategy A
v_c_strA   <- c(H  = c_H, 
                S1 = c_S1 + c_trtA,
                S2 = c_S2 + c_trtA, 
                D  = c_D)
## Vector of state utilities under strategy B
v_u_strB   <- c(H  = u_H, 
                S1 = u_S1, 
                S2 = u_S2, 
                D  = u_D)
## Vector of state costs under strategy B
v_c_strB   <- c(H  = c_H, 
                S1 = c_S1 + c_trtB, 
                S2 = c_S2 + c_trtB, 
                D  = c_D)
## Vector of state utilities under strategy AB
v_u_strAB  <- c(H  = u_H, 
                S1 = u_trtA, 
                S2 = u_S2, 
                D  = u_D)
## Vector of state costs under strategy AB
v_c_strAB  <- c(H  = c_H, 
                S1 = c_S1 + (c_trtA + c_trtB), 
                S2 = c_S2 + (c_trtA + c_trtB), 
                D  = c_D)

## Store the vectors of state utilities for each strategy in a list 
l_u   <- list(SQ = v_u_SoC,
              A  = v_u_strA,
              B  = v_u_strB,
              AB = v_u_strAB)
## Store the vectors of state cost for each strategy in a list 
l_c   <- list(SQ = v_c_SoC,
              A  = v_c_strA,
              B  = v_c_strB,
              AB = v_c_strAB)

# assign strategy names to matching items in the lists
names(l_u) <- names(l_c) <- v_names_str

## create empty vectors to store total utilities and costs 
v_tot_qaly <- v_tot_cost <- vector(mode = "numeric", length = n_str)
names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str

#### Loop through each strategy and calculate total utilities and costs ####
for (i in 1:n_str) {
  v_u_str <- l_u[[i]]   # select the vector of state utilities for the i-th strategy
  v_c_str <- l_c[[i]]   # select the vector of state costs for the i-th strategy
  
  #### Expected QALYs and costs per cycle ####
  ### Vector of QALYs and Costs
  ## Apply state rewards ###
  v_qaly_str <- l_m_M[[i]] %*% v_u_str # sum the utilities of all states for each cycle
  v_cost_str <- l_m_M[[i]] %*% v_c_str # sum the costs of all states for each cycle
  
  #### Discounted total expected QALYs and Costs per strategy and apply half-cycle correction if applicable ####
  ## QALYs
  v_tot_qaly[i] <- t(v_qaly_str) %*% (v_dwe * v_wcc)
  ## Costs
  v_tot_cost[i] <- t(v_cost_str) %*% (v_dwc * v_wcc)
}

########################## Cost-effectiveness analysis #######################
### Calculate incremental cost-effectiveness ratios (ICERs)
df_cea <- calculate_icers(cost       = v_tot_cost, 
                          effect     = v_tot_qaly,
                          strategies = v_names_str)
df_cea

### Create CEA table in proper format
table_cea <- format_table_cea(df_cea)
table_cea

### CEA frontier
plot(df_cea, label = "all") +
  expand_limits(x = max(table_cea$QALYs) + 0.5) 

###################### Probabilistic Sensitivity Analysis (PSA) ##################### 
### Load model, CEA and PSA functions
source("R/Functions STM_01.R")

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

# Test function to compute CE outcomes
calculate_ce_out(l_params_all)

# Test function to generate PSA input dataset
generate_psa_params(10) 

# Number of simulations
n_sim <- 1000

# Generate PSA input dataset
df_psa_input <- generate_psa_params(n_sim = n_sim)
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
v_wtp <- seq(0, 250000, by = 5000)

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