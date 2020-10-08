##############################################################################
### Implementation of cohort state-transition models in R:                  ##
### From conceptualization to implementation 2020                           ##
##############################################################################
# This code forms the basis for the state-transition model of the article: 
# 'Implementation of cohort state-transition models in R' 
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
# Usual care: best available care for the patients with the disease. This scenario reflects the natural history of the disease progressions
# New treatment 1: this new treatment is given to all sick patients, patients in sick and sicker, but does only improves the utility of those being sick.
# New treatment 2: the new treatment reduces disease progression from sick to sicker. However, it is not possible to distinguish those sick from sicker and therefore all individuals in one of the two sick states get the treatment.  
# New treatment 1 & new treatment 2: This strategy combines the new treatment 1 and new treatment 2. The disease progression is reduced and Sick individuals has an improved utility. 

################################ Initial setup ############################### 
rm(list = ls())    # remove any variables in R's memory 

### Install packages
# install.packages("dplyr")      # to manipulate data
# install.packages("data.table") # to manipulate data
# install.packages("tidyr")      # to manipulate data
# install.packages("reshape2")   # to manipulate data
# install.packages("ggplot2")    # to visualize data
# install.packages("scales")     # for dollar signs and commas
# install.packages("boot")       # to handle log odds and log odds ratios
# install.packages("devtools")   # to ensure compatibility among packages
# devtools::install_github("DARTH-git/dampack") # to install dampack form GitHub, for CEA and calculate ICERs

### Load packages
library(dplyr)    
library(data.table)
library(tidyr)
library(reshape2) 
library(ggplot2)   
library(scales)    
library(boot)
library(dampack) 

### Load supplementary functions
source("R/Functions.R")

################################ Model input ################################# 
## General setup
n_age_init  <- 25                      # age at baseline
n_age_max   <- 100                     # maximum age of follow up
n_t         <- n_age_max - n_age_init  # time horizon, number of cycles
v_n         <- c("H", "S1", "S2", "D") # the 4 health states of the model:
# Healthy (H), Sick (S1), Sicker (S2), Dead (D)
v_hcc       <- rep(1, n_t + 1)         # vector of half-cycle correction 
v_hcc[1]    <- v_hcc[n_t + 1] <- 0.5   # half-cycle correction weight 
n_states    <- length(v_n)             # number of health states 
d_c         <- 0.03                    # discount rate for costs 
d_e         <- 0.03                    # discount rate for QALYs
v_names_str <- c("Usual care",         # store the strategy names
                 "New treatment 1", 
                 "New treatment 2", 
                 "New treatments 1 & 2") 
n_str       <- length(v_names_str) # number of strategies

## Transition probabilities (per cycle), hazard ratios and odds ratio
p_HD        <- 0.002                  # constant probability of dying when Healthy (all-cause mortality)
p_HS1       <- 0.15                   # probability to become Sick when Healthy conditional on surviving the cycle
p_S1H       <- 0.5                    # probability to become Healthy when Sick conditional on surviving the cycle
p_S1S2      <- 0.105                  # probability to become Sicker when Sick conditional on surviving the cycle
hr_S1       <- 3                      # hazard ratio of death in Sick vs Healthy
hr_S2       <- 10                     # hazard ratio of death in Sicker vs Healthy 
# For New treatment 2
or_S1S2     <- 0.6                    # odds ratio of becoming Sicker when Sick under New treatment 2

## State rewards
# Costs
c_H    <- 2000  # cost of remaining one cycle Healthy 
c_S1   <- 4000  # cost of remaining one cycle Sick 
c_S2   <- 15000 # cost of remaining one cycle Sicker 
c_D    <- 0     # cost of being dead (per cycle)
c_trt1 <- 12000 # cost of New treatment 1 (per cycle) 
c_trt2 <- 13000 # cost of New treatment 2 (per cycle)
# Utilities
u_H    <- 1     # utility when Healthy 
u_S1   <- 0.75  # utility when Sick 
u_S2   <- 0.5   # utility when Sicker
u_D    <- 0     # utility when Healthy 
u_trt1 <- 0.95  # utility when being treated

# Discount weight (equal discounting is assumed for costs and effects)
v_dwc  <- 1 / ((1 + d_e) ^ (0:(n_t)))
v_dwe  <- 1 / ((1 + d_c) ^ (0:(n_t)))

## Process model inputs
# compute mortality rates
r_HD        <- prob_to_rate(p_HD)     # hazard rate of dying when Healthy
r_S1D       <- r_HD * hr_S1           # hazard rate of dying when Sick
r_S2D       <- r_HD * hr_S2           # hazard rate of dying when Sicker
# transform rates to probabilities
p_S1D       <- rate_to_prob(r_S1D)    # probability of dying when Sick
p_S2D       <- rate_to_prob(r_S2D)    # probability of dying when Sicker
# transform odds ratios to probabilites
lor_S1S2    <- log(or_S1S2)           # log-odds ratio of becoming Sicker when Sick
logit_S1S2  <- logit(p_S1S2)          # log-odds of becoming Sicker when Sick
p_S1S2_trt2 <- inv.logit(logit_S1S2 +
                           lor_S1S2)  # probability to become Sicker when Sick 
# under New treatment 2 conditional on surviving

###################### Construct state-transition models ##################### 
## Initial state vector
# All starting healthy
v_s_init <- c(H = 1, S1 = 0, S2 = 0, D = 0) # initial state vector
v_s_init

## Initialize cohort trace
m_M <- matrix(0, 
              nrow = (n_t + 1), ncol = n_states, 
              dimnames = list(0:n_t, v_n))
# Store the initial state vector in the first row of the cohort trace
m_M[1, ] <- v_s_init
# For treatment 2
m_M_trt2 <- m_M

## Initialize transition probability matrix
m_P <- matrix(0, 
              nrow = n_states, ncol = n_states, 
              dimnames = list(v_n, v_n)) # define row and column names
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

# For New treatment 2
# Only need to update the transition probabilities from S1 involving p_S1S2
m_P_trt2 <- m_P
m_P_trt2["S1", "S1"] <- (1 - p_S1D) * (1 - (p_S1H + p_S1S2_trt2))
m_P_trt2["S1", "S2"] <- (1 - p_S1D) * p_S1S2_trt2

### Check if transition probability matrix is valid (i.e., elements cannot < 0 or > 1) 
check_transition_probability(m_P,      verbose = TRUE)
check_transition_probability(m_P_trt2, verbose = TRUE)
### Check if transition probability matrix sum to 1 (i.e., each row should sum to 1)
check_sum_of_transition_array(m_P,      n_states = n_states, n_t = n_t, verbose = TRUE)
check_sum_of_transition_array(m_P_trt2, n_states = n_states, n_t = n_t, verbose = TRUE)


#### Run Markov model ####
# Iterative solution of time-independent cSTM
for(t in 1:n_t){
  # Fill in cohort trace
  m_M[t + 1, ] <- m_M[t, ] %*% m_P
  m_M_trt2[t + 1, ] <- m_M_trt2[t, ] %*% m_P_trt2

}

## Store the cohort traces in a list
l_m_M <- list(m_M,
              m_M,
              m_M_trt2,
              m_M_trt2)
names(l_m_M) <- v_names_str

#### Plot Outputs ####
### Plot the cohort trace for Usual care and New Treatment 1
plot_trace(m_M)

#### State Rewards ####
## Vector of state utilities under Usual care
v_u_UC     <- c(H  = u_H, 
                S1 = u_S1, 
                S2 = u_S2, 
                D  = u_D)
## Vector of state costs per cycle under Usual care
v_c_UC     <- c(H  = c_H, 
                S1 = c_S1,
                S2 = c_S2, 
                D  = c_D)
## Vector of state utilities under New treatment 1
v_u_trt1   <- c(H  = u_H, 
                S1 = u_trt1, 
                S2 = u_S2, 
                D  = u_D)
## Vector of state costs per cycle under New treatment 1
v_c_trt1   <- c(H  = c_H, 
                S1 = c_S1 + c_trt1,
                S2 = c_S2 + c_trt1, 
                D  = c_D)
## Vector of state utilities under New treatment 2
v_u_trt2   <- c(H  = u_H, 
                S1 = u_S1, 
                S2 = u_S2, 
                D  = u_D)
## Vector of state costs per cycle under New treatment 2
v_c_trt2   <- c(H  = c_H, 
                S1 = c_S1 + c_trt2, 
                S2 = c_S2 + c_trt2, 
                D  = c_D)
## Vector of state utilities under New treatment 1 & 2
v_u_trt1_2 <- c(H  = u_H, 
                S1 = u_trt1, 
                S2 = u_S2, 
                D  = u_D)
## Vector of state costs per cycle under New treatment 1 & 2
v_c_trt1_2 <- c(H  = c_H, 
                S1 = c_S1 + (c_trt1 + c_trt2), 
                S2 = c_S2 + (c_trt1 + c_trt2), 
                D  = c_D)

## Store the vectors of state utilities for each strategy in a list 
l_u   <- list(v_u_UC,
              v_u_trt1,
              v_u_trt2,
              v_u_trt1_2)
## Store the vectors of state cost for each strategy in a list 
l_c   <- list(v_c_UC,
              v_c_trt1,
              v_c_trt2,
              v_c_trt1_2)
## Store the transition array for each strategy in a list
l_m_M <- list(m_M,
              m_M,
              m_M_trt2,
              m_M_trt2)
# assign strategy names to matching items in the lists
names(l_u) <- names(l_c) <- names(l_m_M) <- v_names_str

## create empty vectors to store total utilities and costs 
v_tot_qaly <- v_tot_cost <- vector(mode = "numeric", length = n_str)
names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str 

#### Loop through each strategy and calculate total utilities and costs ####
for (i in 1:n_str) {
  v_u     <- l_u[[i]]   # select the vector of state utilities for the ith strategy
  v_c     <- l_c[[i]]   # select the vector of state costs for the ith strategy
  m_M_str <- l_m_M[[i]] # select the cohort trace for the ith strategy
  
  #### Expected QALYs and Costs per cycle ####
  ## Vector of QALYs under Usual Care
  v_qaly <- m_M_str %*% v_u
  ## Vector of costs under Usual Care
  v_cost <- m_M_str %*% v_c
  
  #### Discounted total expected QALYs and Costs per strategy ####
  ### Apply half-cycle correction if applicable 
  ## QALYs
  v_tot_qaly[i] <- t(v_qaly) %*% (v_dwe * v_hcc)
  ## Costs
  v_tot_cost[i] <- t(v_cost) %*% (v_dwc * v_hcc)
}

######################### Cost-effectiveness analysis ######################## 
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
