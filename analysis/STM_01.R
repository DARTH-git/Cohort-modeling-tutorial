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
n_age_init  <- 25                       # age at baseline
n_age_max   <- 100                      # maximum age of follow up
n_t         <- n_age_max - n_age_init   # time horizon, number of cycles
v_names_states <- c("H", "S1", "S2", "D")  # the 4 health states of the model:
                                           # Healthy (H), Sick (S1), Sicker (S2), Dead (D)
v_hcc       <- rep(1, n_t + 1)          # vector of half-cycle correction 
v_hcc[1]    <- v_hcc[n_t + 1] <- 0.5    # half-cycle correction weight 
n_states    <- length(v_names_states)   # number of health states 
d_c         <- 0.03                     # discount rate for costs 
d_e         <- 0.03                     # discount rate for QALYs
v_names_str <- c("SoC", "A", "B", "AB") # store the strategy names
n_str       <- length(v_names_str)      # number of strategies

## Transition probabilities (per cycle), hazard ratios and odds ratio
p_HD        <- 0.002 # constant probability of dying when Healthy (all-cause mortality)
p_HS1       <- 0.15  # probability to become Sick when Healthy conditional on surviving
p_S1H       <- 0.5   # probability to become Healthy when Sick conditional on surviving
p_S1S2      <- 0.105 # probability to become Sicker when Sick conditional on surviving
hr_S1       <- 3     # hazard ratio of death in Sick vs Healthy 
hr_S2       <- 10    # hazard ratio of death in Sicker vs Healthy 
# For treatment B 
or_S1S2     <- 0.6   # odds ratio of becoming Sicker when Sick under B 

## State rewards
# Costs
c_H    <- 2000  # cost of remaining one cycle Healthy 
c_S1   <- 4000  # cost of remaining one cycle Sick 
c_S2   <- 15000 # cost of remaining one cycle Sicker 
c_D    <- 0     # cost of being dead (per cycle)
c_trtA <- 12000 # cost of treatment A
c_trtB <- 13000 # cost of treatment B
# Utilities
u_H    <- 1     # utility when Healthy 
u_S1   <- 0.75  # utility when Sick 
u_S2   <- 0.5   # utility when Sicker
u_D    <- 0     # utility when Healthy 
u_trtA <- 0.95  # utility when being treated

## Transition rewards
du_HS1 <- 0.01  # disutility when transitioning from Healthy to Sick
ic_HS1 <- 1000  # increase in cost when transitioning from Healthy to Sick
ic_D   <- 2000  # increase in cost when dying

# Discount weight (equal discounting is assumed for costs and effects)
v_dwc <- 1 / ((1 + d_e) ^ (0:n_t))
v_dwe <- 1 / ((1 + d_c) ^ (0:n_t))

### Process model inputs
## Age-specific transition probabilities to the Dead state
# compute mortality rates
r_HD        <- prob_to_rate(p_HD)     # hazard rate of dying when Healthy
r_S1D       <- r_HD * hr_S1           # hazard rate of dying when Sick
r_S2D       <- r_HD * hr_S2           # hazard rate of dying when Sicker
# transform rates to probabilities
p_S1D       <- rate_to_prob(r_S1D)    # probability of dying when Sick
p_S2D       <- rate_to_prob(r_S2D)    # probability of dying when Sicker

## Transition probability of becoming Sicker when Sick for treatment B
# transform odds ratios to probabilites 
lor_S1S2    <- log(or_S1S2)               # log-odds ratio of becoming Sicker when Sick
logit_S1S2  <- boot::logit(p_S1S2)        # log-odds of becoming Sicker when Sick
p_S1S2_trtB <- boot::inv.logit(logit_S1S2 +
                               lor_S1S2)  # probability to become Sicker when Sick 
# under B conditional on surviving

###################### Construct state-transition models ##################### 
## Initial state vector
# All starting healthy
v_s_init <- c(H = 1, S1 = 0, S2 = 0, D = 0) # initial state vector
v_s_init

## Initialize cohort trace for cSTM for strategies SoC and A
m_M <- matrix(0, 
              nrow = (n_t + 1), ncol = n_states, 
              dimnames = list(0:n_t, v_names_states))
# Store the initial state vector in the first row of the cohort trace
m_M[1, ] <- v_s_init
## Initialize cohort trace for strategies B and AB
m_M_strB <- m_M # structure and initial states remain the same.

## Initialize transition probability matrix 
# all transitions to a non-death state are assumed to be conditional on survival 
m_P <- matrix(0, 
              nrow = n_states, ncol = n_states, 
              dimnames = list(v_names_states, v_names_states)) # define row and column names
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
check_sum_of_transition_array(m_P,      n_states = n_states, n_t = n_t, verbose = TRUE)
check_sum_of_transition_array(m_P_strB, n_states = n_states, n_t = n_t, verbose = TRUE)

## Initialize transition array which will capture transitions from each state to another over time 
# for strategies SoC and A
a_A <- array(0,
             dim = c(n_states, n_states, n_t + 1),
             dimnames = list(v_names_states, v_names_states, 0:n_t))
# Set first slice of A with the initial state vector in its diagonal
diag(a_A[, , 1]) <- v_s_init
# For strategies B and AB, the array structure and initial state are identical 
a_A_strB <- a_A

#### Run Markov model ####
# Iterative solution of time-independent cSTM
for(t in 1:n_t){
  ## Fill in cohort trace
  # For strategies SoC and A
  m_M[t + 1, ]      <- m_M[t, ]      %*% m_P
  # For strategies B and AB
  m_M_strB[t + 1, ] <- m_M_strB[t, ] %*% m_P_strB
  
  ## Fill in transition dynamics array
  # For strategies SoC and A
  a_A[, , t + 1]      <- m_M[t, ]      * m_P
  # For strategies B and AB
  a_A_strB[, , t + 1] <- m_M_strB[t, ] * m_P_strB
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
              A =  v_u_strA,
              B =  v_u_strB,
              AB = v_u_strAB)
## Store the vectors of state cost for each strategy in a list 
l_c   <- list(SQ = v_c_SoC,
              A =  v_c_strA,
              B =  v_c_strB,
              AB = v_c_strAB)
## Store the transition array for each strategy in a list
l_a_A <- list(SQ = a_A,
              A =  a_A,
              B =  a_A_strB,
              AB = a_A_strB)

# assign strategy names to matching items in the lists
names(l_u) <- names(l_c) <- names(l_a_A) <- v_names_str

## create empty vectors to store total utilities and costs 
v_tot_qaly <- v_tot_cost <- vector(mode = "numeric", length = n_str)
names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str

#### Loop through each strategy and calculate total utilities and costs ####
for (i in 1:n_str) {
  v_u_str <- l_u[[i]]   # select the vector of state utilities for the ith strategy
  v_c_str <- l_c[[i]]   # select the vector of state costs for the ith strategy
  a_A_str <- l_a_A[[i]] # select the transition array for the ith strategy
  
  #### Array of state rewards ####
  # Create transition matrices of state utilities and state costs for the ith strategy 
  m_u_str   <- matrix(v_u_str, nrow = n_states, ncol = n_states, byrow = T)
  m_c_str   <- matrix(v_c_str, nrow = n_states, ncol = n_states, byrow = T)
  # Expand the transition matrix of state utilities across cycles to form a transition array of state utilities
  a_R_u_str <- array(m_u_str, 
                     dim      = c(n_states, n_states, n_t + 1),
                     dimnames = list(v_names_states, v_names_states, 0:n_t))
  # Expand the transition matrix of state costs across cycles to form a transition array of state costs
  a_R_c_str <- array(m_c_str, 
                     dim      = c(n_states, n_states, n_t + 1),
                     dimnames = list(v_names_states, v_names_states, 0:n_t))
  
  #### Apply transition rewards ####  
  # Apply disutility due to transition from H to S1
  a_R_u_str["H", "S1", ]      <- a_R_u_str["H", "S1", ]       - du_HS1
  # Add transition cost per cycle due to transition from H to S1
  a_R_c_str["H", "S1", ]      <- a_R_c_str["H", "S1", ]       + ic_HS1
  # Add transition cost  per cycle of dying from all non-dead states
  a_R_c_str[-n_states, "D", ] <- a_R_c_str[- n_states, "D", ] + ic_D
  
  #### Expected QALYs and Costs for all transitions per cycle ####
  # QALYs = life years x QoL
  # Note: all parameters are annual in our example. In case your own case example is different make sure you correctly apply .
  a_Y_c_str <- a_A_str * a_R_c_str
  a_Y_u_str <- a_A_str * a_R_u_str 
  
  #### Expected QALYs and costs per cycle ####
  ## Vector of QALYs and Costs
  v_qaly_str <- apply(a_Y_u_str, 3, sum) # sum the proportion of the cohort across transitions 
  v_cost_str <- apply(a_Y_c_str, 3, sum) # sum the proportion of the cohort across transitions
  
  #### Discounted total expected QALYs and Costs per strategy and apply half-cycle correction if applicable ####
  ## QALYs
  v_tot_qaly[i] <- t(v_qaly_str) %*% (v_dwe * v_hcc)
  ## Costs
  v_tot_cost[i] <- t(v_cost_str) %*% (v_dwc * v_hcc)
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

