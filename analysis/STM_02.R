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
# This model incorporates time-dependent transition probabilities 

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
# devtools::install_github("DARTH-git/dampack") # to install dampack form GitHub, for CEA and calculate ICERs

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
p_HS1       <- 0.15  # probability to become Sick when Healthy conditional on surviving
p_S1H       <- 0.5   # probability to become Healthy when Sick conditional on surviving
p_S1S2      <- 0.105 # probability to become Sicker when Sick conditional on surviving
hr_S1       <- 3     # hazard ratio of death in Sick vs Healthy 
hr_S2       <- 10    # hazard ratio of death in Sicker vs Healthy 
# For treatment B 
or_S1S2     <- 0.6   # odds ratio of becoming Sicker when Sick 

## Age-dependent mortality rates
lt_usa_2015 <- read.csv("data/LifeTable_USA_Mx_2015.csv")
v_r_mort_by_age <- lt_usa_2015 %>% 
                   select(Total) %>%
                   as.matrix() # anyone above 100 have the same mortality
names(v_r_mort_by_age) <- lt_usa_2015$Age

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
u_trtA  <- 0.95  # utility when being treated

## Transition rewards
du_HS1 <- 0.01  # disutility when transitioning from Healthy to Sick
ic_HS1 <- 1000  # increase in cost when transitioning from Healthy to Sick
ic_D   <- 2000  # increase in cost when dying

# Discount weight (equal discounting is assumed for costs and effects)
v_dwc <- 1 / ((1 + d_e) ^ (0:n_t))
v_dwe <- 1 / ((1 + d_c) ^ (0:n_t))

### Process model inputs
## Age-specific transition probabilities to the Dead state
# extract age-specific all-cause mortality for ages in model time horizon
v_r_HDage <- v_r_mort_by_age[(n_age_init + 1) + 0:(n_t - 1)]
# compute mortality rates
v_r_S1Dage <- v_r_HDage * hr_S1        # Age-specific mortality rate in the Sick state 
v_r_S2Dage <- v_r_HDage * hr_S2        # Age-specific mortality rate in the Sicker state 
# transform rates to probabilities
v_p_HDage  <- rate_to_prob(v_r_HDage)  # Age-specific mortality risk in the Healthy state 
v_p_S1Dage <- rate_to_prob(v_r_S1Dage) # Age-specific mortality risk in the Sick state
v_p_S2Dage <- rate_to_prob(v_r_S2Dage) # Age-specific mortality risk in the Sicker state

## Transition probability of becoming Sicker when Sick for treatment B
# transform odds ratios to probabilites 
lor_S1S2    <- log(or_S1S2)               # log-odds ratio of becoming Sicker when Sick
logit_S1S2  <- boot::logit(p_S1S2)        # log-odds of becoming Sicker when Sick
p_S1S2_trtB <- boot::inv.logit(logit_S1S2 +
                               lor_S1S2)  # probability to become Sicker when Sick 
# under B conditional on surviving

###################### Construct state-transition models #####################
#### Create transition arrays ####
# Initialize 3-D array
a_P <- array(0, dim      = c(n_states, n_states, n_t),
                dimnames = list(v_names_states, v_names_states, 0:(n_t - 1)))
### Fill in array
## From H
a_P["H", "H", ]   <- (1 - v_p_HDage) * (1 - p_HS1)
a_P["H", "S1", ]  <- (1 - v_p_HDage) * p_HS1
a_P["H", "D", ]   <- v_p_HDage
## From S1
a_P["S1", "H", ]  <- (1 - v_p_S1Dage) * p_S1H
a_P["S1", "S1", ] <- (1 - v_p_S1Dage) * (1 - (p_S1H + p_S1S2))
a_P["S1", "S2", ] <- (1 - v_p_S1Dage) * p_S1S2
a_P["S1", "D", ]  <- v_p_S1Dage
## From S2
a_P["S2", "S2", ] <- 1 - v_p_S2Dage
a_P["S2", "D", ]  <- v_p_S2Dage
## From D
a_P["D", "D", ]   <- 1

### For strategies B and AB
## Initialize transition probability array for strategies B and AB
a_P_strB <- a_P
## Only need to update the probabilities involving the transition from Sick to Sicker, p_S1S2
# From S1
a_P_strB["S1", "S1", ] <- (1 - v_p_S1Dage) * (1 - (p_S1H + p_S1S2_trtB))
a_P_strB["S1", "S2", ] <- (1 - v_p_S1Dage) * p_S1S2_trtB

### Check if transition probability matrix is valid (i.e., elements cannot < 0 or > 1) 
check_transition_probability(a_P,      verbose = TRUE)
check_transition_probability(a_P_strB, verbose = TRUE)
### Check if transition probability matrix sum to 1 (i.e., each row should sum to 1)
check_sum_of_transition_array(a_P,      n_states = n_states, n_t = n_t, verbose = TRUE)
check_sum_of_transition_array(a_P_strB, n_states = n_states, n_t = n_t, verbose = TRUE)

#### Run Markov model ####
## Initial state vector
# All starting healthy
v_s_init <- c(H = 1, S1 = 0, S2 = 0, D = 0) # initial state vector
v_s_init

## Initialize cohort trace for age-dependent (ad) cSTM for srategies SoC and A
m_M_ad <- matrix(0, 
                 nrow     = (n_t + 1), ncol = n_states, 
                 dimnames = list(0:n_t, v_names_states))
# Store the initial state vector in the first row of the cohort trace
m_M_ad[1, ] <- v_s_init
## Initialize cohort trace for srategies B and AB
m_M_ad_strB <- m_M_ad # structure and initial states remain the same.

## Initialize transition array which will capture transitions from each state to another over time 
# for srategies SoC and A
a_A <- array(0,
             dim      = c(n_states, n_states, n_t + 1),
             dimnames = list(v_names_states, v_names_states, 0:n_t))
# Set first slice of a_A with the initial state vector in its diagonal
diag(a_A[, , 1]) <- v_s_init
# For srategies B and AB, the array structure and initial state are identical 
a_A_strB <- a_A

## Iterative solution of age-dependent cSTM
for(t in 1:n_t){
  ## Fill in cohort trace
  # For srategies SoC and A
  m_M_ad[t + 1, ]      <- m_M_ad[t, ]      %*% a_P[, , t]
  # For srategies B and AB
  m_M_ad_strB[t + 1, ] <- m_M_ad_strB[t, ] %*% a_P_strB[, , t]
  
  ## Fill in transition-dynamics array
  # For srategies SoC and A
  a_A[, , t + 1]      <- m_M_ad[t, ]      * a_P[, , t]
  # For srategies B and AB
  a_A_strB[, , t + 1] <- m_M_ad_strB[t, ] * a_P_strB[, , t]
}

## Store the cohort traces in a list
l_m_M <- list(m_M_ad,
              m_M_ad,
              m_M_ad_strB,
              m_M_ad_strB)
names(l_m_M) <- v_names_str

#### Plot Outputs ####
## Plot the cohort traces 
plot_trace(m_M_ad)
plot_trace_strategy(l_m_M)
## Plot the epidemiology outcomes
survival_plot        <- plot_surv(l_m_M, v_names_death_states = "D") +
                        theme(legend.position = "")
prevalence_S1_plot   <- plot_prevalence(l_m_M, 
                                        v_names_sick_states = c("S1"), 
                                        v_names_dead_states = "D")  +
                        theme(legend.position = "")
prevalence_S2_plot   <- plot_prevalence(l_m_M, 
                                        v_names_sick_states = c("S2"), 
                                        v_names_dead_states = "D")  +
                        theme(legend.position = "")
prevalence_S1S2_plot <- plot_prevalence(l_m_M, 
                                        v_names_sick_states = c("S1", "S2"), 
                                        v_names_dead_states = "D") +
                        theme(legend.position = "")
prop_sicker_plot     <- plot_proportion_sicker(l_m_M, 
                                               v_names_sick_states = c("S1", "S2"), 
                                               v_names_sicker_states = c("S2")) +
                        theme(legend.position = "bottom")
grid.arrange(survival_plot, 
             prevalence_S1_plot, 
             prevalence_S2_plot, 
             prevalence_S1S2_plot, 
             prop_sicker_plot, 
             ncol = 1, heights = c(0.75, 0.75, 0.75, 0.75, 1))

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
  
  #### Expected QALYs and Costs per cycle ####
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
plot(df_cea, label = "all", txtsize = 16) +
     expand_limits(x = max(table_cea$QALYs) + 0.1) +
  theme(legend.position = c(0.8, 0.2))

###################### Probabalistic Sensitivty Analysis #####################
# Source functions that contain the model and CEA output
source('R/Functions STM_02.R')

# Function to generate PSA input dataset
generate_psa_params <- function(n_sim = 1000, seed = 071818){
  set.seed(seed) # set a seed to be able to reproduce the same results
  df_psa <- data.frame(
    # Transition probabilities (per cycle)
    p_HS1    = rbeta(n_sim, 30, 170),        # probability to become sick when healthy conditional on surviving
    p_S1H    = rbeta(n_sim, 60, 60) ,        # probability to become healthy when sick conditional on surviving
    p_S1S2   = rbeta(n_sim, 84, 716),        # probability to become Sicker when Sick conditional on surviving
    hr_S1    = rlnorm(n_sim, log(3), 0.01),  # rate ratio of death in S1 vs healthy
    hr_S2    = rlnorm(n_sim, log(10), 0.02), # rate ratio of death in S2 vs healthy 
    lor_S1S2 = rnorm(n_sim, log(0.6), 0.1),  # log-odds ratio of becoming Sicker when Sick under B
    
    # State rewards
    # Costs
    c_H    = rgamma(n_sim, shape = 100,   scale = 20),   # cost of remaining one cycle in state H
    c_S1   = rgamma(n_sim, shape = 177.8, scale = 22.5), # cost of remaining one cycle in state S1
    c_S2   = rgamma(n_sim, shape = 225,   scale = 66.7), # cost of remaining one cycle in state S2
    c_trtA = rgamma(n_sim, shape = 576,   scale = 20.8), # cost of treatment A (per cycle)
    c_trtB = rgamma(n_sim, shape = 676,   scale = 19.2), # cost of treatment B (per cycle)
    c_D    = 0,                                          # cost of being in the death state
    # Utilities
    u_H    = rbeta(n_sim, shape1 = 200, shape2 = 3),    # utility when healthy
    u_S1   = rbeta(n_sim, shape1 = 130, shape2 = 45),   # utility when sick
    u_S2   = rbeta(n_sim, shape1 = 230,  shape2 = 230), # utility when sicker
    u_D    = 0,                                         # utility when dead
    u_trtA = rbeta(n_sim, shape1 = 300, shape2 = 15),   # utility when being treated
    
    # Transition rewards
    du_HS1 = rbeta(n_sim, shape1 = 11,  shape2 = 1088), # disutility when transitioning from Healthy to Sick
    ic_HS1 = rgamma(n_sim, shape = 25,  scale = 40),    # increase in cost when transitioning from Healthy to Sick
    ic_D   = rgamma(n_sim, shape = 100, scale = 20)     # increase in cost when dying
  )
  return(df_psa)
}

# Number of PSA samples
n_sim <- 1000

# Generate PSA input dataset
df_psa_input <- generate_psa_params(n_sim = n_sim)
# First six observations
head(df_psa_input)

# Histogram of parameters
ggplot(melt(df_psa_input, variable.name = "Parameter"), 
       aes(x = value)) +
  facet_wrap(~Parameter, scales = "free") +
  geom_histogram(aes(y = ..density..)) +
  theme_bw(base_size = 16)

# Initialize matrices with PSA output 
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
  df_c[i, ]  <- l_out_temp$Cost  
  df_e[i, ]  <- l_out_temp$Effect
  # Display simulation progress
  if(i/(n_sim/100) == round(i/(n_sim/100), 0)) { # display progress every 5%
    cat('\r', paste(i/n_sim * 100, "% done", sep = " "))
  }
}

# Create PSA object for dampack
l_psa <- make_psa_obj(cost          = df_c, 
                      effectiveness = df_e, 
                      parameters    = df_psa_input, 
                      strategies    = v_names_str)
l_psa$strategies <- v_names_str
colnames(l_psa$effectiveness)<- v_names_str
colnames(l_psa$cost)<- v_names_str

# Vector with willingness-to-pay (WTP) thresholds.
v_wtp <- seq(0, 200000, by = 5000)

# Cost-Effectiveness Scatter plot
plot(l_psa)

## Conduct CEA with probabilistic output
# Compute expected costs and effects for each strategy from the PSA
df_out_ce_psa <- summary(l_psa)

# Calculate incremental cost-effectiveness ratios (ICERs)
df_cea_psa <- calculate_icers(cost       = df_out_ce_psa$meanCost, 
                              effect     = df_out_ce_psa$meanEffect,
                              strategies = df_out_ce_psa$Strategy)
df_cea_psa

## Plot cost-effectiveness frontier
plot(df_cea_psa, label = "all") +
  expand_limits(x = 20.8)

## Cost-effectiveness acceptability curves (CEACs) and frontier (CEAF)
ceac_obj <- ceac(wtp = v_wtp, psa = l_psa)
# Regions of highest probability of cost-effectiveness for each strategy
summary(ceac_obj)
# ceac_obj$Strategy <- ordered(ceac_obj$Strategy, levels = v_names_str)

# CEAC & CEAF plot
plot(ceac_obj, txtsize = 16, xlim = c(0, NA)) +
  theme(legend.position = "bottom")

##  Expected Loss Curves (ELCs)
elc_obj <- calc_exp_loss(wtp = v_wtp, psa = l_psa)
elc_obj

# ELC plot
plot(elc_obj, log_y = FALSE, txtsize = 16, xlim = c(0, NA)) +
  theme(legend.position = "bottom")

## Expected value of perfect information (EVPI)
evpi <- calc_evpi(wtp = v_wtp, psa = l_psa)
# EVPI plot
plot(evpi, effect_units = "QALY")
