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

## Transition probabilities (per cycle) and hazard ratios
p_HS1       <- 0.15                   # probability to become Sick when Healthy conditional on surviving
p_S1H       <- 0.5                    # probability to become Healthy when Sick conditional on surviving
p_S1S2      <- 0.105                  # probability to become Sicker when Sick conditional on surviving
hr_S1       <- 3                      # hazard ratio of death in Sick vs Healthy 
hr_S2       <- 10                     # hazard ratio of death in Sicker vs Healthy 
# For New treatment 2
or_S1S2     <- 0.6                    # odds ratio of becoming Sicker when Sick under New treatment 2
lor_S1S2    <- log(or_S1S2)           # log-odds ratio of becoming Sicker when Sick
logit_S1S2  <- logit(p_S1S2)          # log-odds of becoming Sicker when Sick
p_S1S2_trt2 <- inv.logit(logit_S1S2 +
                           lor_S1S2)  # probability to become Sicker when Sick under New treatment 2

## Age-dependent mortality rates
lt_usa_2015 <- read.csv("data/LifeTable_USA_Mx_2015.csv")
v_r_mort_by_age <- lt_usa_2015 %>% 
                   select(Total) %>%
                   as.matrix() # anyone above 100 have the same mortality, ALTERNATIVELY: conditional probs

## Age-specific transition probabilities
# extract age-specific all-cause mortality for ages in model time horizon
v_r_HDage  <- v_r_mort_by_age[(n_age_init) + 1:n_t]
v_r_S1Dage <- v_r_HDage * hr_S1        # Age-specific mortality rate in the Sick state 
v_r_S2Dage <- v_r_HDage * hr_S2        # Age-specific mortality rate in the Sicker state 
v_p_HDage  <- rate_to_prob(v_r_HDage)  # Age-specific mortality risk in the Healthy state 
v_p_S1Dage <- rate_to_prob(v_r_S1Dage) # Age-specific mortality risk in the Sick state
v_p_S2Dage <- rate_to_prob(v_r_S2Dage) # Age-specific mortality risk in the Sicker state

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
u_trt1 <- 0.95  # utility when being treate

## Transition rewards
du_HS1 <- 0.01  # disutility when transitioning from Healthy to Sick
ic_HS1 <- 1000  # increase in cost when transitioning from Healthy to Sick
ic_D   <- 2000  # increase in cost when dying

# Discount weight (equal discounting is assumed for costs and effects)
v_dwc <- 1 / ((1 + d_e) ^ (0:(n_t)))
v_dwe <- 1 / ((1 + d_c) ^ (0:(n_t)))

###################### Construct state-transition models #####################
#### Create transition arrays ####
# Initialize 3-D array
a_P <- array(0, dim      = c(n_states, n_states, n_t),
                dimnames = list(v_n, v_n, 0:(n_t - 1)))
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

### For New treatment 2
## Initialize transition probability array for new treatment 2
a_P_trt2 <- a_P
## Only need to update the probabilities involving the transition from Sick to Sicker, p_S1S2
# From S1
a_P_trt2["S1", "S1", ] <- (1 - v_p_S1Dage) * (1 - (p_S1H + p_S1S2_trt2))
a_P_trt2["S1", "S2", ] <- (1 - v_p_S1Dage) * p_S1S2_trt2

### Check if transition probability matrix is valid (i.e., elements cannot < 0 or > 1) 
check_transition_probability(a_P,      verbose = TRUE)
check_transition_probability(a_P_trt2, verbose = TRUE)
### Check if transition probability matrix sum to 1 (i.e., each row should sum to 1)
check_sum_of_transition_array(a_P,      n_states = n_states, n_t = n_t, verbose = TRUE)
check_sum_of_transition_array(a_P_trt2, n_states = n_states, n_t = n_t, verbose = TRUE)

#### Run Markov model ####
## Initial state vector
# All starting healthy
v_s_init <- c(H = 1, S1 = 0, S2 = 0, D = 0) # initial state vector
v_s_init

## Initialize cohort trace for age-dependent (ad) cSTM for usual care and new treatment 1
m_M_ad <- matrix(0, 
                 nrow     = (n_t + 1), ncol = n_states, 
                 dimnames = list(0:n_t, v_n))
# Store the initial state vector in the first row of the cohort trace
m_M_ad[1, ] <- v_s_init
## Initialize cohort trace new treatment 1 and combination of both new treatments
m_M_ad_trt2 <- m_M_ad # structure and initial states remain the same.

## Initialize transition array which will capture transitions from each state to another over time 
a_A <- array(0,
             dim      = c(n_states, n_states, n_t + 1),
             dimnames = list(v_n, v_n, 0:n_t))

# Set first slice of A with the initial state vector in its diagonal
diag(a_A[, , 1]) <- v_s_init
# For New treatment 2
a_A_trt2 <- a_A

## Iterative solution of age-dependent cSTM
for(t in 1:n_t){
  ## Fill in cohort trace
  # For usual care and new treatment 1
  m_M_ad[t + 1, ]      <- m_M_ad[t, ]      %*% a_P[, , t]
  # For new treatment 2 and combination of both new treatments
  m_M_ad_trt2[t + 1, ] <- m_M_ad_trt2[t, ] %*% a_P_trt2[, , t]
  
  ## Fill in transition dynamics array
  # For usual care and new treatment 1
  a_A[, , t + 1]       <- m_M_ad[t, ]       * a_P[, , t]
  # For new treatment 2 and combination of both new treatments
  a_A_trt2[, , t + 1]  <- m_M_ad_trt2[t, ]  * a_P_trt2[, , t]
}

## Store the cohort traces in a list
l_m_M <- list(m_M_ad,
              m_M_ad,
              m_M_ad_trt2,
              m_M_ad_trt2)
names(l_m_M) <- v_names_str

#### Plot Outputs ####
## Plot the cohort traces 
plot_trace(m_M_ad)
plot_trace_strategy(l_m_M)
## Plot the epidemiology outcomes
survival_plot        <- plot_surv(l_m_M, v_names_death_states = "D") +
                        theme(legend.position = "")
prevalence_S1_plot   <- plot_prevalence(l_m_M, v_names_sick_states = c("S1"), v_names_dead_states = "D")  +
                        theme(legend.position = "")
prevalence_S2_plot   <- plot_prevalence(l_m_M, v_names_sick_states = c("S2"), v_names_dead_states = "D")  +
                        theme(legend.position = "")
prevalence_S1S2_plot <- plot_prevalence(l_m_M, v_names_sick_states = c("S1", "S2"), v_names_dead_states = "D") +
                        theme(legend.position = "")
prop_sicker_plot     <- plot_proportion_sicker(l_m_M, v_names_sick_states = c("S1", "S2"), v_names_sicker_states = c("S2")) +
                        theme(legend.position = "bottom")
grid.arrange(survival_plot, prevalence_S1_plot, prevalence_S2_plot, prevalence_S1S2_plot, prop_sicker_plot, 
             ncol = 1, heights = c(0.75, 0.75, 0.75, 0.75, 1))

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
## Vector of state costs per cycle under new treatment 1
v_c_trt1   <- c(H  = c_H, 
                S1 = c_S1 + c_trt1,
                S2 = c_S2 + c_trt1, 
                D  = c_D)
## Vector of state utilities under new treatment 2
v_u_trt2   <- c(H  = u_H, 
                S1 = u_S1, 
                S2 = u_S2, 
                D  = u_D)
## Vector of state costs per cycle under new treatment 2
v_c_trt2   <- c(H  = c_H, 
                S1 = c_S1 + c_trt2, 
                S2 = c_S2 + c_trt2, 
                D  = c_D)
## Vector of state utilities under new treatments 1 & 2
v_u_trt1_2 <- c(H  = u_H, 
                S1 = u_trt1, 
                S2 = u_S2, 
                D  = u_D)
## Vector of state costs per cycle under new treatments 1 & 2
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
l_a_A <- list(a_A,
              a_A,
              a_A_trt2,
              a_A_trt2)

# assign strategy names to matching items in the lists
names(l_u) <- names(l_c) <- names(l_a_A) <- v_names_str

## create empty vectors to store total utilities and costs 
v_totqaly <- v_totcost <- vector(mode = "numeric", length = n_str)
names(v_totqaly) <- names(v_totcost) <- v_names_str

#### Loop through each strategy and calculate total utilities and costs ####
for (i in 1:n_str) {
  v_u <- l_u[[i]]    # select the vector of state utilities for the ith strategy
  v_c <- l_c[[i]]    # select the vector of state costs for the ith strategy
  a_A <- l_a_A[[i]]  # select the transition array for the ith strategy
  
  #### Array of state rewards ####
  # Create transition matrices of state utilities and state costs for the ith strategy 
  m_u   <- matrix(v_u, nrow = n_states, ncol = n_states, byrow = T)
  m_c   <- matrix(v_c, nrow = n_states, ncol = n_states, byrow = T)
  # Expand the transition matrix of state utilities across cycles to form a transition array of state utilities
  a_R_u <- array(m_u, 
                 dim      = c(n_states, n_states, n_t + 1),
                 dimnames = list(v_n, v_n, 0:n_t))
  # Expand the transition matrix of state costs across cycles to form a transition array of state costs
  a_R_c <- array(m_c, 
                 dim      = c(n_states, n_states, n_t + 1),
                 dimnames = list(v_n, v_n, 0:n_t))
  
  #### Apply transition rewards ####  
  # Apply disutility due to transition from H to S1
  a_R_u["H", "S1", ]      <- a_R_u["H", "S1", ]       - du_HS1
  # Add transition cost per cycle due to transition from H to S1
  a_R_c["H", "S1", ]      <- a_R_c["H", "S1", ]       + ic_HS1
  # Add transition cost  per cycle of dying from all non-dead states
  a_R_c[-n_states, "D", ] <- a_R_c[- n_states, "D", ] + ic_D
  
  #### Expected QALYs and Costs for all transitions per cycle ####
  # QALYs = life years x QoL
  # Note: all parameters are annual in our example. In case your own case example is different make sure you correctly apply .
  a_Y_c <- a_A * a_R_c
  a_Y_u <- a_A * a_R_u 
  
  #### Expected QALYs and Costs per cycle ####
  ## Vector of QALYs under Usual Care
  v_qaly <- apply(a_Y_u, 3, sum) # sum the proportion of the cohort across transitions 
  v_cost <- apply(a_Y_c, 3, sum) # sum the proportion of the cohort across transitions
  
  #### Discounted total expected QALYs and Costs per strategy and apply half-cycle correction if applicable ####
  ## QALYs
  v_totqaly[i] <- t(v_qaly) %*% (v_dwe * v_hcc)
  ## Costs
  v_totcost[i] <- t(v_cost) %*% (v_dwc * v_hcc)
}

########################## Cost-effectiveness analysis #######################
### Calculate incremental cost-effectiveness ratios (ICERs)
df_cea <- calculate_icers(cost       = v_totcost, 
                          effect     = v_totqaly,
                          strategies = v_names_str)
df_cea

### Create CEA table in proper format
table_cea <- format_table_cea(df_cea)
table_cea

### CEA frontier
plot(df_cea, label = "all") +
     expand_limits(x = max(table_cea$QALYs) + 0.5) 


