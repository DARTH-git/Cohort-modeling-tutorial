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
# This model incorporates time-dependent and history- / state-dependent transition probabilities 

################################ Initial setup ############################### 
rm(list = ls())    # remove any variables in R's memory 

### Install packages
# install.packages("dplyr")    # to manipulate data
# install.packages("reshape2") # to transform data
# install.packages("ggplot2")  # to visualize data
# install.packages("scales")   # for dollar signs and commas
# install.packages("boot")     # to handle log odds and log odds ratios
# install.packages("devtools") # to ensure compatibility among packages
# devtools::install_github("DARTH-git/dampack") # to install dampack form GitHub, for CEA and calculate ICERs

### Load packages
library(dplyr)    
library(reshape2) 
library(ggplot2)   
library(scales)    
library(boot)
library(dampack) 

### Load functions
source("R/Functions.R")

################################ Model input ################################# 
## General setup
n_age_init <- 25                      # age at baseline
n_age_max  <- 110                     # maximum age of follow up
n_t        <- n_age_max - n_age_init  # time horizon, number of cycles
v_n        <- c("H", "S1", "S2", "D") # the 4 health states of the model: # change to v_names_states

## Tunnel inputs
# Number of tunnels
n_tunnel_size <- n_t 
# Name for tunnels states of Sick state
v_Sick_tunnel <- paste("S1_", seq(1, n_tunnel_size), "Yr", sep = "")
# Create variables for model with tunnels
v_n_tunnels <- c("H", v_Sick_tunnel, "S2", "D") # health state names
n_states_tunnels <- length(v_n_tunnels)         # number of health states

v_hcc    <- rep(1, n_t + 1)          # vector of half-cycle correction 
v_hcc[1] <- v_hcc[n_t + 1] <- 0.5
d_c <- 0.03                          # discount rate for costs 
d_e <- 0.03                          # discount rate for QALYs
v_names_str <- c("Usual care",       # store the strategy names
                 "New treatment 1", 
                 "New treatment 2",
                 "New treatments 1 & 2") 
n_str <- length(v_names_str) # number of strategies

## Transition probabilities (per cycle) and hazard ratios
p_HS1 <- 0.15  # probability to become Sick when Healthy conditional on surviving
p_S1H <- 0.50  # probability to become Healthy when Sick conditional on surviving
hr_S1 <- 3     # hazard ratio of death in Sick vs Healthy
hr_S2 <- 10    # hazard ratio of death in Sicker vs Healthy 

## History-dependent transition probability from S1 to S2 conditional on surviving
# Weibull parameters
n_lambda <- 0.08 # scale
n_gamma  <- 1.1  # shape
# Weibull hazard
v_p_S1S2_tunnels <- n_lambda * n_gamma * (1:n_tunnel_size)^{n_gamma-1}
# For new treatment 2
or_S1S2  <- 0.6                       # odds ratio of becoming Sicker when Sick under New treatment 2
lor_S1S2 <- log(or_S1S2)              # log-odd ratio of becoming Sicker when Sick
logit_S1S2 <- logit(v_p_S1S2_tunnels) # log-odds of becoming Sicker when Sick
# probability to become Sicker when Sick under New treatment 2 conditional on surviving
v_p_S1S2_tunnels_trt2 <- inv.logit(logit_S1S2 + lor_S1S2) 

## Age-dependent mortality rates
lt_usa_2015 <- read.csv("data/LifeTable_USA_Mx_2015.csv")
v_r_mort_by_age <- lt_usa_2015 %>% 
                   select(Total) %>%
                   as.matrix()

## Age-specific transition probabilities
# extract age-specific all-cause mortality for ages in model time horizon
v_r_HDage  <- v_r_mort_by_age[(n_age_init + 1) + 0:(n_t - 1)]
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
u_trt1 <- 0.95  # utility when being treated

## Transition rewards
du_HS1 <- 0.01  # disutility when transitioning from Healthy to Sick
ic_HS1 <- 1000  # increase in cost when transitioning from Healthy to Sick
ic_D   <- 2000  # increase in cost when dying

# Discount weight (equal discounting is assumed for costs and effects)
v_dwc <- 1 / ((1 + d_e) ^ (0:(n_t)))
v_dwe <- 1 / ((1 + d_c) ^ (0:(n_t)))

###################### Construct state-transition models #####################
#### Create transition matrix ####
# Initialize 3-D array
a_P_tunnels <- array(0, dim   = c(n_states_tunnels, n_states_tunnels, n_t),
                     dimnames = list(v_n_tunnels, v_n_tunnels, 0:(n_t - 1)))
### Fill in array
## From H
a_P_tunnels["H", "H", ]              <- (1 - v_p_HDage) * (1 - p_HS1)
a_P_tunnels["H", v_Sick_tunnel[1], ] <- (1 - v_p_HDage) * p_HS1
a_P_tunnels["H", "D", ]              <- v_p_HDage
## From S1
for(i in 1:(n_tunnel_size - 1)){
  a_P_tunnels[v_Sick_tunnel[i], "H", ]  <- (1 - v_p_S1Dage) * p_S1H
  a_P_tunnels[v_Sick_tunnel[i], 
              v_Sick_tunnel[i + 1], ]   <- (1 - v_p_S1Dage) * 
                                           (1 - (p_S1H + v_p_S1S2_tunnels[i]))
  a_P_tunnels[v_Sick_tunnel[i], "S2", ] <- (1 - v_p_S1Dage) * v_p_S1S2_tunnels[i]
  a_P_tunnels[v_Sick_tunnel[i], "D", ]  <- v_p_S1Dage
}
# repeat code for the last cycle to force the cohort stay in the last tunnel state of Sick
a_P_tunnels[v_Sick_tunnel[n_tunnel_size], "H", ]  <- (1 - v_p_S1Dage) * p_S1H
a_P_tunnels[v_Sick_tunnel[n_tunnel_size],
            v_Sick_tunnel[n_tunnel_size], ] <- (1 - v_p_S1Dage) *
                                               (1 - (p_S1H + v_p_S1S2_tunnels[n_tunnel_size]))
a_P_tunnels[v_Sick_tunnel[n_tunnel_size], "S2", ] <- (1 - v_p_S1Dage) * 
                                                      v_p_S1S2_tunnels[n_tunnel_size]
a_P_tunnels[v_Sick_tunnel[n_tunnel_size], "D", ]  <- v_p_S1Dage
## From S2
a_P_tunnels["S2", "S2", ] <- 1 - v_p_S2Dage
a_P_tunnels["S2", "D", ]  <- v_p_S2Dage
## From D
a_P_tunnels["D", "D", ] <- 1

# For New treatment 2
# Only need to update the transition probabilities from S1 involving p_S1S2
## From S1
a_P_tunnels_trt2 <- a_P_tunnels
for(i in 1:(n_tunnel_size - 1)){
  a_P_tunnels_trt2[v_Sick_tunnel[i], "H", ]  <- (1 - v_p_S1Dage) * p_S1H
  a_P_tunnels_trt2[v_Sick_tunnel[i], 
              v_Sick_tunnel[i + 1], ]        <- (1 - v_p_S1Dage) * 
                                                (1 - (p_S1H + v_p_S1S2_tunnels_trt2[i]))
  a_P_tunnels_trt2[v_Sick_tunnel[i], "S2", ] <- (1 - v_p_S1Dage) * v_p_S1S2_tunnels_trt2[i]
  a_P_tunnels_trt2[v_Sick_tunnel[i], "D", ]  <- v_p_S1Dage
}
# repeat code for the last cycle to force the cohort stay in the last tunnel state of Sick
a_P_tunnels_trt2[v_Sick_tunnel[n_tunnel_size], "H", ] <- (1 - v_p_S1Dage) * p_S1H
a_P_tunnels_trt2[v_Sick_tunnel[n_tunnel_size],
            v_Sick_tunnel[n_tunnel_size], ] <- (1 - v_p_S1Dage) * 
                                               (1 - (p_S1H +v_p_S1S2_tunnels_trt2[n_tunnel_size]))
a_P_tunnels_trt2[v_Sick_tunnel[n_tunnel_size], "S2", ] <- (1 - v_p_S1Dage) *
                                                           v_p_S1S2_tunnels_trt2[n_tunnel_size]
a_P_tunnels_trt2[v_Sick_tunnel[n_tunnel_size], "D", ]  <- v_p_S1Dage

### Check if transition probability matrix is valid (i.e., elements cannot < 0 or > 1) 
check_transition_probability(a_P_tunnels,      verbose = TRUE)
check_transition_probability(a_P_tunnels_trt2, verbose = TRUE)
### Check if transition probability matrix sum to 1 (i.e., each row should sum to 1)
check_sum_of_transition_array(a_P_tunnels,      n_states = n_states_tunnels, n_t = n_t, verbose = TRUE)
check_sum_of_transition_array(a_P_tunnels_trt2, n_states = n_states_tunnels, n_t = n_t, verbose = TRUE)

#### Run Markov model ####
## Initial state vector
# All starting healthy
v_s_init_tunnels <- c(1, rep(0, n_tunnel_size), 0, 0) 

## Initialize cohort trace for history-dependent cSTM
m_M_tunnels <- matrix(0, 
                      nrow = (n_t + 1), ncol = n_states_tunnels, 
                      dimnames = list(0:n_t, v_n_tunnels))
# Store the initial state vector in the first row of the cohort trace
m_M_tunnels[1, ] <- v_s_init_tunnels
# For New treatment 2
m_M_tunnels_trt2 <- m_M_tunnels

## Initialize transition array
a_A_tunnels <- array(0,
                     dim = c(n_states_tunnels, n_states_tunnels, n_t + 1),
                     dimnames = list(v_n_tunnels, v_n_tunnels, 0:n_t))
# Set first slice of A with the initial state vector in its diagonal
diag(a_A_tunnels[, , 1]) <- v_s_init_tunnels
# For New treatment 2
a_A_tunnels_trt2 <- a_A_tunnels

## Iterative solution of age-dependent cSTM
for(t in 1:n_t){
  # Fill in cohort trace
  m_M_tunnels[t + 1, ]     <- m_M_tunnels[t, ]      %*% a_P_tunnels[, , t]
  m_M_tunnels_trt2[t + 1,] <- m_M_tunnels_trt2[t, ] %*% a_P_tunnels_trt2[, , t]
  # Fill in transition dynamics array
  a_A_tunnels[, , t + 1]       <- m_M_tunnels[t, ]      * a_P_tunnels[, , t]
  a_A_tunnels_trt2[, , t + 1]  <- m_M_tunnels_trt2[t, ] * a_P_tunnels_trt2[, , t]
}
# Create aggregated trace
m_M_tunnels_sum <- cbind(H  = m_M_tunnels[, "H"], 
                         S1 = rowSums(m_M_tunnels[, 2:(n_tunnel_size +1)]), 
                         S2 = m_M_tunnels[, "S2"],
                         D  = m_M_tunnels[, "D"])
m_M_tunnels_sum_trt2 <- cbind(H  = m_M_tunnels_trt2[, "H"], 
                              S1 = rowSums(m_M_tunnels_trt2[, 2:(n_tunnel_size +1)]), 
                              S2 = m_M_tunnels_trt2[, "S2"],
                              D  = m_M_tunnels_trt2[, "D"])

## Store the cohort traces in a list
l_m_M <- list(m_M_tunnels_sum,      # Cohort trace for Usual Care
              m_M_tunnels_sum,      # Cohort trace for New Treatment 1, same as Usual Care
              m_M_tunnels_sum_trt2, # Cohort trace for New Treatment 2
              m_M_tunnels_sum_trt2) # Cohort trace for New Treatment 1 & 2, same as New Treatment 1
names(l_m_M) <- v_names_str

#### Plot Outputs ####
## Plot the cohort trace for Usual care and New Treatment 1
plot_trace(m_M_tunnels_sum)
plot_trace_strategy(l_m_M)

#### State Rewards ####
## Vector of utilities for S1 under Usual care
v_u_S1_UC <- rep(u_S1, n_tunnel_size)
names(v_u_S1_UC) <- v_Sick_tunnel
## Vector of state utilities under Usual care
v_u_UC <- c(H  = u_H, 
            v_u_S1_UC, 
            S2 = u_S2,
            D  = u_D)
## Vector of costs per cycle for S1 under Usual care
v_c_S1_UC <- rep(c_S1, n_tunnel_size)
names(v_c_S1_UC) <- v_Sick_tunnel
## Vector of state costs per cycle under Usual care
v_c_UC <- c(H  = c_H,
            v_c_S1_UC, 
            S2 = c_S2,
            D  = c_D)

## Vector of utilities for S1 under New treatment 1
v_u_S1_trt1 <- rep(u_trt1, n_tunnel_size)
names(v_u_S1_trt1) <- v_Sick_tunnel
## Vector of state utilities under New treatment 1
v_u_trt1 <- c(H  = u_H, 
              v_u_S1_trt1, 
              S2 = u_S2, 
              D  = u_D)
## Vector of costs per cycle for S1 under New treatment 1
v_c_S1_trt1 <- rep(c_S1 + c_trt1, n_tunnel_size)
names(v_c_S1_trt1) <- v_Sick_tunnel
## Vector of state costs per cycle under New treatment 1
v_c_trt1 <- c(H  = c_H, 
              v_c_S1_trt1,
              S2 = c_S2 + c_trt1,
              D  = c_D)

## Vector of utilities for S1 under New treatment 2
v_u_S1_trt2 <- rep(u_S1, n_tunnel_size)
names(v_u_S1_trt2) <- v_Sick_tunnel
## Vector of state utilities under New treatment 2
v_u_trt2 <- c(H  = u_H, 
              v_u_S1_trt2, 
              S2 = u_S2,
              D  = u_D)
## Vector of costs per cycle for S1 under New treatment 2
v_c_S1_trt2 <- rep(c_S1 + c_trt2, n_tunnel_size)
names(v_c_S1_trt2) <- v_Sick_tunnel
## Vector of state costs per cycle under New treatment 2
v_c_trt2 <- c(H  = c_H, 
              v_c_S1_trt2, 
              S2 = c_S2 + c_trt2,
              D  = c_D)

## Vector of utilities for S1 under New treatments 1 & 2
v_u_S1_trt1_2 <- rep(u_trt1, n_tunnel_size)
names(v_u_S1_trt1_2) <- v_Sick_tunnel
## Vector of state utilities under New treatment 2
v_u_trt1_2 <- c(H  = u_H, 
                v_u_S1_trt1_2, 
                S2 = u_S2, 
                D  = u_D)
## Vector of costs per cycle for S1 under New treatment 2
v_c_S1_trt1_2 <- rep(c_S1 + (c_trt1 + c_trt2), n_tunnel_size)
names(v_c_S1_trt1_2) <- v_Sick_tunnel
## Vector of state costs per cycle under New treatment 2
v_c_trt1_2 <- c(H  = c_H, 
                v_c_S1_trt1_2, 
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
l_a_A <- list(a_A_tunnels,
              a_A_tunnels,
              a_A_tunnels_trt2,
              a_A_tunnels_trt2)

# assign strategy names to matching items in the lists
names(l_u) <- names(l_c) <- names(l_a_A) <- v_names_str

## create empty vectors to store total utilities and costs 
v_totqaly <- v_totcost <- vector(mode = "numeric", length = n_str)
names(v_totqaly) <- names(v_totcost) <- v_names_str

#### Loop through each strategy and calculate total utilities and costs ####
for (i in 1:n_str) {
  v_u <- l_u[[i]]    # select the vector of state utilities for the ith strategy
  v_c <- l_c[[i]]    # select the vector of state costs for the ith strategy
  a_A_tunnels <- l_a_A[[i]]  # select the transition array for the ith strategy
  
  #### Array of state utilities and costs ####
  # Create transition matrices of state utilities and state costs for the ith strategy 
  m_u   <- matrix(v_u, nrow = n_states_tunnels, ncol = n_states_tunnels, byrow = T)
  m_c   <- matrix(v_c, nrow = n_states_tunnels, ncol = n_states_tunnels, byrow = T)
  # Expand the transition matrix of state utilities across cycles to form a transition array of state utilities
  a_R_u <- array(m_u, 
                 dim      = c(n_states_tunnels, n_states_tunnels, n_t + 1),
                 dimnames = list(v_n_tunnels, v_n_tunnels, 0:n_t))
  # Expand the transition matrix of state costs across cycles to form a transition array of state costs
  a_R_c <- array(m_c, 
                 dim      = c(n_states_tunnels, n_states_tunnels, n_t + 1),
                 dimnames = list(v_n_tunnels, v_n_tunnels, 0:n_t))
  
  #### Apply transition rewards ####  
  # Add disutility due to transition from H to S1
  a_R_u["H", "S1_1Yr", ]          <- a_R_u["H", "S1_1Yr", ]          - du_HS1
  # Add transition cost per cycle due to transition from H to S1
  a_R_c["H", "S1_1Yr", ]          <- a_R_c["H", "S1_1Yr", ]          + ic_HS1
  # Add transition cost per cycle of dying from all non-dead states
  a_R_c[-n_states_tunnels, "D", ] <- a_R_c[-n_states_tunnels, "D", ] + ic_D
  
  #### Expected QALYs and Costs for all transitions per cycle ####
  # QALYs = life years x QoL
  # Note: all parameters are annual in our example. In case your own case example is different make sure you correctly apply .
  a_Y_c <- a_A_tunnels * a_R_c
  a_Y_u <- a_A_tunnels * a_R_u 
  
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
     expand_limits(x = max(table_cea$QALYs) + 0.3, 
                   y = 250000) 

###################### Probabalistic Sensitivty Analysis #####################
# Source functions that contain the model and CEA output
source('R/Functions STM_03.R')

# Function to generate PSA input dataset
generate_psa_params <- function(n_sim = 1000, seed = 071818){
  set.seed(seed) # set a seed to be able to reproduce the same results
  df_psa <- data.frame(
    # Transition probabilities (per cycle)
    p_HS1    = rbeta(n_sim, 30, 170),          # probability to become sick when healthy
    p_S1H    = rbeta(n_sim, 312, 312) ,        # probability to become healthy when sick
    hr_S1    = rlnorm(n_sim, log(3),  1),      # rate ratio of death in S1 vs healthy
    hr_S2    = rlnorm(n_sim, log(10), 1),      # rate ratio of death in S2 vs healthy 
    n_lambda = rlnorm(n_sim, log(0.08), 0.02), # transition from S1 to S2 - Weibull scale parameter
    n_gamma  = rlnorm(n_sim, log(1.1), 0.05),  # transition from S1 to S2 - Weibull shape parameter
    lor_S1S2 = rnorm(n_sim, log(0.6), 0.1),    # log-odds ratio of becoming Sicker whe 
    # State rewards
    # Costs
    c_H    = rgamma(n_sim, shape = 100,   scale = 20),   # cost of remaining one cycle in state H
    c_S1   = rgamma(n_sim, shape = 177.8, scale = 22.5), # cost of remaining one cycle in state S1
    c_S2   = rgamma(n_sim, shape = 225,   scale = 66.7), # cost of remaining one cycle in state S2
    c_trt1 = rgamma(n_sim, shape = 576,   scale = 20.8), # cost of treatment (per cycle)
    c_trt2 = rgamma(n_sim, shape = 676,   scale = 19.2), # cost of treatment (per cycle)
    c_D    = 0,                                          # cost of being in the death state
    
    # Utilities
    u_H    = rbeta(n_sim, shape1 = 200, shape2 = 3),  # utility when healthy
    u_S1   = rbeta(n_sim, shape1 = 130, shape2 = 45), # utility when sick
    u_S2   = rbeta(n_sim, shape1 = 50,  shape2 = 50), # utility when sicker
    u_D    = 0,                                       # utility when dead
    u_trt1 = rbeta(n_sim, shape1 = 300, shape2 = 15), # utility when being treated
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
  df_c[i, ]  <- l_out_temp$Cost  # HUGE PROBLEM HERE
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
v_wtp <- seq(0, 200000, by = 10000)

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

