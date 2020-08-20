###############################################################################
### Cohort state-transition models in R:                                     ##
### From conceptualization to implementation 2019                            ##
###############################################################################
# This code forms the basis for the state-transition model of the article: 
# 'Cohort state-transition models in R: From conceptualization to implementation' 
# Authors: 
# - Fernando Alarid-Escudero <fernando.alarid@cide.edu>
# - Eline Krijkamp
# - Eva A. Enns
# - Alan Yang
# - Myriam G.M. Hunink
# - Petros Pechlivanoglou
# - Hawre Jalal
# Please cite the article when using this code
#
# To program this tutorial we made use of 
# R version 3.5.0 (2018-04-23)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.14.5
# RStudio: Version 1.1.453 2009-2018 RStudio, Inc

###############################################################################
################# Code of Appendix   ##########################################
###############################################################################
# Implements an age- and history-dependent Sick-Sicker cSTM model                          
# + include code for a probabilistic sensitivity analysis (PSA) 
##################################### Initial setup ###########################
# rm(list = ls())  # remove any variables in R's memory 
library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(scales)   # for dollar signs and commas
library(truncnorm)
library(pryr)
# devtools::install_github("DARTH-git/dampack") # to install dampack form GitHub
library(dampack)  # for CEA and calculate ICERs

################################## DARTH colors  ###############################

# code for the DARTH colors for the figures
DARTHgreen      <- '#009999'  
DARTHyellow     <- '#FDAD1E'  
DARTHblue       <- '#006699' 
DARTHlightgreen <- '#00adad'
DARTHgray       <- '#666666'

##################################### Model input ##############################
## General setup
n_age_init <- 25 # age at baseline
n_age_max <- 110 # maximum age of follow up
n_t <- n_age_max - n_age_init # time horizon, number of cycles


## Tunnel inputs
# Number of tunnels
n_tunnel_size <- n_t 
# Name for tunnels states of Sick state
v_Sick_tunnel <- paste("S1_", seq(1, n_tunnel_size), "Yr", sep = "")
# Create variables for model with tunnels
v_n_tunnels <- c("H", v_Sick_tunnel, "S2", "D") # health state names
n_states_tunnels <- length(v_n_tunnels)         # number of health states

d_c <- 0.03 # discount rate for costs 
d_e <- 0.03 # discount rate for QALYs
v_names_str <- c("Usual care", "New treatment") # store the strategy names
n_str <- length(v_names_str) # number of strategies
v_hcc    <- rep(1, n_t + 1)     # vector of half-cycle correction 
v_hcc[1] <- v_hcc[n_t + 1] <- 0.5


## Transition probabilities (per cycle) and hazard ratios
p_HS1   <- 0.15 # probability to become Sick when Healthy
p_S1H   <- 0.5  # probability to become Healthy when Sick
hr_S1   <- 3    # hazard ratio of death in Sick vs Healthy
hr_S2   <- 10   # hazard ratio of death in Sicker vs Healthy 

## History-dependent transition from S1 to S2
# Weibull parameters
n_lambda <- 0.08 # scale
n_gamma  <- 1.1  # shape
# Weibull function
p_S1S2_tunnels <- n_lambda * n_gamma * (1:n_tunnel_size)^{n_gamma-1}

## Age-dependent mortality rates
lt_usa_2005 <- read.csv("data/LifeTable_USA_Mx_2015.csv")
v_r_mort_by_age <- lt_usa_2005 %>% 
  # filter(Age >= age & Age <= n_age_max) %>%
  select(Total) %>%
  as.matrix()

## State rewards
# Costs
c_H   <- 2000  # cost of remaining one cycle Healthy 
c_S1  <- 4000  # cost of remaining one cycle Sick 
c_S2  <- 15000 # cost of remaining one cycle Sicker 
c_D   <- 0     # cost of being dead (per cycle)
c_Trt <- 12000 # cost of treatment (per cycle) 
# Utilities
u_H   <- 1     # utility when Healthy 
u_S1  <- 0.75  # utility when Sick 
u_S2  <- 0.5   # utility when Sicker
u_D   <- 0     # utility when Healthy 
u_Trt <- 0.95  # utility when being treated

## Transition rewards
du_HS1 <- 0.01  # disutility when transitioning from Healthy to Sick
ic_HS1 <- 1000  # increase in cost when transitioning from Healthy to Sick
ic_D   <- 2000  # increase in cost when dying

# Discount weight (equal discounting is assumed for costs and effects)
v_dwc <- 1 / ((1 + d_e) ^ (0:(n_t)))
v_dwe <- 1 / ((1 + d_c) ^ (0:(n_t)))

## Age-specific transition probabilities
# Age-specific probability of dying when Healthy (all-cause mortality)
v_p_HDage  <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_t - 1)]) 
# Age-specific mortality risk in the Sick state
v_p_S1Dage <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_t - 1)] * hr_S1)
# Age-specific mortality risk in the Sicker state
v_p_S2Dage <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_t - 1)] * hr_S2)

######################## Construct state-transition models ####################
#### Create transition matrix ####
# Initialize 3-D array
a_P_tunnels <- array(0, dim = c(n_states_tunnels, n_states_tunnels, n_t),
                     dimnames = list(v_n_tunnels, v_n_tunnels, 0:(n_t - 1)))

### Fill in array
## From H
a_P_tunnels["H", "H", ]              <- 1 - (p_HS1 + v_p_HDage)
a_P_tunnels["H", v_Sick_tunnel[1], ] <- p_HS1
a_P_tunnels["H", "D", ]              <- v_p_HDage
## From S1
for(i in 1:(n_tunnel_size - 1)){
  a_P_tunnels[v_Sick_tunnel[i], "H", ]  <- p_S1H
  a_P_tunnels[v_Sick_tunnel[i], 
              v_Sick_tunnel[i + 1], ]   <- 1 - (p_S1H + p_S1S2_tunnels[i] + v_p_S1Dage)
  a_P_tunnels[v_Sick_tunnel[i], "S2", ] <- p_S1S2_tunnels[i]
  a_P_tunnels[v_Sick_tunnel[i], "D", ]  <- v_p_S1Dage
}
# repeat code for the last cycle to force the cohort stay in the last tunnel stateof Sick
a_P_tunnels[v_Sick_tunnel[n_tunnel_size], "H", ]  <- p_S1H
a_P_tunnels[v_Sick_tunnel[n_tunnel_size],
            v_Sick_tunnel[n_tunnel_size], ] <- 1 - (p_S1H +
                                                      p_S1S2_tunnels[n_tunnel_size] + 
                                                      v_p_S1Dage)
a_P_tunnels[v_Sick_tunnel[n_tunnel_size], "S2", ] <- p_S1S2_tunnels[n_tunnel_size]
a_P_tunnels[v_Sick_tunnel[n_tunnel_size], "D", ]  <- v_p_S1Dage
## From S2
a_P_tunnels["S2", "S2", ] <- 1 - v_p_S2Dage
a_P_tunnels["S2", "D", ]  <- v_p_S2Dage
## From D
a_P_tunnels["D", "D", ] <- 1

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

## Initialize transition array
a_A_tunnels <- array(0,
                     dim = c(n_states_tunnels, n_states_tunnels, n_t + 1),
                     dimnames = list(v_n_tunnels, v_n_tunnels, 0:n_t))
# Set first slice of A with the initial state vector in its diagonal
diag(a_A_tunnels[, , 1]) <- v_s_init_tunnels

## Iterative solution of age-dependent cSTM
for(t in 1:n_t){
  # Fill in cohort trace
  m_M_tunnels[t + 1, ] <- m_M_tunnels[t, ] %*% a_P_tunnels[, , t]
  # Fill in transition dynamics array
  a_A_tunnels[, , t + 1]  <- m_M_tunnels[t, ] * a_P_tunnels[, , t]
}
# Create aggregated trace
m_M_tunnels_sum <- cbind(H = m_M_tunnels[, "H"], 
                         S1 = rowSums(m_M_tunnels[, 2:(n_tunnel_size +1)]), 
                         S2 = m_M_tunnels[, "S2"],
                         D = m_M_tunnels[, "D"])
#### Plot Outputs ####
### Cohort trace
## Define colors and line types
cols <- c("H" = DARTHgreen, "S1" = DARTHblue, 
          "S2" = DARTHyellow, "D" = DARTHgray)
lty <-  c("H" = 1, "S1" = 2, "S2" = 4, "D" = 3)
## Plot the cohort trace
ggplot(melt(m_M_tunnels_sum), aes(x = Var1, y = value, 
                         color = Var2, linetype = Var2)) +
  geom_line(size = 1) +
  scale_colour_manual(name = "Health state", 
                      values = cols) +
  scale_linetype_manual(name = "Health state",
                        values = lty) +
  xlab("Cycle") +
  ylab("Proportion of the cohort") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", 
        legend.background = element_rect(fill = NA))

#### State and Transition Rewards ####
### State rewards
## Vector of utilities for S1 under Usual Care
v_u_S1_UC <- rep(u_S1, n_tunnel_size)
names(v_u_S1_UC) <- v_Sick_tunnel
## Vector of state utilities under Usual Care
v_u_UC <- c(H = u_H, v_u_S1_UC, S2 = u_S2, D = u_D)

## Vector of costs for S1 under Usual Care
v_c_S1_UC <- rep(c_S1, n_tunnel_size)
names(v_c_S1_UC) <- v_Sick_tunnel
## Vector of state costs under Usual Care
v_c_UC <- c(H = c_H, v_c_S1_UC, S2 = c_S2, D = c_D)

## Vector of utilities for S1 under New Treatment
v_u_S1_Trt <- rep(u_Trt, n_tunnel_size)
names(v_u_S1_Trt) <- v_Sick_tunnel
## Vector of state utilities under New Treatment
v_u_Trt <- c(H = u_H, v_u_S1_Trt, S2 = u_S2, D = u_D)

## Vector of costs for S1 under Usual Care
v_c_S1_Trt <- rep(c_S1 + c_Trt, n_tunnel_size)
names(v_c_S1_Trt) <- v_Sick_tunnel
## Vector of state costs under New Treatment
v_c_Trt <- c(H = c_H, v_c_S1_Trt, S2 = c_S2 + c_Trt, D = c_D)

### Arrays of rewards
## Array of state and transition utilities under Usual Care
a_R_u_UC <- aperm(array(v_u_UC,
                        dim = c(n_states_tunnels, n_states_tunnels, n_t + 1),
                        dimnames = list(v_n_tunnels, v_n_tunnels, 0:n_t)),
                  perm = c(2, 1, 3))
## Array of state and transition costs under Usual Care
a_R_c_UC <- aperm(array(v_c_UC,
                        dim = c(n_states_tunnels, n_states_tunnels, n_t + 1),
                        dimnames = list(v_n_tunnels, v_n_tunnels, 0:n_t)),
                  perm = c(2, 1, 3))

## Array of utilities under New Treatment
a_R_u_Trt <- aperm(array(v_u_Trt,
                        dim = c(n_states_tunnels, n_states_tunnels, n_t + 1),
                        dimnames = list(v_n_tunnels, v_n_tunnels, 0:n_t)),
                  perm = c(2, 1, 3))
## Array of costs per cycle under New Treatment
a_R_c_Trt <- aperm(array(v_c_Trt,
                        dim = c(n_states_tunnels, n_states_tunnels, n_t + 1),
                        dimnames = list(v_n_tunnels, v_n_tunnels, 0:n_t)),
                  perm = c(2, 1, 3))

### Transition rewards
## For Usual Care
# Add disutility due to transition from H to S1
a_R_u_UC["H", "S1_1Yr", ] <- a_R_u_UC["H", "S1_1Yr", ] - du_HS1
# Add transition cost due to transition from H to S1
a_R_c_UC["H", "S1_1Yr", ] <- a_R_c_UC["H", "S1_1Yr", ] + ic_HS1
# Add transition cost of dying from all non-dead states
a_R_c_UC[-n_states_tunnels, "D", ] <- a_R_c_UC[-n_states_tunnels, "D", ] + ic_D

## For New Treatment
# Add disutility due to transition from Healthy to Sick
a_R_u_Trt["H", "S1_1Yr", ] <- a_R_u_Tr["H", "S1_1Yr", ] - du_HS1
# Add transition cost due to transition from Healthy to Sick
a_R_c_Trt["H", "S1_1Yr", ] <- a_R_c_Tr["H", "S1_1Yr", ] + ic_HS1
# Add transition cost of dying from all non-dead states
a_R_c_Trt[-n_states_tunnels, "D", ] <- a_R_c_Trt[-n_states_tunnels, "D", ] + ic_D

#### Expected QALYs and Costs for all transitions per cycle ####
# QALYs = life years x QoL
# life years are markov trace * cycle length  
n_cpy <- 1   # the number of cycles per year (n_cpy <- 0.25 refers to cycles of 3 months)

### For Usual Care
a_Y_c_UC <- a_A_tunnels * a_R_c_UC
a_Y_u_UC <- a_A_tunnels * a_R_u_UC * n_cpy
### For New Treatment
a_Y_c_Trt <- a_A_tunnels * a_R_c_Trt
a_Y_u_Trt <- a_A_tunnels * a_R_u_Trt * n_cpy

#### Expected QALYs and Costs per cycle ####
## Vector of qalys under Usual Care
v_qaly_UC <- rowSums(t(colSums(a_Y_u_UC)))
## Vector of costs under Usual Care
v_cost_UC <- rowSums(t(colSums(a_Y_c_UC)))
## Vector of qalys under New Treatment
v_qaly_Trt <- rowSums(t(colSums(a_Y_u_Trt)))
## Vector of costs under New Treatment
v_cost_Trt <- rowSums(t(colSums(a_Y_c_Trt)))

#### Discounted total expected QALYs and Costs per strategy ####
### For Usual Care
## QALYs
n_totqaly_UC <- sum(v_qaly_UC * v_dwe * v_hcc)
## Costs
n_totcost_UC <- sum(v_cost_UC * v_dwc * v_hcc)
### For New Treatment
## QALYs
n_totqaly_Trt <- sum(v_qaly_Trt * v_dwe * v_hcc)
## Costs
n_totcost_Trt <- sum(v_cost_Trt * v_dwc * v_hcc)

########################### Cost-effectiveness analysis #######################
### Vector of total costs for both strategies
v_ted_cost <- c(n_totcost_UC, n_totcost_Trt)
### Vector of effectiveness for both strategies
v_ted_qaly <- c(n_totqaly_UC, n_totqaly_Trt)

### Calculate incremental cost-effectiveness ratios (ICERs)
df_cea <- calculate_icers(cost = v_ted_cost,
                          effect = v_ted_qaly,
                          strategies = v_names_str)
df_cea
### Create CEA table
table_cea <- df_cea
## Format column names
colnames(table_cea)[2:6] <- c("Costs ($)", "QALYs",
                              "Incremental Costs ($)", "Incremental QALYs",
                              "ICER ($/QALY)") # name the columns
## Format rows
table_cea$`Costs ($)` <- comma(round(table_cea$`Costs ($)`, 0))
table_cea$`Incremental Costs ($)`[2] <- comma(round(table_cea$`Incremental Costs ($)`[2], 0))
table_cea$QALYs <- round(table_cea$QALYs, 3)
table_cea$`Incremental QALYs` <- round(table_cea$`Incremental QALYs`, 3)
table_cea$`ICER ($/QALY)`[2] <- comma(round(table_cea$`ICER ($/QALY)`[2], 0))
table_cea
### CEA frontier
plot(df_cea) +
  expand_limits(x = 21.2)



########################### Probabalistic Sensitivty Analysis #######################

# Function to generate PSA input dataset
generate_psa_params <- function(n_sim = 1000, seed = 071818){
  set.seed(seed) # set a seed to be able to reproduce the same results
  df_psa <- data.frame(
    # Transition probabilities (per cycle)
    p_HS1   = rbeta(n_sim, 30, 170),                          # probability to become sick when healthy
    p_S1H   = rbeta(n_sim, 60, 60) ,                          # probability to become healthy when sick
    hr_S1   = rlnorm(n_sim, log(3),  0.01),                   # rate ratio of death in S1 vs healthy
    hr_S2   = rlnorm(n_sim, log(10), 0.02),                   # rate ratio of death in S2 vs healthy 
    
    # State rewards
    # Costs
    c_H   = rgamma(n_sim, shape = 100,   scale = 20),         # cost of remaining one cycle in state H
    c_S1  = rgamma(n_sim, shape = 177.8, scale = 22.5),       # cost of remaining one cycle in state S1
    c_S2  = rgamma(n_sim, shape = 225,   scale = 66.7),       # cost of remaining one cycle in state S2
    c_Trt = rgamma(n_sim, shape = 73.5,  scale = 163.3),      # cost of treatment (per cycle)
    c_D   = 0                                         ,       # cost of being in the death state
    
    # Utilities
    u_H   = rtruncnorm(n_sim, mean =    1, sd = 0.01, b = 1), # utility when healthy
    u_S1  = rtruncnorm(n_sim, mean = 0.75, sd = 0.02, b = 1), # utility when sick
    u_S2  = rtruncnorm(n_sim, mean = 0.50, sd = 0.03, b = 1), # utility when sicker
    u_D   = 0                                               , # utility when dead
    u_Trt = rtruncnorm(n_sim, mean = 0.95, sd = 0.02, b = 1)  # utility when being treated
  )
  return(df_psa)
}

# Number of simulations
n_sim <- 10000

# Generate PSA input dataset
df_psa_input <- generate_psa_params(n_sim = n_sim)
# First six observations
head(df_psa_input)

# Histogram of parameters
ggplot(melt(df_psa_input, variable.name = "Parameter"), aes(x = value)) +
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
# Source functions that contain the model and CEA output
source('functions/Functions STM_03.R')

# Run Markov model on each parameter set of PSA input dataset
for(i in 1:n_sim){
  l_out_temp <- calculate_ce_out(df_psa_input[i, ])
  df_c[i, ]  <- l_out_temp$Cost
  df_e[i, ]  <- l_out_temp$Effect
  # Display simulation progress
  if(i/(n_sim/10) == round(i/(n_sim/10),0)) { # display progress every 10%
    cat('\r', paste(i/n_sim * 100, "% done", sep = " "))
  }
}

# Create PSA object for dampack
l_psa <- make_psa_obj(cost          = df_c, 
                      effectiveness = df_e, 
                      parameters    = df_psa_input, 
                      strategies    = v_names_str)

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
plot(df_cea_psa)

## Cost-effectiveness acceptability curves (CEACs) and frontier (CEAF)
ceac_obj <- ceac(wtp = v_wtp, psa = l_psa)
# Regions of highest probability of cost-effectiveness for each strategy
summary(ceac_obj)
# CEAC & CEAF plot
plot(ceac_obj)

##  Expected Loss Curves (ELCs)
elc_obj <- calc_exp_loss(wtp = v_wtp, psa = l_psa)
elc_obj
# ELC plot
plot(elc_obj, log_y = FALSE)

## Expected value of perfect information (EVPI)
evpi <- calc_evpi(wtp = v_wtp, psa = l_psa)
# EVPI plot
plot(evpi, effect_units = "QALY")

