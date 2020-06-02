#################################### Analysis on how much memory PSA for STM_03 uses ###############################################

# Implements an age- and history-dependent Sick-Sicker cSTM model                          

##################################### Initial setup ###########################
rm(list = ls())  # remove any variables in R's memory 
library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(scales)   # for dollar signs and commas
library(truncnorm)
# devtools::install_github("DARTH-git/dampack") # to install dampack form GitHub
library(dampack)  # for CEA and calculate ICERs
library(pryr)

################################## DARTH colors  ###############################

# code for the DARTH colors for the figures
DARTHgreen      <- '#009999'  
DARTHyellow     <- '#FDAD1E'  
DARTHblue       <- '#006699' 
DARTHlightgreen <- '#00adad'
DARTHgray       <- '#666666'

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
    c_H   = rgamma(n_sim, shape = 100, scale = 20)    ,       # cost of remaining one cycle in state H
    c_S1  = rgamma(n_sim, shape = 177.8, scale = 22.5),       # cost of remaining one cycle in state S1
    c_S2  = rgamma(n_sim, shape = 225, scale = 66.7)  ,       # cost of remaining one cycle in state S2
    c_Trt = rgamma(n_sim, shape = 73.5, scale = 163.3),       # cost of treatment (per cycle)
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


## General setup
n_age_init <- 25 # age at baseline
n_age_max <- 110 # maximum age of follow up
n_t <- n_age_max - n_age_init # time horizon, number of cycles

n_tuns <- seq(from=n_t, to=(n_t*2-2), by=20) # loop through number of tunnel states
n_tuns[length(n_tuns)] <- 168
n_tuns <- n_t

n_sims <- c(1, 1000, 5000, 10000, 100000)  # loop through number of PSA simulations


############################### loop through number of tunnels and number of PSA simulations #########################################

df_mem <- as.data.frame(matrix(0, nrow=length(n_tuns), ncol=length(n_sims)))
rownames(df_mem) <- n_tuns
colnames(df_mem) <- n_sims

for (i in 1:length(n_tuns)) {

##################################### Model input #########################################

## Tunnel inputs
# Number of tunnels
n_tunnel_size <- n_tuns[[i]]
# Name for tunnels states of Sick state
v_Sick_tunnel <- paste("S1_", seq(1, n_tunnel_size), "Yr", sep = "")
# Create variables for model with tunnels
v_n_tunnels <- c("H", v_Sick_tunnel, "S2", "D") # health state names
n_states_tunnels <- length(v_n_tunnels)         # number of health states

d_c <- 0.03 # discount rate for costs 
d_e <- 0.03 # discount rate for QALYs
v_names_str <- c("Usual care", "New treatment") # store the strategy names
n_str <- length(v_names_str) # number of strategies
v_hcc    <- rep(1, n_t+1)      # vector of half-cycle correction 
v_hcc[1] <- v_hcc[n_t+1] <- 0.5

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


########################### Probabalistic Sensitivty Analysis #######################

for (j in 1:length(n_sims)) {

# Number of simulations
n_sim <- n_sims[j]

# Generate PSA input dataset
df_psa_input <- generate_psa_params(n_sim = n_sim)
# First six observations
# head(df_psa_input)


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

# Vector of memory usage
memory_ussed <- c()

# Run Markov model on each parameter set of PSA input dataset
for(k in 1:n_sim){
  l_out_temp <- calculate_ce_out(df_psa_input[k, ])
  df_c[k, ]  <- l_out_temp$Cost
  df_e[k, ]  <- l_out_temp$Effect
  # Display simulation progress
  if(k/(n_sim/10) == round(k/(n_sim/10),0)) { # display progress every 10%
    cat('\r', paste(k/n_sim * 100, "% done", sep = " "))
  }
}

#df_mem[i,j] <- sum(memory_ussed)
df_mem[i,j] <- mem_used()

} # end sim loop
} # end state loop

# df_mem <- df_mem / 1000000000

write.table(df_mem, 'PSA_sim.txt')

df_mem/1000000000















