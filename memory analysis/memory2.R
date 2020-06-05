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

## General setup
n_age_init <- 25 # age at baseline
n_age_max <- 640 # maximum age of follow up
n_t <- n_age_max - n_age_init # time horizon, number of cycles

n_tuns <- 2:(n_t*2-2)

################################################### loop through number of states ################################################

df_mem <- as.data.frame(matrix(0, nrow=length(n_tuns), ncol=1))
rownames(df_mem) <- n_tuns
colnames(df_mem) <- 'memory'

# for (k in 1:length(n_tuns)) {
  
k <- length(n_tuns)

##################################### Model input #########################################

## Tunnel inputs
# Number of tunnels
n_tunnel_size <- n_tuns[k]
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

# expand v_r_mort_by_age to fit the tunnel state explosion
v_r_mort_by_age <- c(v_r_mort_by_age, rep(v_r_mort_by_age[length(v_r_mort_by_age)], 2*n_t))

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

# size of the 3-D array in GB
as.numeric(object.size(a_P_tunnels))/1000000000

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
v_u_S1_Tr <- rep(u_Trt, n_tunnel_size)
names(v_u_S1_Tr) <- v_Sick_tunnel
## Vector of state utilities under New Treatment
v_u_Tr <- c(H = u_H, v_u_S1_Tr, S2 = u_S2, D = u_D)

## Vector of costs for S1 under Usual Care
v_c_S1_Tr <- rep(c_S1 + c_Trt, n_tunnel_size)
names(v_c_S1_Tr) <- v_Sick_tunnel
## Vector of state costs under New Treatment
v_c_Tr <- c(H = c_H, v_c_S1_Tr, S2 = c_S2 + c_Trt, D = c_D)

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
a_R_u_Tr <- aperm(array(v_u_Tr,
                        dim = c(n_states_tunnels, n_states_tunnels, n_t + 1),
                        dimnames = list(v_n_tunnels, v_n_tunnels, 0:n_t)),
                  perm = c(2, 1, 3))
## Array of costs under New Treatment
a_R_c_Tr <- aperm(array(v_c_Tr,
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
a_R_u_Tr["H", "S1_1Yr", ] <- a_R_u_Tr["H", "S1_1Yr", ] - du_HS1
# Add transition cost due to transition from Healthy to Sick
a_R_c_Tr["H", "S1_1Yr", ] <- a_R_c_Tr["H", "S1_1Yr", ] + ic_HS1
# Add transition cost of dying from all non-dead states
a_R_c_Tr[-n_states_tunnels, "D", ] <- a_R_c_Tr[-n_states_tunnels, "D", ] + ic_D

#### Expected QALYs and Costs for all transitions per cycle ####
### For Usual Care
a_Y_c_UC <- a_A_tunnels * a_R_c_UC
a_Y_u_UC <- a_A_tunnels * a_R_u_UC
### For New Treatment
a_Y_c_Tr <- a_A_tunnels * a_R_c_Tr
a_Y_u_Tr <- a_A_tunnels * a_R_u_Tr

#### Expected QALYs and Costs per cycle ####
## Vector of qalys under Usual Care
v_qaly_UC <- rowSums(t(colSums(a_Y_u_UC)))
## Vector of costs under Usual Care
v_cost_UC <- rowSums(t(colSums(a_Y_c_UC)))
## Vector of qalys under New Treatment
v_qaly_Tr <- rowSums(t(colSums(a_Y_u_Tr)))
## Vector of costs under New Treatment
v_cost_Tr <- rowSums(t(colSums(a_Y_c_Tr)))

#### Discounted total expected QALYs and Costs per strategy ####
### For Usual Care
## QALYs
n_totqaly_UC <- sum(v_qaly_UC * v_dwe * v_hcc)
## Costs
n_totcost_UC <- sum(v_cost_UC * v_dwc * v_hcc)
### For New Treatment
## QALYs
n_totqaly_Tr <- sum(v_qaly_Tr * v_dwe * v_hcc)
## Costs
n_totcost_Tr <- sum(v_cost_Tr * v_dwc * v_hcc)

########################### Cost-effectiveness analysis #######################
### Vector of total costs for both strategies
v_ted_cost <- c(n_totcost_UC, n_totcost_Tr)
### Vector of effectiveness for both strategies
v_ted_qaly <- c(n_totqaly_UC, n_totqaly_Tr)

### Calculate incremental cost-effectiveness ratios (ICERs)
df_cea <- calculate_icers(cost = v_ted_cost,
                          effect = v_ted_qaly,
                          strategies = v_names_str)

df_mem[k,] <- mem_used()

# } # end state loop


df_mem <- df_mem / 1000000000

write.table(df_mem, 'model_sim.txt')


##### EXTRAPOLATION
df_mem <- read.csv('memory analysis/model_sim.txt', header=T)

library(fitdistrplus)
x <- as.numeric(df_mem$memory)
descdist(x, discrete = FALSE)


fit.exp <- fitdist(x, "exp")
plot(fit.exp)

fit.weibull <- fitdist(x, "weibull")
plot(fit.weibull)

fit.gamma <- fitdist(x, "gamma")
plot(fit.gamma)

fit.lognormal <- fitdist(x, "lnorm")
plot(fit.lognormal)

fit.logistic <- fitdist(x, "logis")
plot(fit.logistic)

fit.normal <- fitdist(x, "norm")
plot(fit.normal)







