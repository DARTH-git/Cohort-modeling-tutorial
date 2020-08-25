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
# Implements a time-independent Sick-Sicker cSTM model                        #

##################################### Initial setup ###########################
rm(list = ls())    # remove any variables in R's memory 
library(dplyr)     # to manipulate data
library(reshape2)  # to transform data
library(ggplot2)   # for nice looking plots
library(scales)    # for dollar signs and commas
# devtools::install_github("DARTH-git/dampack") # to install dampack form GitHub
library(dampack)   # for CEA and calculate ICERs

# Define required functions
calc_logit <- function(p) {
  logit <- log(p/(1-p))
  return(logit)
} 

calc_invlogit <- function(x) {
  invlogit <- exp(x)/(1+exp(x))
  return(invlogit)
}

################################## DARTH colors  ###############################

# code for the DARTH colors for the figures
DARTHgreen      <- '#009999'  
DARTHyellow     <- '#FDAD1E'  
DARTHblue       <- '#006699' 
DARTHlightgreen <- '#00adad'
DARTHgray       <- '#666666'

##################################### Model input #########################################
## General setup
n_age_init <- 25 # age at baseline
n_age_max <- 110 # maximum age of follow up
n_t <- n_age_max - n_age_init  # time horizon, number of cycles

v_n <- c("H", "S1", "S2", "D") # the 4 health states of the model:
                               # Healthy (H), Sick (S1), Sicker (S2), Dead (D)
n_states <- length(v_n)        # number of health states 
v_hcc    <- rep(1, n_t + 1)    # vector of half-cycle correction 
v_hcc[1] <- v_hcc[n_t + 1] <- 0.5 # half-cycle correction weight 
d_c <- 0.03 # discount rate for costs 
d_e <- 0.03 # discount rate for QALYs
v_names_str <- c("Usual care", "New treatment 1", "New treatment 2", "New treatments 1 & 2") # store the strategy names

## Transition probabilities (per cycle) and hazard ratios
p_HD    <- 0.002 # constant probability of dying when Healthy (all-cause mortality)
p_HS1   <- 0.15  # probability to become Sick when Healthy
p_S1H   <- 0.5   # probability to become Healthy when Sick
p_S1S2  <- 0.105 # probability to become Sicker when Sick
hr_S1   <- 3     # hazard ratio of death in Sick vs Healthy
hr_S2   <- 10    # hazard ratio of death in Sicker vs Healthy 
p_S1D   <- 1 - exp(log(1 - p_HD) * hr_S1) # probability of dying in Sick
p_S2D   <- 1 - exp(log(1 - p_HD) * hr_S2) # probability of dying in Sicker
# For New treatment 2
or_S1S2  <- 0.7              # odds ratio of becoming Sicker when Sick under New treatment 2 compared to Usual care
lor_S1S2 <- log(or_S1S2)     # log-odds ratio of becoming Sicker when Sick
logitp_S1S2 <- calc_logit(p_S1S2) # log-odds of becoming Sicker when Sick
p_S1S2_trt2 <- calc_invlogit(logitp_S1S2 + lor_S1S2) # probability to become Sicker when Sick under New treatment 2

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
du_HS1 <- 0.01 # disutility when transitioning from Healthy to Sick
ic_HS1 <- 1000 # increase in cost when transitioning from Healthy to Sick
ic_D   <- 2000 # increase in cost when dying

# Discount weight (equal discounting is assumed for costs and effects)
v_dwc <- 1 / ((1 + d_e) ^ (0:(n_t)))
v_dwe <- 1 / ((1 + d_c) ^ (0:(n_t)))

############################# Construct state-transition models ################
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
m_P["H", "H"]   <- 1 - (p_HS1 + p_HD)
m_P["H", "S1"]  <- p_HS1
m_P["H", "D"]   <- p_HD
# From S1
m_P["S1", "H"]  <- p_S1H
m_P["S1", "S1"] <- 1 - (p_S1H + p_S1S2 + p_S1D)
m_P["S1", "S2"] <- p_S1S2
m_P["S1", "D"]  <- p_S1D
# From S2
m_P["S2", "S2"] <- 1 - p_S2D
m_P["S2", "D"]  <- p_S2D
# From D
m_P["D", "D"]   <- 1

# For New treatment 2
# Only need to update the probabilities involving p_S1S2
m_P_trt2 <- m_P
m_P_trt2["S1", "S1"] <- 1 - (p_S1H + p_S1S2_trt2 + p_S1D)
m_P_trt2["S1", "S2"] <- p_S1S2_trt2

### Check if transition matrix is valid (i.e., each row should add up to 1)
valid <- rowSums(m_P) # sum the rows 
if (!isTRUE(all.equal(as.numeric((valid)), as.numeric(rep(1, n_states))))) { #check if the rows are all equal to one 
  stop("This is not a valid transition Matrix")
}
valid2 <- rowSums(m_P_trt2) # sum the rows 
if (!isTRUE(all.equal(as.numeric((valid2)), as.numeric(rep(1, n_states))))) { #check if the rows are all equal to one 
  stop("This is not a valid transition Matrix")
}

## Initialize 3D transition dynamics array
a_A <- array(0,
             dim = c(n_states, n_states, n_t + 1),
             dimnames = list(v_n, v_n, 0:n_t))
# Set first slice of A with the initial state vector in its diagonal
diag(a_A[, , 1]) <- v_s_init
# For New treatment 2
a_A_trt2 <- a_A

#### Run Markov model ####
# Iterative solution of time-independent cSTM
for(t in 1:n_t){
  # Fill in cohort trace
  m_M[t + 1, ] <- m_M[t, ] %*% m_P
  m_M_trt2[t + 1, ] <- m_M_trt2[t, ] %*% m_P_trt2
  # Fill in transition dynamics array
  a_A[, , t + 1]  <- m_M[t, ] * m_P
  a_A_trt2[, , t + 1]  <- m_M_trt2[t, ] * m_P_trt2
}

#### Plot Outputs ####
## Define colors and line types
cols <- c("H"  = DARTHgreen, "S1" = DARTHblue, 
          "S2" = DARTHyellow, "D" = DARTHgray)
lty <-  c("H" = 1, "S1" = 2, "S2" = 4, "D" = 3)
## Plot the cohort trace
ggplot(melt(m_M), aes(x = Var1, y = value, 
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

# For New Treatment 2
ggplot(melt(m_M_trt2), aes(x = Var1, y = value, 
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
## Vector of state utilities under Usual Care
v_u_UC <- c(H = u_H, S1 = u_S1, S2 = u_S2, D = u_D)
## Vector of state costs per cycle under Usual Care
v_c_UC <- c(H = c_H, S1 = c_S1, S2 = c_S2, D = c_D)

## Vector of state utilities under New treatment 1
v_u_trt1 <- c(H = u_H, S1 = u_trt1, S2 = u_S2, D = u_D)
## Vector of state costs per cycle under New treatment 1 
v_c_trt1 <- c(H = c_H, S1 = c_S1 + c_trt1, S2 = c_S2 + c_trt1, D = c_D)

## Vector of state utilities under New treatment 2
v_u_trt2 <- c(H = u_H, S1 = u_S1, S2 = u_S2, D = u_D)
## Vector of state costs per cycle under New treatment 2
v_c_trt2 <- c(H = c_H, S1 = c_S1 + c_trt2, S2 = c_S2 + c_trt2, D = c_D)

## Vector of state utilities under New treatment 1 & 2
v_u_trt1_2 <- c(H = u_H, S1 = u_trt1, S2 = u_S2, D = u_D)
## Vector of state costs per cycle under New treatment 1 & 2
v_c_trt1_2 <- c(H = c_H, S1 = c_S1 + (c_trt1 + c_trt2), S2 = c_S2 + (c_trt1 + c_trt2), D = c_D)

### Arrays of rewards
## Array of state and transition utilities under Usual Care
a_R_u_UC <- aperm(array(v_u_UC, 
                        dim = c(n_states, n_states, n_t + 1),
                        dimnames = list(v_n, v_n, 0:n_t)),
                  perm = c(2, 1, 3))
## Array of state and transition costs per cycle under Usual Care
a_R_c_UC <- aperm(array(v_c_UC, 
                        dim = c(n_states, n_states, n_t + 1),
                        dimnames = list(v_n, v_n, 0:n_t)),
                  perm = c(2, 1, 3))

## Array of utilities under New treatment 1 
a_R_u_trt1 <- aperm(array(v_u_trt1,
                         dim = c(n_states, n_states, n_t + 1),
                         dimnames = list(v_n, v_n, 0:n_t)),
                         perm = c(2, 1, 3))
## Array of costs per cycle under New treatment 1 
a_R_c_trt1 <- aperm(array(v_c_trt1,
                         dim = c(n_states, n_states, n_t + 1),
                         dimnames = list(v_n, v_n, 0:n_t)),
                         perm = c(2, 1, 3))

## Array of utilities under New treatment 2
a_R_u_trt2 <- aperm(array(v_u_trt2,
                         dim = c(n_states, n_states, n_t + 1),
                         dimnames = list(v_n, v_n, 0:n_t)),
                    perm = c(2, 1, 3))
## Array of costs under New treatment 2 
a_R_c_trt2 <- aperm(array(v_c_trt2,
                         dim = c(n_states, n_states, n_t + 1),
                         dimnames = list(v_n, v_n, 0:n_t)),
                    perm = c(2, 1, 3))

## Array of utilities under New treatment 1 & 2
a_R_u_trt1_2 <- aperm(array(v_u_trt1_2,
                          dim = c(n_states, n_states, n_t + 1),
                          dimnames = list(v_n, v_n, 0:n_t)),
                    perm = c(2, 1, 3))
## Array of costs under New treatment 1 & 2 
a_R_c_trt1_2 <- aperm(array(v_c_trt1_2,
                          dim = c(n_states, n_states, n_t + 1),
                          dimnames = list(v_n, v_n, 0:n_t)),
                    perm = c(2, 1, 3))

### Transition rewards
## For Usual Care
# Add disutility due to transition from H to S1
a_R_u_UC["H", "S1", ] <- a_R_u_UC["H", "S1", ] - du_HS1
# Add transition cost per cycle due to transition from H to S1
a_R_c_UC["H", "S1", ] <- a_R_c_UC["H", "S1", ] + ic_HS1
# Add transition cost  per cycle of dying from all non-dead states
a_R_c_UC[-n_states, "D", ] <- a_R_c_UC[-n_states, "D", ] + ic_D

## For New treatment 1
# Add disutility due to transition from Healthy to Sick
a_R_u_trt1["H", "S1", ] <- a_R_u_trt1["H", "S1", ] - du_HS1
# Add transition cost per cycle due to transition from Healthy to Sick
a_R_c_trt1["H", "S1", ] <- a_R_c_trt1["H", "S1", ] + ic_HS1
# Add transition cost per cycle of dying from all non-dead states
a_R_c_trt1[-n_states, "D", ] <- a_R_c_trt1[-n_states, "D", ] + ic_D

# For New treatment 2
# Add disutility due to transition from Healthy to Sick
a_R_u_trt2["H", "S1", ] <- a_R_u_trt2["H", "S1", ] - du_HS1
# Add transition cost per cycle due to transition from Healthy to Sick
a_R_c_trt2["H", "S1", ] <- a_R_c_trt2["H", "S1", ] + ic_HS1
# Add transition cost per cycle of dying from all non-dead states
a_R_c_trt2[-n_states, "D", ] <- a_R_c_trt2[-n_states, "D", ] + ic_D

# For New treatment 1 & 2
# Add disutility due to transition from Healthy to Sick
a_R_u_trt1_2["H", "S1", ] <- a_R_u_trt1_2["H", "S1", ] - du_HS1
# Add transition cost per cycle due to transition from Healthy to Sick
a_R_c_trt1_2["H", "S1", ] <- a_R_c_trt1_2["H", "S1", ] + ic_HS1
# Add transition cost per cycle of dying from all non-dead states
a_R_c_trt1_2[-n_states, "D", ] <- a_R_c_trt1_2[-n_states, "D", ] + ic_D

#### Expected QALYs and Costs for all transitions per cycle ####
# QALYs = life years x QoL
# NOTE: all parameters are annual

### For Usual Care
a_Y_c_UC <- a_A * a_R_c_UC
a_Y_u_UC <- a_A * a_R_u_UC 

### For New treatment 1 
a_Y_c_trt1 <- a_A * a_R_c_trt1
a_Y_u_trt1 <- a_A * a_R_u_trt1 

### For New treatment 2
a_Y_c_trt2 <- a_A_trt2 * a_R_c_trt2
a_Y_u_trt2 <- a_A_trt2 * a_R_u_trt2 

### For New treatment 1 & 2
a_Y_c_trt1_2 <- a_A_trt2 * a_R_c_trt1_2
a_Y_u_trt1_2 <- a_A_trt2 * a_R_u_trt1_2 

#### Expected QALYs and Costs per cycle ####
## Vector of qalys under Usual Care
v_qaly_UC <- rowSums(t(colSums(a_Y_u_UC)))
## Vector of costs under Usual Care
v_cost_UC <- rowSums(t(colSums(a_Y_c_UC)))

## Vector of qalys under New Treatment 1
v_qaly_trt1 <- rowSums(t(colSums(a_Y_u_trt1)))
## Vector of costs under New Treatment 1
v_cost_trt1 <- rowSums(t(colSums(a_Y_c_trt1)))

## Vector of qalys under New Treatment 2
v_qaly_trt2 <- rowSums(t(colSums(a_Y_u_trt2)))
## Vector of costs under New Treatment 2
v_cost_trt2 <- rowSums(t(colSums(a_Y_c_trt2)))

## Vector of qalys under New Treatment 1 & 2
v_qaly_trt1_2 <- rowSums(t(colSums(a_Y_u_trt1_2)))
## Vector of costs under New Treatment 1 & 2
v_cost_trt1_2 <- rowSums(t(colSums(a_Y_c_trt1_2)))

#### Discounted total expected QALYs and Costs per strategy ####
### For Usual Care
## QALYs
n_totqaly_UC <- t(v_qaly_UC) %*% (v_dwe * v_hcc)
## Costs
n_totcost_UC <- t(v_cost_UC) %*% (v_dwc * v_hcc)

### For New treatment 1
## QALYs
n_totqaly_trt1 <- t(v_qaly_trt1) %*% (v_dwe * v_hcc)
## Costs
n_totcost_trt1 <- t(v_cost_trt1) %*% (v_dwc * v_hcc)

### For New treatment 2
## QALYs
n_totqaly_trt2 <- t(v_qaly_trt2) %*% (v_dwe * v_hcc)
## Costs
n_totcost_trt2 <- t(v_cost_trt2) %*% (v_dwc * v_hcc)

### For New treatment 1 & 2
## QALYs
n_totqaly_trt1_2 <- t(v_qaly_trt1_2) %*% (v_dwe * v_hcc)
## Costs
n_totcost_trt1_2 <- t(v_cost_trt1_2) %*% (v_dwc * v_hcc)


########################### Cost-effectiveness analysis #######################
### Vector of total costs for all strategies
v_ted_cost <- c(n_totcost_UC, n_totcost_trt1, n_totcost_trt2, n_totcost_trt1_2)
### Vector of effectiveness for all strategies
v_ted_qaly <- c(n_totqaly_UC, n_totqaly_trt1, n_totqaly_trt2, n_totqaly_trt1_2)

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
table_cea$`Incremental Costs ($)` <- comma(round(table_cea$`Incremental Costs ($)`, 0))
table_cea$QALYs <- round(table_cea$QALYs, 2)
table_cea$`Incremental QALYs` <- round(table_cea$`Incremental QALYs`, 2)
table_cea$`ICER ($/QALY)` <- comma(round(table_cea$`ICER ($/QALY)`, 0))
table_cea
### CEA frontier
plot(df_cea, label="all") +
     expand_limits(x = 23) # change this

