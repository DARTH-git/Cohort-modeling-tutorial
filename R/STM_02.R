###############################################################################
### Implementation of cohort state-transition models in R:                   ##
### From conceptualization to implementation 2020                            ##
###############################################################################
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

###############################################################################
################# Code of Appendix   ##########################################
###############################################################################
# Implements a time-independent Sick-Sicker cSTM model                        #
# Usual care: best available care for the patients with the disease. This scenario reflects the natural history of the disease progressions
# New treatment 1: this new treatment is given to all sick patients, patients in sick and sicker, but does only improves the utility of those being sick.
# New treatment 2: the new treatment reduces disease progression from sick to sicker. However, it is not possible to distinguish those sick from sicker and therefore all individuals in one of the two sick states get the treatment.  
# New treatment 1 & new treatment 2”: This strategy combines the new treatment 1 and new treatment 2. The disease progression is reduced and Sick individuals has an improved utility. 
# This model incorporates time-dependent transition probabilities 

##################################### Initial setup ###########################
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
source("functions/Functions.R")

##################################### Model input ##############################
## General setup
n_age_init <- 25 # age at baseline
n_age_max <- 110 # maximum age of follow up
n_t <- n_age_max - n_age_init  # time horizon, number of cycles

v_n <- c("H", "S1", "S2", "D") # the 4 health states of the model:
                               # Healthy (H), Sick (S1), Sicker (S2), Dead (D)
v_hcc    <- rep(1, n_t+1)      # vector of half-cycle correction 
v_hcc[1] <- v_hcc[n_t+1] <- 0.5 # half-cycle correction weight 
n_states <- length(v_n) # number of health states 
d_c <- 0.03 # discount rate for costs 
d_e <- 0.03 # discount rate for QALYs
v_names_str <- c("Usual care", "New treatment 1", "New treatment 2", "New treatments 1 & 2") # store the strategy names
n_str <- length(v_names_str) # number of strategies

## Transition probabilities (per cycle) and hazard ratios
p_HS1   <- 0.15  # probability to become Sick when Healthy
p_S1H   <- 0.5   # probability to become Healthy when Sick
p_S1S2  <- 0.105 # probability to become Sicker when Sick
hr_S1   <- 3     # hazard ratio of death in Sick vs Healthy
hr_S2   <- 10    # hazard ratio of death in Sicker vs Healthy 
# For New treatment 2
or_S1S2  <- 0.7              # odds ratio of becoming Sicker when Sick under New treatment 2
lor_S1S2 <- log(or_S1S2)     # log-odd ratio of becoming Sicker when Sick
logit_S1S2 <- logit(p_S1S2)  # log-odds of becoming Sicker when Sick
p_S1S2_trt2 <- inv.logit(logit_S1S2 + lor_S1S2) # probability to become Sicker when Sick under New treatment 2

## Age-dependent mortality rates
lt_usa_2005 <- read.csv("data/LifeTable_USA_Mx_2015.csv")
v_r_mort_by_age <- lt_usa_2005 %>% 
  # filter(Age >= age & Age <= n_age_max) %>%  # WE DONT NEED THIS RIGHT?
  select(Total) %>%
  as.matrix()

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

## Age-specific transition probabilities
# Age-specific probability of dying when Healthy (all-cause mortality)
p_HDage  <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_t - 1)])  # NOTE two steps + rate prob function
# Age-specific mortality risk in the Sick state
p_S1Dage <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_t - 1)] * hr_S1)
# Age-specific mortality risk in the Sicker state
p_S2Dage <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_t - 1)] * hr_S2)

######################## Construct state-transition models ####################
#### Create transition arrays ####
# Initialize 3-D array
a_P <- array(0, dim = c(n_states, n_states, n_t),
             dimnames = list(v_n, v_n, 0:(n_t - 1)))
### Fill in array
## From H
a_P["H", "H", ]   <- 1 - (p_HS1 + p_HDage)
a_P["H", "S1", ]  <- p_HS1
a_P["H", "D", ]   <- p_HDage
## From S1
a_P["S1", "H", ]  <- p_S1H
a_P["S1", "S1", ] <- 1 - (p_S1H + p_S1S2 + p_S1Dage)
a_P["S1", "S2", ] <- p_S1S2
a_P["S1", "D", ]  <- p_S1Dage
## From S2
a_P["S2", "S2", ] <- 1 - p_S2Dage
a_P["S2", "D", ]  <- p_S2Dage
## From D
a_P["D", "D", ]   <- 1

# For New treatment 2
# Only need to update the probabilities involving p_S1S2
a_P_trt2 <- a_P
a_P_trt2["S1", "S1", ] <- 1 - (p_S1H + p_S1S2_trt2 + p_S1Dage)
a_P_trt2["S1", "S2", ] <- p_S1S2_trt2

# ### Check if transition matrix is valid (i.e., each row should add up to 1) ## NOTE: put in functions.R
valid <- apply(a_P, 3, function(x) sum(rowSums(x))==n_states)
if (!isTRUE(all.equal(as.numeric(sum(valid)), as.numeric(n_t)))) {
  stop("This is not a valid transition Matrix")
}
valid2 <- apply(a_P_trt2, 3, function(x) sum(rowSums(x))==n_states)
if (!isTRUE(all.equal(as.numeric(sum(valid)), as.numeric(n_t)))) {
  stop("This is not a valid transition Matrix")
}

#### Run Markov model ####
## Initial state vector
# All starting healthy
v_s_init <- c(H = 1, S1 = 0, S2 = 0, D = 0) # initial state vector
v_s_init

## Initialize cohort trace for age-dependent cSTM
m_M_ad <- matrix(0, 
                 nrow = (n_t + 1), ncol = n_states, 
                 dimnames = list(0:n_t, v_n))
# Store the initial state vector in the first row of the cohort trace
m_M_ad[1, ] <- v_s_init
# For New treatment 2
m_M_ad_trt2 <- m_M_ad

## Initialize transition array
a_A <- array(0,
             dim = c(n_states, n_states, n_t + 1),
             dimnames = list(v_n, v_n, 0:n_t))
# Set first slice of A with the initial state vector in its diagonal
diag(a_A[, , 1]) <- v_s_init
# For New treatment 2
a_A_trt2 <- a_A

## Iterative solution of age-dependent cSTM
for(t in 1:n_t){
  # Fill in cohort trace
  m_M_ad[t + 1, ]      <- m_M_ad[t, ]      %*% a_P[, , t]
  m_M_ad_trt2[t + 1, ] <- m_M_ad_trt2[t, ] %*% a_P_trt2[, , t]
  # Fill in transition dynamics array
  a_A[, , t + 1]       <- m_M_ad[t, ]       * a_P[, , t]
  a_A_trt2[, , t + 1]  <- m_M_ad_trt2[t, ]  * a_P_trt2[, , t]
}

#### Plot Outputs ####
### Cohort trace
## Plot the cohort trace for scenarios Usual care and New treatment 1 
plot_trace(m_M_ad)
## Plot the cohort trace for scenarios New treatment 2 and New treatments 1 & 2
plot_trace(m_M_ad_trt2)

### DO YOU GUYS THINK THESE EPI OUTCOMES ARE NECESSARY?
### Survival
v_S_ad <- rowSums(m_M_ad[, -4]) # vector with survival curve
## ggplot
ggplot(data.frame(Cycle = 0:n_t, Survival = v_S_ad), 
          aes(x = Cycle, y = Survival)) +
  geom_line(size = 1.3) +
  xlab("Cycle") +
  ylab("Proportion alive") +
  theme_bw(base_size = 14) +
  theme()

# For New treatment 2
v_S_ad_trt2 <- rowSums(m_M_ad_trt2[, -4]) # vector with survival curve
## ggplot
ggplot(data.frame(Cycle = 0:n_t, Survival = v_S_ad_trt2), 
          aes(x = Cycle, y = Survival)) +
  geom_line(size = 1.3) +
  xlab("Cycle") +
  ylab("Proportion alive") +
  theme_bw(base_size = 14) +
  theme()

### Prevalence
v_prev_S1   <- m_M_ad[, "S1"] / v_S_ad # vector with prevalence of Sick
v_prev_S2   <- m_M_ad[, "S2"] / v_S_ad # vector with prevalence of Sicker
v_prev_S1S2 <- rowSums(m_M_ad[, c("S1", "S2")]) / v_S_ad # vector with prevalence of Sick and Sicker
## Data.frame with all prevalence
df_prev_states <- data.frame(Cycle = 0:n_t, 
                             States  = ordered(rep(c("S1", "S2", "S1 and S2"),
                                                   each = (n_t + 1)), 
                                               levels = c("S1", "S2", "S1 and S2")), 
                             Prevalence = c(v_prev_S1, 
                                            v_prev_S2, 
                                            v_prev_S1S2))
## ggplot
ggplot(df_prev_states, 
       aes(x = Cycle, y = Prevalence, 
           color = States, linetype = States)) +
  geom_line(size = 1) +
  scale_y_continuous(labels = scales::percent) + 
  scale_color_discrete(name = "Health State", l = 50) +
  scale_linetype(name = "Health State") +
  xlab("Cycle") +
  ylab("Prevalence (%)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

# For New treatment 2
v_prev_S1_trt2   <- m_M_ad_trt2[, "S1"] / v_S_ad_trt2 # vector with prevalence of Sick
v_prev_S2_trt2   <- m_M_ad_trt2[, "S2"] / v_S_ad_trt2 # vector with prevalence of Sicker
v_prev_S1S2_trt2 <- rowSums(m_M_ad_trt2[, c("S1", "S2")]) / v_S_ad_trt2 # vector with prevalence of Sick and Sicker
## Data.frame with all prevalence
df_prev_states_trt2 <- data.frame(Cycle = 0:n_t, 
                                  States  = ordered(rep(c("S1", "S2", "S1 and S2"),
                                                   each = (n_t + 1)), 
                                               levels = c("S1", "S2", "S1 and S2")), 
                                  Prevalence = c(v_prev_S1_trt2, 
                                                 v_prev_S2_trt2, 
                                                 v_prev_S1S2_trt2))
## ggplot
ggplot(df_prev_states_trt2, 
       aes(x = Cycle, y = Prevalence, 
           color = States, linetype = States)) +
  geom_line(size = 1) +
  scale_y_continuous(labels = scales::percent) + 
  scale_color_discrete(name = "Health State", l = 50) +
  scale_linetype(name = "Health State") +
  xlab("Cycle") +
  ylab("Prevalence (%)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

### Proportion of Sicker among sick individuals
v_prop_S2 <- m_M_ad[-1, "S2"] / v_prev_S1S2[-1] 
## ggplot
ggplot(data.frame(Cycle = 1:n_t, 
                  Proportion = v_prop_S2), 
       aes(x = Cycle, y = Proportion)) +
  geom_line(size = 1) +
  # scale_y_continuous(labels = scales::percent) + 
  xlab("Cycle") +
  ylab("Proportion") +
  theme_bw(base_size = 14) +
  theme()

# For New treatment2 
v_prop_S2_trt2 <- m_M_ad_trt2[-1, "S2"] / v_prev_S1S2_trt2[-1] 
## ggplot
ggplot(data.frame(Cycle = 1:n_t, 
                  Proportion = v_prop_S2_trt2), 
       aes(x = Cycle, y = Proportion)) +
  geom_line(size = 1) +
  # scale_y_continuous(labels = scales::percent) + 
  xlab("Cycle") +
  ylab("Proportion") +
  theme_bw(base_size = 14) +
  theme()

### Life expectancy
le_ad <- sum(v_S_ad)
le_ad_trt2 <- sum(v_S_ad_trt2)

#### State and Transition Rewards ####
### State rewards
## Vector of state utilities under Usual care
v_u_UC <- c(H = u_H, S1 = u_S1, S2 = u_S2, D = u_D)
## Vector of state costs per cycle under Usual care
v_c_UC <- c(H = c_H, S1 = c_S1, S2 = c_S2, D = c_D)

## Vector of state utilities under New treatment 1
v_u_trt1 <- c(H = u_H, S1 = u_trt1, S2 = u_S2, D = u_D)
## Vector of state costs per cycle under New Treatment
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

v_qaly_UC1 <- apply(a_Y_u_UC, 3, sum) # use this

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
df_cea <- calculate_icers(cost = v_totcost, 
                          effect = v_totqaly,
                          strategies = v_names_str)
df_cea
### Create CEA table
table_cea <- df_cea
## Format column names
colnames(table_cea)[colnames(table_cea) %in% c("Cost", "Effect", "Inc_Cost", "Inc_Effect", "ICER")] <- 
            c("Costs ($)", "QALYs", "Incremental Costs ($)", "Incremental QALYs", "ICER ($/QALY)") 
## Format rows
table_cea$`Costs ($)` <- comma(round(table_cea$`Costs ($)`, 0))
table_cea$`Incremental Costs ($)` <- comma(round(table_cea$`Incremental Costs ($)`, 0))
table_cea$QALYs <- round(table_cea$QALYs, 2)
table_cea$`Incremental QALYs` <- round(table_cea$`Incremental QALYs`, 2)
table_cea$`ICER ($/QALY)` <- comma(round(table_cea$`ICER ($/QALY)`, 0))
table_cea
### CEA frontier
plot(df_cea, label = "all") +
  expand_limits(x = max(table_cea$QALYs + 0.5)) # change this


