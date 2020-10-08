#-----------------------------------------#
####          Decision Model           ####
#-----------------------------------------#
#' Decision Model
#'
#' \code{decision_model} implements the decision model used.
#'
#' @param l_params_all List with all parameters of decision model
#' @param verbose Logical variable to indicate print out of messages
#' @return The transition probability array and the cohort trace matrix.
#' 
decision_model <- function(l_params_all, verbose = FALSE) {
  with(as.list(l_params_all), {
    
    ###################### Process model inputs ######################
    ## Age-specific transiition probability to the Dead state
    # compute mortality rates
    v_r_S1Dage <- v_r_HDage * hr_S1        # Age-specific mortality rate in the Sick state 
    v_r_S2Dage <- v_r_HDage * hr_S2        # Age-specific mortality rate in the Sicker state 
    # transform rates to probabilities
    v_p_S1Dage <- rate_to_prob(v_r_S1Dage) # Age-specific mortality risk in the Sick state
    v_p_S2Dage <- rate_to_prob(v_r_S2Dage) # Age-specific mortality risk in the Sicker state
    
    ## History-dependent transition probability of becoming Sicker when Sick
    # conditional on surviving
    # Weibull hazard
    v_p_S1S2_tunnels <- n_lambda * n_gamma * (1:n_tunnel_size)^{n_gamma-1}
    # transform odds ratios to probabilites for New treatment 2
    # vector of log-odds of becoming Sicker when Sick
    v_logit_S1S2_tunnels <- boot::logit(v_p_S1S2_tunnels) 
    # vector with probabilities of becoming Sicker when Sick under New treatment 2 
    # conditional on surviving
    v_p_S1S2_tunnels_trt2 <- boot::inv.logit(v_logit_S1S2_tunnels + lor_S1S2) 
    
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
                       v_Sick_tunnel[i + 1], ]   <- (1 - v_p_S1Dage) * 
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
    l_m_M <- list(m_M_tunnels_sum,
                  m_M_tunnels_sum,
                  m_M_tunnels_sum_trt2,
                  m_M_tunnels_sum_trt2)
    names(l_m_M) <- v_names_str
    
    ## Store the transition array for each strategy in a list
    l_a_A <- list(a_A_tunnels,
                  a_A_tunnels,
                  a_A_tunnels_trt2,
                  a_A_tunnels_trt2)
    names(l_m_M) <- names(l_a_A) <- v_names_str
    
    ########################################## RETURN OUTPUT  ##########################################
    out <- list(l_m_M = l_m_M,
                l_a_A = l_a_A)
    
    return(out)
  }
  )
}

#---------------------------------------------#
#### Calculate cost-effectiveness outcomes ####
#---------------------------------------------#
#' Calculate cost-effectiveness outcomes
#'
#' \code{calculate_ce_out} calculates costs and effects for a given vector of parameters using a simulation model.
#' @param l_params_all List with all parameters of decision model
#' @param n_wtp Willingness-to-pay threshold to compute net benefits
#' @return A data frame with discounted costs, effectiveness and NMB.
#' 
calculate_ce_out <- function(l_params_all, n_wtp = 100000){ # User defined
  with(as.list(l_params_all), {
    
    ### Run decision model to get transition dynamics array
    model <- decision_model(l_params_all = l_params_all)
    l_a_A <- model[["l_a_A"]]
    
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
    
    # assign strategy names to matching items in the lists
    names(l_u) <- names(l_c) <- v_names_str
    
    ## create empty vectors to store total utilities and costs 
    v_totqaly <- v_totcost <- vector(mode = "numeric", length = n_str)
    names(v_totqaly) <- names(v_totcost) <- v_names_str
    
    #### Loop through each strategy and calculate total utilities and costs ####
    for (i in 1:n_str) {
      v_u             <- l_u[[i]]   # select the vector of state utilities for the ith strategy
      v_c             <- l_c[[i]]   # select the vector of state costs for the ith strategy
      a_A_tunnels_str <- l_a_A[[i]] # select the transition array for the ith strategy
      
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
      a_Y_c <- a_A_tunnels_str * a_R_c
      a_Y_u <- a_A_tunnels_str * a_R_u 
      
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
    df_cea <- df_cea[match(v_names_str, df_cea$Strategy), ]
    return(df_cea)
  }
  )
}
