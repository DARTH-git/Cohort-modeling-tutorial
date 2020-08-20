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
    
    ## History-dependent transition from S1 to S2
    # Create tunnel states
    p_S1S2_tunnels <- n_lambda * n_gamma * (1:n_tunnel_size)^{n_gamma-1}
    logitp_S1S2 <- logit(p_S1S2_tunnels) # log-odds of becoming Sicker when Sick
    p_S1S2_tunnels_trt2 <- invlogit(logitp_S1S2 + lor_S1S2) # probability to become Sicker when Sick under New treatment 2
    
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
    # repeat code for the last cycle to force the cohort stay in the last tunnel state of Sick
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
    
    # For New treatment 2
    # Only need to update the probabilities involving p_S1S2
    a_P_tunnels_trt2 <- a_P_tunnels
    for(i in 1:(n_tunnel_size - 1)){
      a_P_tunnels_trt2[v_Sick_tunnel[i], "H", ]  <- p_S1H
      a_P_tunnels_trt2[v_Sick_tunnel[i], 
                       v_Sick_tunnel[i + 1], ]   <- 1 - (p_S1H + p_S1S2_tunnels_trt2[i] + v_p_S1Dage)
      a_P_tunnels_trt2[v_Sick_tunnel[i], "S2", ] <- p_S1S2_tunnels_trt2[i]
      a_P_tunnels_trt2[v_Sick_tunnel[i], "D", ]  <- v_p_S1Dage
    }
    a_P_tunnels_trt2[v_Sick_tunnel[n_tunnel_size],
                     v_Sick_tunnel[n_tunnel_size], ] <- 1 - (p_S1H +
                                                             p_S1S2_tunnels_trt2[n_tunnel_size] + 
                                                             v_p_S1Dage)
    a_P_tunnels_trt2[v_Sick_tunnel[n_tunnel_size], "S2", ] <- p_S1S2_tunnels_trt2[n_tunnel_size]
    
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
      m_M_tunnels[t + 1, ] <- m_M_tunnels[t, ] %*% a_P_tunnels[, , t]
      m_M_tunnels_trt2[t + 1,] <- m_M_tunnels_trt2[t, ] %*% a_P_tunnels_trt2[, , t]
      # Fill in transition dynamics array
      a_A_tunnels[, , t + 1]  <- m_M_tunnels[t, ] * a_P_tunnels[, , t]
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
    
    ########################################## RETURN OUTPUT  ##########################################
    out <- list(m_M_tunnels = m_M_tunnels,
                a_A_tunnels = a_A_tunnels,
                m_M_tunnels_sum = m_M_tunnels_sum,
                m_M_tunnels_trt2 = m_M_tunnels_trt2,
                a_A_tunnels_trt2 = a_A_tunnels_trt2,
                m_M_tunnels_sum_trt2 = m_M_tunnels_sum_trt2)
    
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
    a_A_tunnels <- model[["a_A_tunnels"]]
    a_A_tunnels_trt2 <- model[["a_A_tunnels_trt2"]]
    
    #### State and Transition Rewards ####
    ### State rewards
    ## Vector of utilities for S1 under Usual care
    v_u_S1_UC <- rep(u_S1, n_tunnel_size)
    names(v_u_S1_UC) <- v_Sick_tunnel
    ## Vector of state utilities under Usual care
    v_u_UC <- c(H = u_H, v_u_S1_UC, S2 = u_S2, D = u_D)
    ## Vector of costs for S1 under Usual care
    v_c_S1_UC <- rep(c_S1, n_tunnel_size)
    names(v_c_S1_UC) <- v_Sick_tunnel
    ## Vector of state costs under Usual care
    v_c_UC <- c(H = c_H, v_c_S1_UC, S2 = c_S2, D = c_D)
    
    ## Vector of utilities for S1 under New treatment 1
    v_u_S1_Trt <- rep(u_Trt, n_tunnel_size)
    names(v_u_S1_Trt) <- v_Sick_tunnel
    ## Vector of state utilities under New treatment 1
    v_u_Trt <- c(H = u_H, v_u_S1_Trt, S2 = u_S2, D = u_D)
    ## Vector of costs for S1 under New treatment 1
    v_c_S1_Trt <- rep(c_S1 + c_Trt, n_tunnel_size)
    names(v_c_S1_Trt) <- v_Sick_tunnel
    ## Vector of state costs under New treatment 1
    v_c_Trt <- c(H = c_H, v_c_S1_Trt, S2 = c_S2 + c_Trt, D = c_D)
    
    ## Vector of utilities for S1 under New treatment 2
    v_u_S1_Trt2 <- rep(u_Trt, n_tunnel_size)
    names(v_u_S1_Trt2) <- v_Sick_tunnel
    ## Vector of state utilities under New treatment 2
    v_u_Trt2 <- c(H = u_H, v_u_S1_Trt2, S2 = u_S2, D = u_D)
    ## Vector of costs for S1 under New treatment 2
    v_c_S1_Trt2 <- rep(c_S1 + c_Trt, n_tunnel_size)
    names(v_c_S1_Trt2) <- v_Sick_tunnel
    ## Vector of state costs under New treatment 2
    v_c_Trt2 <- c(H = c_H, v_c_S1_Trt2, S2 = c_S2 + c_Trt, D = c_D)
    
    ### Arrays of rewards
    ## Array of state and transition utilities under Usual care
    a_R_u_UC <- aperm(array(v_u_UC,
                            dim = c(n_states_tunnels, n_states_tunnels, n_t + 1),
                            dimnames = list(v_n_tunnels, v_n_tunnels, 0:n_t)),
                      perm = c(2, 1, 3))
    ## Array of state and transition costs under Usual care
    a_R_c_UC <- aperm(array(v_c_UC,
                            dim = c(n_states_tunnels, n_states_tunnels, n_t + 1),
                            dimnames = list(v_n_tunnels, v_n_tunnels, 0:n_t)),
                      perm = c(2, 1, 3))
    
    ## Array of utilities under New treatment 1
    a_R_u_Trt <- aperm(array(v_u_Trt,
                             dim = c(n_states_tunnels, n_states_tunnels, n_t + 1),
                             dimnames = list(v_n_tunnels, v_n_tunnels, 0:n_t)),
                       perm = c(2, 1, 3))
    ## Array of costs under New treatment 1
    a_R_c_Trt <- aperm(array(v_c_Trt,
                             dim = c(n_states_tunnels, n_states_tunnels, n_t + 1),
                             dimnames = list(v_n_tunnels, v_n_tunnels, 0:n_t)),
                       perm = c(2, 1, 3))
    
    ## Array of utilities under New treatment 2
    a_R_u_Trt2 <- aperm(array(v_u_Trt2,
                              dim = c(n_states_tunnels, n_states_tunnels, n_t + 1),
                              dimnames = list(v_n_tunnels, v_n_tunnels, 0:n_t)),
                        perm = c(2, 1, 3))
    ## Array of costs under New treatment 2
    a_R_c_Trt2 <- aperm(array(v_c_Trt2,
                              dim = c(n_states_tunnels, n_states_tunnels, n_t + 1),
                              dimnames = list(v_n_tunnels, v_n_tunnels, 0:n_t)),
                        perm = c(2, 1, 3))
    
    ### Transition rewards
    ## For Usual care
    # Add disutility due to transition from H to S1
    a_R_u_UC["H", "S1_1Yr", ] <- a_R_u_UC["H", "S1_1Yr", ] - du_HS1
    # Add transition cost due to transition from H to S1
    a_R_c_UC["H", "S1_1Yr", ] <- a_R_c_UC["H", "S1_1Yr", ] + ic_HS1
    # Add transition cost of dying from all non-dead states
    a_R_c_UC[-n_states_tunnels, "D", ] <- a_R_c_UC[-n_states_tunnels, "D", ] + ic_D
    
    ## For New treatment 1
    # Add disutility due to transition from Healthy to Sick
    a_R_u_Trt["H", "S1_1Yr", ] <- a_R_u_Trt["H", "S1_1Yr", ] - du_HS1
    # Add transition cost due to transition from Healthy to Sick
    a_R_c_Trt["H", "S1_1Yr", ] <- a_R_c_Trt["H", "S1_1Yr", ] + ic_HS1
    # Add transition cost of dying from all non-dead states
    a_R_c_Trt[-n_states_tunnels, "D", ] <- a_R_c_Trt[-n_states_tunnels, "D", ] + ic_D
    
    ## For New treatment 2
    # Add disutility due to transition from Healthy to Sick
    a_R_u_Trt2["H", "S1_1Yr", ] <- a_R_u_Trt2["H", "S1_1Yr", ] - du_HS1
    # Add transition cost due to transition from Healthy to Sick
    a_R_c_Trt2["H", "S1_1Yr", ] <- a_R_c_Trt2["H", "S1_1Yr", ] + ic_HS1
    # Add transition cost of dying from all non-dead states
    a_R_c_Trt2[-n_states_tunnels, "D", ] <- a_R_c_Trt2[-n_states_tunnels, "D", ] + ic_D
    
    #### Expected QALYs and Costs for all transitions per cycle ####
    ### For Usual care
    a_Y_c_UC <- a_A_tunnels * a_R_c_UC
    a_Y_u_UC <- a_A_tunnels * a_R_u_UC
    
    ### For New treatment 1
    a_Y_c_Trt <- a_A_tunnels * a_R_c_Trt
    a_Y_u_Trt <- a_A_tunnels * a_R_u_Trt
    
    ### For New treatment 2
    a_Y_c_Trt2 <- a_A_tunnels_trt2 * a_R_c_Trt2
    a_Y_u_Trt2 <- a_A_tunnels_trt2 * a_R_u_Trt2
    
    #### Expected QALYs and Costs per cycle ####
    ## Vector of qalys under Usual care
    v_qaly_UC <- rowSums(t(colSums(a_Y_u_UC)))
    ## Vector of costs under Usual care
    v_cost_UC <- rowSums(t(colSums(a_Y_c_UC)))
    
    ## Vector of qalys under New treatment 1
    v_qaly_Trt <- rowSums(t(colSums(a_Y_u_Trt)))
    ## Vector of costs under New treatment 1
    v_cost_Trt <- rowSums(t(colSums(a_Y_c_Trt)))
    
    ## Vector of qalys under New treatment 2
    v_qaly_Trt2 <- rowSums(t(colSums(a_Y_u_Trt2)))
    ## Vector of costs under New treatment 2
    v_cost_Trt2 <- rowSums(t(colSums(a_Y_c_Trt2)))
    
    #### Discounted total expected QALYs and Costs per strategy ####
    ### For Usual care
    ## QALYs
    n_totqaly_UC <- sum(v_qaly_UC * v_dwe * v_hcc)
    ## Costs
    n_totcost_UC <- sum(v_cost_UC * v_dwc * v_hcc)
    
    ### For New treatment 1
    ## QALYs
    n_totqaly_Trt <- sum(v_qaly_Trt * v_dwe * v_hcc)
    ## Costs
    n_totcost_Trt <- sum(v_cost_Trt * v_dwc * v_hcc)
    
    ### For New treatment 2
    ## QALYs
    n_totqaly_Trt2 <- sum(v_qaly_Trt2 * v_dwe * v_hcc)
    ## Costs
    n_totcost_Trt2 <- sum(v_cost_Trt2 * v_dwc * v_hcc)
    
    ########################### Cost-effectiveness analysis #######################
    ### Vector of total costs for all strategies
    v_ted_cost <- c(n_totcost_UC, n_totcost_Trt, n_totcost_Trt2)
    ### Vector of effectiveness for all strategies
    v_ted_qaly <- c(n_totqaly_UC, n_totqaly_Trt, n_totqaly_Trt2)
    
    ### Calculate incremental cost-effectiveness ratios (ICERs)
    df_ce <- data.frame(Cost       = v_ted_cost,
                        Effect     = v_ted_qaly,
                        Strategies = v_names_str)
    return(df_ce)
  }
  )
}
