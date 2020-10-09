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
    ## Age-specific transition probabilities to the Dead state
    # compute mortality rates
    v_r_S1Dage  <- v_r_HDage * hr_S1        # Age-specific mortality rate in the Sick state 
    v_r_S2Dage  <- v_r_HDage * hr_S2        # Age-specific mortality rate in the Sicker state 
    # transform rates to probabilities
    v_p_S1Dage  <- rate_to_prob(v_r_S1Dage) # Age-specific mortality risk in the Sick state
    v_p_S2Dage  <- rate_to_prob(v_r_S2Dage) # Age-specific mortality risk in the Sicker state
    
    ## Transition probability of becoming Sicker when Sick for B
    # transform odds ratios to probabilites 
    logit_S1S2  <- boot::logit(p_S1S2)        # log-odds of becoming Sicker when Sick
    p_S1S2_trtB <- boot::inv.logit(logit_S1S2 +
                                     lor_S1S2)  # probability to become Sicker when Sick 
    # under B conditional on surviving
    
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
                     dimnames = list(0:n_t, v_n))
    # Store the initial state vector in the first row of the cohort trace
    m_M_ad[1, ] <- v_s_init
    ## Initialize cohort trace for srategies B and AB
    m_M_ad_strB <- m_M_ad # structure and initial states remain the same.
    
    ## Initialize transition array which will capture transitions from each state to another over time 
    # for srategies SoC and A
    a_A <- array(0,
                 dim      = c(n_states, n_states, n_t + 1),
                 dimnames = list(v_n, v_n, 0:n_t))
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
    
    ## Store the transition array for each strategy in a list
    l_a_A <- list(a_A,
                  a_A,
                  a_A_strB,
                  a_A_strB)
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
                         dimnames = list(v_n, v_n, 0:n_t))
      # Expand the transition matrix of state costs across cycles to form a transition array of state costs
      a_R_c_str <- array(m_c_str, 
                         dim      = c(n_states, n_states, n_t + 1),
                         dimnames = list(v_n, v_n, 0:n_t))
      
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
    df_cea <- df_cea[match(v_names_str, df_cea$Strategy), ]
    return(df_cea)
  }
  )
}
