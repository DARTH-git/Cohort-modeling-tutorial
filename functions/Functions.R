#----------------------------------------------------------------------------#
####              Function to convert probabilities to rates              ####
#----------------------------------------------------------------------------#
#' Convert a probability to a rate
#'
#' \code{prob_to_rate} checks if a probability is between 0 and 1 and convert it to a rate.
#'
#' @param p probability
#' @param t time/ frequency
#' @return a number - converted rate
#' 
prob_to_rate <- function(p, t = 1){
  if (sum(p > 1) >0 | sum(p < 0) > 0 ){
    print("probability not between 0 and 1")
  }
  r <- -log(1-p)/ t
  return(r)
}

#----------------------------------------------------------------------------#
####              Function to convert rates to probabilities              ####
#----------------------------------------------------------------------------#
#' Convert a rate to a probability
#'
#' \code{rate_to_prob} convert a rate to a probability.
#'
#' @param r rate
#' @param t time/ frequency
#' @return a number - converted probability
#' 
# Function to convert rates to probabilities
rate_to_prob <- function(r, t = 1){
  p <- 1 - exp(- r * t)
  return(p)
}

#----------------------------------------------------------------------------#
#### Function to check if transition probability matrix or array is valid ####
#----------------------------------------------------------------------------#
#' Check if transition array is valid
#'
#' \code{check_transition_probability} checks if transition probabilities are in \[0, 1\].
#'
#' @param a_P A transition probability array.
#' @param array Logical variable to specify if it is a transition probability array instead of a matrix. Default = TRUE
#' @param err_stop Logical variable to stop model run if set up as TRUE. Default = FALSE.
#' @param verbose Logical variable to indicate print out of messages. Default = FALSE
#' @return
#' This function stops if transition probability array/matrix is not valid and shows 
#' what are the entries that are not valid
#' @import utils
#' @export
check_transition_probability <- function(a_P, 
                                         array = TRUE, 
                                         err_stop = FALSE, 
                                         verbose = FALSE) {
  ### If it is a transition probability array   # 1) is_array or check dim
                                                # 2) if not an array make it an array (don't repeat code)
  if (array) {
    m_indices_notvalid <- arrayInd(which(a_P < 0 | a_P > 1), 
                                   dim(a_P))
    if(dim(m_indices_notvalid)[1] != 0){
      v_rows_notval   <- rownames(a_P)[m_indices_notvalid[, 1]]
      v_cols_notval   <- colnames(a_P)[m_indices_notvalid[, 2]]
      v_cycles_notval <- dimnames(a_P)[[3]][m_indices_notvalid[, 3]]
      
      df_notvalid <- data.frame(`Transition probabilities not valid:` = 
                                  matrix(paste0(paste(v_rows_notval, v_cols_notval, sep = "->"),
                                                "; at cycle ",
                                                v_cycles_notval), ncol = 1), 
                                check.names = FALSE)
      if(err_stop) {
        stop("Not valid transition probabilities\n",
             paste(capture.output(df_notvalid), collapse = "\n"))
      }
      if(verbose){
        warning("Not valid transition probabilities\n",
                paste(capture.output(df_notvalid), collapse = "\n"))
      } 
    }
    
    ### If it is a transition probability matrix
  } else { 
    a_P <- as.array(a_P)
    m_indices_notvalid <- arrayInd(which(a_P < 0 | a_P > 1), 
                                   dim(a_P))
    if(dim(m_indices_notvalid)[1] != 0){
      v_rows_notval   <- rownames(a_P)[m_indices_notvalid[, 1]]
      v_cols_notval   <- colnames(a_P)[m_indices_notvalid[, 2]]
      
      df_notvalid <- data.frame(`Transition probabilities not valid:` = 
                                  matrix(paste0(paste(v_rows_notval, v_cols_notval, sep = "->")),
                                         ncol = 1), 
                                check.names = FALSE)
      
      if(err_stop) {
        stop("Not valid transition probabilities\n",
             paste(capture.output(df_notvalid), collapse = "\n"))
      }
      
      if(verbose){
        warning("Not valid transition probabilities\n",
                paste(capture.output(df_notvalid), collapse = "\n"))
      } 
    }
  }
}

#----------------------------------------------------------------------------#
####   Function to check if sum of transition probabilities equal to one  ####
#----------------------------------------------------------------------------#
#' Check if the sum of transition probabilities equal to one. 
#'
#' \code{check_sum_of_transition_array} checks if each of the rows of the 
#' transition probability array/matrix sum to one. 
#' 
#' @param a_P A transition probability array
#' @param array Logical variable to specify if it is a transition probability array instead of a matrix. Default = TRUE 
#' @param n_states Number of health states
#' @param n_t Number of cycles
#' @param err_stop Logical variable to stop model run if set up as TRUE. Default = FALSE
#' @param verbose Logical variable to indicate print out of messages. Default = FALSE
#' @return 
#' The transition probability array/matrix
#' @import dplyr
#' @export
check_sum_of_transition_array <- function(a_P,
                                          array = TRUE, 
                                          n_states,
                                          n_t,  
                                          err_stop = FALSE, 
                                          verbose = FALSE) {
  ### If it is a transition probability array
  if (array) {
    valid <- (apply(a_P, 3, function(x) sum(rowSums(x))) == n_states)
    if (!isTRUE(all.equal(as.numeric(sum(valid)), as.numeric(n_t)))) {
      if(err_stop) {
        stop("This is not a valid transition probability array")
      }
      
      if(verbose){
        warning("This is not a valid transition probability array")
      } 
    }
  ### If it is a transition probability matrix  
  } else {
    a_P <- as.array(a_P)
    valid <- (apply(a_P, 1, function(x) sum(x) == 1))
    if (!isTRUE(all.equal(as.numeric(sum(valid)), n_states))) {
      if(err_stop) {
        stop("This is not a valid transition probability matrix")
      }
      
      if(verbose){
        warning("This is not a valid transition probability matrix")
      } 
    }
  }
} 

#----------------------------------------------------------------------------#
####                    Function to plot cohort trace                     ####
#----------------------------------------------------------------------------#
#' Plot cohort trace
#'
#' \code{plot_trace} plots the cohort trace.
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a ggplot object - plot of the cohort trace
#' 
get_DARTH_cols <- function() {
  # DARTH colors
  DARTHgreen      <- '#009999'  
  DARTHyellow     <- '#FDAD1E'  
  DARTHblue       <- '#006699' 
  DARTHlightgreen <- '#00adad'
  DARTHgray       <- '#666666'
  DARTHcols <- c("H"  = DARTHgreen,  "S1" = DARTHblue, 
                 "S2" = DARTHyellow, "D"  = DARTHgray)
  
  return(DARTHcols) 
}

#----------------------------------------------------------------------------#
####                    Function to plot cohort trace                     ####
#----------------------------------------------------------------------------#
#' Plot cohort trace
#'
#' \code{plot_trace} plots the cohort trace.
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a ggplot object - plot of the cohort trace
#' 
plot_trace <- function(m_M) {
  df_M <- data.frame(Cycle = 0:n_t, m_M, check.names = F)
  df_M_long <- tidyr::gather(df_M, key = `Health State`, value, 2:ncol(df_M))

  ## Plot the cohort trace for scenarios Usual care and New treatment 1 
  p <- ggplot(df_M_long, aes(x = Cycle, y = value, 
                            color = `Health State`, linetype = `Health State`)) +
       geom_line(size = 1) +
       xlab("Cycle") +
       ylab("Proportion of the cohort") +
       theme_bw(base_size = 14) +
       theme(legend.position  = "bottom", 
             legend.background = element_rect(fill = NA)) 
       # just one trace udner one strategy, 
 
  return(p) 
}

#----------------------------------------------------------------------------#
####                    Function to plot cohort trace per strategy        ####
#----------------------------------------------------------------------------#
#' Plot cohort trace
#'
#' \code{plot_trace} plots the cohort trace.
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a ggplot object - plot of the cohort trace
#' 
plot_trace_strategy <- function(l_m_M) {
  ## Pre-process the list of cohort traces for plotting
  n_str <- length(l_m_M) 
  l_df_M <- lapply(l_m_M, as.data.frame)
  df_M_strategies <- data.table::rbindlist(l_df_M, use.names = T, 
                                           idcol = "Strategy")
  df_M_strategies$Cycle <- rep(0:n_t, n_str)
  
  m_M_plot <- tidyr::gather(df_M_strategies, key = `Health State`, value, 
                            2:(ncol(df_M_strategies)-1))
  # order the strategy and the states
  m_M_plot$`Health State`    <- factor(m_M_plot$`Health State`, levels = v_n)
  m_M_plot$Strategy <- factor(m_M_plot$Strategy, levels = v_names_str)
  
  ## Plot the cohort trace for scenarios Usual care and New treatment 1 
  p <- ggplot(m_M_plot, aes(x = Cycle, y = value, 
                            color = Strategy, linetype = Strategy)) +
    geom_line(size = 1) +
    xlab("Cycle") +
    ylab("Proportion of the cohort") +
    theme_bw(base_size = 14) +
    theme(legend.position  = "bottom", 
          legend.background = element_rect(fill = NA)) +
    facet_wrap(~ state)
  
  return(p) 
}

#----------------------------------------------------------------------------#
####                   Function to calculate survival curve                    ####
#----------------------------------------------------------------------------#
#' Plot survival curve
#'
#' \code{calc_surv} plots the survival probability curve.
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a ggplot object - plot of the survival curve
#' 
calc_surv <- function(l_m_M, v_names_death_states) {
  df_surv <- as.data.frame(lapply(l_m_M, 
                                  function(x) {
                                    rowSums(x[, !colnames(x) %in% v_names_death_states])
                                    }
                                  ))
  colnames(df_surv) <- v_names_str
  df_surv$Cycle     <- 0:n_t
  df_surv_long      <- tidyr::gather(df_surv, key = Strategy, Survival, 1:n_str)
  # order the strategy and the states
  df_surv_long$Strategy <- ordered(df_surv_long$Strategy, levels = v_names_str)
  df_surv_long <- df_surv_long %>% 
    select(Strategy, Cycle, Survival)
  
  return(df_surv_long)
}

#----------------------------------------------------------------------------#
####                   Function to calculate sick                    ####
#----------------------------------------------------------------------------#
#' Plot survival curve
#'
#' \code{calc_surv} plots the survival probability curve.
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a ggplot object - plot of the survival curve
#' 
calc_sick <- function(l_m_M, v_names_sick_states) {
  df_sick <- as.data.frame(lapply(l_m_M, 
                                  function(x) {
                                    rowSums(x[, colnames(x) %in% v_names_sick_states])
                                  }
  ))
  colnames(df_sick) <- v_names_str
  df_sick$Cycle     <- 0:n_t
  df_sick_long      <- tidyr::gather(df_sick, key = Strategy, Sick, 1:n_str)
  # order the strategy and the states
  df_sick_long$Strategy <- ordered(df_surv_long$Strategy, levels = v_names_str)
  df_sick_long <- df_sick_long %>% 
    select(Strategy, Cycle, Sick)
  
  return(df_sick_long)
}

#----------------------------------------------------------------------------#
####                   Function to calculate prevalence curve             ####
#----------------------------------------------------------------------------#
#' Calculate prevalence curve
#'
#' \code{plot_prevalence} plots the prevalence curve for Sick/Sicker states.
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a ggplot object - plot of the prevalence curve
#' 
calc_prevalence <- function(l_m_M, v_names_sick_states, v_names_dead_states) {
  df_alive      <- calc_surv(l_m_M, v_names_dead_states)
  df_prop_sick  <- calc_sick(l_m_M, v_names_sick_states)
  df_prevalence <- data.frame(Strategy   = df_alive$Strategy, 
                              Cycle      = df_alive$Cycle,
                              Prevalence = df_prop_sick$Sick / df_alive$Survival)
  return(df_prevalence) 
}

#----------------------------------------------------------------------------#
####                   Function to calculate prevalence curve                   ####
#----------------------------------------------------------------------------#
#' Calculate prevalence curve
#'
#' \code{plot_prevalence} plots the prevalence curve for Sick/Sicker states.
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a ggplot object - plot of the prevalence curve
#' 
calc_prop_sicker <- function(l_m_M, v_names_sick_states, v_names_sicker_states) {
  df_prop_sick   <- calc_sick(l_m_M, v_names_sick_states)
  df_prop_sicker <- calc_sick(l_m_M, v_names_sicker_states)
  df_prop_sick_sicker <- data.frame(Strategy   = df_M_strategies$Strategy, 
                                    Cycle      = df_M_strategies$Cycle,
                                    `Proporiton Sicker` = 
                                      df_prop_sicker$Sick / (df_prop_sick$Sick + df_prop_sicker$Sick))
  
  return(df_prop_sick_sicker) 
}

#----------------------------------------------------------------------------#
####                   Function to plot survival curve                    ####
#----------------------------------------------------------------------------#
#' Plot survival curve
#'
#' \code{plot_surv} plots the survival probability curve.
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a ggplot object - plot of the survival curve
#' 
plot_surv <- function(l_m_M, v_names_death_states) {
  # Calculate survival probabilities
  df_surv <- calc_surv(l_m_M, v_names_death_states)
  # Order the strategy and the states
  df_surv$Strategy <- factor(df_surv$Strategy, levels = v_names_str)
  
  ## Plot the cohort trace for scenarios Usual care and New treatment 1 
  p <- ggplot(df_surv, 
              aes(x = Cycle, y = Survival, group = Strategy)) +
    geom_line(aes(linetype = Strategy, col = Strategy), size = 1.2) +
    xlab("Cycle") +
    ylab("Proportion alive") +
    ggtitle("Survival probabilities") + 
    theme_bw(base_size = 14) +
    theme()
  
  return(p) 
}

#----------------------------------------------------------------------------#
####                   Function to plot prevalence curve                   ####
#----------------------------------------------------------------------------#
#' Plot prevalence curve
#'
#' \code{plot_prevalence} plots the prevalence curve for Sick/Sicker states.
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a ggplot object - plot of the prevalence curve
#' 
plot_prevalence <- function(l_m_M, v_names_sick_states, v_names_dead_states) {
  # Calculate prevalence
  df_prevalence <- calc_prevalence(l_m_M, v_names_sick_states, v_names_dead_states)
  # order the strategy and the states
  df_prevalence$Strategy <- factor(df_prevalence$Strategy, levels = v_names_str)
  
  ## Plot the cohort trace for scenarios Usual care and New treatment 1 
  p <- ggplot(df_prevalence, 
              aes(x = Cycle, y = Prevalence, group = Strategy)) +
    geom_line(aes(linetype = Strategy, col = Strategy), size = 1.2) +
    xlab("Cycle") +
    ylab("Prevalence") +
    ggtitle("Prevalence") + 
    theme_bw(base_size = 14) +
    theme()
  
  return(p) 
}

#----------------------------------------------------------------------------#
####                   Function to plot proportion Sicker curve                   ####
#----------------------------------------------------------------------------#
#' Plot proportion Sicker curve
#'
#' \code{plot_prevalence} plots the prevalence curve for Sick/Sicker states.
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a ggplot object - plot of the prevalence curve
#' 
plot_proportion_sicker <- function(l_m_M, v_names_sick_states, v_names_dead_states) {
  # order the strategy and the states
  df_prevalence <- calc_prop_sicker(l_m_M, v_names_sick_states, v_names_sicker_states)
  df_prevalence$Strategy <- factor(df_prevalence$Strategy, levels = v_names_str)
  
  ## Plot the cohort trace for scenarios Usual care and New treatment 1 
  p <- ggplot(df_prevalence, 
              aes(x = Cycle, y = Prevalence, group = Strategy)) +
    geom_line(aes(linetype = Strategy, col = Strategy), size = 1.2) +
    xlab("Cycle") +
    ylab("Proportion sicker") +
    ggtitle("Proportion of sicker among sick") + 
    theme_bw(base_size = 14) +
    theme()
  
  return(p) 
}

#----------------------------------------------------------------------------#
####                     Function to format CEA table                     ####
#----------------------------------------------------------------------------#
#' Format CEA table
#'
#' \code{format_table_cea} formats the CEA table.
#'
#' @param table_cea a dataframe object - table with CEA results
#' @return a dataframe object - formatted CEA table
#' 
format_table_cea <- function(table_cea) {
  ## Format column names
  colnames(table_cea)[colnames(table_cea) %in% c("Cost", "Effect", "Inc_Cost", "Inc_Effect", "ICER")] <- 
    c("Costs ($)", "QALYs", "Incremental Costs ($)", "Incremental QALYs", "ICER ($/QALY)") 
  
  ## Format rows
  table_cea$`Costs ($)` <- comma(round(table_cea$`Costs ($)`, 0))
  table_cea$`Incremental Costs ($)` <- comma(round(table_cea$`Incremental Costs ($)`, 0))
  table_cea$QALYs <- round(table_cea$QALYs, 2)
  table_cea$`Incremental QALYs` <- round(table_cea$`Incremental QALYs`, 2)
  table_cea$`ICER ($/QALY)` <- comma(round(table_cea$`ICER ($/QALY)`, 0))
  return(table_cea)
}
