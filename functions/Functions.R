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
  ### If it is a transition probability array
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
plot_trace <- function(l_m_M) {
  library(data.table)
  library(tidyr)
  
  ## Define colors and line types
  cols      <- c("H"  = "green",   "S1" = "orange", 
                 "S2" = "red",     "D"  = "blue")
  line_type <- c("H"  = "solid",   "S1" = "dashed", 
                 "S2" = "dotdash", "D"  = "dotted")
  
  # DARTH colors
  DARTHgreen      <- '#009999'  
  DARTHyellow     <- '#FDAD1E'  
  DARTHblue       <- '#006699' 
  DARTHlightgreen <- '#00adad'
  DARTHgray       <- '#666666'
  DARTHcols <- c("H"  = DARTHgreen,  "S1" = DARTHblue, 
                 "S2" = DARTHyellow, "D"  = DARTHgray)
  
  ## Pre-process the list of cohort traces for plotting
  # append column indicating which strategy the each trace was constructed under
  l_m_M <- lapply(seq_along(l_m_M), 
                  function(y, n, i) {cbind(as.data.frame(y[[i]]), 
                                           cycle = 0:n_t,
                                           strategy = rep(n[[i]], nrow(y[[i]]))
                  )}, y = l_m_M, n = names(l_m_M))
  # combine the list of traces into one
  m_M_plot <- do.call(rbind, l_m_M)
  m_M_plot <- gather(m_M_plot, key = state, value, H, S1, S2, D)
  # order the strategy and the states
  m_M_plot$state    <- factor(m_M_plot$state, levels = v_n)
  m_M_plot$strategy <- factor(m_M_plot$strategy, levels = v_names_str)
  
  ## Plot the cohort trace for scenarios Usual care and New treatment 1 
  p <- ggplot(m_M_plot, aes(x = cycle, y = value, 
                            color = state, linetype = state)) +
       geom_line(size = 1) +
       scale_colour_manual(name = "Health state", 
                           values = cols) +
       scale_linetype_manual(name = "Health state",
                             values = line_type) +
       xlab("Cycle") +
       ylab("Proportion of the cohort") +
       theme_bw(base_size = 14) +
       theme(legend.position  = "bottom", 
             legend.background = element_rect(fill = NA)) +
       facet_wrap(~strategy)
 
   return(p) 
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
plot_surv <- function(l_m_M) {

  ## Define colors and line types
  cols      <- c("orange", "blue", "red", "green")
  line_type <- c("solid", "dashed", "dotdash", "dotted")
  
  ## Calculate survival probabilities 
  surv <- as.data.frame(lapply(l_m_M, function(x) {rowSums(x[, !colnames(x) == "D"])}))
  colnames(surv) <- v_names_str
  surv$Cycle <- 0:n_t
  surv <- gather(surv, Strategy, Survival, 1:4)
  
  # order the strategy and the states
  surv$Strategy <- factor(surv$Strategy, levels = v_names_str)
  
  ## Plot the cohort trace for scenarios Usual care and New treatment 1 
  p <- ggplot(surv, 
              aes(x = Cycle, y = Survival, group = Strategy)) +
              geom_line(aes(linetype = Strategy, col = Strategy), size = 1.2) +
              xlab("Cycle") +
              ylab("Proportion alive") +
              scale_color_manual(values = as.character(cols)) +
              scale_linetype_manual(values = as.character(line_type)) +
              ggtitle("Survival probabilities") + 
              theme_bw(base_size = 14) +
              theme()
  
  return(p) 
}

#----------------------------------------------------------------------------#
####     Function to plot proportion of Sicker among sick individuals     ####
#----------------------------------------------------------------------------#
#' Plot proportion of Sicker among Sick individuals
#'
#' \code{plot_surv} plots the proportion of individuals in the Sicker state 
#' among those in the Sick or the Sicker state.
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a ggplot object - plot of the proportion of Sicker among sick
#' 
plot_sicker <- function(l_m_M) {
  
  ## Define colors and line types
  cols      <- c("orange", "blue", "red", "green")
  line_type <- c("solid", "dashed", "dotdash", "dotted")
  
  ## Calculate proportion of Sicker among the sick
  prop_S2 <- data.frame(Cycle = 1:n_t)
  for (i in 1:n_str) {
    m_M <- l_m_M[[i]]
    v_S_ad      <- rowSums(m_M[, !colnames(m_M) == "D"])  # vector with survival curve
    v_prev_S1S2 <- rowSums(m_M[, c("S1", "S2")]) / v_S_ad # vector with prevalence of Sick and Sicker
    v_prop_S2   <- m_M[-1, "S2"] / v_prev_S1S2[-1] 
    prop_S2 <- cbind(prop_S2, v_prop_S2)
  }
  colnames(prop_S2) <- c("Cycle", v_names_str)
  prop_S2 <- gather(prop_S2, Strategy, Proportion, 2:5)
  
  # order the strategy and the states
  prop_S2$Strategy <- factor(prop_S2$Strategy, levels = v_names_str)
  
  ## Plot the cohort trace for scenarios Usual care and New treatment 1 
  p <- ggplot(prop_S2, 
              aes(x = Cycle, y = Proportion, group = Strategy)) +
    geom_line(aes(linetype = Strategy, col = Strategy), size = 1.2) +
    xlab("Cycle") +
    ylab("Proportion in Sicker") +
    scale_color_manual(values = as.character(cols)) +
    scale_linetype_manual(values = as.character(line_type)) +
    ggtitle("Proportion of Sicker among sick individuals") + 
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
plot_prevalence <- function(l_m_M) {
  
  ## Define colors and line types
  cols      <- c("H"  = "green",   "S1" = "orange", 
                 "S2" = "red",     "D"  = "blue")
  line_type <- c("H"  = "solid",   "S1" = "dashed", 
                 "S2" = "dotdash", "D"  = "dotted")
  
  ## Calculate the prevalences
  for (i in 1:n_str) {
    m_M <- l_m_M[[i]]
    v_S_ad <- rowSums(m_M[, !colnames(m_M) == "D"])
    v_prev_S1   <- m_M[, "S1"] / v_S_ad # vector with prevalence of Sick
    v_prev_S2   <- m_M[, "S2"] / v_S_ad # vector with prevalence of Sicker
    v_prev_S1S2 <- rowSums(m_M[, c("S1", "S2")]) / v_S_ad # vector with prevalence of Sick and Sicker
    ## Data.frame with all prevalence
    df_prev_states <- data.frame(Cycle = 0:n_t, 
                                 States = ordered(rep(c("S1", "S2", "S1 and S2"),
                                                  each = (n_t + 1)), 
                                                  levels = c("S1", "S2", "S1 and S2")), 
                                 Prevalence = c(v_prev_S1, 
                                                v_prev_S2, 
                                                v_prev_S1S2))
  }
  
  colnames(prev) <- v_names_str
  
  
  surv$cycle <- 0:n_t
  surv <- gather(surv, strategy, survival, 1:4)
  
  # order the strategy and the states
  surv$strategy <- factor(surv$strategy, levels = v_names_str)
  
  ## Plot the cohort trace for scenarios Usual care and New treatment 1 
  p <- ggplot(surv, 
              aes(x = cycle, y = survival, group = strategy)) +
    geom_line(aes(linetype = strategy, col = strategy), size = 1.2) +
    xlab("Cycle") +
    ylab("Proportion alive") +
    scale_color_manual(values = as.character(cols)) +
    scale_linetype_manual(values = as.character(line_type)) +
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
