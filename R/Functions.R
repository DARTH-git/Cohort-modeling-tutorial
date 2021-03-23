#----------------------------------------------------------------------------#
####                    Function to plot cohort trace                     ####
#----------------------------------------------------------------------------#
#' Plot cohort trace
#'
#' \code{plot_trace} plots the cohort trace.
#'
#' @param m_M a cohort trace matrix
#' @return a ggplot object - plot of the cohort trace
#' 
plot_trace <- function(m_M) {
  df_M      <- data.frame(Cycle = 0:n_cycles, m_M, check.names = F)
  df_M_long <- tidyr::gather(df_M, key = `Health State`, value, 2:ncol(df_M))
  df_M_long$`Health State` <- factor(df_M_long$`Health State`, levels = v_names_states)
  gg_trace <- ggplot(df_M_long, aes(x = Cycle, y = value, 
                            color = `Health State`, linetype = `Health State`)) +
    geom_line(size = 1) +
    xlab("Cycle") +
    ylab("Proportion of the cohort") +
    scale_x_continuous(breaks = number_ticks(8)) + 
    theme_bw(base_size = 14) +
    theme(legend.position  = "bottom", 
          legend.background = element_rect(fill = NA)) 

  return(gg_trace) 
}

#----------------------------------------------------------------------------#
####                    Function to plot cohort trace per strategy        ####
#----------------------------------------------------------------------------#
#' Plot cohort trace per strategy
#'
#' \code{plot_trace} plots the cohort trace for each strategy, split by health state.
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a ggplot object - plot of the cohort trace for each strategy split by health state.
#' 
plot_trace_strategy <- function(l_m_M) {
  n_str <- length(l_m_M) 
  l_df_M <- lapply(l_m_M, as.data.frame)
  df_M_strategies <- data.table::rbindlist(l_df_M, use.names = T, 
                                           idcol = "Strategy")
  df_M_strategies$Cycle <- rep(0:n_cycles, n_str)
  m_M_plot <- tidyr::gather(df_M_strategies, key = `Health State`, value, 
                            2:(ncol(df_M_strategies)-1))
  m_M_plot$`Health State`    <- factor(m_M_plot$`Health State`, levels = v_names_states)
  m_M_plot$Strategy <- factor(m_M_plot$Strategy, levels = v_names_str)

  p <- ggplot(m_M_plot, aes(x = Cycle, y = value, 
                            color = Strategy, linetype = Strategy)) +
    geom_line(size = 1) +
    scale_color_brewer(palette="RdBu") +
    xlab("Cycle") +
    ylab("Proportion of the cohort") +
    theme_bw(base_size = 14) +
    theme(legend.position  = "bottom", 
          legend.background = element_rect(fill = NA)) +
    facet_wrap(~ `Health State`)
  
  return(p) 
}

#----------------------------------------------------------------------------#
####             Function to calculate survival probabilities             ####
#----------------------------------------------------------------------------#
#' Calculate survival probabilities
#'
#' \code{calc_surv} calculates the survival probabilities.
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a dataframe containing survival probabilities for each strategy
#' 
calc_surv <- function(l_m_M, v_names_death_states) {
  df_surv <- as.data.frame(lapply(l_m_M, 
                                  function(x) {
                                    rowSums(x[, !colnames(x) %in% v_names_death_states])
                                    }
                                  ))
  colnames(df_surv) <- v_names_str
  df_surv$Cycle     <- 0:n_cycles
  df_surv_long      <- tidyr::gather(df_surv, key = Strategy, Survival, 1:n_str)
  df_surv_long$Strategy <- ordered(df_surv_long$Strategy, levels = v_names_str)
  df_surv_long <- df_surv_long %>% 
    select(Strategy, Cycle, Survival)
  
  return(df_surv_long)
}

#----------------------------------------------------------------------------#
####                Function to calculate state proportions               ####
#----------------------------------------------------------------------------#
#' Calculate state proportions
#'
#' \code{calc_surv} calculates the proportions of the cohort in specified states
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a dataframe containing proportions in specified states for each strategy
#' 
calc_sick <- function(l_m_M, v_names_sick_states) {
  n_sick_states <- length(v_names_sick_states)
  df_sick <- as.data.frame(lapply(l_m_M,
                                  function(x) {
                                    if (n_sick_states == 1) {
                                      x[, colnames(x) %in% v_names_sick_states]
                                    } else {
                                      rowSums(x[, colnames(x) %in% v_names_sick_states])
                                    }
                                  }
  ))
  colnames(df_sick) <- v_names_str
  df_sick$Cycle     <- 0:n_cycles
  df_sick_long      <- tidyr::gather(df_sick, key = Strategy, Sick, 1:n_str)
  df_sick_long$Strategy <- ordered(df_sick_long$Strategy, levels = v_names_str)
  df_sick_long <- df_sick_long %>% 
    select(Strategy, Cycle, Sick)
  
  return(df_sick_long)
}

#----------------------------------------------------------------------------#
####                   Function to calculate prevalence                   ####
#----------------------------------------------------------------------------#
#' Calculate prevalence
#'
#' \code{plot_prevalence} calculate the prevalence for different health states.
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a dataframe containing prevalence of specified health states for each strategy
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
####           Function to calculate state-in-state proportions           ####
#----------------------------------------------------------------------------#
#' Calculate state-in-state proportions
#'
#' \code{plot_prevalence} calculates the proportion of a speciefied subset of states among a set of specified states
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a dataframe containing state-in-state proportions of specified health states for each strategy
#' 
calc_prop_sicker <- function(l_m_M, v_names_sick_states, v_names_sicker_states) {
  df_prop_sick   <- calc_sick(l_m_M, v_names_sick_states)
  df_prop_sicker <- calc_sick(l_m_M, v_names_sicker_states)
  df_prop_sick_sicker <- data.frame(Strategy   = df_prop_sick$Strategy, 
                                    Cycle      = df_prop_sick$Cycle,
                                    `Proportion Sicker` = 
                                     df_prop_sicker$Sick / 
                                    (df_prop_sick$Sick + df_prop_sicker$Sick))
  
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
  df_surv <- calc_surv(l_m_M, v_names_death_states)
  df_surv$Strategy <- factor(df_surv$Strategy, levels = v_names_str)
  df_surv$Survival <- round(df_surv$Survival, 2)
  
  p <- ggplot(df_surv, 
              aes(x = Cycle, y = Survival, group = Strategy)) +
    geom_line(aes(linetype = Strategy, col = Strategy), size = 1.2) +
    scale_color_brewer(palette="RdBu") +
    xlab("Cycle") +
    ylab("Proportion") +
    ggtitle("Survival probabilities") +
    theme_bw(base_size = 14) +
    theme()
  
  return(p) 
}

#----------------------------------------------------------------------------#
####                   Function to plot prevalence curve                  ####
#----------------------------------------------------------------------------#
#' Plot prevalence curve
#'
#' \code{plot_prevalence} plots the prevalence curve for specified health states.
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a ggplot object - plot of the prevalence curve
#' 
plot_prevalence <- function(l_m_M, v_names_sick_states, v_names_dead_states) {
  df_prevalence <- calc_prevalence(l_m_M, v_names_sick_states, v_names_dead_states)
  df_prevalence$Strategy <- factor(df_prevalence$Strategy, levels = v_names_str)
  df_prevalence$Proportion.Sicker <- round(df_prevalence$Prevalence, 2)
  
  p <- ggplot(df_prevalence, 
              aes(x = Cycle, y = Prevalence, group = Strategy)) +
    geom_line(aes(linetype = Strategy, col = Strategy), size = 1.2) +
    scale_color_brewer(palette = "RdBu") +
    xlab("Cycle") +
    ylab("Proportion") + 
    ggtitle(paste("Prevalence", "of", paste(v_names_sick_states, collapse = " "))) + 
    theme_bw(base_size = 14) +
    theme()
  
  return(p) 
}

#----------------------------------------------------------------------------#
####           Function to plot state-in-state proportion curve           ####
#----------------------------------------------------------------------------#
#' Plot state-in-state proportion curve
#'
#' \code{plot_prevalence} plots the 
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a ggplot object - plot of state-in-state proportion curve
#' 
plot_proportion_sicker <- function(l_m_M, v_names_sick_states, v_names_sicker_states) {
  df_proportion_sicker <- calc_prop_sicker(l_m_M, v_names_sick_states, v_names_sicker_states)
  df_proportion_sicker$Strategy <- factor(df_proportion_sicker$Strategy, levels = v_names_str)
  df_proportion_sicker$Proportion.Sicker <- round(df_proportion_sicker$Proportion.Sicker, 2)
  
  p <- ggplot(df_proportion_sicker, 
              aes(x = Cycle, y = Proportion.Sicker, group = Strategy)) +
    geom_line(aes(linetype = Strategy, col = Strategy), size = 1.2, na.rm = T) +
    scale_color_brewer(palette = "RdBu") +
    xlab("Cycle") +
    ylab("Proportion") +
    ggtitle(paste(paste("Proportion of", v_names_sicker_states), 
                  paste(c("among", v_names_sick_states), collapse = " "))) + 
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
  colnames(table_cea)[colnames(table_cea) 
                      %in% c("Cost", 
                             "Effect", 
                             "Inc_Cost", 
                             "Inc_Effect",
                             "ICER")] <- 
                                               
                           c("Costs ($)", 
                             "QALYs", 
                             "Incremental Costs ($)", 
                             "Incremental QALYs", 
                             "ICER ($/QALY)") 
  
  table_cea$`Costs ($)` <- comma(round(table_cea$`Costs ($)`, 0))
  table_cea$`Incremental Costs ($)` <- comma(round(table_cea$`Incremental Costs ($)`, 0))
  table_cea$QALYs <- round(table_cea$QALYs, 2)
  table_cea$`Incremental QALYs` <- round(table_cea$`Incremental QALYs`, 2)
  table_cea$`ICER ($/QALY)` <- comma(round(table_cea$`ICER ($/QALY)`, 0))
  return(table_cea)
}
