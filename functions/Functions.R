#---------------------------------------------------#
####         Function to plot cohort trace       ####
#---------------------------------------------------#
#' Plot cohort trace
#'
#' \code{plot_trace} plots the cohort trace.
#'
#' @param trace cohort trace matrix
#' @return a ggplot object - plot of the cohort trace
#' 
plot_trace <- function(trace) {
  ## Define colors and line types
  cols <- c("H"  = "green", "S1" = "orange", 
            "S2" = "red", "D" = "blue")
  
  # DARTH colors
  DARTHgreen      <- '#009999'  
  DARTHyellow     <- '#FDAD1E'  
  DARTHblue       <- '#006699' 
  DARTHlightgreen <- '#00adad'
  DARTHgray       <- '#666666'
  DARTHcols <- c("H"  = DARTHgreen, "S1" = DARTHblue, 
                 "S2" = DARTHyellow, "D" = DARTHgray)
  
  line_type <-  c("H" = "solid", "S1" = "dashed", "S2" = "dotdash", "D" = "dotted")
  
  ## Plot the cohort trace for scenarios Usual care and New treatment 1 
  p <- ggplot(melt(trace), aes(x = Var1, y = value, 
                               color = Var2, linetype = Var2)) +
       geom_line(size = 1) +
       scale_colour_manual(name = "Health state", 
                           values = cols) +
       scale_linetype_manual(name = "Health state",
                             values = line_type) +
       xlab("Cycle") +
       ylab("Proportion of the cohort") +
       theme_bw(base_size = 14) +
       theme(legend.position = "bottom", 
             legend.background = element_rect(fill = NA))
  
  return(p) 
}

#---------------------------------------------------#
####         Function to format CEA table        ####
#---------------------------------------------------#
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

#------------------------------------------------------------------------------------------#
####         Function to check if transition probability matrix or array is valid       ####
#------------------------------------------------------------------------------------------#
#' Check if transition probabilty matrix or array is valid (rows sum to 1)
#'
#' \code{check_transition} checks if the transition probability matrix or array is valid.
#'
#' @param transition_object transition probability matrix or array
#' @param array indicates whether your transition object is a matrix (= F) or an array (= T)
#' @return a ggplot object - plot of the cohort trace
#' 
check_transition <- function(transition_object, array = F) {
  if (array) {
    valid <- apply(transition_object, 3, function(x) sum(rowSums(x))==n_states)
    if (!isTRUE(all.equal(as.numeric(sum(valid)), as.numeric(n_t)))) {
        stop("This is not a valid transition array")
    }
  }
  else {
    valid <- rowSums(transition_object) # sum the rows 
    if (!isTRUE(all.equal(as.numeric((valid)), as.numeric(rep(1, n_states))))) { #check if the rows are all equal to one 
        stop("This is not a valid transition matrix")
    }
  }
}

#----------------------------------------------------------------#
####         Function to convert probabilities to rates       ####
#----------------------------------------------------------------#
#' Convert a probability to a rate
#'
#' \code{ProbRate} checks if a probability is between 0 and 1 and convert it to a rate.
#'
#' @param p probability
#' @param t time/ frequency
#' @return a number - converted rate
#' 
ProbRate <- function(p, t = 1){
  if (sum(p > 1) >0 | sum(p < 0) > 0 ){
    print("probability not between 0 and 1")
  }
  r <- -log(1-p)/ t
  return(r)
}

#----------------------------------------------------------------#
####         Function to convert rates to probabilities       ####
#----------------------------------------------------------------#
#' Convert a rate to a probability
#'
#' \code{ProbRate} convert a rate to a probability.
#'
#' @param r rate
#' @param t time/ frequency
#' @return a number - converted probability
#' 
# Function to convert rates to probabilities
RateProb <- function(r, t = 1){
  p <- 1 - exp(- r * t)
  return(p)
}

#-------------------------------------------------------------------------------------------------#
####         Function to convert probabilities and probabilities - different frequency         ####
#-------------------------------------------------------------------------------------------------#
#' Convert a probability to a probability with a different frequency
#'
#' \code{ProbRate} convert a probability to a probability with a different frequency.
#'
#' @param p probability
#' @param t time/ frequency
#' @return a number - converted probability
#' 
ProbProb <- function(p, t = 1){
  p_new <- RateProb(ProbRate(p, t))
  return(p_new)
}

