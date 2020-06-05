
## General setup
n_age_init <- 25 # age at baseline
n_age_max <- 640 # maximum age of follow up
n_t <- n_age_max - n_age_init # time horizon, number of cycles

n_tuns <- 2:(n_t*2-2)

## If you increase n_t, then n_tuns would be doubled, it will go over 8GB
## If you decrease n_t, then we are not at the max number of tunnels (n_tuns), we could do better
## We cannot decrease n_t and increase the number of tunnels infinitely, in other words, if we decrease n_t, n_tuns would have to be decreased as well (then we are not at max capability)
## We cannot have many combinations (decrease n_t in exchange for more n_tun)

### We can say, capacity at 8GB is at roughly 615 cycles, if we max the number of states it can have at each total number of cycles, in this case it's roughly 1231.

k <- length(n_tuns)

##################################### Model input #########################################

## Tunnel inputs
# Number of tunnels
n_tunnel_size <- n_tuns[k]
# Name for tunnels states of Sick state
v_Sick_tunnel <- paste("S1_", seq(1, n_tunnel_size), "Yr", sep = "")
# Create variables for model with tunnels
v_n_tunnels <- c("H", v_Sick_tunnel, "S2", "D") # health state names
n_states_tunnels <- length(v_n_tunnels)         # number of health statess



a_P_tunnels <- array(0, dim = c(n_states_tunnels, n_states_tunnels, n_t),
                     dimnames = list(v_n_tunnels, v_n_tunnels, 0:(n_t - 1)))


