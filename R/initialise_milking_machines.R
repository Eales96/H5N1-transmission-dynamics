
# Function to initialise the state of milking machines

initialise_milking_machines <- function(N_milkers){
  
  data.frame(status = rep(0, N_milkers),
             time_inf = rep(NA, N_milkers),
             uses_since_inf = rep(NA, N_milkers),
             cleans_since_inf = rep(NA, N_milkers))
  
}