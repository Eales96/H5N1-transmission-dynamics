
# Function to initalise susceptible cows in fields

initialise_fields <- function(N_cows,
                              Milking_station = rep(1, length(N_cows)),
                              Area = rep(1, length(N_cows))){
  
  data.frame(N_S = N_cows,
             N_E = rep(0, length(N_cows)),
             N_I = rep(0, length(N_cows)),
             N_R = rep(0, length(N_cows)),
             ID = seq(1,length(N_cows)),
             M = Milking_station,
             A = Area)
}