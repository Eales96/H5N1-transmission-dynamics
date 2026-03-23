
# Function for running multiple simulations with the same initial conditions (returns dataframe of all simulations)

run_multi_sims <- function(df_fields = df_fields,
                                    df_milk = df_milk,
                                    params = params,
                                    n_days = 100,
                                    milking_frequency = 2,
                                    N_sims = 10,
                                    cleaning_time = NA,
                                    N_fields=5){
  
  
  list_of_results <- lapply(1:N_sims, function(i) {
    # Your simulation code here that returns a data frame for one run
    print(i)
    sim_temp <- simulation(df_fields = df_fields,
                                    df_milk = df_milk,
                                    params = params,
                                    n_days = n_days,
                                    milking_frequency = milking_frequency,
                                    cleaning_time = cleaning_time,
                                    N_fields = N_fields)
    sim_temp$sim <- factor(i)
    return(sim_temp)
  })
  
  library(data.table)
  sims <- rbindlist(list_of_results)
  
  return(sims)
}

