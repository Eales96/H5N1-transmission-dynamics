
# Function for running multiple simulations with the same initial conditions (returns dataframe of all simulations)

run_multi_sims_reactive <- function(df_fields = df_fields,
                           df_milk = df_milk,
                           params = params,
                           n_days = 100,
                           milking_frequency = 2,
                           N_sims = 10,
                           cleaning_time = NA,
                           reaction_threshold=10,
                           N_fields=5){
  
  
  list_of_results <- lapply(1:N_sims, function(i) {
    # Your simulation code here that returns a data frame for one run
    print(i)
    sim_temp <- simulation_reactive(df_fields = df_fields,
                           df_milk = df_milk,
                           params = params,
                           n_days = n_days,
                           milking_frequency = milking_frequency,
                           cleaning_time = cleaning_time,
                           reaction_threshold=reaction_threshold,
                           N_fields = N_fields)
    sim_temp$sim <- factor(i)
    return(sim_temp)
  })
  
  library(data.table)
  sims <- rbindlist(list_of_results)
  
  
  #sims_tot <- data.frame()
  #for(i in unique(sims$time)){
  #  for(j in unique(sims$sim)){
  #    df <- sims[sims$time==i &sims$sim==j,]
  #    row_df <- data.frame(time=i,
  #                         N_S = sum(df$N_S),
  #                         N_E = sum(df$N_E),
  #                         N_I = sum(df$N_I),
  #                         N_R = sum(df$N_R),
  #                         ID = "Total",
  #                         incidence = sum(df$incidence),
  #                         sim = j)
  #    
  #    sims_tot <- rbind(sims_tot, row_df)
  #  }
  #  
  #  
  #}
  
  #sims <- rbind(sims, sims_tot)
  return(sims)
}
