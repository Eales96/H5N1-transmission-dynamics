
# Code for running a single simulation and returning a dataframe of incidence in cattle per unit time (1/milking frequency)
# and the event log detailing all infection events
simulation_detailed_infections_only <- function(df_fields, df_milk,  params, time=0, n_days=100, milking_time=6/24, milking_frequency=2, cleaning_time=NA){
  
  inc <- data.frame()
  
  event_log <- df_fields
  event_log$time <- time
  event_log$event <- "initial"
  
  for(i in 0:(n_days*milking_frequency) ){
    #print(i)
    g_events <- g_simulation_period_detailed_infections_only(df_fields=df_fields,  params=params, time=time, next_milking_time = milking_time)
    df_fields <- g_events[[1]]
    event_log_g <- g_events[[2]]
    time=milking_time
    
    
    ### Perform milking stuff here
    
    post_milking_outcomes <- milking_period_detailed_infections_only(df_fields=df_fields, df_milk=df_milk,  params=params, time=time, cleaning_time=cleaning_time)
    
    
    df_fields <- post_milking_outcomes[[1]]
    df_milk <- post_milking_outcomes[[2]]
    time <- post_milking_outcomes[[3]]
    event_log_m <- post_milking_outcomes[[4]]
    
    row <- data.frame(time=i/milking_frequency,
                      N_S = df_fields$N_S,
                      N_E = df_fields$N_E,
                      N_I = df_fields$N_I,
                      N_R = df_fields$N_R,
                      ID = df_fields$ID)
    
    inc <- rbind(inc, row)
    ###############################
    milking_time <- milking_time + 1/milking_frequency
    event_log <- rbind(event_log, event_log_g, event_log_m)
  }
  
  
  #inc$inc <- NA
  #inc[2: (nrow(inc)),]$inc <- inc[1: (nrow(inc)-1),]$N_S -inc[2: (nrow(inc)),]$N_S
  
  
  return(event_log)
  
}

