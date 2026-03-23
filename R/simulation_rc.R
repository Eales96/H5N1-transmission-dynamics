#Simulation for reactive cohorting OLD


simulation_rc <- function(df_fields, df_milk,  params, start_time=0, n_days=100, milking_time=6/24, milking_frequency=1){
  time <- start_time
  milking_time <- time+milking_time
  
  init_inc <- df_fields
  
  inc <- data.frame()
  for(i in 0:(n_days*milking_frequency) ){
    
    df_fields <- g_simulation_period(df_fields=df_fields,  params=params, time=time, next_milking_time = milking_time)
    time=milking_time
    
    
    ### Perform milking stuff here
    
    post_milking_outcomes <- milking_period(df_fields=df_fields, df_milk=df_milk,  params=params, time=time)
    df_fields <- post_milking_outcomes[[1]]
    df_milk <- post_milking_outcomes[[2]]
    time <- post_milking_outcomes[[3]]
    row <- data.frame(time=start_time+i/milking_frequency,
                      N_S = df_fields$N_S,
                      N_E = df_fields$N_E,
                      N_I = df_fields$N_I,
                      N_R = df_fields$N_R,
                      ID = df_fields$ID)
    inc <- rbind(inc, row)
    ###############################
    milking_time <- milking_time + 1/milking_frequency
    
  }
  
  inc$incidence <- 0
  for(i in unique(inc$ID)){
    inc[inc$ID == i,]$incidence <- c(init_inc[init_inc$ID == i,]$N_S[1] -inc[inc$ID == i,]$N_S[1] ,inc[inc$ID == i,]$N_S[1:(nrow(inc[inc$ID == i,])-1)] - inc[inc$ID == i,]$N_S[2:(nrow(inc[inc$ID == i,]))])
  }
  return(list(inc, df_fields, df_milk))
  
}
