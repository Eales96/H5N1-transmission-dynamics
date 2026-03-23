# Code for running a single simulation and returning a dataframe of incidence in cattle per unit time (1/milking frequency)

simulation_old <- function(df_fields, df_milk,  params, time=0, n_days=100, milking_time=6/24, milking_frequency=2, cleaning_time=NA){
  
  inc <- data.frame()
  for(i in 0:(n_days*milking_frequency) ){
    
    df_fields <- g_simulation_period(df_fields=df_fields,  params=params, time=time, next_milking_time = milking_time)
    time=milking_time
    
    
    ### Perform milking stuff here
    
    post_milking_outcomes <- milking_period(df_fields=df_fields, df_milk=df_milk,  params=params, time=time, cleaning_time=cleaning_time)
    df_fields <- post_milking_outcomes[[1]]
    df_milk <- post_milking_outcomes[[2]]
    time <- post_milking_outcomes[[3]]
    row <- data.frame(time=i/milking_frequency,
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
    inc[inc$ID == i,]$incidence <- c(0,inc[inc$ID == i,]$N_S[1:(nrow(inc[inc$ID == i,])-1)] - inc[inc$ID == i,]$N_S[2:(nrow(inc[inc$ID == i,]))])
  }
  return(inc)
  
}


# Code for running a single simulation and returning a dataframe of incidence in cattle per unit time (1/milking frequency)

simulation <- function(df_fields, df_milk,  params, time=0, n_days=100,
                                milking_time=6/24, milking_frequency=2,
                                cleaning_time=NA,
                                N_fields=5){

  inc <- data.frame()
  
  while_con <- sum(df_fields$N_E+df_fields$N_I)>0 | any((time-df_milk$time_inf)<5, na.rm=TRUE)
  i=1
  
  while(while_con){
    
    df_fields <- g_simulation_period(df_fields=df_fields,  params=params, time=time, next_milking_time = milking_time)
    time=milking_time
    
    
    ### Perform milking stuff here
    
    post_milking_outcomes <- milking_period(df_fields=df_fields, df_milk=df_milk,  params=params, time=time, cleaning_time=cleaning_time)
    df_fields <- post_milking_outcomes[[1]]
    df_milk <- post_milking_outcomes[[2]]
    time <- post_milking_outcomes[[3]]
    row <- data.frame(time=i/milking_frequency,
                      N_S = df_fields$N_S,
                      N_E = df_fields$N_E,
                      N_I = df_fields$N_I,
                      N_R = df_fields$N_R,
                      ID = df_fields$ID)
    inc <- rbind(inc, row)
    ###############################
    milking_time <- milking_time + 1/milking_frequency
    
    i = i+1
    while_con <- sum(df_fields$N_E+df_fields$N_I)>0 | any((time-df_milk$time_inf)<5, na.rm=TRUE)
    
  }
  
  inc$incidence <- 0
  for(i in unique(inc$ID)){
    if(length(inc[inc$ID == i,]$incidence )>1){
      inc[inc$ID == i,]$incidence <- c(0,inc[inc$ID == i,]$N_S[1:(nrow(inc[inc$ID == i,])-1)] - inc[inc$ID == i,]$N_S[2:(nrow(inc[inc$ID == i,]))])
    }
  }
  
  inc_end <- inc[inc$time==max(inc$time),]
  inc_end$time <- 999999999
  inc <- rbind(inc, inc_end)
  print(paste("time:", time))
  
  sims_tot <- data.frame()
  for(i in unique(inc$time)){
    df <- inc[inc$time==i,]
    row_df <- data.frame(time=i,
                         N_S = sum(df$N_S),
                         N_E = sum(df$N_E),
                         N_I = sum(df$N_I),
                         N_R = sum(df$N_R),
                         ID = "Total",
                         incidence = sum(df$incidence))
    
    sims_tot <- rbind(sims_tot, row_df)
  }
  
  inc <- rbind(inc, sims_tot)
  return(inc)
  
}

