
# Run an entire milking period

milking_period <- function(df_fields, df_milk,  params, time, cleaning_time=NA){
  N_milk <- nrow(df_milk)
  N_fields <- nrow(df_fields)
  
  for(i in 1:N_fields){
    
    if(i%in%cleaning_time){ ######################################################
      df_milk <- perform_cleaning(df_milk)######################################
    } ##########################################################################
    
    df_fields_temp <- rbind(df_fields, df_fields[i,])
    df_fields_temp[N_fields+1,1:4] <- c(0,0,0,0)
    while(sum(df_fields_temp[i,1:4])>0){
      
      
      all_cow_states <- c(rep(1,df_fields_temp$N_S[i]), 
                          rep(2,df_fields_temp$N_E[i]),
                          rep(3,df_fields_temp$N_I[i]),
                          rep(4,df_fields_temp$N_R[i]))
      
      if(length(all_cow_states)>=N_milk){
        milking_cow_states <- sample(all_cow_states, N_milk, replace=FALSE)
      } else{
        milking_cow_states <- c(all_cow_states, rep(100,N_milk-length(all_cow_states)) )
      }
      
      
      delta_N_S <- sum(milking_cow_states==1)
      delta_N_E <- sum(milking_cow_states==2)
      delta_N_I <- sum(milking_cow_states==3)
      delta_N_R <- sum(milking_cow_states==4)
      
      df_fields_temp$N_S[i] <- df_fields_temp$N_S[i] - delta_N_S
      df_fields_temp$N_E[i] <- df_fields_temp$N_E[i] - delta_N_E
      df_fields_temp$N_I[i] <- df_fields_temp$N_I[i] - delta_N_I
      df_fields_temp$N_R[i] <- df_fields_temp$N_R[i] - delta_N_R
      
      ###########################################################################################
      new_states <- update_milking_machines(df_milk, cow_status = milking_cow_states, params, time)
      df_fields_temp <- g_simulation_period(df_fields=df_fields_temp,  params=params, time=time, next_milking_time = time+params$milk_time_cow)
      time = time + params$milk_time_cow
      ###########################################################################################
      
      
      df_milk <- new_states[[1]]
      milking_cows_infected <- new_states[[2]]
      event_occur <- new_states[[3]]
      
      delta_N_S <- delta_N_S - milking_cows_infected
      delta_N_E <- delta_N_E + milking_cows_infected
      delta_N_I <- delta_N_I 
      delta_N_R <- delta_N_R
      
      df_fields_temp$N_S[N_fields+1]  <- df_fields_temp$N_S[N_fields+1] + delta_N_S
      df_fields_temp$N_E[N_fields+1]  <- df_fields_temp$N_E[N_fields+1] + delta_N_E
      df_fields_temp$N_I[N_fields+1]  <- df_fields_temp$N_I[N_fields+1] + delta_N_I
      df_fields_temp$N_R[N_fields+1]  <- df_fields_temp$N_R[N_fields+1] + delta_N_R
      
      
      #time <- time + params$milk_time_cow
    }
    
    
    df_fields <- df_fields_temp[-i,]
    
    df_fields <- df_fields[order(df_fields$ID),]
    
    df_fields <- g_simulation_period(df_fields=df_fields,  params=params, time=time, next_milking_time = time+params$milk_time_field)
    time = time + params$milk_time_field
  }
  
  post_milking_outcomes <- list(df_fields, df_milk, time)
  return(post_milking_outcomes)
  
}
