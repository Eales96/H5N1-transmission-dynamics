
# Run an entire milking period

milking_period_detailed_infections_only <- function(df_fields, df_milk,  params, time, cleaning_time=NA){
  N_milk <- nrow(df_milk)
  N_fields <- nrow(df_fields)
  
  event_log <- data.frame()
  
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
      g_events <- g_simulation_period_detailed_infections_only(df_fields=df_fields_temp,  params=params, time=time, next_milking_time = time+params$milk_time_cow)
      df_fields_temp <- g_events[[1]]
      
      if(nrow(g_events[[2]])>0){
        event_log_temp <- g_events[[2]]
        if(nrow(event_log_temp)>0){
          event_log_temp <- event_log_temp %>% group_by(ID,time, M, A) %>% summarise(N_S = sum(N_S),
                                                                                     N_I = sum(N_I),
                                                                                     N_R = sum(N_R),
                                                                                     N_E = sum(N_E))
          event_log_temp$event <- "cow2"
          event_log_temp[event_log_temp$ID==i,]$N_S <- event_log_temp[event_log_temp$ID==i,]$N_S + delta_N_S
          event_log_temp[event_log_temp$ID==i,]$N_E <- event_log_temp[event_log_temp$ID==i,]$N_E + delta_N_E
          event_log_temp[event_log_temp$ID==i,]$N_I <- event_log_temp[event_log_temp$ID==i,]$N_I + delta_N_I
          event_log_temp[event_log_temp$ID==i,]$N_R <- event_log_temp[event_log_temp$ID==i,]$N_R + delta_N_R
        }
      } else{
        event_log_temp <- data.frame()
      }
      
      
      
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
      
      if(milking_cows_infected>0){
        event_log_m <- df_fields_temp
        event_log_m$time <- time
        
        
        if(nrow(event_log_m)>0){
          event_log_m <- event_log_m %>% group_by(ID,time, M, A) %>% summarise(N_S = sum(N_S),
                                                                               N_I = sum(N_I),
                                                                               N_R = sum(N_R),
                                                                               N_E = sum(N_E))
          event_log_m$event <- "milk"
        }
        
        
      } else{
        event_log_m <- data.frame()
      }
      

      event_log <- rbind(event_log, event_log_temp, event_log_m)
      #time <- time + params$milk_time_cow
    }
    
    
    df_fields <- df_fields_temp[-i,]
    
    df_fields <- df_fields[order(df_fields$ID),]
    
    g_events <- g_simulation_period_detailed_infections_only(df_fields=df_fields,  params=params, time=time, next_milking_time = time+params$milk_time_field)
    df_fields <- g_events[[1]]
    
    if(nrow(g_events[[2]])>0){
      event_log_temp <- g_events[[2]]
      if(nrow(event_log_temp)>0){
        event_log_temp <- event_log_temp %>% group_by(ID,time, M, A) %>% summarise(N_S = sum(N_S),
                                                                                   N_I = sum(N_I),
                                                                                   N_R = sum(N_R),
                                                                                   N_E = sum(N_E))
        event_log_temp$event <- "cow3"
        event_log_temp[event_log_temp$ID==i,]$N_S <- event_log_temp[event_log_temp$ID==i,]$N_S + delta_N_S
        event_log_temp[event_log_temp$ID==i,]$N_E <- event_log_temp[event_log_temp$ID==i,]$N_E + delta_N_E
        event_log_temp[event_log_temp$ID==i,]$N_I <- event_log_temp[event_log_temp$ID==i,]$N_I + delta_N_I
        event_log_temp[event_log_temp$ID==i,]$N_R <- event_log_temp[event_log_temp$ID==i,]$N_R + delta_N_R
      }
    } else{
      event_log_temp <- data.frame()
    }
    
    
    event_log <- rbind(event_log, event_log_temp)
    time = time + params$milk_time_field
  }
  
  post_milking_outcomes <- list(df_fields, df_milk, time, event_log)
  return(post_milking_outcomes)
  
}




