# Get plottable data on infectiousness of milking units over time
# from milk_log output of simulation_detailed

milk_infectious_profile <- function(milk_log, times=seq(0,100,0.1), params){
  unq_times <- unique(milk_log$time)
  
  inf_profile_df <- data.frame()
  
  for(i in times){
    last_time <- unq_times[tail(which(unq_times < i),1)]
    milk_log_temp <- milk_log[milk_log$time == last_time,]
    milk_log_temp
    
    inf_prof <- sum(params$prob_m2c * exp(- (i - milk_log_temp$time_inf) *params$decay_time) * exp(- (milk_log_temp$uses_since_inf-1)*params$decay_use ) * exp(- (milk_log_temp$cleans_since_inf)*params$decay_clean), na.rm=T)
    
    row <- data.frame(time=i,inf_prof=inf_prof)
    inf_profile_df <- rbind(inf_profile_df, row)
  }
  
  for(i in unq_times){
    last_time <- unq_times[tail(which(unq_times < i),1)]
    milk_log_temp <- milk_log[milk_log$time == last_time,]
    milk_log_temp
    
    inf_prof <- sum(params$prob_m2c * exp(- (i - milk_log_temp$time_inf) *params$decay_time) * exp(- (milk_log_temp$uses_since_inf-1)*params$decay_use ) * exp(- (milk_log_temp$cleans_since_inf)*params$decay_clean), na.rm=T)
    
    row <- data.frame(time=i,inf_prof=inf_prof)
    inf_profile_df <- rbind(inf_profile_df, row)
  }
  
  for(i in (unq_times-0.00000001) ){
    last_time <- unq_times[tail(which(unq_times < i),1)]
    milk_log_temp <- milk_log[milk_log$time == last_time,]
    milk_log_temp
    
    inf_prof <- sum(params$prob_m2c * exp(- (i - milk_log_temp$time_inf) *params$decay_time) * exp(- (milk_log_temp$uses_since_inf-1)*params$decay_use ) * exp(- (milk_log_temp$cleans_since_inf)*params$decay_clean), na.rm=T)
    
    row <- data.frame(time=i,inf_prof=inf_prof)
    inf_profile_df <- rbind(inf_profile_df, row)
  }
  
  return(inf_profile_df)
  
}
