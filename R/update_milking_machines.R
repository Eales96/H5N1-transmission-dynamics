
# Update cow and milking unit status for one timestep

update_milking_machines <- function(df_milk, cow_status, params, time){
  
  df_milk_m2c <- df_milk[df_milk$status==1 & cow_status==1,]
  pos_m2c <- nrow(df_milk_m2c)
  inf_m2c <- runif(pos_m2c) < params$prob_m2c * exp(- (time - df_milk_m2c$time_inf) *params$decay_time) * exp(- (df_milk_m2c$uses_since_inf-1)*params$decay_use ) * exp(- (df_milk_m2c$cleans_since_inf)*params$decay_clean)
  cows_infected <- sum(inf_m2c)
  
  # Performing cow to milking machine transmission events
  pos_c2m <- nrow(df_milk[cow_status==3,])
  inf_c2m <- runif(pos_c2m)<params$prob_c2m
  
  if(sum(inf_c2m)>0){
    df_milk[cow_status==3,][inf_c2m,]$status <- 1
    df_milk[cow_status==3,][inf_c2m,]$time_inf <- time
    df_milk[cow_status==3,][inf_c2m,]$uses_since_inf <- 0
    df_milk[cow_status==3,][inf_c2m,]$cleans_since_inf <- 0
    
  }
  
  
  df_milk[df_milk$status==1,]$uses_since_inf <- df_milk[df_milk$status==1,]$uses_since_inf +1
  
  event_occur <- cows_infected+sum(inf_c2m) > 0
  
  return( list(df_milk, cows_infected, event_occur) )
  
}

