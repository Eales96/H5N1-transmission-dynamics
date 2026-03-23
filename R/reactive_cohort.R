
reactive_cohort <- function(df_fields, N_fields=5){
  
  new_fields <- data.frame()
  N_cows <- df_fields$N_S + df_fields$N_E + df_fields$N_I + df_fields$N_R
  field_size <- N_cows / N_fields
  new_area <- df_fields$A/N_fields
  
  
  for(i in 1:N_fields){
    all_cow_states <- c(rep(1,df_fields$N_S), 
                        rep(2,df_fields$N_E),
                        rep(3,df_fields$N_I),
                        rep(4,df_fields$N_R))
    
    new_field_states <- sample(all_cow_states,N_cows/N_fields, replace = FALSE)
    
    delta_N_S <- sum(new_field_states==1)
    delta_N_E <- sum(new_field_states==2)
    delta_N_I <- sum(new_field_states==3)
    delta_N_R <- sum(new_field_states==4)
    
    delta_N_S <- sum(new_field_states==1)
    delta_N_E <- sum(new_field_states==2)
    delta_N_I <- sum(new_field_states==3)
    delta_N_R <- sum(new_field_states==4)
    
    row_df <- data.frame(N_S = delta_N_S,
                         N_E = delta_N_E,
                         N_I = delta_N_I,
                         N_R = delta_N_R,
                         ID = i+1,
                         M  = 1,
                         A  = new_area)
    
    new_fields <- rbind(row_df, new_fields)
    
    df_fields$N_S <- df_fields$N_S - delta_N_S
    df_fields$N_E <- df_fields$N_E - delta_N_E
    df_fields$N_I <- df_fields$N_I - delta_N_I
    df_fields$N_R <- df_fields$N_R - delta_N_R
    
  }
  
  
  return(new_fields)
  
  
  
}
