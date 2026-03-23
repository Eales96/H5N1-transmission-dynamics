
# Run an entire gillespie period

g_simulation_period <- function(df_fields,  params, time, next_milking_time){
  
  while(time < next_milking_time){
    event_details <-  get_g_step(df_fields = df_fields, params = params)
    
    if(time+event_details[1]< next_milking_time){
      time <- time+event_details[1]
      df_fields <- perform_g_event(df_fields = df_fields, event_type = event_details[2])
    } else{
      time <- next_milking_time
    }
    
  }
  
  return(df_fields)
  
}