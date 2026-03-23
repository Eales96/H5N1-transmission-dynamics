
# Run an entire gillespie period returning state of df_fields after period
# and infection event log

g_simulation_period_detailed_infections_only <- function(df_fields,  params, time, next_milking_time){
  event_log <- data.frame()
  while(time < next_milking_time){
    event_details <-  get_g_step(df_fields = df_fields, params = params)
    
    if(time+event_details[1]< next_milking_time){
      time <- time+event_details[1]
      df_fields <- perform_g_event(df_fields = df_fields, event_type = event_details[2])
      
      if(event_details[2]%in%c(3,4)){
        new_event <- df_fields
        new_event$time <- time
        new_event$event <- "cow1"
        event_log <- rbind(event_log, new_event)
      }
      
    } else{
      time <- next_milking_time
    }
    
  }
  
  return(list(df_fields, event_log))
  
}
