
# Perform event decided by Gillespie step

perform_g_event <- function(df_fields, event_type){
  
  N_fields <- nrow(df_fields)
  
  if(event_type==1){
    probs <- cumsum(df_fields$N_E/ sum(df_fields$N_E) )
    field <- which(runif(1)<probs)[1]
    
    df_fields$N_E[field] = df_fields$N_E[field] - 1
    df_fields$N_I[field] = df_fields$N_I[field] + 1
    
  } else if(event_type==2){
    
    probs <- cumsum(df_fields$N_I/ sum(df_fields$N_I) )
    field <- which(runif(1)<probs)[1]
    
    df_fields$N_I[field] = df_fields$N_I[field] - 1
    df_fields$N_R[field] = df_fields$N_R[field] + 1
    
  } else if(event_type==3){
    
    probs <- cumsum(df_fields$N_S/ sum(df_fields$N_S) )
    field <- which(runif(1)<probs)[1]
    
    df_fields$N_S[field] = df_fields$N_S[field] - 1
    df_fields$N_E[field] = df_fields$N_E[field] + 1
    
  } else{
    field <- event_type-3
    
    df_fields$N_S[field] = df_fields$N_S[field] - 1
    df_fields$N_E[field] = df_fields$N_E[field] + 1
    
  }
  
  df_fields
  
}
