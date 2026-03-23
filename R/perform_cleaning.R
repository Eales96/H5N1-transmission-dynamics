# Update milking unit status to document additional cleaning

perform_cleaning <- function(df_milk, params, time){
  N_milk <- nrow(df_milk)
  
  df_milk[df_milk$status==1,]$cleans_since_inf <- df_milk[df_milk$status==1,]$cleans_since_inf + 1
  
  return(df_milk)
}