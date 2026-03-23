
# Get timing and event type of possible Gillespie step

get_g_step <- function(df_fields, params){
  
  rates <- get_rates(df_fields = df_fields,
                     params = params)
  
  R <- sum(rates)
  
  rates_norm <- cumsum(rates/R)
  
  t_step <- log(1/runif(1))*(1/R)
  event_type <- which(runif(1)<rates_norm)[1]
  
  c(t_step, event_type)
}