
# Get rates of transitions

get_rates <- function(df_fields, params){
  
  R_sigma <- sum(df_fields$N_E)*params$sigma
  R_gamma <- sum(df_fields$N_I)*params$gamma
  R_alpha0 <- sum(df_fields$N_S)*sum(df_fields$N_I)*params$alpha0/sum(df_fields$A) ############################################################################
  R_alpha1_f <- (df_fields$N_S * df_fields$N_I)*params$alpha1/df_fields$A ####################################################################
  
  c(R_sigma, R_gamma, R_alpha0, R_alpha1_f)
  
}