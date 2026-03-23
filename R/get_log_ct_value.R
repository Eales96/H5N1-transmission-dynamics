
get_log_ct_value <- function(pars,
                             time,
                             t_inf,
                             log_params){
  
  past_peak <- which(time>=t_inf+pars$t_peak)
  before_peak <- which(time>t_inf & time<t_inf+pars$t_peak)
  before_inf <- which(time<=t_inf)
  
  Ct1 <- pars$peak[past_peak] + (time-t_inf[past_peak]-pars$t_peak[past_peak])*pars$decay[past_peak]
  Ct2 <- pars$peak[before_peak] - (time-t_inf[before_peak]-pars$t_peak[before_peak])*pars$rise[before_peak]
  log_Ct1 <- ( Ct1-log_params$intercept)/log_params$grad
  log_Ct2 <- ( Ct2-log_params$intercept)/log_params$grad
  
  
  log_Ct <- sum(c(exp(log_Ct1),exp(log_Ct2)))
  return(log_Ct)
}


get_ct_value_over_time <- function(inf_log,
                                   ct_params,
                                   log_params,
                                   times = seq(0.5,50,1)){
  
  n_inf <- nrow(inf_log)
  first_row <- inf_log[1,]
  N_cows <- first_row$N_S + first_row$N_E + first_row$N_I + first_row$N_R
  ct_draws <- data.frame(peak = rnorm(n_inf, ct_params$peak_mn, ct_params$peak_sd),
                         t_peak = rnorm(n_inf, ct_params$t_peak_mn, ct_params$t_peak_sd),
                         decay = exp(rnorm(n_inf, log(ct_params$decay_mn), ct_params$decay_sd)),
                         rise = exp(rnorm(n_inf, log(ct_params$rise_mn), ct_params$rise_sd)))
  inf_log <- cbind(inf_log, ct_draws)
  
  
  final <- data.frame()
  
  for(i in times){
    
    log_titre <- get_log_ct_value(t_inf = inf_log$time,
                                  pars = inf_log,
                                  log_params = log_params, 
                                  time=i)
    N_infected <- N_cows - inf_log[tail(which(inf_log$time<i),1),]$N_S
    
    row <- data.frame(time = i,
                      N_inf = N_infected,
                      N_cow = N_cows,
                      log_titre = log_titre,
                      sim=1)
    final <- rbind(final,row)
    
  }
  
  return(final)
}


get_ct_draws <- function(inf_log,
                         ct_params,
                         log_params,
                         times = seq(0.5,50,1),
                         N_draws = 10){
  
  N_sims <- length(unique(inf_log$sim))
  
  output_final <- lapply(1:N_sims, function(j){
    print(j)
    output <- lapply(1:N_draws, function(i){
      draw <- get_ct_value_over_time(inf_log = inf_log[inf_log$sim==j,],
                                     ct_params,
                                     log_params,
                                     times = times)
      draw$sim_Ct <- i
      return(draw)
    })
    
    
    
    output <- rbindlist(output)
    output$sim <- j
    return(output)
  })
  
  output_final <- rbindlist(output_final)
  return(output_final)
}


get_ct_value <- function(final, m=1, log_params, ct_threshold=38){
  
  final$avg_titre <- m*final$log_titre/(m*final$N_inf+final$N_cow)
  final$Ct <- log(final$avg_titre) *log_params$grad + log_params$intercept
  final$prob <- pnorm(ct_threshold, mean=final$Ct, sd=ct_params$sigma, lower.tail=TRUE)
  
  return(final)
}



get_detection_point <- function(ct_sims_prob_ex,
                                first_test=1,
                                freq_test){
  
  
  test_times <- seq(first_test, max(ct_sims_prob_ex$time), by=freq_test)
  ct_tests <- ct_sims_prob_ex[test_times,]
  N_trials <- nrow(ct_tests)
  test_outcomes <- rbinom(N_trials, 1, ct_tests$prob)
  if(1 %in% test_outcomes){
    pos_test <- ct_tests[which(test_outcomes==1)[1],]
    pos_test$detected <- "Yes"
  } else{
    pos_test <- ct_sims_prob_ex[nrow(ct_sims_prob_ex),]
    pos_test$detected <- "No"
  }
  return(pos_test)
}


get_time_until_detection <- function(ct_sims_prob,
                                     freq_test=7,
                                     sims_per_comb=2){
  
  N_sims <- max(ct_sims_prob$sim)
  N_ct <- max(ct_sims_prob$sim_Ct)
  
  detection_points <- data.frame()
  
  for(i in 1:N_sims){
    for(j in 1:N_ct){
      for(k in 1:freq_test){
        
        detect_points <- lapply(1:sims_per_comb, function(l){
          detect_point <- get_detection_point(ct_sims_prob[ct_sims_prob$sim==i&
                                                             ct_sims_prob$sim_Ct==j,],
                                              first_test = k, freq_test = freq_test)
          detect_point$sim_prob <- l
          return(detect_point)
        })
        
        detect_points <- rbindlist(detect_points)
        detection_points <- rbind(detection_points, detect_points)
        
      }
    }
  }
  
  return(detection_points)
}


extract_all_detect_times <- function(ct_sims_prob,
                                     freq_test,
                                     sims_per_comb = 1){
  
  detection_times <- data.frame()
  for(i in freq_test){
    print(i)
    detection_times_temp <- get_time_until_detection(ct_sims_prob = ct_sims_prob,
                                                     freq_test = i,
                                                     sims_per_comb = sims_per_comb)
    detection_times_temp$freq_test <- i
    detection_times <- rbind(detection_times,detection_times_temp)
  }
  return(detection_times)
}