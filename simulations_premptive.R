setwd("C:/Users/EALESO/OneDrive - The University of Melbourne/Projects/H5N1 cow transmission")


# 1. Simulations to investigate effect of pre-emptive cohorting

################################################################################

library(patchwork)
library(data.table)
library(arrow)
################################################################################

# Functions for initialising cows and milking machines (with no infection)
source("R/initialise_fields.R")
source("R/initialise_milking_machines.R")

# Functions for simulating cow-to-cow transmission
source("R/get_rates.R")
source("R/get_g_step.R")
source("R/perform_g_event.R")
source("R/g_simulation_period.R")

# Functions for simulating cow-to-milking unit transmission
source("R/update_milking_machines.R")
source("R/perform_cleaning.R")
source("R/milking_period.R")

# Functions for performing simulations
source("R/simulation.R")
source("R/run_multi_sims.R")

##########################################################################################################
# Setting up main parameter sets for simulations
##########################################################################################################

params <- data.frame(sigma=1/1.1, # rate E->I, Eales et al 2026 (time until minimum Ct reached)
                     gamma=1/7.8, # rate I->R (recovery), Eales et al 2026 (average duration of infectiousness)
                     alpha0 = 0.0, # between cohort cow-cow transmission, set to zero for main analysis
                     alpha1 = 0.0, # within cohort cow-cow transmission, varied for main analysis
                     prob_c2m = 0.0, # probability of infected cow contaminating milking unit, varied for main analysis
                     prob_m2c = 0.0, # probability of contaminated milking unit infecting cow, varied for main analysis
                     decay_use = 0.0, # Set viral decay with each use to 0 (could set virus concentration to also decay exponentially with number of uses)
                     decay_time = 1.15*24, # Assuming viral decay rate is approximately 1.15 hr^-1 (Le Sage et al 2024 Figure 1B (approximate half life of virus on rubber liners)
                     decay_clean = 0.0, # This will need to be changed when investigating the timing of cleaning
                     milk_time_cow = 5/(60*24), # Assume average time to milk one cow is 5 minutes
                     milk_time_field = 0/(60*24)) # Assume no time between milking cohorts


params_scen1 <- params
params_scen2 <- params
params_scen3 <- params
params_scen4 <- params
params_scen5 <- params

################################################################################
# Only cow-to-cow transmission
params_scen1$alpha1 <- 0.34
params_scen1$prob_c2m  <- 0.0
params_scen1$prob_m2c  <- 0.0

################################################################################
# Predominantly cow-to-milking unit transmission with some cow-to-cow transnmission
params_scen3$alpha1 <- 0.07
params_scen3$prob_c2m  <- 0.16
params_scen3$prob_m2c  <- 0.16
################################################################################
# Only cow-to-milking unit transmission
params_scen5$alpha1 <- 0.0
params_scen5$prob_c2m  <- 0.18
params_scen5$prob_m2c  <- 0.18

##########################################################################################################
# Setting up initial conditions for fields and milking units
##########################################################################################################
# Mote that when initialising fields we set the area to be the same as the number of animals in each field

################################################################################
# Single cohort of 1000, with 10 initial infections
df_fields1 <- initialise_fields(1000, A=1000) 
df_fields1$N_S[1] <- df_fields1$N_S[1] - 10
df_fields1$N_E[1] <- 10

################################################################################
# Two cohorts of 500, with 5 initial infections per cohort
df_fields2 <- initialise_fields(rep(500,2), A=500) 
df_fields2$N_S[1] <- df_fields2$N_S[1] - 5
df_fields2$N_S[2] <- df_fields2$N_S[2] - 5
df_fields2$N_E[1] <- 5
df_fields2$N_E[2] <- 5

################################################################################
# Three cohorts of 333/334, with 3/4 initial infections per cohort
df_fields3 <- initialise_fields(c(333,333,334), A=c(333,333,334)) 
df_fields3$N_S[1] <- df_fields3$N_S[1] - 3
df_fields3$N_S[2] <- df_fields3$N_S[2] - 4
df_fields3$N_S[3] <- df_fields3$N_S[3] - 3
df_fields3$N_E[1] <- 3
df_fields3$N_E[2] <- 4
df_fields3$N_E[3] <- 3

################################################################################
# Four cohorts of 250, with 3/4 initial infections per cohort
df_fields4 <- initialise_fields(rep(250,4), A=250) 
df_fields4$N_S[1] <- df_fields4$N_S[1] - 3
df_fields4$N_S[2] <- df_fields4$N_S[2] - 2
df_fields4$N_S[3] <- df_fields4$N_S[3] - 3
df_fields4$N_S[4] <- df_fields4$N_S[4] - 2
df_fields4$N_E[1] <- 3
df_fields4$N_E[2] <- 2
df_fields4$N_E[3] <- 3
df_fields4$N_E[4] <- 2

################################################################################
# Five cohorts of 200, with 2 initial infections per cohort
df_fields5 <- initialise_fields(rep(200,5), A=200) 
df_fields5$N_S[1] <- df_fields5$N_S[1] - 2
df_fields5$N_S[2] <- df_fields5$N_S[2] - 2
df_fields5$N_S[3] <- df_fields5$N_S[3] - 2
df_fields5$N_S[4] <- df_fields5$N_S[4] - 2
df_fields5$N_S[5] <- df_fields5$N_S[5] - 2
df_fields5$N_E[1] <- 2
df_fields5$N_E[2] <- 2
df_fields5$N_E[3] <- 2
df_fields5$N_E[4] <- 2
df_fields5$N_E[5] <- 2

################################################################################
# Ten cohorts of 100, with 1 initial infections per cohort
df_fields10 <- initialise_fields(rep(100,10), A=100) 
df_fields10$N_S[1] <- df_fields10$N_S[1] - 1
df_fields10$N_S[2] <- df_fields10$N_S[2] - 1
df_fields10$N_S[3] <- df_fields10$N_S[3] - 1
df_fields10$N_S[4] <- df_fields10$N_S[4] - 1
df_fields10$N_S[5] <- df_fields10$N_S[5] - 1
df_fields10$N_S[6] <- df_fields10$N_S[6] - 1
df_fields10$N_S[7] <- df_fields10$N_S[7] - 1
df_fields10$N_S[8] <- df_fields10$N_S[8] - 1
df_fields10$N_S[9] <- df_fields10$N_S[9] - 1
df_fields10$N_S[10] <- df_fields10$N_S[10] - 1
df_fields10$N_E[1] <- 1
df_fields10$N_E[2] <- 1
df_fields10$N_E[3] <- 1
df_fields10$N_E[4] <- 1
df_fields10$N_E[5] <- 1
df_fields10$N_E[6] <- 1
df_fields10$N_E[7] <- 1
df_fields10$N_E[8] <- 1
df_fields10$N_E[9] <- 1
df_fields10$N_E[10] <- 1

################################################################################
# Five cohorts of 200, with initial infections in first cohort
df_fields5_1 <- initialise_fields(rep(200,5), A=200) 
df_fields5_1$N_S[1] <- df_fields5_1$N_S[1] - 10
df_fields5_1$N_E[1] <- 10

################################################################################
# Five cohorts of 200, with initial infections in second cohort
df_fields5_2 <- initialise_fields(rep(200,5), A=200) 
df_fields5_2$N_S[2] <- df_fields5_2$N_S[2] - 10
df_fields5_2$N_E[2] <- 10

################################################################################
# Five cohorts of 200, with initial infections in third cohort
df_fields5_3 <- initialise_fields(rep(200,5), A=200) 
df_fields5_3$N_S[3] <- df_fields5_3$N_S[3] - 10
df_fields5_3$N_E[3] <- 10

################################################################################
# Five cohorts of 200, with initial infections in fourth cohort
df_fields5_4 <- initialise_fields(rep(200,5), A=200) 
df_fields5_4$N_S[4] <- df_fields5_4$N_S[4] - 10
df_fields5_4$N_E[4] <- 10

################################################################################
# Five cohorts of 200, with initial infections in fifth cohort
df_fields5_5 <- initialise_fields(rep(200,5), A=200) 
df_fields5_5$N_S[5] <- df_fields5_5$N_S[5] - 10
df_fields5_5$N_E[5] <- 10

################################################################################
# Intialise milking units for main scenarios
df_milk <- initialise_milking_machines(50)

##########################################################################################################
# Running simulations for scenarios comparing number of cohorts
##########################################################################################################

################################################################################
# Simulations for single cohort
set.seed(12345)
sims_p1_c1 <- run_multi_sims(df_fields = df_fields1,
                             df_milk = df_milk,
                             params = params_scen1,
                             n_days=200,
                             N_sims = 100)
set.seed(12345)
sims_p3_c1 <- run_multi_sims(df_fields = df_fields1,
                             df_milk = df_milk,
                             params = params_scen3,
                             n_days=200,
                             N_sims = 100)
set.seed(12345)
sims_p5_c1 <- run_multi_sims(df_fields = df_fields1,
                             df_milk = df_milk,
                             params = params_scen5,
                             n_days=200,
                             N_sims = 100)

################################################################################
# Confirm parameters give a median of approximately 80%
temp1 <- sims_p1_c1[sims_p1_c1$time==200,]
quantile(temp1$N_R, c(0.05,0.5,0.95))/1000
ggplot(sims_p1_c1, aes(x=time, y=N_I, group=sim))+
  geom_line()

temp3 <- sims_p3_c1[sims_p3_c1$time==200,]
quantile(temp3$N_R, c(0.05,0.5,0.95))/1000
ggplot(sims_p5_c1, aes(x=time, y=N_I, group=sim))+
  geom_line()

temp5 <- sims_p5_c1[sims_p5_c1$time==200,]
quantile(temp5$N_R, c(0.05,0.5,0.95))/1000
ggplot(sims_p5_c1, aes(x=time, y=N_R, group=sim))+
  geom_line()

################################################################################
# Simulations for two cohorts
set.seed(12345)
sims_p1_c2 <- run_multi_sims(df_fields = df_fields2,
                             df_milk = df_milk,
                             params = params_scen1,
                             n_days=200,
                             N_sims = 100)
set.seed(12345)
sims_p3_c2 <- run_multi_sims(df_fields = df_fields2,
                             df_milk = df_milk,
                             params = params_scen3,
                             n_days=200,
                             N_sims = 100)
set.seed(12345)
sims_p5_c2 <- run_multi_sims(df_fields = df_fields2,
                             df_milk = df_milk,
                             params = params_scen5,
                             n_days=200,
                             N_sims = 100)

################################################################################
# Simulations for three cohorts
set.seed(12345)
sims_p1_c3 <- run_multi_sims(df_fields = df_fields3,
                             df_milk = df_milk,
                             params = params_scen1,
                             n_days=200,
                             N_sims = 100)
set.seed(12345)
sims_p3_c3 <- run_multi_sims(df_fields = df_fields3,
                             df_milk = df_milk,
                             params = params_scen3,
                             n_days=200,
                             N_sims = 100)
set.seed(12345)
sims_p5_c3 <- run_multi_sims(df_fields = df_fields3,
                             df_milk = df_milk,
                             params = params_scen5,
                             n_days=200,
                             N_sims = 100)

################################################################################
# Simulations for four cohorts
set.seed(12345)
sims_p1_c4 <- run_multi_sims(df_fields = df_fields4,
                             df_milk = df_milk,
                             params = params_scen1,
                             n_days=200,
                             N_sims = 100)
set.seed(12345)
sims_p3_c4 <- run_multi_sims(df_fields = df_fields4,
                             df_milk = df_milk,
                             params = params_scen3,
                             n_days=200,
                             N_sims = 100)
set.seed(12345)
sims_p5_c4 <- run_multi_sims(df_fields = df_fields4,
                             df_milk = df_milk,
                             params = params_scen5,
                             n_days=200,
                             N_sims = 100)


################################################################################
# Simulations for five cohorts
set.seed(12345)
sims_p1_c5 <- run_multi_sims(df_fields = df_fields5,
                             df_milk = df_milk,
                             params = params_scen1,
                             n_days=200,
                             N_sims = 100)
set.seed(12345)
sims_p3_c5 <- run_multi_sims(df_fields = df_fields5,
                             df_milk = df_milk,
                             params = params_scen3,
                             n_days=200,
                             N_sims = 100)
set.seed(12345)
sims_p5_c5 <- run_multi_sims(df_fields = df_fields5,
                             df_milk = df_milk,
                             params = params_scen5,
                             n_days=200,
                             N_sims = 100)

################################################################################
# Simulations for ten cohorts
set.seed(12345)
sims_p1_c10 <- run_multi_sims(df_fields = df_fields10,
                             df_milk = df_milk,
                             params = params_scen1,
                             n_days=200,
                             N_sims = 100)
set.seed(12345)
sims_p3_c10 <- run_multi_sims(df_fields = df_fields10,
                             df_milk = df_milk,
                             params = params_scen3,
                             n_days=200,
                             N_sims = 100)
set.seed(12345)
sims_p5_c10 <- run_multi_sims(df_fields = df_fields10,
                             df_milk = df_milk,
                             params = params_scen5,
                             n_days=200,
                             N_sims = 100)


################################################################################
# Saving outcomes of simulations
# 1 cohort
sims_p1_c1$n_cohorts <- 1
sims_p1_c1$param_set <- 1

sims_p3_c1$n_cohorts <- 1
sims_p3_c1$param_set <- 3

sims_p5_c1$n_cohorts <- 1
sims_p5_c1$param_set <- 5

# 2 cohorts
sims_p1_c2$n_cohorts <- 2
sims_p1_c2$param_set <- 1

sims_p3_c2$n_cohorts <- 2
sims_p3_c2$param_set <- 3

sims_p5_c2$n_cohorts <- 2
sims_p5_c2$param_set <- 5

# 3 cohorts
sims_p1_c3$n_cohorts <- 3
sims_p1_c3$param_set <- 1

sims_p3_c3$n_cohorts <- 3
sims_p3_c3$param_set <- 3

sims_p5_c3$n_cohorts <- 3
sims_p5_c3$param_set <- 5

# 4 cohorts
sims_p1_c4$n_cohorts <- 4
sims_p1_c4$param_set <- 1

sims_p3_c4$n_cohorts <- 4
sims_p3_c4$param_set <- 3

sims_p5_c4$n_cohorts <- 4
sims_p5_c4$param_set <- 5
# 5 cohorts
sims_p1_c5$n_cohorts <- 5
sims_p1_c5$param_set <- 1

sims_p3_c5$n_cohorts <- 5
sims_p3_c5$param_set <- 3

sims_p5_c5$n_cohorts <- 5
sims_p5_c5$param_set <- 5

# 10 cohorts
sims_p1_c10$n_cohorts <- 10
sims_p1_c10$param_set <- 1

sims_p3_c10$n_cohorts <- 10
sims_p3_c10$param_set <- 3

sims_p5_c10$n_cohorts <- 10
sims_p5_c10$param_set <- 5

sims_nChorts <- rbind(sims_p1_c1, sims_p3_c1, sims_p5_c1,
                      sims_p1_c2, sims_p3_c2, sims_p5_c2,
                      sims_p1_c3, sims_p3_c3, sims_p5_c3,
                      sims_p1_c4, sims_p3_c4, sims_p5_c4,
                      sims_p1_c5, sims_p3_c5, sims_p5_c5,
                      sims_p1_c10, sims_p3_c10, sims_p5_c10)
write_parquet(sims_nChorts, "simulations/sims_nCohortsFINAL.parquet")

##########################################################################################################
# Running simulations for scenarios comparing cohort of seeding
##########################################################################################################

################################################################################
# Simulations for seeding in first cohort
set.seed(12345)
sims_p1_c5_s1 <- run_multi_sims(df_fields = df_fields5_1,
                                df_milk = df_milk,
                                params = params_scen1,
                                n_days=200,
                                N_sims = 100)
set.seed(12345)
sims_p3_c5_s1 <- run_multi_sims(df_fields = df_fields5_1,
                                df_milk = df_milk,
                                params = params_scen3,
                                n_days=200,
                                N_sims = 100)
set.seed(12345)
sims_p5_c5_s1 <- run_multi_sims(df_fields = df_fields5_1,
                                df_milk = df_milk,
                                params = params_scen5,
                                n_days=200,
                                N_sims = 100)

################################################################################
# Simulations for seeding in second cohort
set.seed(12345)
sims_p1_c5_s2 <- run_multi_sims(df_fields = df_fields5_2,
                                df_milk = df_milk,
                                params = params_scen1,
                                n_days=200,
                                N_sims = 100)
set.seed(12345)
sims_p3_c5_s2 <- run_multi_sims(df_fields = df_fields5_2,
                                df_milk = df_milk,
                                params = params_scen3,
                                n_days=200,
                                N_sims = 100)
set.seed(12345)
sims_p5_c5_s2 <- run_multi_sims(df_fields = df_fields5_2,
                                df_milk = df_milk,
                                params = params_scen5,
                                n_days=200,
                                N_sims = 100)

################################################################################
# Simulations for seeding in third cohort
set.seed(12345)
sims_p1_c5_s3 <- run_multi_sims(df_fields = df_fields5_3,
                                df_milk = df_milk,
                                params = params_scen1,
                                n_days=200,
                                N_sims = 100)
set.seed(12345)
sims_p3_c5_s3 <- run_multi_sims(df_fields = df_fields5_3,
                                df_milk = df_milk,
                                params = params_scen3,
                                n_days=200,
                                N_sims = 100)
set.seed(12345)
sims_p5_c5_s3 <- run_multi_sims(df_fields = df_fields5_3,
                                df_milk = df_milk,
                                params = params_scen5,
                                n_days=200,
                                N_sims = 100)

################################################################################
# Simulations for seeding in fourth cohort
set.seed(12345)
sims_p1_c5_s4 <- run_multi_sims(df_fields = df_fields5_4,
                                df_milk = df_milk,
                                params = params_scen1,
                                n_days=200,
                                N_sims = 100)
set.seed(12345)
sims_p3_c5_s4 <- run_multi_sims(df_fields = df_fields5_4,
                                df_milk = df_milk,
                                params = params_scen3,
                                n_days=200,
                                N_sims = 100)
set.seed(12345)
sims_p5_c5_s4 <- run_multi_sims(df_fields = df_fields5_4,
                                df_milk = df_milk,
                                params = params_scen5,
                                n_days=200,
                                N_sims = 100)

################################################################################
# Simulations for seeding in fifth cohort
set.seed(12345)
sims_p1_c5_s5 <- run_multi_sims(df_fields = df_fields5_5,
                                df_milk = df_milk,
                                params = params_scen1,
                                n_days=200,
                                N_sims = 100)
set.seed(12345)
sims_p3_c5_s5 <- run_multi_sims(df_fields = df_fields5_5,
                                df_milk = df_milk,
                                params = params_scen3,
                                n_days=200,
                                N_sims = 100)
set.seed(12345)
sims_p5_c5_s5 <- run_multi_sims(df_fields = df_fields5_5,
                                df_milk = df_milk,
                                params = params_scen5,
                                n_days=200,
                                N_sims = 100)

################################################################################
# Saving outcomes of simulations
# 1st cohort
sims_p1_c5_s1$seed_cohorts <- 1
sims_p1_c5_s1$param_set <- 1

sims_p3_c5_s1$seed_cohorts <- 1
sims_p3_c5_s1$param_set <- 3

sims_p5_c5_s1$seed_cohorts <- 1
sims_p5_c5_s1$param_set <- 5

# 2nd cohort
sims_p1_c5_s2$seed_cohorts <- 2
sims_p1_c5_s2$param_set <- 1

sims_p3_c5_s2$seed_cohorts <- 2
sims_p3_c5_s2$param_set <- 3

sims_p5_c5_s2$seed_cohorts <- 2
sims_p5_c5_s2$param_set <- 5

# 3rd cohort
sims_p1_c5_s3$seed_cohorts <- 3
sims_p1_c5_s3$param_set <- 1

sims_p3_c5_s3$seed_cohorts <- 3
sims_p3_c5_s3$param_set <- 3

sims_p5_c5_s3$seed_cohorts <- 3
sims_p5_c5_s3$param_set <- 5

# 4th cohort
sims_p1_c5_s4$seed_cohorts <- 4
sims_p1_c5_s4$param_set <- 1

sims_p3_c5_s4$seed_cohorts <- 4
sims_p3_c5_s4$param_set <- 3

sims_p5_c5_s4$seed_cohorts <- 4
sims_p5_c5_s4$param_set <- 5

# 5th cohort
sims_p1_c5_s5$seed_cohorts <- 5
sims_p1_c5_s5$param_set <- 1

sims_p3_c5_s5$seed_cohorts <- 5
sims_p3_c5_s5$param_set <- 3

sims_p5_c5_s5$seed_cohorts <- 5
sims_p5_c5_s5$param_set <- 5

sims_seedCohorts <- rbind(sims_p1_c5_s1, sims_p1_c5_s2, sims_p1_c5_s3, sims_p1_c5_s4, sims_p1_c5_s5,
                          sims_p3_c5_s1, sims_p3_c5_s2, sims_p3_c5_s3, sims_p3_c5_s4, sims_p3_c5_s5,
                          sims_p5_c5_s1, sims_p5_c5_s2, sims_p5_c5_s3, sims_p5_c5_s4, sims_p5_c5_s5)
write_parquet(sims_seedCohorts, "simulations/sims_seedCohortsFINAL.parquet")
##########################################################################################################
# Running simulations for comparing cleaning_times 90%
##########################################################################################################
params_scen1$decay_clean <- -log(0.1) # Cleaning reduces infectiousness of all units by 90%
params_scen3$decay_clean <- -log(0.1)
params_scen5$decay_clean <- -log(0.1)
################################################################################
# Simulations for cleaning between fifth and first unit


set.seed(12345)
sims_p3_c5_clb1_90 <- run_multi_sims(df_fields = df_fields5,
                                  df_milk = df_milk,
                                  params = params_scen3,
                                  n_days=200,
                                  N_sims = 100,
                                  cleaning_time = 1)
set.seed(12345)
sims_p5_c5_clb1_90 <- run_multi_sims(df_fields = df_fields5,
                                  df_milk = df_milk,
                                  params = params_scen5,
                                  n_days=200,
                                  N_sims = 100,
                                  cleaning_time = 1)

################################################################################
# Simulations for cleaning between first and second unit


set.seed(12345)
sims_p3_c5_clb2_90 <- run_multi_sims(df_fields = df_fields5,
                                  df_milk = df_milk,
                                  params = params_scen3,
                                  n_days=200,
                                  N_sims = 100,
                                  cleaning_time = 2)
set.seed(12345)
sims_p5_c5_clb2_90 <- run_multi_sims(df_fields = df_fields5,
                                  df_milk = df_milk,
                                  params = params_scen5,
                                  n_days=200,
                                  N_sims = 100,
                                  cleaning_time = 2)

################################################################################
# Simulations for cleaning between second and third unit

set.seed(12345)
sims_p3_c5_clb3_90 <- run_multi_sims(df_fields = df_fields5,
                                  df_milk = df_milk,
                                  params = params_scen3,
                                  n_days=200,
                                  N_sims = 100,
                                  cleaning_time = 3)
set.seed(12345)
sims_p5_c5_clb3_90 <- run_multi_sims(df_fields = df_fields5,
                                  df_milk = df_milk,
                                  params = params_scen5,
                                  n_days=200,
                                  N_sims = 100,
                                  cleaning_time = 3)

################################################################################
# Simulations for cleaning between fifth and first unit

set.seed(12345)
sims_p3_c5_clb4_90 <- run_multi_sims(df_fields = df_fields5,
                                  df_milk = df_milk,
                                  params = params_scen3,
                                  n_days=200,
                                  N_sims = 100,
                                  cleaning_time = 4)
set.seed(12345)
sims_p5_c5_clb4_90 <- run_multi_sims(df_fields = df_fields5,
                                  df_milk = df_milk,
                                  params = params_scen5,
                                  n_days=200,
                                  N_sims = 100,
                                  cleaning_time = 4)

################################################################################
# Simulations for cleaning between fourth and fifth unit

set.seed(12345)
sims_p3_c5_clb5_90 <- run_multi_sims(df_fields = df_fields5,
                                  df_milk = df_milk,
                                  params = params_scen3,
                                  n_days=200,
                                  N_sims = 100,
                                  cleaning_time = 5)
set.seed(12345)
sims_p5_c5_clb5_90 <- run_multi_sims(df_fields = df_fields5,
                                  df_milk = df_milk,
                                  params = params_scen5,
                                  n_days=200,
                                  N_sims = 100,
                                  cleaning_time = 5)


################################################################################
##########################################################################################################
# Running simulations for comparing cleaning_times
##########################################################################################################
params_scen1$decay_clean <- -log(0.5) # Cleaning reduces infectiousness of all units by 90%
params_scen3$decay_clean <- -log(0.5)
params_scen5$decay_clean <- -log(0.5)
################################################################################
# Simulations for cleaning between fifth and first unit


set.seed(12345)
sims_p3_c5_clb1_50 <- run_multi_sims(df_fields = df_fields5,
                                  df_milk = df_milk,
                                  params = params_scen3,
                                  n_days=200,
                                  N_sims = 100,
                                  cleaning_time = 1)
set.seed(12345)
sims_p5_c5_clb1_50 <- run_multi_sims(df_fields = df_fields5,
                                  df_milk = df_milk,
                                  params = params_scen5,
                                  n_days=200,
                                  N_sims = 100,
                                  cleaning_time = 1)

################################################################################
# Simulations for cleaning between first and second unit


set.seed(12345)
sims_p3_c5_clb2_50 <- run_multi_sims(df_fields = df_fields5,
                                     df_milk = df_milk,
                                     params = params_scen3,
                                     n_days=200,
                                     N_sims = 100,
                                     cleaning_time = 1)
set.seed(12345)
sims_p5_c5_clb2_50 <- run_multi_sims(df_fields = df_fields5,
                                     df_milk = df_milk,
                                     params = params_scen5,
                                     n_days=200,
                                     N_sims = 100,
                                     cleaning_time = 2)

################################################################################
# Simulations for cleaning between second and third unit

set.seed(12345)
sims_p3_c5_clb3_50 <- run_multi_sims(df_fields = df_fields5,
                                     df_milk = df_milk,
                                     params = params_scen3,
                                     n_days=200,
                                     N_sims = 100,
                                     cleaning_time = 3)
set.seed(12345)
sims_p5_c5_clb3_50 <- run_multi_sims(df_fields = df_fields5,
                                     df_milk = df_milk,
                                     params = params_scen5,
                                     n_days=200,
                                     N_sims = 100,
                                     cleaning_time = 3)

################################################################################
# Simulations for cleaning between fifth and first unit

set.seed(12345)
sims_p3_c5_clb4_50 <- run_multi_sims(df_fields = df_fields5,
                                     df_milk = df_milk,
                                     params = params_scen3,
                                     n_days=200,
                                     N_sims = 100,
                                     cleaning_time = 4)
set.seed(12345)
sims_p5_c5_clb4_50 <- run_multi_sims(df_fields = df_fields5,
                                     df_milk = df_milk,
                                     params = params_scen5,
                                     n_days=200,
                                     N_sims = 100,
                                     cleaning_time = 4)

################################################################################
# Simulations for cleaning between fourth and fifth unit

set.seed(12345)
sims_p3_c5_clb5_50 <- run_multi_sims(df_fields = df_fields5,
                                     df_milk = df_milk,
                                     params = params_scen3,
                                     n_days=200,
                                     N_sims = 100,
                                     cleaning_time = 5)
set.seed(12345)
sims_p5_c5_clb5_50 <- run_multi_sims(df_fields = df_fields5,
                                     df_milk = df_milk,
                                     params = params_scen5,
                                     n_days=200,
                                     N_sims = 100,
                                     cleaning_time = 5)


sims_p3_c5_clb1_90$cleaning_before <- 1
sims_p3_c5_clb2_90$cleaning_before <- 2
sims_p3_c5_clb3_90$cleaning_before <- 3
sims_p3_c5_clb4_90$cleaning_before <- 4
sims_p3_c5_clb5_90$cleaning_before <- 5

sims_p3_c5_clb1_50$cleaning_before <- 1
sims_p3_c5_clb2_50$cleaning_before <- 2
sims_p3_c5_clb3_50$cleaning_before <- 3
sims_p3_c5_clb4_50$cleaning_before <- 4
sims_p3_c5_clb5_50$cleaning_before <- 5

sims_p3_c5_clb1_90$cleaning_effect <- 90
sims_p3_c5_clb2_90$cleaning_effect <- 90
sims_p3_c5_clb3_90$cleaning_effect <- 90
sims_p3_c5_clb4_90$cleaning_effect <- 90
sims_p3_c5_clb5_90$cleaning_effect <- 90

sims_p3_c5_clb1_50$cleaning_effect <- 50
sims_p3_c5_clb2_50$cleaning_effect <- 50
sims_p3_c5_clb3_50$cleaning_effect <- 50
sims_p3_c5_clb4_50$cleaning_effect <- 50
sims_p3_c5_clb5_50$cleaning_effect <- 50

sims_p5_c5_clb1_90$cleaning_before <- 1
sims_p5_c5_clb2_90$cleaning_before <- 2
sims_p5_c5_clb3_90$cleaning_before <- 3
sims_p5_c5_clb4_90$cleaning_before <- 4
sims_p5_c5_clb5_90$cleaning_before <- 5

sims_p5_c5_clb1_50$cleaning_before <- 1
sims_p5_c5_clb2_50$cleaning_before <- 2
sims_p5_c5_clb3_50$cleaning_before <- 3
sims_p5_c5_clb4_50$cleaning_before <- 4
sims_p5_c5_clb5_50$cleaning_before <- 5

sims_p5_c5_clb1_90$cleaning_effect <- 90
sims_p5_c5_clb2_90$cleaning_effect <- 90
sims_p5_c5_clb3_90$cleaning_effect <- 90
sims_p5_c5_clb4_90$cleaning_effect <- 90
sims_p5_c5_clb5_90$cleaning_effect <- 90

sims_p5_c5_clb1_50$cleaning_effect <- 50
sims_p5_c5_clb2_50$cleaning_effect <- 50
sims_p5_c5_clb3_50$cleaning_effect <- 50
sims_p5_c5_clb4_50$cleaning_effect <- 50
sims_p5_c5_clb5_50$cleaning_effect <- 50



sims_p3_c5_clb1_90$param_set <- 3
sims_p3_c5_clb2_90$param_set <- 3
sims_p3_c5_clb3_90$param_set <- 3
sims_p3_c5_clb4_90$param_set <- 3
sims_p3_c5_clb5_90$param_set <- 3

sims_p3_c5_clb1_50$param_set <- 3
sims_p3_c5_clb2_50$param_set <- 3
sims_p3_c5_clb3_50$param_set <- 3
sims_p3_c5_clb4_50$param_set <- 3
sims_p3_c5_clb5_50$param_set <- 3

sims_p5_c5_clb1_90$param_set <- 5
sims_p5_c5_clb2_90$param_set <- 5
sims_p5_c5_clb3_90$param_set <- 5
sims_p5_c5_clb4_90$param_set <- 5
sims_p5_c5_clb5_90$param_set <- 5

sims_p5_c5_clb1_50$param_set <- 5
sims_p5_c5_clb2_50$param_set <- 5
sims_p5_c5_clb3_50$param_set <- 5
sims_p5_c5_clb4_50$param_set <- 5
sims_p5_c5_clb5_50$param_set <- 5


sims_clean <- rbind(sims_p3_c5_clb1_90, sims_p3_c5_clb2_90, sims_p3_c5_clb3_90, sims_p3_c5_clb4_90, sims_p3_c5_clb5_90,
                    sims_p3_c5_clb1_50, sims_p3_c5_clb2_50, sims_p3_c5_clb3_50, sims_p3_c5_clb4_50, sims_p3_c5_clb5_50,
                    sims_p5_c5_clb1_90, sims_p5_c5_clb2_90, sims_p5_c5_clb3_90, sims_p5_c5_clb4_90, sims_p5_c5_clb5_90,
                    sims_p5_c5_clb1_50, sims_p5_c5_clb2_50, sims_p5_c5_clb3_50, sims_p5_c5_clb4_50, sims_p5_c5_clb5_50)
write_parquet(sims_clean, "simulations/sims_cleansFINAL.parquet")


params_scen1$decay_clean <- 0
params_scen3$decay_clean <- 0
params_scen5$decay_clean <- 0


##########################################################################################################
# Running simulations for scenarios comparing cohort of seeding
# when only cow-cow transmission 
# with some between cohort transmission
##########################################################################################################

################################################################################
# Simulations for seeding in first cohort

params_scen1_bct1 <- params_scen1
params_scen1_bct2 <- params_scen1
params_scen1_bct3 <- params_scen1
params_scen1_bct4 <- params_scen1


params_scen1_bct1$alpha0 <- params_scen1_bct1$alpha1*0.001
params_scen1_bct1$alpha1 <- params_scen1_bct1$alpha1*0.999

params_scen1_bct2$alpha0 <- params_scen1_bct2$alpha1*0.005
params_scen1_bct2$alpha1 <- params_scen1_bct2$alpha1*0.995

params_scen1_bct3$alpha0 <- params_scen1_bct3$alpha1*0.01
params_scen1_bct3$alpha1 <- params_scen1_bct3$alpha1*0.99

params_scen1_bct4$alpha0 <- params_scen1_bct4$alpha1*0.02
params_scen1_bct4$alpha1 <- params_scen1_bct4$alpha1*0.98


params_scen1_bct1$alpha0 + params_scen1_bct1$alpha1
params_scen1_bct2$alpha0 + params_scen1_bct2$alpha1
params_scen1_bct3$alpha0 + params_scen1_bct3$alpha1
params_scen1_bct4$alpha0 + params_scen1_bct4$alpha1

set.seed(12345)
sims_p1_c5_s1_bct1 <- run_multi_sims(df_fields = df_fields5_1,
                                     df_milk = df_milk,
                                     params = params_scen1_bct1,
                                     n_days=200,
                                     N_sims = 100)

set.seed(12345)
sims_p1_c5_s1_bct2 <- run_multi_sims(df_fields = df_fields5_1,
                                     df_milk = df_milk,
                                     params = params_scen1_bct2,
                                     n_days=200,
                                     N_sims = 100)

set.seed(12345)
sims_p1_c5_s1_bct3 <- run_multi_sims(df_fields = df_fields5_1,
                                     df_milk = df_milk,
                                     params = params_scen1_bct3,
                                     n_days=200,
                                     N_sims = 100)

set.seed(12345)
sims_p1_c5_s1_bct4 <- run_multi_sims(df_fields = df_fields5_1,
                                     df_milk = df_milk,
                                     params = params_scen1_bct4,
                                     n_days=200,
                                     N_sims = 100)


sims_p1_c5_s1$seed_cohorts <- 1
sims_p1_c5_s1$param_set <- 1
sims_p1_c5_s1$between_cohort <- 0


sims_p1_c5_s1_bct1$seed_cohorts <- 1
sims_p1_c5_s1_bct1$param_set <- 1
sims_p1_c5_s1_bct1$between_cohort <- 0.001

sims_p1_c5_s1_bct2$seed_cohorts <- 1
sims_p1_c5_s1_bct2$param_set <- 1
sims_p1_c5_s1_bct2$between_cohort <- 0.005

sims_p1_c5_s1_bct3$seed_cohorts <- 1
sims_p1_c5_s1_bct3$param_set <- 1
sims_p1_c5_s1_bct3$between_cohort <- 0.01

sims_p1_c5_s1_bct4$seed_cohorts <- 1
sims_p1_c5_s1_bct4$param_set <- 1
sims_p1_c5_s1_bct4$between_cohort <- 0.02

sims_bct <- rbind(sims_p1_c5_s1,
                    sims_p1_c5_s1_bct1,
                    sims_p1_c5_s1_bct2,
                    sims_p1_c5_s1_bct3,
                    sims_p1_c5_s1_bct4)
write_parquet(sims_bct, "simulations/sims_bctFINAL.parquet")

##########################################################################################################
# Running simulations for scenarios comparing cohort of seeding
# when only cow-milking unit transmissions
# but different decay rates

##########################################################################################################

params_scen5$decay_time

params_scen5_dt1 <- params_scen5
params_scen5_dt2 <- params_scen5
params_scen5_dt3 <- params_scen5

params_scen5_dt1$decay_time <- params_scen5$decay_time/2
params_scen5_dt2$decay_time <- params_scen5$decay_time/4
params_scen5_dt3$decay_time <- params_scen5$decay_time/3

################################################################################
# Simulations for decay time half original

set.seed(12345)
sims_p5_c5_s1_dt1 <- run_multi_sims(df_fields = df_fields5_1,
                                    df_milk = df_milk,
                                    params = params_scen5_dt1,
                                    n_days=200,
                                    N_sims = 100)

set.seed(12345)
sims_p5_c5_s2_dt1 <- run_multi_sims(df_fields = df_fields5_2,
                                    df_milk = df_milk,
                                    params = params_scen5_dt1,
                                    n_days=200,
                                    N_sims = 100)

set.seed(12345)
sims_p5_c5_s3_dt1 <- run_multi_sims(df_fields = df_fields5_3,
                                    df_milk = df_milk,
                                    params = params_scen5_dt1,
                                    n_days=200,
                                    N_sims = 100)

set.seed(12345)
sims_p5_c5_s4_dt1 <- run_multi_sims(df_fields = df_fields5_4,
                                    df_milk = df_milk,
                                    params = params_scen5_dt1,
                                    n_days=200,
                                    N_sims = 100)

set.seed(12345)
sims_p5_c5_s5_dt1 <- run_multi_sims(df_fields = df_fields5_5,
                                    df_milk = df_milk,
                                    params = params_scen5_dt1,
                                    n_days=200,
                                    N_sims = 100)


##################################################################
# Simulations for decay time for quarter original
set.seed(12345)
sims_p5_c5_s1_dt2 <- run_multi_sims(df_fields = df_fields5_1,
                                    df_milk = df_milk,
                                    params = params_scen5_dt2,
                                    n_days=200,
                                    N_sims = 100)

set.seed(12345)
sims_p5_c5_s2_dt2 <- run_multi_sims(df_fields = df_fields5_2,
                                    df_milk = df_milk,
                                    params = params_scen5_dt2,
                                    n_days=200,
                                    N_sims = 100)

set.seed(12345)
sims_p5_c5_s3_dt2 <- run_multi_sims(df_fields = df_fields5_3,
                                    df_milk = df_milk,
                                    params = params_scen5_dt2,
                                    n_days=200,
                                    N_sims = 100)

set.seed(12345)
sims_p5_c5_s4_dt2 <- run_multi_sims(df_fields = df_fields5_4,
                                    df_milk = df_milk,
                                    params = params_scen5_dt2,
                                    n_days=200,
                                    N_sims = 100)

set.seed(12345)
sims_p5_c5_s5_dt2 <- run_multi_sims(df_fields = df_fields5_5,
                                    df_milk = df_milk,
                                    params = params_scen5_dt2,
                                    n_days=200,
                                    N_sims = 100)

##################################################################
# Simulations for decay time for quarter original
set.seed(12345)
sims_p5_c5_s1_dt3 <- run_multi_sims(df_fields = df_fields5_1,
                                    df_milk = df_milk,
                                    params = params_scen5_dt3,
                                    n_days=200,
                                    N_sims = 100)

set.seed(12345)
sims_p5_c5_s2_dt3 <- run_multi_sims(df_fields = df_fields5_2,
                                    df_milk = df_milk,
                                    params = params_scen5_dt3,
                                    n_days=200,
                                    N_sims = 100)

set.seed(12345)
sims_p5_c5_s3_dt3 <- run_multi_sims(df_fields = df_fields5_3,
                                    df_milk = df_milk,
                                    params = params_scen5_dt3,
                                    n_days=200,
                                    N_sims = 100)

set.seed(12345)
sims_p5_c5_s4_dt3 <- run_multi_sims(df_fields = df_fields5_4,
                                    df_milk = df_milk,
                                    params = params_scen5_dt3,
                                    n_days=200,
                                    N_sims = 100)

set.seed(12345)
sims_p5_c5_s5_dt3 <- run_multi_sims(df_fields = df_fields5_5,
                                    df_milk = df_milk,
                                    params = params_scen5_dt3,
                                    n_days=200,
                                    N_sims = 100)

################################################################################

sims_p5_c5_s1$seed_cohorts <- 1
sims_p5_c5_s1$param_set <- 5
sims_p5_c5_s1$decay_time<- "Original"

sims_p5_c5_s2$seed_cohorts <- 2
sims_p5_c5_s2$param_set <- 5
sims_p5_c5_s2$decay_time<- "Original"

sims_p5_c5_s3$seed_cohorts <- 3
sims_p5_c5_s3$param_set <- 5
sims_p5_c5_s3$decay_time<- "Original"

sims_p5_c5_s4$seed_cohorts <- 4
sims_p5_c5_s4$param_set <- 5
sims_p5_c5_s4$decay_time<- "Original"

sims_p5_c5_s5$seed_cohorts <- 5
sims_p5_c5_s5$param_set <- 5
sims_p5_c5_s5$decay_time<- "Original"


sims_p5_c5_s1_dt1$seed_cohorts <- 1
sims_p5_c5_s1_dt1$param_set <- 5
sims_p5_c5_s1_dt1$decay_time<- 1

sims_p5_c5_s2_dt1$seed_cohorts <- 2
sims_p5_c5_s2_dt1$param_set <- 5
sims_p5_c5_s2_dt1$decay_time<- 1

sims_p5_c5_s3_dt1$seed_cohorts <- 3
sims_p5_c5_s3_dt1$param_set <- 5
sims_p5_c5_s3_dt1$decay_time<- 1

sims_p5_c5_s4_dt1$seed_cohorts <- 4
sims_p5_c5_s4_dt1$param_set <- 5
sims_p5_c5_s4_dt1$decay_time<- 1

sims_p5_c5_s5_dt1$seed_cohorts <- 5
sims_p5_c5_s5_dt1$param_set <- 5
sims_p5_c5_s5_dt1$decay_time<- 1


sims_p5_c5_s1_dt2$seed_cohorts <- 1
sims_p5_c5_s1_dt2$param_set <- 5
sims_p5_c5_s1_dt2$decay_time<- 2

sims_p5_c5_s2_dt2$seed_cohorts <- 2
sims_p5_c5_s2_dt2$param_set <- 5
sims_p5_c5_s2_dt2$decay_time<- 2

sims_p5_c5_s3_dt2$seed_cohorts <- 3
sims_p5_c5_s3_dt2$param_set <- 5
sims_p5_c5_s3_dt2$decay_time<- 2

sims_p5_c5_s4_dt2$seed_cohorts <- 4
sims_p5_c5_s4_dt2$param_set <- 5
sims_p5_c5_s4_dt2$decay_time<- 2

sims_p5_c5_s5_dt2$seed_cohorts <- 5
sims_p5_c5_s5_dt2$param_set <- 5
sims_p5_c5_s5_dt2$decay_time<- 2



sims_p5_c5_s1_dt3$seed_cohorts <- 1
sims_p5_c5_s1_dt3$param_set <- 5
sims_p5_c5_s1_dt3$decay_time<- 3

sims_p5_c5_s2_dt3$seed_cohorts <- 2
sims_p5_c5_s2_dt3$param_set <- 5
sims_p5_c5_s2_dt3$decay_time<- 3

sims_p5_c5_s3_dt3$seed_cohorts <- 3
sims_p5_c5_s3_dt3$param_set <- 5
sims_p5_c5_s3_dt3$decay_time<- 3

sims_p5_c5_s4_dt3$seed_cohorts <- 4
sims_p5_c5_s4_dt3$param_set <- 5
sims_p5_c5_s4_dt3$decay_time<- 3

sims_p5_c5_s5_dt3$seed_cohorts <- 5
sims_p5_c5_s5_dt3$param_set <- 5
sims_p5_c5_s5_dt3$decay_time<- 3


sims_dt <- rbind(sims_p5_c5_s1, sims_p5_c5_s2, sims_p5_c5_s3, sims_p5_c5_s4, sims_p5_c5_s5,
                 sims_p5_c5_s1_dt1, sims_p5_c5_s2_dt1, sims_p5_c5_s3_dt1, sims_p5_c5_s4_dt1, sims_p5_c5_s5_dt1,
                 sims_p5_c5_s1_dt2, sims_p5_c5_s2_dt2, sims_p5_c5_s3_dt2, sims_p5_c5_s4_dt2, sims_p5_c5_s5_dt2,
                 sims_p5_c5_s1_dt3, sims_p5_c5_s2_dt3, sims_p5_c5_s3_dt3, sims_p5_c5_s4_dt3, sims_p5_c5_s5_dt3)
write_parquet(sims_dt, "simulations/sims_dtFINAL.parquet")


################################################################################################################
# Cleaning between milking periods when milking units remain infectious for longer
################################################################################################################

################################################################################
params_scen5_dt2$decay_clean <- -log(1) 
# Cleaning reduces infectiousness of all units by 0%
# i.e. no effect
################################################################################
# Simulations for decay time for quarter original 
# with cleaning between 5th and 1st cohort at 0%

test1_0 <- run_multi_sims(df_fields = df_fields5_1,
                          df_milk = df_milk,
                          params = params_scen5_dt2,
                          n_days=200,
                          N_sims = 100,
                          cleaning_time = 1)
test2_0 <- run_multi_sims(df_fields = df_fields5_2,
                          df_milk = df_milk,
                          params = params_scen5_dt2,
                          n_days=200,
                          N_sims = 100,
                          cleaning_time = 1)
test3_0 <- run_multi_sims(df_fields = df_fields5_3,
                          df_milk = df_milk,
                          params = params_scen5_dt2,
                          n_days=200,
                          N_sims = 100,
                          cleaning_time = 1)

test4_0 <- run_multi_sims(df_fields = df_fields5_4,
                          df_milk = df_milk,
                          params = params_scen5_dt2,
                          n_days=200,
                          N_sims = 100,
                          cleaning_time = 1)

test5_0 <- run_multi_sims(df_fields = df_fields5_5,
                          df_milk = df_milk,
                          params = params_scen5_dt2,
                          n_days=200,
                          N_sims = 100,
                          cleaning_time = 1)


################################################################################
params_scen5_dt2$decay_clean <- -log(0.5) 
# Cleaning reduces infectiousness of all units by 50%
##########################################################################################
# Simulations for decay time for quarter original 
# with cleaning between 5th and 1st cohort at 50%

test1_50 <- run_multi_sims(df_fields = df_fields5_1,
                          df_milk = df_milk,
                          params = params_scen5_dt2,
                          n_days=200,
                          N_sims = 100,
                          cleaning_time = 1)
test2_50 <- run_multi_sims(df_fields = df_fields5_2,
                          df_milk = df_milk,
                          params = params_scen5_dt2,
                          n_days=200,
                          N_sims = 100,
                          cleaning_time = 1)
test3_50 <- run_multi_sims(df_fields = df_fields5_3,
                          df_milk = df_milk,
                          params = params_scen5_dt2,
                          n_days=200,
                          N_sims = 100,
                          cleaning_time = 1)

test4_50 <- run_multi_sims(df_fields = df_fields5_4,
                          df_milk = df_milk,
                          params = params_scen5_dt2,
                          n_days=200,
                          N_sims = 100,
                          cleaning_time = 1)

test5_50 <- run_multi_sims(df_fields = df_fields5_5,
                          df_milk = df_milk,
                          params = params_scen5_dt2,
                          n_days=200,
                          N_sims = 100,
                          cleaning_time = 1)

################################################################################
params_scen5_dt2$decay_clean <- -log(0.1) 
# Cleaning reduces infectiousness of all units by 90%
##########################################################################################
# Simulations for decay time for quarter original 
# with cleaning between 5th and 1st cohort at 90%

test1_90 <- run_multi_sims(df_fields = df_fields5_1,
                           df_milk = df_milk,
                           params = params_scen5_dt2,
                           n_days=200,
                           N_sims = 100,
                           cleaning_time = 1)
test2_90 <- run_multi_sims(df_fields = df_fields5_2,
                           df_milk = df_milk,
                           params = params_scen5_dt2,
                           n_days=200,
                           N_sims = 100,
                           cleaning_time = 1)
test3_90 <- run_multi_sims(df_fields = df_fields5_3,
                           df_milk = df_milk,
                           params = params_scen5_dt2,
                           n_days=200,
                           N_sims = 100,
                           cleaning_time = 1)

test4_90 <- run_multi_sims(df_fields = df_fields5_4,
                           df_milk = df_milk,
                           params = params_scen5_dt2,
                           n_days=200,
                           N_sims = 100,
                           cleaning_time = 1)

test5_90 <- run_multi_sims(df_fields = df_fields5_5,
                           df_milk = df_milk,
                           params = params_scen5_dt2,
                           n_days=200,
                           N_sims = 100,
                           cleaning_time = 1)

################################################################################
params_scen5_dt2$decay_clean <- -log(0.01) 
# Cleaning reduces infectiousness of all units by 99%
##########################################################################################
# Simulations for decay time for quarter original 
# with cleaning between 5th and 1st cohort at 99%

test1_99 <- run_multi_sims(df_fields = df_fields5_1,
                            df_milk = df_milk,
                            params = params_scen5_dt2,
                            n_days=200,
                            N_sims = 100,
                            cleaning_time = 1)
test2_99 <- run_multi_sims(df_fields = df_fields5_2,
                            df_milk = df_milk,
                            params = params_scen5_dt2,
                            n_days=200,
                            N_sims = 100,
                            cleaning_time = 1)
test3_99 <- run_multi_sims(df_fields = df_fields5_3,
                            df_milk = df_milk,
                            params = params_scen5_dt2,
                            n_days=200,
                            N_sims = 100,
                            cleaning_time = 1)

test4_99 <- run_multi_sims(df_fields = df_fields5_4,
                            df_milk = df_milk,
                            params = params_scen5_dt2,
                            n_days=200,
                            N_sims = 100,
                            cleaning_time = 1)

test5_99 <- run_multi_sims(df_fields = df_fields5_5,
                            df_milk = df_milk,
                            params = params_scen5_dt2,
                            n_days=200,
                            N_sims = 100,
                            cleaning_time = 1)

################################################################################
params_scen5_dt2$decay_clean <- -log(0.001) 
# Cleaning reduces infectiousness of all units by 99.9%
##########################################################################################
# Simulations for decay time for quarter original 
# with cleaning between 5th and 1st cohort at 99.9%

test1_999 <- run_multi_sims(df_fields = df_fields5_1,
                           df_milk = df_milk,
                           params = params_scen5_dt2,
                           n_days=200,
                           N_sims = 100,
                           cleaning_time = 1)
test2_999 <- run_multi_sims(df_fields = df_fields5_2,
                           df_milk = df_milk,
                           params = params_scen5_dt2,
                           n_days=200,
                           N_sims = 100,
                           cleaning_time = 1)
test3_999 <- run_multi_sims(df_fields = df_fields5_3,
                           df_milk = df_milk,
                           params = params_scen5_dt2,
                           n_days=200,
                           N_sims = 100,
                           cleaning_time = 1)

test4_999 <- run_multi_sims(df_fields = df_fields5_4,
                           df_milk = df_milk,
                           params = params_scen5_dt2,
                           n_days=200,
                           N_sims = 100,
                           cleaning_time = 1)

test5_999 <- run_multi_sims(df_fields = df_fields5_5,
                           df_milk = df_milk,
                           params = params_scen5_dt2,
                           n_days=200,
                           N_sims = 100,
                           cleaning_time = 1)

################################################################################
test1_0$seed_cohorts <- 1
test2_0$seed_cohorts <- 2
test3_0$seed_cohorts <- 3
test4_0$seed_cohorts <- 4
test5_0$seed_cohorts <- 5

test1_50$seed_cohorts <- 1
test2_50$seed_cohorts <- 2
test3_50$seed_cohorts <- 3
test4_50$seed_cohorts <- 4
test5_50$seed_cohorts <- 5

test1_90$seed_cohorts <- 1
test2_90$seed_cohorts <- 2
test3_90$seed_cohorts <- 3
test4_90$seed_cohorts <- 4
test5_90$seed_cohorts <- 5

test1_99$seed_cohorts <- 1
test2_99$seed_cohorts <- 2
test3_99$seed_cohorts <- 3
test4_99$seed_cohorts <- 4
test5_99$seed_cohorts <- 5

test1_999$seed_cohorts <- 1
test2_999$seed_cohorts <- 2
test3_999$seed_cohorts <- 3
test4_999$seed_cohorts <- 4
test5_999$seed_cohorts <- 5

test0 <- rbind(test1_0,
               test2_0,
               test3_0,
               test4_0,
               test5_0)

test0$cleaning_effect <- "0%"

test50 <- rbind(test1_50,
               test2_50,
               test3_50,
               test4_50,
               test5_50)

test50$cleaning_effect <- "50%"

test90 <- rbind(test1_90,
                test2_90,
                test3_90,
                test4_90,
                test5_90)

test90$cleaning_effect <- "90%"
test99 <- rbind(test1_99,
                test2_99,
                test3_99,
                test4_99,
                test5_99)
test99$cleaning_effect <- "99%"

test999 <- rbind(test1_999,
                test2_999,
                test3_999,
                test4_999,
                test5_999)
test999$cleaning_effect <- "99.9%"

test <- rbind(test0, test50, test90, test99, test999)
################################################################################

sims_p5_c5_s1_dt2_clb1_99$seed_cohorts <- 1
sims_p5_c5_s2_dt2_clb1_99$seed_cohorts <- 2
sims_p5_c5_s3_dt2_clb1_99$seed_cohorts <- 3
sims_p5_c5_s4_dt2_clb1_99$seed_cohorts <- 4
sims_p5_c5_s5_dt2_clb1_99$seed_cohorts <- 5

sims_p5_c5_s1_dt2_clb1_95$seed_cohorts <- 1
sims_p5_c5_s2_dt2_clb1_95$seed_cohorts <- 2
sims_p5_c5_s3_dt2_clb1_95$seed_cohorts <- 3
sims_p5_c5_s4_dt2_clb1_95$seed_cohorts <- 4
sims_p5_c5_s5_dt2_clb1_95$seed_cohorts <- 5

sims_p5_c5_s1_dt2_clb1_90$seed_cohorts <- 1
sims_p5_c5_s2_dt2_clb1_90$seed_cohorts <- 2
sims_p5_c5_s3_dt2_clb1_90$seed_cohorts <- 3
sims_p5_c5_s4_dt2_clb1_90$seed_cohorts <- 4
sims_p5_c5_s5_dt2_clb1_90$seed_cohorts <- 5

sims_p5_c5_dt2_clb1_90 <- rbind(sims_p5_c5_s1_dt2_clb1_90,
                                sims_p5_c5_s2_dt2_clb1_90,
                                sims_p5_c5_s3_dt2_clb1_90,
                                sims_p5_c5_s4_dt2_clb1_90,
                                sims_p5_c5_s5_dt2_clb1_90)

sims_p5_c5_dt2_clb1_90$cleaning_effect <-"90%"

sims_p5_c5_dt2_clb1_95 <- rbind(sims_p5_c5_s1_dt2_clb1_95,
                                sims_p5_c5_s2_dt2_clb1_95,
                                sims_p5_c5_s3_dt2_clb1_95,
                                sims_p5_c5_s4_dt2_clb1_95,
                                sims_p5_c5_s5_dt2_clb1_95)

sims_p5_c5_dt2_clb1_95$cleaning_effect <-"95%"

sims_p5_c5_dt2_clb1_99 <- rbind(sims_p5_c5_s1_dt2_clb1_99,
                                sims_p5_c5_s2_dt2_clb1_99,
                                sims_p5_c5_s3_dt2_clb1_99,
                                sims_p5_c5_s4_dt2_clb1_99,
                                sims_p5_c5_s5_dt2_clb1_99)

sims_p5_c5_dt2_clb1_99$cleaning_effect <-"99%"

sims_p5_c5_dt2_clb1 <- rbind(sims_p5_c5_dt2_clb1_90,
                             sims_p5_c5_dt2_clb1_95,
                             sims_p5_c5_dt2_clb1_99)

sims_p5_c5_dt2_clb1$param_set <- 5
sims_p5_c5_dt2_clb1$decay_time<- 2

write_parquet(test, "simulations/sims_dt_cleaningFINAL.parquet")
