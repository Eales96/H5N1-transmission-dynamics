setwd("C:/Users/EALESO/OneDrive - The University of Melbourne/Projects/H5N1 cow transmission")


# 1. Simulations to investigate effect of reactive cohorting

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
source("R/reactive_cohort.R")
source("R/simulation_reactive.R")
source("R/run_multi_sims_reactive.R")

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
params_scen3 <- params
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
# Single cohort of 1000, with 1 initial infections
df_fields1 <- initialise_fields(1000, A=1000) 
df_fields1$N_S[1] <- df_fields1$N_S[1] - 1
df_fields1$N_E[1] <- 1

################################################################################
# Single cohort of 100, with 1 initial infections
df_fields2 <- initialise_fields(100, A=100) 
df_fields2$N_S[1] <- df_fields2$N_S[1] - 1
df_fields2$N_E[1] <- 1

################################################################################
# Single cohort of 10000, with 1 initial infections
df_fields3 <- initialise_fields(10000, A=10000) 
df_fields3$N_S[1] <- df_fields3$N_S[1] - 1
df_fields3$N_E[1] <- 1


################################################################################
# Intialise milking units for main scenarios
df_milk1 <- initialise_milking_machines(50)
df_milk2 <- initialise_milking_machines(5)
df_milk3 <- initialise_milking_machines(500)

##########################################################################################################
# Running simulations for scenarios comparing timing of reactive cohorting
##########################################################################################################

################ Simulations for field with 1000 cattle ########################
set.seed(12345)
sims_fm1_p5_r1 <- run_multi_sims_reactive(df_fields = df_fields1,
                                          df_milk = df_milk1,
                                          params = params_scen5,
                                          n_days=400,
                                          N_sims = 200,
                                          reaction_threshold = 1,
                                          N_fields=5)

set.seed(12345)
sims_fm1_p5_r5 <- run_multi_sims_reactive(df_fields = df_fields1,
                                          df_milk = df_milk1,
                                          params = params_scen5,
                                          n_days=400,
                                          N_sims = 200,
                                          reaction_threshold = 5,
                                          N_fields=5)

set.seed(12345)
sims_fm1_p5_r10 <- run_multi_sims_reactive(df_fields = df_fields1,
                                           df_milk = df_milk1,
                                           params = params_scen5,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 10,
                                           N_fields=5)

set.seed(12345)
sims_fm1_p5_r20 <- run_multi_sims_reactive(df_fields = df_fields1,
                                           df_milk = df_milk1,
                                           params = params_scen5,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 20,
                                           N_fields=5)

set.seed(12345)
sims_fm1_p5_r50 <- run_multi_sims_reactive(df_fields = df_fields1,
                                           df_milk = df_milk1,
                                           params = params_scen5,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 50,
                                           N_fields=5)

set.seed(12345)
sims_fm1_p5_r100000 <- run_multi_sims_reactive(df_fields = df_fields1,
                                           df_milk = df_milk1,
                                           params = params_scen5,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 100000,
                                           N_fields=5)

################ Simulations for field with 100 cattle ########################
set.seed(12345)
sims_fm2_p5_r1 <- run_multi_sims_reactive(df_fields = df_fields2,
                                          df_milk = df_milk2,
                                          params = params_scen5,
                                          n_days=200,
                                          N_sims = 200,
                                          reaction_threshold = 1,
                                          N_fields=5)

set.seed(12345)
sims_fm2_p5_r5 <- run_multi_sims_reactive(df_fields = df_fields2,
                                          df_milk = df_milk2,
                                          params = params_scen5,
                                          n_days=200,
                                          N_sims = 200,
                                          reaction_threshold = 5,
                                          N_fields=5)

set.seed(12345)
sims_fm2_p5_r10 <- run_multi_sims_reactive(df_fields = df_fields2,
                                           df_milk = df_milk2,
                                           params = params_scen5,
                                           n_days=200,
                                           N_sims = 200,
                                           reaction_threshold = 10,
                                           N_fields=5)

set.seed(12345)
sims_fm2_p5_r20 <- run_multi_sims_reactive(df_fields = df_fields2,
                                           df_milk = df_milk2,
                                           params = params_scen5,
                                           n_days=200,
                                           N_sims = 200,
                                           reaction_threshold = 20,
                                           N_fields=5)

set.seed(12345)
sims_fm2_p5_r50 <- run_multi_sims_reactive(df_fields = df_fields2,
                                           df_milk = df_milk2,
                                           params = params_scen5,
                                           n_days=200,
                                           N_sims = 200,
                                           reaction_threshold = 50,
                                           N_fields=5)

set.seed(12345)
sims_fm2_p5_r100000 <- run_multi_sims_reactive(df_fields = df_fields2,
                                               df_milk = df_milk2,
                                               params = params_scen5,
                                               n_days=200,
                                               N_sims = 200,
                                               reaction_threshold = 100000,
                                               N_fields=5)



################ Simulations for field with 10000 cattle ########################
set.seed(12345)
sims_fm3_p5_r1 <- run_multi_sims_reactive(df_fields = df_fields3,
                                          df_milk = df_milk3,
                                          params = params_scen5,
                                          n_days=400,
                                          N_sims = 200,
                                          reaction_threshold = 1,
                                          N_fields=5)

set.seed(12345)
sims_fm3_p5_r5 <- run_multi_sims_reactive(df_fields = df_fields3,
                                          df_milk = df_milk3,
                                          params = params_scen5,
                                          n_days=400,
                                          N_sims = 200,
                                          reaction_threshold = 5,
                                          N_fields=5)

set.seed(12345)
sims_fm3_p5_r10 <- run_multi_sims_reactive(df_fields = df_fields3,
                                           df_milk = df_milk3,
                                           params = params_scen5,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 10,
                                           N_fields=5)

set.seed(12345)
sims_fm3_p5_r20 <- run_multi_sims_reactive(df_fields = df_fields3,
                                           df_milk = df_milk3,
                                           params = params_scen5,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 20,
                                           N_fields=5)

set.seed(12345)
sims_fm3_p5_r50 <- run_multi_sims_reactive(df_fields = df_fields3,
                                           df_milk = df_milk3,
                                           params = params_scen5,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 50,
                                           N_fields=5)

set.seed(12345)
sims_fm3_p5_r100000 <- run_multi_sims_reactive(df_fields = df_fields3,
                                               df_milk = df_milk3,
                                               params = params_scen5,
                                               n_days=400,
                                               N_sims = 200,
                                               reaction_threshold = 100000,
                                               N_fields=5)




###########################################################################################################
##################### Saving outcomes of simulations ###########################
###########################################################################################################

################################################################################
# 1000 cows in field
sims_fm1_p5_r1$reaction <- 1
sims_fm1_p5_r5$reaction <- 5
sims_fm1_p5_r10$reaction <- 10
sims_fm1_p5_r20$reaction <- 20
sims_fm1_p5_r50$reaction <- 50
sims_fm1_p5_r100000$reaction <- 100000

sims_fm1_p5 <- rbind(sims_fm1_p5_r1,
                     sims_fm1_p5_r5,
                     sims_fm1_p5_r10,
                     sims_fm1_p5_r20,
                     sims_fm1_p5_r50,
                     sims_fm1_p5_r100000)

sims_fm1_p5$N_cows <- 1000
sims_fm1_p5$param_set <- 5


################################################################################
# 100 cows in field

sims_fm2_p5_r1$reaction <- 1
sims_fm2_p5_r5$reaction <- 5
sims_fm2_p5_r10$reaction <- 10
sims_fm2_p5_r20$reaction <- 20
sims_fm2_p5_r50$reaction <- 50
sims_fm2_p5_r100000$reaction <- 100000

sims_fm2_p5 <- rbind(sims_fm2_p5_r1,
                     sims_fm2_p5_r5,
                     sims_fm2_p5_r10,
                     sims_fm2_p5_r20,
                     sims_fm2_p5_r50,
                     sims_fm2_p5_r100000)

sims_fm2_p5$N_cows <- 100
sims_fm2_p5$param_set <- 5


################################################################################
# 10000 cows in field
sims_fm3_p5_r1$reaction <- 1
sims_fm3_p5_r5$reaction <- 5
sims_fm3_p5_r10$reaction <- 10
sims_fm3_p5_r20$reaction <- 20
sims_fm3_p5_r50$reaction <- 50
sims_fm3_p5_r100000$reaction <- 100000

sims_fm3_p5 <- rbind(sims_fm3_p5_r1,
                     sims_fm3_p5_r5,
                     sims_fm3_p5_r10,
                     sims_fm3_p5_r20,
                     sims_fm3_p5_r50,
                     sims_fm3_p5_r100000)

sims_fm3_p5$N_cows <- 10000
sims_fm3_p5$param_set <- 5

################################################################################

sims_p5 <- rbind(sims_fm1_p5,
                 sims_fm2_p5,
                 sims_fm3_p5)

write_parquet(sims_p5, "simulations/sims_reactive_p5.parquet")


##########################################################################################################
# Running simulations (Parameter scenario 3)
##########################################################################################################

################ Simulations for field with 1000 cattle ########################
set.seed(12345)
sims_fm1_p3_r1 <- run_multi_sims_reactive(df_fields = df_fields1,
                                          df_milk = df_milk1,
                                          params = params_scen3,
                                          n_days=400,
                                          N_sims = 200,
                                          reaction_threshold = 1,
                                          N_fields=5)

set.seed(12345)
sims_fm1_p3_r5 <- run_multi_sims_reactive(df_fields = df_fields1,
                                          df_milk = df_milk1,
                                          params = params_scen3,
                                          n_days=400,
                                          N_sims = 200,
                                          reaction_threshold = 5,
                                          N_fields=5)

set.seed(12345)
sims_fm1_p3_r10 <- run_multi_sims_reactive(df_fields = df_fields1,
                                           df_milk = df_milk1,
                                           params = params_scen3,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 10,
                                           N_fields=5)

set.seed(12345)
sims_fm1_p3_r20 <- run_multi_sims_reactive(df_fields = df_fields1,
                                           df_milk = df_milk1,
                                           params = params_scen3,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 20,
                                           N_fields=5)

set.seed(12345)
sims_fm1_p3_r50 <- run_multi_sims_reactive(df_fields = df_fields1,
                                           df_milk = df_milk1,
                                           params = params_scen3,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 50,
                                           N_fields=5)

set.seed(12345)
sims_fm1_p3_r100000 <- run_multi_sims_reactive(df_fields = df_fields1,
                                               df_milk = df_milk1,
                                               params = params_scen3,
                                               n_days=400,
                                               N_sims = 200,
                                               reaction_threshold = 100000,
                                               N_fields=5)

################ Simulations for field with 100 cattle ########################
set.seed(12345)
sims_fm2_p3_r1 <- run_multi_sims_reactive(df_fields = df_fields2,
                                          df_milk = df_milk2,
                                          params = params_scen3,
                                          n_days=200,
                                          N_sims = 200,
                                          reaction_threshold = 1,
                                          N_fields=5)

set.seed(12345)
sims_fm2_p3_r5 <- run_multi_sims_reactive(df_fields = df_fields2,
                                          df_milk = df_milk2,
                                          params = params_scen3,
                                          n_days=200,
                                          N_sims = 200,
                                          reaction_threshold = 5,
                                          N_fields=5)

set.seed(12345)
sims_fm2_p3_r10 <- run_multi_sims_reactive(df_fields = df_fields2,
                                           df_milk = df_milk2,
                                           params = params_scen3,
                                           n_days=200,
                                           N_sims = 200,
                                           reaction_threshold = 10,
                                           N_fields=5)

set.seed(12345)
sims_fm2_p3_r20 <- run_multi_sims_reactive(df_fields = df_fields2,
                                           df_milk = df_milk2,
                                           params = params_scen3,
                                           n_days=200,
                                           N_sims = 200,
                                           reaction_threshold = 20,
                                           N_fields=5)

set.seed(12345)
sims_fm2_p3_r50 <- run_multi_sims_reactive(df_fields = df_fields2,
                                           df_milk = df_milk2,
                                           params = params_scen3,
                                           n_days=200,
                                           N_sims = 200,
                                           reaction_threshold = 50,
                                           N_fields=5)

set.seed(12345)
sims_fm2_p3_r100000 <- run_multi_sims_reactive(df_fields = df_fields2,
                                               df_milk = df_milk2,
                                               params = params_scen3,
                                               n_days=200,
                                               N_sims = 200,
                                               reaction_threshold = 100000,
                                               N_fields=5)



################ Simulations for field with 10000 cattle ########################
set.seed(12345)
sims_fm3_p3_r1 <- run_multi_sims_reactive(df_fields = df_fields3,
                                          df_milk = df_milk3,
                                          params = params_scen3,
                                          n_days=400,
                                          N_sims = 200,
                                          reaction_threshold = 1,
                                          N_fields=5)

set.seed(12345)
sims_fm3_p3_r5 <- run_multi_sims_reactive(df_fields = df_fields3,
                                          df_milk = df_milk3,
                                          params = params_scen3,
                                          n_days=200,
                                          N_sims = 200,
                                          reaction_threshold = 5,
                                          N_fields=5)

set.seed(12345)
sims_fm3_p3_r10 <- run_multi_sims_reactive(df_fields = df_fields3,
                                           df_milk = df_milk3,
                                           params = params_scen3,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 10,
                                           N_fields=5)

set.seed(12345)
sims_fm3_p3_r20 <- run_multi_sims_reactive(df_fields = df_fields3,
                                           df_milk = df_milk3,
                                           params = params_scen3,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 20,
                                           N_fields=5)

set.seed(12345)
sims_fm3_p3_r50 <- run_multi_sims_reactive(df_fields = df_fields3,
                                           df_milk = df_milk3,
                                           params = params_scen3,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 50,
                                           N_fields=5)

set.seed(12345)
sims_fm3_p3_r100000 <- run_multi_sims_reactive(df_fields = df_fields3,
                                               df_milk = df_milk3,
                                               params = params_scen3,
                                               n_days=400,
                                               N_sims = 200,
                                               reaction_threshold = 100000,
                                               N_fields=5)




###########################################################################################################
##################### Saving outcomes of simulations ###########################
###########################################################################################################

################################################################################
# 1000 cows in field
sims_fm1_p3_r1$reaction <- 1
sims_fm1_p3_r5$reaction <- 5
sims_fm1_p3_r10$reaction <- 10
sims_fm1_p3_r20$reaction <- 20
sims_fm1_p3_r50$reaction <- 50
sims_fm1_p3_r100000$reaction <- 100000

sims_fm1_p3 <- rbind(sims_fm1_p3_r1,
                     sims_fm1_p3_r5,
                     sims_fm1_p3_r10,
                     sims_fm1_p3_r20,
                     sims_fm1_p3_r50,
                     sims_fm1_p3_r100000)

sims_fm1_p3$N_cows <- 1000
sims_fm1_p3$param_set <- 3


################################################################################
# 100 cows in field

sims_fm2_p3_r1$reaction <- 1
sims_fm2_p3_r5$reaction <- 5
sims_fm2_p3_r10$reaction <- 10
sims_fm2_p3_r20$reaction <- 20
sims_fm2_p3_r50$reaction <- 50
sims_fm2_p3_r100000$reaction <- 100000

sims_fm2_p3 <- rbind(sims_fm2_p3_r1,
                     sims_fm2_p3_r5,
                     sims_fm2_p3_r10,
                     sims_fm2_p3_r20,
                     sims_fm2_p3_r50,
                     sims_fm2_p3_r100000)

sims_fm2_p3$N_cows <- 100
sims_fm2_p3$param_set <- 3


################################################################################
# 10000 cows in field
sims_fm3_p3_r1$reaction <- 1
sims_fm3_p3_r5$reaction <- 5
sims_fm3_p3_r10$reaction <- 10
sims_fm3_p3_r20$reaction <- 20
sims_fm3_p3_r50$reaction <- 50
sims_fm3_p3_r100000$reaction <- 100000

sims_fm3_p3 <- rbind(sims_fm3_p3_r1,
                     sims_fm3_p3_r5,
                     sims_fm3_p3_r10,
                     sims_fm3_p3_r20,
                     sims_fm3_p3_r50,
                     sims_fm3_p3_r100000)

sims_fm3_p3$N_cows <- 10000
sims_fm3_p3$param_set <- 3

################################################################################

sims_p3 <- rbind(sims_fm1_p3,
                 sims_fm2_p3,
                 sims_fm3_p3)

write_parquet(sims_p3, "simulations/sims_reactive_p3.parquet")


##########################################################################################################
# Running simulations (parameter scenario 1)
##########################################################################################################

################ Simulations for field with 1000 cattle ########################
set.seed(12345)
sims_fm1_p1_r1 <- run_multi_sims_reactive(df_fields = df_fields1,
                                          df_milk = df_milk1,
                                          params = params_scen1,
                                          n_days=400,
                                          N_sims = 200,
                                          reaction_threshold = 1,
                                          N_fields=5)

set.seed(12345)
sims_fm1_p1_r5 <- run_multi_sims_reactive(df_fields = df_fields1,
                                          df_milk = df_milk1,
                                          params = params_scen1,
                                          n_days=400,
                                          N_sims = 200,
                                          reaction_threshold = 5,
                                          N_fields=5)

set.seed(12345)
sims_fm1_p1_r10 <- run_multi_sims_reactive(df_fields = df_fields1,
                                           df_milk = df_milk1,
                                           params = params_scen1,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 10,
                                           N_fields=5)

set.seed(12345)
sims_fm1_p1_r20 <- run_multi_sims_reactive(df_fields = df_fields1,
                                           df_milk = df_milk1,
                                           params = params_scen1,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 20,
                                           N_fields=5)

set.seed(12345)
sims_fm1_p1_r50 <- run_multi_sims_reactive(df_fields = df_fields1,
                                           df_milk = df_milk1,
                                           params = params_scen1,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 50,
                                           N_fields=5)

set.seed(12345)
sims_fm1_p1_r100000 <- run_multi_sims_reactive(df_fields = df_fields1,
                                               df_milk = df_milk1,
                                               params = params_scen1,
                                               n_days=400,
                                               N_sims = 200,
                                               reaction_threshold = 100000,
                                               N_fields=5)

################ Simulations for field with 100 cattle ########################
set.seed(12345)
sims_fm2_p1_r1 <- run_multi_sims_reactive(df_fields = df_fields2,
                                          df_milk = df_milk2,
                                          params = params_scen1,
                                          n_days=400,
                                          N_sims = 200,
                                          reaction_threshold = 1,
                                          N_fields=5)

set.seed(12345)
sims_fm2_p1_r5 <- run_multi_sims_reactive(df_fields = df_fields2,
                                          df_milk = df_milk2,
                                          params = params_scen1,
                                          n_days=400,
                                          N_sims = 200,
                                          reaction_threshold = 5,
                                          N_fields=5)

set.seed(12345)
sims_fm2_p1_r10 <- run_multi_sims_reactive(df_fields = df_fields2,
                                           df_milk = df_milk2,
                                           params = params_scen1,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 10,
                                           N_fields=5)

set.seed(12345)
sims_fm2_p1_r20 <- run_multi_sims_reactive(df_fields = df_fields2,
                                           df_milk = df_milk2,
                                           params = params_scen1,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 20,
                                           N_fields=5)

set.seed(12345)
sims_fm2_p1_r50 <- run_multi_sims_reactive(df_fields = df_fields2,
                                           df_milk = df_milk2,
                                           params = params_scen1,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 50,
                                           N_fields=5)

set.seed(12345)
sims_fm2_p1_r100000 <- run_multi_sims_reactive(df_fields = df_fields2,
                                               df_milk = df_milk2,
                                               params = params_scen1,
                                               n_days=400,
                                               N_sims = 200,
                                               reaction_threshold = 100000,
                                               N_fields=5)



################ Simulations for field with 10000 cattle ########################
set.seed(12345)
sims_fm3_p1_r1 <- run_multi_sims_reactive(df_fields = df_fields3,
                                          df_milk = df_milk3,
                                          params = params_scen1,
                                          n_days=400,
                                          N_sims = 200,
                                          reaction_threshold = 1,
                                          N_fields=5)

set.seed(12345)
sims_fm3_p1_r5 <- run_multi_sims_reactive(df_fields = df_fields3,
                                          df_milk = df_milk3,
                                          params = params_scen1,
                                          n_days=400,
                                          N_sims = 200,
                                          reaction_threshold = 5,
                                          N_fields=5)

set.seed(12345)
sims_fm3_p1_r10 <- run_multi_sims_reactive(df_fields = df_fields3,
                                           df_milk = df_milk3,
                                           params = params_scen1,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 10,
                                           N_fields=5)

set.seed(12345)
sims_fm3_p1_r20 <- run_multi_sims_reactive(df_fields = df_fields3,
                                           df_milk = df_milk3,
                                           params = params_scen1,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 20,
                                           N_fields=5)

set.seed(12345)
sims_fm3_p1_r50 <- run_multi_sims_reactive(df_fields = df_fields3,
                                           df_milk = df_milk3,
                                           params = params_scen1,
                                           n_days=400,
                                           N_sims = 200,
                                           reaction_threshold = 50,
                                           N_fields=5)

set.seed(12345)
sims_fm3_p1_r100000 <- run_multi_sims_reactive(df_fields = df_fields3,
                                               df_milk = df_milk3,
                                               params = params_scen1,
                                               n_days=400,
                                               N_sims = 200,
                                               reaction_threshold = 100000,
                                               N_fields=5)




###########################################################################################################
##################### Saving outcomes of simulations ###########################
###########################################################################################################

################################################################################
# 1000 cows in field
sims_fm1_p1_r1$reaction <- 1
sims_fm1_p1_r5$reaction <- 5
sims_fm1_p1_r10$reaction <- 10
sims_fm1_p1_r20$reaction <- 20
sims_fm1_p1_r50$reaction <- 50
sims_fm1_p1_r100000$reaction <- 100000

sims_fm1_p1 <- rbind(sims_fm1_p1_r1,
                     sims_fm1_p1_r5,
                     sims_fm1_p1_r10,
                     sims_fm1_p1_r20,
                     sims_fm1_p1_r50,
                     sims_fm1_p1_r100000)

sims_fm1_p1$N_cows <- 1000
sims_fm1_p1$param_set <- 1


################################################################################
# 100 cows in field

sims_fm2_p1_r1$reaction <- 1
sims_fm2_p1_r5$reaction <- 5
sims_fm2_p1_r10$reaction <- 10
sims_fm2_p1_r20$reaction <- 20
sims_fm2_p1_r50$reaction <- 50
sims_fm2_p1_r100000$reaction <- 100000

sims_fm2_p1 <- rbind(sims_fm2_p1_r1,
                     sims_fm2_p1_r5,
                     sims_fm2_p1_r10,
                     sims_fm2_p1_r20,
                     sims_fm2_p1_r50,
                     sims_fm2_p1_r100000)

sims_fm2_p1$N_cows <- 100
sims_fm2_p1$param_set <- 1


################################################################################
# 10000 cows in field
sims_fm3_p1_r1$reaction <- 1
sims_fm3_p1_r5$reaction <- 5
sims_fm3_p1_r10$reaction <- 10
sims_fm3_p1_r20$reaction <- 20
sims_fm3_p1_r50$reaction <- 50
sims_fm3_p1_r100000$reaction <- 100000

sims_fm3_p1 <- rbind(sims_fm3_p1_r1,
                     sims_fm3_p1_r5,
                     sims_fm3_p1_r10,
                     sims_fm3_p1_r20,
                     sims_fm3_p1_r50,
                     sims_fm3_p1_r100000)

sims_fm3_p1$N_cows <- 10000
sims_fm3_p1$param_set <- 1

################################################################################

sims_p1 <- rbind(sims_fm1_p1,
                 sims_fm2_p1,
                 sims_fm3_p1)

write_parquet(sims_p1, "simulations/sims_reactive_p1.parquet")
