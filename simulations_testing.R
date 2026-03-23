setwd("C:/Users/EALESO/OneDrive - The University of Melbourne/Projects/H5N1 cow transmission")


# 1. Simulations to investigate effect of testing

################################################################################

library(patchwork)
library(data.table)
library(arrow)
library(rstan)
library(dplyr)
library(ggplot2)
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

# Functions for getting detailed outputs of simulations (timing of all infections)
source("R/g_simulation_period_detailed_infections_only.R")
source("R/milking_period_detailed_infections_only.R")
source("R/simulation_detailed_infections_only.R")
source("R/run_multi_sims_infections_only.R")




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

################################################################################
# Intialise milking units for scenarios
df_milk1 <- initialise_milking_machines(50)
df_milk2 <- initialise_milking_machines(5)
df_milk3 <- initialise_milking_machines(500)


##########################################################################################################
# Running 200 days of simulations getting exact timing of infections
##########################################################################################################

set.seed(12345)
options(dplyr.summarise.inform = FALSE)
inf_sims1 <- run_multi_sims_infections_only(df_fields = df_fields1,
                                            df_milk = df_milk1,
                                            params = params_scen5,
                                            n_days=200,
                                            N_sims = 100)
write_parquet(inf_sims1, "simulations/sims_testing1.parquet")

set.seed(12345)
options(dplyr.summarise.inform = FALSE)
inf_sims2 <- run_multi_sims_infections_only(df_fields = df_fields2,
                                            df_milk = df_milk2,
                                            params = params_scen5,
                                            n_days=200,
                                            N_sims = 100)
write_parquet(inf_sims2, "simulations/sims_testing2.parquet")


set.seed(12345)
options(dplyr.summarise.inform = FALSE)
inf_sims3 <- run_multi_sims_infections_only(df_fields = df_fields3,
                                            df_milk = df_milk3,
                                            params = params_scen5,
                                            n_days=200,
                                            N_sims = 100)
write_parquet(inf_sims3, "simulations/sims_testing3.parquet")

