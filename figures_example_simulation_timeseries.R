setwd("C:/Users/EALESO/OneDrive - The University of Melbourne/Projects/H5N1 cow transmission")


# 1. Code for plotting figure of example simulations for the transmission regimes

################################################################################
library(ggplot2)
library(patchwork)
library(data.table)
library(arrow)
library(dplyr)
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

# Functions for getting detailed outputs of simulations (timing of all events)
source("R/g_simulation_period_detailed.R")
source("R/milking_period_detailed.R")
source("R/simulation_detailed.R")
source('R/milk_infectious_profile.R')

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
# Single cohort of 1000, with 10 initial infections
df_fields1 <- initialise_fields(1000, A=1000) 
df_fields1$N_S[1] <- df_fields1$N_S[1] - 10
df_fields1$N_E[1] <- 10

################################################################################
# Intialise milking units for main scenarios
df_milk <- initialise_milking_machines(50)

##########################################################################################################
# Run simulations to get example outcomes
##########################################################################################################

set.seed(1234)
options(dplyr.summarise.inform = FALSE)
det_sim1 <- simulation_detailed(df_fields=df_fields1,
                                df_milk = df_milk,
                                params = params_scen1,
                                n_days=100, 
                                milking_frequency = 2)
saveRDS(det_sim1, "simulations/example_traj1.rds")

set.seed(1234)
options(dplyr.summarise.inform = FALSE)
det_sim3 <- simulation_detailed(df_fields=df_fields1,
                                df_milk = df_milk,
                                params = params_scen3,
                                n_days=100, 
                                milking_frequency = 2)
saveRDS(det_sim3, "simulations/example_traj3.rds")

set.seed(1234)
options(dplyr.summarise.inform = FALSE)
det_sim5 <- simulation_detailed(df_fields=df_fields1,
                                df_milk = df_milk,
                                params = params_scen5,
                                n_days=100, 
                                milking_frequency = 2)
saveRDS(det_sim5, "simulations/example_traj5.rds")


##########################################################################################################
# Reformat outputs ready for plotting
##########################################################################################################
# Reload if previously ran code
det_sim1 <- readRDS("simulations/example_traj1.rds")
det_sim3 <- readRDS("simulations/example_traj3.rds")
det_sim5 <- readRDS("simulations/example_traj5.rds")

################################################################################
# Assign objects from simulation_detailed function call

event_log1 <- det_sim1[[1]]
event_log3 <- det_sim3[[1]]
event_log5 <- det_sim5[[1]]

milk_log1 <- det_sim1[[3]]
milk_log3 <- det_sim3[[3]]
milk_log5 <- det_sim5[[3]]

milk_log1 <- milk_infectious_profile(milk_log = milk_log1, times=seq(0,100,0.01), params=params_scen1)
milk_log3 <- milk_infectious_profile(milk_log = milk_log3, times=seq(0,100,0.01), params=params_scen3)
milk_log5 <- milk_infectious_profile(milk_log = milk_log5, times=seq(0,100,0.01), params=params_scen5)

################################################################################
# Reformat event_log1
event_log1_duplicate <- event_log1
event_log1_duplicate$time <- c(0, event_log1$time[1:(nrow(event_log1)-1)])
event_log1 <- rbind(event_log1,event_log1_duplicate)
event_log1 <- event_log1[order(event_log1$time),]

################################################################################
# Reformat event_log3
event_log3_duplicate <- event_log3
event_log3_duplicate$time <- c(0, event_log3$time[1:(nrow(event_log3)-1)])
event_log3 <- rbind(event_log3,event_log3_duplicate)
event_log3 <- event_log3[order(event_log3$time),]

################################################################################
# Reformat event_log5
event_log5_duplicate <- event_log5
event_log5_duplicate$time <- c(0, event_log5$time[1:(nrow(event_log5)-1)])
event_log5 <- rbind(event_log5,event_log5_duplicate)
event_log5 <- event_log5[order(event_log5$time),]


################################################################################
# Merge for plotting
event_log1$param <- "Cow\u2013cow transmission only"
event_log3$param <- "Mixed mode transmission"
event_log5$param <- "Cow\u2013milking unit transmission only"

event_log <- rbind(event_log1, event_log3, event_log5)
event_log$param <- factor(event_log$param, levels = c("Cow\u2013cow transmission only",
                                                      "Mixed mode transmission",
                                                      "Cow\u2013milking unit transmission only"))

milk_log1$param <- "Cow\u2013cow transmission only"
milk_log3$param <- "Mixed mode transmission"
milk_log5$param <- "Cow\u2013milking unit transmission only"

milk_log <- rbind(milk_log1, milk_log3, milk_log5)
milk_log$param <- factor(milk_log$param, levels = c("Cow\u2013cow transmission only",
                                                      "Mixed mode transmission",
                                                      "Cow\u2013milking unit transmission only"))

milk_times <- data.frame(xmin=(6/24)+seq(0,100,0.5),
                         xmax=(6/24)+seq(0,100,0.5)+(20*5/(60*24) ) )
##########################################################################################################
# Plot figures
##########################################################################################################
cols <- RColorBrewer::brewer.pal(8, "Dark2")


plt1<-ggplot(event_log)+
  geom_line(aes(x=time, y=N_E+N_I, color=param))+
  ylab("Infection prevalence")+
  xlab("Time (days)")+
  scale_color_brewer("Transmission regime",palette="Dark2")+
  annotate("rect", xmin=23, xmax=26,ymin=170,ymax=242, color="black", alpha=0.,linetype="dashed")+
  scale_x_continuous(limits=c(0,100),expand=expansion(0))+
  theme_bw(base_size = 14)+
  theme(legend.background = element_rect(color="black"),
        legend.position = c(0.69,0.88))


plt2 <- ggplot(event_log[event_log$time>23&event_log$time <26,])+
  geom_line(aes(x=time, y=N_I+N_E, color=param))+
  ylab("Infection prevalence")+
  xlab("Time (days)")+
  scale_x_continuous(limits=c(23,26),expand=expansion(0))+
  scale_color_brewer("Tranmission regime",palette="Dark2")+
  geom_rect(data=milk_times, aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf), fill="grey",alpha=0.4)+
  #scale_y_continuous(limits=c(0,250),expand=expansion(0))+
  theme_bw(base_size = 14)+
  theme(legend.position = "none")


plt3 <- ggplot(milk_log[milk_log$time>23&milk_log$time <26, ])+
  geom_line(aes(x=time, y=inf_prof/50, color=param))+
  ylab("Average infectiousness\nof milking units")+
  xlab("Time (days)")+
  scale_x_continuous(limits=c(23,26),expand=expansion(0))+
  scale_color_brewer("Tranmission regime",palette="Dark2")+
  geom_rect(data=milk_times, aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf), fill="grey",alpha=0.4)+
  #scale_y_continuous(limits=c(0,250),expand=expansion(0))+
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

plt1 <- plt1+labs(tag="A")+theme(plot.tag.position = c(0.01,0.98))
plt2 <- plt2+labs(tag="B")+theme(plot.tag.position = c(0.01,0.98))
plt3 <- plt3+labs(tag="C")+theme(plot.tag.position = c(0.01,0.98))

plt2 <- plt2+theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank())

plt23 <- plt2+plt3 +plot_layout(nrow=2)
plt1+plt23 +plot_layout(widths=c(0.6,1))

ggsave('figures/example_traj.png', width=15.5, height=6.5)
