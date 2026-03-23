setwd("C:/Users/EALESO/OneDrive - The University of Melbourne/Projects/H5N1 cow transmission")

################################################################################

library(patchwork)
library(data.table)
library(arrow)
library(rstan)
library(dplyr)
library(ggplot2)
################################################################################


# Functions for testing
source('R/get_log_ct_value.R')




################################################################################
# Load fitted stan model for Ct value trajectories (from Eales et al 2026)
mod1 <- readRDS('fitted_stan_models/mod_ft_ADC.rds')
post1 <- rstan::extract(mod1)

# Extract median parameters for modelling test positivity
ct_params <- data.frame(t_peak_mn = median(post1$t_peak_mn),
                        t_peak_sd = median(post1$t_peak_sd),
                        rise_mn = median(post1$rise_mn),
                        rise_sd = median(post1$rise_sd),
                        peak_mn = median(post1$peak_mn),
                        peak_sd = median(post1$peak_sd),
                        decay_mn = median(post1$decay_mn),
                        decay_sd = median(post1$decay_sd),
                        sigma = median(post1$sigma_Ct))


ct_params <- data.frame(t_peak_mn = 1.100091, #(post1$t_peak_mn),
                        t_peak_sd = 0.05049329, #median(post1$t_peak_sd),
                        rise_mn = 77.59923, #median(post1$rise_mn),
                        rise_sd = 0.102828, #median(post1$rise_sd),
                        peak_mn = 15.67295, #median(post1$peak_mn),
                        peak_sd = 2.745692, #median(post1$peak_sd),
                        decay_mn = 0.8849212, #median(post1$decay_mn),
                        decay_sd = 0.4387093, #median(post1$decay_sd),
                        sigma = 2.490441)#median(post1$sigma_Ct))

################################################################################
# Define standard curve for Ct value to log-titre
# Using standard curve for Ct value to log-titre in Facciuolo et al (same used in Eales et al 2026)

log_params <- data.frame(intercept=44.83, 
                         grad=-1.736)


################################################################################
# Read in simulations
inf_sims1 <- read_parquet("simulations/sims_testing1.parquet")
inf_sims2 <- read_parquet("simulations/sims_testing2.parquet")
inf_sims3 <- read_parquet("simulations/sims_testing3.parquet")


# TEMPORARY
inf_sims1 <- inf_sims1[inf_sims1$event%in%c("milk","initial"),]
inf_sims2 <- inf_sims2[inf_sims2$event%in%c("milk","initial"),]
inf_sims3 <- inf_sims3[inf_sims3$event%in%c("milk","initial"),]


organise_duplicates <- function(inf_sims1){
  
  new_df <- data.frame()
  for(i in unique(inf_sims1$sim)){
    print(i)
    temp <- inf_sims1[inf_sims1$sim==i,]
    temp$inc <- 1
    if(nrow(temp)>1){
      temp[2:nrow(temp),]$inc <- temp[1:(nrow(temp)-1),]$N_S-temp[2:nrow(temp),]$N_S
    }
    
    for(j in 1:nrow(temp)){
      for(k in temp[j,]$inc:1){
        new_row <- temp[j,]
        new_row$N_S <- new_row$N_S+k-1
        new_row$N_E <- new_row$N_E-k+1
        new_row$inc <- 1
        new_df <- rbind(new_df, new_row)
      }
    }
    
  }
  
  return(new_df)
  
}

inf_sims1 <- organise_duplicates(inf_sims1)
inf_sims2 <- organise_duplicates(inf_sims2)
inf_sims3 <- organise_duplicates(inf_sims3)

################################################################################
# Get Ct draws for each infection trajectory

set.seed(12345)
ct_sims1 <- get_ct_draws(inf_log = inf_sims1,
                         ct_params = ct_params,
                         log_params = log_params,
                         N_draws = 10,
                         times = seq(0.5,200,1))

ct_sims2 <- get_ct_draws(inf_log = inf_sims2,
                         ct_params = ct_params,
                         log_params = log_params,
                         N_draws = 10,
                         times = seq(0.5,200,1))

ct_sims3 <- get_ct_draws(inf_log = inf_sims3,
                         ct_params = ct_params,
                         log_params = log_params,
                         N_draws = 10,
                         times = seq(0.5,200,1))

################################################################################
# Get the probability for specific Ct and milk parameter sets

ct_sims_prob1 <- get_ct_value(final=ct_sims1,
                              m=0.25,
                              log_params = log_params,
                              ct_threshold = 36)

ct_sims_prob2 <- get_ct_value(final=ct_sims2,
                              m=0.25,
                              log_params = log_params,
                              ct_threshold = 36)

ct_sims_prob3 <- get_ct_value(final=ct_sims3,
                              m=0.25,
                              log_params = log_params,
                              ct_threshold = 36)

################################################################################
# Get the detection times for each scenario

detection_times1 <- extract_all_detect_times(ct_sims_prob = ct_sims_prob1,
                                             freq_test = c(1,2,3,4,5,6,7,10,14,21,28),
                                             sims_per_comb = 10)

detection_times2 <- extract_all_detect_times(ct_sims_prob = ct_sims_prob2,
                                             freq_test = c(1,2,3,4,5,6,7,10,14,21,28),
                                             sims_per_comb = 10)

detection_times3 <- extract_all_detect_times(ct_sims_prob = ct_sims_prob3,
                                             freq_test = c(1,2,3,4,5,6,7,10,14,21,28),
                                             sims_per_comb = 10)



take_off1 <- ct_sims_prob1[ct_sims_prob1$time==max(ct_sims_prob1$time) & ct_sims_prob1$N_inf>10,]$sim
take_off2 <- ct_sims_prob2[ct_sims_prob2$time==max(ct_sims_prob2$time) & ct_sims_prob2$N_inf>10,]$sim
take_off3 <- ct_sims_prob3[ct_sims_prob3$time==max(ct_sims_prob3$time) & ct_sims_prob3$N_inf>10,]$sim

detection_times_to1 <- detection_times1[detection_times1$sim %in% take_off1,]
detection_times_to2 <- detection_times2[detection_times2$sim %in% take_off2,]
detection_times_to3 <- detection_times3[detection_times3$sim %in% take_off3,]

detection_times_to <- rbind(detection_times_to1,
                            detection_times_to2,
                            detection_times_to3)

#table(detection_times_to[detection_times_to$detected=="Yes",]$N_inf)

###################################################################################
# Plot figures

detection_times_to$N_cow2 <- factor(detection_times_to$N_cow)
levels(detection_times_to$N_cow2 ) <- c("100 cattle",
                                        "1000 cattle",
                                        "10,000 cattle")


panel1 <- ggplot(, aes(x=time ,color=factor(freq_test)))+
  stat_ecdf(linewidth=1)+
  facet_wrap(.~factor(N_cow2), scales = "free_y", nrow=3)+
  theme_bw(base_size=14)+
  xlab("Time outbreak is detected\n(days since first infection)")+
  ylab("Cumulative probability")+
  scale_x_continuous(limits=c(0,NA), expand=expansion(0.0), breaks=c(0,7,14,21,28,35,42,49))+
  scale_y_continuous(limits=c(0,1), expand=expansion(0.01))+
  scale_color_brewer("Testing interval (days)",palette = "Paired")+
  #geom_hline(yintercept = 0.5, linetype="dashed")+
  #geom_hline(yintercept = 0.95, linetype="dotted")+
  coord_cartesian(xlim=c(0,50))+
  theme(legend.position = c(0.75,0.83),
        legend.background = element_rect(color="black"),
        strip.background = element_rect(color="black", fill="white"))+
  guides(color=guide_legend(ncol=3))


panel2 <- ggplot(detection_times_to, aes(x=N_inf ,color=factor(freq_test)))+
  stat_ecdf(linewidth=1)+
  facet_wrap(.~factor(N_cow2), scales = "free_y", nrow=3)+
  theme_bw(base_size=14)+
  xlab("Outbreak size at detection\n(number of cattle infected)")+
  ylab("Cumulative probability")+
  scale_x_continuous(limits=c(0,NA), expand=expansion(0.0),
                     breaks=c(0,10,20,30,40,50,60,70,80,90,100))+#, breaks=c(0,1,5,10,20,50)
  scale_y_continuous(limits=c(0,1), expand=expansion(0.01))+
  scale_color_brewer("Testing interval (days)",palette = "Paired")+
  #geom_hline(yintercept = 0.5, linetype="dashed")+
  #geom_hline(yintercept = 0.95, linetype="dotted")+
  coord_cartesian(xlim=c(0,100))+
  geom_vline(xintercept = 1, linetype="dashed")+
  geom_vline(xintercept = 5, linetype="dashed")+
  geom_vline(xintercept = 10, linetype="dashed")+
  geom_vline(xintercept = 20, linetype="dashed")+
  geom_vline(xintercept = 50, linetype="dashed")+
  theme(legend.position = "none",
        strip.background = element_rect(color="black", fill="white"))


panel2 <- panel2+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank())

panel1 <- panel1+
  theme(panel.grid.minor = element_blank())

panel1 + panel2

ggsave('figures/testing_main.png', width=10, height=8)

################################################################################
# Stats
dt7 <- detection_times_to[detection_times_to$freq_test==28,]
quantile(dt7[dt7$N_cow==100,]$N_inf, 0.5)
quantile(dt7[dt7$N_cow==1000,]$N_inf, 0.5)
quantile(dt7[dt7$N_cow==10000,]$N_inf, 0.5)

quantile(dt7[dt7$N_cow==100,]$time, 0.5)
quantile(dt7[dt7$N_cow==1000,]$time, 0.5)
quantile(dt7[dt7$N_cow==10000,]$time, 0.5)


dt1 <- detection_times_to[detection_times_to$freq_test==7,]
nrow(dt1[dt1$N_cow==100 &dt1$N_inf<=1,])/nrow(dt1[dt1$N_cow==100,])
nrow(dt1[dt1$N_cow==1000 &dt1$N_inf<=1,])/nrow(dt1[dt1$N_cow==1000,])
nrow(dt1[dt1$N_cow==10000 &dt1$N_inf<=1,])/nrow(dt1[dt1$N_cow==10000,])

###########################################################################
# Explanation figure

inf_sims_ex <- inf_sims3[inf_sims3$sim %in% c(21,22),]
df_inc <- data.frame()
for(i in 0:200){
  for(j in unique(inf_sims_ex$sim)){
    temp <- inf_sims_ex[inf_sims_ex$sim==j & inf_sims_ex$time>=i &inf_sims_ex$time<(i+1),]
    df_inc <- rbind(df_inc,
                    data.frame(time=i,sim=j,inc=nrow(temp)),
                    data.frame(time=i+0.999,sim=j,inc=nrow(temp)) )
    
  }
  
}

cols <- RColorBrewer::brewer.pal(6, "Dark2")

methods1 <- ggplot(df_inc, aes(x=time, y=inc ,color=factor(sim)))+
  geom_line()+
  ylab("Daily infection\nincidence")+
  xlab("Time (days)")+
  scale_color_manual(values=cols[4:5])+
  coord_cartesian(xlim=c(0,100))+
  theme_bw(base_size = 14)+
  theme(legend.position = "none")



df_ct <- ct_sims_prob3[ct_sims_prob3$sim%in%c(21,22),]
df_ct1 <- df_ct
df_ct2 <- df_ct
df_ct1$time <- df_ct1$time-0.49999
df_ct2$time <- df_ct2$time+0.49999
df_ct <- rbind(df_ct1)

methods2 <- ggplot(df_ct, aes(x=time, y=Ct, color=factor(sim), group=interaction(sim,sim_Ct)))+
  geom_line()+
  #geom_line(aes(y=prob*35+10), linetype="dashed")+
  ylab("Ct value")+
  xlab("Time (days)")+
  geom_hline(yintercept = 36, linetype="dashed")+
  coord_cartesian(ylim=c(60,20),xlim=c(0,100), expand = expansion(0.))+
  #scale_y_continuous(sec.axis = sec_axis(~(.-10)/35, "test"))+
  scale_color_manual(values=cols[4:5])+
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

methods3 <- ggplot(df_ct, aes(x=time, y=prob, color=factor(sim), group=interaction(sim,sim_Ct)))+
  geom_line()+
  #geom_line(aes(y=prob*35+10), linetype="dashed")+
  ylab("Probability of\npositive test")+
  xlab("Time (days)")+
  geom_hline(yintercept = 36, linetype="dashed")+
  coord_cartesian(ylim=c(-0.01,1.01),xlim=c(0,100), expand = expansion(0.))+
  #scale_y_continuous(sec.axis = sec_axis(~(.-10)/35, "test"))+
  scale_color_manual(values=cols[4:5])+
  theme_bw(base_size = 14)+
  theme(legend.position = "none")


dt_ex <- detection_times3[detection_times3$sim %in% c(21,22),]

dt_ex$lab <- factor(dt_ex$freq_test)
levels(dt_ex$lab) <- c("Testing every 1 day",
                       "Testing every 2 days",
                       "Testing every 3 days",
                       "Testing every 4 days",
                       "Testing every 5 days",
                       "Testing every 6 days",
                       "Testing every 7 days",
                       "Testing every 10 days",
                       "Testing every 14 days",
                       "Testing every 21 days",
                       "Testing every 28 days")

methods4 <- ggplot(dt_ex[dt_ex$freq_test%in% c(1,3,5,7,14,28),], aes(x=time, fill=interaction(sim) ))+
  geom_histogram(stat="bin",position = position_dodge())+
  facet_wrap(.~lab, nrow=3, scale="free_y")+
  scale_fill_manual(values=cols[4:5])+
  theme_bw(base_size=14)+
  coord_cartesian(xlim=c(0,100), expand = expansion(1))+
  ylab("Probability outbreak detected")+
  xlab("Time (days)")+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_rect(fill="white"))


methods1 <- methods1+
  labs(tag="A")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  annotate("label", x=85,y=330, label="Simulated outbreaks", size=5, fill="white")
  #ggtitle("Simulated outbreaks")

methods2 <- methods2+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  annotate("label", x=70,y=50, label="Simulated Ct value in milk vat", size=5, fill="white")
  #ggtitle("Simulated Ct value in milk vat")

methods3 <- methods3+
  annotate("label", x=70,y=0.75, label="Probability milk sample tests\npositive (if tested)", size=5, fill="white")
  #theme(plot.title = element_text(size=14, hjust=0.5))+
  #ggtitle("Probability milk sample tests positive")


methods4 <- methods4+
  labs(tag="B")+
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size=14, hjust=0.5))+
  geom_hline(yintercept = 0)+
  ggtitle("Timing of outbreak detection under different testing frequencies")

#methods1 + methods2 + methods3+methods4 +plot_layout(nrow=4, heights=c(0.5,0.5,0.5))


methods_123<-methods1 + methods2 + methods3+ plot_layout(nrow=3, heights=c(0.5,0.5,0.5))
plot_grid(methods_123, methods4)

ggsave('figures/testing_method.png', width=14, height=8)




