setwd("C:/Users/EALESO/OneDrive - The University of Melbourne/Projects/H5N1 cow transmission")

################################################################################
# Load R packages

library(ggplot2)
library(patchwork)
library(arrow)
library(cowplot)

################################################################################
# Read in the simulation outcomes

df1 <- read_parquet('simulations/sims_nCohortsFINAL.parquet')
df2 <- read_parquet('simulations/sims_seedCohortsFINAL.parquet')
df3 <- read_parquet('simulations/sims_cleansFINAL.parquet')
df4 <- read_parquet('simulations/sims_bctFINAL.parquet')
df5 <- read_parquet('simulations/sims_dtFINAL.parquet')
df6 <- read_parquet('simulations/sims_dt_cleaningFINAL.parquet')

################################################################################
# Relabel parameter set names

df1$param_set <- factor(df1$param_set)
levels(df1$param_set) <- c("Cow\u2013cow\ntransmission only",
                           "Mixed mode\ntransmission",
                           "Cow\u2013milking unit\ntransmission only")

df2$param_set <- factor(df2$param_set)
levels(df2$param_set) <- c("Cow\u2013cow\ntransmission only",
                           "Mixed mode\ntransmission",
                           "Cow\u2013milking unit\ntransmission only")

df3$param_set <- factor(df3$param_set)
levels(df3$param_set) <- c("Mixed mode\ntransmission",
                           "Cow\u2013milking unit\ntransmission only")


df4$param_set <- factor(df4$param_set)
levels(df4$param_set) <- c("Cow\u2013cow\ntransmission only")

df5$param_set <- factor(df5$param_set)
levels(df5$param_set) <- c("Cow\u2013milking unit\ntransmission only")

df6$param_set <- factor(df6$param_set)
levels(df5$param_set) <- c("Cow\u2013milking unit\ntransmission only")

#######################################################################################################################
# Code to reformat and plot df1
# Analyses to investigate the effect of milking cohorts
#######################################################################################################################

################################################################################
# Reformat df1 to 
# - Get median summary lines
# - Get example trajectories
# - Get final size for all outbreaks by cohorts

df1_total <- df1[df1$time==max(df1$time) & df1$ID=="Total",]

df1_summary <- data.frame()
for(i in unique(df1_total$n_cohorts)){
  for(j in unique(df1_total$param_set)){
    
    df1_temp <- df1_total[df1_total$n_cohorts==i & df1_total$param_set==j, ]
    
    row_df <- data.frame(n_cohorts = i,
                         param_set = j,
                         median = median(df1_temp$N_R),
                         mean = mean(df1_temp$N_R),
                         lwr = quantile(df1_temp$N_R, 0.05),
                         upr = quantile(df1_temp$N_R, 0.95))
    df1_summary <- rbind(df1_summary, row_df)
  }
}

df1_summary$param_set <- factor(df1_summary$param_set, levels=c("Cow\u2013cow\ntransmission only",
                                                                "Mixed mode\ntransmission",
                                                                "Cow\u2013milking unit\ntransmission only"))

# Example trajectories
df1_example <- df1[df1$sim==1 & df1$n_cohorts==5,]
#12

# Get final size for all paraem sets
df1_final <- df1[df1$time==999999999 & df1$n_cohorts==5,]
df1_final$new_sim = -99
for(i in unique(df1_final$param_set)){
  
  total_order <- order(df1_final[df1_final$param_set==i & df1_final$ID=="Total",]$N_R)
  for(j in unique(as.numeric(df1_final$sim))){
    df1_final[df1_final$param_set==i & df1_final$sim==j,]$new_sim <- which(total_order==j) 
  }
}


################################################################################
# Create figure 1
plt1a <- ggplot(df1_total, aes(x=factor(n_cohorts), y=N_R, color=param_set ))+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0.75),shape=16, alpha=0.4)+
  geom_crossbar(data=df1_summary, aes(x=factor(n_cohorts),fill=factor(param_set), ymax=median, ymin=median, y=median),color="black", position=position_dodge(width = 0.75), width=0.6, size=0.3)+
  ylab("Total infections")+
  xlab("Number of ordered milking cohorts")+
  scale_y_continuous(limits=c(0,1000), expand=expansion(0))+
  #geom_point(data=df1_total[df1_total$sim==12 &df1_total$n_cohorts==5
  #                          & df1_total$param_set!="Mixed mode\ntransmission",], #!!!!!!!!!!!!!!!!!!!!!
  #           shape=4, stroke=2,size=5,
  #           position=position_jitterdodge(dodge.width = 1))+
  coord_cartesian(ylim=c(0,1000))+
  scale_color_brewer("Transmission\nregime",palette="Dark2")+
  scale_fill_brewer("Transmission\nregime",palette="Dark2")+
  theme_bw(base_size = 14)+
  theme(legend.position=c(0.2,0.15),#legend.position=c(0.105,0.3),
        legend.background = element_rect(color="black"),
        panel.grid.major.x = element_blank())+
  guides(fill="none")


plt1b<-ggplot(df1_example[df1_example$ID!="Total" &df1_example$param_set!="Mixed mode\ntransmission",], aes(x=time, y=N_I+N_E, color=ID))+
  facet_wrap(.~param_set, nrow=1)+
  ylab("Infection prevalence")+
  xlab("Time (days)")+
  scale_color_brewer("Milking\ncohort",palette = "Set1")+
  geom_line()+
  scale_x_continuous(limits=c(0,110), expand = expansion(0))+
  theme_bw(base_size=14)+
  theme(strip.background = element_rect(color="black", fill="white"),
        legend.position = c(0.875,0.65),
        legend.background = element_rect(color="black"))



plt1c<-ggplot(df1_final[df1_final$ID!="Total" & df1_final$param_set!="Mixed mode\ntransmission",], aes(x=new_sim,y=N_R, fill=ID, group=new_sim))+
  geom_col(position = "fill")+
  geom_point(data=df1_final[df1_final$ID=="Total" & df1_final$param_set!="Mixed mode\ntransmission",], aes(x=new_sim, y=N_R/1000), size=0.1)+
  scale_fill_brewer(palette = "Set1")+
  xlab("Simulations (ordered by final epidemic size)")+
  ylab("Proportion of infections")+
  scale_x_continuous(expand=expansion(0))+
  scale_y_continuous(expand=expansion(0),
                     sec.axis = sec_axis(~., name="Total infections",breaks=c(0,0.250,0.5,0.75,1.),labels=c(0,250,500,750,1000)))+
  facet_wrap(.~param_set)+
  theme_bw(base_size=14)+
  theme(strip.background = element_rect(color="black", fill="white"),
        legend.position = "none",
        axis.text.x = element_blank())

plt1c

plt1a <- plt1a+labs(tag="A")
plt1b <- plt1b+labs(tag="B")
plt1c <- plt1c+labs(tag="C")

plt1a <- plt1a + theme(plot.tag.position = c(0.01, 0.97))
plt1b <- plt1b + theme(plot.tag.position = c(0.01, 0.97),
                       panel.spacing = unit(0.1, "lines"))
plt1c <- plt1c + theme(plot.tag.position = c(0.01, 0.97),
                       axis.ticks.x = element_blank(),
                       panel.spacing = unit(0.1, "lines"))


#plt1bc <- plt1b+plt1c+plot_layout(nrow=1, widths=c(1,1)) 
#plt1a/plt1bc


#ggsave('figures/n_cohorts.png', width=11, height=8)

plt1bc <- plt1b+plt1c+plot_layout(nrow=2) 
plot_grid(plt1a, plt1bc, rel_widths = c(1,0.9))

ggsave('figures/n_cohortsFINAL.png', width=11, height=8)
#ggsave('figures/n_cohorts_presentation.png', width=11, height=8)

#######################################################################################################################
# Code to reformat and plot df2
# Analyses investigating effect of cohort of seeding
#######################################################################################################################

################################################################################
# Reformat df2 to 
# - Get median summary lines
# - Get example trajectories
# - Get final size for all outbreaks by cohorts

df2_total <- df2[df2$time==max(df2$time) & df2$ID=="Total",]

df2_summary <- data.frame()
for(i in unique(df2_total$seed_cohorts)){
  for(j in unique(df2_total$param_set)){
    
    df2_temp <- df2_total[df2_total$seed_cohorts==i & df2_total$param_set==j, ]
    
    row_df <- data.frame(seed_cohorts = i,
                         param_set = j,
                         median = median(df2_temp$N_R),
                         mean = mean(df2_temp$N_R),
                         lwr = quantile(df2_temp$N_R, 0.05),
                         upr = quantile(df2_temp$N_R, 0.95))
    df2_summary <- rbind(df2_summary, row_df)
  }
}

df2_summary$param_set <- factor(df2_summary$param_set, levels=c("Cow\u2013cow\ntransmission only",
                                                                "Mixed mode\ntransmission",
                                                                "Cow\u2013milking unit\ntransmission only"))

# Example trajectories
df2_example <- df2[df2$sim==17 & df2$param_set=="Cow\u2013milking unit\ntransmission only",]

# Get final size for all paraem sets
df2_final <- df2[df2$time==999999999 & df2$param_set=="Cow\u2013milking unit\ntransmission only",]
df2_final$new_sim = -99
for(i in unique(df2_final$seed_cohorts)){
  
  total_order <- order(df2_final[df2_final$seed_cohorts==i & df2_final$ID=="Total",]$N_R)
  for(j in unique(as.numeric(df2_final$sim))){
    df2_final[df2_final$seed_cohorts==i & df2_final$sim==j,]$new_sim <- which(total_order==j) 
  }
}


################################################################################
# Create figure 2
plt2a <- ggplot(df2_total, aes(x=factor(seed_cohorts), y=N_R, color=param_set))+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0.75),shape=16, alpha=0.4)+
  geom_crossbar(data=df2_summary, aes(x=factor(seed_cohorts),fill=factor(param_set), ymax=median, ymin=median, y=median),color="black", position=position_dodge(width = 0.75), width=0.6, size=0.3)+
  ylab("Total infections")+
  xlab("Cohort infections are introduced to")+
  scale_y_continuous(limits=c(0,1000), expand=expansion(0))+
  #geom_point(data=df2_total[df2_total$sim==17 &df2_total$param_set=="Cow\u2013milking unit\ntransmission only",], #!!!!!!!!!!!!!!!!!!!!!
  #           shape=4, stroke=2,size=2,
  #           position=position_jitterdodge(dodge.width = 1))+
  #coord_cartesian(ylim=c(0,1000))+
  scale_color_brewer("Transmission regime",palette="Dark2")+
  scale_fill_brewer("Transmission regime",palette="Dark2")+
  theme_bw(base_size = 14)+
  theme(legend.position=c(0.7,0.8),
        legend.background = element_rect(color="black"),
        panel.grid.major.x = element_blank())+
  guides(fill="none")


df2_example$seed_cohorts <- factor(df2_example$seed_cohorts)
levels(df2_example$seed_cohorts) <- c("Infections seeding\ncohort 1",
                                    "Infections seeding\ncohort 2",
                                    "Infections seeding\ncohort 3",
                                    "Infections seeding\ncohort 4",
                                    "Infections seeding\ncohort 5")


plt2b<-ggplot(df2_example[df2_example$ID!="Total" 
                          &df2_example$seed_cohorts%in%c("Infections seeding\ncohort 1",
                                                         "Infections seeding\ncohort 3",
                                                         "Infections seeding\ncohort 5") ,], aes(x=time, y=N_I+N_E, color=ID))+
  facet_wrap(.~seed_cohorts, nrow=1)+
  ylab("Infection prevalence")+
  xlab("Time (days)")+
  scale_color_brewer("Milking\ncohort",palette = "Set1")+
  geom_line()+
  scale_x_continuous(limits=c(0,150), expand = expansion(0), breaks=c(0,50,100))+
  theme_bw(base_size=14)+
  theme(strip.background = element_rect(color="black", fill="white"),
        legend.position = c(0.88,0.57),
        legend.background = element_rect(color="black"))

df2_final$seed_cohorts <- factor(df2_final$seed_cohorts)
levels(df2_final$seed_cohorts) <- c("Infections seeding\ncohort 1",
                                    "Infections seeding\ncohort 2",
                                    "Infections seeding\ncohort 3",
                                    "Infections seeding\ncohort 4",
                                    "Infections seeding\ncohort 5")

plt2c<-ggplot(df2_final[df2_final$ID!="Total"
                        & df2_final$seed_cohorts%in%c("Infections seeding\ncohort 1",
                                                      "Infections seeding\ncohort 3",
                                                      "Infections seeding\ncohort 5"),], aes(x=new_sim,y=N_R, fill=ID, group=new_sim))+
  geom_col(position = "fill")+
  geom_point(data=df2_final[df2_final$ID=="Total" & df2_final$param_set!="Mixed mode\ntransmission"
                            & df2_final$seed_cohorts%in%c("Infections seeding\ncohort 1",
                                                          "Infections seeding\ncohort 3",
                                                          "Infections seeding\ncohort 5"),], aes(x=new_sim, y=N_R/1000), size=0.1)+
  scale_fill_brewer(palette = "Set1")+
  xlab("Simulations (ordered by final epidemic size)")+
  ylab("Proportion of infections")+
  scale_x_continuous(expand=expansion(0.01))+
  scale_y_continuous(expand=expansion(0),
                     sec.axis = sec_axis(~., name="Total infections",breaks=c(0,0.250,0.5,0.75,1.),labels=c(0,250,500,750,1000)))+
  facet_wrap(.~seed_cohorts, nrow=1)+
  theme_bw(base_size=14)+
  theme(strip.background = element_rect(color="black", fill="white"),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        title = element_text(hjust=1))

plt2c
plt2a <- plt2a+labs(tag="A")
plt2b <- plt2b+labs(tag="B")
plt2c <- plt2c+labs(tag="C")

plt2a <- plt2a + theme(plot.tag.position = c(0.01, 0.97))
plt2b <- plt2b + theme(plot.tag.position = c(0.01, 0.97),
                       panel.spacing = unit(0.1, "lines"))
plt2c <- plt2c + theme(plot.tag.position = c(0.01, 0.97),
                       panel.spacing = unit(0.1, "lines"))

plt2bc <- plt2b+ plt2c +plot_layout(nrow=2, heights=c(1,1))
plot_grid(plt2a,plt2bc, rel_widths = c(1,1.2))
ggsave('figures/seed_cohortsFINAL.png', width=11, height=7)

#######################################################################################################################
# Code to reformat and plot df3 
# Sensitivity analyses assuming enhanced cleaning once per milking period
#######################################################################################################################

################################################################################
# Reformat df3 to 
# - Get median summary lines
# - Get example trajectories
# - Get final size for all outbreaks by cohorts

df3_total <- df3[df3$time==max(df3$time) & df3$ID=="Total",]

df3_summary <- data.frame()

for(i in unique(df3_total$cleaning_before)){
  for(j in unique(df3_total$param_set)){
    for(k in unique(df3_total$cleaning_effect)){
      
      df3_temp <- df3_total[df3_total$cleaning_before==i
                            & df3_total$param_set==j
                            & df3_total$cleaning_effect==k, ]
      
      row_df <- data.frame(cleaning_before = i,
                           param_set = j,
                           cleaning_effect=k,
                           median = median(df3_temp$N_R),
                           mean = mean(df3_temp$N_R),
                           lwr = quantile(df3_temp$N_R, 0.05),
                           upr = quantile(df3_temp$N_R, 0.95))
      df3_summary <- rbind(df3_summary, row_df)
      
    }
  }
}

df3_summary$param_set <- factor(df3_summary$param_set, levels=c("Mixed mode\ntransmission",
                                                                "Cow\u2013milking unit\ntransmission only"))


df3_final <- df3[df3$time==999999999 & df3$param_set%in%c("Cow\u2013milking unit\ntransmission only","Mixed mode\ntransmission"),]
df3_final$new_sim = -99
for(i in unique(df3_final$cleaning_before)){
  
  for(k in unique(df3_final$cleaning_effect)){
    for(l in unique(df3_final$param_set)){
      total_order <- order(df3_final[df3_final$cleaning_before==i 
                                     & df3_final$ID=="Total"
                                     & df3_final$cleaning_effect==k
                                     & df3_final$param_set==l,]$N_R)
      for(j in unique(as.numeric(df3_final$sim))){
        df3_final[df3_final$cleaning_effect==k 
                  &df3_final$cleaning_before==i 
                  & df3_final$sim==j
                  & df3_final$param_set==l,]$new_sim <- which(total_order==j) 
      }
    }
    
    
  }
  
  
}


################################################################################
# Some cleaning of variables and merging with the baseline from df1

df3_total$cleaning_before <- factor(df3_total$cleaning_before)
levels(df3_total$cleaning_before) <- c("5\u20131","1\u20132","2\u20133","3\u20134","4\u20135")

df3_noclean <- df1_total[df1_total$param_set%in%c("Cow\u2013milking unit\ntransmission only", "Mixed mode\ntransmission") &df1_total$n_cohorts==5,]
df3_noclean$cleaning_before <- "No cleaning"
df3_noclean$cleaning_effect <- "No effect"

df3_total$n_cohorts=5
df3_total <- rbind(df3_total, df3_noclean)


df3_summary$cleaning_before <- factor(df3_summary$cleaning_before)
levels(df3_summary$cleaning_before) <- c("5\u20131","1\u20132","2\u20133","3\u20134","4\u20135")


df3_noclean_summary <- df1_summary[df1_summary$param_set%in%c("Cow\u2013milking unit\ntransmission only", "Mixed mode\ntransmission") &df1_summary$n_cohorts==5,]
df3_noclean_summary$cleaning_before <- "No cleaning"
df3_noclean_summary$cleaning_effect <- "No effect"

df3_summary$n_cohorts = 5
df3_summary <- rbind(df3_summary, df3_noclean_summary)

df3_total$cleaning_before <- factor(df3_total$cleaning_before,
                                    levels=c("No cleaning","5\u20131","1\u20132","2\u20133","3\u20134","4\u20135"))



df3_summary$cleaning_effect <- factor(df3_summary$cleaning_effect)
levels(df3_summary$cleaning_effect) <- c("50% reduction in infectiousness",
                                         "90% reduction in infectiousness",
                                         "Not applicable")

df3_total$cleaning_effect <- factor(df3_total$cleaning_effect)
levels(df3_total$cleaning_effect) <- c("50% reduction in infectiousness",
                                         "90% reduction in infectiousness",
                                         "Not applicable")

################################################################################
# Create cleaning sensitivity figure 1

ggplot(df3_total, aes(x=factor(cleaning_before), y=N_R, color=factor(cleaning_effect) ))+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0.75),shape=16, alpha=0.4)+
  geom_crossbar(data=df3_summary, aes(x=factor(cleaning_before),
                                      fill=factor(cleaning_effect),
                                      ymax=median, ymin=median, y=median),
                color="black",
                position=position_dodge(width = 0.75), width=0.6, size=0.3)+
  ylab("Total infections")+
  xlab("Cleaning between milking cohorts")+
  scale_y_continuous(limits=c(0,NA), expand=expansion(0.02))+
  facet_wrap(.~param_set, nrow=2, scales="free_y")+
  scale_color_brewer("Effect of cleaning\nmilking units",palette="Dark2")+
  scale_fill_brewer("Effect of cleaning\nmilking units",palette="Dark2")+
  theme_bw(base_size = 14)+
  theme(strip.background = element_rect(color="black", fill="white"),
        legend.position="right",
        legend.background = element_rect(color="black"),
        panel.grid.major.x = element_blank())+
  guides(fill="none")

ggsave('figures/sens_cleaning1FINAL.png', width=9, height=9)


################################################################################
# Creating cleaning sensitivity figure 2

df3_final$cleaning_before <- factor(df3_final$cleaning_before)
levels(df3_final$cleaning_before) <- c("Cleaning between\ncohorts 5\u20131",
                                       "Cleaning between\ncohorts 1\u20132",
                                       "Cleaning between\ncohorts 2\u20133",
                                       "Cleaning between\ncohorts 3\u20134",
                                       "Cleaning between\ncohorts 4\u20135")




df3_final$cleaning_effect <- factor(df3_final$cleaning_effect)
levels(df3_final$cleaning_effect) <- c("Cleaning causes 50%\nreduction in infectiousness",
                                         "Cleaning causes 90%\nreduction in infectiousness",
                                         "Not applicable")

panelA <- ggplot(df3_final[df3_final$ID!="Total" & df3_final$param_set=="Mixed mode\ntransmission",], aes(x=new_sim,y=N_R, fill=ID, group=new_sim))+
  geom_col(position = "fill")+
  geom_point(data=df3_final[df3_final$ID=="Total" & df3_final$param_set=="Mixed mode\ntransmission",], aes(x=new_sim, y=N_R/1000), size=0.1)+
  scale_fill_brewer(palette = "Set1")+
  xlab("Simulations (ordered by final epidemic size)")+
  ylab("Proportion of infections")+
  scale_x_continuous(expand=expansion(0))+
  scale_y_continuous(expand=expansion(0),
                     sec.axis = sec_axis(~., name="Total infections",breaks=c(0,0.250,0.5,0.75,1.),labels=c(0,250,500,750,1000)))+
  facet_grid(rows=vars(cleaning_effect), cols=vars(cleaning_before) )+
  theme_bw(base_size=14)+
  theme(strip.background = element_rect(color="black", fill="white"),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

panelB <- ggplot(df3_final[df3_final$ID!="Total" & df3_final$param_set=="Cow\u2013milking unit\ntransmission only",], aes(x=new_sim,y=N_R, fill=ID, group=new_sim))+
  geom_col(position = "fill")+
  geom_point(data=df3_final[df3_final$ID=="Total" & df3_final$param_set=="Cow\u2013milking unit\ntransmission only",], aes(x=new_sim, y=N_R/1000),fill="black", size=0.1)+
  scale_fill_brewer("Milking cohort",palette = "Set1")+
  xlab("Simulations (ordered by final epidemic size)")+
  ylab("Proportion of infections")+
  scale_x_continuous(expand=expansion(0))+
  scale_y_continuous(expand=expansion(0),
                     sec.axis = sec_axis(~., name="Total infections",breaks=c(0,0.250,0.5,0.75,1.),labels=c(0,250,500,750,1000)))+
  facet_grid(rows=vars(cleaning_effect), cols=vars(cleaning_before) )+
  theme_bw(base_size=14)+
  theme(strip.background = element_rect(color="black", fill="white"),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

panelA <- panelA +ggtitle("Mixed mode transmission")+
  theme(plot.title=element_text(hjust=0.5))
panelB <- panelB +ggtitle("Cow\u2013milking unit transmission only")+
  theme(plot.title=element_text(hjust=0.5))
panelA +panelB + plot_layout(nrow=2)

ggsave('figures/sens_cleaning2FINAL.png', width=10, height=12)



#######################################################################################################################
# Code to reformat and plot df4 
# Sensitivity analysis assuming between cohort cow-cow transmission
#######################################################################################################################

################################################################################
# Reformat df4 to 
# - Get median summary lines
# - Get example trajectories
# - Get final size for all outbreaks by cohorts

df4_total <- df4[df4$time==max(df4$time) & df4$ID=="Total",]


#!!!!!!!!!!!!!!!!!
df4_summary <- data.frame()
for(i in unique(df4_total$between_cohort)){
  for(j in unique(df4_total$param_set)){
    
    df4_temp <- df4_total[df4_total$between_cohort==i & df4_total$param_set==j, ]
    
    row_df <- data.frame(between_cohort = i,
                         param_set = j,
                         median = median(df4_temp$N_R),
                         mean = mean(df4_temp$N_R),
                         lwr = quantile(df4_temp$N_R, 0.05),
                         upr = quantile(df4_temp$N_R, 0.95))
    df4_summary <- rbind(df4_summary, row_df)
  }
}


# Get final size for all paraem sets
df4_final <- df4[df4$time==999999999,]
df4_final$new_sim = -99
for(i in unique(df4_final$between_cohort)){
  
  total_order <- order(df4_final[df4_final$between_cohort==i & df4_final$ID=="Total",]$N_R)
  for(j in unique(as.numeric(df4_final$sim))){
    df4_final[df4_final$between_cohort==i & df4_final$sim==j,]$new_sim <- which(total_order==j) 
  }
}




################################################################################
# Plotting figure 

plt4a <- ggplot(df4_total, aes(x=factor(between_cohort), y=N_R))+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0.75), color="grey32",shape=16, alpha=0.4)+
  geom_crossbar(data=df4_summary, aes(x=factor(between_cohort), ymax=median, ymin=median, y=median), position=position_dodge(width = 0.75), width=0.6, size=0.3)+
  ylab("Total infections")+
  xlab("Proportion of cow\u2013cow\ntransmission rate between cohorts")+
  scale_y_continuous(limits=c(0,1000), expand=expansion(0))+
  coord_cartesian(ylim=c(0,1000))+
  theme_bw(base_size = 14)+
  theme(legend.position=c(0.2,0.9),#legend.position=c(0.105,0.3),
        legend.background = element_rect(color="black"),
        panel.grid.major.x = element_blank())+
  guides(color = guide_legend(override.aes = list(shape = 19)))



plt4b <- ggplot(df4_final[df4_final$ID!="Total",], aes(x=new_sim,y=N_R, fill=ID, group=new_sim))+
  geom_col(position = "fill")+
  geom_point(data=df4_final[df4_final$ID=="Total",], aes(x=new_sim, y=N_R/1000),fill="black", size=0.1)+
  scale_fill_brewer("Milking cohort",palette = "Set1")+
  xlab("Simulations (ordered by final epidemic size)")+
  ylab("Proportion of infections")+
  scale_x_continuous(expand=expansion(0))+
  scale_y_continuous(expand=expansion(0),
                     sec.axis = sec_axis(~., name="Total infections",breaks=c(0,0.250,0.5,0.75,1.),labels=c(0,250,500,750,1000)))+
  facet_wrap(.~between_cohort, nrow=5)+
  theme_bw(base_size=14)+
  theme(strip.background = element_rect(color="black", fill="white"),
        legend.position = "bottom",
        axis.text.x = element_blank())



plt4a <- plt4a+labs(tag="A")
plt4b <- plt4b+labs(tag="B")

plt4a <- plt4a + theme(plot.tag.position = c(0.01, 0.99))
plt4b <- plt4b + ggtitle("Proportion of cow\u2013cow\ntransmission rate between cohorts")+
  theme(plot.tag.position = c(0.01, 0.99),
                       panel.spacing = unit(0.1, "lines"),
                       plot.title=element_text(hjust=0.5, size=14))


plot_grid(plt4a,plt4b)

ggsave('figures/sens_bctFINAL.png', width=10, height=12)


#######################################################################################################################
# Code to reformat and plot df5
# Sensitivity assuming longer half life of milking unit infectiousness
#######################################################################################################################

################################################################################
# Reformat df5 to 
# - Get median summary lines
# - Get example trajectories
# - Get final size for all outbreaks by cohorts

df5_total <- df5[df5$time==max(df5$time) & df5$ID=="Total",]


df5_summary <- data.frame()
for(i in unique(df5_total$decay_time)){
  for(j in unique(df5_total$seed_cohorts)){
    
    df5_temp <- df5_total[df5_total$decay_time==i & df5_total$seed_cohorts==j, ]
    
    row_df <- data.frame(decay_time = i,
                         seed_cohorts = j,
                         median = median(df5_temp$N_R),
                         mean = mean(df5_temp$N_R),
                         lwr = quantile(df5_temp$N_R, 0.05),
                         upr = quantile(df5_temp$N_R, 0.95))
    df5_summary <- rbind(df5_summary, row_df)
  }
}


# Get final size for all paraem sets
df5_final <- df5[df5$time==999999999,]
df5_final$new_sim = -99
for(i in unique(df5_final$decay_time)){
  
  for(k in unique(df5_final$seed_cohorts)){
    total_order <- order(df5_final[df5_final$decay_time==i
                                   & df5_final$seed_cohorts==k
                                   & df5_final$ID=="Total",]$N_R)
    for(j in unique(as.numeric(df5_final$sim))){
      df5_final[df5_final$decay_time==i & df5_final$seed_cohorts==k & df5_final$sim==j,]$new_sim <- which(total_order==j) 
    }
  }
  
}


################################################################################
# Minor formatting

df5_total$decay_time <- factor(df5_total$decay_time)
levels(df5_total$decay_time) <- c("1.74 hours\n(2 times original)",
                                  "3.48 hours\n(4 times original)",
                                  "2.61 hours\n(3 times original)",
                                  "0.87 hours\n(original)")

df5_total$decay_time <- factor(df5_total$decay_time, 
                               levels = c("0.87 hours\n(original)",
                                          "1.74 hours\n(2 times original)",
                                          "2.61 hours\n(3 times original)",
                                          "3.48 hours\n(4 times original)"))


df5_summary$decay_time <- factor(df5_summary$decay_time)
levels(df5_summary$decay_time) <- c("1.74 hours\n(2 times original)",
                                  "3.48 hours\n(4 times original)",
                                  "2.61 hours\n(3 times original)",
                                  "0.87 hours\n(original)")

df5_summary$decay_time <- factor(df5_summary$decay_time, 
                               levels = c("0.87 hours\n(original)",
                                          "1.74 hours\n(2 times original)",
                                          "2.61 hours\n(3 times original)",
                                          "3.48 hours\n(4 times original)"))


df5_final$decay_time <- factor(df5_final$decay_time)
levels(df5_final$decay_time) <- c("1.74 hours\n(2 times original)",
                                  "3.48 hours\n(4 times original)",
                                  "2.61 hours\n(3 times original)",
                                  "0.87 hours\n(original)")

df5_final$decay_time <- factor(df5_final$decay_time, 
                               levels = c("0.87 hours\n(original)",
                                          "1.74 hours\n(2 times original)",
                                          "2.61 hours\n(3 times original)",
                                          "3.48 hours\n(4 times original)"))


################################################################################
# Plotting figure 

ggplot(df5_total, aes(x=factor(seed_cohorts), y=N_R, color=decay_time))+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0.75),shape=16, alpha=0.4)+
  geom_crossbar(data=df5_summary, aes(x=factor(seed_cohorts), ymax=median, ymin=median, y=median, fill=decay_time),color="black", position=position_dodge(width = 0.75), width=0.6, size=0.3)+
  ylab("Total infections")+
  xlab("Cohort infections are introduced to")+
  scale_y_continuous(limits=c(0,1000), expand=expansion(0))+
  scale_color_brewer("Half life of\nmilking unit\ninfectiousness",palette="Dark2")+
  scale_fill_brewer("Half life of\nmilking unit\ninfectiousness",palette="Dark2")+
  coord_cartesian(ylim=c(0,1000))+
  theme_bw(base_size = 14)+
  theme(legend.position="right",#legend.position=c(0.105,0.3),
        legend.background = element_rect(color="black"),
        panel.grid.major.x = element_blank(),
        legend.key.size = unit(1, "cm"))+
  guides(fill="none")

ggsave('figures/sens_dt1FINAL.png', width=8, height=6)

df5_final$seed_cohorts <- factor(df5_final$seed_cohorts)
levels(df5_final$seed_cohorts) <- c("Infections seeded\nin cohort 1",
                                    "Infections seeded\nin cohort 2",
                                    "Infections seeded\nin cohort 3",
                                    "Infections seeded\nin cohort 4",
                                    "Infections seeded\nin cohort 5")

ggplot(df5_final[df5_final$ID!="Total",], aes(x=new_sim,y=N_R, fill=ID, group=new_sim))+
  geom_col(position = "fill")+
  geom_point(data=df5_final[df5_final$ID=="Total",], aes(x=new_sim, y=N_R/1000),fill="black", size=0.1)+
  scale_fill_brewer("Milking cohort",palette = "Set1")+
  xlab("Simulations (ordered by final epidemic size)")+
  ylab("Proportion of infections")+
  scale_x_continuous(expand=expansion(0))+
  scale_y_continuous(expand=expansion(0),
                     sec.axis = sec_axis(~., name="Total infections",breaks=c(0,0.250,0.5,0.75,1.),labels=c(0,250,500,750,1000)))+
  facet_grid(cols=vars(decay_time), rows=vars(seed_cohorts))+
  theme_bw(base_size=14)+
  ggtitle("Half life of milking unit infectiousness")+
  theme(strip.background = element_rect(color="black", fill="white"),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        plot.title=element_text(hjust=0.5, size=14))

ggsave('figures/sens_dt2FINAL.png', width=10, height=12)


#######################################################################################################
# Investigating effect of cleaning between milking when decay rate slower
# (Milking units infectious for longer)
#######################################################################################################
df6$cleaning_effect <- factor(df6$cleaning_effect)
levels(df6$cleaning_effect) <- c("No effect",
                                 "50% reduction\nin infectiousness",
                                 "90% reduction\nin infectiousness",
                                 "99% reduction\nin infectiousness",
                                 "99.9% reduction\nin infectiousness")

df6_total <- df6[df6$time==max(df6$time) & df6$ID=="Total",]


#!!!!!!!!!!!!!!!!!
df6_summary <- data.frame()
for(i in unique(df6_total$seed_cohorts)){
  for(j in unique(df6_total$cleaning_effect)){
    
    df6_temp <- df6_total[df6_total$seed_cohorts==i & df6_total$cleaning_effect==j, ]
    
    row_df <- data.frame(seed_cohorts = i,
                         cleaning_effect = j,
                         median = median(df6_temp$N_R),
                         mean = mean(df6_temp$N_R),
                         lwr = quantile(df6_temp$N_R, 0.05),
                         upr = quantile(df6_temp$N_R, 0.95))
    df6_summary <- rbind(df6_summary, row_df)
  }
}


# Get final size for all paraem sets
df6_final <- df6[df6$time==999999999,]
df6_final$new_sim = -99
for(i in unique(df6_final$cleaning_effect)){
  
  for(k in unique(df6_final$seed_cohorts)){
    total_order <- order(df6_final[df6_final$cleaning_effect==i
                                   & df6_final$seed_cohorts==k
                                   & df6_final$ID=="Total",]$N_R)
    for(j in unique(as.numeric(df6_final$sim))){
      df6_final[df6_final$cleaning_effect==i & df6_final$seed_cohorts==k & df6_final$sim==j,]$new_sim <- which(total_order==j) 
    }
  }
  
}


df6_summary$cleaning_effect <- factor(df6_summary$cleaning_effect, levels =c("No effect",
                                                                             "50% reduction\nin infectiousness",
                                                                             "90% reduction\nin infectiousness",
                                                                             "99% reduction\nin infectiousness",
                                                                             "99.9% reduction\nin infectiousness") )


################################################################################
#
################################################################################
# Create cleaning with long decay time

ggplot(df6_total, aes(x=factor(seed_cohorts), y=N_R, color=cleaning_effect ))+
  geom_point(position=position_jitterdodge(jitter.width=0.25, dodge.width = 0.5),shape=16, alpha=0.4)+
  geom_crossbar(data=df6_summary, aes(x=factor(seed_cohorts),fill=cleaning_effect, ymax=median, ymin=median, y=median),color="black", position=position_dodge(width = 0.5), width=0.6, size=0.3)+
  ylab("Total infections")+
  xlab("Cohort infections are introduced to")+
  scale_y_continuous(limits=c(0,1000), expand=expansion(0))+
  #geom_point(data=df1_total[df1_total$sim==12 &df1_total$n_cohorts==5
  #                          & df1_total$param_set!="Mixed mode\ntransmission",], #!!!!!!!!!!!!!!!!!!!!!
  #           shape=4, stroke=2,size=5,
  #           position=position_jitterdodge(dodge.width = 1))+
  coord_cartesian(ylim=c(0,1000))+
  scale_color_brewer("Effect of cleaning\nmilking units",palette="Dark2")+
  scale_fill_brewer("Effect of cleaning\nmilking units",palette="Dark2")+
  theme_bw(base_size = 14)+
  theme(legend.position=c(0.15,0.2),#legend.position=c(0.105,0.3),
        legend.background = element_rect(color="black"),
        panel.grid.major.x = element_blank())+
  guides(fill="none")


ggsave('figures/seed_cohorts_with_cleaning1FINAL.png', width=11, height=8)


df6_final$seed_cohorts <- factor(df6_final$seed_cohorts)
levels(df6_final$seed_cohorts) <- c("Infections seeded\nin cohort 1",
                                    "Infections seeded\nin cohort 2",
                                    "Infections seeded\nin cohort 3",
                                    "Infections seeded\nin cohort 4",
                                    "Infections seeded\nin cohort 5")

ggplot(df6_final[df6_final$ID!="Total"], aes(x=new_sim,y=N_R, fill=ID, group=new_sim))+
  geom_col(position = "fill")+
  geom_point(data=df6_final[df6_final$ID=="Total",], aes(x=new_sim, y=N_R/1000),fill="black", size=0.1)+
  scale_fill_brewer("Milking cohort",palette = "Set1")+
  xlab("Simulations (ordered by final epidemic size)")+
  ylab("Proportion of infections")+
  scale_x_continuous(expand=expansion(0))+
  scale_y_continuous(expand=expansion(0),
                     sec.axis = sec_axis(~., name="Total infections",breaks=c(0,0.250,0.5,0.75,1.),labels=c(0,250,500,750,1000)))+
  facet_grid(rows=vars(cleaning_effect), cols=vars(seed_cohorts) )+
  theme_bw(base_size=14)+
  theme(strip.background = element_rect(color="black", fill="white"),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing.y = unit(0.8, 'lines'))

ggsave('figures/seed_cohorts_with_cleaning2FINAL.png', width=9, height=9)

