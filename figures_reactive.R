setwd("C:/Users/EALESO/OneDrive - The University of Melbourne/Projects/H5N1 cow transmission")

################################################################################
# Load R packages

library(ggplot2)
library(patchwork)
library(arrow)
library(cowplot)

################################################################################
# Read in the simulation outcomes

df1 <- read_parquet("simulations/sims_reactive_p1.parquet")
df2 <- read_parquet("simulations/sims_reactive_p3.parquet")
df3 <- read_parquet("simulations/sims_reactive_p5.parquet")


################################################################################
# Relabel parameter set names

df1$param_set <- as.character(df1$param_set)
df2$param_set <- as.character(df2$param_set)
df3$param_set <- as.character(df3$param_set)

df1$param_set <- "Cow\u2013cow\ntransmission only"
df2$param_set <- "Mixed mode\ntransmission"
df3$param_set <- "Cow\u2013milking unit\ntransmission only"

df <- rbind(df1, df2, df3)

################################################################################
# New label for farm sizes

df$N_cows2 <- factor(df$N_cows)

levels(df$N_cows2 ) <- c("100 cattle",
                         "1000 cattle",
                         "10,000 cattle")

################################################################################
# Relevel reaction 

df$reaction <- factor(df$reaction)
levels(df$reaction) <- c("1","5","10","20","50","No intervention")
#df$reaction <- factor(df$reaction, levels= rev(c("1","5","10","20","50","No intervention")))


################################################################################
# Formatting for plotting

df_total <- df[df$time==max(df$time) & df$ID=="Total",]

df_prop <- data.frame()

for(i in unique(df$N_cows)){
  print(i)
  for(j in unique(df$reaction)){
    print(j)
    for(k in unique(df$param_set)){
      
      df_temp <- df_total[df_total$N_cows==i & df_total$reaction==j & df_total$param_set==k,]
      
      
      row <- data.frame(prop_react = sum(df_temp$N_R+df_temp$N_E+df_temp$N_I>=as.numeric(j) ),
                        prop_unfinished = sum(df_temp$N_E+df_temp$N_I>0),
                        prop_5 = sum(df_temp$N_R+df_temp$N_E+df_temp$N_I < 0.05*i),
                        N_cows = i,
                        reaction = j,
                        param_set = k)
      if(is.na(row$prop_react)){
        row$prop_react <- 0
      }
      
      df_prop <- rbind(df_prop, row)
      
    }
  }
}


df_prop$reaction <- factor(df_prop$reaction, levels= c("1","5","10","20","50","No intervention"))


df_prop[df_prop$N_cows==10000,]

df_prop[df_prop$param_set == "Mixed mode\ntransmission",]
df_prop[df_prop$param_set == "Cow\u2013milking unit\ntransmission only",]
df_prop[df_prop$param_set == "Cow\u2013cow\ntransmission only",]


ggplot(df[df$ID=="Total" & df$N_cows==10000 &df$param_set=="Mixed mode\ntransmission",], aes(x=time, y=(N_E+N_I)/N_cows,color=N_cows2 ,group=factor(sim)) )+
  geom_line()+
  facet_wrap(.~reaction, nrow=3)+
  coord_cartesian(ylim=c(0,0.1))+
  theme(legend.position = "none")
#levels(df_prop$reaction) <- as.numeric(c(1,2,3,4,5,6))

df_prop$N_cows2 <- factor(df_prop$N_cows)

levels(df_prop$N_cows2 ) <- c("100 cattle",
                         "1000 cattle",
                         "10,000 cattle")

################################################################################
cols <- RColorBrewer::brewer.pal(8,"Dark2")

################################################################################

df_prop1 <- df_prop[df_prop$param_set == "Cow\u2013milking unit\ntransmission only",]

panel1 <- ggplot(df_total[df_total$param_set=="Cow\u2013milking unit\ntransmission only",], aes(x=reaction, y=N_R/N_cows, color=N_cows2))+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0.75),shape=16, alpha=0.4)+
  #geom_crossbar(data=df1_summary, aes(x=factor(n_cohorts),fill=factor(param_set), ymax=median, ymin=median, y=median),color="black", position=position_dodge(width = 0.75), width=0.6, size=0.3)+
  ylab("Total infections (proportion of cattle)")+
  xlab("Outbreak size when reactive cohorting introduced\n(cumulative number of infected cattle)")+
  scale_y_continuous(limits=c(0,1), expand=expansion(0), breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
  scale_color_manual("Herd size", values = cols[1:3])+
  scale_fill_manual("Herd size", values = cols[1:3])+
  theme_bw(base_size = 14)+
  theme(legend.position=c(0.15,0.85),#legend.position=c(0.105,0.3),
        legend.background = element_rect(color="black"),
        panel.grid.major.x = element_blank())+
  geom_hline(yintercept = 0.05, linetype="dashed")+
  guides(fill="none",
         color=guide_legend(nrow=3))

panel2 <- ggplot(df_prop1, aes(y=1-prop_5/200, x=reaction, color=factor(N_cows), group=factor(N_cows) ))+
  geom_line(position=position_jitterdodge(jitter.width=0.0, dodge.width = 0.))+
  geom_point(position=position_jitterdodge(jitter.width=0.0, dodge.width = 0.0))+
  ylab("Proportion of outbreaks\n with >5% attack rate")+
  xlab("Outbreak size when reactive cohorting introduced\n(cumulative number of infected cattle)")+
  scale_y_continuous(limits=c(0,1), expand=expansion(0), breaks=c(0,0.5,1.0))+
  #scale_x_discrete(breaks=c(1,2,3,4,5,6), labels=rev(c("1","5","10","20","50","No intervention")))+
  scale_color_manual("Herd size", values = cols[1:3])+
  scale_fill_manual("Herd size", values = cols[1:3])+
  theme_bw(base_size = 14)+
  theme(legend.position="none",#legend.position=c(0.105,0.3),
        legend.background = element_rect(color="black"),
        panel.grid = element_blank())+
  guides(fill="none")


panel1 <- panel1+theme(axis.title.x = element_blank(),
                       axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       legend.position = c(0.92,0.3))

panel1 + panel2 +plot_layout(nrow=2, heights=c(1,0.4))

ggsave('figures/reactive_cohorts_param5.png', width=11, height=8)


################################################################################################################
# Cow-cow transmission

df_prop2 <- df_prop[df_prop$param_set == "Cow\u2013cow\ntransmission only",]

panel1 <- ggplot(df_total[df_total$param_set=="Cow\u2013cow\ntransmission only",], aes(x=reaction, y=N_R/N_cows, color=N_cows2))+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0.75),shape=16, alpha=0.4)+
  #geom_crossbar(data=df1_summary, aes(x=factor(n_cohorts),fill=factor(param_set), ymax=median, ymin=median, y=median),color="black", position=position_dodge(width = 0.75), width=0.6, size=0.3)+
  ylab("Total infections (proportion of cattle)")+
  xlab("Outbreak size when reactive cohorting introduced\n(cumulative number of infected cattle)")+
  scale_y_continuous(limits=c(0,1), expand=expansion(0), breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
  scale_color_manual("Herd size", values = cols[1:3])+
  scale_fill_manual("Herd size", values = cols[1:3])+
  theme_bw(base_size = 14)+
  theme(legend.position=c(0.15,0.85),#legend.position=c(0.105,0.3),
        legend.background = element_rect(color="black"),
        panel.grid.major.x = element_blank())+
  geom_hline(yintercept = 0.05, linetype="dashed")+
  guides(fill="none",
         color=guide_legend(nrow=3))

panel2 <- ggplot(df_prop2, aes(y=1-prop_5/200, x=reaction, color=factor(N_cows), group=factor(N_cows) ))+
  geom_line(position=position_jitterdodge(jitter.width=0.0, dodge.width = 0.))+
  geom_point(position=position_jitterdodge(jitter.width=0.0, dodge.width = 0.0))+
  ylab("Proportion of outbreaks\n with >5% attack rate")+
  xlab("Outbreak size when reactive cohorting introduced\n(cumulative number of infected cattle)")+
  scale_y_continuous(limits=c(0,1), expand=expansion(0), breaks=c(0,0.5,1.0))+
  #scale_x_discrete(breaks=c(1,2,3,4,5,6), labels=rev(c("1","5","10","20","50","No intervention")))+
  scale_color_manual("Herd size", values = cols[1:3])+
  scale_fill_manual("Herd size", values = cols[1:3])+
  theme_bw(base_size = 14)+
  theme(legend.position="none",#legend.position=c(0.105,0.3),
        legend.background = element_rect(color="black"),
        panel.grid = element_blank())+
  guides(fill="none")


panel1 <- panel1+theme(axis.title.x = element_blank(),
                       axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       legend.position = c(0.92,0.3))

panel1 + panel2 +plot_layout(nrow=2, heights=c(1,0.4))

ggsave('figures/reactive_cohorts_param1.png', width=11, height=8)


################################################################################################################
# Cow-cow transmission

df_prop3 <- df_prop[df_prop$param_set == "Mixed mode\ntransmission",]

panel1 <- ggplot(df_total[df_total$param_set=="Mixed mode\ntransmission",], aes(x=reaction, y=N_R/N_cows, color=N_cows2))+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0.75),shape=16, alpha=0.4)+
  #geom_crossbar(data=df1_summary, aes(x=factor(n_cohorts),fill=factor(param_set), ymax=median, ymin=median, y=median),color="black", position=position_dodge(width = 0.75), width=0.6, size=0.3)+
  ylab("Total infections (proportion of cattle)")+
  xlab("Outbreak size when reactive cohorting introduced\n(cumulative number of infected cattle)")+
  scale_y_continuous(limits=c(0,1), expand=expansion(0), breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
  scale_color_manual("Herd size", values = cols[1:3])+
  scale_fill_manual("Herd size", values = cols[1:3])+
  theme_bw(base_size = 14)+
  theme(legend.position=c(0.15,0.85),#legend.position=c(0.105,0.3),
        legend.background = element_rect(color="black"),
        panel.grid.major.x = element_blank())+
  geom_hline(yintercept = 0.05, linetype="dashed")+
  guides(fill="none",
         color=guide_legend(nrow=3))

panel2 <- ggplot(df_prop3, aes(y=1-prop_5/200, x=reaction, color=factor(N_cows), group=factor(N_cows) ))+
  geom_line(position=position_jitterdodge(jitter.width=0.0, dodge.width = 0.))+
  geom_point(position=position_jitterdodge(jitter.width=0.0, dodge.width = 0.0))+
  ylab("Proportion of outbreaks\n with >5% attack rate")+
  xlab("Outbreak size when reactive cohorting introduced\n(cumulative number of infected cattle)")+
  scale_y_continuous(limits=c(0,1), expand=expansion(0), breaks=c(0,0.5,1.0))+
  #scale_x_discrete(breaks=c(1,2,3,4,5,6), labels=rev(c("1","5","10","20","50","No intervention")))+
  scale_color_manual("Herd size", values = cols[1:3])+
  scale_fill_manual("Herd size", values = cols[1:3])+
  theme_bw(base_size = 14)+
  theme(legend.position="none",#legend.position=c(0.105,0.3),
        legend.background = element_rect(color="black"),
        panel.grid = element_blank())+
  guides(fill="none")


panel1 <- panel1+theme(axis.title.x = element_blank(),
                       axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       legend.position = c(0.92,0.3))

panel1 + panel2 +plot_layout(nrow=2, heights=c(1,0.4))

ggsave('figures/reactive_cohorts_param3.png', width=11, height=8)
