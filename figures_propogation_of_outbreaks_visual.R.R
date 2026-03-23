
setwd("C:/Users/EALESO/OneDrive - The University of Melbourne/Projects/H5N1 cow transmission")
library(patchwork)
library(egg)

cols3 <- c("black", RColorBrewer::brewer.pal(8, 'Dark2')[4:5])

###############################################################################################
# Propagation of outbreaks explanatory figure
##############################################################################################

basic_SIR <- function(params, init, tmax=100){
  
  
  S = init$S
  I = init$I
  R = init$R
  N = S+I+R
  
  df <- data.frame(time=0, S=S, I=I, R=R)
  
  beta = params$beta
  gamma = params$gamma
  
  for(i in 1:tmax){
    if(I>0){
      S_I <- rpois(1, lambda = beta*S*I/N)
      I_R <- rpois(1, lambda = gamma*I)
      
      S <- S - S_I
      I <- I + S_I -I_R
      R <- R + I_R
      
      row_df <- data.frame(time=i, S=S, I=I, R=R)
      df <- rbind(df, row_df)
    }
    print(i)
    
    
  }
  
  df
  
}

reformat_sim <- function(sim1){
  sim1$inc <- 0
  sim1$inc[2:nrow(sim1)] <- sim1$S[1:(nrow(sim1)-1)]-sim1$S[2:nrow(sim1)]
  
  
  sim1a <- sim1
  sim1b <- sim1
  sim1b$time <- sim1b$time +0.9999999
  
  sim1 <- rbind(sim1a, sim1b)
  sim1
}

set.seed(12345)


init <- data.frame(S=1995,I=5, R=0)

set.seed(12345)
params <- data.frame(beta=0.3, gamma=0.16)
sim1 <- basic_SIR(params, init, tmax=100)

set.seed(12345)
params <- data.frame(beta=0.23, gamma=0.16)
sim2 <- basic_SIR(params, init, tmax=100)


set.seed(12345)
init <- sim1[30,]
params <- data.frame(beta=0.16, gamma=0.16)
sim3 <- basic_SIR(params, init, tmax=70)
sim3$time <- sim3$time+29
sim3 <-rbind(sim1[1:29,],sim3)


sim1 <- reformat_sim(sim1)
sim2 <- reformat_sim(sim2)
sim3 <- reformat_sim(sim3)


sim1$label <- "Unmitigated\noutbreak"
sim2$label <- "Pre-emptive\ninterventions"
sim3$label <- "Reactive\ninterventions"

sim <- rbind(sim1, sim2, sim3)

p1 <- ggplot(sim)+
  geom_line(aes(x=time, y=inc, color=label), size=1)+
  geom_line(data=sim1,aes(x=time, y=inc, color=label),size=1)+
  theme_bw(base_size = 14)+
  xlab("Time")+
  ylab(expression(i[j](t)))+
  scale_color_manual("Source farm", values=cols3, breaks = c("Unmitigated\noutbreak","Pre-emptive\ninterventions","Reactive\ninterventions"))+
  geom_vline(xintercept = 30, linetype="dashed")+
  ggtitle("Infection incidence on source farm")+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(angle=0, vjust=0.5, size=16),
        legend.position = c(0.15,0.72),
        legend.background = element_rect(color="black"),
        plot.title = element_text(hjust = 0.5, color = "black", size=13))

p1
r1 <- data.frame(time = seq(0,100),
                 rt   = rep(1, 101),
                 label2 = "baseline")
r2 <- data.frame(time = c(seq(0,30), seq(30,100)) ,
                 rt   = c(rep(1, 31),rep(0.05,71) ),
                 label2 = "Reactive export controls")
r3 <- data.frame(time = seq(0,100),
                 rt   = rep(0.3, 101),
                 label2 = "Targeted export controls")

r <- rbind(r1, r2, r3)


p2 <- ggplot(r)+
  geom_line(aes(x=time, y=rt, color=label2), size=1)+
  geom_line(aes(x=time, y=rt, color=label2), size=1)+
  theme_bw(base_size=14)+
  xlab("Time")+
  ylab(expression(r[j%->%k](t)))+
  #geom_vline(xintercept = 30, linetype="dashed")+
  scale_color_manual(values=c(cols3[1] ,cols3[3],cols3[2]) )+
  ggtitle("Export rate of infected cattle from source to recipient farm")+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle=0, vjust=0.5, size=16),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, color = "black", size=13))


get_mu <- function(inc_t, r_t, phi_tau, scale){
  
  inc_t <- c(rep(0,40), inc_t)
  r_t <- c(rep(0,40), r_t)
  
  int1 <- c()
  int2 <- rep(0,40)
  
  for(i in 41:length(inc_t)){
    
    int2 <- c(int2, sum(inc_t[(i-39):i]*rev(phi_tau)))
    int1 <- c(int1, sum(r_t[1:i]*int2[1:i]))
    
  }
  mu_t <- int1
  p_t <- exp(-mu_t[2:length(mu_t)])*(mu_t[2:length(mu_t)]-mu_t[1:(length(mu_t)-1)])
  
  data.frame(time=seq(1,100),
             mu_t = int1,
             p_t_pdf = c(0, p_t),
             p_t_cdf = 1- exp(-mu_t))
}
sigmoid <- function(x, a){
  1/(1+exp(-a*x))
}

1-sigmoid(seq(-4,36), a=1)
plot(1-sigmoid(seq(-5,35), a=1))

phi_tau1 <- (1-sigmoid(seq(-7,32), a=1))*0.5
phi_tau2 <- (1-sigmoid(seq(-7,32), a=1))*0.2


mu_t1a <- get_mu(inc_t = sim1$inc[1:100], r_t = r1$rt[1:100]*0.001, phi_tau = phi_tau1)
mu_t2a <- get_mu(inc_t = sim2$inc[1:100], r_t = r3$rt[1:100]*0.001, phi_tau = phi_tau1)
mu_t3a <- get_mu(sim3$inc[1:100], r2$rt[1:100]*0.001, phi_tau = phi_tau1)

mu_t1b <- get_mu(sim1$inc[1:100], r1$rt[1:100]*0.001, phi_tau = phi_tau2)
mu_t2b <- get_mu(sim2$inc[1:100], r3$rt[1:100]*0.001, phi_tau = phi_tau2)
mu_t3b <- get_mu(sim3$inc[1:100], r2$rt[1:100]*0.001, phi_tau = phi_tau2)

mu_t1a$label = "Unmitigated outbreak"
mu_t1a$label2 = "No pre-emptive\ninterventions"

mu_t2a$label = "Pre-emptive interventions"
mu_t2a$label2 = "No pre-emptive\ninterventions"

mu_t3a$label = "Reactive interventions"
mu_t3a$label2 = "No pre-emptive\ninterventions"

mu_t1b$label = "Unmitigated outbreak"
mu_t1b$label2 = "Pre-emptive\ninterventions"

mu_t2b$label = "Pre-emptive interventions"
mu_t2b$label2 = "Pre-emptive\ninterventions"

mu_t3b$label = "Reactive interventions"
mu_t3b$label2 = "Pre-emptive\ninterventions"

square_times <- function(mu_t){
  mu_ta <- mu_t
  mu_tb <- mu_t
  mu_tb$time <- mu_t$time+0.9999999
  
  rbind(mu_ta)
}

df_mu <- rbind(square_times(mu_t1a), square_times(mu_t1b),
               square_times(mu_t2a), square_times(mu_t2b),
               square_times(mu_t3a), square_times(mu_t3b) )


p3 <- ggplot(df_mu)+
  geom_line(aes(x=time, y=p_t_cdf, color=label, linetype=label2), size=1)+
  theme_bw(base_size=14)+
  xlab("Time")+
  ylab(expression(P[j%->%k](t)))+
  scale_color_manual("Source farm",values=cols3, breaks = c("Unmitigated outbreak", "Pre-emptive interventions", "Reactive interventions"))+
  scale_linetype_manual("Recipient farm",values=c("solid", "dotted"))+
  ggtitle("Probability of outbreak on recipient farm")+
  #geom_vline(xintercept = 30, linetype="dashed")+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(angle=0, vjust=0.5, size=16),
        legend.position = c(0.18,0.78),
        legend.background = element_rect(color="black"),
        plot.title = element_text(hjust = 0.5, color = "black", size=13))+
  guides(linetype=guide_legend(order=1),
         color=guide_none())

p3
df_phi1 <- data.frame(time=seq(0,39),
                      phi = phi_tau1,
                      label2 = "No pre-emptive interventions" )

df_phi2 <- data.frame(time=seq(0,39),
                      phi = phi_tau2,
                      label2 = "Pre-emptive interventions" )

df_phi <- rbind(df_phi1, df_phi2)

p4 <- ggplot(df_phi)+
  geom_line(aes(x=time, y=phi,linetype=label2), size=1)+
  theme_bw(base_size=14)+
  xlab(expression("Time since infection, "~tau) )+
  ylab(expression(phi[k](tau)))+
  #geom_vline(xintercept = 30, linetype="dashed")+
  scale_linetype_manual("Recipient farm",values=c("solid", "dotted"))+
  ggtitle("Probability an imported infected cow\ncauses outbreak on recipient farm")+
  scale_y_continuous(position='left')+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(angle=0, vjust=0.5),
        legend.position = c(0.67,0.8),
        legend.background = element_rect(color="black"),
        plot.title = element_text(hjust = 0, color = "black", size=14))+
  coord_cartesian(xlim=c(2,15)) + theme(plot.background = element_rect(color="black"))

p4
ggsave(paste('figure/',  'figure1-SUP.png', sep=""), width=5, height=5)





p1 <- p1+theme(axis.title.x = element_blank())+coord_cartesian(xlim=c(2,95))
p2 <- p2+coord_cartesian(xlim=c(2,95))
p3 <- p3+coord_cartesian(xlim=c(2,95))





p1 + p2 + p3 +plot_layout(nrow = 3, ncol=1, heights=c(2,1, 2))


ggsave(paste('figure/',  'figure1.png', sep=""), width=6, height=9)


