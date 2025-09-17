# Experiment on phases
rm(list=ls()) # remove all variables
#graphics.off() # remove all figures

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))

# change the subplots of figure a,b,c, and d into two columns, one column x with high frequency
# mimicing flood, another x with low frequency mimicing drought, then show the phase shift figb,
# then do the variance transformation figc show the variace distribution with transformed timeseirs
# last, do the scatter plot figd

# load library ----
library(WASP2.0)
library(data.table)
library(ggplot2)
library(egg)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(ggpmisc)

#display.brewer.all(colorblindFriendly = T, type=c("div","qual","seq","all")[2])
col_len <- 8
col_qual <- brewer.pal(n=col_len, name="Dark2")
col_div <- brewer.pal(n=col_len, name="RdYlBu")
col_seq <- brewer.pal(n=col_len, name="Blues")

font_size <- 10

# synthetic data ----
set.seed(101)
#frequency, sampled from a given range
fp <- 25
n=nobs = 300
sd.y = sd.x = 0
fs = 1/nobs
t <- seq(0, 1, length.out = nobs)
y <- sin(2 * pi * fp * t ) # + rnorm(nobs, 0, sd.y)

phi1 = pi/2
x1 <- sin(2 * pi * fp * t + phi1) + rnorm(nobs, 0, sd.x)

data <- list(x = y, dp = as.matrix(x1))
#plot.ts(cbind(data$x,data$dp), main="")

# phase transformation of raw data
dp.n <- phase.wasp.index(data$dp, data$x, index=1:n, k=7)
#dp.n <- phase.wasp.index(data$dp, data$x, index=1:(n/2),k=5)


### fig-a ----
df_raw <- data.frame(No=1:n, cbind(y=y, X1=x1, X2=dp.n[,1])) %>% gather(Group,value,3:4)

my_labeller <- as_labeller(c(X1="x[1]",
                             X2="x[2]",
                             y="y [Truth]"),
                           default = label_parsed)
df_raw$Group <- factor(df_raw$Group, levels = c("y","X1","X2"))
summary(df_raw)

p1a <- ggplot(df_raw) +
  
  geom_line(aes(x=No, y=y), color="black",linewidth=0.5) +
  geom_line(aes(x=No, y=value, color=Group),linewidth=0.5) +
  
  facet_wrap(Group~., strip.position="left", ncol = 1, labeller = my_labeller) +
  
  #scale_y_continuous(limits=c(-2.5,2.5)) +
  scale_x_continuous(limits=c(0,100)) + 
  scale_color_manual(values=col_qual) +
  labs(x=NULL, y=NULL) +
  
  #egg::theme_article() +
  theme(text = element_text(size = font_size),
        plot.margin = unit(c(0.8,0.1,0.5,0.2), "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        #panel.background = element_rect(fill=alpha(col_qual[8],0.15)),
        plot.background= element_blank(),
        
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0,size = 12,
                                         margin = margin(t = 0, r = 5, b = 0, l = 0)),
        #strip.text = element_text(size = 12),
        
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        
        legend.position = "none",
        legend.title=element_blank(),
        legend.key.width = unit(1,"cm"))

p1a %>% print()

# WASP
mode <- "DWT"
wf <- "d6"

dwt0 <- WASP2.0::dwt.vt(list(x=y, dp=y), wf, verbose = F)
dwt <- WASP2.0::dwt.vt(data, wf, verbose = F)
dwt1 <- WASP2.0::dwt.wasp.val(data, wf, phase=TRUE, k=5, max_iter=3, verbose = F)

# Auto-k + momentum demo ----
dwt1_auto <- WASP2.0::dwt.wasp.val(
  data, wf,
  phase = TRUE,
  k_auto = TRUE, power_thresh = 0.8, k_max = 10,
  momentum = 0.3,
  max_iter = 3,
  verbose = FALSE
)
cat("\n[Auto-k demo] k_used per predictor:\n")
print(sapply(dwt1_auto$diag$freq_idx, length))
cat("Top freq_idx of first predictor (head):\n")
print(head(dwt1_auto$diag$freq_idx[[1]]))

## fig-b ----
# before normalization
data.frame(No=1:nrow(dwt$S), S0=dwt0$S, S1=dwt$S[,1], 
           S2=dwt1$S[,1], S3=dwt1_auto$S[,1]) %>% print()

df_cov <- data.frame(No=1:nrow(dwt$S), S0=dwt0$S/norm(dwt0$S, type = "2"), 
                     S1=dwt$S[,1]/norm(dwt$S[,1], type = "2"),
                     S2=dwt1$S[,1]/norm(dwt1$S[,1], type = "2"),
                     S3=dwt1_auto$S[,1]/norm(dwt1_auto$S[,1], type = "2")) %>% 
  tidyr::gather(group, S, 3:5) %>% mutate(S=round(S,3), S0=round(S0,3))
# after normalization
df_cov %>% spread(group, S) %>% print()

summary(df_cov)
p1b <- ggplot(df_cov,aes(x=No, y=S,fill=group)) +
  geom_bar(position="dodge", stat="identity") +
  
  geom_line(aes(x=No, y=S0), color="black",linewidth=1) +  
  geom_point(aes(x=No, y=S0)) + 
  
  facet_wrap(group~., scales="free",ncol=1) +
  
  #scale_y_continuous(limits = c(0,5)) +
  scale_fill_manual(values=col_qual) +
  #coord_fixed(ratio=6) +
  
  labs(x="Frequency",y="Covariance") +
  
  #egg::theme_article() +
  theme(text = element_text(size = font_size),
        plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        #panel.background = element_rect(fill=alpha(col_qual[8],0.15)),
        plot.background= element_blank(),
        
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0,
                                         margin = margin(t = 0, r = 12, b = 0, l = 0)),
        strip.text = element_blank(),
        
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        
        axis.title.x = element_text(size=10, margin=margin(t=2)), #add margin to x-axis title
        axis.title.y = element_text(size=10, margin=margin(r=2)), #add margin to y-axis title
        
        legend.position = "none",
        legend.title=element_blank(),
        legend.key.width = unit(1,"cm"))

p1b #%>% print()


## fig-c ----
df_x_n <- data.frame(No=1:n, y=dwt$x, X1=dwt$dp.n[,1], X2=dwt1$dp.n[,1]) %>% 
  tidyr::gather(group, x, 3:4)

my_labeller <- as_labeller(c(X1="x[1]^\"'\"",
                             X2="x[2]^\"'\"",
                             y="y [Truth]"),
                           default = label_parsed)


summary(df_x_n)
p1c <- ggplot(df_x_n) +
  
  geom_line(aes(x=No, y=y), color="black") + 
  geom_line(aes(x=No, y=x, color=group)) + 
  
  facet_wrap(group~., scale="free",strip.position="left", ncol = 1, labeller = my_labeller) +
  
  #scale_y_continuous(limits = c(0,5)) +
  scale_x_continuous(limits=c(0,100)) + 
  scale_color_manual(values=col_qual) +
  labs(x=NULL, y=NULL) +
  
  theme(text = element_text(size = font_size),
        plot.margin = unit(c(0,0.1,0,0.2), "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        #panel.background = element_rect(fill=alpha(col_qual[8],0.15)),
        plot.background= element_blank(),
        
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, size=12, 
                                         margin = margin(t = 0, r = 5, b = 0, l = 0)),
        #strip.text = element_text(size = font_size),
        
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        
        legend.position = "none",
        legend.title=element_blank(),
        legend.key.width = unit(1,"cm"))

p1c #%>% print()

## fig-d----
fit1 <- stats::lm(dwt$x ~., data= dwt$dp.n[,1]%>% data.frame()) # baseline
fit2 <- stats::lm(dwt1$x ~., data= dwt1$dp.n[,1]%>% data.frame()) # variance only 

pred_cal_lm <- data.frame(No=1:n, obs=dwt$x, 
                          mod1=fit1$fitted.values, mod2=fit2$fitted.values) %>% gather(Group,pred,3:4)

summary(pred_cal_lm)

p1d <- ggplot(data = pred_cal_lm, 
              aes(x = pred, y = obs, color=Group, fill=Group)) +
  #stat_poly_line(size = 1, alpha = 1, se=FALSE, color="grey") +
  stat_poly_eq(label.y = 0.9, size=3, color="black") + 
  geom_point(size = 2, alpha = 1, stroke=1, shape=21, color="grey4") + 
  geom_abline(intercept=0, slope=1, color="grey4") + 
  
  coord_fixed() + 
  facet_wrap(.~Group, scale="fixed", ncol=1) + 
  
  scale_x_continuous(limits = c(-1,1), exp=c(0.2,0)) +
  scale_y_continuous(limits = c(-1,1), exp=c(0.2,0)) +
  scale_fill_manual(values=col_qual) + 
  
  labs(x="Predicted",y="Observed") +
  #egg::theme_article() + 
  egg::theme_presentation() + 
  theme(text = element_text(size = 7, family="sans"),
        plot.margin = unit(c(0.1,0,0.1, 0), "cm"),
        plot.background = element_blank(),
        
        strip.background = element_blank(),
        strip.placement = "outside",
        #strip.text.y = element_text(angle = 0), 
        strip.text.x = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        
        axis.title.x = element_text(size=10, margin=margin(t=2)), #add margin to x-axis title
        axis.title.y = element_text(size=10, margin=margin(r=10,l=-5)), #add margin to y-axis title
        
        legend.position = "none",
        legend.background = element_rect(fill="transparent"),
        legend.title=element_blank(),
        legend.text=element_text(face = "bold"), 
        legend.key.width = unit(1,"cm")
  )

p1d #%>% print()

fig <- cowplot::plot_grid(p1a,  p1b, p1c, p1d, ncol=2, #labels=letters,
                          rel_widths = c(1,0.8),
                          #label_size= font_size + 2, hjust=1,
                          label_fontfamily = "sans",
                          label_fontface = "bold",
                          label_colour = "black")

fig %>% print()

# Figure 1 ----
if(F){
  filen <- "Figure_conept"
  export::graph2pdf(x=fig, file=filen, aspectr=2, font = "Arial",
                    height = 16/2.54, width = 12/2.54, bg = "transparent")
  
  png(filen, height=16, width=12, units="cm", bg = "transparent", res=600)
  fig %>% print()
  dev.off()
}
