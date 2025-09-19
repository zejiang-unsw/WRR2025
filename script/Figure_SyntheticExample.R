# Experiment on phases
rm(list=ls()) # remove all variables
graphics.off() # remove all figures

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))

flag.save <- F

# load library ----
library(data.table)
library(ggplot2)
library(egg)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(doParallel)
library(foreach)

#display.brewer.all(colorblindFriendly = T, type=c("div","qual","seq","all")[2])
col_len <- 8
col_qual <- brewer.pal(n=col_len, name="Dark2")
col_div <- brewer.pal(n=col_len, name="RdYlBu")
col_seq <- brewer.pal(n=col_len, name="Blues")

# global variables
mode <- switch(1, "DWT","MODWT")
wf <- switch(3, "haar","d4","d16")

flag.v <- switch(1, "","_v1")
# synthetic data ----
set.seed(2025-07-08)
ensemble <- 100

fd <- c("x","y")
n <- 5000
noise <- 0.1

ts.r <- WASP2.0::data.gen.Rossler(a = 0.2, b = 0.2, w = 5.7, start = c(-2, -10, 0.2), n, s=0)
data <- list(x=ts.r$z, dp=cbind(ts.r$x, ts.r$y))
plot.ts(cbind(data$x,data$dp), main="",plot.type = c("multiple", "single")[2], col=col_qual)


df.data <- data.frame(No=1:n, x=ts.r$x, y=ts.r$y, z=ts.r$z) %>% gather(group,value,2:3)

font_size = 8
p <- ggplot(data=df.data) +

  geom_line(aes(x=No, y=value, color=group),linewidth=0.5) +
  geom_line(aes(x=No, y=z-20, color="z"),linewidth=0.5) +

  #facet_wrap(group~., ncol=2, scales = "fixed", strip.position = "top") +
  
  #scale_x_continuous(limits = c(1, 60), breaks = scales::pretty_breaks(n=6)) +
  scale_y_continuous( breaks = scales::pretty_breaks(n=6)) +
  
  scale_color_manual(values = col_qual) +
  scale_fill_manual(values = col_qual) +
  
  #coord_cartesian(clip = "off") + # label outside margin
  
  labs(x =NULL, y= NULL, #y = expression(bold("Forecast Skill"~(rho))),
       color=NULL, fill=NULL) +
  #guides(color = "none") +
  egg::theme_article() +
  theme(text = element_text(size = font_size),
        plot.margin = unit(c(0.5,0.1,0.1, 0.1), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        
        #panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0,size = 12,
                                         margin = margin(t = 0, r = 5, b = 0, l = 0)),
        #strip.text = element_text(size = 12),
        
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        
        legend.position = "top",
        legend.title=element_blank(),
        legend.key.width = unit(1,"cm"))

p %>% print()

if(flag.save){
filen <- paste0("Figure_Rossler_",wf,"_n",n,"_sd",noise,"_0813.png")
png(filen, height=10, width=12, units="cm", bg = "white", res=600)
p %>% print()
dev.off()
}


# data.list <- NULL
# for(r in 1:ensemble){
#   if(flag.v=="_v1"){
#     z0 <- c(runif(1,-3,-1), runif(1,-11,-9), runif(1,0,0.5))
#   } else {
#     z0 <- c(-2, -10, 0.2)
#   }
#   ts.r <- WASP2.0::data.gen.Rossler(a = 0.2, b = 0.2, w = 5.7, start = z0,
#                                     n, s=noise)
#   
#   data.list[[r]] <- list(x=ts.r$z, dp=cbind(ts.r$x, ts.r$y))
#   
# }
#sapply(data.list, function(ls) is.na(ls$x)) %>% sum()
#data1 <- data.list[[sample(1:100,1)]]
#plot.ts(cbind(data1$x,data1$dp), main="")

# load results ----
# k.folds <- 2
# folds <- cut(seq(1,n),breaks=k.folds,labels=FALSE)

file_dat <- paste0("../result/Rossler_",wf,"_n",n,"_sd",noise,
                   "_r",ensemble,flag.v,".Rdata")
load(file_dat)


summary(pred_df_list)
# metrics ----
pred_dt <- pred_df_list %>% gather(model, pred, which(colnames(pred_df_list) %in% paste0("mod",1:3))) %>% 
  mutate_if(is.character, as.factor) %>% 
  data.table()
summary(pred_dt)

metric_df <- pred_dt[,.(RMSE=Metrics::rmse(obs,pred),
                        R = cor(obs,pred)), by=c("method", "group", "model","r")]
metric_df$group <- factor(metric_df$group, labels = c("Calibration","Validation"))
metric_df$model <- factor(metric_df$model, labels=c("x","x'","x'(φ)"))
summary(metric_df)

# plot ----
metric_df1 <- metric_df %>% gather(metric,value,5:6)
p_val <- ggplot(metric_df1 %>% subset(group=="Validation"&method=="lm"), 
                 aes(x=model, y=value,fill=model)) + 
  #geom_bar(position="dodge", stat="identity") + 
  geom_boxplot(linewidth = 0.3,       # thinner box lines
               outlier.size = 0.8) +  # smaller outlier points
  
  facet_wrap(metric~., scales="free") + 
  
  scale_y_continuous(breaks=scales::pretty_breaks(n=5)) + 
  scale_fill_manual(values=col_qual) + 
  
  labs(x=NULL, y=NULL, title=NULL) + 
  egg::theme_article() +
  #egg::theme_presentation() +
  theme(text = element_text(family = "sans", face='bold',size =font_size),
        plot.margin = unit(c(0.5,0.1,0.1, 0.1), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        
        strip.background.x = element_rect(fill="#DEDEDE", colour="black"),
        strip.background.y = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = font_size,family="sans", face='plain'),
        axis.text = element_text(size=font_size, family="sans", face="plain",color="black"),
        axis.text.x = element_text(size=font_size, family="sans", face="bold",color="black"),
        
        legend.position = "none",
        legend.title = element_text(family = "sans", face='bold',size = font_size, hjust=5),
        legend.background = element_rect(fill = "transparent"),
        legend.spacing.x = unit(0., 'cm'),
        legend.key.width = unit(1.5, "cm"),
        legend.text = element_text(size=font_size, margin = margin(l = -10))
  )

p_val %>% print()

if(flag.save){
filen <- paste0("Figure_Rossler_",wf,"_n",n,"_sd",noise,"_r",ensemble,
                flag.v,"_0813.png")
png(filen, height=10, width=12, units="cm", bg = "white", res=600)
p_val %>% print()
dev.off()

# export::graph2pdf(x=fig, file="Figure1_demo", aspectr=2, font = "Arial",
#                   width = 18/2.54, height = 14/2.54, bg = "transparent")
}

