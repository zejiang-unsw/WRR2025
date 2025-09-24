# stationary model as the baseline model
rm(list=ls())
#graphics.off() # clear environment and graphics

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))

flag.save <- F

# Load required packages
#library(GPST)
library(WASP2.0)
library(hydroTSM) #time2season
library(sf)
library(ggplot2)
library(extRemes)
library(dplyr)
library(tidyr)
library(data.table)
library(export)
library(RColorBrewer)

library(scoringRules)

library(doParallel)
library(parallel)
library(foreach)
# Global parameters----
#display.brewer.all(colorblindFriendly = T, type=c("div","qual","seq","all")[2])
col_len <- 8
col_qual <- brewer.pal(n=col_len, name="Paired")
col_div <- brewer.pal(n=col_len, name="RdYlBu")
col_seq <- brewer.pal(n=col_len, name="YlOrRd")

labels_mod_exp <- c(expression(italic(mu~" ")), expression(italic(mu(x))),
                    expression(italic(mu("x'"))), expression(italic(mu("x'(φ)"))))

lead_all <- c(1:3)
#season_all <- c("DJF", "MAM", "JJA", "SON")
flag.par <- T
flag.v <- switch(1, "","_k5", "_k3", "_v2","_v3")
log.tran = 1

flag.v1 <- 1
flag.var <- switch(flag.v1, "wet","dry")
#flag.acc <- switch(flag.v1, 7, 30)

# v0 use SSI_k as the covariate
# v1 use Disc as the covariate 
flag.s <- switch(2,"","v1")

## WASP ----
mode <- switch(1, "DWT","MODWT") ## must use MODWT
if(mode=="DWT") wf <- switch(3, "haar","d4","d16") else wf <- "d4"

if(mode=="DWT") {method = "modwt"; pad = "sym"} # !!symmetric padding very important!!
boundary = switch(1,"periodic","reflection"); sign = "pos" # !! pos matters!!

## GEV -----
sig_lev <- 0.05 # likelihood ratio test
method_fevd <- switch(1, "MLE","GMLE") # !!! reliable parameter estimation is important !!!
#optim_args <- list(method=c("BFGS", "Nelder-Mead")[1])

# Load data---------------------------------------------------------------------
# Processed dataset
load(paste0("../data/Stn_SSI_AU_wet_",flag.s,".Rdat"))
load(paste0("../data/Stn_SSI_AU_dry_",flag.s,".Rdat"))

stn_list <- intersect(stn_list_dry, stn_list_wet)

if(flag.var=="wet"){
  #data("Qmax_SSI_AU_wet")
  load(paste0("../data/Qmax_SSI_AU_wet_",flag.s,".Rdat"))
  Qmax_SSI_AU <- Qmax_SSI_AU_wet[stn_list]
} else {
  #data("Qmax_SSI_AU_dry")
  load(paste0("../data/Qmax_SSI_AU_dry_",flag.s,".Rdat"))
  Qmax_SSI_AU <- Qmax_SSI_AU_dry[stn_list]
}
stn_nos <- names(Qmax_SSI_AU)
summary(Qmax_SSI_AU[[1]])

CI_names <- c('Nino3', 'Nino4', 'Nino3.4', 'EMI', 'DMI','EPI', 'WPI', 'II', 
              'TSI', 'TBV', 'Disc')
CI_names1 <- c('Niño 3', 'Niño 4', 'Niño 3.4', 'EMI', 'DMI','EPI', 'WPI', 'II', 
              'TSI', 'TBV', 'Disc')
CI_longname <-c('Niño 3', 'Niño 4', 'Niño 3.4', 'El Niño Modoki Index', 
                'Dipole Mode Index','Indian Ocean East Pole Index', 
                'Indian Ocean West Pole Index', 'Indonesian Index', 
                'Tasman Sea Index', 'Tropical Trans-basin Variability Index')

stn_check <- stn_nos#[1:3]
stn_samp <- sapply(stn_check, function(v) which(v==stn_nos))

# start parallel-------------------------------------------------
system.time({
  if(flag.par){
    no_cores <- detectCores() - 4
    #cl <- makeCluster(no_cores, type="FORK")
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
    
    # Load libraries and data outside the loop whenever possible.
    clusterExport(cl, varlist=c("Qmax_SSI_AU","method", "pad", "boundary", "sign", "stn_nos"), envir=environment())
    clusterEvalQ(cl, {
      library(WASP2.0)
      library(extRemes)
      library(dplyr)
      library(readr)
    }) # Load once on all workers
    
    # on.exit({
    #   stopCluster(cl)
    #   registerDoSEQ()
    # })
  } else {
    registerDoSEQ() # test setup
  }
})



# Fit GEV-------------------------------------------------
file_name <- paste0("../result/SSI_Phase_",mode,"_",wf,flag.v,"_ALL_AU_",flag.var,".Rdata")

if(!file.exists(file_name)){
  phase_df_all <- NULL

  out_list <- list()

  system.time({out_list <- foreach(i_lead=lead_all, .combine=c) %:%
    foreach(stn_id=iter(stn_samp, by=11), .combine = c) %dopar% {

      phase.diff <- function(dp, y, index_cal, k = 1) {
        dp <- as.matrix(dp)
        n <- NROW(dp)
        if(length(y)!=n) stop("dp and y length are different!")
        if(k>(floor(n/2)-1)) stop("k must be smaller than half of data length!")
        out <- matrix(NA, nrow = n, ncol = NCOL(dp))
        
        # Calibration + Validation indices
        index_val <- setdiff(1:n, index_cal)
        
        # Demean y
        y <- y - mean(y)
        y_cal <- y
        y_cal[index_val] <- 0
        
        # FFT of calibration-phase y
        Y <- fft(y_cal)
        phase_y <- Arg(Y)
        power_y <- Mod(Y)^2
        
        # Top-k dominant frequencies (exclude DC)
        freq_idx <- order(power_y[2:floor(n / 2)], decreasing = TRUE)[1:k] + 1
        # if (n %% 2 == 0) {
        #   nyquist_idx <- n / 2 + 1
        #   if (!(nyquist_idx %in% freq_idx)) {
        #     freq_idx <- c(freq_idx, nyquist_idx)
        #   }
        # }
        #cat("\n",freq_idx)
        #cat("\n Y_fft: ",Y[c(1, freq_idx)])
        # Loop over predictors
        phase_diff_out <- vector(length=ncol(dp))
        for (i in 1:ncol(dp)) {
          # Step 1: CALIBRATION ----
          x_cal <- rep(0, n)
          x_cal[index_cal] <- dp[index_cal, i]
          X_cal <- fft(x_cal)
          
          power_x_cal <- Mod(X_cal)^2
          phase_x_cal <- Arg(X_cal)
          #cat("\n",which(power_x_cal[1:floor(n / 2)]==max(power_x_cal[1:floor(n / 2)])))
          
          phase_diff <- ((phase_y - phase_x_cal + pi) %% (2 * pi)) - pi
          #cat("\n phase diff: ",phase_diff[c(1,freq_idx)])
          #cat("\n X_fft: ",X_cal[c(1, freq_idx)])
          
          phase_diff_out[i] <- phase_diff[freq_idx]
        }
        
        return(list(phase_diff=phase_diff_out, freq_idx=rep(freq_idx,ncol(dp))))
      }
      
      Qmax_lead <- lapply(stn_id, function(i) Qmax_SSI_AU[[i]] %>% mutate(Lead=as.numeric(Lead)) %>% 
        subset(Lead==i_lead) %>% arrange(Date) %>% data.frame() %>% mutate(Disc=log(Disc+log.tran)))
      
      Qmax <- lapply(Qmax_lead, function(ls) list(Date = ls$Date,
                        x = ls$Qmax,
                        dp = ls[CI_names]
                        ))
      summary(Qmax[[1]]$dp)
      #head(Qmax[[1]]$dp)

      # k-fold cross validation ------------------------------------------------
      k.folds <- 2
      folds <- lapply(Qmax, function(ls) cut(seq(1,length(ls$x)),breaks=k.folds,labels=FALSE))

      lr_lead <- NULL
      for(k_fold in 1:k.folds){
        # k_fold <- 1
        ind_sub <- lapply(folds, function(x) which(x!=k_fold, arr.ind=TRUE)) # calibration index
        ind_val <- lapply(folds, function(x) which(x==k_fold, arr.ind=TRUE)) # validation index
        

        phase.df <- lapply(1:length(Qmax), function(i) phase.diff(Qmax[[i]]$dp, Qmax[[i]]$x, 
                                                                      ind_sub[[i]]))
       
        lr_lead <- rbind(lr_lead, data.frame(fold=k_fold, cpy=CI_names, lead=i_lead, Station=stn_nos[stn_id], 
                                             dif=sapply(phase.df, function(ls) ls$phase_diff), 
                                             freq=sapply(phase.df, function(ls) ls$freq)))
        
        
      }
    
      out <- list(lr_lead=lr_lead)
      list(out)
  }
  })

  # Combine each component from the list of lists
  lr_lead   <- do.call(rbind, lapply(out_list, function(x) x$lr_lead))

  phase_df_all <- rbind(phase_df_all, data.frame(lr_lead))
  
  ### save----
  save(phase_df_all, file=file_name)
} else {
  load(file_name)
}

if(flag.par) stopCluster(cl)
registerDoSEQ()

foreach::getDoParName() %>% print()  # (e.g., "doSEQ", "doParallel", "doSNOW").
foreach::getDoParWorkers() %>% print()  # number of workers currently active

summary(phase_df_all)
phase_df_all$fold <- factor(phase_df_all$fold)

phase_df_all1 <- phase_df_all %>% select(!freq) %>% spread(fold,dif)

phase_df_all2 <- phase_df_all %>% select(!dif) %>% spread(fold,freq)
summary(phase_df_all2)

phase_df_all1$cpy <- factor(phase_df_all1$cpy, levels=CI_names, labels = CI_names1)
phase_df_all1$lead <- factor(phase_df_all1$lead, labels = paste0("Lead ",0:2))

font_size <- 6

p_scatter_all <- ggplot(data=phase_df_all1, aes(x=`1`, y=`2`, fill=lead)) +
  
  geom_point(size = 1, alpha = 1, stroke=1, shape=21, color="grey4")+
  
  geom_abline(aes(intercept = 0, slope = 1), color="red") + 
  
  scale_color_manual(values=col_qual) +
  coord_fixed() + 
  facet_grid(lead~cpy, scales ="fixed") +
  #facet_wrap(.~lead, ncol=3, scales ="fixed") +
  
  scale_x_continuous(limits=c(-pi,pi), breaks = scales::pretty_breaks(n = 4), expand=c(0.2,0)) +
  scale_y_continuous(limits=c(-pi,pi), breaks = scales::pretty_breaks(n = 4), expand=c(0.2,0)) +
  labs(x="Δφ in Fold 1",y="Δφ in Fold 2", color="GEV model") +
  # guides(fill = guide_legend(title = "Lead",
  #                            label.position = "bottom",
  #                            title.position = "left", title.vjust = 1)) +
  
  guides(color = guide_legend(reverse = TRUE)) +
  
  egg::theme_article() +
  theme(text = element_text(size = font_size,family="sans"),
        plot.margin = unit(c(0.1,0.5,0.1,0.1), "cm"),
        
        #panel.spacing.y = unit(1, "lines"), # space between subplot
        #strip.background = element_blank(),
        strip.background = element_rect(fill="#DDDDDD", colour="grey4"),
        strip.placement = "outside",
        strip.text = element_text(size = font_size,family="sans", face='bold'),
        #strip.text.y = element_text(face='bold'),
        #strip.text.x.top = element_blank(),
        axis.text = element_text(size=font_size, family="sans", face="plain",color="black"),
        axis.title= element_text(size=font_size, family="sans", face="plain"),
        
        axis.title.x = element_text(margin = margin(t = 5)), 
        
        axis.ticks.length  = unit(0.1, "cm"),
        axis.ticks = element_line(linewidth = 0.2),
        
        legend.position = "none",
        #legend.position = c(-0.35,0.15),
        legend.title = element_text(hjust = 0.8),
        legend.direction = c("horizontal","vertical")[2],
        legend.background = element_rect(fill="transparent"),
        legend.key.spacing.x =  unit(0.1, "cm"),
        legend.key.width =  unit(0.8, "cm"),
        legend.text = element_text(size=font_size, margin = margin(r = 1))
        
  )
p_scatter_all %>% print()

## scater ----
p_scatter <- ggplot(data=phase_df_all1 %>% subset(cpy=="Disc"), aes(x=`1`, y=`2`, fill=lead)) +
  
  geom_point(size = 2, alpha = 1, stroke=1, shape=21, color="grey4")+
  
  geom_abline(aes(intercept = 0, slope = 1), color="red") + 
  
  scale_color_manual(values=col_qual) +
  coord_fixed() + 
  #facet_grid(cpy~lead, scales ="fixed") +
  facet_wrap(.~lead, ncol=3, scales ="fixed") +
  
  #scale_x_discrete(labels=CI_longname) +
  #scale_y_continuous(limits=c(-0.04,0.1), breaks = scales::pretty_breaks(n = 4), expand=c(0.1,0)) +
  labs(x="Δφ in Fold 1",y="Δφ in Fold 2", color="GEV model") +
  # guides(fill = guide_legend(title = "Lead",
  #                            label.position = "bottom",
  #                            title.position = "left", title.vjust = 1)) +
  
  guides(color = guide_legend(reverse = TRUE)) +
  
  egg::theme_article() +
  theme(text = element_text(size = 7,family="sans"),
        plot.margin = unit(c(0.1,0.5,0.1,0.5), "cm"),
        
        #panel.spacing.y = unit(1, "lines"), # space between subplot
        #strip.background = element_blank(),
        strip.background = element_rect(fill="#DDDDDD", colour="grey4"),
        strip.placement = "outside",
        strip.text.x = element_text(size = 7,family="sans", face='bold'),
        #strip.text.y = element_text(face='bold'),
        #strip.text.x.top = element_blank(),
        axis.text = element_text(size=7, family="sans", face="plain",color="black"),
        axis.title= element_text(size=7, family="sans", face="plain"),
        
        axis.title.x = element_text(margin = margin(t = 5)), 
        
        axis.ticks.length  = unit(0.1, "cm"),
        axis.ticks = element_line(linewidth = 0.2),
        
        legend.position = "none",
        #legend.position = c(-0.35,0.15),
        legend.title = element_text(hjust = 0.8),
        legend.direction = c("horizontal","vertical")[2],
        legend.background = element_rect(fill="transparent"),
        legend.key.spacing.x =  unit(0.1, "cm"),
        legend.key.width =  unit(0.8, "cm"),
        legend.text = element_text(size=7, margin = margin(r = 1))
        
  )
p_scatter %>% print()

phase_df_all$lead <- factor(phase_df_all$lead, labels = paste0("Lead ",0:2))
phase_df_all$cpy <- factor(phase_df_all$cpy, levels=CI_names, labels = CI_names1)

p_box_all <- ggplot(data=phase_df_all, aes(x=fold, y=dif, fill=lead)) +
  
  geom_boxplot(linewidth = 0.3,       # thinner box lines
               outlier.size = 0.8) +  # smaller outlier points
  
  #geom_hline(data=RPSS_df_mu, aes(yintercept = mu, color=Model),linetype="longdash", lwd=0.2) +
  
  scale_color_manual(values=col_qual) +
  
  facet_grid(lead~cpy, scales ="fixed") +
  #facet_wrap(.~lead, ncol=3, scales ="fixed") +
  
  scale_x_discrete(labels=c("Fold 1", "Fold 2")) +
  scale_y_continuous(limits=c(-pi,pi), breaks = scales::pretty_breaks(n = 4), expand=c(0.2,0)) +
  labs(x=NULL,y=expression("Δφ"), color=NULL) +
  # guides(fill = guide_legend(title = "Lead",
  #                            label.position = "bottom",
  #                            title.position = "left", title.vjust = 1)) +
  
  guides(color = guide_legend(reverse = TRUE)) +
  
  egg::theme_article() +
  theme(text = element_text(size = font_size,family="sans"),
        plot.margin = unit(c(0.1,0.5,0.1,0.1), "cm"),
        
        #panel.spacing.y = unit(1, "lines"), # space between subplot
        #strip.background = element_blank(),
        strip.background = element_rect(fill="#DDDDDD", colour="grey4"),
        strip.placement = "outside",
        strip.text = element_text(size = font_size,family="sans", face='bold'),
        #strip.text.y = element_text(face='bold'),
        #strip.text.x.top = element_blank(),
        axis.text = element_text(size=font_size, family="sans", face="plain",color="black"),
        axis.title= element_text(size=font_size, family="sans", face="plain"),
        
        axis.title.x = element_text(margin = margin(t = 5)), 
        
        axis.ticks.length  = unit(0.1, "cm"),
        axis.ticks = element_line(linewidth = 0.2),
        
        legend.position = "none",
        #legend.position = c(-0.35,0.15),
        legend.title = element_text(hjust = 0.8),
        legend.direction = c("horizontal","vertical")[2],
        legend.background = element_rect(fill="transparent"),
        legend.key.spacing.x =  unit(0.1, "cm"),
        legend.key.width =  unit(0.8, "cm"),
        legend.text = element_text(size=font_size, margin = margin(r = 1))
        
  )
p_box_all %>% print()

summary(phase_df_all)

p_box <- ggplot(data=phase_df_all %>% subset(cpy=="Disc"), aes(x=fold, y=dif, fill=lead)) +
  
  geom_boxplot()+

  #geom_hline(data=RPSS_df_mu, aes(yintercept = mu, color=Model),linetype="longdash", lwd=0.2) +
  
  scale_color_manual(values=col_qual) +
  
  #facet_grid(cpy~lead, scales ="fixed") +
  facet_wrap(.~lead, ncol=3, scales ="fixed") +

  scale_x_discrete(labels=c("Fold 1", "Fold 2")) +
  scale_y_continuous(limits=c(-pi,pi), breaks = scales::pretty_breaks(n = 4), expand=c(0.1,0)) +
  labs(x=NULL,y=expression("Δφ"), color=NULL) +
  # guides(fill = guide_legend(title = "Lead",
  #                            label.position = "bottom",
  #                            title.position = "left", title.vjust = 1)) +
  
  guides(color = guide_legend(reverse = TRUE)) +
  
  egg::theme_article() +
  theme(text = element_text(size = 7,family="sans"),
        plot.margin = unit(c(0.1,0.5,0.1,0.4), "cm"),
        
        #panel.spacing.y = unit(1, "lines"), # space between subplot
        #strip.background = element_blank(),
        strip.background = element_rect(fill="#DDDDDD", colour="grey4"),
        strip.placement = "outside",
        strip.text.x = element_text(size = 7,family="sans", face='bold'),
        #strip.text.y = element_text(face='bold'),
        #strip.text.x.top = element_blank(),
        axis.text = element_text(size=7, family="sans", face="plain",color="black"),
        axis.title= element_text(size=8, family="sans", face="plain"),
        
        axis.title.x = element_text(margin = margin(t = 5)), 
        
        axis.ticks.length  = unit(0.1, "cm"),
        axis.ticks = element_line(linewidth = 0.2),
        
        legend.position = "none",
        #legend.position = c(-0.35,0.15),
        legend.title = element_text(hjust = 0.8),
        legend.direction = c("horizontal","vertical")[2],
        legend.background = element_rect(fill="transparent"),
        legend.key.spacing.x =  unit(0.1, "cm"),
        legend.key.width =  unit(0.8, "cm"),
        legend.text = element_text(size=7, margin = margin(r = 1))
        
  )
p_box %>% print()

### Figure 6----
fig1 <- cowplot::plot_grid(p_scatter_all, p_box_all, ncol=1, #labels=letters,
                          rel_heights = c(1,0.8),
                          #label_size= font_size + 2, hjust=1,
                          label_fontfamily = "sans",
                          label_fontface = "bold",
                          label_colour = "black")

fig1 %>% print()


fig <- cowplot::plot_grid(p_scatter, p_box, ncol=1, #labels=letters,
                          rel_widths = c(1,0.8),
                          #label_size= font_size + 2, hjust=1,
                          label_fontfamily = "sans",
                          label_fontface = "bold",
                          label_colour = "black")

fig %>% print()


if(flag.save){
# export::graph2pdf(x=fig, file="Figure1_demo", aspectr=2, font = "Arial",
#                   width = 18/2.54, height = 14/2.54, bg = "transparent")
filen <- paste0("Figure6_AU_SSI_",flag.var,".png")
png(filen, height=14, width=16, units="cm", bg = "white", res=600)
fig %>% print()
dev.off()

filen1 <- paste0("Figure6_ALL_AU_SSI_",flag.var,".png")
png(filen1, height=12, width=18, units="cm", bg = "white", res=600)
fig1 %>% print()
dev.off()
}

#shell.exec(filen1)
