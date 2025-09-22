# stationary model as the baseline model
rm(list=ls())
#graphics.off() # clear environment and graphics

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))

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

flag.save <- F

# Global parameters----
#display.brewer.all(colorblindFriendly = T, type=c("div","qual","seq","all")[2])
col_len <- 8
col_qual <- brewer.pal(n=col_len, name="Dark2")
col_pair <- brewer.pal(n=col_len+4, name="Paired")
col_div <- brewer.pal(n=col_len, name="RdYlBu")
col_seq <- brewer.pal(n=col_len, name="YlOrRd")

#col_qual <- viridisLite::viridis(n=4, option = "D")

#labels_mod <- c("Stationary","Non-Stationary","Non-stationary (WASP)")
# labels_mod_exp <- c(expression(italic(mu~" ")), expression(italic(mu(x))),
#                     expression(italic(mu("x'"))), expression(italic(mu("x'(φ)"))))
labels_mod_exp <- c(expression(mu~" "), expression(mu(x)),
                    expression(mu("x'")), expression(mu("x'(φ)")))
lead_all <- c(1:3)
season_all <- c("DJF", "MAM", "JJA", "SON")

#flag.k <- 5; flag.iter= 3
flag.k <- 3; flag.iter= 1
log.tran <- 1
flag.v <- switch(1, paste0("_k",flag.k,"_log",log.tran), "_v1", "_v2","_v3") # k=5, iter=1

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
CI_names <- c('Nino3', 'Nino4', 'Nino3.4', 'EMI', 'DMI','EPI', 'WPI', 'II', 
              'TSI', 'TBV', 'Disc')

CI_longname <-c('Niño 3', 'Niño 4', 'Niño 3.4', 'El Niño Modoki Index', 
                'Dipole Mode Index','Indian Ocean East Pole Index', 
                'Indian Ocean West Pole Index', 'Indonesian Index', 
                'Tasman Sea Index', 'Tropical Trans-basin Variability Index')

CI_levels_plot <- c('Nino3', 'Nino4', 'Nino3.4', 'EMI', 'DMI','EPI', 'WPI', 'II', 'TSI', 'TBV')
CI_labels_plot <- CI_longname[match(CI_levels_plot, CI_names)]
CI_levels_plot_rev <- rev(CI_levels_plot)

model_levels_keep <- c("Mod3", "Mod1")

create_case_plot <- function(flag.var, panel_label) {
  if(flag.var=="wet"){
    data("Qmax_SSI_AU_wet")
    Qmax_SSI_AU <- Qmax_SSI_AU_wet
  } else {
    data("Qmax_SSI_AU_dry")
    Qmax_SSI_AU <- Qmax_SSI_AU_dry
  }
  data_q <- Qmax_SSI_AU
  stn_nos <- names(data_q)
  stn_check <- stn_nos
  stn_samp <- sapply(stn_check, function(v) which(v==stn_nos))
  
  # Fit GEV-------------------------------------------------
  file_save <- paste0("../result/SSI_CRPSS_",mode,"_",wf,"_Location_",
                      method_fevd,flag.v,"_ALL_AU_",flag.var,".Rdata")
  load(file_save)
  
  ## RPSS ----
  summary(RPSS_df_all)
  RPSS_df_all1 <- RPSS_df_all %>% spread(Model, RPSS) %>% 
    mutate(Mod1=(Mod1-Mod0)/(1-Mod0)) %>% 
    mutate(Mod2=(Mod2-Mod0)/(1-Mod0)) %>% 
    mutate(Mod3=(Mod3-Mod0)/(1-Mod0)) %>%
    dplyr::select(!Mod0) %>%
    dplyr::filter(CI %in% CI_levels_plot)
  
  RPSS_df_all1$CI <- factor(RPSS_df_all1$CI, levels=CI_levels_plot)
  
  # average across two folds
  RPSS_df <- RPSS_df_all1 %>% dplyr::group_by(CI, Lead, Station) %>% 
    summarise(Mod1=mean(Mod1), Mod2=mean(Mod2), Mod3=mean(Mod3))
  summary(RPSS_df)
  
  RPSS_dt <- RPSS_df %>%  data.table()
  tol_no <- length(stn_nos)
  RPSS_pct <- RPSS_dt[,.(Mod1=sum(Mod1>0, na.rm=T)/tol_no,
                         Mod2=sum(Mod2>0, na.rm=T)/tol_no,
                         Mod3=sum(Mod3>0, na.rm=T)/tol_no), by=c("Lead","CI")]
  
  ## percentage of improvements ----
  RPSS_pct %>% print()
  
  ## fig-box ----
  RPSS_df_all1 <- RPSS_df %>% gather(Model, RPSS, 4:6) %>% 
    dplyr::filter(Model %in% model_levels_keep,
                  CI %in% CI_levels_plot,
                  RPSS > 0)
  #RPSS_df_all1 <- RPSS_df_all %>% subset(RPSS > 0)
  #RPSS_df_all1 <- RPSS_df_all 
  summary(RPSS_df_all1)
  
  RPSS_df_all1$CI <- factor(RPSS_df_all1$CI, levels=CI_levels_plot)
  RPSS_df_all1$Model <- factor(RPSS_df_all1$Model, levels=model_levels_keep)
  RPSS_df_all1$Lead <- factor(RPSS_df_all1$Lead, labels = paste0("Lead ",0:2))
  
  RPSS_df_mu <- RPSS_df_all1 %>% group_by(Lead, Model) %>% summarise(mu=mean(RPSS))
  
  ## absolute improvements ----
  #RPSS_df_mu %>% spread(Model, mu) %>% mutate(Dif=Mod2-Mod1) %>% print()
  RPSS_df_mu %>% spread(Model, mu) %>% print()
  
  # Step 1: summarise your data
  RPSS_summary <- RPSS_df_all1 %>%
    dplyr::filter(!is.na(CI)) %>%
    group_by(CI, Model, Lead) %>%
    summarise(
      mean_RPSS = mean(RPSS, na.rm = TRUE),
      sd_RPSS   = sd(RPSS, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    mutate(
      ymin = mean_RPSS - sd_RPSS,
      ymax = mean_RPSS + sd_RPSS
    )
  summary(RPSS_summary) %>% print()
  
  # Step 2: plot
  RPSS_summary <- RPSS_summary %>%
    mutate(CI = factor(CI, levels = CI_levels_plot)) %>% 
    mutate(Model = factor(Model, levels = c("Mod3", "Mod1")))

  axis_position <- ifelse(flag.var == "wet", "right", "left")
  dodge_width <- 0.2
  pd <- position_dodge2(width = dodge_width, padding = 0.15, preserve = "single")
  
  if(flag.var=="wet"){
    scale_case <- scale_x_reverse(
      #limits = c(-0.1, 0.45),
      limits = c(0.4, -0.),
      breaks = scales::pretty_breaks(n = 6),
      expand = c(0.05, 0)
    ) 
  } else {
    scale_case <- scale_x_continuous(
      limits = c(-0., 0.4),
      breaks = scales::pretty_breaks(n = 6),
      expand = c(0.05, 0)
    ) 
  }
  
  color_val = c("Mod3" = col_pair[2], "Mod1" = alpha(col_pair[1],0.7))
  plot_case <- ggplot() +
    geom_errorbar(data = RPSS_summary %>% subset(Model=="Mod3"),
                  aes(x = mean_RPSS, y = CI, xmin = mean_RPSS, xmax = ymax),
                  width = 0.6) +
    geom_col(data = RPSS_summary %>% subset(Model=="Mod3"),
             aes(x = mean_RPSS, y = CI, color = Model, 
                 fill = Model),width=0.8) +
    geom_col(data = RPSS_summary %>% subset(Model=="Mod1"),
             aes(x = mean_RPSS, y = CI, color = Model, 
                 fill = Model),width=0.5) +
    # geom_point(data = RPSS_summary, aes(y = CI, x = ymax, fill = Model),
    #            #position = pd,
    #            pch = 23,
    #            size = 2.1,
    #            stroke = 0.3) +
    geom_vline(xintercept = 0, color = "black", lwd = 0.5) +
    
    scale_fill_manual(
      values = color_val,
      labels = c(expression(mu("x'(φ)")), expression(mu(x)))
    ) + 
    scale_color_manual(
      values = color_val,
      labels = c(expression(mu("x'(φ)")), expression(mu(x)))
    ) +
    
    facet_wrap(.~Lead, ncol = 1, scales = "fixed") +
    scale_y_discrete(
      limits = CI_levels_plot_rev,
      labels = rev(CI_labels_plot),
      position = axis_position
    ) +
    scale_case + 
    labs(
      x = "Continuous Ranked Probability Skill Score (CRPSS)",
      y = NULL,
      color = NULL, fill=NULL,
      title = panel_label
    ) +
   guides(color = guide_legend(reverse = TRUE),
          fill = guide_legend(reverse = TRUE)) +
    
    egg::theme_article() +
    theme(
      text = element_text(size = 7, family = "sans"),
      plot.margin = unit(c(0.1, 0.5, 0.1, 0.5), "cm"),
      plot.title = element_text(hjust = 0, face = "bold"),
      strip.background = element_rect(fill = "#DDDDDD", colour = "grey4"),
      strip.placement = "outside",
      strip.text.x = element_text(size = 7, family = "sans", face = 'bold'),
      axis.text.x = element_text(size = 7, family = "sans", face = "plain", color = "black"),
      axis.text.y = element_text(size = 7, family = "sans", face = "plain", color = "black", hjust=0.5,
                                 margin = margin(r = 20)),
      #axis.text.y.right = element_text(size = 7, family = "sans", face = "plain", color = "black",hjust=0.5),
      axis.text.y.right = element_blank(),
      axis.title.x = element_text(margin = margin(t = 5)),
      axis.ticks.length  = unit(0.1, "cm"),
      axis.ticks = element_line(linewidth = 0.2),
      legend.position = c(1.5, -0.02),
      legend.title = element_text(hjust = 0.8),
      legend.direction = "horizontal",
      legend.background = element_rect(fill = "transparent"),
      legend.key.spacing.x = unit(0.1, "cm"),
      legend.key.width = unit(0.8, "cm"),
      legend.text = element_text(size = 7, margin = margin(r = 1))
    )

  if(flag.var == "wet"){
    plot_case <- plot_case +
      theme(
        axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.title.y = element_blank()
      )
  } else {
    plot_case <- plot_case +
      theme(
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
      )
  }

  plot_case
}

p_wet <- create_case_plot("wet", "a) SSI 7-day maximum")
p_dry <- create_case_plot("dry", "b) SSI 30-day minimum") +
  theme(legend.position = "none")

fig <- egg::ggarrange(p_wet, p_dry, ncol = 2)

fig %>% print()


### Figure 3 ----
if(flag.save){
  filen <- paste0("Figure_CRPSS_ALL_",mode,"_",wf,"_",method_fevd,flag.v,"_AU_SSI_wet_dry.png")
  
  export::graph2pdf(x = fig, file = filen, aspectr = 2, font = "Arial",
                    height = 18 / 2.54, width = 18 / 2.54, bg = "white")
  
  png(filen, height=10,width=16,units="cm", bg = "white", res=600)
  grid::grid.newpage()
  grid::grid.draw(fig)
  dev.off()
}
