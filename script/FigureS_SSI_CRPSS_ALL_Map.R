# stationary model as the baseline model
rm(list=ls())
#graphics.off() # clear environment and graphics

# Getting the path of your current open file
if(requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()){
  current_path <- rstudioapi::getActiveDocumentContext()$path
  if(!is.null(current_path) && nzchar(current_path)){
    setwd(dirname(current_path))
  }
}

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

flag.save.supp <- F

# Global parameters----
#display.brewer.all(colorblindFriendly = T, type=c("div","qual","seq","all")[2])
col_len <- 8
col_qual <- brewer.pal(n=col_len, name="Dark2")
col_pair <- brewer.pal(n=col_len+4, name="Paired")
col_div <- brewer.pal(n=col_len, name="RdYlBu")
col_seq <- brewer.pal(n=col_len, name="YlOrRd")

col_qual <- viridisLite::viridis(n=4, option = "D")

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


target_group <- 3

create_condition_plot <- function(flag.var,
                                  panel_letters,
                                  inset_params = list(x = 0.035, y = 0.055, width = 0.38, height = 0.3)) {
  condition_title <- ifelse(flag.var == "wet", "Wet", "Dry")

  if(flag.var == "wet"){
    data("Qmax_SSI_AU_wet")
    Qmax_SSI_AU <- Qmax_SSI_AU_wet
  } else {
    data("Qmax_SSI_AU_dry")
    Qmax_SSI_AU <- Qmax_SSI_AU_dry
  }

  data_q <- Qmax_SSI_AU
  stn_nos <- names(data_q)
  CI_names <- c("Nino3", "Nino4", "Nino3.4", "EMI", "DMI","EPI", "WPI", "II",
                "TSI", "TBV", "Disc")

  file_save <- paste0("../result/SSI_CRPSS_", mode, "_", wf, "_Location_",
                      method_fevd, flag.v, "_ALL_AU_", flag.var, ".Rdata")
  load(file_save)

  RPSS_df_all1 <- RPSS_df_all %>% tidyr::spread(Model, RPSS) %>%
    mutate(Mod1 = (Mod1 - Mod0)/(1 - Mod0)) %>%
    mutate(Mod2 = (Mod2 - Mod0)/(1 - Mod0)) %>%
    mutate(Mod3 = (Mod3 - Mod0)/(1 - Mod0)) %>%
    dplyr::select(!Mod0)

  RPSS_df_all1$CI <- factor(RPSS_df_all1$CI, levels = CI_names)

  RPSS_df <- RPSS_df_all1 %>%
    dplyr::group_by(CI, Lead, Station) %>%
    summarise(Mod1 = mean(Mod1),
              Mod2 = mean(Mod2),
              Mod3 = mean(Mod3),
              .groups = "drop")

  RPSS_dt <- RPSS_df %>% data.table()

  filename <- "GRDC_Stations_Australia.csv"
  excel_GRDC <- read.csv(paste0("../data/", filename), header = TRUE,
                         stringsAsFactors = FALSE)
  ind_stn <- sapply(stn_nos, function(x) which(x == excel_GRDC[, 1]))

  loc_stn <- data.frame(Station = stn_nos,
                        excel_GRDC[ind_stn, c("LONG_NEW", "LAT_NEW")])

  world <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")

  labels_pct <- c("C1: 0~5%", "C2: 5~10%", "C3: 10~20%", "C4: >20%")
  rpss_levels <- c("(0,0.05]", "(0.05,0.1]", "(0.1,0.2]", "(0.2,Inf]")
  tol_no <- 219
  legend_override <- list(alpha = 1,shape = 21,
                          stroke = 0.25)

  RPSS_dt_pct1 <- NULL
  map_sf_by_group <- vector("list", 3)

  for(i_mod in 1:3){
    if(i_mod == 1){
      RPSS_dt1 <- RPSS_dt %>%
        subset(Mod2 > Mod1 & Mod1 > 0) %>%
        mutate(Dif = (Mod2 - Mod1))
    } else if(i_mod == 2){
      RPSS_dt1 <- RPSS_dt %>%
        subset(Mod3 > Mod1 & Mod1 > 0) %>%
        mutate(Dif = (Mod3 - Mod1))
    } else {
      RPSS_dt1 <- RPSS_dt %>%
        subset(Mod3 > Mod2 & Mod2 > 0) %>%
        mutate(Dif = (Mod3 - Mod2))
    }

    RPSS_dt_avg <- RPSS_dt1[, .(RPSS = max(Dif)), by = c("Lead", "Station")]

    RPSS_dt_sub <- RPSS_dt_avg %>%
      merge(loc_stn, by = "Station") %>%
      data.table() %>%
      mutate(flag = 1) %>%
      mutate(RPSS_cut = factor(cut(RPSS, breaks = c(0, 0.05, 0.1, 0.2, Inf)),
                               levels = rpss_levels))
    
    #RPSS_dt_sub %>% subset(Lead==3) %>% summary() %>% print()

    loc_stn_sf <- sf::st_as_sf(x = RPSS_dt_sub,
                               coords = c("LONG_NEW", "LAT_NEW"),
                               crs = sf::st_crs(4326))
    loc_stn_sf$Lead <- as.integer(as.character(RPSS_dt_sub$Lead))

    map_sf_by_group[[i_mod]] <- loc_stn_sf

    # Ensure all Lead x RPSS_cut combinations are present, fill missing with pct=0
    all_leads <- sort(unique(RPSS_dt_sub$Lead))
    all_cuts <- rpss_levels
    all_combos <- expand.grid(Lead = all_leads, RPSS_cut = all_cuts, stringsAsFactors = FALSE)
    pct_by_cut <- RPSS_dt_sub[, .(pct = sum(flag, na.rm = TRUE)/tol_no), by = .(Lead, RPSS_cut)]
    pct_by_cut <- merge(all_combos, pct_by_cut, by = c("Lead", "RPSS_cut"), all.x = TRUE)
    pct_by_cut$pct[is.na(pct_by_cut$pct)] <- 0
    pct_by_cut$Group <- i_mod
    # Ensure RPSS_cut is a factor with all levels, so C4 is always present in legend
    pct_by_cut$RPSS_cut <- factor(pct_by_cut$RPSS_cut, levels = all_cuts)
    RPSS_dt_pct1 <- rbind(RPSS_dt_pct1, pct_by_cut)
  }
  
  #RPSS_dt_pct1 %>% spread(Lead, pct) %>% print()
  #RPSS_dt_pct1 %>% summary() %>% print()


  target_map_sf <- map_sf_by_group[[target_group]]
  target_pct <- RPSS_dt_pct1 %>%
    dplyr::filter(Group == target_group) %>%
    mutate(Lead = as.integer(as.character(Lead)))
  

  lead_levels <- sort(unique(target_map_sf$Lead))
  lead_labels_all <- paste0("Lead ", seq_along(lead_levels) - 1)
  lead_labeller <- stats::setNames(lead_labels_all, lead_levels)
  panel_letters <- rep(panel_letters, length.out = length(lead_levels))
  lead_indices <- seq_along(lead_levels) - 1
  panel_titles <- sprintf("%s) Lead %s (%s)", panel_letters, lead_indices, condition_title)

  target_map_sf$RPSS_cut <- factor(target_map_sf$RPSS_cut, levels = rpss_levels)
  target_map_sf$Lead <- factor(target_map_sf$Lead, levels = lead_levels)
  target_pct$Lead <- factor(target_pct$Lead, levels = lead_levels)
  target_pct$RPSS_cut <- factor(target_pct$RPSS_cut, levels = rpss_levels)

  bar_y_max <- max(target_pct$pct) * 100
  if(bar_y_max == 0){
    bar_y_max <- 1
  }

  inset_x <- 0.035
  inset_y <- 0.055
  inset_w <- 0.38
  inset_h <- 0.3

  panel_list <- vector("list", length(lead_levels))
  map_list <- vector("list", length(lead_levels))
  bar_list <- vector("list", length(lead_levels))

  for(i in seq_along(lead_levels)) {
    lead_level <- lead_levels[i]
    map_data <- target_map_sf %>% dplyr::filter(Lead == lead_level)
    bar_data <- target_pct %>% dplyr::filter(Lead == lead_level)
    map_data$RPSS_cut <- factor(map_data$RPSS_cut, levels = rpss_levels)
    missing_cuts <- setdiff(rpss_levels, as.character(unique(map_data$RPSS_cut)))

    total_pct <- sum(bar_data$pct) * 100
    title_text <- sprintf("Total Percent: %.1f%%", total_pct)

    #map_data %>% summary() %>% print()
    ## map ----
    map_plot <- ggplot(data = map_data) +
      geom_sf(aes(fill = RPSS_cut, size = RPSS_cut), shape = 21, alpha = 0.8, stroke = 0.1, 
              show.legend = TRUE) +
      geom_sf(data = world, aes(group = continent), color = "black", alpha = 0.2) +
      coord_sf(
        crs = sf::st_crs(3577),
        xlim = c(-2400000, 2550000),
        ylim = c(-6000000, -1000000),
        label_axes = list(
          bottom = "E", top = NA,
          left = "N", right = "N"
        )) +
      facet_wrap(Lead~., labeller = labeller(Lead = lead_labeller)) +
      scale_size_manual(
        values = stats::setNames(c(0.5, 1, 1.5, 2), rpss_levels),
        labels = labels_pct,
        drop = FALSE,
        breaks = rpss_levels,
        limits = rpss_levels
      ) +
      scale_fill_manual(
        values = stats::setNames(col_pair[3:6], rpss_levels),
        labels = labels_pct,
        drop = FALSE,
        breaks = rpss_levels,
        limits = rpss_levels
      ) +
      guides(fill = guide_legend(ncol = 2,
                                 drop = FALSE,
                                 #override.aes = legend_override
                                 ),
             size = guide_legend(ncol = 2,
                                 drop = FALSE,
                                 override.aes = legend_override)) +
      labs(title = NULL,
           x = NULL, y = NULL,
           fill = "Improve in CRPSS",
           size = "Improve in CRPSS") +
      egg::theme_article(base_size = 7) +
      theme(text = element_text(size = 7, family = "sans"),
            plot.margin = unit(c(0.5, 0.2, 0.15, 0.2), "cm"),
            panel.background = element_rect(fill = NA),
            panel.grid.major = element_line(linewidth = 0.05,
                                            color = alpha("grey", 0.8)),
            strip.background = element_rect(fill="#DDDDDD", colour="grey4"),
            strip.placement = "outside",
            strip.text.x = element_text(size = 7,family="sans", face='bold'),

            axis.ticks.length = unit(0.1, "cm"),
            axis.ticks = element_line(linewidth = 0.2),
            axis.text = element_text(size = 7, family = "sans", face = "plain", color = "black"),
            plot.title = element_text(hjust = 0, face = "bold"))

    map_plot <- map_plot +
      theme(legend.position = c(0.7,0.13),
            legend.direction = "vertical",
            legend.title = element_text(size=6, family = "sans", face="bold"),
            legend.text = element_text(size = 5, margin = margin(l = -5)),
            legend.key.spacing.y = unit(-1,"pt"),
            legend.spacing.x = unit(-1,"pt"),
            legend.box.spacing = unit(1, "pt"))


    if(target_group==1) limits.y <-  c(0, 50) else if(target_group==3) limits.y <- c(0,80)
    else if(target_group==2) limits.y <- c(0, 35)

    ## bar ----
    # Ensure all RPSS_cut levels (C1 to C4) are present, even if pct=0 for some
    # Create complete data with all levels
    all_cuts_complete <- rpss_levels
    
    #all_cuts_complete %>% print()
    
    bar_data_complete <- expand.grid(
      RPSS_cut = all_cuts_complete,
      stringsAsFactors = FALSE
    ) %>%
      left_join(bar_data, by = "RPSS_cut") %>%
      mutate(pct = ifelse(is.na(pct), 0, pct)) %>%
      mutate(RPSS_cut = factor(RPSS_cut, levels = rpss_levels))
    
    #bar_data_complete %>% summary() %>% print()
    
    bar_plot <- ggplot(data = bar_data_complete,
                       aes(x = RPSS_cut, y = pct * 100, fill = RPSS_cut)) +
      geom_col(width = 0.6, color = NA) +
      scale_fill_manual(
        values = stats::setNames(col_pair[3:6], rpss_levels),
        labels = labels_pct,
        guide = "none",
        drop = FALSE,
        limits = rpss_levels
      ) +
      scale_x_discrete(
        labels = paste0("C", 1:4),
        drop = FALSE
      ) +
      scale_y_continuous(limits = limits.y,
                         breaks = scales::pretty_breaks(n = 3),
                         expand = c(0, 0)) +
      labs(title = title_text,
           x = NULL, y = "%") +
      egg::theme_article(base_size = 6) +
      theme(plot.background = element_rect(fill = "transparent", color = "transparent"),
            panel.background = element_rect(fill = "transparent"),
            panel.grid.major = element_blank(),
            plot.margin = unit(c(0.15, 0.2, 0.1, 0.2), "cm"),
            
            plot.title = element_text(size = 6, hjust = 0.4, face = "bold"),
            axis.text.x = element_text(size = 5, family = "sans", color = "black", 
                                       angle = 0),
            axis.text.y = element_text(size = 5, family = "sans", color = "black"),
            axis.ticks.length = unit(0.05, "cm"),
            axis.ticks = element_line(linewidth = 0.2))

    composite <- cowplot::ggdraw(map_plot) +
      cowplot::draw_plot(bar_plot,
                         x = inset_x,
                         y = inset_y,
                         width = inset_w*0.8,
                         height = inset_h*1,
                         hjust = -0.3,
                         vjust = -0.1)

    panel_list[[i]] <- composite
    map_list[[i]] <- map_plot
    bar_list[[i]] <- bar_plot
  }

  list(panels = panel_list,
       maps = map_list,
       bars = bar_list,
       bar_data = target_pct,
       panel_titles = panel_titles)
}

panel_letters_list <- list(
  wet = letters[1:3],
  dry = letters[4:6]
)

condition_plots <- lapply(names(panel_letters_list), function(flag.var) {
  create_condition_plot(flag.var = flag.var,
                        panel_letters = panel_letters_list[[flag.var]])
})
names(condition_plots) <- names(panel_letters_list)

wet_column <- cowplot::plot_grid(plotlist = condition_plots[["wet"]]$panels,
                                 ncol = 1,
                                 align = "v",
                                 rel_heights = rep(1, length(condition_plots[["wet"]]$panels)))


dry_column <- cowplot::plot_grid(plotlist = condition_plots[["dry"]]$panels,
                                  ncol = 1,
                                  align = "v",
                                  rel_heights = rep(1, length(condition_plots[["dry"]]$panels)))

figS_map <- cowplot::plot_grid(wet_column,
                           dry_column,
                           ncol = 2,
                           labels = c("a) SSI 7-day maximum","b) SSI 30-day minimum"),
                           hjust=-0.25, vjust=1.1,
                           label_size = 10, 
                           align = "hv",
                           rel_widths = c(1, 1))

figS_map %>% print()

if(flag.save.supp){
  filen_s <- paste0("../figure/FigureS_CRPSS_ALL_", mode, "_", wf, "_", method_fevd, flag.v, 
                    "_AU_SSI_map_group", target_group, ".png")
  
  png(filen_s, height=21.6, width=15, units="cm", bg = "white", res=600)
  figS_map %>% print()
  dev.off()
}
