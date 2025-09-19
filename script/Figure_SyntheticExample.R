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
library(cowplot)

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
n <- 5000

system_specs <- list(
  Duffing = list(
    label = "Duffing",
    prefix = "Duffing",
    sd_label = "0.1",
    n = n,
    generator = function() {
      ts <- WASP2.0::data.gen.Duffing(nobs = n, s = 0.1, do.plot = FALSE)
      data.frame(No = 1:n, x = ts$x, y = ts$y)
    }
  ),
  Rossler = list(
    label = "Rossler",
    prefix = "Rossler",
    sd_label = "0.1",
    n = n,
    generator = function() {
      ts <- WASP2.0::data.gen.Rossler(a = 0.2, b = 0.2, w = 5.7, start = c(-2, -10, 0.2),
                                      n = n, s = 0.1)
      data.frame(No = 1:n, x = ts$x, y = ts$y, z = ts$z)
    }
  ),
  Lorenz = list(
    label = "Lorenz",
    prefix = "Lorenz",
    sd_label = "5",
    n = n,
    generator = function() {
      time <- seq(0, 50, length.out = n)
      ts <- WASP2.0::data.gen.Lorenz(time = time, s = 5)
      data.frame(No = 1:n, x = ts$x, y = ts$y, z = ts$z)
    }
  )
)

font_size <- 8
time_series_plots <- list()
metric_plots <- list()
system_levels <- sapply(system_specs, function(spec) spec$label)

for (sys_name in names(system_specs)) {
  spec <- system_specs[[sys_name]]
  label <- spec$label
  ts_wide <- spec$generator()
  ts_long <- ts_wide %>%
    gather(group, value, -No) %>%
    mutate(system = label)

  p_ts <- ggplot(ts_long, aes(x = No, y = value, color = group)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = col_qual) +
    labs(x = NULL, y = NULL, title = label, color = NULL) +
    egg::theme_article() +
    theme(text = element_text(size = font_size),
          plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "cm"),
          plot.background = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "top",
          legend.title = element_blank(),
          legend.key.width = unit(1, "cm"))

  time_series_plots[[label]] <- p_ts

  file_dat <- file.path("..", "result",
                        sprintf("%s_%s_n%d_sd%s_r%d%s.Rdata",
                                spec$prefix, wf, spec$n, spec$sd_label,
                                ensemble, flag.v))
  stopifnot(file.exists(file_dat))
  load(file_dat)

  pred_dt <- pred_df_list %>%
    gather(model, pred, which(colnames(pred_df_list) %in% paste0("mod", 1:3))) %>%
    mutate(system = label) %>%
    mutate_if(is.character, as.factor) %>%
    data.table()

  pred_dt$system <- factor(pred_dt$system, levels = system_levels)

  metric_df <- pred_dt[, .(RMSE = Metrics::rmse(obs, pred),
                           R = cor(obs, pred)),
                       by = c("system", "method", "group", "model", "r")]

  metric_df$group <- factor(as.character(metric_df$group),
                            levels = c("cal", "val"),
                            labels = c("Calibration", "Validation"))
  metric_df$model <- factor(as.character(metric_df$model),
                            levels = paste0("mod", 1:3),
                            labels = c("x", "x'", "x'(φ)"))
  metric_df$method <- factor(as.character(metric_df$method))

  metric_df1 <- metric_df %>%
    gather(metric, value, R, RMSE)

  metric_df1$metric <- factor(metric_df1$metric, levels = c("R", "RMSE"))

  metric_subset <- metric_df1 %>%
    filter(system == label, group == "Validation", method == "lm") %>%
    droplevels()

  p_metric <- ggplot(metric_subset, aes(x = model, y = value, fill = model)) +
    geom_boxplot(linewidth = 0.3, outlier.size = 0.8) +
    facet_wrap(metric ~ ., scales = "free") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_fill_manual(values = col_qual) +
    labs(x = NULL, y = NULL, title = label) +
    egg::theme_article() +
    theme(text = element_text(family = "sans", face = "bold", size = font_size),
          plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "cm"),
          plot.background = element_blank(),
          panel.background = element_blank(),
          strip.background.x = element_rect(fill = "#DEDEDE", colour = "black"),
          strip.background.y = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(size = font_size, family = "sans", face = "plain"),
          axis.text = element_text(size = font_size, family = "sans", face = "plain", color = "black"),
          axis.text.x = element_text(size = font_size, family = "sans", face = "bold", color = "black"),
          legend.position = "none")

  metric_plots[[label]] <- p_metric
  rm(pred_df_list)
}

figure_grobs <- c(time_series_plots[system_levels], metric_plots[system_levels])
fig <- cowplot::plot_grid(plotlist = figure_grobs, ncol = 3, align = "hv")
fig %>% print()

if(flag.save){
  filen <- paste0("Figure_SyntheticExample_", wf, "_n", n, "_0813.png")
  png(filen, height = 14, width = 18, units = "cm", bg = "white", res = 600)
  fig %>% print()
  dev.off()
}
