# Experiment on phases
rm(list=ls()) # remove all variables
graphics.off() # remove all figures

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))

flag.save <- T

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
library(rgl)

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

aux_dir <- file.path(".","insets")
if(!dir.exists(aux_dir)) dir.create(aux_dir, recursive = TRUE)

system_specs <- list(
  Duffing = list(
    label = "Duffing",
    prefix = "Duffing",
    sd_label = "0.1",
    n = n,
    var_labels = c(x = "x", y = "y"),
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
    var_labels = c(x = "x1", y = "x2", z = "y"),
    generator = function() {
      ts <- WASP2.0::data.gen.Rossler(a = 0.2, b = 0.2, w = 5.7, start = c(-2, -10, 0.2),
                                      n = n, s = 0.1)
      data.frame(No = 1:n, x = ts$x, y = ts$y, z = ts$z)
    }
  ),
  Lorenz = list(
    label = "Lorenz",
    prefix = "Lorenz",
    sd_label = "1",
    n = n,
    var_labels = c(x = "x1", y = "x2", z = "y"),
    generator = function() {
      time <- seq(0, 50, length.out = n)
      ts <- WASP2.0::data.gen.Lorenz(time = time, s = 1)
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
    mutate(group = factor(group,
                          levels = names(spec$var_labels),
                          labels = spec$var_labels),
           system = label)

  x_limit <- if (sys_name == "Duffing") c(0, 100) else c(0, 500)

  base_theme <- egg::theme_article() +
    theme(text = element_text(size = font_size),
          plot.margin = grid::unit(c(0.35, 0.1, 0.1, 0.5), "cm"),
          plot.background = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          legend.key = element_blank(),
          legend.key.width = grid::unit(0.6, "cm"),
          axis.title = element_text(size = font_size),
          axis.text = element_text(size = font_size))

  p_ts_base <- ggplot(ts_long, aes(x = No, y = value, color = group)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = col_qual) +
    scale_x_continuous(limits = x_limit, breaks = scales::pretty_breaks(n = 6)) +
    labs(x = "Time step", y = NULL, title = NULL, color = NULL) +
    base_theme +
    theme(legend.position = c(0.8,0.8))
  p_ts_base
  if (label == "Lorenz") {
    p_ts_base <- p_ts_base +
      scale_y_continuous(limits = c(-30, 80), breaks = scales::pretty_breaks(n = 6))
  }

  legend_plot <- ggplot(ts_long, aes(x = No, y = value, color = group)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = col_qual) +
    scale_x_continuous(limits = x_limit) +
    labs(color = NULL) +
    base_theme +
    theme(legend.position = "top",
          legend.justification = "left",
          legend.box.just = "left",
          legend.background = element_rect(fill = scales::alpha("white", 0.6), colour = NA),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())

  legend_grob <- cowplot::get_legend(legend_plot)

  if (label == "Duffing") {
    ts <- WASP2.0::data.gen.Duffing(nobs = n, s=0, do.plot=F)
    ts_wide <- data.frame(No = 1:n, x = ts$x, y = ts$y)
    scatter_df <- ts_wide %>% select(x, y)
    scatter_plot <- ggplot(scatter_df, aes(x = x, y = y)) +
      geom_point(alpha = 1, size = 0.5, color = "black") +
      coord_equal() +
      theme_void()
    export::graph2pdf(x = scatter_plot,
                      file = file.path(aux_dir, paste0(label, "_phase.pdf")),
                      aspectr = 2,
                      font = "Arial",
                      height = 4 / 2.54,
                      width = 4 / 2.54,
                      bg = "transparent")
  } else if (label %in% c("Rossler", "Lorenz")) {
    if(label=="Rossler"){
      scatter_df <- WASP2.0::data.gen.Rossler(a = 0.2, b = 0.2, w = 5.7, 
                                        start = c(-2, -10, 0.2), n, s=0)
    } else {
      scatter_df <- WASP2.0::data.gen.Lorenz(sigma = 10, beta = 8/3, rho = 28, start = c(-13, -14, 47),
                                     time = seq(0, by=0.05, length.out = n), s=0)
    }
    
    pdf_file <- file.path(aux_dir, paste0(label, "_phase.pdf"))
    rgl::bg3d(color = "transparent")
    rgl::par3d(windowRect = c(0, 0, 500, 500))
    fig <- rgl::plot3d(scatter_df$x,
                scatter_df$y,
                scatter_df$z,
                type = "l",
                size = 3,
                col = "black", #scales::alpha(col_qual[2], 0.6),
                box = FALSE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                zlab = "")
    fig %>% print()
    rgl.postscript(pdf_file, fmt="pdf")

  }

  p_ts <- cowplot::ggdraw(p_ts_base) +
    cowplot::draw_grob(legend_grob, x = 0.98, y = 0.98, width = 0.32,
                       height = 0.22, hjust = 1, vjust = 1)

  time_series_plots[[label]] <- p_ts

  file_dat <- file.path("..", "result",
                        sprintf("%s_%s_n%d_sd%s_r%d%s.Rdata",
                                spec$prefix, wf, spec$n, spec$sd_label,
                                ensemble, flag.v))
  stopifnot(file.exists(file_dat))
  load(file_dat)

  pred_dt <- pred_df_list %>%
    gather(model, pred, mod1, mod3) %>%
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
                            levels = c("mod1", "mod3"),
                            labels = c("x", "x'(φ)"))
  metric_df$method <- factor(as.character(metric_df$method))

  metric_df1 <- metric_df %>%
    gather(metric, value, R, RMSE)

  metric_df1$metric <- factor(metric_df1$metric, levels = c("R", "RMSE"))

  metric_subset <- metric_df1 %>%
    filter(system == label, group == "Validation", method == "lm") %>%
    droplevels() %>%
    group_by(metric) %>%
    mutate(q_low = quantile(value, 0.05),
           q_high = quantile(value, 0.95),
           pad = (q_high - q_low) * 0.25 + 1e-6,
           plot_value = pmin(pmax(value, q_low - pad), q_high + pad)) %>%
    ungroup()
  ## metric ----
  p_metric <- ggplot(metric_subset, aes(x = model, y = plot_value, fill = model)) +
    geom_violin(trim = FALSE, linewidth = 0.3, colour = NA) +
    stat_summary(fun = mean, geom = "point", shape = 21, size = 1.3, fill = "red") +
    facet_wrap(metric ~ ., scales = "free") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
    scale_fill_manual(values = col_qual) +
    labs(x = NULL, y = NULL, title = NULL) +
    egg::theme_article() +
    theme(text = element_text(family = "sans", face = "bold", size = font_size),
          plot.margin = grid::unit(c(0.35, 0.1, 0.1, 0.5), "cm"),
          plot.background = element_blank(),
          panel.background = element_blank(),
          strip.background.x = element_rect(fill = "#DEDEDE", colour = "black"),
          strip.background.y = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(size = font_size, family = "sans", face = "plain"),
          axis.text = element_text(size = font_size, family = "sans", face = "plain", color = "black"),
          axis.text.x = element_text(size = font_size, family = "sans", face = "bold", color = "black"),
          legend.position = "none")
  p_metric
  metric_plots[[label]] <- p_metric
  rm(pred_df_list)
}

figure_grobs <- c(time_series_plots[system_levels], metric_plots[system_levels])
fig <- egg::ggarrange(plots = figure_grobs,
                      ncol = 3,
                      labels = letters[seq_along(figure_grobs)],
                      label.args = list(gp = grid::gpar(fontsize = 10, fontface = "bold"),
                                        hjust = -0.1,
                                        vjust = 1.2))
fig <- cowplot::plot_grid(plotlist = figure_grobs, ncol=3,
                          labels = letters[seq_along(figure_grobs)],
                          label_size = 10)
fig %>% print()

if(flag.save){
  filen <- paste0("Figure_SyntheticExample_", wf, "_n", n, ".png")
  
  export::graph2pdf(x = fig, file = filen, aspectr = 2, font = "Arial",
                    height = 12 / 2.54, width = 18 / 2.54, bg = "white")
  
  # png(filen, height = 16, width = 18, units = "cm", bg = "white", res = 600)
  # fig %>% print()
  # dev.off()
}
