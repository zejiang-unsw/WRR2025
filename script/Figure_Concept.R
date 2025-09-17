# Experiment on phases
rm(list = ls()) # remove all variables
#graphics.off() # remove all figures

# Getting the path of your current open file
current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))

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
col_qual <- brewer.pal(n = col_len, name = "Dark2")
col_div <- brewer.pal(n = col_len, name = "RdYlBu")
col_seq <- brewer.pal(n = col_len, name = "Blues")

font_size <- 10

# synthetic data ----
set.seed(101)
n <- nobs <- 512/2
sd.y <- sd.x <- 0
fs <- 1 / nobs
t <- seq(0, 1, length.out = nobs)
phi1 <- pi / 2

freq_x <- 12 # common predictor frequency
x_common <- sin(2 * pi * freq_x * t + phi1) + rnorm(nobs, 0, sd.x)

wf <- "d6"

l2_norm <- function(x) {
  sqrt(sum(x^2))
}

build_case <- function(freq_y) {
  y <- sin(2 * pi * freq_y * t)
  data <- list(x = y, dp = as.matrix(x_common))
  dp.n <- phase.wasp.index(data$dp, data$x, index = 1:n, k = n/2-1)

  dwt0 <- WASP2.0::dwt.vt(list(x = y, dp = y), wf, verbose = FALSE)
  dwt <- WASP2.0::dwt.vt(data, wf, verbose = FALSE)
  dwt1 <- WASP2.0::dwt.wasp.val(data, wf, phase = TRUE, index_c = 1:n, 
                                k = n/2-1, max_iter = 3, 
                                #k_auto = TRUE, power_thresh = 0.8, k_max = n/2-1,
                                #momentum = 0.3,
                                verbose = FALSE)

  list(
    y = y,
    dp.n = dp.n,
    dwt0 = dwt0,
    dwt = dwt,
    dwt1 = dwt1
  )
}

freq_high <- 24
freq_low <- 6

case_high <- build_case(freq_high)
case_low <- build_case(freq_low)

# convenience themes ----
line_theme <- function(base_size = font_size) {
  theme(text = element_text(size = base_size),
        plot.margin = unit(c(0.4, 0.2, 0.2, 0.2), "cm"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
}

cov_theme <- theme(text = element_text(size = font_size),
                   plot.margin = unit(c(0.6, 0.3, 0.3, 0.3), "cm"),
                   panel.grid = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank(),
                   plot.background = element_blank(),
                   axis.title = element_blank(),
                   axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   legend.position = "none")

scatter_theme <- theme(text = element_text(size = font_size),
                       plot.margin = unit(c(-0.5, 0.2, 0, 0.2), "cm"),
                       panel.grid = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank(),
                       plot.background = element_blank(),
                       axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       legend.position = "none")

# Fig-a ----
df_y_high <- tibble::tibble(No = 1:n, value = case_high$y)
df_y_low <- tibble::tibble(No = 1:n, value = case_low$y)
df_x_common <- tibble::tibble(No = 1:n, value = x_common)

p1a_y_high <- ggplot(df_y_high, aes(No, value)) +
  geom_line(color = "black", linewidth = 0.6) +
  scale_x_continuous(limits = c(0, 100)) +
  labs(title = "High-frequency y (Flood)") +
  line_theme()

p1a_y_low <- ggplot(df_y_low, aes(No, value)) +
  geom_line(color = "black", linewidth = 0.6) +
  scale_x_continuous(limits = c(0, 100)) +
  labs(title = "Low-frequency y (Drought)") +
  line_theme()

p1a_x <- ggplot(df_x_common, aes(No, value)) +
  geom_line(color = col_qual[2], linewidth = 0.6) +
  scale_x_continuous(limits = c(0, 100)) +
  labs(title = "Predictor x") +
  line_theme()

p1a_high <- cowplot::plot_grid(p1a_y_high, p1a_x, ncol = 1, align = "v")
p1a_low <- cowplot::plot_grid(p1a_y_low, p1a_x, ncol = 1, align = "v")


figa <- cowplot::plot_grid(p1a_high, p1a_low, ncol=2)
figa

# Fig-b ----
df_phase_high <- tibble::tibble(
  No = 1:n,
  y = case_high$y,
  x_phase = case_high$dp.n[, 1],
  x_raw = x_common
)

df_phase_low <- tibble::tibble(
  No = 1:n,
  y = case_low$y,
  x_phase = case_low$dp.n[, 1],
  x_raw = x_common
)

p1b_high_y <- ggplot(df_phase_high) +
  geom_line(aes(No, y), color = "black", linewidth = 0.6) +
  #geom_line(aes(No, x_phase), color = col_qual[3], linewidth = 0.6, alpha = 0.8) +
  scale_x_continuous(limits = c(0, 100)) +
  labs(title = "High-frequency y") +
  line_theme()

p1b_low_y <- ggplot(df_phase_low) +
  geom_line(aes(No, y), color = "black", linewidth = 0.6) +
  #geom_line(aes(No, x_phase), color = col_qual[3], linewidth = 0.6, alpha = 0.8) +
  scale_x_continuous(limits = c(0, 100)) +
  labs(title = "Low-frequency y") +
  line_theme()

p1b_high_shift <- df_phase_high %>%
  tidyr::pivot_longer(cols = c("x_raw", "x_phase"), names_to = "Series", values_to = "value") %>%
  ggplot(aes(No, value, colour = Series, linetype = Series)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c(x_raw = col_qual[2], x_phase = col_qual[3])) +
  scale_linetype_manual(values = c(x_raw = "dashed", x_phase = "solid")) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(-1.3, 1.3)) +
  labs(title = "Phase-aligned x") +
  line_theme()

p1b_low_shift <- df_phase_low %>%
  tidyr::pivot_longer(cols = c("x_raw", "x_phase"), names_to = "Series", values_to = "value") %>%
  ggplot(aes(No, value, colour = Series, linetype = Series)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c(x_raw = col_qual[2], x_phase = col_qual[3])) +
  scale_linetype_manual(values = c(x_raw = "dashed", x_phase = "solid")) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(-1.3, 1.3)) +
  labs(title = "Phase-aligned x") +
  line_theme()

p1b_high <- cowplot::plot_grid(p1b_high_y, p1b_high_shift, ncol = 1, align = "v")
p1b_low <- cowplot::plot_grid(p1b_low_y, p1b_low_shift, ncol = 1, align = "v")


figb <- cowplot::plot_grid(p1b_high_shift, p1b_low_shift, ncol=2)
figb

# Fig-c ----
df_cov_high <- data.frame(No = 1:nrow(case_high$dwt$S),
                          S0 = case_high$dwt0$S / l2_norm(case_high$dwt0$S),
                          S1 = case_high$dwt$S[, 1] / l2_norm(case_high$dwt$S[, 1]),
                          S2 = case_high$dwt1$S[, 1] / l2_norm(case_high$dwt1$S[, 1])) %>%
  tidyr::pivot_longer(cols = c("S1", "S2"), names_to = "Series", values_to = "S") %>%
  mutate(S = round(S, 3), S0 = round(S0, 3))

df_cov_low <- data.frame(No = 1:nrow(case_low$dwt$S),
                         S0 = case_low$dwt0$S / l2_norm(case_low$dwt0$S),
                         S1 = case_low$dwt$S[, 1] / l2_norm(case_low$dwt$S[, 1]),
                         S2 = case_low$dwt1$S[, 1] / l2_norm(case_low$dwt1$S[, 1])) %>%
  tidyr::pivot_longer(cols = c("S1", "S2"), names_to = "Series", values_to = "S") %>%
  mutate(S = round(S, 3), S0 = round(S0, 3))

p1c_cov_high <- ggplot(df_cov_high, aes(No, S, fill = Series)) +
  geom_col(position = "dodge") +
  geom_line(aes(No, S0), color = "black", linewidth = 0.6) +
  geom_point(aes(No, S0), color = "black", size = 1.5) +
  scale_fill_manual(values = col_qual[c(1, 3)]) +
  labs(title = "Variance-modulated x") +
  cov_theme

p1c_cov_low <- ggplot(df_cov_low, aes(No, S, fill = Series)) +
  geom_col(position = "dodge") +
  geom_line(aes(No, S0), color = "black", linewidth = 0.6) +
  geom_point(aes(No, S0), color = "black", size = 1.5) +
  scale_fill_manual(values = col_qual[c(1, 3)]) +
  labs(title = "Variance-modulated x") +
  cov_theme

df_x_trans_high <- tibble::tibble(No = 1:n,
                                  x_raw = case_high$dwt$dp[, 1],
                                  x_trans = case_high$dwt1$dp.n[, 1]) %>%
  tidyr::pivot_longer(cols = c("x_raw", "x_trans"), names_to = "Series", values_to = "value")

df_x_trans_low <- tibble::tibble(No = 1:n,
                                 x_raw = case_low$dwt$dp[, 1],
                                 x_trans = case_low$dwt1$dp.n[, 1]) %>%
  tidyr::pivot_longer(cols = c("x_raw", "x_trans"), names_to = "Series", values_to = "value")

p1c_x_high <- ggplot(df_x_trans_high, aes(No, value, colour = Series, linetype = Series)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c(x_raw = col_qual[2], x_trans = col_qual[4])) +
  scale_linetype_manual(values = c(x_raw = "dashed", x_trans = "solid")) +
  scale_x_continuous(limits = c(0, 100)) +
  labs(title = "Transformed x") +
  line_theme()

p1c_x_low <- ggplot(df_x_trans_low, aes(No, value, colour = Series, linetype = Series)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c(x_raw = col_qual[2], x_trans = col_qual[4])) +
  scale_linetype_manual(values = c(x_raw = "dashed", x_trans = "solid")) +
  scale_x_continuous(limits = c(0, 100)) +
  labs(title = NULL) +
  line_theme()

p1c_high <- cowplot::plot_grid(p1c_cov_high, p1c_x_high, ncol = 1, align = "v")
p1c_low <- cowplot::plot_grid(p1c_cov_low, p1c_x_low, ncol = 1, align = "v")

figc <- cowplot::plot_grid(p1c_high, p1c_low, ncol=2)
figc

# Fig-d ----
fit_raw_high <- stats::lm(case_high$dwt$x ~ ., data = case_high$dwt$dp[, 1] %>% data.frame())
fit_trans_high <- stats::lm(case_high$dwt1$x ~ ., data = case_high$dwt1$dp.n[, 1] %>% data.frame())

scatter_high <- tibble::tibble(
  obs = case_high$dwt$x,
  `Raw predictor` = fit_raw_high$fitted.values,
  `Transformed predictor` = fit_trans_high$fitted.values
) %>%
  tidyr::pivot_longer(cols = c("Raw predictor", "Transformed predictor"),
                      names_to = "Series", values_to = "pred")

fit_raw_low <- stats::lm(case_low$dwt$x ~ ., data = case_low$dwt$dp[, 1] %>% data.frame())
fit_trans_low <- stats::lm(case_low$dwt1$x ~ ., data = case_low$dwt1$dp.n[, 1] %>% data.frame())

scatter_low <- tibble::tibble(
  obs = case_low$dwt$x,
  `Raw predictor` = fit_raw_low$fitted.values,
  `Transformed predictor` = fit_trans_low$fitted.values
) %>%
  tidyr::pivot_longer(cols = c("Raw predictor", "Transformed predictor"),
                      names_to = "Series", values_to = "pred")

plot_scatter <- function(df, title_text = NULL) {
  ggplot(df, aes(pred, obs)) +
    geom_point(color = col_qual[2], fill = scales::alpha(col_qual[2], 0.4), shape = 21, size = 2) +
    geom_abline(intercept = 0, slope = 1, color = "grey40") +
    labs(title = title_text, x = "Predicted", y = "Observed") +
    coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
    scatter_theme
}

plot_scatter_trans <- function(df, title_text = NULL) {
  ggplot(df, aes(pred, obs)) +
    geom_point(color = col_qual[4], fill = scales::alpha(col_qual[4], 0.4), shape = 21, size = 2) +
    geom_abline(intercept = 0, slope = 1, color = "grey40") +
    labs(title = title_text, x = "Predicted", y = "Observed") +
    coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
    scatter_theme
}

p1d_high <- cowplot::plot_grid(
  plot_scatter(scatter_high %>% filter(Series == "Raw predictor"), ""),
  plot_scatter_trans(scatter_high %>% filter(Series == "Transformed predictor"), ""),
  ncol = 2,
  align = "hv"
)

p1d_low <- cowplot::plot_grid(
  plot_scatter(scatter_low %>% filter(Series == "Raw predictor"), ""),
  plot_scatter_trans(scatter_low %>% filter(Series == "Transformed predictor"), ""),
  ncol = 2,
  align = "hv"
)

figd <- cowplot::plot_grid(p1d_high, p1d_low, ncol=2)
figd

# Final figure ----
fig <- cowplot::plot_grid(figa, figb,
                          figc, figd,
                          rel_heights = c(1,0.5,1.2,1),
                          ncol = 1,
                          label_fontfamily = "sans",
                          label_fontface = "bold",
                          label_colour = "black")

fig %>% print()

# Figure export template ----
if (FALSE) {
  filen <- "Figure_concept"
  export::graph2pdf(x = fig, file = filen, aspectr = 2, font = "Arial",
                    height = 16 / 2.54, width = 12 / 2.54, bg = "transparent")

  png(filen, height = 16, width = 12, units = "cm", bg = "transparent", res = 600)
  fig %>% print()
  dev.off()
}
