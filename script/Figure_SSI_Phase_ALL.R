rm(list = ls())

current_path <- ""
# Align working directory with active document when available
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  current_path <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                           error = function(e) "")
  if (nzchar(current_path)) {
    setwd(dirname(current_path))
  }
}

flag.save <- TRUE

# Packages --------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(cowplot)
library(ggpubr)

# Shared styling ---------------------------------------------------------------
col_qual <- brewer.pal(n = 8, name = "Paired")
CI_names <- c("Nino3", "Nino4", "Nino3.4", "EMI", "DMI", "EPI", "WPI", "II",
              "TSI", "TBV", "Disc")
CI_names1 <- c("Niño 3", "Niño 4", "Niño 3.4", "EMI", "DMI", "EPI", "WPI", "II",
               "TSI", "TBV", "Disc")

mode <- "DWT"
wf <- "d16"

case_lookup <- data.frame(
  flag = c("wet", "dry"),
  label = c("SSI 7-day maximum", "SSI 30-day minimum"),
  stringsAsFactors = FALSE
)
case_palette <- stats::setNames(col_qual[c(1, 2)], case_lookup$label)

load_phase_case <- function(case_flag, case_label) {
  file_name <- sprintf("../result/SSI_Phase_%s_%s_ALL_AU_%s.Rdata", mode, wf, case_flag)
  if (!file.exists(file_name)) {
    stop(sprintf("Result file not found: %s", file_name))
  }
  loaded_objs <- load(file_name)
  if (!"phase_df_all" %in% loaded_objs) {
    stop(sprintf("Object 'phase_df_all' missing in %s", file_name))
  }
  df <- phase_df_all
  df$case <- case_label
  df$fold <- factor(df$fold, levels = c(1, 2), labels = c("Fold 1", "Fold 2"))
  df$lead <- factor(df$lead, labels = paste0("Lead ", 0:2))
  df$cpy <- factor(df$cpy, levels = CI_names, labels = CI_names1)
  df
}

phase_list <- vector("list", nrow(case_lookup))
for (i in seq_len(nrow(case_lookup))) {
  phase_list[[i]] <- load_phase_case(case_lookup$flag[i], case_lookup$label[i])
}
phase_df <- dplyr::bind_rows(phase_list)
rm(phase_list)

phase_df$case <- factor(phase_df$case, levels = case_lookup$label)

phase_df_wide <- phase_df %>%
  select(-freq) %>%
  pivot_wider(names_from = fold, values_from = dif) %>%
  rename(fold1 = `Fold 1`, fold2 = `Fold 2`)

disc_df <- phase_df_wide %>%
  filter(cpy == "Disc") %>%
  mutate(phase_gap = abs(fold1 - fold2))

case_levels <- c("Zero reference", case_lookup$label)
case_palette_full <- c("Zero reference" = "grey70", case_palette)[case_levels]

box_data <- dplyr::bind_rows(
  disc_df %>% transmute(case = "Zero reference", lead, phase_gap = 0),
  disc_df %>% transmute(case, lead, phase_gap)
) %>%
  mutate(case = factor(case, levels = case_levels))

scatter_disc <- ggplot(disc_df, aes(x = fold1, y = fold2, fill = lead)) +
  geom_point(size = 1.6, shape = 21, colour = "grey35", alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, colour = "firebrick", linewidth = 0.3) +
  coord_fixed() +
  facet_grid(lead ~case, switch="y") +
  scale_fill_manual(values = col_qual, name = "Lead") +
  scale_x_continuous(limits = c(-pi, pi),
                     breaks = scales::pretty_breaks(n = 4),
                     expand = c(0.1, 0)) +
  scale_y_continuous(limits = c(-pi, pi),
                     breaks = scales::pretty_breaks(n = 4),
                     expand = c(0.1, 0)) +
  labs( #title = "Phase comparison for Disc covariate",
       x = expression(Delta * phi["Fold 1"]),
       y = expression(Delta * phi["Fold 2"])) +
  #theme_minimal(base_size = 7, base_family = "sans") +
  egg::theme_article(base_size = 10, base_family = "sans") + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    strip.background = element_rect(fill = "#DDDDDD", colour = "grey40"),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(colour = "black", size = 8),
    axis.title = element_text(size = 8),
    legend.position = "none",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    plot.margin = margin(4, 6, 4, 4)
  )

scatter_disc

max_phase_gap <- max(box_data$phase_gap, na.rm = TRUE)

ggpubr::compare_means(phase_gap ~ case, group.by = "lead", #method="t.test",
                      ref.group = "Zero reference",
                      data = box_data) %>%
  print()


ggpubr::compare_means(phase_gap ~ case, group.by = "lead", method="t.test",
                      #ref.group = "Zero reference",
                      data = box_data %>% subset(case!="Zero reference")) %>%
  print()

summary(box_data)

box_disc <- ggplot(box_data %>% subset(case!="Zero reference"), aes(x = case, y = phase_gap, fill = case)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.7, linewidth = 0.3, outliers =F) +
  facet_wrap(~lead, nrow = 1) +
  scale_fill_manual(values = case_palette_full, guide = "none") +
  labs(x = NULL,
       y = expression("|" * Delta * phi["Fold 1"] - Delta * phi["Fold 2"] * "| (rad)")) +
  ggpubr::stat_compare_means(#ref.group = "Zero reference",
                             label = "p.signif",
                             method = "wilcox.test",
                             size = 3,
                             label.y = max_phase_gap + 0.1) +
  geom_hline(yintercept = 0, colour = "grey20", linewidth = 0.2) +
  theme_minimal(base_size = 7, base_family = "sans") +
  theme(
    axis.text.x = element_text(colour = "black", size = 6, face = "bold", angle = 15, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 6),
    axis.title = element_text(size = 7),
    strip.background = element_rect(fill = "#DDDDDD", colour = "grey40"),
    strip.text = element_text(face = "bold"),
    plot.margin = margin(4, 6, 4, 4)
  )

box_disc

figure_s3 <- cowplot::plot_grid(
  scatter_disc,
  box_disc,
  labels = c("a", "b"),
  label_size = 8,
  label_fontfamily = "sans",
  label_fontface = "bold",
  ncol = 1,
  rel_heights = c(1.2, 1)
)

figure_s3 %>% print()

if (flag.save) {
  dev.size(units = c("in", "cm", "px")[2])
  dev.size(units = c("in", "cm", "px")[1])
  
  filen_s = "../figure/Figure_S3_SSI_Phase_Disc.png"
  
  export::graph2pdf(x = scatter_disc, file = filen_s, aspectr = 2, font = "Arial",
                    height = 5, width = 4, 
                    #height = 12 / 2.54, width = 9 / 2.54, 
                    bg = "white")
  
  png(filen_s, height=5, width=4, units="in", bg = "white", res=600)
  scatter_disc %>% print()
  dev.off()
  
}
