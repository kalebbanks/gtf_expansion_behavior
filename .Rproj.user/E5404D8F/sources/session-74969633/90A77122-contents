# Script: 03_figures.R
# Figures for "Behavioral divergence during rapid range expansion of a native amphibian"
# Purpose: Create Figures 2, 4, 5 and Figure S1
# Inputs: data/emergence_long.csv, data/movement_long.csv, data/acoustic_results.csv,
#         data/agg.wav, data/adv.wav, data/playback.wav
# Outputs: outputs/figures/
# Author: XXXXXXXXXXXXXX
# Date: 10/12/2026

#### Packages ####
#install.packages(c("ggplot2", "tidyr", "dplyr", "ggpubr", "tuneR", "seewave", "survival"))
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(tuneR)
library(seewave)
library(survival)

#### Data ####
emergence <- read.csv("data/emergence_long.csv")   # latency to exit shelter (max 1200 s); event = 0 if censored
movement  <- read.csv("data/movement_long.csv")    # arena movement measures (max 1800 s)
acoustics <- read.csv("data/acoustic_results.csv")

agg      <- readWave("data/agg.wav")
adv      <- readWave("data/adv.wav")
playback <- readWave("data/playback.wav")


dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

# Figure 1 was created in QGIS.
# Figure 3 was created in Adobe Illustrator.


#### Figure 2: Waveforms of aggressive call, advertisement call and playback stimulus ####
# Edited in Adobe Illustrator for presentation
par(mfrow = c(3, 1), mar = c(4, 4, 2, 2))
oscillo(agg,      colwave = "black", labels = FALSE, title = "Aggressive call")
oscillo(adv,      colwave = "black", labels = FALSE, title = "Advertisement call")
oscillo(playback, colwave = "black", alab = "",      title = "Playback stimulus")


#### Figure 4: Range comparisons of repeatable behaviors ####

# Kaplan-Meier curves: cumulative proportion that exited/moved, censored trials shown as crosses
km_plot <- function(dat, time_col, xlab, max_time) {
  fit <- survfit(Surv(dat[[time_col]], dat$event) ~ Range, data = dat)
  km <- data.frame(
    time       = fit$time,
    prop_event = 1 - fit$surv,
    n.censor   = fit$n.censor,
    Range      = sub("Range=", "", rep(names(fit$strata), fit$strata))
  )
  # start each curve at time 0
  km <- rbind(data.frame(time = 0, prop_event = 0, n.censor = 0,
                         Range = unique(km$Range)), km)
  km <- km[order(km$Range, km$time), ]

  ggplot(km, aes(x = time, y = prop_event, color = Range, linetype = Range)) +
    geom_step(linewidth = 1.2) +
    geom_point(data = subset(km, n.censor > 0), shape = 3, size = 4, stroke = 1.2,
               show.legend = FALSE) +
    scale_color_manual(values = c("Expanded" = "orange", "Historical" = "grey40")) +
    scale_linetype_manual(values = c("Expanded" = "solid", "Historical" = "dashed")) +
    scale_x_continuous(limits = c(0, max_time * 1.04),
                       breaks = seq(0, max_time, by = max_time / 4),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(x = xlab, y = "Cumulative proportion", color = NULL, linetype = NULL) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = c(0.97, 0.05),
      legend.justification = c(1, 0),
      legend.background = element_blank(),
      legend.key.width = unit(1.5, "cm"),
      legend.text = element_text(size = 13),
      axis.title = element_text(face = "bold", size = 15),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.margin = margin(10, 15, 10, 10)
    )
}

## Figure 4A: Latency to exit shelter
fig4a <- km_plot(emergence, "latency", "Latency to exit shelter (s)", 1200) +
  annotate("text", x = 1200 * 0.97, y = 0.30, hjust = 1,
           label = "italic(P) == 0.031", parse = TRUE, size = 5)

## Figure 4B: Latency to first movement
latency_move <- movement %>%
  filter(measure == "Latency", !is.na(value)) %>%
  mutate(event = ifelse(value == 1800, 0, 1))     # 0 = did not move within 1800 s

fig4b <- km_plot(latency_move, "value", "Latency to first movement (s)", 1800) +
  annotate("text", x = 1800 * 0.97, y = 0.30, hjust = 1,
           label = "italic(P) == 0.77", parse = TRUE, size = 5)

## Figure 4C: Number of jumps (mean ± SE of individual means)
summary_jumps <- movement %>%
  filter(measure == "Jumps") %>%
  group_by(ID, Range) %>%
  summarise(mean_jumps = mean(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(Range) %>%
  summarise(Mean = mean(mean_jumps, na.rm = TRUE),
            SE   = sd(mean_jumps, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

fig4c <- ggplot(summary_jumps, aes(x = Range, y = Mean, fill = Range)) +
  geom_col(width = 0.5, alpha = 0.9, color = "black") +
  geom_point(size = 5, shape = 18, color = "black") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.15, linewidth = 1) +
  scale_fill_manual(values = c("Historical" = "grey", "Expanded" = "orange")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Mean number of jumps") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 16),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  annotate("text", x = 1.5, y = max(summary_jumps$Mean + summary_jumps$SE) * 1.05,
           label = "italic(P) == 0.040", parse = TRUE, size = 5)

figure4 <- ggarrange(fig4a, fig4b, fig4c,
                     labels = c("A", "B", "C"), hjust = -0.5, vjust = 2,
                     ncol = 3, align = "h",
                     font.label = list(size = 18, color = "black", face = "bold", family = "Arial"))
figure4

#ggsave("outputs/figures/figure4.pdf", plot = figure4,height = 11, width = 25, units = "cm", device = cairo_pdf)


#### Figure 5: Proportion of aggressive calls ####
acoustics$population <- factor(acoustics$population,
                               levels = c("expanded", "historical"),
                               labels = c("Expanded", "Historical"))

figure5 <- ggplot(acoustics, aes(x = population, y = proportion_aggressive, fill = population)) +
  geom_boxplot(width = 0.5, alpha = 0.9, color = "black", linewidth = 0.8,
               outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 123),
             size = 3, color = "black", alpha = 0.6) +
  scale_fill_manual(values = c("Historical" = "grey", "Expanded" = "orange")) +
  labs(x = NULL, y = "Proportion of aggressive calls") +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 16),
    axis.title = element_text(face = "bold", size = 15),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  annotate("text", x = 1.5, y = max(acoustics$proportion_aggressive, na.rm = TRUE) * 1.05,
           label = "italic(P) < 0.001", parse = TRUE, size = 5)
figure5

#ggsave("outputs/figures/figure5.pdf", plot = figure5,width = 10, height = 12, units = "cm", device = cairo_pdf)


#### Figure S1: Individual means of repeatable behaviors ####
ind_means <- function(dat, value_col) {
  dat %>%
    filter(!is.na(.data[[value_col]])) %>%
    group_by(ID, Range) %>%
    summarise(ind_mean = mean(.data[[value_col]]), .groups = "drop")
}

exit_ind <- ind_means(emergence, "latency")
move_ind <- ind_means(filter(movement, measure == "Latency"), "value")
jump_ind <- ind_means(filter(movement, measure == "Jumps"), "value")

## Figure S1A: Latency to exit shelter
S1_exit <- ggplot(exit_ind, aes(x = Range, y = ind_mean, fill = Range)) +
  geom_boxplot(width = 0.5, alpha = 0.9, color = "black", linewidth = 0.8,
               outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 123),
             size = 2.5, color = "black", alpha = 0.6) +
  scale_fill_manual(values = c("Expanded" = "orange", "Historical" = "grey")) +
  labs(x = NULL, y = "Mean latency to exit shelter (s)") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  annotate("text", x = 1.5, y = max(exit_ind$ind_mean) * 1.08,
           label = "italic(P) == 0.031", parse = TRUE, size = 5)

## Figure S1B: Latency to first movement
S1_move <- ggplot(move_ind, aes(x = Range, y = ind_mean, fill = Range)) +
  geom_boxplot(width = 0.5, alpha = 0.9, color = "black", linewidth = 0.8,
               outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 123),
             size = 2.5, color = "black", alpha = 0.6) +
  scale_fill_manual(values = c("Expanded" = "orange", "Historical" = "grey")) +
  labs(x = NULL, y = "Mean latency to first movement (s)") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  annotate("text", x = 1.5, y = max(move_ind$ind_mean) * 1.08,
           label = "italic(P) == 0.77", parse = TRUE, size = 5)

## Figure S1C: Number of jumps
S1_jumps <- ggplot(jump_ind, aes(x = Range, y = ind_mean, fill = Range)) +
  geom_boxplot(width = 0.5, alpha = 0.9, color = "black", linewidth = 0.8,
               outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 123),
             size = 2.5, color = "black", alpha = 0.6) +
  scale_fill_manual(values = c("Expanded" = "orange", "Historical" = "grey")) +
  labs(x = NULL, y = "Mean number of jumps") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  annotate("text", x = 1.5, y = max(jump_ind$ind_mean) * 1.08,
           label = "italic(P) == 0.040", parse = TRUE, size = 5)

figureS1 <- ggarrange(S1_exit, S1_move, S1_jumps,
                      labels = c("A", "B", "C"), ncol = 3, align = "h",
                      font.label = list(size = 18, color = "black", face = "bold"))
figureS1

#ggsave("outputs/figures/figureS1.pdf", plot = figureS1,height = 11, width = 25, units = "cm", device = cairo_pdf)
