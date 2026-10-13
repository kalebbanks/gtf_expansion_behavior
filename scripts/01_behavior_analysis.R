# Script: 01_behavior_analysis.R
# Statistical analyses for "Behavioral divergence during rapid range expansion of a native amphibian"
# Purpose: Repeatability (Table 1, Table S3), range comparisons (Tables 2-3), trial timing (Table S2)
# Inputs: data/emergence_long.csv, data/movement_long.csv, data/acoustic_results.csv,
#         data/GTF_Behavior_Master_Sheet.csv
# Author: XXXXXXXX
# Date: 10/12/2026
# Recorded outputs are shown as comments below each result.

#### Packages ####
# install.packages(c("survival", "coxme", "lme4", "glmmTMB", "rptR","tidyr", "dplyr", "performance"))
library(survival)
library(coxme)
library(lme4)
library(glmmTMB)
library(rptR)
library(tidyr)
library(dplyr)
library(performance)

#### Data ####
emergence  <- read.csv("data/emergence_long.csv")   # latency to exit shelter (max 1200 s)
movement <- read.csv("data/movement_long.csv")       # arena movement measures (max 1800 s)
acoustics  <- read.csv("data/acoustic_results.csv")
sched <- read.csv("data/GTF_Behavior_Master_Sheet.csv", colClasses = "character",check.names = FALSE) %>%
             rename_with(trimws) %>%
             mutate(County = trimws(County),Range = ifelse(County %in% c("Vanderburgh", "Crawford", "Hopkins", "Webster"),"Expanded", "Historical"))

#### Functions ####
# Repeatability for Cox models (residual variance = pi^2/6)
icc  <- function(v, p) v / (v + p + pi^2/6)          # Rv2: individual only
icc3 <- function(v, p) (v + p) / (v + p + pi^2/6)    # Rv3: multi-cluster (Nakagawa et al. 2026)

# Profile-likelihood grid for the ID variance in a Cox model with Population and ID
profile_id <- function(formula, data) {
  grid <- seq(0.001, 5, length.out = 200)
  out <- lapply(grid, function(v) {
    tryCatch({
      fit <- coxme(formula, data = data, vfixed = list(NULL, v))   # Population estimated, ID fixed at v
      data.frame(var_id = v,
                 var_pop = VarCorr(fit)$Population[[1]],
                 loglik = fit$loglik["Integrated"])
    }, error = function(e) NULL)
  })
  do.call(rbind, out)
}

# CI from the profile: values of the ID variance within the 95% likelihood cutoff
profile_ci <- function(prof, model, fun) {
  inside <- prof[2 * (model$loglik["Integrated"] - prof$loglik) <= qchisq(0.95, 1), ]
  c(lower = fun(min(inside$var_id), inside$var_pop[which.min(inside$var_id)]),
    upper = fun(max(inside$var_id), inside$var_pop[which.max(inside$var_id)]))
}

#Formating time functions
to_hour <- function(x) {
  x <- as.numeric(gsub(":", "", x))      # "1622" or "14:56" -> 1622 / 1456
  floor(x / 100) + (x %% 100) / 60       # HHMM -> decimal hours
}

fmt_dates <- function(d) {
  d <- sort(unique(na.omit(d)))
  if (length(d) == 1) format(d, "%d %b") else paste(format(min(d), "%d %b"), "-", format(max(d), "%d %b"))
}

#### Table 1 and Table S3: Repeatability ####
# Each behavior: model with ID and Population, variance components (Table S3),
# Rv2 (Table S3), Rv3 (Table 1), and LRT for the ID variance (P in Table 1).

## Latency to exit shelter (Cox mixed model)
exit_rpt <- coxme(Surv(latency, event) ~ 1 + (1 | Population) + (1 | ID), data = emergence)
var_id_exit  <- VarCorr(exit_rpt)$ID[[1]]
var_pop_exit <- VarCorr(exit_rpt)$Population[[1]]
c(var_id = var_id_exit, var_pop = var_pop_exit, residual = pi^2/6)
# var_id 0.983, var_pop 0.143, residual 1.645

exit_noID <- coxme(Surv(latency, event) ~ 1 + (1 | Population), data = emergence)
anova(exit_noID, exit_rpt)
# Chisq = 14.64, df = 1, P = 0.00013

exit_prof <- profile_id(Surv(latency, event) ~ 1 + (1 | Population) + (1 | ID), emergence)
c(Rv2 = icc(var_id_exit, var_pop_exit),  profile_ci(exit_prof, exit_rpt, icc))
# Rv2 = 0.35 (0.17-0.54)
c(Rv3 = icc3(var_id_exit, var_pop_exit), profile_ci(exit_prof, exit_rpt, icc3))
# Rv3 = 0.41 (0.22-0.58)

# Latency to first movement; event = 0 if the frog did not move within 1800 s
lat <- movement %>%
  filter(measure == "Latency", !is.na(value)) %>%
  mutate(event = ifelse(value == 1800, 0, 1))

## Latency to first movement (Cox mixed model)
move_rpt <- coxme(Surv(value, event) ~ 1 + (1 | Population) + (1 | ID), data = lat)
var_id_move  <- VarCorr(move_rpt)$ID[[1]]
var_pop_move <- VarCorr(move_rpt)$Population[[1]]
c(var_id = var_id_move, var_pop = var_pop_move, residual = pi^2/6)
# var_id 0.371, var_pop 0.048, residual 1.645

move_noID <- coxme(Surv(value, event) ~ 1 + (1 | Population), data = lat)
anova(move_noID, move_rpt)
# Chisq = 4.93, df = 1, P = 0.026

move_prof <- profile_id(Surv(value, event) ~ 1 + (1 | Population) + (1 | ID), lat)
c(Rv2 = icc(var_id_move, var_pop_move),  profile_ci(move_prof, move_rpt, icc))
# Rv2 = 0.18 (0.03-0.36)
c(Rv3 = icc3(var_id_move, var_pop_move), profile_ci(move_prof, move_rpt, icc3))
# Rv3 = 0.20 (0.05-0.38)

## Movement events (Poisson GLMM)
set.seed(123)
events_rpt <- rpt(value ~ 1 + (1 | ID) + (1 | Population),
                  grname = c("ID", "Population"),
                  data = filter(movement, measure == "Numb_exploration_occurances"),
                  datatype = "Poisson")
events_rpt
# Rv2 (R for ID, original scale) = 0.04 (0.00-0.13), LRT P = 0.158

events_rpt$R["R_org", "ID"] + events_rpt$R["R_org", "Population"]
quantile(events_rpt$R_boot_org$ID + events_rpt$R_boot_org$Population,
         c(0.025, 0.975), na.rm = TRUE)
# Rv3 = 0.06 (0.00-0.15)

events_var <- rpt(value ~ 1 + (1 | ID) + (1 | Population),
                  grname = c("ID", "Population", "Residual"),
                  data = filter(movement, measure == "Numb_exploration_occurances"),
                  datatype = "Poisson", ratio = FALSE, nboot = 0, npermut = 0)
events_var$R["R_link", ]
# link scale: ID 0.227, Population 0.134, Residual 1.734

## Jumps (Poisson GLMM)
set.seed(123)
jumps_rpt <- rpt(value ~ 1 + (1 | ID) + (1 | Population),
                 grname = c("ID", "Population"),
                 data = filter(movement, measure == "Jumps"),
                 datatype = "Poisson")
jumps_rpt
# Rv2 = 0.21 (0.05-0.43), LRT P = 0.004

jumps_rpt$R["R_org", "ID"] + jumps_rpt$R["R_org", "Population"]
quantile(jumps_rpt$R_boot_org$ID + jumps_rpt$R_boot_org$Population,
         c(0.025, 0.975), na.rm = TRUE)
# Rv3 = 0.36 (0.16-0.51)

jumps_var <- rpt(value ~ 1 + (1 | ID) + (1 | Population),
                 grname = c("ID", "Population", "Residual"),
                 data = filter(movement, measure == "Jumps"),
                 datatype = "Poisson", ratio = FALSE, nboot = 0, npermut = 0)
jumps_var$R["R_link", ]
# link scale: ID 0.267, Population 0.202, Residual 0.514

## Time spent moving (Gaussian LMM)
set.seed(123)
moving_rpt <- rpt(value ~ 1 + (1 | ID) + (1 | Population),
                  grname = c("ID", "Population"),
                  data = filter(movement, measure == "Exploration"),
                  datatype = "Gaussian")
moving_rpt
# Rv2 = 0.00 (0.00-0.18), LRT P = 1.000
# Singular fit: ID variance estimated at zero

moving_rpt$R$ID + moving_rpt$R$Population
quantile(moving_rpt$R_boot$ID + moving_rpt$R_boot$Population,
         c(0.025, 0.975), na.rm = TRUE)
# Rv3 = 0.25 (0.06-0.50)

moving_var <- rpt(value ~ 1 + (1 | ID) + (1 | Population),
                  grname = c("ID", "Population", "Residual"),
                  data = filter(movement, measure == "Exploration"),
                  datatype = "Gaussian", ratio = FALSE, nboot = 0, npermut = 0)
moving_var$R
# ID 0, Population 3076.3, Residual 9339.1

## Time on floor (Gaussian LMM)
set.seed(123)
floor_rpt <- rpt(value ~ 1 + (1 | ID) + (1 | Population),
                 grname = c("ID", "Population"),
                 data = filter(movement, measure == "Floor"),
                 datatype = "Gaussian")
floor_rpt
# Rv2 = 0.12 (0.00-0.32), LRT P = 0.137

floor_rpt$R$ID + floor_rpt$R$Population
quantile(floor_rpt$R_boot$ID + floor_rpt$R_boot$Population,
         c(0.025, 0.975), na.rm = TRUE)
# Rv3 = 0.14 (0.00-0.34)

floor_var <- rpt(value ~ 1 + (1 | ID) + (1 | Population),
                 grname = c("ID", "Population", "Residual"),
                 data = filter(movement, measure == "Floor"),
                 datatype = "Gaussian", ratio = FALSE, nboot = 0, npermut = 0)
floor_var$R
# ID 37433.4, Population 6627.1, Residual 264231.2

# Repeatable behaviors (LRT P < 0.05): latency to exit, latency to first movement, jumps.
# Only these were compared between ranges (Tables 2-3).


#### Table 2: Range comparisons for repeatable behaviors ####

## Jumps (negative binomial GLMM)
jumps <- filter(movement, measure == "Jumps")

# Poisson model first, to check overdispersion
jumps_pois <- glmer(value ~ Range + Mass + (1 | Population) + (1 | ID),
                    data = jumps, family = poisson)
sum(residuals(jumps_pois, type = "pearson")^2) / df.residual(jumps_pois)
# 1.59: overdispersed, so a negative binomial model was used

jumps_nb <- glmmTMB(value ~ Range + Mass + (1 | Population) + (1 | ID),
                    data = jumps, family = nbinom2)
summary(jumps_nb)
#                 Estimate Std. Error z value Pr(>|z|)
# (Intercept)      2.60422    0.38991   6.679  < 0.001
# Rangehistorical -0.60433    0.29344  -2.059   0.0395
# Mass            -0.08706    0.04713  -1.847   0.0647
# Random-effect variances: Population 0.118, ID 0.180

## Latency to first movement (Cox mixed model)
move_range <- coxme(Surv(value, event) ~ Range + Mass + (1 | Population) + (1 | ID),
                    data = lat)
move_range
#                    coef exp(coef) se(coef)     z    p
# Rangehistorical -0.0798    0.9233   0.2783 -0.29 0.77
# Mass            -0.0518    0.9495   0.0560 -0.93 0.35
# Random-effect variances: Population 0.060, ID 0.364

## Latency to exit shelter (Cox mixed model)
exit_range <- coxme(Surv(latency, event) ~ Range + Mass + (1 | Population) + (1 | ID),
                    data = emergence)
exit_range
#                    coef exp(coef) se(coef)     z      p
# RangeHistorical -0.6594    0.5172   0.3053 -2.16 0.0308
# Mass            -0.1340    0.8746   0.0724 -1.85 0.0643
# Random-effect variances: Population 0.029, ID 0.923

## Trials with no response (frog did not exit or did not move), by range
emergence %>%
  group_by(Range) %>%
  summarise(n_trials = n(), n_no_response = sum(event == 0),
            pct_no_response = round(100 * mean(event == 0), 1))

lat %>%
  group_by(Range) %>%
  summarise(n_trials = n(), n_no_response = sum(event == 0),
            pct_no_response = round(100 * mean(event == 0), 1))


#### Table 3: Proportion of aggressive calls ####

## Chorus SPL did not differ between ranges (Welch's t-test)
t.test(chorus_spl ~ population, data = acoustics)

acoustics$non_aggressive <- acoustics$total_calls - acoustics$aggressive
acoustics$OLRE <- seq_len(nrow(acoustics))   # observation-level random effect

# Binomial GLMM without OLRE, to check overdispersion
agg_binom <- glmer(cbind(aggressive, non_aggressive) ~
                     population + chorus_spl + body_temp + mass + (1 | county),
                   data = acoustics, family = binomial)
check_overdispersion(agg_binom)
# Overdispersion detected, so an observation-level random effect was added

agg_model <- glmer(cbind(aggressive, non_aggressive) ~
                     population + scale(chorus_spl) + scale(body_temp) + scale(mass) +
                     (1 | county) + (1 | OLRE),
                   data = acoustics, family = binomial,
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(agg_model)

#                      Estimate Std. Error z value Pr(>|z|)
# (Intercept)           -0.3927     0.6971  -0.563 0.573
# populationhistorical  -4.2578     1.1147  -3.820 < 0.001
# scale(chorus_spl)     -1.2702     0.4696  -2.705 0.007
# scale(body_temp)       0.3356     0.5615   0.598 0.550
# scale(mass)           -0.4447     0.4616  -0.963 0.335
# Random-effect variances: OLRE 15.36, county 0 (singular fit); n = 83 frogs, 8 populations

# Sensitivity check: same model without the population random effect
agg_noPop <- glmer(cbind(aggressive, non_aggressive) ~
                     population + scale(chorus_spl) + scale(body_temp) + scale(mass) +
                     (1 | OLRE),
                   data = acoustics, family = binomial,
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
cbind(with_pop = fixef(agg_model), without_pop = fixef(agg_noPop))
# Fixed-effect estimates are nearly identical


#### Table S2: Collection dates and trial start times by population ####

# One row per frog x trial, with date and start time
sched_long <- bind_rows(
  sched %>% transmute(County, Range, Trial = "Arena movement 1",
                      date = `Exploration trial 1 date`, start = `Exploration Trial 1 start`),
  sched %>% transmute(County, Range, Trial = "Shelter emergence 1",
                      date = `Boldness Trial 1 date`, start = `Boldness Trial 1 start`),
  sched %>% transmute(County, Range, Trial = "Arena movement 2",
                      date = `Exploration trial 2 date`, start = `Exploration Trial 2 start`),
  sched %>% transmute(County, Range, Trial = "Shelter emergence 2",
                      date = `Boldness Trial 2 date`, start = `Boldness Trial 2 start`)
) %>%
  mutate(date = as.Date(date, format = "%m/%d/%y"),
         hour = to_hour(start),
         Trial = factor(Trial, levels = c("Arena movement 1", "Shelter emergence 1","Arena movement 2", "Shelter emergence 2")))

sched_table <- sched_long %>%
  group_by(Range, County, Trial) %>%
  summarise(cell = sprintf("%.1f (%.1f–%.1f)",
                           mean(hour, na.rm = TRUE),
                           min(hour, na.rm = TRUE),
                           max(hour, na.rm = TRUE)),
            .groups = "drop") %>%
  pivot_wider(names_from = Trial, values_from = cell) %>%
  left_join(
    sched %>%
      mutate(collected = as.Date(`Date Collected`, format = "%m/%d/%y")) %>%
      group_by(County) %>%
      summarise(first_collected = min(collected, na.rm = TRUE),
                `Date collected` = fmt_dates(collected), .groups = "drop"),
    by = "County"
  ) %>%
  mutate(Range = factor(Range, levels = c("Historical", "Expanded"))) %>%
  arrange(Range, first_collected) %>%
  select(Range, County, `Date collected`,
         `Arena movement 1`, `Shelter emergence 1`,
         `Arena movement 2`, `Shelter emergence 2`)

sched_table
