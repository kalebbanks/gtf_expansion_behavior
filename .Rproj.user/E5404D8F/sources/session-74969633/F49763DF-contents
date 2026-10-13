####################
#Statistical Analysis for "Behavioral divergence during rapid range expansion of a native amphibian" 
#Script: 02_trait_correlations
#Purpose: To see if repeatable traits are correlated with one another
#Inputs: data/
#Author: XXXXXXXX
#Date: 7/4/2026
####################

#####Packages
#install.packages("Hmisc")
library(Hmisc)


#####Load Data
#Traits together in one spreadsheet
All_Traits <- read.csv("data/all_traits.csv")


#####Correlation analysis

Trait_means <- All_Traits %>%
  mutate(
    LatencyExit = rowMeans(select(., LatencyExit1, LatencyExit2), na.rm = TRUE),
    LatencyMovement = rowMeans(select(., LatencyMovement1, LatencyMovement2), na.rm = TRUE),
    Jump = rowMeans(select(., Jump1, Jump2), na.rm = TRUE)
  )

traits <- Trait_means[, c("LatencyExit", "LatencyMovement", "Jump", "Aggressive")]

cor(traits, method = "spearman", use = "complete.obs")

rcorr(as.matrix(traits), type = "spearman")

