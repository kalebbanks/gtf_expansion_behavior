# Behavioral divergence during rapid range expansion of a native amphibian

Authors: XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## Overview

This repository contains the data and R code used to analyze behavioral differences between historical and recently expanded populations of the green treefrog (*Hyla cinerea*).

### Associated study

Paper currently in review, will be updated later

---

## Repository structure

```
├── data/
│   ├── emergence_long.csv
│   ├── movement_long.csv
│   ├── acoustic_results.csv
│   ├── All_Traits.csv
│   ├── GTF_Behavior_Master_Sheet.csv
│   ├── agg.wav
│   ├── adv.wav
│   └── playback.wav
│
├── scripts/
│   ├── 01_behavior_analysis.R
│   └── 03_figures.R
│
├── outputs/
│   └── figures/        
│
└── README.md
```

---

## Data files

### gtf_behavior_master_sheet.csv 
Contains all meta data associated with trials and field collection

### acoustic_results.csv

Contains the results of the acoustic playback experiment. Each observation represents an individual frog and includes the proportion of aggressive calls produced during playback trials.

### all_traits.csv

Contains summarized behavioral trait data for each individual used in repeatability analyses, correlation analyses, and comparisons between historical and expanded populations.

### emergence_long.csv

Contains latency-to-exit measurements.

### movement_long.csv

Contains repeated movement behavior measurements, including:

* Latency to first movement
* Latency to exit
* Time spent on the floor
* Time spent moving
* Number of jumping events
* Number of movement events

---

## Analysis scripts

### 01_behavior_analysis.R

Performs the primary behavioral analyses, including:

* Repeatability analyses
* Comparisons between historical and expanded populations

### 02_trait_correlations.R

Tests correlations among behavioral traits to evaluate relationships between repeatable behaviors.

### 03_figures.R

Generates most figures presented in the manuscript using the processed data and statistical results.

---

## Software

Analyses were conducted in R version 4.5.2 (R Core Team 2024)


---

## Reproducing the analyses

1. Download this repository.
2. Create R project and put in folder with data, ouputs, and scripts folders
3. Run the scripts in numerical order:

```
01_behavior_analysis.R
02_trait_correlations.R
03_figures.R
```

---

## Contact

For questions regarding the data or analyses, please contact Owen M. Edwards (owenedwards@ou.edu)
