#####################
#Script: 03_NSDM_modeling.R
#Purpose:  
#Inputs: 
#Outputs: 
#Author:Kaleb M. Banks
#Date: 2025-09-03
#####################
#packages
library(dplyr)
library(sabinaNSDM)
library(terra)
library(covsel)
library(biomod2)
library(ecospat)
library(fs)
library(sgsR)

#load datasets
glob_occ <- read.csv("data/occ/glob_occ.csv")
glob_abs <- read.csv("data/occ/glob_abs.csv")
reg_occ <- read.csv("data/occ/reg_occ.csv")
reg_abs <- read.csv("data/occ/reg_abs.csv")


#Load rasters
scenario_1951_1959 <- rast("data/env_var/scenarios/1951_1959.tif")
scenario_1960_1969 <- rast("data/env_var/scenarios/1960_1969.tif")
scenario_1970_1979 <- rast("data/env_var/scenarios/1970_1979.tif")
scenario_1980_1989 <- rast("data/env_var/scenarios/1980_1989.tif")
scenario_1990_1999 <- rast("data/env_var/scenarios/1990_1999.tif")
scenario_2000_2009 <- rast("data/env_var/scenarios/2000_2009.tif")
scenario_2010_2019 <- rast("data/env_var/scenarios/2010_2019.tif")
scenario_2020_2024 <- rast("data/env_var/scenarios/2020_2024.tif")

reg_training_var <- rast("data/env_var/scenarios/reg_training_var.tif")
glob_training_var <- rast("data/env_var/scenarios/glob_training_var.tif")

#Model 

nsdm_input <- NSDM.InputData(
  SpeciesName = "h.cinerea",
  spp.data.global = glob_occ,
  spp.data.regional = reg_occ,
  expl.var.global = glob_training_var,
  expl.var.regional = reg_training_var,
  new.env = list(scenario_1951_1959, scenario_1960_1969, scenario_1970_1979, scenario_1980_1989,scenario_1990_1999, scenario_2000_2009, scenario_2010_2019, scenario_2020_2024),
  new.env.names = c("1951_1959", "1960_1969", "1970_1979", "1980_1989", "1990_1999", "2000_2009", "2010_2019", "2020_2024"),
  Background.Global = NULL,
  Background.Regional = NULL,
  Absences.Global = glob_abs,
  Absences.Regional = reg_abs
)


nsdm_finput <- NSDM.FormattingData(nsdm_input,
                                   nPoints = NULL, # number of background points
                                   Min.Dist.Global = "resolution",
                                   Min.Dist.Regional = "resolution",
                                   Background.method = "random", # method “random" or "stratified” to generate background points 
                                   save.output = TRUE) #save outputs locally



nsdm_selvars <- NSDM.SelectCovariates(nsdm_finput,
                                      maxncov.Global = 3,   # Max number of covariates to be selected at the global scale
                                      maxncov.Regional = 4, # Max number of covariates to be selected at the regional scale
                                      corcut = 0.99, #  correlation threshold
                                      algorithms = c("glm","gam","rf"),
                                      ClimaticVariablesBands = NULL, # covariate bands to be excluded in the covariate selection at the regional scale
                                      save.output = TRUE)

nsdm_global <- NSDM.Global(nsdm_selvars,
                           algorithms = c("MARS"),# Statistical algorithms used for modelling
                           CV.nb.rep = 10, # number of cross-validation repetitions
                           CV.perc = 0.8, # percentage of the data will be used for training in each cross-validation fold
                           metric.select.thresh = 0.9, #  AUC threshold to include replicates in the final ensemble model
                           CustomModelOptions = NULL, # Allows users to apply custom modelling options. 
                           save.output = TRUE, 
                           rm.biomod.folder = FALSE) # Remove the temporary folders created by `biomod2` 


######Eyeball test
#Current predictions
plot(terra::unwrap(nsdm_global$current.projections$Pred))
plot(terra::unwrap(nsdm_global$current.projections$Pred.bin.TSS))

#CV metrics 
summary(nsdm_global)


nsdm_regional <- NSDM.Regional(nsdm_selvars,
                               algorithms = c("MARS"),
                               CV.nb.rep = 10,
                               CV.perc = 0.8,
                               # metric.select.thresh = 0.8,
                               CustomModelOptions = NULL, 
                               save.output = TRUE,
                               rm.biomod.folder = FALSE)

#####Eyeball test
#Current predictions
plot(terra::unwrap(nsdm_regional$current.projections$Pred))
plot(terra::unwrap(nsdm_regional$current.projections$Pred.bin.TSS))

#CV metrics 
summary(nsdm_regional)

#Future global predictions 
plot(terra::unwrap(nsdm_regional$new.projections$Pred.Scenario[[1]]))
plot(terra::unwrap(nsdm_regional$new.projections$Pred.bin.TSS.Scenario[[1]]))
plot(terra::unwrap(nsdm_regional$new.projections$Pred.Scenario[[2]]))
plot(terra::unwrap(nsdm_regional$new.projections$Pred.bin.TSS.Scenario[[2]]))
plot(terra::unwrap(nsdm_regional$new.projections$Pred.Scenario[[3]]))
plot(terra::unwrap(nsdm_regional$new.projections$Pred.bin.TSS.Scenario[[3]]))
plot(terra::unwrap(nsdm_regional$new.projections$Pred.Scenario[[4]]))
plot(terra::unwrap(nsdm_regional$new.projections$Pred.bin.TSS.Scenario[[4]]))
plot(terra::unwrap(nsdm_regional$new.projections$Pred.Scenario[[5]]))
plot(terra::unwrap(nsdm_regional$new.projections$Pred.bin.TSS.Scenario[[5]]))
plot(terra::unwrap(nsdm_regional$new.projections$Pred.Scenario[[6]]))
plot(terra::unwrap(nsdm_regional$new.projections$Pred.bin.TSS.Scenario[[6]]))
plot(terra::unwrap(nsdm_regional$new.projections$Pred.Scenario[[7]]))
plot(terra::unwrap(nsdm_regional$new.projections$Pred.bin.TSS.Scenario[[7]]))
plot(terra::unwrap(nsdm_regional$new.projections$Pred.Scenario[[8]]))
plot(terra::unwrap(nsdm_regional$new.projections$Pred.bin.TSS.Scenario[[8]]))


nsdm_covariate <- NSDM.Covariate(nsdm_global,
                                 algorithms = c("MARS", "GBM", "GAM", "GLM"),
                                 rm.corr=TRUE,
                                 CV.nb.rep = 5,
                                 CV.perc = 0.8,
                                 # metric.select.thresh = 0.8,
                                 CustomModelOptions = NULL,
                                 save.output = TRUE,
                                 rm.biomod.folder = FALSE)

plot(terra::unwrap(nsdm_covariate$current.projections$Pred))
plot(terra::unwrap(nsdm_covariate$current.projections$Pred.bin.TSS), 
     main = "Binary Prediction, TSS Threshold (2014-2024)")
#CV metrics 
summary(nsdm_covariate)


plot(terra::unwrap(nsdm_covariate$new.projections$Pred.Scenario[[1]]))
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.bin.TSS.Scenario[[1]]), main = "Binary Prediction, TSS Threshold (1950-1959)")
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.Scenario[[2]]))
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.bin.TSS.Scenario[[2]]), main = "Binary Prediction, TSS Threshold (1960-1969)")
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.Scenario[[3]]))
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.bin.TSS.Scenario[[3]]), main = "Binary Prediction, TSS Threshold (1970-1979)")
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.Scenario[[4]]))
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.bin.TSS.Scenario[[4]]), main = "Binary Prediction, TSS Threshold (1980-1989)")
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.Scenario[[5]]))
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.bin.TSS.Scenario[[5]]), main = "Binary Prediction, TSS Threshold (1990-1999)")
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.Scenario[[6]]))
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.bin.TSS.Scenario[[6]]), main = "Binary Prediction, TSS Threshold (2000-2010)")
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.Scenario[[7]]))
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.bin.TSS.Scenario[[7]]), main = "Binary Prediction, TSS Threshold (2010-2020)")
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.Scenario[[8]]))
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.bin.TSS.Scenario[[8]]), main = "Binary Prediction, TSS Threshold (2020-2024)")






######################################################################
#####Demo to show nate and Owen
#covariate model
plot(terra::unwrap(nsdm_covariate$current.projections$Pred))
plot(terra::unwrap(nsdm_covariate$current.projections$Pred.bin.TSS), 
     main = "Binary Prediction, TSS Threshold (2014-2024)")
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.bin.TSS.Scenario[[1]]), main = "Binary Prediction, TSS Threshold (1950-1959)")
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.bin.TSS.Scenario[[2]]), main = "Binary Prediction, TSS Threshold (1960-1969)")
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.bin.TSS.Scenario[[3]]), main = "Binary Prediction, TSS Threshold (1970-1979)")
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.bin.TSS.Scenario[[4]]), main = "Binary Prediction, TSS Threshold (1980-1989)")
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.bin.TSS.Scenario[[5]]), main = "Binary Prediction, TSS Threshold (1990-1999)")
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.bin.TSS.Scenario[[6]]), main = "Binary Prediction, TSS Threshold (2000-2010)")
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.bin.TSS.Scenario[[7]]), main = "Binary Prediction, TSS Threshold (2010-2020)")
plot(terra::unwrap(nsdm_covariate$new.projections$Pred.bin.TSS.Scenario[[8]]), main = "Binary Prediction, TSS Threshold (2020-2024)")


covar_models_var_import <- read.csv("Results/Covariate/Values/h.cinerea_indvar.csv")
covar_indiv_models <- read.csv("Results/Covariate/Values/h.cinerea_replica.csv")
covar_indiv_models <- covar_indiv_models[covar_indiv_models$metric.eval == "ROC", ]
covar_indiv_models <- covar_indiv_models[covar_indiv_models$validation >= 0.80, ]
covar_models_var_import <- covar_models_var_import[
  covar_models_var_import$full.name %in% covar_indiv_models$full.name, 
]

covar_ensemble_var_import_mean <- aggregate(var.imp ~ expl.var, data = covar_models_var_import, FUN = mean, na.rm = TRUE)
covar_ensemble_var_import_sd <- aggregate(var.imp ~ expl.var, data = covar_models_var_import, FUN = sd, na.rm = TRUE)

covar_ensemble_var_import <- merge(covar_ensemble_var_import_mean, covar_ensemble_var_import_sd, by = "expl.var")
colnames(covar_ensemble_var_import) <- c("Variable", "MeanVarImp", "SDVarImp")

covar_ensemble_var_import


library(ggplot2)
library(rprojroot)
library(biomod2)

#########RUN CODE HERE TO SEE VARIABLE IMPORTANCE
ggplot(covar_ensemble_var_import, aes(x = reorder(Variable, -MeanVarImp), y = MeanVarImp)) +
  geom_col(fill = "grey40", width = 0.6) +
  geom_errorbar(aes(ymin = MeanVarImp - SDVarImp, ymax = MeanVarImp + SDVarImp), 
                width = 0.2) +
  labs(x = NULL, y = "Mean Variable Importance") +
  theme_classic(base_size = 13)

load("Results/Covariate/Models/h.cinerea/h.cinerea.AllModels.ensemble.models.out")
covar_ensemble_BIOMOD_output <- h.cinerea.AllModels.ensemble.models.out
rm(h.cinerea.AllModels.ensemble.models.out)

setwd("Results/Covariate/Models/")

covar_response_curves <- biomod2::bm_PlotResponseCurves(
  bm.out = covar_ensemble_BIOMOD_output
  , models.chosen = biomod2::get_built_models(covar_ensemble_BIOMOD_output)
  , fixed.var = 'median')

covar_response_curves <- covar_response_curves$plot$data
covar_response_curves <- covar_response_curves[covar_response_curves$pred.name != "R.areolata_EMcvByROC_mergedData_mergedRun_mergedAlgo", ]
covar_response_curves <- covar_response_curves[, !(names(covar_response_curves) %in% c("pred.name", "id"))]
covar_response_curves$model <- "Covariate"
setwd(find_rstudio_root_file())

load("Results/Global/Models/h.cinerea/h.cinerea.AllModels.ensemble.models.out")
glob_ensemble_BIOMOD_output <- h.cinerea.AllModels.ensemble.models.out
rm(h.cinerea.AllModels.ensemble.models.out)

setwd("Results/Global/Models/")

glob_response_curves <- biomod2::bm_PlotResponseCurves(
  bm.out = glob_ensemble_BIOMOD_output
  , models.chosen = biomod2::get_built_models(glob_ensemble_BIOMOD_output)
  , fixed.var = 'median')

glob_response_curves <- glob_response_curves$plot$data
glob_response_curves <- glob_response_curves[glob_response_curves$pred.name != "h.cinerea_EMcvByROC_mergedData_mergedRun_mergedAlgo", ]
glob_response_curves <- glob_response_curves[, !(names(glob_response_curves) %in% c("pred.name", "id"))]
glob_response_curves$model <- "Global"
#reset working directory back to project root
setwd(find_rstudio_root_file())

response_curves <- rbind(glob_response_curves, covar_response_curves)
response_curves <- response_curves %>%
  group_by(expl.name) %>%
  filter(pred.val >= 0, pred.val <= 1) %>%
  ungroup()

library(ggplot2)

##################Go here for response curves
ggplot(response_curves, aes(x = expl.val, y = pred.val, color = model)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ expl.name, scales = "free_x") +
  labs(x = "Predictor value", y = "Predicted response", color = "Model") +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")



############## SEE rangewide Global model projection
load("Results/Global/Models/h.cinerea/h.cinerea.AllModels.ensemble.models.out")
glob_ensemble_BIOMOD_output <- h.cinerea.AllModels.ensemble.models.out

global_ensemble_BIOMOD_output <- h.cinerea.AllModels.ensemble.models.out

BIOMOD_Projection(
  global_ensemble_BIOMOD_output,
  proj.name,
  new.env.xy = glob_training_var,
  models.chosen = "all",
  metric.binary = NULL,
  metric.filter = NULL,
  compress = TRUE,
  build.clamping.mask = TRUE,
  nb.cpu = 1,
  seed.val = NULL,
)
setwd("Results/Global/Models")
range_wide <- BIOMOD_EnsembleForecasting(
  global_ensemble_BIOMOD_output,
  proj.name = 'CurrentEM',
  new.env = glob_training_var,
  new.env.xy = NULL,
  models.chosen = "all",
  metric.binary = NULL,
  metric.filter = NULL,
  compress = TRUE,
  nb.cpu = 1,
  na.rm = TRUE,
)

biomod2::plot(range_wide,
              main="Crawfish Frog Species Distribution Model Map")







