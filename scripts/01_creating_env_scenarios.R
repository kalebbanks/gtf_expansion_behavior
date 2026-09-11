#####################
#Script: 02_creating_env_scenarios.R
#Purpose: Edit, mask, and crop raw raster files and then assemble spatrasters for modeling. 
#Inputs: data/env_var/scenarios_raw
#Outputs: data/env_var/scenarios_spatraster
#Author:Kaleb M. Banks
#Date: 2025-09-03
#####################

#####NOTE TO FUTURE SELF (8/26/2026): THERE IS A LOT OF INEFFICIENIES IN THIS CODE SUCH AS SAVING FILES TO DRIVE WHEN YOU DIDNT HAVE TO. CAREFULLY REVIEW THIS CODE AND FIND A WAY TO BETTER OPTIMIZE IT BEFORE YOU PUBLISH IT



#####Packages
#install.packages("terra")
install.packages("fastbioclim")
library(terra)
library(geodata)
library(fastbioclim)

#####Functions
make_bioclim <- function(years,
                         bios = 1:19,
                         monthly_dir = "data/env_var/monthly_raw",
                         out_dir = "data/env_var/bioclim",
                         mask_poly = NULL) {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  period_label <- paste0(min(years), "-", max(years))
  
  # ---- Helper: locate the exact folder for one variable + period, and read its files ----
  get_var_stack <- function(var) {
    
    folder_name <- paste0("wc2.1_cruts4.09_2.5m_", var, "_", period_label)
    folder_path <- file.path(monthly_dir, folder_name)
    
    if (!dir.exists(folder_path)) {
      stop(var, ": expected folder not found: ", folder_path,
           "\nCheck that a folder matching this exact year range exists.")
    }
    
    files <- list.files(folder_path, pattern = "\\.tif$", full.names = TRUE)
    
    if (length(files) == 0) {
      stop(var, ": folder exists but contains no .tif files: ", folder_path)
    }
    
    # extract the month from filenames like ..._2014-01.tif
    file_months <- as.integer(sub(".*-(\\d{2})\\.tif$", "\\1", files))
    
    if (any(is.na(file_months))) {
      stop(var, ": could not parse month from one or more filenames in ", folder_path)
    }
    
    # sort chronologically by filename so month order is predictable
    ord <- order(files)
    
    list(stack = terra::rast(files[ord]), months = file_months[ord])
  }
  
  tmin_data <- get_var_stack("tmin")
  tmax_data <- get_var_stack("tmax")
  prec_data <- get_var_stack("prec")
  
  # ---- Flag incomplete data rather than silently proceeding ----
  n_expected <- length(years) * 12
  
  check_complete <- function(data, var) {
    if (length(data$months) != n_expected) {
      stop(var, ": found ", length(data$months), " files, expected ", n_expected,
           " for ", period_label, ". Check the folder contents.")
    }
  }
  check_complete(tmin_data, "tmin")
  check_complete(tmax_data, "tmax")
  check_complete(prec_data, "prec")
  
  # ---- Flag invalid bios up front rather than letting derive_bioclim fail deep in the call ----
  if (any(!bios %in% 1:19)) {
    stop("bios must be a subset of 1:19. Received: ", paste(bios, collapse = ", "))
  }
  
  # ---- Step 1: average monthly values across years into a 12-month climatology ----
  tmin_avg <- calculate_average(
    x = tmin_data$stack, index = tmin_data$months,
    output_dir = file.path(out_dir, paste0("tmin_avg_", period_label)),
    overwrite = TRUE, output_names = "tmin_avg"
  )
  
  tmax_avg <- calculate_average(
    x = tmax_data$stack, index = tmax_data$months,
    output_dir = file.path(out_dir, paste0("tmax_avg_", period_label)),
    overwrite = TRUE, output_names = "tmax_avg"
  )
  
  prec_avg <- calculate_average(
    x = prec_data$stack, index = prec_data$months,
    output_dir = file.path(out_dir, paste0("prec_avg_", period_label)),
    overwrite = TRUE, output_names = "prec_avg"
  )
  
  # ---- Step 2: derive only the requested bioclim variables from that climatology ----
  bio_out_dir <- file.path(out_dir, paste0("bioclim_", period_label))
  
  bio_vars <- derive_bioclim(
    bios = bios,
    tmin = tmin_avg,
    tmax = tmax_avg,
    prcp = prec_avg,
    output_dir = bio_out_dir,
    overwrite = TRUE
  )
  
  # ---- Step 3 (Optional): Crop and Mask to Shapefile ----
  if (!is.null(mask_poly)) {
    
    # Load vector if passed as a filepath string
    if (is.character(mask_poly)) {
      mask_poly <- terra::vect(mask_poly)
    }
    
    # Reproject mask to match raster CRS if they differ
    if (terra::crs(bio_vars) != terra::crs(mask_poly)) {
      mask_poly <- terra::project(mask_poly, terra::crs(bio_vars))
    }
    
    # Crop bounding box and mask exact shape (in correct sequence)
    bio_vars <- terra::crop(bio_vars, mask_poly)
    bio_vars <- terra::mask(bio_vars, mask_poly)
    
    # Overwrite on disk with cropped/masked versions
    terra::writeRaster(bio_vars, 
                       filename = file.path(bio_out_dir, paste0(names(bio_vars), ".tif")), 
                       overwrite = TRUE)
  }
  
  message("Saved bio", paste(bios, collapse = ", "), " for ", period_label,
          " to: ", bio_out_dir)
  return(bio_vars)
}


#####Load shapefiles
regional_mask <- vect('data/shapefiles/ok_selected_counties.shp')
global_mask <- vect('data/shapefiles/usa_selected_states.shp')

#####Load Raster
dist_rivers <- rast('data/env_var/dis_rivers.tif')


######Load all historical monthly worldclim data
#Also an additional warning, this will use about 37 gb of space on your local drive. Later in the code after the bioclim variables are made, we will delete it. 
#might be an easier way to do this, but the code works. 
out_dir <- "data/env_var/monthly_raw"
base_url <- "https://geodata.ucdavis.edu/climate/worldclim/2_1/hist/cts4.09"
variables <- c("tmin", "tmax", "prec")
periods <- c("1950-1959","1960-1969","1970-1979","1980-1989","1990-1999","2000-2009","2010-2019","2020-2024")

#build full list of files to download
file_grid <- expand.grid(variable = variables, period = periods,
                         stringsAsFactors = FALSE)

file_grid$filename <- paste0("wc2.1_cruts4.09_2.5m_", file_grid$variable,
                             "_", file_grid$period, ".zip")
file_grid$url  <- file.path(base_url, file_grid$filename)
file_grid$dest <- file.path(out_dir, file_grid$filename)

#download loop
for (i in seq_len(nrow(file_grid))) {
  
  url  <- file_grid$url[i]
  dest <- file_grid$dest[i]
  
  if (file.exists(dest)) {
    message("Already downloaded, skipping: ", basename(dest))
    next
  }
  
  message("Downloading (", i, "/", nrow(file_grid), "): ", basename(dest))
  
  result <- tryCatch({
    download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
    TRUE
  }, error = function(e) {
    message("  FAILED: ", basename(dest), " -- ", conditionMessage(e))
    FALSE
  })
  
  if (!result) next
}

#unzip each into its own folder
for (i in seq_len(nrow(file_grid))) {
  
  dest <- file_grid$dest[i]
  if (!file.exists(dest)) next
  
  unzip_dir <- file.path(out_dir, tools::file_path_sans_ext(basename(dest)))
  dir.create(unzip_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("Unzipping: ", basename(dest))
  unzip(dest, exdir = unzip_dir)

  extracted_files <- list.files(unzip_dir, pattern = "\\.tif$")
  if (length(extracted_files) > 0) {
    file.remove(dest)
    message("Deleted zip: ", basename(dest))
  } else {
    message("WARNING: no .tif files found after unzipping, keeping zip: ",
            basename(dest))
  }
}


######Make folder group 2014-2024, which will be the years to train the model
#We already downloaded these, just need to get in proper folders
base_dir <- "data/env_var/monthly_raw"
if (!dir.exists(base_dir)) {
  dir.create(base_dir, recursive = TRUE)
}

variables <- c("tmin", "tmax", "prec")

for (var in variables) {
  
  # Source folders (already unzipped, named to match the zip files)
  src_2010s <- file.path(base_dir,
                         paste0("wc2.1_cruts4.09_2.5m_", var, "_2010-2019"))
  src_2020s <- file.path(base_dir,
                         paste0("wc2.1_cruts4.09_2.5m_", var, "_2020-2024"))
  
  # Destination folder for the new combined period
  dest_dir <- file.path(base_dir,
                        paste0("wc2.1_cruts4.09_2.5m_", var, "_2014-2024"))
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- Pull 2014-2019 files out of the 2010-2019 folder ----
  files_2010s <- list.files(src_2010s, pattern = "\\.tif$", full.names = TRUE)
  
  # Extract the 4-digit year from filenames like ..._2014-06.tif
  years_2010s <- as.integer(sub(".*_(\\d{4})-\\d{2}\\.tif$", "\\1", files_2010s))
  
  keep_2010s <- files_2010s[years_2010s >= 2014 & years_2010s <= 2019]
  
  message(var, ": found ", length(keep_2010s),
          " files from 2014-2019 in the 2010-2019 folder (expect 72)")
  
  file.copy(keep_2010s, dest_dir, overwrite = FALSE)
  
  # ---- Copy all of 2020-2024 into the same folder ----
  files_2020s <- list.files(src_2020s, pattern = "\\.tif$", full.names = TRUE)
  
  message(var, ": found ", length(files_2020s),
          " files in the 2020-2024 folder (expect 60)")
  
  file.copy(files_2020s, dest_dir, overwrite = FALSE)
  
  #check to make sure its right
  n_final <- length(list.files(dest_dir, pattern = "\\.tif$"))
  message(var, ": ", n_final, " total files now in ", dest_dir,
          " (expect 132)\n")
}

#####Crop all TIFF files to global extent so data is easier to work with and saves time. 
mask_sf <- vect('data/shapefiles/usa_selected_states.shp')
base_dir <- "data/env_var/monthly_raw"

tif_files <- list.files(base_dir, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)

for (f in tif_files) {
  r <- rast(f)
  
  if (crs(r) != crs(mask_sf)) {
    mask_sf_proj <- project(mask_sf, crs(r))
  } else {
    mask_sf_proj <- mask_sf
  }
  
  r_cropped <- crop(r, mask_sf_proj)
  r_masked  <- mask(r_cropped, mask_sf_proj)
  
  # Temporary filename to safely handle Windows file-locking before overwriting
  temp_out <- paste0(f, "_temp.tif")
  writeRaster(r_masked, filename = temp_out, overwrite = TRUE)
  
  # Remove R pointer to unlock file
  rm(r, r_cropped, r_masked)
  tmp <- gc(verbose = FALSE)
  
  # Replace original file with temporary cropped file
  file.remove(f)
  file.rename(temp_out, f)
}

#####Create bioclim variables for 2014-2024, which will be used to train the models. Then test for correlations. 

glob_bio_2014_2024  <- make_bioclim(years = 2014:2024)
regional_mask <- project(regional_mask, crs(glob_bio_2014_2024))
reg_bio_2014_2024 <- crop(glob_bio_2014_2024, regional_mask)
reg_bio_2014_2024 <- mask(reg_bio_2014_2024, regional_mask)

dist_rivers <- crop(dist_rivers, reg_bio_2014_2024[[1]])
dist_rivers <- mask(dist_rivers, reg_bio_2014_2024[[1]])

reg_bio_2014_2024 <- c(reg_bio_2014_2024, dist_rivers)

plot(reg_bio_2014_2024[[20]])
names(reg_bio_2014_2024)[20] <- "dist_to_river"

#####Do correlation tests:
reg_env <- as.data.frame(reg_bio_2014_2024, na.rm = TRUE)
glob_env <- as.data.frame(glob_bio_2014_2024, na.rm = TRUE)

reg_spearman_matrix <- cor(reg_env, method = "spearman")
glob_spearman_matrix <- cor(glob_env, method = "spearman")

#use variables Bio2: Mean Diurnal Range (Mean of monthly (max temp - min temp)), Bio6: Min Temperature of Coldest Month, Bio12: Annual Precipitation, Dist_to_rivers (only for the regional model)

######Env var to train the models 
#Global
names(glob_bio_2014_2024)
glob_training_var <- glob_bio_2014_2024[[c("bio02", "bio06", "bio12")]]
names(glob_training_var)

names(reg_bio_2014_2024)
reg_training_var <- reg_bio_2014_2024[[c("bio02", "bio06", "bio12", "dist_to_river")]]

######Make all scenarios
#1950-1959, 1960-1969, 1970-1979, 1980-1989, 1990-1999, 2000-2009, 2010-2019, 2020-2024
decade_list <- list(
  "1950_1959" = 1951:1959,
  "1960_1969" = 1960:1969,
  "1970_1979" = 1970:1979,
  "1980_1989" = 1980:1989,
  "1990_1999" = 1990:1999,
  "2000_2009" = 2000:2009,
  "2010_2019" = 2010:2019,
  "2020_2024" = 2020:2024
)

# 2. Container list to store the output SpatRaster objects
bioclim_stacks <- list()

# 3. Loop through each decade
for (decade_label in names(decade_list)) {
  
  years_vec <- decade_list[[decade_label]]
  
  message("--- Processing decade: ", decade_label, " ---")
  
  # Run make_bioclim for the current decade
  bioclim_stacks[[decade_label]] <- make_bioclim(
    years = years_vec,
    bios = c(2, 6, 12),
    mask_poly = regional_mask
  )
}

#now stack all tiffs into terra objects
base_dir <- "data/env_var/bioclim"
names(dist_rivers)
names(dist_rivers) <- "dist_to_river"

# Step 1: Find scenario subfolders
scenario_folders <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
scenario_folders <- scenario_folders[grep("bioclim_", scenario_folders)]

bioclim_scenarios <- list()

# Step 2: Loop, load, align dist_rivers, and assign custom names
for (folder in scenario_folders) {
  
  stack_name <- gsub("-", "_", sub("^bioclim_", "scenario_", basename(folder)))
  tif_files  <- list.files(folder, pattern = "\\.tif$", full.names = TRUE)
  
  if (length(tif_files) > 0) {
    bio_stack <- terra::rast(tif_files)
    
    # Check if geometry matches; if not, resample/crop dist_layer to match bio_stack
    dist_aligned <- dist_rivers
    if (!compareGeom(bio_stack, dist_aligned, stopOnError = FALSE)) {
      dist_aligned <- crop(dist_aligned, bio_stack[[1]])
      dist_aligned <- mask(dist_aligned, bio_stack[[1]])
    }
    
    # Combine bioclim stack + aligned dist_rivers layer
    full_stack <- c(bio_stack, dist_aligned)
    
    # Save using formatted scenario name
    bioclim_scenarios[[stack_name]] <- full_stack
    
    message("Loaded ", nlyr(full_stack), " layers into: ", stack_name)
    
  } else {
    warning("No .tif files found in: ", folder)
  }
}

#I have no idea why 2020-2024 did not work on the loop above, but i did it manually below
bio2_2020_2024 <- rast("data/env_var/bioclim/bioclim_2020-2024/bio02.tif")
bio6_2020_2024 <- rast("data/env_var/bioclim/bioclim_2020-2024/bio06.tif")
bio12_2020_2024 <- rast("data/env_var/bioclim/bioclim_2020-2024/bio12.tif")

bio2_2020_2024 <-crop(bio2_2020_2024, reg_bio_2014_2024[[1]])
bio2_2020_2024 <- mask(bio2_2020_2024, reg_bio_2014_2024[[1]])
bio6_2020_2024 <-crop(bio6_2020_2024, reg_bio_2014_2024[[1]])
bio6_2020_2024 <- mask(bio6_2020_2024, reg_bio_2014_2024[[1]])
bio12_2020_2024 <-crop(bio12_2020_2024, reg_bio_2014_2024[[1]])
bio12_2020_2024 <- mask(bio12_2020_2024, reg_bio_2014_2024[[1]])

scenario_2020_2024 <- c(bio2_2020_2024, bio6_2020_2024, bio12_2020_2024, dist_rivers)

#make folder for scenarios
dir.create("data/env_var/scenarios/", recursive = TRUE, showWarnings = FALSE)

writeRaster(bioclim_scenarios[["scenario_1951_1959"]], "data/env_var/scenarios/1951_1959.tif")
writeRaster(bioclim_scenarios[["scenario_1960_1969"]], "data/env_var/scenarios/1960_1969.tif")
writeRaster(bioclim_scenarios[["scenario_1970_1979"]], "data/env_var/scenarios/1970_1979.tif")
writeRaster(bioclim_scenarios[["scenario_1980_1989"]], "data/env_var/scenarios/1980_1989.tif")
writeRaster(bioclim_scenarios[["scenario_1990_1999"]], "data/env_var/scenarios/1990_1999.tif")
writeRaster(bioclim_scenarios[["scenario_2000_2009"]], "data/env_var/scenarios/2000_2009.tif")
writeRaster(bioclim_scenarios[["scenario_2010_2019"]], "data/env_var/scenarios/2010_2019.tif")
writeRaster(scenario_2020_2024, "data/env_var/scenarios/2020_2024.tif")
writeRaster(glob_training_var, "data/env_var/scenarios/glob_training_var.tif")
writeRaster(reg_training_var, "data/env_var/scenarios/reg_training_var.tif")

#####Try to see if you can visualize differences between scenarios, lots of this below is scrap
library(dplyr)
library(tidyr)
library(ggplot2)

# Convert all scenario stacks into a single long-format data frame
scenario_df <- lapply(names(bioclim_scenarios), function(scen_name) {
  df <- as.data.frame(bioclim_scenarios[[scen_name]], xy = TRUE)
  df$scenario <- scen_name
  return(df)
}) %>% bind_rows()

# Clean up scenario labels for nice plot facets (e.g., "1951-1959")
scenario_df <- scenario_df %>%
  mutate(decade = gsub("scenario_", "", scenario),
         decade = gsub("_", "-", decade))

ggplot(scenario_df, aes(x = x, y = y, fill = bio12)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Bio12 (mm)", option = "magma") +
  facet_wrap(~ decade, ncol = 4) +
  coord_equal() +
  theme_minimal() +
  labs(title = "Annual Precipitation (Bio12) Trajectory Across Decades",
       x = "Longitude", y = "Latitude") +
  theme(strip.text = element_text(face = "bold", size = 11))

ggplot(scenario_df, aes(x = x, y = y, fill = bio06)) +
  geom_raster() +
  scale_fill_viridis_c(name = "minimum temp", option = "magma") +
  facet_wrap(~ decade, ncol = 4) +
  coord_equal() +
  theme_minimal() +
  labs(title = "Minimum temp of coldest month (bio6) Trajectory Across Decades",
       x = "Longitude", y = "Latitude") +
  theme(strip.text = element_text(face = "bold", size = 11))


long_df <- scenario_df %>%
  pivot_longer(cols = c("bio02", "bio06", "bio12", "dist_to_river"),
               names_to = "variable", 
               values_to = "value")

# Boxplots showing variable shifts across decades
ggplot(long_df, aes(x = decade, y = value, fill = decade)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  facet_wrap(~ variable, scales = "free_y", ncol = 2) +
  scale_fill_viridis_d(option = "mako", guide = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Environmental Variable Distribution Across Decades",
       x = "Decade Scenario", y = "Value")


# Subtract baseline raster from current raster
delta_bio06 <- bioclim_scenarios[["scenario_2010_2019"]][["bio02"]] - 
  bioclim_scenarios[["scenario_1951_1959"]][["bio02"]]

# Plot difference (Diverging color scale: Red = Decrease, Blue = Increase)
plot(delta_bio02, 
     col = hcl.colors(100, "PuOr"), 
     main = "Bio06 Delta: (2010-2019) vs (1951-1959)")
