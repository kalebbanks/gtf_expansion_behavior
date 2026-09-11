#####################
#Script: 01_organizing_occurences.R
#Purpose: organize and thin occurences for modeling
#Inputs: data/occ
#Outputs: 
#Author: Kaleb M. Banks
#Date: 2026-
#####################

########8/25,NOTE:ONCE MODEL FRAMEWORK DONE, USE NEWER GBIF DATASET#########


######Packages
install.packages("spThin")
library(dplyr)
library(ggplot2)
library(maps)
library(terra)
library(sf)
library(SDMtune)
library(ecospat)
library(spThin)

######Functions

######Load survey data, GBIF, and shapefile
raw_surveys <- read.csv("data/occ/raw_surveys.csv")
raw_surveys <- raw_surveys %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y"))
GBIF <- read.csv("data/occ/GBIF_9_2026_dirty.csv")
expanded_range <- st_read("data/shapefiles/gtf_expanded_range.shp")
reg_low_sample_range <- st_read("data/shapefiles/ok_low_sample_range.shp")

#####load rasters
reg_training_var <- rast("data/env_var/scenarios/1951_1959.tif")
reg_training_var_1 <- reg_training_var[[1]]
glob_training_var <- rast("data/env_var/scenarios/glob_training_var.tif")
glob_training_var_1 <- glob_training_var[[1]]

#####GLOBAL OCCURRENCE DATASET
GBIF <- GBIF[!is.na(GBIF$decimalLatitude) & !is.na(GBIF$decimalLongitude),]
GBIF <- GBIF[!is.na(GBIF$decimalLatitude) & !is.na(GBIF$decimalLongitude) &
               !is.na(GBIF$year) & GBIF$year >= 2014 & GBIF$year <= 2024, ]
GBIF <- GBIF[is.na(GBIF$coordinateUncertaintyInMeters) | 
               GBIF$coordinateUncertaintyInMeters <= 20000, ]


#Many "hitch hiker" observations, need to get rid of those with a custom shapefile i made in arc gis
expanded_range <- st_transform(expanded_range, crs = 4326)
ggplot() +
  geom_polygon(data = us_states, aes(x = long, y = lat, group = group),
               fill = "grey95", color = "grey50", linewidth = 0.2) +
  geom_sf(data = expanded_range, fill = "orange", alpha = 0.2, color = "darkorange", 
          linewidth = 0.6, inherit.aes = FALSE) +
  geom_point(data = GBIF, aes(x = decimalLongitude, y = decimalLatitude),
             color = "forestgreen", size = 1.5, alpha = 0.6) +
  coord_sf() +
  labs(title = "Hyla cinerea occurrences (2014–2024)",
       x = "Longitude", y = "Latitude") +
  theme_minimal()

GBIF_sf <- st_as_sf(GBIF, coords = c("decimalLongitude", "decimalLatitude"), 
                    crs = 4326, remove = FALSE)

inside <- st_intersects(GBIF_sf, expanded_range, sparse = FALSE)[, 1]
GBIF <- GBIF[inside, ]

glob_occ <- setNames(GBIF[, c("decimalLongitude", "decimalLatitude")], c("lon", "lat"))

#####REGIONAL OCCURRENCE DATASETS


#####Convert raw survey data to presence and absence list
#A few sites were surveyed more than once, we need to see which ones 
raw_surveys$site[duplicated(raw_surveys$site)]

raw_surveys[raw_surveys$lon > 0, ]

#This thins out sites that have been surveyed twice, retaining the one with the highest h.cinerea call index value. Discards surveys where we visited a site, didn't hear anything, then surveyed it again and heard stuff. 
raw_surveys_clean <- raw_surveys %>%
  group_by(site) %>%
  group_modify(~ {
    if (nrow(.x) == 1) return(.x)  # not a duplicate, leave alone
    
    max_val <- max(.x$h.cinerea)
    tied <- .x[.x$h.cinerea == max_val, ]
    
    if (nrow(tied) == 1) {
      tied
    } else {
      # tied on max value — keep the most recent one
      tied[which.max(tied$date), ]
    }
  }) %>%
  ungroup()

raw_surveys_clean <- raw_surveys_clean %>%
  mutate(presence = if_else(h.cinerea > 0, 1, 0)) %>%
  select(site, lat, lon, date, presence)

#divide into survey occ and survey abs
survey_occ <- raw_surveys_clean[!is.na(raw_surveys_clean$presence) & raw_surveys_clean$presence == 1, ]

survey_abs <- raw_surveys_clean[!is.na(raw_surveys_clean$presence) & raw_surveys_clean$presence == 0, ]


######Filter GBIF data to only include records from 2014-2024
ok_GBIF <- GBIF[!is.na(GBIF$stateProvince) & GBIF$stateProvince == "Oklahoma", ]

ok_GBIF <- ok_GBIF[!is.na(ok_GBIF$year) & ok_GBIF$year >= 2014 & ok_GBIF$year <= 2026, ]

ok_gbif_occ <- ok_GBIF[, c("decimalLatitude", "decimalLongitude")]
colnames(ok_gbif_occ) <- c("lat", "lon")

#####Combine Regional occurrence datasets
coords1 <- survey_occ[, c("lat", "lon")]
coords1$source <- "surveys"

coords2 <- ok_gbif_occ[, c("lat", "lon")]
coords2$source <- "GBIF"

reg_occ <- rbind(coords1, coords2)
reg_abs <- survey_abs[, c("lat", "lon")]

reg_low_sample_range <- st_transform(reg_low_sample_range, crs = 4326)
ok_state <- subset(us_states, region == "oklahoma")
ggplot() +
  geom_polygon(data = ok_state, aes(x = long, y = lat, group = group),
               fill = "grey95", color = "grey50", linewidth = 0.2) +
  geom_sf(data = reg_low_sample_range, fill = "orange", alpha = 0.2, color = "darkorange", 
          linewidth = 0.6, inherit.aes = FALSE) +
  geom_point(data = reg_occ, aes(x = lon, y = lat, color = source),
             size = 1.5, alpha = 0.6) +
  geom_point(data = reg_abs, aes(x = lon, y = lat),
             color = "black", shape = 4, size = 1.5, alpha = 0.6) +
  coord_sf() +
  labs(title = "Regional occurrences and absences (Oklahoma)",
       x = "Longitude", y = "Latitude", color = "Occurrence source") +
  theme_minimal()


#### Thin regional occurrences #####

reg_occ_sf <- st_as_sf(reg_occ, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
in_low_sample <- st_intersects(reg_occ_sf, reg_low_sample_range, sparse = FALSE)[, 1]

occ_low_sample  <- reg_occ[in_low_sample, ]
occ_high_sample <- reg_occ[!in_low_sample, ]

# Step 1: remove duplicate/near-duplicate points in the same raster cell (e.g. state park hotspot)
cell_thin <- function(df) {
  as.data.frame(thinData(
    coords = as.data.frame(df[, c("lon", "lat")]),
    env = reg_training_var_1, x = "lon", y = "lat"
  ))
}

occ_low_sample_cellthin  <- cell_thin(occ_low_sample)
occ_high_sample_cellthin <- cell_thin(occ_high_sample)

# Step 2: nearest-neighbor spacing (km), now meaningful, used to pick a thinning distance
get_nn_dist <- function(df, lon_col = "lon", lat_col = "lat") {
  pts_proj <- st_transform(st_as_sf(df, coords = c(lon_col, lat_col), crs = 4326), 5070)
  d <- st_distance(pts_proj); diag(d) <- NA
  apply(d, 1, min, na.rm = TRUE) / 1000
}

summary(get_nn_dist(occ_low_sample_cellthin))
summary(get_nn_dist(occ_high_sample_cellthin))

# Step 3: pick thin.par based on the summaries above, then thin to a shared distance
thin_occ <- function(df, dist_km) {
  df$species <- "Hcinerea"
  thin(df, lat.col = "lat", long.col = "lon", spec.col = "species",
       thin.par = dist_km, reps = 10,
       locs.thinned.list.return = TRUE, write.files = FALSE)[[1]]
}

occ_low_sample_thinned  <- thin_occ(occ_low_sample_cellthin, dist_km = 5)
occ_high_sample_thinned <- thin_occ(occ_high_sample_cellthin, dist_km = 5)

reg_occ_thinned <- rbind(occ_low_sample_thinned, occ_high_sample_thinned)
names(reg_occ_thinned) <- c("X", "Y")
##### Plot #####
ok_state <- subset(us_states, region == "oklahoma")

ggplot() +
  geom_polygon(data = ok_state, aes(x = long, y = lat, group = group),
               fill = "grey95", color = "grey50", linewidth = 0.2) +
  geom_point(data = reg_abs_thinned, aes(x = Longitude, y = Latitude),
             size = 1.5, alpha = 0.6) +
  coord_sf() +
  labs(title = "Regional occurrences (thinned, Oklahoma)",
       x = "Longitude", y = "Latitude") +
  theme_minimal()

#####Regional absences 
abs_sf <- st_as_sf(reg_abs, coords = c("lon","lat"), crs = 4326, remove = FALSE)
summary(get_nn_dist(reg_abs))


reg_abs$species <- "Hcinerea"
reg_abs_thinned <- thin(reg_abs, lat.col = "lat", long.col = "lon", spec.col = "species",
                        thin.par = 5, reps = 10,
                        locs.thinned.list.return = TRUE, write.files = FALSE)[[1]]

nrow(reg_abs)          # before
nrow(reg_abs_thinned)  # after
names(reg_abs_thinned) <- c("X", "Y")


#####Thin Global occurences 
glob_occ_thinned <- as.data.frame(thinData(
  coords = as.data.frame(glob_occ[, c("lon", "lat")]),
  env = glob_training_var,   # your global-scale terra SpatRaster
  x = "lon", y = "lat"
))

nrow(glob_occ)          # before
nrow(glob_occ_thinned)  # after
plot(glob_occ)
plot(glob_occ_thinned)
names(glob_occ_thinned) <- c("X", "Y")

#Generate global absences 
glob_occ_sf <- st_as_sf(glob_occ_thinned, coords = c("lon", "lat"), crs = 4326)

# Project to a metric CRS for accurate buffering (Albers Equal Area, CONUS)
glob_occ_proj <- st_transform(glob_occ_sf, crs = 5070)

# Buffer each point by 300 km, then dissolve into one polygon
buffer_300km <- st_buffer(glob_occ_proj, dist = 300000)  # 300 km in meters
buffer_300km <- st_union(buffer_300km)

# Reproject back to lat/lon for consistency with your rasters/plots
buffer_300km_ll <- st_transform(buffer_300km, crs = 4326)
plot(glob_training_var[[1]])
plot(buffer_300km_ll)
plot(st_geometry(buffer_300km_ll))
points(glob_occ_thinned$lon, glob_occ_thinned$lat, pch = 20, col = "forestgreen", cex = 0.3)

# 1. Mask the global climate raster to the 300 km accessible area
glob_var_masked <- mask(crop(glob_training_var, vect(buffer_300km_ll)), vect(buffer_300km_ll))
plot(glob_var_masked)

# 2. Build a small exclusion buffer around each presence point (e.g. 10 km)
#    to prevent background points from landing right on top of known occurrences
occ_sf <- st_as_sf(glob_occ_thinned, coords = c("lon","lat"), crs = 4326)
exclusion_buffer <- st_buffer(occ_sf, dist = 10000)   # 10 km, in a projected CRS ideally — see note below
exclusion_union <- st_union(exclusion_buffer)

# 3. Mask out the exclusion zone from the accessible area raster
glob_var_for_bg <- mask(glob_var_masked, vect(exclusion_union), inverse = TRUE)
plot(glob_var_for_bg)

library(terra)

set.seed(123)  # makes the random PA selection reproducible

# Find cells in glob_var_for_bg that contain valid data
valid_cells <- which(!is.na(values(glob_var_for_bg[[1]])))

# Make sure there are enough valid cells
length(valid_cells)

# Randomly select 10,000 unique raster cells
pa_cells <- sample(valid_cells, size = 10000, replace = FALSE)

# Convert selected cells to point coordinates
PA <- xyFromCell(glob_var_for_bg, pa_cells)

# Convert to a data frame
PA <- as.data.frame(PA)
names(PA) <- c("Longitude", "Latitude")

# Check
head(PA)
nrow(PA)

plot(PA)

names(PA) <- c("X", "Y")

write.csv(reg_abs_thinned, "data/occ/reg_abs.csv")
write.csv(reg_occ_thinned, "data/occ/reg_occ.csv")
write.csv(glob_occ_thinned, "data/occ/glob_occ.csv")
write.csv(PA, "data/occ/glob_abs.csv")





