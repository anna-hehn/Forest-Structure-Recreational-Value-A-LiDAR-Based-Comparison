##STATISTSIC COMPACT: Forest Structure Parameters

library(lidR)
library(terra)
library(sf)

#Load Data
setwd("D:/LiDen/Daten")

las <- readLAS(
  "603_5547_seg.laz",
  filter = "-drop_withheld"
)
crs(las) <- 25832

#Create DTM (20 cm resolution) and normalize tree heights 
dtm <- rasterize_terrain(
  las,
  res = 0.2,
  algorithm = knnidw()
)

las_norm <- normalize_height(las, dtm)

#Detect TreeTops
f <- function(x) {
  x * 0.1 + 3
}

ttops <- locate_trees(
  las_norm,
  algorithm = lmf(f),
  uniqueness = "bitmerge"
)


#Create and smooth Canopy Height Model (CHM)
chm <- rasterize_canopy(
  las_norm,
  res = 0.5,
  algorithm = pitfree(subcircle = 0.15)
)

chm_smooth <- terra::focal(
  chm,
  w = matrix(1, 3, 3),
  fun = median,
  na.rm = TRUE
)


#Segment individual trees and calculate crown metrics
las_seg <- segment_trees(
  las_norm,
  silva2016(chm_smooth, ttops)
)

#Calculate standard tree metrics and crown area
trees <- crown_metrics(
  las_seg,
  func = .stdtreemetrics,
  geom = "convex"
)


#Extract tree height and crown area
inventory <- data.frame(
  height = trees$Z,
  crown_area = trees$convhull_area
)


#Calculate stem density (trees per hectare)
fine_raster <- terra::rast(
  ttops,
  res = 1
)

stem_presence <- terra::rasterize(
  ttops,
  fine_raster,
  fun = "length",
  background = 0
)

stem_sum_20m <- terra::aggregate(
  stem_presence,
  fact = 20,
  fun = sum,
  na.rm = TRUE
)

# Convert trees per 400 m² into trees per hectare
# 1 hectare = 10,000 m²
# 10,000 / 400 = 25
stem_density <- stem_sum_20m * 25


#Calculate canopy cover
crown_mask <- chm_smooth > 2

# Calculate the percentage of canopy cover within 20 x 20 m cells
canopy_cover <- terra::aggregate(
  crown_mask,
  fact = 20,
  fun = mean,
  na.rm = TRUE
) * 100


#OVERVIEW!
stats_summary <- data.frame(
  Metric = c(
    "Mean Tree Height",
    "Median Tree Height",
    "Maximum Tree Height",
    "Mean Crown Area",
    "Median Crown Area",
    "Mean Stem Density (trees/ha)",
    "Mean Canopy Cover (%)"
  ),
  Value = c(
    mean(inventory$height, na.rm = TRUE),
    median(inventory$height, na.rm = TRUE),
    max(inventory$height, na.rm = TRUE),
    mean(inventory$crown_area, na.rm = TRUE),
    median(inventory$crown_area, na.rm = TRUE),
    mean(terra::values(stem_density), na.rm = TRUE),
    mean(terra::values(canopy_cover), na.rm = TRUE)
  )
)

print(stats_summary)