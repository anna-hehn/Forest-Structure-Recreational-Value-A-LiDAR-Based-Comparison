# Forest-Structure-Recreational-Value-A-LiDAR-Based-Comparison
A standardized LiDAR workflow in R analyzing and comparing forest structure parameters between a designated healing forest ('Wald für die Seele', Bad Kissingen) and a managed forest ('University Forest Sailershausen') using the lidR and TreeCompR packages.

## Study Motivation
The phrase "forests are good for you" ("Wald tut gut")  reflects a widely held belief. Forests serve as essential environments for health and recreation, allowing individuals to escape everyday stress, noise, and urban hectic. However, the growing demand for recreational spaces introduces new expectations for forest management (Wagner et al. 2022).
### This study adresses the folowing questions:
1. How do recreation-relevant forest structure parameters derived from airborne LiDAR data (ALS) differ between a designated recreational forest and a production-oriented commercial forest?
2. To what extent do these different forest types offer accessible, barrier-free terrain?

## Study Area
The Area of Interest (AOI) was manually defined via point cloud segmentation within the open-source software 'CloudCompare'. The AOI is divided into two distinct forest plots, each covering an extent of approximately 20 hectares.

### Forest Plots
#### Wald für die Seele, Bad Kissingen
The "Wald für die Seele" (Forest for the Soul) is a 14.5 hectare nature experience and art project located within the Klauswald in Bad Kissingen, managed by the Foundation for Consciousness Sciences since 2015. The project aims to demonstrate how nature and biodiversity promote mental, physical, and social well-being. It serves as a dedicated space for ecotherapy, meditation, reflection, and the experience of art in nature (Wald für die Seele 2023; Galuska 2023).

#### University Forest, Sailershausen
The "University Forest Sailershausen" covers approximately 2,300 hectares near Haßfurt and is managed by the forestry office of the University of Würzburg. Operating as a sustainable commercial forest, its tree composition is dominated by deciduous species. In addition to fulfilling ecological functions like a large bird sanctuary, the woodland serves a dual purpose as a production forest and a scientific research area for academic institutions (Schmidt/Täufer 2023; DFV 2017).

## Data Source
The raw datasets utilized in this study were obtained from the OpenData portal of the Bavarian Agency for Digitisation (Bayerische Vermessungsverwaltung (BVV)) and consist of Airborne Laser Scanning (ALS) point clouds. The flight campaign for the University Forest Sailershausen took place between January 18 and February 1, 2025, while the "Wald für die Seele" was surveyed between January 10 and January 28, 2024. 
Both datasets feature a minimum point density of 4 points/m². The vertical accuracy is approximately 0.12 m in flat, open terrain, while the horizontal positioning accuracy is around 0.30 m (BVV 2026).

---

# Workflow

## Derivation of Forest Structure Parameters using the [lidR](https://github.com/r-lidar/lidR) package

### Environment Setup & Data Loading

```R
library(lidR)
library(terra)  
library(sf)     
library(dplyr)
library(ggplot2)

# Set working directory and load segmented point cloud
setwd("D:/LiDen/Daten")
las <- readLAS("603_5547_seg.laz", filter = "-drop_withheld")

# Assign Coordinate Reference System (CRS: ETRS89 / UTM zone 32N)
crs(las) <- 25832

# Set local export path and inspect raw point cloud
setwd("D:/LiDen/Daten/Sailershausen/lidR_Sailershausen")
plot(las)
```

### DTM Generation, Canopy Modeling and Tree Segmentation

#### DTM Creation & Height Normalization
A Digital Terrain Model (DTM) is generated at a 0.2 m resolution using a k-nearest neighbor interpolation with inverse distance weighting (knnidw). This DTM is used to normalize the point cloud, converting absolute elevations into true tree heights above the ground.

```R
## DTM Creation & Height Normalization
dtm <- rasterize_terrain(las, res = 0.2, algorithm = knnidw())
las_norm <- normalize_height(las, dtm)
```

#### Treetop Detection
Individual treetops are located using a Local Maximum Filter (lmf) with a dynamic window size defined by a custom height function to adjust for varying crown widths.

```R
## Treetop Detection
f <- function(x) { x * 0.1 + 3 }
ttops <- locate_trees(las = las_norm, algorithm = lmf(f), uniqueness = "bitmerge")
plot(ttops["Z"], cex = 0.5, pch = 19, pal = height.colors, nbreaks = 30, main = "Detected Treetops")
```

#### Canopy Height Model (CHM)
A pit-free Canopy Height Model (CHM) is rasterized at a 0.5 m resolution and subsequently smoothed using a 3x3 median filter (terra::focal) to eliminate spatial artifacts and data pits.

```R
## Canopy Height Model (CHM)
chm <- rasterize_canopy(las = las_norm, res = 0.5, algorithm = pitfree(subcircle = 0.15))
# Smooth CHM to reduce artifacts
kernel <- matrix(1, 3, 3)
chm_smooth <- terra::focal(chm, w = kernel, fun = median, na.rm = TRUE)
plot(chm_smooth, main = "Smoothed Canopy Height Model")
```

#### Tree Segmentation and Crown Metrics
* **Tree Segmentation:** Individual trees are isolated using the `silva2016` algorithm based on the smoothed CHM and detected treetops.
* **Metrics Extraction:** Standard forest metrics are derived via `crown_metrics()` using two geometric approaches:
  * `geom = "point"` → Extracts precise spatial locations and individual tree heights ($Z$).
  * `geom = "convex"` → Computes crown areas using calculated convex hulls.

```R
## Tree Segmentation & Crown Metrics
# Segment trees using Silva et al. (2016) algorithm
algo <- silva2016(chm_smooth, ttops)
las_seg <- segment_trees(las_norm, algo)

## Derive crown metrics (Convex Hull for area and Points for locations)
crowns <- crown_metrics(las_seg, func = .stdtreemetrics, geom = "convex")
trees_lidar <- crown_metrics(las_seg, func = .stdtreemetrics, geom = "point")
plot(crowns["convhull_area"], main = "Crown area (convex hull)")
plot(trees_lidar["Z"], main = "Tree heights", pch = 16)
summary(crowns$convhull_area)
```

### Inventory Data Export for TreeCompare

```R
## Export Inventory Data for TreeCompare
inventory_table <- as.data.frame(sf::st_coordinates(trees_lidar))
names(inventory_table) <- c("x", "y")
inventory_table$hoehe <- trees_lidar$Z
inventory_table$kronenflaeche <- trees_lidar$convhull_area

## Save dataset
write.csv2(inventory_table, "Baumparameter_603_5547.csv", row.names = FALSE)
```

---

### Density Metrics: Stem Density and Canopy Cover
* **Stem Density ($N/\text{ha}$):** Treetop point geometries are rasterized onto a fine 1 m grid to log tree presence. 
  * Aggregation → Pixels are grouped into a **20x20 m forest inventory grid** (400 $\text{m}^2$).
  * Hectare Extrapolation → To scale the local tree count up to a standard hectare (10,000 $\text{m}^2$), values are multiplied by a scaling factor of 25 ($10,000 / 400 = 25$).
* **Canopy Cover (%):** The smoothed Canopy Height Model (CHM) is thresholded at 2 meters to isolate the upper canopy layer from low vegetation and ground noise.
  * Aggregation → The binary mask is aggregated to the same 20x20 m grid to compute the percentage of crown coverage.
* **Geospatial Export:** The resulting stem density matrix is exported as a GeoTIFF raster (`stammzahl_dichte_sail.tif`) for spatial mapping and statistical comparison across the plots.

```R
## Stem Density (Trees per Hectare)
# Rasterize treetop locations at 1m resolution
fine_raster <- terra::rast(ttops, res = 1)
stem_presence <- terra::rasterize(ttops, fine_raster, fun = "length", background = 0)

# Aggregate to 20x20m grid cells (400 sqm)
stem_sum_20m <- terra::aggregate(stem_presence, fact = 20, fun = sum, na.rm = TRUE)

# Extrapolate to 1 hectare scale (factor: 10000 / 400 = 25)
stammzahl_dichte <- stem_sum_20m * 25

## Canopy Cover (%)
# Threshold CHM at 2 meters height and aggregate to 20x20m grid cells
crown_mask <- chm_smooth > 2
kronenschluss <- terra::aggregate(crown_mask, fact = 20, fun = mean, na.rm = TRUE) * 100

# Plot Canopy Cover
plot(kronenschluss, 
     main = "Canopy Cover (%)", 
     col = colorRampPalette(c("lightgoldenrod1", "yellow4", "darkgreen"))(255),
     range = c(0, 100),
     plg = list(title = "Cover (%)"))

## Export stem density raster
terra::writeRaster(stammzahl_dichte, "D:/LiDen/Daten/Sailershausen/lidR_Sailershausen/stammzahl_dichte_sail.tif", overwrite = TRUE)
```
---

## Terrain Analysis: Slope and Accessibility
