# Unsupervised k-means classification of Sentinel-2 imagery for Eyre Bird
# Observatory, Western Australia (32°14'47"S 126°18'06"E).
#
# The observatory sits on sandy coastal scrub with patches of mallee
# woodland, which true-colour or plain NIR/Red/Blue bands don't separate
# well since both cover types are green vegetation. Instead this stacks:
#   - B04  Red        (10m) baseline reflectance / NDVI input
#   - B08  NIR         (10m) vegetation vigour and canopy density
#   - B11  SWIR1       (20m) sand/soil exposure vs. vegetation moisture
#   - B05  Red Edge 1  (20m) canopy structure / chlorophyll
#   - B06  Red Edge 2  (20m) canopy structure / chlorophyll
# The 20m bands are resampled to the 10m grid before stacking. B03 (Green)
# is loaded only to compute MNDWI so ocean pixels can be masked out before
# clustering, rather than let a water class eat one of the 6 clusters.
#
# The stack is clipped to a 1km buffer around the observatory and
# classified into 6 clusters with k-means. Results are plotted with
# ggplot2 and saved as a PNG.
#
# Adjust file paths and band file names below for your dataset.

#install.packages(c("terra", "ggplot2"))
library(terra)
library(ggplot2)

# set working directory to where your Sentinel-2 bands are stored
setwd("/home/si/gdrive/landsat/Working")

# ensure you are working from the correct directory
print(getwd())

# --- 1. SETUP PATHS ---
# Replace these with your actual file paths
s2_red_path <- "2026-05-18-00:00_2026-05-18-23:59_Sentinel-2_L2A_B04_(Raw).tiff" # Red
s2_nir_path <- "2026-05-18-00:00_2026-05-18-23:59_Sentinel-2_L2A_B08_(Raw).tiff" # NIR
s2_sw1_path <- "2026-05-18-00:00_2026-05-18-23:59_Sentinel-2_L2A_B11_(Raw).tiff" # SWIR1
s2_re1_path <- "2026-05-18-00:00_2026-05-18-23:59_Sentinel-2_L2A_B05_(Raw).tiff" # Red Edge 1
s2_re2_path <- "2026-05-18-00:00_2026-05-18-23:59_Sentinel-2_L2A_B06_(Raw).tiff" # Red Edge 2
s2_grn_path <- "2026-05-18-00:00_2026-05-18-23:59_Sentinel-2_L2A_B03_(Raw).tiff" # Green (masking only)

# --- 2. LOAD BANDS ---
s2_red <- rast(s2_red_path)
s2_nir <- rast(s2_nir_path)
s2_sw1 <- rast(s2_sw1_path)
s2_re1 <- rast(s2_re1_path)
s2_re2 <- rast(s2_re2_path)
s2_grn <- rast(s2_grn_path)

print("Sentinel-2 bands loaded.")

# --- 3. CLIP TO EYRE BIRD OBSERVATORY (1KM ZONE) ---
# Coordinates: 32°14'47"S 126°18'06"E
eyre_coords <- matrix(c(126.3017, -32.2464), ncol = 2)
eyre_point <- vect(eyre_coords, crs = "EPSG:4326")

# Project point to match imagery CRS, then buffer by 1km
eyre_proj <- project(eyre_point, crs(s2_nir))
eyre_buffer <- buffer(eyre_proj, width = 1000)

# Crop each band to the buffer first (cheap), then resample the 20m bands
red_c <- crop(s2_red, eyre_buffer, mask = TRUE)
nir_c <- crop(s2_nir, eyre_buffer, mask = TRUE)
grn_c <- crop(s2_grn, eyre_buffer, mask = TRUE)
sw1_c <- crop(s2_sw1, eyre_buffer, mask = TRUE)
re1_c <- crop(s2_re1, eyre_buffer, mask = TRUE)
re2_c <- crop(s2_re2, eyre_buffer, mask = TRUE)

# --- 4. RESAMPLE 20M BANDS TO THE 10M GRID ---
sw1_10m <- resample(sw1_c, nir_c, method = "bilinear")
re1_10m <- resample(re1_c, nir_c, method = "bilinear")
re2_10m <- resample(re2_c, nir_c, method = "bilinear")

# --- 5. BUILD RED / NIR / SWIR1 / RED-EDGE COMPOSITE ---
composite <- c(red_c, nir_c, sw1_10m, re1_10m, re2_10m)
names(composite) <- c("Red", "NIR", "SWIR1", "RedEdge1", "RedEdge2")

# --- 6. MASK OUT THE OCEAN ---
# MNDWI (Xu, 2006) = (Green - SWIR1) / (Green + SWIR1); water is typically > 0.
mndwi <- (grn_c - sw1_10m) / (grn_c + sw1_10m)
water_threshold <- 0
land_mask <- ifel(mndwi > water_threshold, NA, 1)

img_clipped <- mask(composite, land_mask)

# --- 7. K-MEANS CLASSIFICATION ---
# Convert to dataframe, keeping cell numbers, dropping NAs (ocean + outside the 1km circle)
img_df <- as.data.frame(img_clipped, cells = TRUE, na.rm = TRUE)

# Run k-means with 6 clusters
set.seed(123) # for reproducibility
k_results <- kmeans(img_df[, -1], centers = 6, iter.max = 100, nstart = 10)

# Create a blank raster and fill with cluster results
km_raster <- rast(img_clipped, nlyr = 1)
names(km_raster) <- "cluster"
km_raster[img_df$cell] <- k_results$cluster

print("kmeans completed.")

# --- 8. PLOT RESULTS WITH GGPLOT2 ---
num_clusters <- 6

km_df <- as.data.frame(km_raster, xy = TRUE, na.rm = TRUE)
km_df$cluster <- factor(km_df$cluster)

p <- ggplot(km_df, aes(x = x, y = y, fill = cluster)) +
  geom_raster() +
  scale_fill_brewer(palette = "Set2", name = "Class") +
  coord_equal() +
  labs(
    title = "Sentinel-2 6-Class K-means: Eyre (Air) Bird Observatory (1km, ocean masked)",
    x = NULL, y = NULL
  ) +
  theme_minimal()

print(p)

# Save the plot to a PNG file
png_filename <- "Sentinel-2-airbird-kmeans.png"
ggsave(png_filename, plot = p, width = 6, height = 6, units = "in", dpi = 300)

print(paste("Saved plot to", png_filename))
