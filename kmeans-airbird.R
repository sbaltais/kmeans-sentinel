# Unsupervised k-means classification of Sentinel-2 imagery for Eyre Bird
# Observatory, Western Australia (32°14'47"S 126°18'06"E).
#
# Stacks the Infrared (B08), Red (B04), and Blue (B02) bands into a
# composite raster, clips it to a 1km buffer around the observatory, and
# classifies the clipped image into 6 clusters with k-means. Results are
# plotted with ggplot2 and saved as a PNG.
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
s2_nir_path <- "2026-05-18-00:00_2026-05-18-23:59_Sentinel-2_L2A_B08_(Raw).tiff" # Sentinel-2 Infrared (NIR)
s2_red_path <- "2026-05-18-00:00_2026-05-18-23:59_Sentinel-2_L2A_B04_(Raw).tiff" # Sentinel-2 Red
s2_blu_path <- "2026-05-18-00:00_2026-05-18-23:59_Sentinel-2_L2A_B02_(Raw).tiff" # Sentinel-2 Blue

# --- 2. LOAD BANDS ---
s2_nir <- rast(s2_nir_path)
s2_red <- rast(s2_red_path)
s2_blu <- rast(s2_blu_path)

print("Sentinel-2 bands loaded.")

# --- 3. CLIP TO EYRE BIRD OBSERVATORY (1KM ZONE) ---
# Coordinates: 32°14'47"S 126°18'06"E
eyre_coords <- matrix(c(126.3017, -32.2464), ncol = 2)
eyre_point <- vect(eyre_coords, crs = "EPSG:4326")

# Project point to match imagery CRS, then buffer by 1km
eyre_proj <- project(eyre_point, crs(s2_nir))
eyre_buffer <- buffer(eyre_proj, width = 1000)

# --- 4. BUILD INFRARED / RED / BLUE COMPOSITE AND CLIP ---
composite <- c(s2_nir, s2_red, s2_blu)
names(composite) <- c("Infrared", "Red", "Blue")

img_clipped <- crop(composite, eyre_buffer, mask = TRUE)

# --- 5. K-MEANS CLASSIFICATION ---
# Convert to dataframe, keeping cell numbers, dropping NAs (pixels outside the 1km circle)
img_df <- as.data.frame(img_clipped, cells = TRUE, na.rm = TRUE)

# Run k-means with 6 clusters
set.seed(123) # for reproducibility
k_results <- kmeans(img_df[, -1], centers = 6, iter.max = 100, nstart = 10)

# Create a blank raster and fill with cluster results
km_raster <- rast(img_clipped, nlyr = 1)
names(km_raster) <- "cluster"
km_raster[img_df$cell] <- k_results$cluster

print("kmeans completed.")

# --- 6. PLOT RESULTS WITH GGPLOT2 ---
num_clusters <- 6

km_df <- as.data.frame(km_raster, xy = TRUE, na.rm = TRUE)
km_df$cluster <- factor(km_df$cluster)

p <- ggplot(km_df, aes(x = x, y = y, fill = cluster)) +
  geom_raster() +
  scale_fill_brewer(palette = "Set2", name = "Class") +
  coord_equal() +
  labs(
    title = "Sentinel-2 6-Class K-means: Eyre (Air) Bird Observatory (1km)",
    x = NULL, y = NULL
  ) +
  theme_minimal()

print(p)

# Save the plot to a PNG file
png_filename <- "Sentinel-2-airbird-kmeans.png"
ggsave(png_filename, plot = p, width = 6, height = 6, units = "in", dpi = 300)

print(paste("Saved plot to", png_filename))
