# kmeans-sentinel
using kmeans to classify sentinel images
# The purpose of this R script is to perform an unsupervised k-means classification on Sentinel-2 imagery for the Eyre Bird Observatory area. The script loads three bands (Red, Green, Blue), processes them to create a composite raster, clips the raster to a 3km buffer around the Eyre Bird Observatory coordinates, and then applies k-means clustering to classify the land cover into 6 distinct classes. Finally, it visualizes the results using ggplot2 and saves the plot as a PNG file.
# R script uses three inputed Sentinel-2 bands, converts them to a raster, clips them based on the Eyre Bird
# Observatory coordinates, and then runs a k-means classification on the clipped raster. Finally, it plots the results and
# saves the plot as a PNG file. Make sure to adjust file paths and band names as needed for your specific dataset and setup.

