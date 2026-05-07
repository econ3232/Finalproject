library(sf)
library(terra)
library(dplyr)
library(ggplot2)

# ---------------------------
# 1. POINTS + BUFFERS
# ---------------------------
sites  <- c("Serpent Mound", "Nobles Pond", "Sheridan Cave")
coords <- data.frame(
  labels = sites,
  longt = c(-83.4302, -81.4819, -83.444763),
  latit = c(39.0252, 40.8552, 40.980496)
)

arch_points <- st_as_sf(coords, coords = c("longt", "latit"), crs = 4326)
arch_points$site <- c("Serpent", "Nobles", "Sheridan")

crs_proj <- 26917
arch_points <- st_transform(arch_points, crs_proj)

buffers_sf <- st_buffer(arch_points, 1000)

buffer_serp <- buffers_sf %>% filter(site == "Serpent")
buffer_np   <- buffers_sf %>% filter(site == "Nobles")
buffer_sc   <- buffers_sf %>% filter(site == "Sheridan")

# ---------------------------
# 2. DEM + SLOPE
# ---------------------------
serpdem <- rast("./Ohio/USGS_13_n40w084_20250915 (1).tif")
NPdem   <- rast("./Ohio/USGS_13_n41w082_20230911 (1).tif")
SCdem   <- rast("./Ohio/USGS_13_n41w084_20250915.tif")

serpdem <- project(serpdem, paste0("EPSG:", crs_proj))
NPdem   <- project(NPdem, paste0("EPSG:", crs_proj))
SCdem   <- project(SCdem, paste0("EPSG:", crs_proj))


# ---------------------------
# 2. SLOPE SEVERITY (ROBUST FIX)
# ---------------------------
# ---------------------------
# 2. MEAN SLOPE + CLASSIFICATION
# ---------------------------

# --- SERPENT ---
serp_crop <- crop(serpdem, vect(buffer_serp)) |> mask(vect(buffer_serp))
slope1 <- terrain(serp_crop, v="slope", unit="degrees")

buffer_serp$slope_mean <- terra::extract(slope1, vect(buffer_serp), fun=mean, na.rm=TRUE)[,2]


# --- NOBLES ---
np_crop <- crop(NPdem, vect(buffer_np)) |> mask(vect(buffer_np))
slope2 <- terrain(np_crop, v="slope", unit="degrees")

buffer_np$slope_mean <- terra::extract(slope2, vect(buffer_np), fun=mean, na.rm=TRUE)[,2]


# --- SHERIDAN ---
sc_crop <- crop(SCdem, vect(buffer_sc)) |> mask(vect(buffer_sc))
slope3 <- terrain(sc_crop, v="slope", unit="degrees")

buffer_sc$slope_mean <- terra::extract(slope3, vect(buffer_sc), fun=mean, na.rm=TRUE)[,2]


# recombine
buffers_sf <- rbind(buffer_serp, buffer_np, buffer_sc)

buffers_sf$slope_class <- case_when(
  buffers_sf$slope_mean < 3  ~ 1,   # low
  buffers_sf$slope_mean < 6  ~ 2,   # moderate
  TRUE                        ~ 3    # high
)

# ---------------------------
# 3. HYDROGRAPHY
# ---------------------------
serphydro <- rbind(
  st_read("./Ohio/Shape/NHDFlowline.shp"),
  st_read("./Ohio/Shape2/Shape/NHDFlowline.shp"),
  st_read("./Ohio/Shape3/Shape/NHDFlowline.shp"),
  st_read("./Ohio/Shape/Shape/NHDFlowline.shp")
)

NPhydro <- st_read("./Ohio/Shape4/Shape/Shape/NHDFlowline.shp")

SChydro <- rbind(
  st_read("./Ohio/Shape5/Shape/Shape/NHDFlowline.shp"),
  st_read("./Ohio/Shape6/Shape/NHDFlowline.shp")
)

serphydro <- st_transform(serphydro, crs_proj) |> st_zm(drop = TRUE)
NPhydro   <- st_transform(NPhydro, crs_proj)   |> st_zm(drop = TRUE)
SChydro   <- st_transform(SChydro, crs_proj)   |> st_zm(drop = TRUE)

intersects1 <- st_intersection(serphydro, buffer_serp)
intersects2 <- st_intersection(NPhydro, buffer_np)
intersects3 <- st_intersection(SChydro, buffer_sc)

# HYDROGRAPHY DISTANCE FACTOR
# ---------------------------

# combine all streams
streams_all <- rbind(serphydro, NPhydro, SChydro)
streams_all <- st_union(streams_all)

# distance from each buffer to nearest stream
buffers_sf$stream_dist <- apply(
  st_distance(buffers_sf, streams_all),
  1,
  min
)

# convert to factor (closer = higher erosion)
buffers_sf$stream_factor <- 1 / (as.numeric(buffers_sf$stream_dist) + 1)




# ---------------------------
# 4. SOILS (CLEAN VERSION)
# ---------------------------
clay <- rast("./Ohio/out (3).tif")
sand <- rast("./Ohio/out (2).tif")

clay <- project(clay, "EPSG:26917")
sand <- project(sand, "EPSG:26917")

buffers_sf$clay_mean <- terra::extract(clay, vect(buffers_sf), fun = mean, na.rm = TRUE)[,2]

buffers_sf$sand_mean <- terra::extract(sand, vect(buffers_sf), fun = mean, na.rm = TRUE)[,2]


buffers_sf$clay_mean <- buffers_sf$clay_mean / 10
buffers_sf$sand_mean <- buffers_sf$sand_mean / 10

buffers_sf$soil_risk <- 
  (buffers_sf$clay_mean * 0.6) +
  ((100 - buffers_sf$sand_mean) * 0.4)

buffers_sf$slope_class <- case_when(
  buffers_sf$slope_mean < 3 ~ "Low",
  buffers_sf$slope_mean < 6 ~ "Moderate",
  TRUE ~ "High"
)
buffers_sf$stream_class <- case_when(
  buffers_sf$stream_dist < 500  ~ "High",
  buffers_sf$stream_dist < 1500 ~ "Moderate",
  TRUE ~ "Low"
)
buffers_sf <- buffers_sf %>%
  mutate(
    clay_n = (clay_mean - min(clay_mean, na.rm = TRUE)) /
      (max(clay_mean, na.rm = TRUE) - min(clay_mean, na.rm = TRUE)),
    
    sand_n = (sand_mean - min(sand_mean, na.rm = TRUE)) /
      (max(sand_mean, na.rm = TRUE) - min(sand_mean, na.rm = TRUE))
  )
buffers_sf$soil_class <- case_when(
  buffers_sf$sand_n > 0.66 & buffers_sf$clay_n < 0.33 ~ "High",
  buffers_sf$sand_n < 0.33 & buffers_sf$clay_n > 0.66 ~ "Low",
  TRUE ~ "Moderate"
)

class_score <- function(x){
  case_when(
    x == "Low" ~ 1,
    x == "Moderate" ~ 2,
    x == "High" ~ 3
  )
}

buffers_sf$erosion_risk <- 
  class_score(buffers_sf$slope_class) +
  class_score(buffers_sf$stream_class) +
  class_score(buffers_sf$soil_class)


buffers_sf$erosion_class <- case_when(
  buffers_sf$erosion_risk <= 5 ~ "Low",
  buffers_sf$erosion_risk <= 7 ~ "Moderate",
  TRUE ~ "High"
)



# ---------------------------
# 5. OHIO BASE MAP
# ---------------------------
ohio <- st_read("./Ohio/Shape/GU_StateorTerritory.shp") |>
  st_transform(crs_proj)

# ---------------------------
# 6. PLOT
# ---------------------------
ggplot() +
  geom_sf(data = ohio, fill = "gray95", color = "black") +
  geom_sf(data = intersects1, color = "blue") +
  geom_sf(data = intersects2, color = "blue") +
  geom_sf(data = intersects3, color = "blue") +
  geom_sf(data = buffers_sf, fill = NA, color = "black") +
  geom_sf(data = arch_points, aes(color = labels), size = 3) +
  theme_minimal()



library(leaflet)

pal_slope <- colorNumeric("YlOrBr", buffers_sf$slope_mean)
pal_stream <- colorNumeric("Blues", buffers_sf$stream_factor)
pal_soil  <- colorNumeric("Greens", buffers_sf$soil_risk)
pal_eros  <- colorNumeric("Reds", buffers_sf$erosion_risk)


buffers_sf  <- st_transform(buffers_sf, 26917)
arch_points <- st_transform(arch_points, 26917)
ohio        <- st_transform(ohio, 26917)

library(ggplot2)

risk_colors <- c(
  "Low" = "green3",
  "Moderate" = "orange",
  "High" = "red3"
)

ggplot() +
  geom_sf(data = ohio, fill = "gray95", color = "black") +
  
  geom_sf(
    data = buffers_sf,
    aes(color = erosion_class),
    size = 4
  ) +
  
  scale_color_manual(values = risk_colors) +
  
  theme_minimal() +
  ggtitle("Erosion Risk Class by Site")