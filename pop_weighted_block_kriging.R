# === 0) Install needed packages (if you haven’t already) ===
# install.packages(c("sf","sp","raster","gstat","exactextractr","tidyverse"))

# === 1) Load libraries ===
library(sf)
library(sp)
library(raster)
library(gstat)
library(exactextractr)
library(tidyverse)

# === 2) Unzip & read the PHU boundary shapefile ===
unzip("PHU_BOUNDARY.zip", exdir = "PHU_BOUNDARY")
shp_path <- list.files("PHU_BOUNDARY", pattern = "\\.shp$", full.names = TRUE)[1]
phu <- st_read(shp_path)

# === 3) Read the 2013 population‐density raster ===
pop_raster <- raster("can_pd_2013.tif")

# === 4) Read the station CSV (annual PM2.5 means) ===
pm_data <- read_csv("pm25_2013.csv") %>%
  filter(Pollutant == "PM2.5")  # keeps only PM2.5 records

# === 5) Convert station data to sf & reproject to match the raster CRS ===
stations <- st_as_sf(pm_data,
                     coords = c("Longitude","Latitude"),
                     crs    = 4326)
common_crs <- crs(pop_raster)
phu       <- st_transform(phu,       crs = common_crs)
stations  <- st_transform(stations,  crs = common_crs)

# === 6) Crop & mask the population raster to the PHU footprints ===
pop_masked <- pop_raster %>%
  crop(extent(phu)) %>%
  mask(phu)

# === 7) Build the kriging grid as the centers of the population raster cells ===
grid_pts <- rasterToPoints(pop_masked, spatial = TRUE)
# grid_pts is already a SpatialPointsDataFrame, suitable for gstat

# === 8) Define the two PM metrics to krige ===
metrics <- c("PM25_annu_mean_3hMax","PM25_annu_mean_8hMax")

# === 9) Loop over each metric: ordinary kriging → population‐weighted mean per PHU ===
for (met in metrics) {
  # a) attach the metric to the stations & drop any NA values
  stations_sf2 <- stations %>%
    mutate(val = .data[[met]]) %>%
    filter(!is.na(val))
  
  # convert sf → SpatialPointsDataFrame
  stations_spdf <- sf::as_Spatial(stations_sf2)  
  
  # b) compute empirical variogram & fit a Gaussian model
  v_emp <- variogram(val ~ 1, data = stations_spdf)
  v_mod <- fit.variogram(v_emp,
                         vgm(psill  = var(stations_spdf$val, na.rm = TRUE),
                             model  = "Gau",
                             range  = max(dist(coordinates(stations_spdf))) / 3,
                             nugget = var(stations_spdf$val, na.rm = TRUE) * 0.1))
  
  # c) ordinary kriging onto the grid
  kriged_sp <- krige(val ~ 1,
                     locations = stations_spdf,
                     newdata   = grid_pts,
                     model     = v_mod)
  
  #can do cross-validation if needed
  # leave-one-out cross-validation (do k-fold validation if more data)
  cv_sp <- krige.cv(val ~ 1,
                    locations = stations_spdf,
                    model     = v_mod)
  observed  <- stations_spdf$val
  predicted <- cv_sp@data$var1.pred
  rmse <- sqrt(mean((observed - predicted)^2, na.rm = TRUE))
  mae  <- mean(abs(observed - predicted), na.rm = TRUE)
  cat(met, "CV → RMSE:", round(rmse,3),
      " MAE:", round(mae,3), "\n")
  
  # d) convert kriged points → raster
  # get coordinates and prediction into a plain data.frame
  coords <- coordinates(kriged_sp)
  df_pred <- data.frame(
    x         = coords[,1],
    y         = coords[,2],
    var1.pred = kriged_sp$var1.pred
  )
  pred_rast <- rasterFromXYZ(df_pred, crs = common_crs)
  
  # e) mask the kriged raster to PHU boundaries
  pred_masked <- pred_rast %>%
    crop(extent(phu)) %>%
    mask(phu)
  
  # f) extract population‐weighted mean for each PHU
  field_name <- paste0("PopW_", met)
  phu[[field_name]] <- exact_extract(
    pred_masked,
    phu,
    fun           = "weighted_mean",
    weights = pop_masked
  )
}

# === 10) Inspect & plot the results ===
phu %>%
  select(starts_with("PopW_")) %>%
  st_drop_geometry() %>%
  head()

plot(phu["PopW_PM25_annu_mean_3hMax"],
     main = "Population‐Weighted PM2.5 (3hMax) by PHU")
plot(phu["PopW_PM25_annu_mean_8hMax"],
     main = "Population‐Weighted PM2.5 (8hMax) by PHU")


--------------------
library(RColorBrewer)
pal <- brewer.pal(7, "YlOrRd")

plot(
  phu["PopW_PM25_annu_mean_3hMax"],
  pal  = pal,
  main = "Population‐Weighted PM2.5 (3hMax) by PHU"
)

library(ggplot2)
ggplot(phu) +
  geom_sf(aes(fill = PopW_PM25_annu_mean_3hMax)) +
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlOrRd"),
    name   = "PM2.5 (3hMax)"
  )
-------------------------------------


#Model comparison

models <- list(
  Sph = vgm(psill=var(stations_spdf$val), model="Sph", range=300, nugget=0),
  Exp = vgm(psill=var(stations_spdf$val), model="Exp", range=300, nugget=0),
  Gau = vgm(psill=var(stations_spdf$val), model="Gau", range=300, nugget=0)
)

fits <- lapply(models, function(mod0) fit.variogram(v_emp, mod0))
SSErr <- sapply(fits, function(m) attr(m, "SSErr"))
print(SSErr)


plot(v_emp, fits$Sph, main="Variogram fit: Sph (black), Exp (blue), Gau (red)")
plot(v_emp, fits$Exp, add=TRUE, col="blue")
plot(v_emp, fits$Gau, add=TRUE, col="red")
legend("bottomright", legend=c("Sph","Exp","Gau"),
       col=c("black","blue","red"), lwd=2)
