######
#Maps and cluster reallocation etc

source("DM/Map_packages.R")

library(spatstat)
library(tmap)
library(rgdal)
library(sf)
library(dplyr)

load("Data/2020-12-08_work-ts-spat.rdata")


#Setting up colour/shape palette for dengue serotypes
df_col_sero <- data.frame(serotype = c(paste0("DENV", 1:4), "Unknown", "Test-negative"),
                          sero.color = c("#0AA4D1", "#F3C73E", "#EF7822", "#76C1A8", "#003D58", "gray"),
                          sero.shape = c(21:24, 4, 25))


#Recoding variables in dataset and joining with colour/shape palette
work_ts_spat_2 <- work_ts_spat %>% 
  mutate(serotype = case_when(serotype == "denv1" ~ "DENV1",
                              serotype == "denv2" ~ "DENV2",
                              serotype == "denv3" ~ "DENV3",
                              serotype == "denv4" ~ "DENV4",
                              serotype == "unk_serotype"~ "Unknown",
                              serotype == "test-negative control" ~ "Test-negative"),
         intervention = ifelse(intervention == 0, "Untreated", "Intervention")) %>% 
  full_join(df_col_sero) 


#Obtaining shape files for study area
Yogya_Adm_1 <- st_read("Data/yogya-shape-files", "Yk_CaseControl_20160512")

#Transform onto a new projection 
Yogya_Adm_1_UTM <- st_transform(Yogya_Adm_1, crs = 32749)


#change dataset such that it recodes the cluster number and joins with sf data
int_assignment <- work_ts_spat_2 %>% 
  dplyr::select(cluster, Intervention = intervention) %>% 
  distinct() %>% 
  mutate(New_Clustr = ifelse(cluster < 10, paste0("0", cluster), as.character(cluster))) %>% 
  dplyr::select(-cluster)
Yogya_Adm_1_UTM <- left_join(Yogya_Adm_1_UTM, int_assignment, by = "New_Clustr")


##############
#Map stuff
df_col_sero <- data.frame(serotype = c(paste0("DENV", 1:4), "Unknown", "Test-negative"),
                          sero.color = c("#0AA4D1", "#F3C73E", "#EF7822", "#76C1A8", "#003D58", "gray"),
                          sero.shape = c(21:24, 4, 25))


#Obtain the simple map showing each cluster's ID
map_yogya <- 
  tm_shape(Yogya_Adm_1_UTM) +
  tm_borders()

map_yogya + 
  tm_text(text = "New_Clustr")

#Also create a basic map only showing the cluster ID and intervention class

map_yogya + 
  tm_text(text = "New_Clustr") +
  tm_polygons(col = "Intervention",
              alpha = 0.5,
              palette = c("lightgray", "white"))


#Create a spatial points data frame 
#So each observation can be mapped to study area map

CRTND_dengue_data_SPDF <- SpatialPointsDataFrame(coords = work_ts_spat_2[,c("longitude", "latitude")],
                                                 data = work_ts_spat_2[, c("cluster", "intervention", "enrolment_date", "illness_onset", "dengue", "age", "sex", "serotype", "sero.color", "sero.shape")],
                                                 proj4string = CRS("+init=epsg:4326")) # sets the projection to WGS 1984 using lat/long. Optional but good to specify

# Setting up a dataset with just the dengue cases

CRTND_VCD_data_SPDF <- subset(CRTND_dengue_data_SPDF, dengue == 1)

#Map containing all dengue points 

map_yogya +
  tm_polygons(col = "Intervention",
              alpha = 0.5,
              palette = c("lightgray", "white")) +
  tm_shape(subset(CRTND_dengue_data_SPDF, dengue == 1)) + 
  tm_symbols(col = "serotype",
             shape = "serotype",
             legend.col.show = FALSE,
             legend.shape.show = FALSE,
             palette = df_col_sero$sero.color,
             shapes = df_col_sero$sero.shape,
             size = 0.1) + 
  tm_add_legend("symbol", 
                col = df_col_sero$sero.color[1:5], 
                shape = df_col_sero$sero.shape[1:5],
                labels = df_col_sero$serotype[1:5]) +
  tm_layout(legend.text.size = 0.8,
            legend.title.size = 0.1,
            legend.title.color = "white",
            legend.outside = TRUE,
            legend.position= c("right","bottom"),
            legend.width = 0.5,
            inner.margins = 0.1) +
  tm_scale_bar(breaks = c(0,.50,1.00), text.size = 0.7, position = c("left", "top"))

#Try and use st_buffer next 
