
#load the dataset 
load("2020-12-08_work-ts-spat.rdata")

#Do some basic analysis of the variables e.g number of dengue cases and
#test negative controls 
table(work_ts_spat$dengue)
#There are 6313 participants with 392 being dengue cases and 5921 being controls 

#specifically to serotype, we get a table of these
table(work_ts_spat$serotype)
#Most cases were denv2, fewest were denv3 and there were 67 cases that couldn't 
#Be identified 


#Summarise all the variables
summary(work_ts_spat)

#Here, first enrolment date was 08/01/2018 and last was 18/03/2020
#First date of illness onset was 04/01/2020 and last was 17/03/2020
#People in the study were between 3 and 45 with the mean age being 15 

###
###Now will try and produce some nice maps
####

install.packages("spatstat")
install.packages("rgdal")
install.packages("knitr")
install.packages("dplyr")
install.packages("sf")
install.packages("tmap")
install.packages("geepack")
library(geepack)
library(spatstat)
library(tmap)
library(rgdal)
library(dplyr)
library(sf)
# Setting up a color palette for the different DENV serotypes
df_col_sero <- data.frame(serotype = c(paste0("DENV", 1:4), "Unknown", "Test-negative"),
                          sero.color = c("#0AA4D1", "#F3C73E", "#EF7822", "#76C1A8", "#003D58", "gray"),
                          sero.shape = c(21:24, 4, 25))

work_ts_spat_2 <- work_ts_spat %>% 
  mutate(serotype = case_when(serotype == "denv1" ~ "DENV1",
                              serotype == "denv2" ~ "DENV2",
                              serotype == "denv3" ~ "DENV3",
                              serotype == "denv4" ~ "DENV4",
                              serotype == "unk_serotype"~ "Unknown",
                              serotype == "test-negative control" ~ "Test-negative"),
         intervention = ifelse(intervention == 0, "Untreated", "Intervention")) %>% 
  full_join(df_col_sero) 

# Getting the "shape file" for the study area
Yogya_Adm_1 <- readOGR("C:/Users/jerom/OneDrive/Documents/Med Stats 2020/Summer Project/Data/yogya-shape-files", "Yk_CaseControl_20160512")
int_assignment <- work_ts_spat_2 %>% 
  dplyr::select(cluster, Intervention = intervention) %>% 
  distinct() %>% 
  mutate(New_Clustr = ifelse(cluster < 10, paste0("0", cluster), as.character(cluster))) %>% 
  dplyr::select(-cluster)
Yogya_Adm_1@data <- left_join(Yogya_Adm_1@data, int_assignment, by = "New_Clustr")

# Creating Spatial Points Dataframe 
# This is so that the points from our dataset can be mapped onto the study area map
CRTND_dengue_data_SPDF <- SpatialPointsDataFrame(coords = work_ts_spat_2[,c("longitude", "latitude")],
                                                 data = work_ts_spat_2[, c("cluster", "intervention", "enrolment_date", "illness_onset", "dengue", "age", "sex", "serotype", "sero.color", "sero.shape")],
                                                 proj4string = CRS("+init=epsg:4326")) # sets the projection to WGS 1984 using lat/long. Optional but good to specify

# Setting up a dataset with just the dengue cases
CRTND_VCD_data_SPDF <- subset(CRTND_dengue_data_SPDF, dengue == 1)

#####
#Map Preliminaries
#####
map_yogya <- 
  tm_shape(Yogya_Adm_1) +
  tm_borders()

map_yogya + 
  tm_text(text = "New_Clustr")

##########
#Dengue Point Plot
############

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
  tm_layout(legend.text.size = 0.5,
            legend.title.size = 0.1,
            legend.title.color = "white",
            legend.position= c("right","bottom"),
            legend.width = 0.5,
            inner.margins = 0.14) +
  tm_scale_bar(breaks = c(0,.50,1.00), text.size = 0.7, position = c("left", "top"))

##########
#Mixed effect + GEE model
##########
case_cohort_function <- function(dataset, d1, d2, time){
  
  temp <- dataset %>% 
    # CRITICAL - put the individuals in time order according to date of illness onset
    arrange(illness_onset)
  
  # Identify the test-positive dengue cases
  vcd_ids <- temp$participant_id[temp$dengue == 1]
  
  # Setting up distance and time matrices
  dist_temp <- geodist(x = cbind(lon = temp$longitude, lat = temp$latitude), y = cbind(lon = temp$longitude, lat = temp$latitude))
  time_temp <- expand.grid(temp$illness_onset, temp$illness_onset)
  time_temp <- matrix(as.Date(time_temp$Var1) - as.Date(time_temp$Var2), nrow = nrow(temp))
  
  diag(dist_temp) <- diag(time_temp) <- NA
  rownames(dist_temp) <- rownames(time_temp) <- temp$participant_id
  colnames(dist_temp) <- colnames(time_temp) <- temp$participant_id
  
  output <- NULL
  exposures <- NULL
  iter <- 1
  for (i in 1:length(time)){
    for (j in 1:length(d2)){
      time_temp.bin <- 1*(time_temp < time[i] & time_temp >= 0)
      dist_temp.bin <- 1*(dist_temp >= d1[j] & dist_temp < d2[j])
      prod_mat <- time_temp.bin*dist_temp.bin # only indicates if anyone is within time and space of each other, to determine exposure, need to focus on vcds
      exposure_set <- prod_mat[,which(colnames(prod_mat) %in% vcd_ids)]
      t1 <- unlist(apply(exposure_set, 1, function(x){max(x,na.rm = TRUE)}))
      
      temp2 <- data.frame(participant_id = names(t1), t1) %>% 
        full_join(dplyr::select(temp, c(participant_id, dengue))) %>% 
        group_by(dengue) %>% 
        summarise(exposed = sum(t1)) %>% 
        mutate(time = time[i],
               distance1 = d1[j],
               distance2 = d2[j])
      
      # Setting up exposure data for the individual-level analyses (GEE, mixed effects, etc.)
      exposures <- cbind(exposures, t1)
      colnames(exposures)[iter] <- paste0("d", d2[j],".t", time[i])
      
      iter <- iter + 1
      output <- bind_rows(output, temp2)
    }
  }
  
  ## Producing simple contingency tables and contingency table-based estimates
  contingency.table <- temp %>% 
    group_by(dengue) %>% 
    summarise(total = n_distinct(participant_id)) %>% 
    full_join(output) %>% 
    mutate(unexposed = total - exposed) %>% 
    group_by(time, distance1, distance2) %>% 
    mutate(numerator = exposed[dengue == 1]*unexposed[dengue == 0], 
           denominator = exposed[dengue == 0]*unexposed[dengue == 1]) %>% 
    mutate(rr = numerator/denominator) %>% 
    mutate(var.log.rr = 1/exposed[dengue == 1] + 1/exposed[dengue == 0] + 1/unexposed[dengue == 1] + 1/unexposed[dengue == 0]) %>% 
    mutate(CI.l = exp(log(rr) - 1.96*sqrt(var.log.rr)),
           CI.u = exp(log(rr) + 1.96*sqrt(var.log.rr)))
  
  ## More complicated individual-level modeling (GEE, mixed effects, etc.)
  mixed.effects.results <- NULL
  gee.results <- NULL
  for (i in 1:ncol(exposures)){ # loop through the various exposures, as determined by times and distances
    mixed.effects.model <-  temp %>% 
      full_join(data.frame(participant_id = rownames(exposures), exposure = exposures[,i])) %>% 
      arrange(cluster) %>% 
      glmer(dengue ~ (1 | cluster) + exposure,
            family = binomial,
            data = .)
    
    mixed.effects.results <- append(mixed.effects.results,
                                    mixed.effects.model)
    names(mixed.effects.results)[i] <- colnames(exposures)[i]
    
    gee.model <-  temp %>% 
      full_join(data.frame(participant_id = rownames(exposures), exposure = exposures[,i])) %>% 
      arrange(cluster) %>% 
      geeglm(dengue ~ exposure,
             family = binomial ("logit"),
             id = cluster,
             scale.fix = TRUE,
             data = .)
    
    sum.gee.model <- summary(gee.model)$coefficients %>% 
      mutate(distance.time = colnames(exposures)[i])
    gee.results <- rbind(gee.results,
                         sum.gee.model)
  }
  
  overall <- list(contingency.table = contingency.table,
                  mixed.effects.results = mixed.effects.results,
                  gee.results = gee.results)
  
  return(overall)
}

install.packages("here")
install.packages("lme4")
install.packages("tidyverse")
install.packages("geodist")
install.packages("kableExtra")
install.packages("ggpubr")
library(ggpubr)
library(here)
library(tidyverse)
library(geodist)
library(knitr)
library(kableExtra)
library(lme4)

########
#Setting up the data, excluding unknown serotype
#######
work_all_dengue <- work_ts_spat %>% 
  # Removing the unknown serotypes
  filter(serotype != "unk_serotype") %>% 
  dplyr::select(-serotype) %>% 
  distinct()

##########
#Setting up the time and distance variables
#########
# time
t <- c(7,14,30) # MODIFY THIS
# distances 
r.upper <- c(100,200,500) # MODIFY THIS
r.lower <- c(0,0,0) # MODIFY THIS - must be the same length as r.upper

output_all_t <- case_cohort_function(dataset = work_ts_spat,
                                     d1 = r.lower,
                                     d2 = r.upper,
                                     time = t)


#########
#GEE Results
###########
gee.estimates <- output_all_t$gee.results %>% 
  rownames_to_column(var = "Coefficients") %>% 
  mutate(or = exp(Estimate),
         CI.l = exp(Estimate - 1.96*`Std.err`),
         CI.u = exp(Estimate + 1.96*`Std.err`),
         distance = as.numeric(str_match(distance.time, "d\\s*(.*?)\\.t")[,2]),
         days = str_split_fixed(distance.time, pattern = ".t", n = 2)[,2]) %>% 
  mutate(days = factor(days, levels = c("7", "14", "30"))) %>% 
  filter(str_detect(Coefficients, pattern = "exposure"))

gee.estimates %>% 
  ggplot(aes(x = distance, y = or, col = days, fill = days, group = days)) + 
  facet_wrap(~days, ncol = 1) +
  geom_line() +
  geom_ribbon(aes(ymin = CI.l, ymax = CI.u),
              alpha = 0.3,
              col = "white") +
  geom_hline(yintercept = 1, lty = 2) + 
  ggpubr::theme_pubr() +
  xlab("Distance from a VCD (m)") +
  ylab("Case-Cohort Odds Ratio") + 
  theme(legend.position = "none")

gee.estimates %>% 
  dplyr::select(distance, days, or, CI.l, CI.u) %>% 
  kable(digits = 2) %>% 
  kable_styling() %>% 
  collapse_rows()




