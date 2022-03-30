#Packages necessary in order to conduct cluster reallocation analysis
source("DM/Tau_packages.R")

library(spatstat)
library(tmap)
library(rgdal)
library(sf)
library(dplyr)
library(geodist)

library(purrr)
library(kableExtra)
library(raster)
library(ggpubr)
library(permute)
library(furrr)

#calling in the functions needed also
source("Scripts/Functions/permutation-function.R")
source("Scripts/Functions/binary-matrix-function_odds.R")

#load the data
load("Data/2020-12-08_work-ts-spat.rdata")

#Need to filter out any unknown serotypes before starting Tau analysis 
work_ts_known <- work_ts_spat %>%
  filter(serotype != "unk_serotype")

#Need to transform dataset into a spatial dataframe, and transforming the 
#coordinate reference system to the one that it used for Indonesia
work_ts_spdf_known <- st_as_sf(work_ts_known, coords = c("longitude", "latitude"), crs=4326) %>% 
  st_transform(crs = 32749)

#There are now 6246 observations instead of 6313

class(work_ts_spdf_known)
#It is now also a simple feature dataset

st_geometry(work_ts_spdf_known)
#Shows that all of the 6246 observations are read as points in this geometry

#Obtaining shape files for study area
Yogya_Adm_1 <- st_read("Data/yogya-shape-files", "Yk_CaseControl_20160512")

#check the class of the simple feature 
class(Yogya_Adm_1)


#First alter dataset so that the cluster column matches that of the shape file 
#and then merge them together
int_assignment <- work_ts_spat%>% 
  dplyr::select(cluster, Intervention = intervention) %>% 
  distinct() %>% 
  mutate(New_Clustr = ifelse(cluster < 10, paste0("0", cluster), as.character(cluster))) %>% 
  dplyr::select(-cluster)
Yogya_Adm_1 <- left_join(Yogya_Adm_1, int_assignment, by = "New_Clustr")

#Transform the shape file so it is the same as the dataset
Yogya_Adm_1_UTM <- st_transform(Yogya_Adm_1, crs = 32749)

#Print all 24 features 
print(Yogya_Adm_1_UTM, n=24)
#Here, the fields (attributes) are ID (all=0), area(m^2), and cluster assignment

#Check the class of the combined dataset of the original data and shape file
class(Yogya_Adm_1_UTM)

#Change above dataset so it only includes intervention clusters
interventions <- filter(Yogya_Adm_1_UTM, Intervention== 1) 

#Add a buffer of 50m to the filtered intervention dataset above using st_buffer 
#and then join all boundaries together to get one shape using st_union 
#so there are no overlaps 
interventions_50 <- filter(Yogya_Adm_1_UTM, Intervention== 1) %>% 
  st_buffer(dist = 50) %>% 
  st_union()

#check the difference between the two plots
plot(st_geometry(interventions))
plot(st_geometry(interventions_50))


class(interventions)

st_geometry(interventions)
st_crs(interventions)

#plot all observations using st_geom 
plot(st_geometry(work_ts_spdf_known),add= TRUE)
#Including add=TRUE here adds a boundary outline to the plot 

#Create a logical vector such that if observations are within the 50m buffer
#of the intervention clusters, they are TRUE, otherwise it is FALSE
#Transform above vector so it is a binary vector i.e 0 if it is a control since
#it is not in the intervention buffer and 1 if it is in the intervention buffer
#i.e now an intervention 
intervention_dist_50 <- as.numeric(st_within(work_ts_spdf_known, interventions_50) %>%
                                  lengths>0) 

#Add this vector to the spatial dataframe to obtain "new" intervention column 
work_ts_known$New_intervention <- intervention_dist_50

#Then to obtain tau estimate, we partition data by intervention status,
#excluding those with unknown serotype 
work_ts_int <- work_ts_known %>%
  filter(New_intervention == 1)

work_ts_cont <- work_ts_known %>%
  filter(New_intervention ==0)


#Use an arbitrary distance to base tau off e.g 0-100m and time of 7 days
d1 <- 0 # lower bound distance
d2 <- 100 # upper bound distance
t <- 7 # time in days

#Calculate tau for intervention and see if it works 
observed_intervention_50 <- permutation_function(
  df = work_ts_int,
  r1 = d1,
  r2 = d2,
  t = t,
  permute = FALSE
)

observed_intervention_50 <- observed_intervention_50 %>% 
  mutate(intervention = "Intervention") %>%
  mutate(buffer = 50)

#Calculate tau for control and see if works
observed_untreated_50 <- permutation_function(
  df = work_ts_cont,
  r1 = d1,
  r2 = d2,
  t = t,
  permute = FALSE
)

observed_untreated_50 <- observed_untreated_50 %>% 
  mutate(intervention = "Untreated") %>%
  mutate(buffer = 50)



#Run some permutations in order to obtain a null distribution for intervention
perms <- as.vector(1:1000, mode = "list")

#apply the permutation function to each element of the vector i.e apply function
#1000 times to the intervention subset of data
null_intervention_50 <- future_map_dfr(perms, ~permutation_function(df = work_ts_int,
                                                                 r1 = d1,
                                                                 r2 = d2,
                                                                 t = t,
                                                                 permute = TRUE),
                                    .id = "permutation")

null_intervention_50 <- null_intervention_50 %>% 
  mutate(intervention = "Intervention")

#same as above except applying function to the control subset of data
null_untreated_50 <- future_map_dfr(perms, ~permutation_function(df = work_ts_cont,
                                                                 r1 = d1,
                                                                 r2 = d2,
                                                                 t = t,
                                                                 permute = TRUE),
                                    .id = "permutation")

null_untreated_50 <- null_untreated_50 %>% 
  mutate(intervention = "Untreated")


# generate quantile-based estimates of CI for each distance interval 
bounds_intervention_50 <- null_intervention_50 %>% 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))
bounds_intervention_50 <- bounds_intervention_50 %>%
  mutate(buffer = 50)


# generate quantile-based estimates of CI for each distance interval 
bounds_untreated_50 <- null_untreated_50 %>% 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))
bounds_untreated_50 <- bounds_untreated_50 %>%
  mutate(buffer = 50)

write.csv(null_intervention_50, "Data/Null_intervention_50.csv")

write.csv(null_untreated_50, "Data/Null_untreated_50.csv")

###################
#Now add a 100m buffer to intervention
interventions_100 <- filter(Yogya_Adm_1_UTM, Intervention== 1) %>% 
  st_buffer(dist = 100) %>% 
  st_union()

intervention_dist_100 <- as.numeric(st_within(work_ts_spdf_known, interventions_100) %>%
                                     lengths>0) 

work_ts_known$New_intervention <- intervention_dist_100

work_ts_int <- work_ts_known %>%
  filter(New_intervention == 1)

work_ts_cont <- work_ts_known %>%
  filter(New_intervention ==0)

observed_intervention_100 <- permutation_function(
  df = work_ts_int,
  r1 = d1,
  r2 = d2,
  t = t,
  permute = FALSE
)

observed_intervention_100 <- observed_intervention_100 %>% 
  mutate(intervention = "Intervention") %>%
  mutate(buffer = 100)

#Calculate tau for control and see if works
observed_untreated_100 <- permutation_function(
  df = work_ts_cont,
  r1 = d1,
  r2 = d2,
  t = t,
  permute = FALSE
)
observed_untreated_100 <- observed_untreated_100 %>% 
  mutate(intervention = "Untreated") %>%
  mutate(buffer = 100)

null_intervention_100 <- future_map_dfr(perms, ~permutation_function(df = work_ts_int,
                                                                    r1 = d1,
                                                                    r2 = d2,
                                                                    t = t,
                                                                    permute = TRUE),
                                       .id = "permutation")
null_intervention_100 <- null_intervention_100 %>% 
  mutate(intervention = "Intervention")

write.csv(null_intervention_100, "Data/Null_intervention_100.csv")

null_untreated_100 <- future_map_dfr(perms, ~permutation_function(df = work_ts_cont,
                                                                 r1 = d1,
                                                                 r2 = d2,
                                                                 t = t,
                                                                 permute = TRUE),
                                    .id = "permutation")

null_untreated_100 <- null_untreated_100 %>% 
  mutate(intervention = "Untreated")

write.csv(null_untreated_100, "Data/Null_untreated_100.csv")

bounds_intervention_100 <- null_intervention_100 %>% 
  # generate quantile-based estimates of CI for each distance interval 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))
bounds_intervention_100 <- bounds_intervention_100 %>%
  mutate(buffer = 100)


bounds_untreated_100 <- null_untreated_100 %>% 
  # generate quantile-based estimates of CI for each distance interval 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))
bounds_untreated_100 <- bounds_untreated_100 %>%
  mutate(buffer = 100)

#######
#Now 150m
#######

interventions_150 <- filter(Yogya_Adm_1_UTM, Intervention== 1) %>% 
  st_buffer(dist = 150) %>% 
  st_union()

intervention_dist_150 <- as.numeric(st_within(work_ts_spdf_known, interventions_150) %>%
                                      lengths>0) 

work_ts_known$New_intervention <- intervention_dist_150

work_ts_int <- work_ts_known %>%
  filter(New_intervention == 1)

work_ts_cont <- work_ts_known %>%
  filter(New_intervention ==0)

observed_intervention_150 <- permutation_function(
  df = work_ts_int,
  r1 = d1,
  r2 = d2,
  t = t,
  permute = FALSE
)

observed_intervention_150 <- observed_intervention_150 %>% 
  mutate(intervention = "Intervention") %>%
  mutate(buffer = 150)

#Calculate tau for control and see if works
observed_untreated_150 <- permutation_function(
  df = work_ts_cont,
  r1 = d1,
  r2 = d2,
  t = t,
  permute = FALSE
)
observed_untreated_150 <- observed_untreated_150 %>% 
  mutate(intervention = "Untreated") %>%
  mutate(buffer = 150)

null_intervention_150 <- future_map_dfr(perms, ~permutation_function(df = work_ts_int,
                                                                     r1 = d1,
                                                                     r2 = d2,
                                                                     t = t,
                                                                     permute = TRUE),
                                        .id = "permutation")
null_intervention_150 <- null_intervention_150 %>% 
  mutate(intervention = "Intervention")

write.csv(null_intervention_150, "Data/Null_intervention_150.csv")

null_untreated_150 <- future_map_dfr(perms, ~permutation_function(df = work_ts_cont,
                                                                  r1 = d1,
                                                                  r2 = d2,
                                                                  t = t,
                                                                  permute = TRUE),
                                     .id = "permutation")

null_untreated_150 <- null_untreated_150 %>% 
  mutate(intervention = "Untreated")

write.csv(null_untreated_150, "Data/Null_untreated_150.csv")

bounds_intervention_150 <- null_intervention_150 %>% 
  # generate quantile-based estimates of CI for each distance interval 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))
bounds_intervention_150 <- bounds_intervention_150 %>%
  mutate(buffer = 150)


bounds_untreated_150 <- null_untreated_150 %>% 
  # generate quantile-based estimates of CI for each distance interval 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))
bounds_untreated_150 <- bounds_untreated_150 %>%
  mutate(buffer= 150)

########
#200m
########
interventions_200 <- filter(Yogya_Adm_1_UTM, Intervention== 1) %>% 
  st_buffer(dist = 200) %>% 
  st_union()

intervention_dist_200 <- as.numeric(st_within(work_ts_spdf_known, interventions_200) %>%
                                      lengths>0) 

work_ts_known$New_intervention <- intervention_dist_200

work_ts_int <- work_ts_known %>%
  filter(New_intervention == 1)

work_ts_cont <- work_ts_known %>%
  filter(New_intervention ==0)

observed_intervention_200 <- permutation_function(
  df = work_ts_int,
  r1 = d1,
  r2 = d2,
  t = t,
  permute = FALSE
)

observed_intervention_200 <- observed_intervention_200 %>% 
  mutate(intervention = "Intervention") %>%
  mutate(buffer = 200)

#Calculate tau for control and see if works
observed_untreated_200 <- permutation_function(
  df = work_ts_cont,
  r1 = d1,
  r2 = d2,
  t = t,
  permute = FALSE
)
observed_untreated_200 <- observed_untreated_200 %>% 
  mutate(intervention = "Untreated") %>%
  mutate(buffer = 200)

null_intervention_200 <- future_map_dfr(perms, ~permutation_function(df = work_ts_int,
                                                                     r1 = d1,
                                                                     r2 = d2,
                                                                     t = t,
                                                                     permute = TRUE),
                                        .id = "permutation")
null_intervention_200 <- null_intervention_200 %>% 
  mutate(intervention = "Intervention")

write.csv(null_intervention_200, "Data/Null_intervention_200.csv")

null_untreated_200 <- future_map_dfr(perms, ~permutation_function(df = work_ts_cont,
                                                                  r1 = d1,
                                                                  r2 = d2,
                                                                  t = t,
                                                                  permute = TRUE),
                                     .id = "permutation")

null_untreated_200 <- null_untreated_200 %>% 
  mutate(intervention = "Untreated")

write.csv(null_untreated_200, "Data/Null_untreated_200.csv")

bounds_intervention_200 <- null_intervention_200 %>% 
  # generate quantile-based estimates of CI for each distance interval 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))
bounds_intervention_200 <- bounds_intervention_200 %>%
  mutate(buffer = 200)


bounds_untreated_200 <- null_untreated_200 %>% 
  # generate quantile-based estimates of CI for each distance interval 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))
bounds_untreated_200 <- bounds_untreated_200 %>%
  mutate(buffer = 200)

##################
#Now look at expanding controls, first at 50m
##################
controls <- filter(Yogya_Adm_1_UTM, Intervention==0)

controls_50 <- filter(Yogya_Adm_1_UTM, Intervention== 0) %>% 
  st_buffer(dist = 50) %>% 
  st_union()

plot(st_geometry(controls))
plot(st_geometry(controls_50))


control_dist_50 <- as.numeric(st_within(work_ts_spdf_known, controls_50) %>%
                                      lengths>0) 

#Now controls have value 1 and interventions have value 0
work_ts_known$New_intervention <- control_dist_50

#Now need to recode new intervention variable such that control=0, int=1 
work_ts_known <- work_ts_known %>%
  mutate(New_intervention = recode(New_intervention, `0`=1L, `1`=0L))


work_ts_int <- work_ts_known %>%
  filter(New_intervention == 1)

work_ts_cont <- work_ts_known %>%
  filter(New_intervention ==0)

observed_intervention_m50 <- permutation_function(
  df = work_ts_int,
  r1 = d1,
  r2 = d2,
  t = t,
  permute = FALSE
)

observed_intervention_m50 <- observed_intervention_m50 %>% 
  mutate(intervention = "Intervention") %>%
  mutate(buffer = -50)

#Calculate tau for control and see if works
observed_untreated_m50 <- permutation_function(
  df = work_ts_cont,
  r1 = d1,
  r2 = d2,
  t = t,
  permute = FALSE
)
observed_untreated_m50 <- observed_untreated_m50 %>% 
  mutate(intervention = "Untreated") %>%
  mutate(buffer = -50)

null_intervention_m50 <- future_map_dfr(perms, ~permutation_function(df = work_ts_int,
                                                                     r1 = d1,
                                                                     r2 = d2,
                                                                     t = t,
                                                                     permute = TRUE),
                                        .id = "permutation")
null_intervention_m50 <- null_intervention_m50 %>% 
  mutate(intervention = "Intervention")

write.csv(null_intervention_m50, "Data/Null_intervention_m50.csv")

null_untreated_m50 <- future_map_dfr(perms, ~permutation_function(df = work_ts_cont,
                                                                  r1 = d1,
                                                                  r2 = d2,
                                                                  t = t,
                                                                  permute = TRUE),
                                     .id = "permutation")

null_untreated_m50 <- null_untreated_m50 %>% 
  mutate(intervention = "Untreated")

write.csv(null_untreated_m50, "Data/Null_untreated_m50.csv")

bounds_intervention_m50 <- null_intervention_m50 %>% 
  # generate quantile-based estimates of CI for each distance interval 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))
bounds_intervention_m50 <- bounds_intervention_m50 %>%
  mutate(buffer = -50)


bounds_untreated_m50 <- null_untreated_m50 %>% 
  # generate quantile-based estimates of CI for each distance interval 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))
bounds_untreated_m50 <- bounds_untreated_m50 %>%
  mutate(buffer = -50)

#########################
#Now for controls expanded at 100m
###########################
controls_100 <- filter(Yogya_Adm_1_UTM, Intervention== 0) %>% 
  st_buffer(dist = 100) %>% 
  st_union()



control_dist_100 <- as.numeric(st_within(work_ts_spdf_known, controls_100) %>%
                                lengths>0) 

#Now controls have value 1 and interventions have value 0
work_ts_known$New_intervention <- control_dist_100

#Now need to recode new intervention variable such that control=0, int=1 
work_ts_known <- work_ts_known %>%
  mutate(New_intervention = recode(New_intervention, `0`=1L, `1`=0L))


work_ts_int <- work_ts_known %>%
  filter(New_intervention == 1)

work_ts_cont <- work_ts_known %>%
  filter(New_intervention ==0)

observed_intervention_m100 <- permutation_function(
  df = work_ts_int,
  r1 = d1,
  r2 = d2,
  t = t,
  permute = FALSE
)

observed_intervention_m100 <- observed_intervention_m100 %>% 
  mutate(intervention = "Intervention") %>%
  mutate(buffer = -100)

#Calculate tau for control and see if works
observed_untreated_m100 <- permutation_function(
  df = work_ts_cont,
  r1 = d1,
  r2 = d2,
  t = t,
  permute = FALSE
)
observed_untreated_m100 <- observed_untreated_m100 %>% 
  mutate(intervention = "Untreated") %>%
  mutate(buffer = -100)

null_intervention_m100 <- future_map_dfr(perms, ~permutation_function(df = work_ts_int,
                                                                     r1 = d1,
                                                                     r2 = d2,
                                                                     t = t,
                                                                     permute = TRUE),
                                        .id = "permutation")
null_intervention_m100 <- null_intervention_m100 %>% 
  mutate(intervention = "Intervention")

write.csv(null_intervention_m100, "Data/Null_intervention_m100.csv")

null_untreated_m100 <- future_map_dfr(perms, ~permutation_function(df = work_ts_cont,
                                                                  r1 = d1,
                                                                  r2 = d2,
                                                                  t = t,
                                                                  permute = TRUE),
                                     .id = "permutation")

null_untreated_m100 <- null_untreated_m100 %>% 
  mutate(intervention = "Untreated")

write.csv(null_untreated_m100, "Data/Null_untreated_m100.csv")

bounds_intervention_m100 <- null_intervention_m100 %>% 
  # generate quantile-based estimates of CI for each distance interval 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))
bounds_intervention_m100 <- bounds_intervention_m100 %>%
  mutate(buffer = -100)


bounds_untreated_m100 <- null_untreated_m100 %>% 
  # generate quantile-based estimates of CI for each distance interval 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))
bounds_untreated_m100 <- bounds_untreated_m100 %>%
  mutate(buffer = -100)

#######################
#Expand controls for 150m
######################
controls_150 <- filter(Yogya_Adm_1_UTM, Intervention== 0) %>% 
  st_buffer(dist = 150) %>% 
  st_union()



control_dist_150 <- as.numeric(st_within(work_ts_spdf_known, controls_150) %>%
                                 lengths>0) 

#Now controls have value 1 and interventions have value 0
work_ts_known$New_intervention <- control_dist_150

#Now need to recode new intervention variable such that control=0, int=1 
work_ts_known <- work_ts_known %>%
  mutate(New_intervention = recode(New_intervention, `0`=1L, `1`=0L))


work_ts_int <- work_ts_known %>%
  filter(New_intervention == 1)

work_ts_cont <- work_ts_known %>%
  filter(New_intervention ==0)

observed_intervention_m150 <- permutation_function(
  df = work_ts_int,
  r1 = d1,
  r2 = d2,
  t = t,
  permute = FALSE
)

observed_intervention_m150 <- observed_intervention_m150 %>% 
  mutate(intervention = "Intervention") %>%
  mutate(buffer = -150)

#Calculate tau for control and see if works
observed_untreated_m150 <- permutation_function(
  df = work_ts_cont,
  r1 = d1,
  r2 = d2,
  t = t,
  permute = FALSE
)
observed_untreated_m150 <- observed_untreated_m150 %>% 
  mutate(intervention = "Untreated") %>%
  mutate(buffer = -150)

null_intervention_m150 <- future_map_dfr(perms, ~permutation_function(df = work_ts_int,
                                                                      r1 = d1,
                                                                      r2 = d2,
                                                                      t = t,
                                                                      permute = TRUE),
                                         .id = "permutation")
null_intervention_m150 <- null_intervention_m150 %>% 
  mutate(intervention = "Intervention")

write.csv(null_intervention_m150, "Data/Null_intervention_m150.csv")

null_untreated_m150 <- future_map_dfr(perms, ~permutation_function(df = work_ts_cont,
                                                                   r1 = d1,
                                                                   r2 = d2,
                                                                   t = t,
                                                                   permute = TRUE),
                                      .id = "permutation")

null_untreated_m150 <- null_untreated_m150 %>% 
  mutate(intervention = "Untreated")

write.csv(null_untreated_m150, "Data/Null_untreated_m150.csv")

bounds_intervention_m150 <- null_intervention_m150 %>% 
  # generate quantile-based estimates of CI for each distance interval 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))
bounds_intervention_m150 <- bounds_intervention_m150 %>%
  mutate(buffer= -150)

bounds_untreated_m150 <- null_untreated_m150 %>% 
  # generate quantile-based estimates of CI for each distance interval 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))
bounds_untreated_m150 <- bounds_untreated_m150 %>%
  mutate(buffer = -150)

###################
#Finally controls expanded at 200m 
####################
controls_200 <- filter(Yogya_Adm_1_UTM, Intervention== 0) %>% 
  st_buffer(dist = 200) %>% 
  st_union()



control_dist_200 <- as.numeric(st_within(work_ts_spdf_known, controls_200) %>%
                                 lengths>0) 

#Now controls have value 1 and interventions have value 0
work_ts_known$New_intervention <- control_dist_200

#Now need to recode new intervention variable such that control=0, int=1 
work_ts_known <- work_ts_known %>%
  mutate(New_intervention = recode(New_intervention, `0`=1L, `1`=0L))


work_ts_int <- work_ts_known %>%
  filter(New_intervention == 1)

work_ts_cont <- work_ts_known %>%
  filter(New_intervention ==0)

observed_intervention_m200 <- permutation_function(
  df = work_ts_int,
  r1 = d1,
  r2 = d2,
  t = t,
  permute = FALSE
)

observed_intervention_m200 <- observed_intervention_m200 %>% 
  mutate(intervention = "Intervention") %>%
  mutate(buffer = -200)

#Calculate tau for control and see if works
observed_untreated_m200 <- permutation_function(
  df = work_ts_cont,
  r1 = d1,
  r2 = d2,
  t = t,
  permute = FALSE
)
observed_untreated_m200 <- observed_untreated_m200 %>% 
  mutate(intervention = "Untreated") %>%
  mutate(buffer = -200)
 

null_intervention_m200 <- future_map_dfr(perms, ~permutation_function(df = work_ts_int,
                                                                      r1 = d1,
                                                                      r2 = d2,
                                                                      t = t,
                                                                      permute = TRUE),
                                         .id = "permutation")
null_intervention_m200 <- null_intervention_m200 %>% 
  mutate(intervention = "Intervention")

write.csv(null_intervention_m200, "Data/Null_intervention_m200.csv")

null_untreated_m200 <- future_map_dfr(perms, ~permutation_function(df = work_ts_cont,
                                                                   r1 = d1,
                                                                   r2 = d2,
                                                                   t = t,
                                                                   permute = TRUE),
                                      .id = "permutation")

null_untreated_m200 <- null_untreated_m200 %>% 
  mutate(intervention = "Untreated")

write.csv(null_untreated_m200, "Data/Null_untreated_m200.csv")

bounds_intervention_m200 <- null_intervention_m200 %>% 
  # generate quantile-based estimates of CI for each distance interval 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))
bounds_intervention_m200 <- bounds_intervention_m200 %>%
  mutate(buffer= -200)

bounds_untreated_m200 <- null_untreated_m200 %>% 
  # generate quantile-based estimates of CI for each distance interval 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))
bounds_untreated_m200 <- bounds_intervention_m200 %>%
  mutate(buffer = -200)

########
#Obtain dataframe containing all 95% CIs for intervention and untreated group
########
Intervention_95CI <- rbind(bounds_intervention_m200, bounds_intervention_m150,
                           bounds_intervention_m100, bounds_intervention_m50,
                           bounds_intervention_1, bounds_intervention_50,
                           bounds_intervention_100, bounds_intervention_150,
                           bounds_intervention_200)

Untreated_95CI <- rbind(bounds_untreated_m200, bounds_untreated_m150,
                        bounds_untreated_m100, bounds_untreated_m50,
                        bounds_untreated_1, bounds_untreated_50,
                        bounds_untreated_100, bounds_untreated_150,
                        bounds_untreated_200)

#Now need dataframe containing point estimates for intervention and untreated
Full_observed_intervention <- rbind(observed_intervention_m200,
                                    observed_intervention_m150, 
                                    observed_intervention_m100,
                                    observed_intervention_m50,
                                    observed_intervention,
                                    observed_intervention_50,
                                    observed_intervention_100,
                                    observed_intervention_150,
                                    observed_intervention_200)

Full_observed_untreated <- rbind(observed_untreated_m200,
                                 observed_untreated_m150,
                                 observed_untreated_m100,
                                 observed_untreated_m50,
                                 observed_untreated,
                                 observed_untreated_50,
                                 observed_untreated_100,
                                 observed_untreated_150,
                                 observed_untreated_200)

################
#Try to plot the 95CIs for intervention on one graph 
###############
intervention_plot <- Full_observed_intervention %>%
  full_join(Intervention_95CI) %>%
  ggplot(aes(x = buffer, y= tau)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymax= CI.u, ymin= CI.l))+
  theme_minimal()+
  labs(x = "Buffer Distance (metres)",
       y = TeX("Estimated $\\tau"),
       subtitle = "Intervention arm")
intervention_plot

untreated_plot <- Full_observed_untreated %>%
  full_join(Untreated_95CI) %>%
  ggplot(aes(x = buffer, y= tau)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymax= CI.u, ymin= CI.l))+
  theme_minimal()+
  labs(x = "Buffer Distance (metres)",
       y = TeX("Estimated $\\tau"),
       subtitle = "Untreated arm")

untreated_plot
