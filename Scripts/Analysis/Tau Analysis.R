#Tau Analysis

source("DM/Tau_packages.R")

library(spatstat)
library(tmap)
library(rgdal)
library(sf)
library(dplyr)
library(geodist)
library(latex2exp)
library(purrr)
library(kableExtra)
library(raster)
library(ggpubr)
library(permute)
library(furrr)

#calling in the functions
source("Scripts/Functions/permutation-function.R")
source("Scripts/Functions/binary-matrix-function_odds.R")

load("Data/2020-12-08_work-ts-spat.rdata")

#Also need to filter out any unknown serotypes before starting Tau analysis 
work_ts_known <- work_ts_spat %>%
  filter(serotype != "unk_serotype")

#Partition known dataset by intervention status
work_ts_int <- work_ts_known %>%
  filter(intervention == 1)

work_ts_cont <- work_ts_known %>%
  filter(intervention ==0)

#First calculate tau at 7 days 
##Put in the time and distance parameters
d1 <- c(0, 0, 0, 0, 0) # lower bounds distance
d2 <- c(100, 200, 500, 750, 1000) # upper bounds distance
t <- 7 # time in days


# Generate the observed estimates of tau for intervention data
observed_intervention <- permutation_function(
  df = work_ts_int,
  r1 = d1,
  r2 = d2,
  t = t,
  permute = FALSE
)

observed_intervention <- observed_intervention %>% 
  mutate(intervention = "Intervention")

# Generate the observed estimates of tau for untreated data
observed_untreated <- permutation_function(
  df = work_ts_cont,
  r1 = d1,
  r2 = d2,
  t = t,
  permute = FALSE
)

observed_untreated <- observed_untreated %>% 
  mutate(intervention = "Untreated")


perms <- as.vector(1:1000, mode = "list")

#Obtain null distribution for intervention arm using permutations
null_intervention_7days <- future_map_dfr(perms, ~permutation_function(df = work_ts_int,
                                                                 r1 = d1,
                                                                 r2 = d2,
                                                                 t = t,
                                                                 permute = TRUE),
                                    .id = "permutation")

null_intervention_7days <- null_intervention_7days %>% 
  mutate(intervention = "Intervention")

#Look at the results from intervention arm 
null_intervention %>% 
  kable(digits = 3) %>% 
  kable_styling()

#Would expect distribution of tau values to be around 1
null_intervention %>% 
  mutate(permutation = factor(permutation, levels = 1:10)) %>% # for plotting purposes
  ggplot(aes(x = r_upper, y = tau, col = permutation)) + 
  geom_point(position = position_jitter(width = 10, height = 0)) + 
  geom_hline(yintercept = 1,
             lty = 2) +
  coord_cartesian(ylim = c(0,8)) + 
  theme_minimal() + 
  labs(x = TeX("Distance: $d_2"),
       y = TeX("Odds Ratio: $\\tau"),
       title = TeX("Estimated arm-level $\\tau")) +
  theme(legend.title = element_blank())

#Obtain null distribution for control arm using permutations 
null_untreated <- future_map_dfr(perms, ~permutation_function(df = work_ts_cont,
                                                              r1 = d1,
                                                              r2 = d2,
                                                              t = t,
                                                              permute = TRUE),
                                 .id = "permutation")

null_untreated <- null_untreated %>% 
  mutate(intervention = "Untreated")

#Look at the distribution of tau values for untreated arm as well:
null_untreated %>% 
  mutate(permutation = factor(permutation, levels = 1:10)) %>% # for plotting purposes
  ggplot(aes(x = r_upper, y = tau, col = permutation)) + 
  geom_point(position = position_jitter(width = 10, height = 0)) + 
  geom_hline(yintercept = 1,
             lty = 2) +
  coord_cartesian(ylim = c(0,8)) + 
  theme_minimal() + 
  labs(x = TeX("Distance: $d_2"),
       y = TeX("Odds Ratio: $\\tau"),
       title = TeX("Estimated arm-level $\\tau")) +
  theme(legend.title = element_blank())



#####
#Now obtain 95% Confidence Intervals based on percentiles for each arm 
bounds_intervention <- null_intervention %>% 
  # generate quantile-based estimates of CI for each distance interval 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))


bounds_untreated <- null_untreated %>% 
  # generate quantile-based estimates of CI for each distance interval 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))

#Finally obtain a graph showing both observed tau measurements for each arm with 
#the 95% CIs 

p1 <- observed_intervention %>% 
  full_join(bounds_intervention) %>% 
  ggplot(aes(x = r_upper, y = tau)) + 
  geom_line() +
  geom_ribbon(aes(ymin = CI.l,
                  ymax = CI.u),
              col = "lightgray",
              alpha = 0.3) +
  geom_hline(yintercept = 1,
             lty = 2) +
  coord_cartesian(ylim = c(0,8)) + 
  theme_minimal() + 
  labs(x = TeX("Distance: $d_2"),
       y = TeX("Odds Ratio: $\\tau"),
       title = TeX("Estimated $\\tau"),
       subtitle = "Intervention arm") +
  theme(legend.title = element_blank())


p2 <- observed_untreated %>% 
  full_join(bounds_untreated) %>% 
  ggplot(aes(x = r_upper, y = tau)) + 
  geom_line() +
  geom_ribbon(aes(ymin = CI.l,
                  ymax = CI.u),
              col = "lightgray",
              alpha = 0.3) +
  geom_hline(yintercept = 1,
             lty = 2) +
  coord_cartesian(ylim = c(0,6)) + 
  theme_minimal() + 
  labs(x = TeX("Distance: $d_2"),
       y = TeX("Odds Ratio: $\\tau"),
       title = TeX("Estimated $\\tau"),
       subtitle = "Untreated arm") +
  theme(legend.title = element_blank())

ggarrange(p1, p2)

#Save the two data frames since they take forever to load
write.csv(null_intervention, "Data/Null_intervention.csv")

write.csv(null_untreated, "Data/Null_untreated.csv")

#######
#Observe the same but for 14 and 30 days 
d1_14 <- c(0,0,0,0,0) # lower bounds distance
d2_14 <- c(100,200,500,750,1000) # upper bounds distance
t_14 <- 14 # time in days

observed_intervention_14 <- permutation_function(
  df = work_ts_int,
  r1 = d1,
  r2 = d2,
  t = t_14,
  permute = FALSE
)

observed_intervention_14 <- observed_intervention_14 %>% 
  mutate(Days = "14")

observed_untreated_14 <- permutation_function(
  df = work_ts_cont,
  r1 = d1,
  r2 = d2,
  t = t_14,
  permute = FALSE
)

observed_untreated_14 <- observed_untreated_14 %>% 
  mutate(Days = "14")

#Obtain null distribution for intervention arm using permutations
null_intervention_14 <- future_map_dfr(perms, ~permutation_function(df = work_ts_int,
                                                                 r1 = d1_14,
                                                                 r2 = d2_14,
                                                                 t = t_14,
                                                                 permute = TRUE),
                                    .id = "permutation")

null_intervention_14 <- null_intervention_14 %>% 
  mutate(intervention = "Intervention")

#Obtain null distribution for control arm using permutations 
null_untreated_14 <- future_map_dfr(perms, ~permutation_function(df = work_ts_cont,
                                                              r1 = d1_14,
                                                              r2 = d2_14,
                                                              t = t_14,
                                                              permute = TRUE),
                                 .id = "permutation")

null_untreated_14 <- null_untreated_14 %>% 
  mutate(intervention = "Untreated")

#Obtain new 95% CIs 

bounds_intervention_14 <- null_intervention_14 %>% 
  # generate quantile-based estimates of CI for each distance interval 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))


bounds_untreated_14 <- null_untreated_14 %>% 
  # generate quantile-based estimates of CI for each distance interval 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))

p1_14 <- observed_intervention_14 %>% 
  full_join(bounds_intervention_14) %>% 
  ggplot(aes(x = r_upper, y = tau)) + 
  geom_line() +
  geom_ribbon(aes(ymin = CI.l,
                  ymax = CI.u),
              col = "lightgray",
              alpha = 0.3) +
  geom_hline(yintercept = 1,
             lty = 2) +
  coord_cartesian(ylim = c(0,6)) + 
  theme_minimal() + 
  labs(x = TeX("Distance: $d_2"),
       y = TeX("Odds Ratio: $\\tau"),
       title = TeX("Estimated $\\tau"),
       subtitle = "Intervention arm") +
  theme(legend.title = element_blank())


p2_14 <- observed_untreated_14 %>% 
  full_join(bounds_untreated_14) %>% 
  ggplot(aes(x = r_upper, y = tau)) + 
  geom_line() +
  geom_ribbon(aes(ymin = CI.l,
                  ymax = CI.u),
              col = "lightgray",
              alpha = 0.3) +
  geom_hline(yintercept = 1,
             lty = 2) +
  coord_cartesian(ylim = c(0,6)) + 
  theme_minimal() + 
  labs(x = TeX("Distance: $d_2"),
       y = TeX("Odds Ratio: $\\tau"),
       title = TeX("Estimated $\\tau"),
       subtitle = "Untreated arm") +
  theme(legend.title = element_blank())

ggarrange(p1_14, p2_14)

write.csv(null_intervention_14, "Data/Null_intervention_14.csv")

write.csv(null_untreated_14, "Data/Null_untreated_14.csv")

#####
#And now for 30 days
d1_30 <- c(0,0,0,0,0) # lower bounds distance
d2_30 <- c(100,200,500,750,1000) # upper bounds distance
t_30 <- 30 # time in days

observed_intervention_30 <- permutation_function(
  df = work_ts_int,
  r1 = d1,
  r2 = d2,
  t = t_30,
  permute = FALSE
)

observed_intervention_30 <- observed_intervention_30 %>% 
  mutate(Days = "30")

observed_untreated_30 <- permutation_function(
  df = work_ts_cont,
  r1 = d1,
  r2 = d2,
  t = t_30,
  permute = FALSE
)

observed_untreated_30 <- observed_untreated_30 %>% 
  mutate(Days = "30")

perms <- as.vector(1:1000, mode = "list")

#Obtain null distribution for intervention arm using permutations
null_intervention_30 <- future_map_dfr(perms, ~permutation_function(df = work_ts_int,
                                                                    r1 = d1_30,
                                                                    r2 = d2_30,
                                                                    t = t_30,
                                                                    permute = TRUE),
                                       .id = "permutation")

null_intervention_30 <- null_intervention_30 %>% 
  mutate(intervention = "Intervention")

#Obtain null distribution for control arm using permutations 
null_untreated_30 <- future_map_dfr(perms, ~permutation_function(df = work_ts_cont,
                                                                 r1 = d1_30,
                                                                 r2 = d2_30,
                                                                 t = t_30,
                                                                 permute = TRUE),
                                    .id = "permutation")

null_untreated_30 <- null_untreated_30 %>% 
  mutate(intervention = "Untreated")

#Obtain new 95% CIs 

bounds_intervention_30 <- null_intervention_30 %>% 
  # generate quantile-based estimates of CI for each distance interval 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))


bounds_untreated_30 <- null_untreated_30 %>% 
  # generate quantile-based estimates of CI for each distance interval 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025),
            CI.u = quantile(tau, probs = 0.975))

p1_30 <- observed_intervention_30 %>% 
  full_join(bounds_intervention_30) %>% 
  ggplot(aes(x = r_upper, y = tau)) + 
  geom_line() +
  geom_ribbon(aes(ymin = CI.l,
                  ymax = CI.u),
              col = "lightgray",
              alpha = 0.3) +
  geom_hline(yintercept = 1,
             lty = 2) +
  coord_cartesian(ylim = c(0,6)) + 
  theme_minimal() + 
  labs(x = TeX("Distance: $d_2"),
       y = TeX("Odds Ratio: $\\tau"),
       title = TeX("Estimated $\\tau"),
       subtitle = "Intervention arm") +
  theme(legend.title = element_blank())


p2_30 <- observed_untreated_30 %>% 
  full_join(bounds_untreated_30) %>% 
  ggplot(aes(x = r_upper, y = tau)) + 
  geom_line() +
  geom_ribbon(aes(ymin = CI.l,
                  ymax = CI.u),
              col = "lightgray",
              alpha = 0.3) +
  geom_hline(yintercept = 1,
             lty = 2) +
  coord_cartesian(ylim = c(0,6)) + 
  theme_minimal() + 
  labs(x = TeX("Distance: $d_2"),
       y = TeX("Odds Ratio: $\\tau"),
       title = TeX("Estimated $\\tau"),
       subtitle = "Untreated arm") +
  theme(legend.title = element_blank())

ggarrange(p1_30, p2_30)

write.csv(null_intervention_30, "Data/Null_intervention_30.csv")

write.csv(null_untreated_30, "Data/Null_untreated_30.csv")


#Obtain a dataframe with all observed tau values
Full_observed <- rbind(observed_intervention, observed_intervention_14, 
                       observed_intervention_30, observed_untreated,
                       observed_untreated_14, observed_untreated_30)

write.csv(Full_observed, "Data/Full_observed.csv")

