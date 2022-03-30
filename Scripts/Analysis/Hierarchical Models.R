#######
#Analysis results from running GLMM and GEE model
source("DM/Case Cohort Packages.R")

library(dplyr)
library(geepack)
library(ggpubr)
library(here)
library(tidyverse)
library(geodist)
library(knitr)
library(kableExtra)
library(lme4)

load("Data/2020-12-08_work-ts-spat.rdata")

#Load the function required for analysis 
source("Scripts/Functions/Case Cohort Estimator.R")


#Setting up the time and distance constraints
# time
t <- c(7,14,30) # MODIFY THIS
# distances 
r.upper <- c(100,200,500,750,1000) # MODIFY THIS
r.lower <- c(0,0,0,0,0) # MODIFY THIS - must be the same length as r.upper

#Run the function
output_all_t <- case_cohort_function(dataset = work_ts_spat,
                                     d1 = r.lower,
                                     d2 = r.upper,
                                     time = t)

#######
#Contingency Table results
day_label <- c("7 Days", "14 Days", "30 Days")
names(day_label) <- c("7", "14", "30")
output_all_t$contingency.table %>% 
  ungroup() %>% 
  filter(dengue == 1) %>% 
  mutate(days = as.factor(time)) %>% 
  ggplot(aes(x = distance2, y = rr, col = days, fill = days)) +
  facet_wrap(~days, ncol=1, labeller= labeller(days = day_label)) +
  geom_line() + 
  geom_ribbon(aes(ymin = CI.l, ymax = CI.u),
              alpha = 0.3,
              col = "white") +
  geom_hline(yintercept = 1, lty = 2) + 
  ggpubr::theme_pubr() +
  xlab("Distance from a VCD (m)") +
  ylab("Case-Cohort Odds Ratio") + 
  theme(legend.position = "none")

output_all_t$contingency.table %>% 
  ungroup() %>% 
  filter(dengue == 1) %>% 
  mutate(days = as.factor(time)) %>% 
  dplyr::select(distance = distance2, days, or = rr, CI.l, CI.u) %>% 
  kable(digits = 2) %>% 
  kable_styling() %>% 
  collapse_rows(columns = 2)

#######
#GLMM Results
mixed_effects_estimates <- map(output_all_t$mixed.effects.results, ~summary(.)$coefficients) %>% 
  map_dfr(~.["exposure",1:2], # extract the estimated log or (column 1) and estimated Standard Error (column 2) for the exposure variable
          .id = "distance.time") %>% 
  mutate(or = exp(Estimate),
         CI.l = exp(Estimate - 1.96*`Std. Error`),
         CI.u = exp(Estimate + 1.96*`Std. Error`),
         distance = as.numeric(str_match(distance.time, "d\\s*(.*?)\\.t")[,2]),
         days = str_split_fixed(distance.time, pattern = ".t", n = 2)[,2]) %>% 
  mutate(days = factor(days, levels = c("7", "14", "30"))) 

mixed_effects_estimates %>% 
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

mixed_effects_estimates %>% 
  dplyr::select(distance, days, or, CI.l, CI.u) %>% 
  kable(digits = 2) %>% 
  kable_styling() %>% 
  collapse_rows(columns = 2)

##########
#GEE Results
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

###########
#Obtain side by side estimates
cont.results <- output_all_t$contingency.table %>% 
  ungroup() %>% 
  filter(dengue == 1) %>% 
  dplyr::select(distance = distance2,
                days = time,
                or = rr, 
                CI.l, CI.u) %>% 
  mutate(days = as.factor(days),
         type = "contingency table") 

me.results <- mixed_effects_estimates %>% 
  dplyr::select(-distance.time, -Estimate, - `Std. Error`) %>% 
  mutate(type = "mixed effects")

gee.results <- gee.estimates %>% 
  dplyr::select(or, CI.l, CI.u, distance, days) %>% 
  mutate(type = "GEE")


cont.results %>% 
  bind_rows(me.results) %>% 
  bind_rows(gee.results) %>% 
  ggplot(aes(x = distance, y = or, col = type, fill = type)) +
  facet_wrap(~days, ncol = 1) +
  geom_line(aes(lty = type)) + 
  # geom_ribbon(aes(ymin = CI.l, ymax = CI.u),
  #             alpha = 0.3,
  #             col = "white") +
  geom_hline(yintercept = 1, lty = 2) + 
  ggpubr::theme_pubr() +
  xlab("Distance from a VCD (m)") +
  ylab("Case-Cohort Odds Ratio") +
  theme(legend.position = "bottom")


