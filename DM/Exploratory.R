######## 
#Exploratory analysis of the data

#Load the dataset 
load("Data/2020-12-08_work-ts-spat.rdata")

install.packages("dplyr")
library(dplyr)
#Total no. of people with and without dengue in sample
table(work_ts_spat$dengue)

#Specifically serotype
table(work_ts_spat$serotype)

#General summary 
summary(work_ts_spat)

#Identify clusters with the most and least dengue cases
#Need to cumulate case numbers by clusters 
cases_by_cluster<- work_ts_spat %>%
  group_by(cluster) %>%
  tally(dengue)

#Cluster 13 had the most dengue cases (76) and cluster 24 had the least (2)
#Also, clusters 2, 16 and 21 only had 3 cases each, they are clusters that are on the edge
#Cluster 13 is landlocked 
#Cluster with second highest cases is 18 with 49 cases, much lower than 13 


#Identify the number of people in each cluster 
cluster_total<- work_ts_spat %>%
  count(cluster)

#This may be used in a table, but too many observations and clusters aren't that 
#meaningful... but cluster 13 has most inhabitants (794) and cluster 17 has the 
#fewest (89), though clusters 21 and 2 have 90 and 91 inhabitants respectively 
#Cluster 6 had second highest number of inhabitants with 666 people in the study


#Want to find the average age of dengue cases and test-negatives
#And overall average age 
case_age <- work_ts_spat %>%
  filter(dengue==1) %>%
  summarise(mean(age))

control_age<- work_ts_spat %>%
  filter(dengue==0) %>%
  summarise(mean(age))

#Mean age for dengue cases is 13 y/o whereas mean age for controls is ~15 y/o
#Not much difference

#Group age by cluster
cluster_age <- work_ts_spat %>%
  group_by(cluster) %>%
  summarise(mean = mean(age), n = n(), sd = sd(age))

#obtain boxplot of age by cluster 
boxplot(age~cluster, data=work_ts_spat, xlab="Cluster", ylab= "Age (years)")

#Try to stratify mean ages of cases by serotype
sero_age<- work_ts_spat %>%
  group_by(serotype) %>%
  summarise(mean(age))

#Doesn't say much - mean age of DENV1 is 11, DENV2 12.6, DENV3 15.6 and DENV4 14.1

#Identify the sex with more cases 
cases_sex<- work_ts_spat %>%
  group_by(sex) %>%
  tally(dengue)

#Roughly same number of cases between males and females: 201 to 191 respectively(?)

serotype_num <- work_ts_spat %>%
  group_by(serotype) %>%
  tally(dengue)
#Stratify sex by cluster
cluster_sex <- work_ts_spat %>%
  group_by(cluster) %>%
  tally(sex)






