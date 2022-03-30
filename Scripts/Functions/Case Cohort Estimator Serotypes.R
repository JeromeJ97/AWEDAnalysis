#Try to modify case-cohort estimator such that it partitions vcd cases by serotype 

case_cohort_function <- function(dataset, d1, d2, time){
  
  temp <- dataset %>% 
    # CRITICAL - put the individuals in time order according to date of illness onset
    arrange(illness_onset)
  
  # Identify the test-positive dengue cases, by serotype
  vcd_s1 <- temp$participant_id[temp$serotype == "denv1"]
  vcd_s2 <- temp$participant_id[temp$serotype == "denv2"]
  vcd_s3 <- temp$participant_id[temp$serotype == "denv3"]
  vcd_s4 <- temp$participant_id[temp$serotype == "denv4"]
  
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
  
  
  overall <- list(contingency.table = contingency.table,
                  mixed.effects.results = mixed.effects.results,
                  gee.results = gee.results)
  
  return(overall)
}