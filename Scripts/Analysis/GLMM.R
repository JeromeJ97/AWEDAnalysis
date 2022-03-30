mixed.effects.results <- NULL
gee.results <- NULL
for (i in 1:ncol(exposures)){ # loop through the various exposures, as determined by times and distances
  mixed.effects.model <-  temp %>% 
    full_join(data.frame(participant_id = rownames(exposures), exposure = exposures[,i])) %>% 
    arrange(cluster) %>% 
    glmer(dengue ~ (1 | cluster) + exposure + intervention,
          family = binomial,
          data = .)
  
  mixed.effects.results <- append(mixed.effects.results,
                                  mixed.effects.model)
  names(mixed.effects.results)[i] <- colnames(exposures)[i]}