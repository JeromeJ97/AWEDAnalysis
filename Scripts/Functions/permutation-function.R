permutation_function <- function(df, r1, r2, t, permute = FALSE){
  # This function permutes individual locations to break the 
  # relationship between location and relatedness (e.g. serotype, time)
  
  if (permute == TRUE){
    # Pulling out location data
    location <- df %>% 
      dplyr::select(participant_id, latitude, longitude) %>% 
      distinct() 
    
    # Shuffling participant location
    location$shuffled_id <- location$participant_id[shuffle(location$participant_id)] 
    location <- location %>% 
      dplyr::select(-participant_id)
    
    temp <- df %>% 
      dplyr::select(-latitude, -longitude) %>% 
      full_join(location, by = c("participant_id" = "shuffled_id"))  
    
  } else {
    temp <- df
  }
  
  ###############
  # MATRICES
  ###############
  # First, estimate pairwise distances in METERS
  dist.mat <- geodist(x = cbind(longitude = temp$longitude, latitude = temp$latitude))
  # Set up for difference in serotype and time (absolute time!)
  time.pairs <- expand.grid(temp$illness_onset, temp$illness_onset)
  time.dif <- as.Date(time.pairs[,1]) - as.Date(time.pairs[,2])
  time.mat <- abs(matrix(time.dif, nrow = nrow(temp)))
  # Serotype
  sero.pairs <- expand.grid(temp$serotype, temp$serotype)
  sero.related <- 1*(sero.pairs[,1] != "test-negative control")*(sero.pairs[,1] == sero.pairs[,2])
  sero.mat <- matrix(sero.related, nrow = nrow(temp))
  
  # Remove self-comparisons
  diag(sero.mat) <- diag(time.mat) <- diag(dist.mat) <- NA
  
  
  ###############
  # ESTIMATION
  ###############
  
  tau_overall <- map2_dfr(r1, # x 
                          r2, # y
                          ~bin.fun(r1 = .x, 
                                   r2 = .y,
                                   t = t,
                                   A = time.mat,
                                   B = dist.mat,
                                   C = sero.mat))
  
  return(tau_overall)
}
