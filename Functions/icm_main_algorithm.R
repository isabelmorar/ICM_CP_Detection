icm_changepoint_detection <- function(training_set, new_data, threshold = 2.0, non_conformity_measure, betting_function){
  
  # Shuffle the training set to induce exchangeability
  shuffled_training_set <- sample(training_set)
  
  # Vector to store non_conformity measure in every iteration
  alpha_vals <- numeric(length(data))
  
  # Initialize the martingale and modified martingale
  S <- 1
  C <- 0
  C_og <- 0 
  
  # Store results for later analysis
  results <- data.frame(Observation = numeric(0), 
                        Martingale_S = numeric(0), 
                        Modified_C = numeric(0), 
                        Detected_Change = logical(0))
  
  first_detection_point <- 0
  
  for (n in 1:length(new_data)) {
    z_n <- new_data[n] 
    
    # Step 1: Calculate non-conformity score for the new observation z_n
    alpha_n <- non_conformity_measure(z_n, shuffled_training_set, n = n, new_data = new_data)
    alpha_vals[n] <- alpha_n
    
    # Step 2: Calculate p-value for z_n using precomputed non-conformity scores
    U <- runif(1) 
    p_n <- (sum(alpha_vals > alpha_n) + U * sum(alpha_vals == alpha_n)) / n
    
    if (p_n == 0){
      p_n = 0.001
    }
    
    # Step 3: Update the martingale value S_n using the betting function
    g_n <- betting_function(p_n)
    S <- S * g_n
    
    # Step 4: Update the modified martingale value C_n
    C <- max(0, C + log(g_n))
    
    vector_S <- c(results$Martingale_S, S)
    C_og <- log(S) - min(vector_S)
    
    # Step 5: Check if a change-point is detected
    detected_change <- C >= threshold
    
    # Store results
    results <- rbind(results, 
                     data.frame(Observation = z_n, Martingale_S = S, 
                                Modified_C = C, Modified_C_og = C_og, Detected_Change = detected_change))
  }
  return(results)
}