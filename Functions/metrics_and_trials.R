# Delay (time until change-point detection)
delay <- function(results, change_point_index) {
  detected_change_points <- which(results$Detected_Change == TRUE)
  
  if (length(detected_change_points) > 0) {
    first_detection_point <- detected_change_points[1]
    
    # Ensure we are considering the first detection after the change-point (i.e., when T > theta)
    if (first_detection_point > change_point_index) {
      delay <- first_detection_point - change_point_index
      return(delay)  # Return the delay in terms of observations
    } else {
      return(NA)  # No detection after the change-point
    }
  } else {
    return(NA)  # No change-point detected
  }
}


false_alarm <- function(results, change_point_index) {
  # Check if a change-point was detected before the actual change-point
  detected_indices <- which(results$Detected_Change == TRUE)
  any_detected_before_theta <- any(detected_indices < change_point_index)
  
  if (any_detected_before_theta) {
    return(1)  # False alarm detected
  } else {
    return(0)  # No false alarm
  }
}


conduct_trials <- function(theta, mu1, test_size = 200, n = 1000, num_trials = 100, betting_function, 
                           non_conformity_measure, th = 2, precomputed_kernel = FALSE){
  vector_delays = c()
  false_alarm_count = 0 
  
  for (trial in 1:num_trials) {
    data_before <- rnorm(theta + test_size, mean = 0, sd = 1)
    data_after <- rnorm(n - theta, mean = mu1, sd = 1)
    data_norm <- c(data_before, data_after)
    
    upper <- test_size-1
    lower <- test_size
    final <- n + test_size
    
    if (precomputed_kernel == TRUE){
      training_set = data_norm[1:upper]
      betting_function <- make_precomputed_kd_BF(training_set, test_size, non_conformity_measure)}
    
    results <- icm_changepoint_detection(training_set = data_norm[1:upper], new_data = data_norm[lower:final], 
                                         threshold = th, non_conformity_measure, betting_function)
    
    delay_trial <- delay(results, theta)
    prob <- false_alarm(results, theta)
    
    false_alarm_count = false_alarm_count + prob
    vector_delays = c(vector_delays, delay_trial)
  }
  
  E_mean_delay = mean(na.omit(vector_delays))
  log_mean_delay = log10(1+E_mean_delay)
  prob_false_alarm = false_alarm_count/num_trials
  
  results_validation = data.frame("Threshold" = th, "Mean Delay" = E_mean_delay, 
                                  log_E_md = log_mean_delay, "Prob False Alarm" = prob_false_alarm)
  
  return(results_validation)
}


run_multiple_threshold_trials <- function(threshold_values, theta, mu1, test_size = 200, n = 1000, num_trials = 100, 
                                          betting_function, non_conformity_measure, precomputed_kernel = FALSE) {
  
  all_results <- lapply(threshold_values, function(th) {
    conduct_trials(theta = theta, mu1 = mu1, test_size = test_size, n = n, 
                   num_trials = num_trials, betting_function = betting_function, 
                   non_conformity_measure = non_conformity_measure, th = th, precomputed_kernel = precomputed_kernel)
  })
  
  final_results <- do.call(rbind, all_results)
  
  return(final_results)
}


evaluate_tidychangepoint_method <- function(theta, mu1, method, n = 1000, num_trials = 100, max_inter_ga = 5) {
  
  vector_delays <- c()
  false_alarm_count <- 0
  detected_count <- 0
  
  for (trial in 1:num_trials) {
    data_before <- rnorm(theta, mean = 0, sd = 1)
    data_after <- rnorm(n - theta, mean = mu1, sd = 1)
    data_norm <- c(data_before, data_after)
    
    # Run the selected changepoint method
    if (method %in% c("ga", "ga-shi", "ga-coen")){
      result = segment(data_norm, method = method, maxiter = max_inter_ga)
    } else {
      result = segment(data_norm, method = method)
    }
    
    # Extract changepoints
    change_point <- as.vector(changepoints(result)[1])
    
    # Calculate mean delay 
    if (!is.na(change_point)) {
      detected_count <- detected_count + 1
      # Ensure we are considering the first detection after the change-point (i.e., when T > theta)
      if (change_point > theta){
        delay <- change_point - theta
        vector_delays <- c(vector_delays, delay)} 
    }
    
    # Calculate probability of false alarm
    if (!is.na(change_point)) {
      # Check if a change-point was detected before the actual change-point
      if (change_point < theta) {
        false_alarm_count <- false_alarm_count + 1}
    }
  }
  
  E_mean_delay = mean(na.omit(vector_delays))
  log_mean_delay = log10(1+E_mean_delay)
  
  if (detected_count > 0){
    prob_false_alarm = false_alarm_count/num_trials
  } else {prob_false_alarm = NA}
  
  results_validation = data.frame("Mean Delay" = E_mean_delay, 
                                  log_E_md = log_mean_delay, "Prob False Alarm" = prob_false_alarm)
  
  return(results_validation)
}