# Distance to the mean
make_NCM_mean <- function() {
  return(function(z, training_set, ...) {
    abs(z - mean(training_set))
  })
}

# k Nearest Neighbors
make_NCM_knn <- function(k = 5) {
  return(function(z, training_set, ...) {
    
    if (length(training_set) == 0) {
      return(0)}
    
    training_set <- if (is.vector(training_set)) matrix(training_set, ncol = 1) else training_set
    z <- matrix(z, nrow = 1)
    data_with_z <- rbind(training_set, z)
    
    knn_result <- FNN::get.knn(data_with_z, k = k)
    distances_z <- knn_result$nn.dist[nrow(knn_result$nn.dist), ]
    return(mean(distances_z))
  })
}

# Likelihood ratio for normal distribution with fixed mu_r
make_NCN_lr <- function(mu_r = 1, sigma2 = 1, sigma2_r = 1) {
  return(function(z, training_set, ..., n = NULL, new_data = NULL) {
    mu_hat_0 <- mean(training_set)
    
    # Likelihood under the prior predictive
    f1 <- dnorm(z, mean = mu_r, sd = sqrt(sigma2 + sigma2_r))
    
    # Likelihood under the reference (pre-change) distribution
    f0 <- dnorm(z, mean = mu_hat_0, sd = sqrt(sigma2))
    
    return(f1/f0)
  })
}

# Likelihood ratio for normal distribution with rolling window estimation
make_NCN_lr_rolling <- function(window_size = 20) {
  return(function(z, training_set, n, new_data, ...) {
    
    # Estimate reference (pre-change) parameters
    mu_0 <- mean(training_set)
    sigma2 <- var(training_set)
    
    # Compute the NCM of the training set with fixed post-change distribution
    if (n < window_size){
      window_data <- new_data[1:n]
      sigma2_r <- 1
      mu_r <- 1
    } else{
      window_data <- tail(new_data[1:n], window_size) # Extract rolling window of recent observations
      mu_r <- mean(window_data)
      sigma2_r <- var(window_data)
      sigma2_r <- max(sigma2_r, 1e-4)  # avoid zero variance
    }
    
    # Compute likelihoods
    f1 <- dnorm(z, mean = mu_r, sd = sqrt(sigma2 + sigma2_r))
    f0 <- dnorm(z, mean = mu_0, sd = sqrt(sigma2))
    
    return(f1/f0)
  })
}