# Constant Betting Function
BF_constant <- function(p) {
  if (p >= 0 && p < 0.5) {
    val <- 1.5
  } else if(p >=0.5 && p <= 1) {
    val <- 0.5
  } else {
    val <- 0}
}

# Mixture Betting Function
BF_mixture <- function(p) {
  integrand <- function(epsilon) {
    epsilon * p^(epsilon - 1)}
  
  result <- integrate(integrand, lower = 0, upper = 1)
  return(result$value)
}

# Kernel Density Betting Function
make_BF_kde <- function(p_values, n_grid = 512) {
  
  # Create extended sample with reflections
  extended_sample <- c()
  for(p in p_values) {
    extended_sample <- c(extended_sample, -p, p, 2-p)}
  
  # Create grid for density estimation - extended to [-1.5, 2.5] as suggested
  grid <- seq(-1.5, 2.5, length.out = n_grid)
  
  # Calculate bandwidth using Silverman's rule of thumb
  n <- length(p_values)
  sigma <- min(sd(p_values), IQR(p_values)/1.34)
  h <- 0.9 * sigma * n^(-0.2) # Silverman's rule
  
  # Compute KDE using the stats package
  kde <- stats::density(extended_sample, kernel = "gaussian", bw = h, n = n_grid, 
                        from = -1.5, to = 2.5)
  
  # Extract the density values within [0,1]
  in_range <- kde$x >= 0 & kde$x <= 1
  x_vals <- kde$x[in_range]
  y_vals <- kde$y[in_range]
  
  # Normalize to integrate to 1 over [0,1]
  delta <- x_vals[2] - x_vals[1]
  area <- sum(y_vals) * delta
  y_vals <- y_vals / area
  
  # Return a betting function of p 
  betting_function <- function(p) {
    if (any(p < 0 | p > 1)) return(0)
    return(approx(x_vals, y_vals, xout = p, rule = 2)$y)
  }
  
  return(betting_function)
}

# Precomputed Kernel Density Betting Function
make_precomputed_kd_BF <- function(training_set, test_size, non_conformity_measure){
  train_alphas <- numeric(test_size)
  train_base <- training_set[1]
  remaining_train <- training_set[2:test_size]
  
  # Calculate non-conformity measure of the elements in the training set
  for (i in 1:(test_size-1)) {
    train_alphas[i] <- non_conformity_measure(training_set[i], train_base)}
  
  # Calculate the p-values of the elements in the training set
  train_p_values <- numeric(test_size - 1)
  for (i in 1:(test_size-1)) {
    alpha_i <- train_alphas[i]
    U <- runif(1)
    train_p_values[i] <- (sum(train_alphas[1:i] > alpha_i) + U * sum(train_alphas[1:i] == alpha_i)) / i}
  
  # Create the kernel density betting function based on the p-values of the training set
  bf_precomputed_kde <- make_BF_kde(train_p_values)
  
  return(bf_precomputed_kde)
}