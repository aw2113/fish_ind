get_mean_fit_ratio <- function(model_fits, ts_dat) {
  
  ## get number of time series
  num_ts <- dim(ts_dat)[1]
  
  # sum of squared residuals
  ssq_resid <- matrix(nrow = num_ts, ncol = 1)
  for(i in 1:num_ts) {
    ssq_resid[i,] <- sum((model_fits$ex[i,] - ts_dat[i,])^2, na.rm = TRUE)
  }
  
  # sum of squared observations
  ssq_obs <- matrix(nrow = num_ts, ncol = 1)
  for(i in 1:num_ts) {
    ssq_obs[i,] <- sum((ts_dat[i,] - 0)^2, na.rm = TRUE)
  }
  
  # fit ratios
  mean_fit_ratio <- mean(ssq_resid/ssq_obs)
  return(mean_fit_ratio)

}

# testing_123 <- get_mean_fit_ratio(mod_fit_ape_eq, dat)

