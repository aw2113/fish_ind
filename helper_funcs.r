## ----dfa-defn-get-DFA-fits---------------------------------------------------------------------
get_DFA_fits <- function(MLEobj, dd = NULL, alpha = 0.05) {
  ## empty list for results
  fits <- list()
  ## extra stuff for var() calcs
  Ey <- MARSS:::MARSShatyt(MLEobj)
  ## model params
  ZZ <- coef(MLEobj, type="matrix")$Z
  ## number of obs ts
  nn <- dim(Ey$ytT)[1]
  ## number of time steps
  TT <- dim(Ey$ytT)[2]
  ## get the inverse of the rotation matrix
  H_inv <- varimax(ZZ)$rotmat
  ## check for covars
  if(!is.null(dd)) {
    DD <- coef(MLEobj, type = "matrix")$D
    ## model expectation
    fits$ex <- ZZ %*% H_inv %*% MLEobj$states + DD %*% dd
  } else {
    ## model expectation
    fits$ex <- ZZ %*% H_inv %*% MLEobj$states
  }
  ## Var in model fits
  VtT <- MARSSkfss(MLEobj)$VtT
  VV <- NULL
  for(tt in 1:TT) {
    RZVZ <- coef(MLEobj, type = "matrix")$R - ZZ %*% VtT[,,tt] %*% t(ZZ)
    SS <- Ey$yxtT[,,tt] - Ey$ytT[,tt,drop = FALSE] %*% t(MLEobj$states[,tt,drop = FALSE])
    VV <- cbind(VV, diag(RZVZ + SS %*% t(ZZ) + ZZ %*% t(SS)))
  }
  SE <- sqrt(VV)
  ## upper & lower (1-alpha)% CI
  fits$up <- qnorm(1-alpha/2)*SE + fits$ex
  fits$lo <- qnorm(alpha/2)*SE + fits$ex
  return(fits)
}


# get DFA fits with only 1 trend --------------------------------------------- #
get_DFA1_fits <- function(MLEobj, dd = NULL, alpha = 0.05) {
  ## empty list for results
  fits <- list()
  ## extra stuff for var() calcs
  Ey <- MARSS:::MARSShatyt(MLEobj)
  ## model params
  ZZ <- coef(MLEobj, type="matrix")$Z
  ## number of obs ts
  nn <- dim(Ey$ytT)[1]
  ## number of time steps
  TT <- dim(Ey$ytT)[2]
  ## get the inverse of the rotation matrix
  # H_inv <- varimax(ZZ)$rotmat
  ## check for covars
  if(!is.null(dd)) {
    DD <- coef(MLEobj, type = "matrix")$D
    ## model expectation
    # fits$ex <- ZZ %*% H_inv %*% MLEobj$states + DD %*% dd
    fits$ex <- ZZ %*% MLEobj$states + DD %*% dd
  } else {
    ## model expectation
    fits$ex <- ZZ %*% MLEobj$states
  }
  ## Var in model fits
  VtT <- MARSSkfss(MLEobj)$VtT
  VV <- NULL
  for(tt in 1:TT) {
    RZVZ <- coef(MLEobj, type = "matrix")$R - ZZ %*% VtT[,,tt] %*% t(ZZ)
    SS <- Ey$yxtT[,,tt] - Ey$ytT[,tt,drop = FALSE] %*% t(MLEobj$states[,tt,drop = FALSE])
    VV <- cbind(VV, diag(RZVZ + SS %*% t(ZZ) + ZZ %*% t(SS)))
  }
  SE <- sqrt(VV)
  ## upper & lower (1-alpha)% CI
  fits$up <- qnorm(1-alpha/2)*SE + fits$ex
  fits$lo <- qnorm(alpha/2)*SE + fits$ex
  return(fits)
}


# get the mean fit ratio ----------------------------------------------------- #
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
