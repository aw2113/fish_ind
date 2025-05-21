library('MARSS')
library('viridis')
library('Hmisc')
source('helper_funcs.r')
# based on Fish 550 lecture 10: DFA

# Apex predators
apex_dat <- read.csv("data/apex_2023.csv")
# They "demean" in their example as opposed to z-score but I'm pretty sure I need
# divide by the sd in addition to demean-ing, hence scale()
apex_wide <- scale(apex_dat[,2:10])

## ----dfa-trans-data----------------------------------------------------------------------------
## transpose data so time goes across columns
dat <- t(apex_wide)
colnames(dat) <- apex_dat[,1]
## get number of time series
N_ts <- dim(dat)[1]
## get length of time series
TT <- dim(dat)[2] 
rownames(dat)
apex_predators <- rownames(dat)

## ----dfa-plot-phytos, fig.height=9, fig.width=8, fig.cap='Demeaned time series of Lake Washington phytoplankton.'----
spp <- rownames(dat)
clr <- viridis(dim(dat)[1]+1)
cnt <- 1
par(mfrow = c(3,3), mai = c(0.5,0.7,0.2,0.1), omi = c(0,0,0,0))
for(i in spp){
  plot(1982:2023, dat[i,], xlab = "",ylab="biomass", bty = "L", xaxt = "n", 
       pch=16, col=clr[cnt], type="b")
  axis(1, labels=TRUE)#1:dim(dat)[2]+1, 1982:2023)
  minor.tick(nx = 5, tick.ratio=0.3)
  title(i)
  cnt <- cnt + 1
}

## -------------------------------------------------------------------------- ##
## Covariates
# DFA with covariate

# coldpool data
coldpool_data <- read.csv("data/coldpool.csv")
cpz <- scale(coldpool_data[,2:12])
# cp_area <- t(cpz[,"AREA_LTE2_KM2"])
npi <- t(cpz[,"npi"])
ice_ext <- t(cpz[,"ice_ext"])
totc <- t(cpz[,"totc"])
# pdo <- t(cpz[,"pdo"])
# mgt <- t(cpz[,"MEAN_GEAR_TEMPERATURE"])
# mst <- t(cpz[,"MEAN_SURFACE_TEMPERATURE"])
# cp_area <- t(coldpool_data$AREA_LTE2_KM2)
# npi <- t(coldpool_data$npi)
# ice_ext <- t(coldpool_data$ice_ext)
# pdo <- t(coldpool_data$pdo)


## ----dfa-dfa-obs-eqn---------------------------------------------------------------------------
## 'ZZ' is loadings matrix

# ----- ----- ***** 3 trends ***** ----- ----- #
# 9 time series
# Z_vals <- list("z11",  0  ,  0  ,
#                "z21","z22",  0  ,
#                "z31","z32","z33",
#                "z41","z42","z43",
#                "z51","z52","z53",
#                "z61","z62","z63",
#                "z71","z72","z73",
#                "z81","z82","z83",
#                "z91","z92","z93")
# ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 3, byrow = TRUE)
# ZZ
# ## 'aa' is the offset/scaling
# aa <- "zero"
# ## 'DD' and 'd' are for covariates
# DD <- "unconstrained"  # unconstrained when we have covariates
# dd <- rbind(npi, ice_ext, totc)  # covariate time series
# ## 'RR' is var-cov matrix for obs errors
# RR <- "diagonal and equal"
# RRdu <- "diagonal and unequal"
# RReq <- "equalvarcov"
# RRuneq <- "unconstrained"

## ----dfa-dfa-proc-eqn--------------------------------------------------------------------------
## number of processes
# mm <- 3
# ## 'BB' is identity: 1's along the diagonal & 0's elsewhere
# BB <- "identity"  # diag(mm)
# ## 'uu' is a column vector of 0's
# uu <- "zero"  # matrix(0, mm, 1)
# ## 'CC' and 'cc' are for covariates
# CC <- "zero"  # matrix(0, mm, 1)
# cc <- "zero"  # matrix(0, 1, wk_last)
# ## 'QQ' is identity
# QQ <- "identity"  # diag(mm)

## ----dfa-create-model-lists--------------------------------------------------------------------
## list with specifications for model vectors/matrices
# mod_list <- list(Z = ZZ, A = aa, D = DD, d = dd, R = RR,
#                  B = BB, U = uu, C = CC, c = cc, Q = QQ)
# mod_list_du <- list(Z = ZZ, A = aa, D = DD, d = dd, R = RRdu,
#                     B = BB, U = uu, C = CC, c = cc, Q = QQ)
# mod_list_eq <- list(Z = ZZ, A = aa, D = DD, d = dd, R = RReq,
#                     B = BB, U = uu, C = CC, c = cc, Q = QQ)
# mod_list_uneq <- list(Z = ZZ, A = aa, D = DD, d = dd, R = RRuneq,
#                       B = BB, U = uu, C = CC, c = cc, Q = QQ)
## list with model inits
# init_list <- list(x0 = matrix(rep(0, mm), mm, 1))
## list with model control parameters
con_list <- list(maxit = 10000, allow.degen = TRUE)#, abstol=0.001)

## ----dfa-fit-dfa-1, cache=TRUE-----------------------------------------------------------------
## fit MARSS
ptm <- proc.time()
# if you use form="dfa" or MARSS.dfa(), then the mod_list cannot contain Z, e.g.,
# determine # of trends and include all covariates and fix R at default
mod_list = list(m = 4, R = "diagonal and equal")
ape_dfa4_cov123_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                          control = con_list, covariates = rbind(npi, ice_ext, totc))
mod_list = list(m = 3, R = "diagonal and equal")
ape_dfa3_cov123_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                       control = con_list, covariates = rbind(npi, ice_ext, totc))
mod_list = list(m = 2, R = "diagonal and equal")
ape_dfa2_cov123_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                          control = con_list, covariates = rbind(npi, ice_ext, totc))
mod_list = list(m = 1, R = "diagonal and equal")
ape_dfa1_cov123_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                          control = con_list, covariates = rbind(npi, ice_ext, totc))
# MARSSparamCIs(ape_dfa3_cov123_de)
# Compare AICc and logLik
print(cbind(model=c(4, 3, 2, 1),
            logLik=round(c(ape_dfa4_cov123_de$logLik, ape_dfa3_cov123_de$logLik, 
                           ape_dfa2_cov123_de$logLik, ape_dfa1_cov123_de$logLik),3),
            AICc=round(c(ape_dfa4_cov123_de$AICc, ape_dfa3_cov123_de$AICc, 
                         ape_dfa2_cov123_de$AICc, ape_dfa1_cov123_de$AICc),3),
            convergence=c(ape_dfa4_cov123_de$convergence, ape_dfa3_cov123_de$convergence, 
                          ape_dfa2_cov123_de$convergence, ape_dfa1_cov123_de$convergence)),
      quote=FALSE)
# proc.time() - ptm
# The most supported model has 4 trends
# now evaluate the most supported form of R
mod_list = list(m = 4, R = "diagonal and unequal")
ape_dfa4_cov123_du <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                          control = con_list, covariates = rbind(npi, ice_ext, totc))
mod_list = list(m = 4, R = "equalvarcov")
ape_dfa4_cov123_eq <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                          control = con_list, covariates = rbind(npi, ice_ext, totc))
mod_list = list(m = 4, R = "unconstrained")
ape_dfa4_cov123_un <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                          control = con_list, covariates = rbind(npi, ice_ext, totc))
# unconstrained did not converge
# Compare AICc and logLik
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            logLik=round(c(ape_dfa4_cov123_de$logLik, ape_dfa4_cov123_du$logLik, 
                           ape_dfa4_cov123_eq$logLik, ape_dfa4_cov123_un$logLik),3),
            AICc=round(c(ape_dfa4_cov123_de$AICc, ape_dfa4_cov123_du$AICc, 
                         ape_dfa4_cov123_eq$AICc, ape_dfa4_cov123_un$AICc),3),
            convergence=c(ape_dfa4_cov123_de$convergence, ape_dfa4_cov123_du$convergence, 
                          ape_dfa4_cov123_eq$convergence, ape_dfa4_cov123_un$convergence)),
quote=FALSE)
# the trends are in $states and the SE's are in $states.se

proc.time() - ptm

## get model fits & CI's
mod_fit_ape_de <- get_DFA_fits(ape_dfa4_cov123_de)
mod_fit_ape_eq <- get_DFA_fits(ape_dfa4_cov123_eq)
# # fit ratios
ape_dfa4_cov123_de_fit_ratio <- get_mean_fit_ratio(mod_fit_ape_de, dat)
ape_dfa4_cov123_eq_fit_ratio <- get_mean_fit_ratio(mod_fit_ape_eq, dat)

# ---------------------------------------------------------------------------- #


# Select covariates
# m = 3, R = "diagonal and equal"

ptm <- proc.time()
# all three covariates
mod_list = list(m = 4, R = "diagonal and equal")
# ape_dfa4_cov123_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
#                             control = con_list, covariates = rbind(npi, ice_ext, totc))
ape_dfa4_cov12_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                           control = con_list, covariates = rbind(npi, ice_ext))
ape_dfa4_cov13_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                           control = con_list, covariates = rbind(npi, totc))
ape_dfa4_cov23_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                           control = con_list, covariates = rbind(ice_ext, totc))
ape_dfa4_cov1_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                          control = con_list, covariates = npi)
ape_dfa4_cov2_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                          control = con_list, covariates = ice_ext)
ape_dfa4_cov3_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                          control = con_list, covariates = totc)
ape_dfa4_cov0_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                          control = con_list)

proc.time() - ptm

# get model fits
ape_dfa4_cov123_de_mod_fits <- get_DFA_fits(ape_dfa4_cov123_de)
ape_dfa4_cov12_de_mod_fits <- get_DFA_fits(ape_dfa4_cov12_de)
ape_dfa4_cov13_de_mod_fits <- get_DFA_fits(ape_dfa4_cov13_de)
ape_dfa4_cov23_de_mod_fits <- get_DFA_fits(ape_dfa4_cov23_de)
ape_dfa4_cov1_de_mod_fits <- get_DFA_fits(ape_dfa4_cov1_de)
ape_dfa4_cov2_de_mod_fits <- get_DFA_fits(ape_dfa4_cov2_de)
ape_dfa4_cov3_de_mod_fits <- get_DFA_fits(ape_dfa4_cov3_de)
ape_dfa4_cov0_de_mod_fits <- get_DFA_fits(ape_dfa4_cov0_de)
# get mean fit ratios
ape_dfa4_cov123_de_mean_fit <- get_mean_fit_ratio(ape_dfa4_cov123_de_mod_fits, dat)
ape_dfa4_cov12_de_mean_fit <- get_mean_fit_ratio(ape_dfa4_cov12_de_mod_fits, dat)
ape_dfa4_cov13_de_mean_fit <- get_mean_fit_ratio(ape_dfa4_cov13_de_mod_fits, dat)
ape_dfa4_cov23_de_mean_fit <- get_mean_fit_ratio(ape_dfa4_cov23_de_mod_fits, dat)
ape_dfa4_cov1_de_mean_fit <- get_mean_fit_ratio(ape_dfa4_cov1_de_mod_fits, dat)
ape_dfa4_cov2_de_mean_fit <- get_mean_fit_ratio(ape_dfa4_cov2_de_mod_fits, dat)
ape_dfa4_cov3_de_mean_fit <- get_mean_fit_ratio(ape_dfa4_cov3_de_mod_fits, dat)
ape_dfa4_cov0_de_mean_fit <- get_mean_fit_ratio(ape_dfa4_cov0_de_mod_fits, dat)

# Compare AICc, logLik, and mean fit ratio
print(cbind(R_structure=rep("diagonal and equal", 8),
            covars=c(123, 12, 13, 23, 1, 2, 3, 0),
            logLik=round(c(ape_dfa4_cov123_de$logLik, 
                           ape_dfa4_cov12_de$logLik, 
                           ape_dfa4_cov13_de$logLik, 
                           ape_dfa4_cov23_de$logLik,
                           ape_dfa4_cov1_de$logLik,
                           ape_dfa4_cov2_de$logLik,
                           ape_dfa4_cov3_de$logLik,
                           ape_dfa4_cov0_de$logLik),3),
            AICc=round(c(ape_dfa4_cov123_de$AICc, 
                         ape_dfa4_cov12_de$AICc, 
                         ape_dfa4_cov13_de$AICc, 
                         ape_dfa4_cov23_de$AICc,
                         ape_dfa4_cov1_de$AICc,
                         ape_dfa4_cov2_de$AICc,
                         ape_dfa4_cov3_de$AICc,
                         ape_dfa4_cov0_de$AICc),3),
            mean_fit_ratio=round(c(ape_dfa4_cov123_de_mean_fit,
                                   ape_dfa4_cov12_de_mean_fit,
                                   ape_dfa4_cov13_de_mean_fit,
                                   ape_dfa4_cov23_de_mean_fit,
                                   ape_dfa4_cov1_de_mean_fit,
                                   ape_dfa4_cov2_de_mean_fit,
                                   ape_dfa4_cov3_de_mean_fit,
                                   ape_dfa4_cov0_de_mean_fit),3),
            convergence=c(ape_dfa4_cov123_de$convergence, 
                          ape_dfa4_cov12_de$convergence, 
                          ape_dfa4_cov13_de$convergence, 
                          ape_dfa4_cov23_de$convergence,
                          ape_dfa4_cov1_de$convergence,
                          ape_dfa4_cov2_de$convergence,
                          ape_dfa4_cov3_de$convergence,
                          ape_dfa4_cov0_de$convergence)),
      quote=FALSE)

## get model fits & CI's
mod_fit_ape4_cov3_de <- get_DFA_fits(ape_dfa4_cov3_de)
mod_fit_ape4_cov2_de <- get_DFA_fits(ape_dfa4_cov2_de)
# # fit ratios
ape_dfa4_cov3_de_fit_ratio <- get_mean_fit_ratio(mod_fit_ape4_cov3_de, dat)
ape_dfa4_cov2_eq_fit_ratio <- get_mean_fit_ratio(mod_fit_ape4_cov2_de, dat)

# ============================================================================ #

# The best model [model with the most support]
# 4 trends
# R = "diagonal and equal"
# includes ice_ext as covariate (cov2)

# *****Plot***** ------------------------------------------------------------- #
## ----dfa-get-H-inv-----------------------------------------------------------------------------
mm <- 4
##### ----- ***** diagonal and equal ***** ----- #####
## get the estimated ZZ
Z_est <- coef(ape_dfa4_cov2_de, type = "matrix")$Z
## get the inverse of the rotation matrix
H_inv <- varimax(Z_est)$rotmat

## ----dfa-rotate-Z-x----------------------------------------------------------------------------
## rotate factor loadings
Z_rot = Z_est %*% H_inv   
## rotate processes
proc_rot = solve(H_inv) %*% ape_dfa4_cov2_de$states

## ----dfa-plot-dfa1, fig.height=9, fig.width=8, eval=TRUE, fig.cap='Estimated states from the DFA model.'----
# portrait 6.5 x 8
ylbl <- apex_predators
w_ts <- seq(dim(dat)[2])
# the first four trends
layout(matrix(c(1,2,3,4,5,6,7,8), mm, 2), widths = c(2,1))
## par(mfcol=c(mm,2), mai = c(0.5,0.5,0.5,0.1), omi = c(0,0,0,0))
par(mai = c(0.5,0.5,0.5,0.1), omi = c(0,0,0,0))
## plot the processes
for(i in 1:mm) {
  ylm <- c(-1,1)*max(abs(proc_rot[i,]))
  ## set up plot area
  plot(w_ts,proc_rot[i,], type = "n", bty = "L",
       ylim = ylm, xlab = "", ylab = "", xaxt = "n", cex=1)
  ## draw zero-line
  abline(h=0, col="gray")
  ## plot trend line
  lines(w_ts, proc_rot[i,], lwd = 2)
  lines(w_ts, proc_rot[i,], lwd = 2)
  ## add panel labels
  mtext(paste("State",i," apex"," diag_equal"), side = 3, line = 0.5, cex = 0.75)
  # axis(1,12*(0:dim(dat_1980)[2])+1,yr_frst+0:dim(dat_1980)[2])
  axis(1,seq(1,42,10),as.character(seq(1982,2023,10)))
  # axis(1, labels = TRUE)
  # minor.tick(nx = 5, tick.ratio=0.3)
  
}
## plot the loadings
minZ <- 0
ylm <- c(-1,1)*max(abs(Z_rot))
for(i in 1:mm) {
  plot(c(1:N_ts)[abs(Z_rot[,i])>minZ], as.vector(Z_rot[abs(Z_rot[,i])>minZ,i]), 
       type="h", lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylm, 
       xlim=c(0.5,N_ts+0.5), 
       col=clr, cex=1, bty="n")
  abline(h=c(-0.2,0.2), lty=2, col="gray", lwd=1)
  abline(h=0, lwd=1.1, col="gray")
  for(j in 1:N_ts) {
    if(Z_rot[j,i] > minZ) {text(j, -0.03, ylbl[j], srt=90, adj=1, cex=0.75, col=clr[j])}
    if(Z_rot[j,i] < -minZ) {text(j, 0.03, ylbl[j], srt=90, adj=0, cex=0.75, col=clr[j])}
  } 
  mtext(paste("Factor loadings"),side=3,line=0.5, cex=0.75)
}

## ----dfa-xy-states12, height=4, width=5, fig.cap='Cross-correlation plot of the two rotations.'----
par(mai = c(0.9,0.9,0.1,0.1))
ccf(proc_rot[1,],proc_rot[2,], lag.max = 12, main="")
acf(proc_rot[1,], type = "covariance")



## ----dfa-plot-dfa-fits, fig.height=9, fig.width=8, fig.cap='Data and fits from the DFA model.'----
## get model fits & CI's
mod_fit_ape4_cov2_de <- get_DFA_fits(ape_dfa4_cov2_de)
## plot the fits
ylbl <- apex_predators
par(mfrow = c(3,3), mai = c(0.5,0.7,0.1,0.1), omi = c(0,0,0,0))
for(i in 1:N_ts) {
  up <- mod_fit_ape4_cov2_de$up[i,]
  mn <- mod_fit_ape4_cov2_de$ex[i,]
  lo <- mod_fit_ape4_cov2_de$lo[i,]
  plot(w_ts,mn,xlab = "",ylab=ylbl[i],xaxt = "n",type = "n", cex.lab = 1.2,
       ylim=c(min(lo),max(up)), bty = "n")
  axis(1,seq(1,42,10),as.character(seq(1982,2023,10)))
  # minor.tick(nx = 3, tick.ratio=0.3)
  points(w_ts,dat[i,], pch=16, col=clr[i])
  lines(w_ts, up, col="darkgray")
  lines(w_ts, mn, col="gray40", lwd = 2)
  lines(w_ts, lo, col="darkgray")
}


# ---------------------------------------------------------------------------- #
# Covariate effects
ape_dfa4_cov2_de_mle <- MARSSparamCIs(ape_dfa4_cov2_de)
ape_dfa4_cov2_de_mle$parMean[31:39]
ape_dfa4_cov2_de_mle$par.upCI$A
ape_dfa4_cov2_de_mle$par.lowCI$A

install.packages("plotrix")      
library("plotrix")
plotCI(1:9, y=ape_dfa4_cov2_de_mle$parMean[31:39],
       li = ape_dfa4_cov2_de_mle$par.lowCI$A,
       ui = ape_dfa4_cov2_de_mle$par.upCI$A,
       ylab = "covariate effect (sea ice)",
       xaxt= 'n', xlab= "",
       scol = clr[1:9],
       main = "Apex predators")
abline(h=0, col="gray50")



# # sum of squared residuals
# ssq_resid_de <- matrix(nrow = length(apex_predators), ncol = 1)
# ssq_resid_eq <- matrix(nrow = length(apex_predators), ncol = 1)
# for(i in 1:(length(apex_predators))) {
#   ssq_resid_de[i,] <- sum((mod_fit_ape_de$ex[apex_predators[i],] - dat[apex_predators[i],])^2, na.rm = TRUE)
#   ssq_resid_eq[i,] <- sum((mod_fit_ape_eq$ex[apex_predators[i],] - dat[apex_predators[i],])^2, na.rm = TRUE)
# }
# # mean of observations
# rowMeans(dat, na.rm = TRUE)
# # sum of squared observations
# ssq_obs_de <- matrix(nrow = length(apex_predators), ncol = 1)
# # ssq_obs_eq <- matrix(nrow = length(apex_predators), ncol = 1)
# for(i in 1:(length(apex_predators))) {
#   ssq_obs_de[i,] <- sum((dat[apex_predators[i],] - 0)^2, na.rm = TRUE)
#   # ssq_obs_eq[i,] <- sum((mod_fit_ape_eq$ex[apex_predators[i],] - 0)^2, na.rm = TRUE)
# }
# 
# # fit ratios
# mean_fit_ratio_de <- mean(ssq_resid_de/ssq_obs_de)
# mean_fit_ratio_eq <- mean(ssq_resid_eq/ssq_obs_de)
# ape_dfa3_cov3_de_fit_ratio <- get_mean_fit_ratio(mod_fit_ape_de, dat)
# ape_dfa3_cov3_eq_fit_ratio <- get_mean_fit_ratio(mod_fit_ape_eq, dat)

