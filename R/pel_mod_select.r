library('MARSS')
library('viridis')
library('Hmisc')
source('helper_funcs.r')
# based on Fish 550 lecture 10: DFA
# pelagic foragers
pelagic_dat <- read.csv("data/pelag_forg_2023.csv")
# They "demean" in their example as opposed to z-score but I'm pretty sure I need
# divide by the sd in addition to demean-ing, hence scale()
pelagic_wide <- scale(pelagic_dat[,2:12])

## ----dfa-trans-data----------------------------------------------------------------------------
## transpose data so time goes across columns
dat <- t(pelagic_wide)
colnames(dat) <- pelagic_dat[,1]
## get number of time series
N_ts <- dim(dat)[1]
## get length of time series
TT <- dim(dat)[2] 
rownames(dat)
pelagic_foragers <- rownames(dat)

## ----dfa-plot-phytos, fig.height=9, fig.width=8, fig.cap='Demeaned time series of Lake Washington phytoplankton.'----
spp <- rownames(dat)
clr <- viridis(dim(dat)[1]+1)
cnt <- 1
par(mfrow = c(3,4), mai = c(0.5,0.7,0.2,0.1), omi = c(0,0,0,0))
for(i in spp){
  plot(1982:2023, dat[i,], xlab = "", ylab="biomass", bty = "L", xaxt = "n", 
       pch=16, col=clr[cnt], type="b", cex.axis=0.75)
  # axis(1,1:dim(dat)[2],1982:2023)
  axis(1, labels = TRUE, cex.axis=0.75)
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

## ----dfa-dfa-proc-eqn--------------------------------------------------------------------------
## number of processes
# mm <- 3

## list with model control parameters
con_list <- list(maxit = 10000, allow.degen = TRUE)#, abstol=0.001)

## ----dfa-fit-dfa-1, cache=TRUE-----------------------------------------------------------------
## fit MARSS
ptm <- proc.time()
# if you use form="dfa" or MARSS.dfa(), then the mod_list cannot contain Z, e.g.,
# determine # of trends and include all covariates and fix R at default
mod_list = list(m = 4, R = "diagonal and equal")
pel_dfa4_cov123_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                            control = con_list, covariates = rbind(npi, ice_ext, totc))
mod_list = list(m = 3, R = "diagonal and equal")
pel_dfa3_cov123_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                          control = con_list, covariates = rbind(npi, ice_ext, totc))
mod_list = list(m = 2, R = "diagonal and equal")
pel_dfa2_cov123_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                          control = con_list, covariates = rbind(npi, ice_ext, totc))
mod_list = list(m = 1, R = "diagonal and equal")
pel_dfa1_cov123_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                          control = con_list, covariates = rbind(npi, ice_ext, totc))
# MARSSparamCIs(pel_dfa3_cov123_de)
# Compare AICc and logLik
print(cbind(model=c(4, 3, 2, 1),
            logLik=round(c(pel_dfa4_cov123_de$logLik, pel_dfa3_cov123_de$logLik, 
                           pel_dfa2_cov123_de$logLik, pel_dfa1_cov123_de$logLik),3),
            AICc=round(c(pel_dfa4_cov123_de$AICc, pel_dfa3_cov123_de$AICc, 
                         pel_dfa2_cov123_de$AICc, pel_dfa1_cov123_de$AICc),3),
            convergence=c(pel_dfa4_cov123_de$convergence, pel_dfa3_cov123_de$convergence, 
                          pel_dfa2_cov123_de$convergence, pel_dfa1_cov123_de$convergence)),
      quote=FALSE)
# get model fits for 1 and 2 trends
# get model fits
pel_dfa1_cov123_de_mod_fits <- get_DFA1_fits(pel_dfa1_cov123_de)
pel_dfa2_cov123_de_mod_fits <- get_DFA_fits(pel_dfa2_cov123_de)
# get mean fit ratios
pel_dfa1_cov123_de_mean_fit <- get_mean_fit_ratio(pel_dfa1_cov123_de_mod_fits, dat)
pel_dfa2_cov123_de_mean_fit <- get_mean_fit_ratio(pel_dfa2_cov123_de_mod_fits, dat)
# Compare AICc, logLik, and mean fit ratio
print(cbind(R_structure=rep("diagonal and equal", 2),
            trends=c(1, 2),
            logLik=round(c(pel_dfa1_cov123_de$logLik, 
                           pel_dfa2_cov123_de$logLik),3),
            AICc=round(c(pel_dfa1_cov123_de$AICc, 
                         pel_dfa2_cov123_de$AICc),3),
            mean_fit_ratio=round(c(pel_dfa1_cov123_de_mean_fit,
                                   pel_dfa2_cov123_de_mean_fit),3)),
      quote=FALSE)


# The most supported model has 2 trends
# now evaluate the most supported form of R
mod_list = list(m = 2, R = "diagonal and unequal")
pel_dfa2_cov123_du <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                          control = con_list, covariates = rbind(npi, ice_ext, totc))
mod_list = list(m = 2, R = "equalvarcov")
pel_dfa2_cov123_eq <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                          control = con_list, covariates = rbind(npi, ice_ext, totc))
mod_list = list(m = 2, R = "unconstrained")
pel_dfa2_cov123_un <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                          control = con_list, covariates = rbind(npi, ice_ext, totc))
# unconstrained did not converge
# Compare AICc and logLik
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            logLik=round(c(pel_dfa2_cov123_de$logLik, pel_dfa2_cov123_du$logLik,
                           pel_dfa2_cov123_eq$logLik, pel_dfa2_cov123_un$logLik),3),
            AICc=round(c(pel_dfa2_cov123_de$AICc, pel_dfa2_cov123_du$AICc,
                         pel_dfa2_cov123_eq$AICc, pel_dfa2_cov123_un$AICc),3),
            convergence=c(pel_dfa2_cov123_de$convergence, pel_dfa2_cov123_du$convergence,
                          pel_dfa2_cov123_eq$convergence, pel_dfa2_cov123_un$convergence)),
      quote=FALSE)

# get model fits
pel_dfa2_cov123_de_mod_fits <- get_DFA_fits(pel_dfa2_cov123_de)
# get mean fit ratios
pel_dfa2_cov123_de_mean_fit <- get_mean_fit_ratio(pel_dfa2_cov123_de_mod_fits, dat)
# get model fits
pel_dfa2_cov123_eq_mod_fits <- get_DFA_fits(pel_dfa2_cov123_eq)
# get mean fit ratios
pel_dfa2_cov123_eq_mean_fit <- get_mean_fit_ratio(pel_dfa2_cov123_eq_mod_fits, dat)
# Compare AICc, logLik, and mean fit ratio
print(cbind(R_structure=c("diagonal and equal", "equalvarcov"),
            trends=c(2, 2),
            logLik=round(c(pel_dfa2_cov123_de$logLik, 
                           pel_dfa2_cov123_eq$logLik),3),
            AICc=round(c(pel_dfa2_cov123_de$AICc, 
                         pel_dfa2_cov123_eq$AICc),3),
            mean_fit_ratio=round(c(pel_dfa2_cov123_de_mean_fit,
                                   pel_dfa2_cov123_eq_mean_fit),3)),
      quote=FALSE)

proc.time() - ptm


# Select covariates ---------- #
# m = 2, R = "diagonal and equal"

ptm <- proc.time()
# all three covariates
mod_list = list(m = 2, R = "diagonal and equal")
pel_dfa2_cov123_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                            control = con_list, covariates = rbind(npi, ice_ext, totc))
pel_dfa2_cov12_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                           control = con_list, covariates = rbind(npi, ice_ext))
pel_dfa2_cov13_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                           control = con_list, covariates = rbind(npi, totc))
pel_dfa2_cov23_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                           control = con_list, covariates = rbind(ice_ext, totc))
pel_dfa2_cov1_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                          control = con_list, covariates = npi)
pel_dfa2_cov2_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                          control = con_list, covariates = ice_ext)
pel_dfa2_cov3_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                          control = con_list, covariates = totc)
pel_dfa2_cov0_de <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                          control = con_list)

# get model fits
pel_dfa2_cov123_de_mod_fits <- get_DFA_fits(pel_dfa2_cov123_de)
pel_dfa2_cov12_de_mod_fits <- get_DFA_fits(pel_dfa2_cov12_de)
pel_dfa2_cov13_de_mod_fits <- get_DFA_fits(pel_dfa2_cov13_de)
pel_dfa2_cov23_de_mod_fits <- get_DFA_fits(pel_dfa2_cov23_de)
pel_dfa2_cov1_de_mod_fits <- get_DFA_fits(pel_dfa2_cov1_de)
pel_dfa2_cov2_de_mod_fits <- get_DFA_fits(pel_dfa2_cov2_de)
pel_dfa2_cov3_de_mod_fits <- get_DFA_fits(pel_dfa2_cov3_de)
pel_dfa2_cov0_de_mod_fits <- get_DFA_fits(pel_dfa2_cov0_de)
# get mean fit ratios
pel_dfa2_cov123_de_mean_fit <- get_mean_fit_ratio(pel_dfa2_cov123_de_mod_fits, dat)
pel_dfa2_cov12_de_mean_fit <- get_mean_fit_ratio(pel_dfa2_cov12_de_mod_fits, dat)
pel_dfa2_cov13_de_mean_fit <- get_mean_fit_ratio(pel_dfa2_cov13_de_mod_fits, dat)
pel_dfa2_cov23_de_mean_fit <- get_mean_fit_ratio(pel_dfa2_cov23_de_mod_fits, dat)
pel_dfa2_cov1_de_mean_fit <- get_mean_fit_ratio(pel_dfa2_cov1_de_mod_fits, dat)
pel_dfa2_cov2_de_mean_fit <- get_mean_fit_ratio(pel_dfa2_cov2_de_mod_fits, dat)
pel_dfa2_cov3_de_mean_fit <- get_mean_fit_ratio(pel_dfa2_cov3_de_mod_fits, dat)
pel_dfa2_cov0_de_mean_fit <- get_mean_fit_ratio(pel_dfa2_cov0_de_mod_fits, dat)

# Compare AICc, logLik, and mean fit ratio
print(cbind(R_structure=rep("diagonal and equal", 8),
            covars=c(123, 12, 13, 23, 1, 2, 3, 0),
            logLik=round(c(pel_dfa2_cov123_de$logLik, 
                           pel_dfa2_cov12_de$logLik, 
                           pel_dfa2_cov13_de$logLik, 
                           pel_dfa2_cov23_de$logLik,
                           pel_dfa2_cov1_de$logLik,
                           pel_dfa2_cov2_de$logLik,
                           pel_dfa2_cov3_de$logLik,
                           pel_dfa2_cov0_de$logLik),3),
            AICc=round(c(pel_dfa2_cov123_de$AICc, 
                         pel_dfa2_cov12_de$AICc, 
                         pel_dfa2_cov13_de$AICc, 
                         pel_dfa2_cov23_de$AICc,
                         pel_dfa2_cov1_de$AICc,
                         pel_dfa2_cov2_de$AICc,
                         pel_dfa2_cov3_de$AICc,
                         pel_dfa2_cov0_de$AICc),3),
            mean_fit_ratio=round(c(pel_dfa2_cov123_de_mean_fit,
                                   pel_dfa2_cov12_de_mean_fit,
                                   pel_dfa2_cov13_de_mean_fit,
                                   pel_dfa2_cov23_de_mean_fit,
                                   pel_dfa2_cov1_de_mean_fit,
                                   pel_dfa2_cov2_de_mean_fit,
                                   pel_dfa2_cov3_de_mean_fit,
                                   pel_dfa2_cov0_de_mean_fit),3),
            convergence=c(pel_dfa2_cov123_de$convergence, 
                          pel_dfa2_cov12_de$convergence, 
                          pel_dfa2_cov13_de$convergence, 
                          pel_dfa2_cov23_de$convergence,
                          pel_dfa2_cov1_de$convergence,
                          pel_dfa2_cov2_de$convergence,
                          pel_dfa2_cov3_de$convergence,
                          pel_dfa2_cov0_de$convergence)),
      quote=FALSE)

proc.time() - ptm


# ============================================================================ #
# Best model plots

## number of processes
mm <- 2
## get the estimated ZZ
Z_est <- coef(pel_dfa2_cov2_de, type = "matrix")$Z
## get the inverse of the rotation matrix
H_inv <- varimax(Z_est)$rotmat
## ----dfa-rotate-Z-x----------------------------------------------------------------------------
## rotate factor loadings
Z_rot = Z_est %*% H_inv   
## rotate processes
proc_rot = solve(H_inv) %*% pel_dfa2_cov2_de$states

## ----dfa-plot-dfa1, fig.height=9, fig.width=8, eval=TRUE, fig.cap='Estimated states from the DFA model.'----
ylbl <- pelagic_foragers
w_ts <- seq(dim(dat)[2])
# the first four trends
layout(matrix(c(1,2,3,4), mm, 2), widths = c(2,1))
## par(mfcol=c(mm,2), mai = c(0.5,0.5,0.5,0.1), omi = c(0,0,0,0))
par(mai = c(0.5,0.5,0.5,0.1), omi = c(0,0,0,0))
## plot the processes
for(i in 1:mm) {
  ylm <- c(-1,1)*max(abs(proc_rot[i,]))
  ## set up plot area
  plot(w_ts,proc_rot[i,], type = "n", bty = "L",
       ylim = ylm, xlab = "", ylab = "", xaxt = "n", cex=1, cex.axis=0.75)
  ## draw zero-line
  abline(h=0, col="gray")
  ## plot trend line
  lines(w_ts, proc_rot[i,], lwd = 2)
  lines(w_ts, proc_rot[i,], lwd = 2)
  ## add panel labels
  mtext(paste("State",i,"_pelagic_foragers", " diagonal and equal"), side = 3, 
        line = 0.5, cex=0.75)
  # axis(1,12*(0:dim(dat_1980)[2])+1,yr_frst+0:dim(dat_1980)[2])
  axis(1,seq(1,42,10),as.character(seq(1982,2023,10)), cex.axis=0.75)
}
## plot the loadings
minZ <- 0
ylm <- c(-1,1)*max(abs(Z_rot))
for(i in 1:mm) {
  plot(c(1:N_ts)[abs(Z_rot[,i])>minZ], as.vector(Z_rot[abs(Z_rot[,i])>minZ,i]), 
       type="h", lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylm, 
       xlim=c(0.5,N_ts+0.5), 
       col=clr, cex=1, bty="n", cex.axis=0.75)
  abline(h=c(-0.2,0.2), lty=2, col="gray")
  for(j in 1:N_ts) {
    if(Z_rot[j,i] > minZ) {text(j, -0.03, ylbl[j], srt=90, adj=1, cex=0.75, col=clr[j])}
    if(Z_rot[j,i] < -minZ) {text(j, 0.03, ylbl[j], srt=90, adj=0, cex=0.75, col=clr[j])}
    abline(h=0, lwd=1.5, col="gray")
  } 
  mtext(paste("Factor loadings"),side=3,line=0.5, cex=0.75)
}

## ----dfa-xy-states12, height=4, width=5, fig.cap='Cross-correlation plot of the two rotations.'----
par(mai = c(0.9,0.9,0.1,0.1))
ccf(proc_rot[1,],proc_rot[2,], lag.max = 12, main="")
acf(proc_rot[1,], type = "covariance")


## ----dfa-plot-dfa-fits, fig.height=9, fig.width=8, fig.cap='Data and fits from the DFA model.'----
## get model fits & CI's
# pel_dfa2_cov2_de_mod_fits <- get_DFA_fits(pel_dfa2_cov2_de)
## plot the fits
ylbl <- pelagic_foragers
par(mfrow = c(3,4), mai = c(0.5,0.7,0.1,0.1), omi = c(0,0,0,0))
for(i in 1:N_ts) {
  up <- pel_dfa2_cov2_de_mod_fits$up[i,]
  mn <- pel_dfa2_cov2_de_mod_fits$ex[i,]
  lo <- pel_dfa2_cov2_de_mod_fits$lo[i,]
  plot(w_ts,mn,xlab = "",ylab=ylbl[i],xaxt = "n",type = "n", cex.lab = 1.2,
       ylim=c(min(lo),max(up)), bty="n", cex.axis=0.75)
  axis(1,seq(1,42,10),as.character(seq(1982,2023,10)), cex.axis=0.75)
  points(w_ts,dat[i,], pch=16, col=clr[i])
  lines(w_ts, up, col="darkgray")
  lines(w_ts, mn, col="gray40", lwd = 2)
  lines(w_ts, lo, col="darkgray")
}






#==============================================================================#
#==============================================================================#
# comparing plots for 1 and 2 trends (with all 3 covariates and R="diagonal and equal")
## number of processes
mm <- 2

## ----dfa-plot-dfa-fits, fig.height=9, fig.width=8, fig.cap='Data and fits from the DFA model.'----
## get model fits & CI's
pel_dfa2_cov123_de_mod_fits <- get_DFA_fits(pel_dfa2_cov123_de)
## plot the fits
ylbl <- pelagic_foragers
par(mfrow = c(3,4), mai = c(0.5,0.7,0.1,0.1), omi = c(0,0,0,0))
for(i in 1:N_ts) {
  up <- pel_dfa2_cov123_de_mod_fits$up[i,]
  mn <- pel_dfa2_cov123_de_mod_fits$ex[i,]
  lo <- pel_dfa2_cov123_de_mod_fits$lo[i,]
  plot(w_ts,mn,xlab = "",ylab=ylbl[i],xaxt = "n",type = "n", cex.lab = 1.2,
       ylim=c(min(lo),max(up)))
  axis(1,seq(1,42,10),as.character(seq(1982,2023,10)))
  points(w_ts,dat[i,], pch=16, col=clr[i])
  lines(w_ts, up, col="darkgray")
  lines(w_ts, mn, col="black", lwd = 2)
  lines(w_ts, lo, col="darkgray")
}

mm <- 1

## ----dfa-plot-dfa-fits, fig.height=9, fig.width=8, fig.cap='Data and fits from the DFA model.'----
## get model fits & CI's
pel_dfa1_cov123_de_mod_fits <- get_DFA1_fits(pel_dfa1_cov123_de)
## plot the fits
ylbl <- pelagic_foragers
par(mfrow = c(3,4), mai = c(0.5,0.7,0.1,0.1), omi = c(0,0,0,0))
for(i in 1:N_ts) {
  up <- pel_dfa1_cov123_de_mod_fits$up[i,]
  mn <- pel_dfa1_cov123_de_mod_fits$ex[i,]
  lo <- pel_dfa1_cov123_de_mod_fits$lo[i,]
  plot(w_ts,mn,xlab = "",ylab=ylbl[i],xaxt = "n",type = "n", cex.lab = 1.2,
       ylim=c(min(lo),max(up)))
  axis(1,seq(1,42,10),as.character(seq(1982,2023,10)))
  points(w_ts,dat[i,], pch=16, col=clr[i])
  lines(w_ts, up, col="darkgray")
  lines(w_ts, mn, col="black", lwd = 2)
  lines(w_ts, lo, col="darkgray")
}



