# DFA with MARSS
# install.packages("MARSS")
library('MARSS')
library('viridis')
# Fish 550 lecture 10: DFA
# Modified here for use with feeding guilds
setwd("C:/Users/andy.whitehouse/Work/Andy/REEM/Eco_Considerations/time_series_analysis")

guild_dat <- read.csv("ebs_2023_guilds.csv")
# pull out the apex predators
apex <- guild_dat[guild_dat$guild == "Apex predators", ]
apex_names <- unique(apex$ebs_name)

# convert from long to wide
# library('tidyr')
# 
# apex %>% 
#   pivot_wider(names_from = ebs_name, values_from = YEAR)
# 
# apex_wide <- reshape(apex, idvar = "ebs_name", timevar = "YEAR", v.names = "tot_bio_tons",direction = "wide")
# reshape(df, idvar="year", timevar="month", v.names="values", direction="wide", sep="_")

# or just read in a csv
apex_dat <- read.csv("apex_2023.csv")
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
clr <- viridis(dim(dat)[1])
cnt <- 1
par(mfrow = c(3,3), mai = c(0.5,0.7,0.1,0.1), omi = c(0,0,0,0))
for(i in spp){
  plot(dat[i,],xlab = "",ylab="biomass", bty = "L", xaxt = "n", pch=16, col=clr[cnt], type="b")
  axis(1,1:dim(dat)[2],1982:2023)
  title(i)
  cnt <- cnt + 1
}


## ----dfa-dfa-obs-eqn---------------------------------------------------------------------------
## 'ZZ' is loadings matrix
# 3 trends
# 9 time series
Z_vals <- list("z11",  0  ,  0  ,
               "z21","z22",  0  ,
               "z31","z32","z33",
               "z41","z42","z43",
               "z51","z52","z53",
               "z61","z62","z63",
               "z71","z72","z73",
               "z81","z82","z83",
               "z91","z92","z93")
ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 3, byrow = TRUE)
ZZ
## 'aa' is the offset/scaling
aa <- "zero"
## 'DD' and 'd' are for covariates
DD <- "zero"  # matrix(0,mm,1)
dd <- "zero"  # matrix(0,1,wk_last)
## 'RR' is var-cov matrix for obs errors
RR <- "diagonal and equal"


## ----dfa-dfa-proc-eqn--------------------------------------------------------------------------
## number of processes
mm <- 3
## 'BB' is identity: 1's along the diagonal & 0's elsewhere
BB <- "identity"  # diag(mm)
## 'uu' is a column vector of 0's
uu <- "zero"  # matrix(0, mm, 1)
## 'CC' and 'cc' are for covariates
CC <- "zero"  # matrix(0, mm, 1)
cc <- "zero"  # matrix(0, 1, wk_last)
## 'QQ' is identity
QQ <- "identity"  # diag(mm)


## ----dfa-create-model-lists--------------------------------------------------------------------
## list with specifications for model vectors/matrices
mod_list <- list(Z = ZZ, A = aa, D = DD, d = dd, R = RR,
                 B = BB, U = uu, C = CC, c = cc, Q = QQ)
## list with model inits
init_list <- list(x0 = matrix(rep(0, mm), mm, 1))
## list with model control parameters
con_list <- list(maxit = 4000, allow.degen = TRUE)


## ----dfa-fit-dfa-1, cache=TRUE-----------------------------------------------------------------
## fit MARSS
ptm <- proc.time()
dfa_1 <- MARSS(y = dat, model = mod_list, inits = init_list, control = con_list)
proc.time() - ptm
names(dfa_1)
# the trends are in $states and the SE's are in $states.se

## ----dfa-get-H-inv-----------------------------------------------------------------------------
## get the estimated ZZ
Z_est <- coef(dfa_1, type = "matrix")$Z
## get the inverse of the rotation matrix
H_inv <- varimax(Z_est)$rotmat


## ----dfa-rotate-Z-x----------------------------------------------------------------------------
## rotate factor loadings
Z_rot = Z_est %*% H_inv   
## rotate processes
proc_rot = solve(H_inv) %*% dfa_1$states


## ----dfa-plot-dfa1, fig.height=9, fig.width=8, eval=TRUE, fig.cap='Estimated states from the DFA model.'----
ylbl <- apex_predators
w_ts <- seq(dim(dat)[2])
layout(matrix(c(1,2,3,4,5,6), mm, 2), widths = c(2,1))
## par(mfcol=c(mm,2), mai = c(0.5,0.5,0.5,0.1), omi = c(0,0,0,0))
par(mai = c(0.5,0.5,0.5,0.1), omi = c(0,0,0,0))
## plot the processes
for(i in 1:mm) {
  ylm <- c(-1,1)*max(abs(proc_rot[i,]))
  ## set up plot area
  plot(w_ts,proc_rot[i,], type = "n", bty = "L",
       ylim = ylm, xlab = "", ylab = "", xaxt = "n")
  ## draw zero-line
  abline(h=0, col="gray")
  ## plot trend line
  lines(w_ts, proc_rot[i,], lwd = 2)
  lines(w_ts, proc_rot[i,], lwd = 2)
  ## add panel labels
  mtext(paste("State",i,"_apex"), side = 3, line = 0.5)
  # axis(1,12*(0:dim(dat_1980)[2])+1,yr_frst+0:dim(dat_1980)[2])
  axis(1,seq(1,42,10),as.character(seq(1982,2023,10)))
}
## plot the loadings
minZ <- 0
ylm <- c(-1,1)*max(abs(Z_rot))
for(i in 1:mm) {
  plot(c(1:N_ts)[abs(Z_rot[,i])>minZ], as.vector(Z_rot[abs(Z_rot[,i])>minZ,i]), type="h",
       lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylm, xlim=c(0.5,N_ts+0.5), 
       col=clr, cex=1.5)
  for(j in 1:N_ts) {
    if(Z_rot[j,i] > minZ) {text(j, -0.03, ylbl[j], srt=90, adj=1, cex=1.2, col=clr[j])}
    if(Z_rot[j,i] < -minZ) {text(j, 0.03, ylbl[j], srt=90, adj=0, cex=1.2, col=clr[j])}
    abline(h=0, lwd=1.5, col="gray")
  } 
  mtext(paste("Factor loadings on state",i),side=3,line=0.5)
}


## ----dfa-xy-states12, height=4, width=5, fig.cap='Cross-correlation plot of the two rotations.'----
par(mai = c(0.9,0.9,0.1,0.1))
ccf(proc_rot[1,],proc_rot[2,], lag.max = 12, main="")


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


## ----dfa-plot-dfa-fits, fig.height=9, fig.width=8, fig.cap='Data and fits from the DFA model.'----
## get model fits & CI's
mod_fit <- get_DFA_fits(dfa_1)
## plot the fits
ylbl <- apex_predators
par(mfrow = c(3,3), mai = c(0.5,0.7,0.1,0.1), omi = c(0,0,0,0))
for(i in 1:N_ts) {
  up <- mod_fit$up[i,]
  mn <- mod_fit$ex[i,]
  lo <- mod_fit$lo[i,]
  plot(w_ts,mn,xlab = "",ylab=ylbl[i],xaxt = "n",type = "n", cex.lab = 1.2,
       ylim=c(min(lo),max(up)))
  axis(1,seq(1,42,10),as.character(seq(1982,2023,10)))
  points(w_ts,dat[i,], pch=16, col=clr[i])
  lines(w_ts, up, col="darkgray")
  lines(w_ts, mn, col="black", lwd = 2)
  lines(w_ts, lo, col="darkgray")
}

# DFA with covariate

# coldpool data
coldpool_data <- read.csv("coldpool.csv")
cpz <- scale(coldpool_data[,2:11])
cp_area <- t(cpz[,"AREA_LTE2_KM2"])
npi <- t(cpz[,"npi"])
ice_ext <- t(cpz[,"ice_ext"])
pdo <- t(cpz[,"pdo"])
mgt <- t(cpz[,"MEAN_GEAR_TEMPERATURE"])
mst <- t(cpz[,"MEAN_SURFACE_TEMPERATURE"])
# cp_area <- t(coldpool_data$AREA_LTE2_KM2)
# npi <- t(coldpool_data$npi)
# ice_ext <- t(coldpool_data$ice_ext)
# pdo <- t(coldpool_data$pdo)

## ----dfa-fit-DFA-covars, cache=TRUE, results='hide'--------------------------------------------
mod_list = list(m = 3, R = "diagonal and equal")
dfa_apex_cpa <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                      control = con_list, covariates=cp_area)
dfa_apex_npi <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                      control = con_list, covariates=npi)
dfa_apex_ice <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                      control = con_list, covariates=ice_ext)
dfa_apex_pdo <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                      control = con_list, covariates=pdo)
## ----dfa-model-selection-----------------------------------------------------------------------
print(cbind(model=c("no covars", "cpa", "npi", "ice", "pdo"),
            AICc=round(c(dfa_1$AICc, dfa_apex_cpa$AICc, dfa_apex_npi$AICc, 
                         dfa_apex_ice$AICc, dfa_apex_pdo$AICc))),
      quote=FALSE)

print(cbind(model=c("no covars", "cpa", "npi", "ice", "pdo"),
            logLik=round(c(dfa_1$logLik, dfa_apex_cpa$logLik, dfa_apex_npi$logLik, 
                           dfa_apex_ice$logLik, dfa_apex_pdo$logLik))),
      quote=FALSE)

# ============================================================================ #
# pelagic foragers
# or just read in a csv
pelagic_dat <- read.csv("pelag_forg_2023.csv")
# They "demean" in their example as opposed to z-score but I'm pretty sure I need
# divide by the sd in addition to demean-ing, hence scale()
pel_wide <- scale(pelagic_dat[,2:dim(pelagic_dat)[2]])
## ----dfa-trans-data----------------------------------------------------------------------------
## transpose data so time goes across columns
pel_dat <- t(pel_wide)
colnames(pel_dat) <- pelagic_dat[,1]
## get number of time series
N_ts <- dim(pel_dat)[1]
## get length of time series
TT <- dim(pel_dat)[2] 
rownames(pel_dat)
pelagic_foragers <- rownames(pel_dat)

## ----dfa-plot-phytos, fig.height=9, fig.width=8, fig.cap='Demeaned time series of Lake Washington phytoplankton.'----
spp <- rownames(pel_dat)
clr <- viridis(length(pelagic_foragers))
cnt <- 1
par(mfrow = c(4,3), mai = c(0.5,0.7,0.1,0.1), omi = c(0,0,0,0))
for(i in spp){
  plot(pel_dat[i,],xlab = "",ylab="biomass", bty = "L", xaxt = "n", pch=16, col=clr[cnt], type="b")
  axis(1,1:dim(pel_dat)[2],1982:2023)
  title(i)
  cnt <- cnt + 1
}


## ----dfa-dfa-obs-eqn---------------------------------------------------------------------------
## 'ZZ' is loadings matrix
# 3 trends
# 11 time series
Z_vals <- list("z11",  0  ,  0  ,
               "z21","z22",  0  ,
               "z31","z32","z33",
               "z41","z42","z43",
               "z51","z52","z53",
               "z61","z62","z63",
               "z71","z72","z73",
               "z81","z82","z83",
               "z91","z92","z93",
               "z101","z102","z103",
               "z111","z112","z113")
ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 3, byrow = TRUE)
ZZ
## 'aa' is the offset/scaling
aa <- "zero"
## 'DD' and 'd' are for covariates
DD <- "zero"  # matrix(0,mm,1)
dd <- "zero"  # matrix(0,1,wk_last)
## 'RR' is var-cov matrix for obs errors
RR <- "diagonal and equal"


## ----dfa-dfa-proc-eqn--------------------------------------------------------------------------
## number of processes
mm <- 3
## 'BB' is identity: 1's along the diagonal & 0's elsewhere
BB <- "identity"  # diag(mm)
## 'uu' is a column vector of 0's
uu <- "zero"  # matrix(0, mm, 1)
## 'CC' and 'cc' are for covariates
CC <- "zero"  # matrix(0, mm, 1)
cc <- "zero"  # matrix(0, 1, wk_last)
## 'QQ' is identity
QQ <- "identity"  # diag(mm)


## ----dfa-create-model-lists--------------------------------------------------------------------
## list with specifications for model vectors/matrices
mod_list <- list(Z = ZZ, A = aa, D = DD, d = dd, R = RR,
                 B = BB, U = uu, C = CC, c = cc, Q = QQ)
## list with model inits
init_list <- list(x0 = matrix(rep(0, mm), mm, 1))
## list with model control parameters
con_list <- list(maxit = 4000, allow.degen = TRUE)


## ----dfa-fit-dfa-1, cache=TRUE-----------------------------------------------------------------
## fit MARSS
ptm <- proc.time()
dfa_pel <- MARSS(y = pel_dat, model = mod_list, inits = init_list, control = con_list)
proc.time() - ptm
names(dfa_pel)
# the trends are in $states and the SE's are in $states.se

## ----dfa-get-H-inv-----------------------------------------------------------------------------
## get the estimated ZZ
Z_est <- coef(dfa_pel, type = "matrix")$Z
## get the inverse of the rotation matrix
H_inv <- varimax(Z_est)$rotmat


## ----dfa-rotate-Z-x----------------------------------------------------------------------------
## rotate factor loadings
Z_rot = Z_est %*% H_inv   
## rotate processes
proc_rot = solve(H_inv) %*% dfa_pel$states


## ----dfa-plot-dfa1, fig.height=9, fig.width=8, eval=TRUE, fig.cap='Estimated states from the DFA model.'----
ylbl <- pelagic_foragers
w_ts <- seq(dim(pel_dat)[2])
layout(matrix(c(1,2,3,4,5,6), mm, 2), widths = c(2,1))
## par(mfcol=c(mm,2), mai = c(0.5,0.5,0.5,0.1), omi = c(0,0,0,0))
par(mai = c(0.5,0.5,0.5,0.1), omi = c(0,0,0,0))
## plot the processes
for(i in 1:mm) {
  ylm <- c(-1,1)*max(abs(proc_rot[i,]))
  ## set up plot area
  plot(w_ts,proc_rot[i,], type = "n", bty = "L",
       ylim = ylm, xlab = "", ylab = "", xaxt = "n")
  ## draw zero-line
  abline(h=0, col="gray")
  ## plot trend line
  lines(w_ts, proc_rot[i,], lwd = 2)
  lines(w_ts, proc_rot[i,], lwd = 2)
  ## add panel labels
  mtext(paste("State",i,"_pel_forg"), side = 3, line = 0.5)
  # axis(1,12*(0:dim(dat_1980)[2])+1,yr_frst+0:dim(dat_1980)[2])
  axis(1,seq(1,42,10),as.character(seq(1982,2023,10)))
}
## plot the loadings
minZ <- 0
ylm <- c(-1,1)*max(abs(Z_rot))
for(i in 1:mm) {
  plot(c(1:N_ts)[abs(Z_rot[,i])>minZ], as.vector(Z_rot[abs(Z_rot[,i])>minZ,i]), type="h",
       lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylm, xlim=c(0.5,N_ts+0.5), 
       col=clr, cex=1.5)
  for(j in 1:N_ts) {
    if(Z_rot[j,i] > minZ) {text(j, -0.03, ylbl[j], srt=90, adj=1, cex=1.2, col=clr[j])}
    if(Z_rot[j,i] < -minZ) {text(j, 0.03, ylbl[j], srt=90, adj=0, cex=1.2, col=clr[j])}
    abline(h=0, lwd=1.5, col="gray")
  } 
  mtext(paste("Factor loadings on state",i),side=3,line=0.5)
}


## ----dfa-xy-states12, height=4, width=5, fig.cap='Cross-correlation plot of the two rotations.'----
par(mai = c(0.9,0.9,0.1,0.1))
ccf(proc_rot[1,],proc_rot[2,], lag.max = 12, main="")


## ----dfa-plot-dfa-fits, fig.height=9, fig.width=8, fig.cap='Data and fits from the DFA model.'----
## get model fits & CI's
mod_fit <- get_DFA_fits(dfa_pel)
## plot the fits
ylbl <- pelagic_foragers
par(mfrow = c(4,3), mai = c(0.5,0.7,0.1,0.1), omi = c(0,0,0,0))
for(i in 1:N_ts) {
  up <- mod_fit$up[i,]
  mn <- mod_fit$ex[i,]
  lo <- mod_fit$lo[i,]
  plot(w_ts,mn,xlab = "",ylab=ylbl[i],xaxt = "n",type = "n", cex.lab = 1.2,
       ylim=c(min(lo),max(up)))
  axis(1,seq(1,42,10),as.character(seq(1982,2023,10)))
  points(w_ts,pel_dat[i,], pch=16, col=clr[i])
  lines(w_ts, up, col="darkgray")
  lines(w_ts, mn, col="black", lwd = 2)
  lines(w_ts, lo, col="darkgray")
}

## ----dfa-fit-DFA-covars, cache=TRUE, results='hide'--------------------------------------------
mod_list = list(m = 3, R = "diagonal and equal")
dfa_pel_cpa <- MARSS(pel_dat, model = mod_list, form = "dfa", z.score = FALSE,
                     control = con_list, covariates=cp_area)
dfa_pel_npi <- MARSS(pel_dat, model = mod_list, form = "dfa", z.score = FALSE,
                      control = con_list, covariates=npi)
dfa_pel_ice <- MARSS(pel_dat, model = mod_list, form = "dfa", z.score = FALSE,
                      control = con_list, covariates=ice_ext)
dfa_pel_pdo <- MARSS(pel_dat, model = mod_list, form = "dfa", z.score = FALSE,
                      control = con_list, covariates=pdo)
dfa_pel_mgt <- MARSS(pel_dat, model = mod_list, form = "dfa", z.score = FALSE,
                     control = con_list, covariates=mgt)
dfa_pel_mst <- MARSS(pel_dat, model = mod_list, form = "dfa", z.score = FALSE,
                     control = con_list, covariates=mst)
dfa_pel_ice_mst <- MARSS(pel_dat, model = mod_list, form = "dfa", z.score = FALSE,
                     control = con_list, covariates=rbind(ice_ext,mst))
## ----dfa-get-H-inv-----------------------------------------------------------------------------
## get the estimated ZZ
Z_est <- coef(dfa_pel_cpa, type = "matrix")$Z
## get the inverse of the rotation matrix
H_inv <- varimax(Z_est)$rotmat
## ----dfa-rotate-Z-x----------------------------------------------------------------------------
## rotate factor loadings
Z_rot = Z_est %*% H_inv   
## rotate processes
proc_rot = solve(H_inv) %*% dfa_pel_cpa$states
## ----dfa-xy-states12, height=4, width=5, fig.cap='Cross-correlation plot of the two rotations.'----
par(mai = c(0.9,0.9,0.1,0.1))
ccf(proc_rot[1,],proc_rot[2,], lag.max = 12, main="")


## ----dfa-model-selection-----------------------------------------------------------------------
print(cbind(model=c("no covars", "cpa", "npi", "ice", "pdo", "mgt", "mst", "ice_mst"),
            AICc=round(c(dfa_1$AICc, dfa_pel_cpa$AICc, dfa_pel_npi$AICc, dfa_pel_ice$AICc, 
                         dfa_pel_pdo$AICc, dfa_pel_mgt$AICc, dfa_pel_mst$AICc,
                         dfa_pel_ice_mst$AICc))),
      quote=FALSE)

print(cbind(model=c("no covars", "cpa", "npi", "ice", "pdo", "mgt", "mst", "ice_mst"),
            logLik=round(c(dfa_1$logLik, dfa_pel_cpa$logLik, dfa_pel_npi$logLik, 
                           dfa_pel_ice$logLik, dfa_pel_pdo$logLik, 
                           dfa_pel_mgt$logLik, dfa_pel_mst$logLik,
                           dfa_pel_ice_mst$logLik))),
      quote=FALSE)


# ============================================================================ #
# benthic foragers
# or just read in a csv
benthic_dat <- read.csv("ben_forg_2023.csv")
# They "demean" in their example as opposed to z-score but I'm pretty sure I need
# divide by the sd in addition to demean-ing, hence scale()
ben_wide <- scale(benthic_dat[,2:dim(benthic_dat)[2]])
## ----dfa-trans-data----------------------------------------------------------------------------
## transpose data so time goes across columns
ben_dat <- t(ben_wide)
colnames(ben_dat) <- benthic_dat[,1]
## get number of time series
N_ts <- dim(ben_dat)[1]
## get length of time series
TT <- dim(ben_dat)[2] 
rownames(ben_dat)
benthic_foragers <- rownames(ben_dat)

## ----dfa-plot-phytos, fig.height=9, fig.width=8, fig.cap='Demeaned time series of Lake Washington phytoplankton.'----
spp <- rownames(ben_dat)
clr <- viridis(length(benthic_foragers))
cnt <- 1
par(mfrow = c(3,3), mai = c(0.5,0.7,0.1,0.1), omi = c(0,0,0,0))
for(i in spp){
  plot(ben_dat[i,],xlab = "",ylab="biomass", bty = "L", xaxt = "n", pch=16, col=clr[cnt], type="b")
  axis(1,1:dim(ben_dat)[2],1982:2023)
  title(i)
  cnt <- cnt + 1
}


## ----dfa-dfa-obs-eqn---------------------------------------------------------------------------
## 'ZZ' is loadings matrix
# 3 trends
# 11 time series
Z_vals <- list("z11",  0  ,  0  ,
               "z21","z22",  0  ,
               "z31","z32","z33",
               "z41","z42","z43",
               "z51","z52","z53",
               "z61","z62","z63",
               "z71","z72","z73",
               "z81","z82","z83",
               "z91","z92","z93")
ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 3, byrow = TRUE)
ZZ
## 'aa' is the offset/scaling
aa <- "zero"
## 'DD' and 'd' are for covariates
DD <- "zero"  # matrix(0,mm,1)
dd <- "zero"  # matrix(0,1,wk_last)
## 'RR' is var-cov matrix for obs errors
RR <- "diagonal and equal"


## ----dfa-dfa-proc-eqn--------------------------------------------------------------------------
## number of processes
mm <- 3
## 'BB' is identity: 1's along the diagonal & 0's elsewhere
BB <- "identity"  # diag(mm)
## 'uu' is a column vector of 0's
uu <- "zero"  # matrix(0, mm, 1)
## 'CC' and 'cc' are for covariates
CC <- "zero"  # matrix(0, mm, 1)
cc <- "zero"  # matrix(0, 1, wk_last)
## 'QQ' is identity
QQ <- "identity"  # diag(mm)


## ----dfa-create-model-lists--------------------------------------------------------------------
## list with specifications for model vectors/matrices
mod_list <- list(Z = ZZ, A = aa, D = DD, d = dd, R = RR,
                 B = BB, U = uu, C = CC, c = cc, Q = QQ)
## list with model inits
init_list <- list(x0 = matrix(rep(0, mm), mm, 1))
## list with model control parameters
con_list <- list(maxit = 4000, allow.degen = TRUE)


## ----dfa-fit-dfa-1, cache=TRUE-----------------------------------------------------------------
## fit MARSS
ptm <- proc.time()
dfa_ben <- MARSS(y = ben_dat, model = mod_list, inits = init_list, control = con_list)
proc.time() - ptm
names(dfa_ben)
# the trends are in $states and the SE's are in $states.se

## ----dfa-get-H-inv-----------------------------------------------------------------------------
## get the estimated ZZ
Z_est <- coef(dfa_ben, type = "matrix")$Z
## get the inverse of the rotation matrix
H_inv <- varimax(Z_est)$rotmat


## ----dfa-rotate-Z-x----------------------------------------------------------------------------
## rotate factor loadings
Z_rot = Z_est %*% H_inv   
## rotate processes
proc_rot = solve(H_inv) %*% dfa_ben$states


## ----dfa-plot-dfa1, fig.height=9, fig.width=8, eval=TRUE, fig.cap='Estimated states from the DFA model.'----
ylbl <- benthic_foragers
w_ts <- seq(dim(ben_dat)[2])
layout(matrix(c(1,2,3,4,5,6), mm, 2), widths = c(2,1))
## par(mfcol=c(mm,2), mai = c(0.5,0.5,0.5,0.1), omi = c(0,0,0,0))
par(mai = c(0.5,0.5,0.5,0.1), omi = c(0,0,0,0))
## plot the processes
for(i in 1:mm) {
  ylm <- c(-1,1)*max(abs(proc_rot[i,]))
  ## set up plot area
  plot(w_ts,proc_rot[i,], type = "n", bty = "L",
       ylim = ylm, xlab = "", ylab = "", xaxt = "n")
  ## draw zero-line
  abline(h=0, col="gray")
  ## plot trend line
  lines(w_ts, proc_rot[i,], lwd = 2)
  lines(w_ts, proc_rot[i,], lwd = 2)
  ## add panel labels
  mtext(paste("State",i,"_ben_forg"), side = 3, line = 0.5)
  # axis(1,12*(0:dim(dat_1980)[2])+1,yr_frst+0:dim(dat_1980)[2])
  axis(1,seq(1,42,10),as.character(seq(1982,2023,10)))
}
## plot the loadings
minZ <- 0
ylm <- c(-1,1)*max(abs(Z_rot))
for(i in 1:mm) {
  plot(c(1:N_ts)[abs(Z_rot[,i])>minZ], as.vector(Z_rot[abs(Z_rot[,i])>minZ,i]), type="h",
       lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylm, xlim=c(0.5,N_ts+0.5), 
       col=clr, cex=1.5)
  for(j in 1:N_ts) {
    if(Z_rot[j,i] > minZ) {text(j, -0.03, ylbl[j], srt=90, adj=1, cex=1.2, col=clr[j])}
    if(Z_rot[j,i] < -minZ) {text(j, 0.03, ylbl[j], srt=90, adj=0, cex=1.2, col=clr[j])}
    abline(h=0, lwd=1.5, col="gray")
  } 
  mtext(paste("Factor loadings on state",i),side=3,line=0.5)
}


## ----dfa-xy-states12, height=4, width=5, fig.cap='Cross-correlation plot of the two rotations.'----
par(mai = c(0.9,0.9,0.1,0.1))
ccf(proc_rot[1,],proc_rot[2,], lag.max = 12, main="")


## ----dfa-plot-dfa-fits, fig.height=9, fig.width=8, fig.cap='Data and fits from the DFA model.'----
## get model fits & CI's
mod_fit <- get_DFA_fits(dfa_ben)
## plot the fits
ylbl <- benthic_foragers
par(mfrow = c(3,3), mai = c(0.5,0.7,0.1,0.1), omi = c(0,0,0,0))
for(i in 1:N_ts) {
  up <- mod_fit$up[i,]
  mn <- mod_fit$ex[i,]
  lo <- mod_fit$lo[i,]
  plot(w_ts,mn,xlab = "",ylab=ylbl[i],xaxt = "n",type = "n", cex.lab = 1.2,
       ylim=c(min(lo),max(up)))
  axis(1,seq(1,42,10),as.character(seq(1982,2023,10)))
  points(w_ts,ben_dat[i,], pch=16, col=clr[i])
  lines(w_ts, up, col="darkgray")
  lines(w_ts, mn, col="black", lwd = 2)
  lines(w_ts, lo, col="darkgray")
}


# ============================================================================ #
# motile_epifauna
# or just read in a csv
motile_dat <- read.csv("mot_epi_2023.csv")
# They "demean" in their example as opposed to z-score but I'm pretty sure I need
# divide by the sd in addition to demean-ing, hence scale()
mote_wide <- scale(motile_dat[,2:dim(motile_dat)[2]])
## ----dfa-trans-data----------------------------------------------------------------------------
## transpose data so time goes across columns
mote_dat <- t(mote_wide)
colnames(mote_dat) <- motile_dat[,1]
## get number of time series
N_ts <- dim(mote_dat)[1]
## get length of time series
TT <- dim(mote_dat)[2] 
rownames(mote_dat)
motile_epifauna <- rownames(mote_dat)

## ----dfa-plot-phytos, fig.height=9, fig.width=8, fig.cap='Demeaned time series of Lake Washington phytoplankton.'----
spp <- rownames(mote_dat)
clr <- viridis(length(motile_epifauna))
cnt <- 1
par(mfrow = c(4,3), mai = c(0.5,0.7,0.1,0.1), omi = c(0,0,0,0))
for(i in spp){
  plot(mote_dat[i,],xlab = "",ylab="biomass", bty = "L", xaxt = "n", pch=16, col=clr[cnt], type="b")
  axis(1,1:dim(mote_dat)[2],1982:2023)
  title(i)
  cnt <- cnt + 1
}


## ----dfa-dfa-obs-eqn---------------------------------------------------------------------------
## 'ZZ' is loadings matrix
# 3 trends
# 11 time series
Z_vals <- list("z11",  0  ,  0  ,
               "z21","z22",  0  ,
               "z31","z32","z33",
               "z41","z42","z43",
               "z51","z52","z53",
               "z61","z62","z63",
               "z71","z72","z73",
               "z81","z82","z83",
               "z91","z92","z93",
               "z101","z102","z103",
               "z111","z112","z113",
               "z121","z122","z123")
ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 3, byrow = TRUE)
ZZ
## 'aa' is the offset/scaling
aa <- "zero"
## 'DD' and 'd' are for covariates
DD <- "zero"  # matrix(0,mm,1)
dd <- "zero"  # matrix(0,1,wk_last)
## 'RR' is var-cov matrix for obs errors
RR <- "diagonal and equal"


## ----dfa-dfa-proc-eqn--------------------------------------------------------------------------
## number of processes
mm <- 3
## 'BB' is identity: 1's along the diagonal & 0's elsewhere
BB <- "identity"  # diag(mm)
## 'uu' is a column vector of 0's
uu <- "zero"  # matrix(0, mm, 1)
## 'CC' and 'cc' are for covariates
CC <- "zero"  # matrix(0, mm, 1)
cc <- "zero"  # matrix(0, 1, wk_last)
## 'QQ' is identity
QQ <- "identity"  # diag(mm)


## ----dfa-create-model-lists--------------------------------------------------------------------
## list with specifications for model vectors/matrices
mod_list <- list(Z = ZZ, A = aa, D = DD, d = dd, R = RR,
                 B = BB, U = uu, C = CC, c = cc, Q = QQ)
## list with model inits
init_list <- list(x0 = matrix(rep(0, mm), mm, 1))
## list with model control parameters
con_list <- list(maxit = 4000, allow.degen = TRUE)


## ----dfa-fit-dfa-1, cache=TRUE-----------------------------------------------------------------
## fit MARSS
ptm <- proc.time()
dfa_mote <- MARSS(y = mote_dat, model = mod_list, inits = init_list, control = con_list)
proc.time() - ptm
names(dfa_mote)
# the trends are in $states and the SE's are in $states.se

## ----dfa-get-H-inv-----------------------------------------------------------------------------
## get the estimated ZZ
Z_est <- coef(dfa_mote, type = "matrix")$Z
## get the inverse of the rotation matrix
H_inv <- varimax(Z_est)$rotmat


## ----dfa-rotate-Z-x----------------------------------------------------------------------------
## rotate factor loadings
Z_rot = Z_est %*% H_inv   
## rotate processes
proc_rot = solve(H_inv) %*% dfa_mote$states


## ----dfa-plot-dfa1, fig.height=9, fig.width=8, eval=TRUE, fig.cap='Estimated states from the DFA model.'----
ylbl <- motile_epifauna
w_ts <- seq(dim(mote_dat)[2])
layout(matrix(c(1,2,3,4,5,6), mm, 2), widths = c(2,1))
## par(mfcol=c(mm,2), mai = c(0.5,0.5,0.5,0.1), omi = c(0,0,0,0))
par(mai = c(0.5,0.5,0.5,0.1), omi = c(0,0,0,0))
## plot the processes
for(i in 1:mm) {
  ylm <- c(-1,1)*max(abs(proc_rot[i,]))
  ## set up plot area
  plot(w_ts,proc_rot[i,], type = "n", bty = "L",
       ylim = ylm, xlab = "", ylab = "", xaxt = "n")
  ## draw zero-line
  abline(h=0, col="gray")
  ## plot trend line
  lines(w_ts, proc_rot[i,], lwd = 2)
  lines(w_ts, proc_rot[i,], lwd = 2)
  ## add panel labels
  mtext(paste("State",i,"_mot_epi"), side = 3, line = 0.5)
  # axis(1,12*(0:dim(dat_1980)[2])+1,yr_frst+0:dim(dat_1980)[2])
  axis(1,seq(1,42,10),as.character(seq(1982,2023,10)))
}
## plot the loadings
minZ <- 0
ylm <- c(-1,1)*max(abs(Z_rot))
for(i in 1:mm) {
  plot(c(1:N_ts)[abs(Z_rot[,i])>minZ], as.vector(Z_rot[abs(Z_rot[,i])>minZ,i]), type="h",
       lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylm, xlim=c(0.5,N_ts+0.5), 
       col=clr, cex=1.5)
  for(j in 1:N_ts) {
    if(Z_rot[j,i] > minZ) {text(j, -0.03, ylbl[j], srt=90, adj=1, cex=1.2, col=clr[j])}
    if(Z_rot[j,i] < -minZ) {text(j, 0.03, ylbl[j], srt=90, adj=0, cex=1.2, col=clr[j])}
    abline(h=0, lwd=1.5, col="gray")
  } 
  mtext(paste("Factor loadings on state",i),side=3,line=0.5)
}


## ----dfa-xy-states12, height=4, width=5, fig.cap='Cross-correlation plot of the two rotations.'----
par(mai = c(0.9,0.9,0.1,0.1))
ccf(proc_rot[1,],proc_rot[2,], lag.max = 12, main="")


## ----dfa-plot-dfa-fits, fig.height=9, fig.width=8, fig.cap='Data and fits from the DFA model.'----
## get model fits & CI's
mod_fit <- get_DFA_fits(dfa_mote)
## plot the fits
ylbl <- motile_epifauna
par(mfrow = c(4,3), mai = c(0.5,0.7,0.1,0.1), omi = c(0,0,0,0))
for(i in 1:N_ts) {
  up <- mod_fit$up[i,]
  mn <- mod_fit$ex[i,]
  lo <- mod_fit$lo[i,]
  plot(w_ts,mn,xlab = "",ylab=ylbl[i],xaxt = "n",type = "n", cex.lab = 1.2,
       ylim=c(min(lo),max(up)))
  axis(1,seq(1,42,10),as.character(seq(1982,2023,10)))
  points(w_ts,mote_dat[i,], pch=16, col=clr[i])
  lines(w_ts, up, col="darkgray")
  lines(w_ts, mn, col="black", lwd = 2)
  lines(w_ts, lo, col="darkgray")
}



# DFA with covariate

# coldpool data
coldpool_data <- read.csv("C:/Users/andy.whitehouse/Work/Andy/REEM/Eco_Considerations/coldpool.csv")
cp_area <- t(coldpool_data$AREA_LTE2_KM2)
## ----dfa-fit-DFA-covars, cache=TRUE, results='hide'--------------------------------------------
mod_list = list(m = 3, R = "diagonal and unequal")
dfa_apex_cpa <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                  control = con_list, covariates=temp)
