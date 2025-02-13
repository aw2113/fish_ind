library('MARSS')
library('viridis')
source('get_DFA_fits.r')
source('get_DFA1_fits.r')
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
clr <- viridis(dim(dat)[1])
cnt <- 1
par(mfrow = c(3,3), mai = c(0.5,0.7,0.1,0.1), omi = c(0,0,0,0))
for(i in spp){
  plot(dat[i,],xlab = "",ylab="biomass", bty = "L", xaxt = "n", pch=16, col=clr[cnt], type="b")
  axis(1,1:dim(dat)[2],1982:2023)
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
DD <- "unconstrained"  # unconstrained when we have covariates
dd <- rbind(npi, ice_ext, totc)  # covariate time series
## 'RR' is var-cov matrix for obs errors
RR <- "diagonal and equal"
# RRdu <- "diagonal and unequal"
# RReq <- "equalvarcov"
# RRuneq <- "unconstrained"

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
# mod_list_du <- list(Z = ZZ, A = aa, D = DD, d = dd, R = RRdu,
#                     B = BB, U = uu, C = CC, c = cc, Q = QQ)
# mod_list_eq <- list(Z = ZZ, A = aa, D = DD, d = dd, R = RReq,
#                     B = BB, U = uu, C = CC, c = cc, Q = QQ)
# mod_list_uneq <- list(Z = ZZ, A = aa, D = DD, d = dd, R = RRuneq,
#                       B = BB, U = uu, C = CC, c = cc, Q = QQ)
## list with model inits
init_list <- list(x0 = matrix(rep(0, mm), mm, 1))
## list with model control parameters
con_list <- list(maxit = 10000, allow.degen = TRUE)#, abstol=0.001)

## ----dfa-fit-dfa-1, cache=TRUE-----------------------------------------------------------------
## fit MARSS
ptm <- proc.time()
# if you use form="dfa" or MARSS.dfa(), then the mod_list cannot contain Z, e.g.,
mod_list = list(m = 3, R = "diagonal and equal")
ape_dfa3_cov3 <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                       control = con_list, covariates = dd)
MARSSparamCIs(ape_dfa3_cov3)
ape_dfa_3_de_init <- MARSS(y = dat, model = mod_list, inits = init_list,
                      control = con_list)
MARSSparamCIs(ape_dfa_3_de_init)
# ape_dfa_3_de <- MARSS(y = dat, model = mod_list, control = con_list)
# ape_dfa_3_du <- MARSS(y = dat, model = mod_list_du, inits = init_list, control = con_list)
# ape_dfa_3_eq <- MARSS(y = dat, model = mod_list_eq, inits = init_list, control = con_list)
# ape_dfa_3_uneq <- MARSS(y = dat, model = mod_list_uneq, inits = init_list, control = con_list)
# MARSSparamCIs(ape_dfa_3_de)
proc.time() - ptm
# names(ape_dfa_3_de)
# the trends are in $states and the SE's are in $states.se
# print(c(ape_dfa_3_de$logLik, ape_dfa_3_de$AICc))
# 
## ----dfa-fit-DFA-covars, cache=TRUE, results='hide'--------------------------------------------
mod_list = list(m = 3, R = "diagonal and equal")
dfa3_apex_cpa <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                      control = con_list, covariates=cp_area)
dfa_apex_npi <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                      control = con_list, covariates=npi)
dfa_apex_ice <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                      control = con_list, covariates=ice_ext)
dfa_apex_pdo <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                      control = con_list, covariates=pdo)


## ----dfa-model-selection-----------------------------------------------------------------------
# Compare AICc and logLik
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            logLik=round(c(ape_dfa_3_de$logLik, ape_dfa_3_du$logLik, ape_dfa_3_eq$logLik, ape_dfa_3_uneq$logLik),3)),
      quote=FALSE)
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            AICc=round(c(ape_dfa_3_de$AICc, ape_dfa_3_du$AICc, ape_dfa_3_eq$AICc, ape_dfa_3_uneq$AICc),3)),
      quote=FALSE)



# *****Plot***** ------------------------------------------------------------- #
## ----dfa-get-H-inv-----------------------------------------------------------------------------
## get the estimated ZZ
Z_est <- coef(ape_dfa3_cov3, type = "matrix")$Z
## get the inverse of the rotation matrix
H_inv <- varimax(Z_est)$rotmat

## ----dfa-rotate-Z-x----------------------------------------------------------------------------
## rotate factor loadings
Z_rot = Z_est %*% H_inv   
## rotate processes
proc_rot = solve(H_inv) %*% ape_dfa3_cov3$states

## ----dfa-plot-dfa1, fig.height=9, fig.width=8, eval=TRUE, fig.cap='Estimated states from the DFA model.'----
ylbl <- apex_predators
w_ts <- seq(dim(dat)[2])
# the first four trends
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
  abline(h=c(-0.2,0.2), lty=2, col="gray")
  for(j in 1:N_ts) {
    if(Z_rot[j,i] > minZ) {text(j, -0.03, ylbl[j], srt=90, adj=1, cex=1, col=clr[j])}
    if(Z_rot[j,i] < -minZ) {text(j, 0.03, ylbl[j], srt=90, adj=0, cex=1, col=clr[j])}
    abline(h=0, lwd=1.5, col="gray")
  } 
  mtext(paste("Factor loadings on state",i),side=3,line=0.5)
}

## ----dfa-xy-states12, height=4, width=5, fig.cap='Cross-correlation plot of the two rotations.'----
par(mai = c(0.9,0.9,0.1,0.1))
ccf(proc_rot[1,],proc_rot[2,], lag.max = 12, main="")
acf(proc_rot[1,], type = "covariance")

## ----dfa-plot-dfa-fits, fig.height=9, fig.width=8, fig.cap='Data and fits from the DFA model.'----
## get model fits & CI's
mod_fit <- get_DFA_fits(ape_dfa3_cov3)
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
# ---------------------------------------------------------------------------- #






























# ============================================================================ #
# ----- ----- ***** 2 trends ***** ----- ----- #
# 9 time series
Z_vals <- list("z11",  0  ,
               "z21","z22",
               "z31","z32",
               "z41","z42",
               "z51","z52",
               "z61","z62",
               "z71","z72",
               "z81","z82",
               "z91","z92")
ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 2, byrow = TRUE)
ZZ
## 'aa' is the offset/scaling
aa <- "zero"
## 'DD' and 'd' are for covariates
DD <- "zero"  # matrix(0,mm,1)
dd <- "zero"  # matrix(0,1,wk_last)
## 'RR' is var-cov matrix for obs errors
RR <- "diagonal and equal"
RRdu <- "diagonal and unequal"
RReq <- "equalvarcov"
RRuneq <- "unconstrained"

## ----dfa-dfa-proc-eqn--------------------------------------------------------------------------
## number of processes
mm <- 2
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
mod_list_du <- list(Z = ZZ, A = aa, D = DD, d = dd, R = RRdu,
                    B = BB, U = uu, C = CC, c = cc, Q = QQ)
mod_list_eq <- list(Z = ZZ, A = aa, D = DD, d = dd, R = RReq,
                    B = BB, U = uu, C = CC, c = cc, Q = QQ)
mod_list_uneq <- list(Z = ZZ, A = aa, D = DD, d = dd, R = RRuneq,
                      B = BB, U = uu, C = CC, c = cc, Q = QQ)
## list with model inits
init_list <- list(x0 = matrix(rep(0, mm), mm, 1))
## list with model control parameters
con_list <- list(maxit = 10000, allow.degen = TRUE)

## ----dfa-fit-dfa-1, cache=TRUE-----------------------------------------------------------------
## fit MARSS
ptm <- proc.time()
ape_dfa_2_de <- MARSS(y = dat, model = mod_list, inits = init_list, control = con_list)
ape_dfa_2_du <- MARSS(y = dat, model = mod_list_du, inits = init_list, control = con_list)
ape_dfa_2_eq <- MARSS(y = dat, model = mod_list_eq, inits = init_list, control = con_list)
ape_dfa_2_uneq <- MARSS(y = dat, model = mod_list_uneq, inits = init_list, control = con_list)
proc.time() - ptm
# the trends are in $states and the SE's are in $states.se

## ----dfa-model-selection-----------------------------------------------------------------------
# Compare AICc and logLik
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            logLik=round(c(ape_dfa_2_de$logLik, ape_dfa_2_du$logLik, ape_dfa_2_eq$logLik, ape_dfa_2_uneq$logLik),3)),
      quote=FALSE)
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            AICc=round(c(ape_dfa_2_de$AICc, ape_dfa_2_du$AICc, ape_dfa_2_eq$AICc, ape_dfa_2_uneq$AICc),3)),
      quote=FALSE)


# ============================================================================ #
# ----- ----- ***** 1 trends ***** ----- ----- #
# 9 time series
Z_vals <- list("z11",
               "z21",
               "z31",
               "z41",
               "z51",
               "z61",
               "z71",
               "z81",
               "z91")
ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 1, byrow = TRUE)
ZZ
## 'aa' is the offset/scaling
aa <- "zero"
## 'DD' and 'd' are for covariates
DD <- "zero"  # matrix(0,mm,1)
dd <- "zero"  # matrix(0,1,wk_last)
## 'RR' is var-cov matrix for obs errors
RR <- "diagonal and equal"
RRdu <- "diagonal and unequal"
RReq <- "equalvarcov"
RRuneq <- "unconstrained"

## ----dfa-dfa-proc-eqn--------------------------------------------------------------------------
## number of processes
mm <- 1
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
mod_list_du <- list(Z = ZZ, A = aa, D = DD, d = dd, R = RRdu,
                    B = BB, U = uu, C = CC, c = cc, Q = QQ)
mod_list_eq <- list(Z = ZZ, A = aa, D = DD, d = dd, R = RReq,
                    B = BB, U = uu, C = CC, c = cc, Q = QQ)
mod_list_uneq <- list(Z = ZZ, A = aa, D = DD, d = dd, R = RRuneq,
                      B = BB, U = uu, C = CC, c = cc, Q = QQ)
## list with model inits
init_list <- list(x0 = matrix(rep(0, mm), mm, 1))
## list with model control parameters
con_list <- list(maxit = 10000, allow.degen = TRUE)

## ----dfa-fit-dfa-1, cache=TRUE-----------------------------------------------------------------
## fit MARSS
ptm <- proc.time()
ape_dfa_1_de <- MARSS(y = dat, model = mod_list, inits = init_list, control = con_list)
ape_dfa_1_du <- MARSS(y = dat, model = mod_list_du, inits = init_list, control = con_list)
ape_dfa_1_eq <- MARSS(y = dat, model = mod_list_eq, inits = init_list, control = con_list)
ape_dfa_1_uneq <- MARSS(y = dat, model = mod_list_uneq, inits = init_list, control = con_list)
proc.time() - ptm
# the trends are in $states and the SE's are in $states.se

# Trying to figure out how to get AICb --------------------------------------- #
# ptm <- proc.time()
# ape_dfa_1_de_aicb <- MARSSaic(ape_dfa_1_de,
#                           output = c("AICbp", "boot.params"),
#                           Options = list(nboot = 10, silent = TRUE))
# proc.time() - ptm
# residuals(ape_dfa_1_de_aicb)
## ----dfa-model-selection-----------------------------------------------------------------------
# Compare AICc and logLik
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            logLik=round(c(ape_dfa_1_de$logLik, ape_dfa_1_du$logLik, ape_dfa_1_eq$logLik, ape_dfa_1_uneq$logLik),3)),
      quote=FALSE)
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            AICc=round(c(ape_dfa_1_de$AICc, ape_dfa_1_du$AICc, ape_dfa_1_eq$AICc, ape_dfa_1_uneq$AICc),3)),
      quote=FALSE)
