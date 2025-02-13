library('MARSS')
library('viridis')
source('get_DFA_fits.r')
source('get_DFA1_fits.r')
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
clr <- viridis(dim(dat)[1])
cnt <- 1
par(mfrow = c(3,4), mai = c(0.5,0.7,0.1,0.1), omi = c(0,0,0,0))
for(i in spp){
  plot(dat[i,],xlab = "",ylab="biomass", bty = "L", xaxt = "n", pch=16, col=clr[cnt], type="b")
  axis(1,1:dim(dat)[2],1982:2023)
  title(i)
  cnt <- cnt + 1
}


# ============================================================================ #
# ----- ----- ***** 1 trends ***** ----- ----- #
# 11 time series
Z_vals <- list("z11",
               "z21",
               "z31",
               "z41",
               "z51",
               "z61",
               "z71",
               "z81",
               "z91",
               "z101",
               "z111")
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
pel_dfa_1_de <- MARSS(y = dat, model = mod_list, inits = init_list, control = con_list)
pel_dfa_1_du <- MARSS(y = dat, model = mod_list_du, inits = init_list, control = con_list)
pel_dfa_1_eq <- MARSS(y = dat, model = mod_list_eq, inits = init_list, control = con_list)
pel_dfa_1_uneq <- MARSS(y = dat, model = mod_list_uneq, inits = init_list, control = con_list)
proc.time() - ptm
# the trends are in $states and the SE's are in $states.se
MARSSparamCIs(pel_dfa_1_de)

# Trying to figure out how to get AICb --------------------------------------- #
# ptm <- proc.time()
# pel_dfa_1_de_aicb <- MARSSaic(pel_dfa_1_de,
#                           output = c("AICbp", "boot.params"),
#                           Options = list(nboot = 10, silent = TRUE))
# proc.time() - ptm
# residuals(pel_dfa_1_de_aicb)
## ----dfa-model-selection-----------------------------------------------------------------------
# Compare AICc and logLik
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            logLik=round(c(pel_dfa_1_de$logLik, pel_dfa_1_du$logLik, pel_dfa_1_eq$logLik, pel_dfa_1_uneq$logLik),3),
            AICc=round(c(pel_dfa_1_de$AICc, pel_dfa_1_du$AICc, pel_dfa_1_eq$AICc, pel_dfa_1_uneq$AICc),3)),
      quote=FALSE)
# print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
#             AICc=round(c(pel_dfa_1_de$AICc, pel_dfa_1_du$AICc, pel_dfa_1_eq$AICc, pel_dfa_1_uneq$AICc),3)),
#       quote=FALSE)


# ============================================================================ #
# ----- ----- ***** 2 trends ***** ----- ----- #
# 11 time series
Z_vals <- list("z11",  0  ,
               "z21","z22",
               "z31","z32",
               "z41","z42",
               "z51","z52",
               "z61","z62",
               "z71","z72",
               "z81","z82",
               "z91","z92",
               "z101","z102",
               "z111","z112")
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
pel_dfa_2_de <- MARSS(y = dat, model = mod_list, inits = init_list, control = con_list)
pel_dfa_2_du <- MARSS(y = dat, model = mod_list_du, inits = init_list, control = con_list)
pel_dfa_2_eq <- MARSS(y = dat, model = mod_list_eq, inits = init_list, control = con_list)
pel_dfa_2_uneq <- MARSS(y = dat, model = mod_list_uneq, inits = init_list, control = con_list)
proc.time() - ptm
# the trends are in $states and the SE's are in $states.se

## ----dfa-model-selection-----------------------------------------------------------------------
# Compare AICc and logLik
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            logLik=round(c(pel_dfa_2_de$logLik, pel_dfa_2_du$logLik, pel_dfa_2_eq$logLik, pel_dfa_2_uneq$logLik),3)),
      quote=FALSE)
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            AICc=round(c(pel_dfa_2_de$AICc, pel_dfa_2_du$AICc, pel_dfa_2_eq$AICc, pel_dfa_2_uneq$AICc),3)),
      quote=FALSE)


# ============================================================================ #
# ----- ----- ***** 3 trends ***** ----- ----- #
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
RRdu <- "diagonal and unequal"
RReq <- "equalvarcov"
RRuneq <- "unconstrained"

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
mod_list_du <- list(Z = ZZ, A = aa, D = DD, d = dd, R = RRdu,
                    B = BB, U = uu, C = CC, c = cc, Q = QQ)
mod_list_eq <- list(Z = ZZ, A = aa, D = DD, d = dd, R = RReq,
                    B = BB, U = uu, C = CC, c = cc, Q = QQ)
mod_list_uneq <- list(Z = ZZ, A = aa, D = DD, d = dd, R = RRuneq,
                      B = BB, U = uu, C = CC, c = cc, Q = QQ)
## list with model inits
init_list <- list(x0 = matrix(rep(0, mm), mm, 1))
## list with model control parameters
con_list <- list(maxit = 10000, allow.degen = TRUE)#, abstol=0.001)

## ----dfa-fit-dfa-1, cache=TRUE-----------------------------------------------------------------
## fit MARSS
ptm <- proc.time()
pel_dfa_3_de <- MARSS(y = dat, model = mod_list, inits = init_list, control = con_list)
pel_dfa_3_du <- MARSS(y = dat, model = mod_list_du, inits = init_list, control = con_list)
pel_dfa_3_eq <- MARSS(y = dat, model = mod_list_eq, inits = init_list, control = con_list)
pel_dfa_3_uneq <- MARSS(y = dat, model = mod_list_uneq, inits = init_list, control = con_list)
proc.time() - ptm
# the trends are in $states and the SE's are in $states.se
# print(c(pel_dfa_3_de$logLik, pel_dfa_3_de$AICc))

# residuals(pel_dfa_3_eq)
## ----dfa-model-selection-----------------------------------------------------------------------
# Compare AICc and logLik
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            logLik=round(c(pel_dfa_3_de$logLik, pel_dfa_3_du$logLik, pel_dfa_3_eq$logLik, pel_dfa_3_uneq$logLik),3)),
      quote=FALSE)
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            AICc=round(c(pel_dfa_3_de$AICc, pel_dfa_3_du$AICc, pel_dfa_3_eq$AICc, pel_dfa_3_uneq$AICc),3)),
      quote=FALSE)


# ============================================================================ #
# ----- ----- ***** 4 trends ***** ----- ----- #
# 11 time series
Z_vals <- list("z11",  0  ,  0  ,  0  ,
               "z21","z22",  0  ,  0  ,
               "z31","z32","z33",  0  ,
               "z41","z42","z43","z44",
               "z51","z52","z53","z54",
               "z61","z62","z63","z64",
               "z71","z72","z73","z74",
               "z81","z82","z83","z84",
               "z91","z92","z93","z94",
               "z101","z102","z103","z104",
               "z111","z112","z113","z114")
ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 4, byrow = TRUE)
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
mm <- 4
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
pel_dfa_4_de <- MARSS(y = dat, model = mod_list, inits = init_list, control = con_list)
pel_dfa_4_du <- MARSS(y = dat, model = mod_list_du, inits = init_list, control = con_list)
pel_dfa_4_eq <- MARSS(y = dat, model = mod_list_eq, inits = init_list, control = con_list)
pel_dfa_4_uneq <- MARSS(y = dat, model = mod_list_uneq, inits = init_list, control = con_list)
proc.time() - ptm
# the trends are in $states and the SE's are in $states.se

## ----dfa-model-selection-----------------------------------------------------------------------
# Compare AICc and logLik
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            logLik=round(c(pel_dfa_4_de$logLik, pel_dfa_4_du$logLik, pel_dfa_4_eq$logLik, pel_dfa_4_uneq$logLik),3)),
      quote=FALSE)
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            AICc=round(c(pel_dfa_4_de$AICc, pel_dfa_4_du$AICc, pel_dfa_4_eq$AICc, pel_dfa_4_uneq$AICc),3)),
      quote=FALSE)

# ============================================================================ #
# ----- ----- ***** 5 trends ***** ----- ----- #
# 11 time series
Z_vals <- list("z11",  0  ,  0  ,  0  ,  0  ,
               "z21","z22",  0  ,  0  ,  0  ,
               "z31","z32","z33",  0  ,  0  ,
               "z41","z42","z43","z44",  0  ,
               "z51","z52","z53","z54","z55",
               "z61","z62","z63","z64","z65",
               "z71","z72","z73","z74","z75",
               "z81","z82","z83","z84","z85",
               "z91","z92","z93","z94","z95",
               "z101","z102","z103","z104","z105",
               "z111","z112","z113","z114","z115")
ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 5, byrow = TRUE)
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
mm <- 5
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
pel_dfa_5_de <- MARSS(y = dat, model = mod_list, inits = init_list, control = con_list)
pel_dfa_5_du <- MARSS(y = dat, model = mod_list_du, inits = init_list, control = con_list)
pel_dfa_5_eq <- MARSS(y = dat, model = mod_list_eq, inits = init_list, control = con_list)
pel_dfa_5_uneq <- MARSS(y = dat, model = mod_list_uneq, inits = init_list, control = con_list)
proc.time() - ptm
# the trends are in $states and the SE's are in $states.se

## ----dfa-model-selection-----------------------------------------------------------------------
# Compare AICc and logLik
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            logLik=round(c(pel_dfa_5_de$logLik, pel_dfa_5_du$logLik, pel_dfa_5_eq$logLik, pel_dfa_5_uneq$logLik),3)),
      quote=FALSE)
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            AICc=round(c(pel_dfa_5_de$AICc, pel_dfa_5_du$AICc, pel_dfa_5_eq$AICc, pel_dfa_5_uneq$AICc),3)),
      quote=FALSE)

# ============================================================================ #
# ----- ----- ***** 6 trends ***** ----- ----- #
# 11 time series
Z_vals <- list("z11",  0  ,  0  ,  0  ,  0  ,  0  ,
               "z21","z22",  0  ,  0  ,  0  ,  0  ,
               "z31","z32","z33",  0  ,  0  ,  0  ,
               "z41","z42","z43","z44",  0  ,  0  ,
               "z51","z52","z53","z54","z55",  0  ,
               "z61","z62","z63","z64","z65","z66",
               "z71","z72","z73","z74","z75","z76",
               "z81","z82","z83","z84","z85","z86",
               "z91","z92","z93","z94","z95","z96",
               "z101","z102","z103","z104","z105","z106",
               "z111","z112","z113","z114","z115","z116")
ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 6, byrow = TRUE)
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
mm <- 6
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
pel_dfa_6_de <- MARSS(y = dat, model = mod_list, inits = init_list, control = con_list)
pel_dfa_6_du <- MARSS(y = dat, model = mod_list_du, inits = init_list, control = con_list)
pel_dfa_6_eq <- MARSS(y = dat, model = mod_list_eq, inits = init_list, control = con_list)
pel_dfa_6_uneq <- MARSS(y = dat, model = mod_list_uneq, inits = init_list, control = con_list)
proc.time() - ptm
# the trends are in $states and the SE's are in $states.se

## ----dfa-model-selection-----------------------------------------------------------------------
# Compare AICc and logLik
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            logLik=round(c(pel_dfa_6_de$logLik, pel_dfa_6_du$logLik, pel_dfa_6_eq$logLik, pel_dfa_6_uneq$logLik),3)),
      quote=FALSE)
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            AICc=round(c(pel_dfa_6_de$AICc, pel_dfa_6_du$AICc, pel_dfa_6_eq$AICc, pel_dfa_6_uneq$AICc),3)),
      quote=FALSE)

# ============================================================================ #
# ----- ----- ***** 7 trends ***** ----- ----- #
# 11 time series
Z_vals <- list("z11",  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,
               "z21","z22",  0  ,  0  ,  0  ,  0  ,  0  ,
               "z31","z32","z33",  0  ,  0  ,  0  ,  0  ,
               "z41","z42","z43","z44",  0  ,  0  ,  0  ,
               "z51","z52","z53","z54","z55",  0  ,  0  ,
               "z61","z62","z63","z64","z65","z66",  0  ,
               "z71","z72","z73","z74","z75","z76","z77",
               "z81","z82","z83","z84","z85","z86","z87",
               "z91","z92","z93","z94","z95","z96","z97",
               "z101","z102","z103","z104","z105","z106","z107",
               "z111","z112","z113","z114","z115","z116","z117")
ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 7, byrow = TRUE)
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
mm <- 7
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
pel_dfa_7_de <- MARSS(y = dat, model = mod_list, inits = init_list, control = con_list)
pel_dfa_7_du <- MARSS(y = dat, model = mod_list_du, inits = init_list, control = con_list)
pel_dfa_7_eq <- MARSS(y = dat, model = mod_list_eq, inits = init_list, control = con_list)
pel_dfa_7_uneq <- MARSS(y = dat, model = mod_list_uneq, inits = init_list, control = con_list)
proc.time() - ptm
# the trends are in $states and the SE's are in $states.se

# Trying to figure out how to get AICb --------------------------------------- #
# ptm <- proc.time()
# pel_dfa_7_eq_aicb <- MARSSaic(pel_dfa_7_eq,
#                           output = c("AICbp", "boot.params"), 
#                           Options = list(nboot = 10, silent = TRUE))
# proc.time() - ptm


## ----dfa-model-selection-----------------------------------------------------------------------
# Compare AICc and logLik
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            logLik=round(c(pel_dfa_7_de$logLik, pel_dfa_7_du$logLik, pel_dfa_7_eq$logLik, pel_dfa_7_uneq$logLik),3)),
      quote=FALSE)
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            AICc=round(c(pel_dfa_7_de$AICc, pel_dfa_7_du$AICc, pel_dfa_7_eq$AICc, pel_dfa_7_uneq$AICc),3)),
      quote=FALSE)




# ============================================================================ #
# ----- ----- ***** 8 trends ***** ----- ----- #
# 11 time series
Z_vals <- list("z11" ,  0   ,  0   ,  0   ,  0   ,  0   ,  0   ,  0  ,
               "z21" ,"z22" ,  0   ,  0   ,  0   ,  0   ,  0   ,  0  ,
               "z31" ,"z32" ,"z33" ,  0   ,  0   ,  0   ,  0   ,  0  ,
               "z41" ,"z42" ,"z43" ,"z44" ,  0   ,  0   ,  0   ,  0  ,
               "z51" ,"z52" ,"z53" ,"z54" ,"z55" ,  0   ,  0   ,  0  ,
               "z61" ,"z62" ,"z63" ,"z64" ,"z65" ,"z66" ,  0   ,  0  ,
               "z71" ,"z72" ,"z73" ,"z74" ,"z75" ,"z76" ,"z77" ,  0  ,
               "z81" ,"z82" ,"z83" ,"z84" ,"z85" ,"z86" ,"z87" ,"z88",
               "z91" ,"z92" ,"z93" ,"z94" ,"z95" ,"z96" ,"z97" ,"z98",
               "z101","z102","z103","z104","z105","z106","z107","z108",
               "z111","z112","z113","z114","z115","z116","z117","z118")
ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 8, byrow = TRUE)
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
mm <- 8
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
pel_dfa_8_de <- MARSS(y = dat, model = mod_list, inits = init_list, control = con_list)
pel_dfa_8_du <- MARSS(y = dat, model = mod_list_du, inits = init_list, control = con_list)
pel_dfa_8_eq <- MARSS(y = dat, model = mod_list_eq, inits = init_list, control = con_list)
pel_dfa_8_uneq <- MARSS(y = dat, model = mod_list_uneq, inits = init_list, control = con_list)
proc.time() - ptm
# the trends are in $states and the SE's are in $states.se

## ----dfa-model-selection-----------------------------------------------------------------------
# Compare AICc and logLik
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            logLik=round(c(pel_dfa_8_de$logLik, pel_dfa_8_du$logLik, pel_dfa_8_eq$logLik, pel_dfa_8_uneq$logLik),3)),
      quote=FALSE)
print(cbind(model=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            AICc=round(c(pel_dfa_8_de$AICc, pel_dfa_8_du$AICc, pel_dfa_8_eq$AICc, pel_dfa_8_uneq$AICc),3)),
      quote=FALSE)

# ============================================================================ #
# ----- ----- ***** 9 trends ***** ----- ----- #
# 11 time series
Z_vals <- list("z11" ,  0   ,  0   ,  0   ,  0   ,  0   ,  0   ,  0   ,  0  ,
               "z21" ,"z22" ,  0   ,  0   ,  0   ,  0   ,  0   ,  0   ,  0  ,
               "z31" ,"z32" ,"z33" ,  0   ,  0   ,  0   ,  0   ,  0   ,  0  ,
               "z41" ,"z42" ,"z43" ,"z44" ,  0   ,  0   ,  0   ,  0   ,  0  ,
               "z51" ,"z52" ,"z53" ,"z54" ,"z55" ,  0   ,  0   ,  0   ,  0  ,
               "z61" ,"z62" ,"z63" ,"z64" ,"z65" ,"z66" ,  0   ,  0   ,  0  ,
               "z71" ,"z72" ,"z73" ,"z74" ,"z75" ,"z76" ,"z77" ,  0   ,  0  ,
               "z81" ,"z82" ,"z83" ,"z84" ,"z85" ,"z86" ,"z87" ,"z88" ,  0  ,
               "z91" ,"z92" ,"z93" ,"z94" ,"z95" ,"z96" ,"z97" ,"z98" ,"z99",
               "z101","z102","z103","z104","z105","z106","z107","z108","z109",
               "z111","z112","z113","z114","z115","z116","z117","z118","z119")
ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 9, byrow = TRUE)
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
mm <- 9
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
pel_dfa_9_de <- MARSS(y = dat, model = mod_list, inits = init_list, control = con_list)
pel_dfa_9_du <- MARSS(y = dat, model = mod_list_du, inits = init_list, control = con_list)
pel_dfa_9_eq <- MARSS(y = dat, model = mod_list_eq, inits = init_list, control = con_list)
pel_dfa_9_uneq <- MARSS(y = dat, model = mod_list_uneq, inits = init_list, control = con_list)
proc.time() - ptm
# the trends are in $states and the SE's are in $states.se

## ----dfa-model-selection-----------------------------------------------------------------------
# Compare AICc and logLik
print(cbind(trends=rep(9,4),
            varcovar=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            logLik=round(c(pel_dfa_9_de$logLik, pel_dfa_9_du$logLik, pel_dfa_9_eq$logLik, pel_dfa_9_uneq$logLik),3),
            AICc=round(c(pel_dfa_9_de$AICc, pel_dfa_9_du$AICc, pel_dfa_9_eq$AICc, pel_dfa_9_uneq$AICc),3),
            converg=c(pel_dfa_9_de$convergence, pel_dfa_9_du$convergence, pel_dfa_9_eq$convergence, pel_dfa_9_uneq$convergence)),
      quote=FALSE)


# ============================================================================ #
# ----- ----- ***** 10 trends ***** ----- ----- #
# 11 time series
Z_vals <- list("z11" ,  0   ,  0   ,  0   ,  0   ,  0   ,  0   ,  0   ,  0   ,   0   ,
               "z21" ,"z22" ,  0   ,  0   ,  0   ,  0   ,  0   ,  0   ,  0   ,   0   ,
               "z31" ,"z32" ,"z33" ,  0   ,  0   ,  0   ,  0   ,  0   ,  0   ,   0   ,
               "z41" ,"z42" ,"z43" ,"z44" ,  0   ,  0   ,  0   ,  0   ,  0   ,   0   ,
               "z51" ,"z52" ,"z53" ,"z54" ,"z55" ,  0   ,  0   ,  0   ,  0   ,   0   ,
               "z61" ,"z62" ,"z63" ,"z64" ,"z65" ,"z66" ,  0   ,  0   ,  0   ,   0   ,
               "z71" ,"z72" ,"z73" ,"z74" ,"z75" ,"z76" ,"z77" ,  0   ,  0   ,   0   ,
               "z81" ,"z82" ,"z83" ,"z84" ,"z85" ,"z86" ,"z87" ,"z88" ,  0   ,   0   ,
               "z91" ,"z92" ,"z93" ,"z94" ,"z95" ,"z96" ,"z97" ,"z98" ,"z99" ,   0   ,
               "z101","z102","z103","z104","z105","z106","z107","z108","z109","z1010",
               "z111","z112","z113","z114","z115","z116","z117","z118","z119","z1110")
ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 10, byrow = TRUE)
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
mm <- 10
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
pel_dfa_10_de <- MARSS(y = dat, model = mod_list, inits = init_list, control = con_list)
pel_dfa_10_du <- MARSS(y = dat, model = mod_list_du, inits = init_list, control = con_list)
pel_dfa_10_eq <- MARSS(y = dat, model = mod_list_eq, inits = init_list, control = con_list)
pel_dfa_10_uneq <- MARSS(y = dat, model = mod_list_uneq, inits = init_list, control = con_list)
proc.time() - ptm
# the trends are in $states and the SE's are in $states.se

## ----dfa-model-selection-----------------------------------------------------------------------
# Compare AICc and logLik
print(cbind(trends=rep(mm,4),
            varcovar=c("diag equal", "diag unequal", "equalvarcov", "unconstrained"),
            logLik=round(c(pel_dfa_10_de$logLik, pel_dfa_10_du$logLik, pel_dfa_10_eq$logLik, pel_dfa_10_uneq$logLik),3),
            AICc=round(c(pel_dfa_10_de$AICc, pel_dfa_10_du$AICc, pel_dfa_10_eq$AICc, pel_dfa_10_uneq$AICc),3),
            converg=c(pel_dfa_10_de$convergence, pel_dfa_10_du$convergence, pel_dfa_10_eq$convergence, pel_dfa_10_uneq$convergence)),
      quote=FALSE)



# ============================================================================ #
# 1 trend is best with diagonal and equal and equalvarcov is second best
## number of processes
mm <- 1

# *****Plot***** ------------------------------------------------------------- #
# diagonal and equal
## get the estimated ZZ
Z_est <- coef(pel_dfa_1_de, type = "matrix")$Z
# There's no varimax() rotation because one trend
## get the inverse of the rotation matrix
# H_inv <- varimax(Z_est)$rotmat

## ----dfa-rotate-Z-x----------------------------------------------------------------------------
## rotate factor loadings
# Z_rot = Z_est %*% H_inv   
## rotate processes
# proc_rot = solve(H_inv) %*% pel_dfa_1_eq$states
proc_rot = Z_est %*% pel_dfa_1_de$states

## ----dfa-plot-dfa1, fig.height=9, fig.width=8, eval=TRUE, fig.cap='Estimated states from the DFA model.'----
ylbl <- pelagic_foragers
w_ts <- seq(dim(dat)[2])
layout(matrix(c(1,2), mm, 2), widths = c(2,1))
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
  mtext(paste("State",i), side = 3, line = 0.5)
  # axis(1,12*(0:dim(dat_1980)[2])+1,yr_frst+0:dim(dat_1980)[2])
  axis(1,seq(1,42,10),as.character(seq(1982,2023,10)))
}
## plot the loadings
minZ <- 0
ylm <- c(-1,1)*max(abs(proc_rot))
for(i in 1:mm) {
  plot(c(1:N_ts)[abs(Z_est[,i])>minZ], as.vector(Z_est[abs(Z_est[,i])>minZ,i]), type="h",
       lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylm, xlim=c(0.5,N_ts+0.5), 
       col=clr, cex=1.5)
  for(j in 1:N_ts) {
    if(Z_est[j,i] > minZ) {text(j, -0.03, ylbl[j], srt=90, adj=1, cex=1, col=clr[j])}
    if(Z_est[j,i] < -minZ) {text(j, 0.03, ylbl[j], srt=90, adj=0, cex=1, col=clr[j])}
    abline(h=0, lwd=1.5, col="gray")
  } 
  mtext(paste("Factor loadings on state",i),side=3,line=0.5)
}
abline(h=c(-0.2,0.2), lty=2, col="gray")

## ----dfa-xy-states12, height=4, width=5, fig.cap='Cross-correlation plot of the two rotations.'----
# par(mai = c(0.9,0.9,0.1,0.1))
# ccf(proc_rot[1,],proc_rot[2,], lag.max = 12, main="")

## ----dfa-plot-dfa-fits, fig.height=9, fig.width=8, fig.cap='Data and fits from the DFA model.'----
## get model fits & CI's
# get_DFA_fits() does not work with one trend
# mod_fit <- get_DFA_fits(pel_dfa_1_de)
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

mod_fit <- get_DFA1_fits(pel_dfa_1_de)

## plot the fits
ylbl <- pelagic_foragers
par(mfrow = c(3,4), mai = c(0.5,0.7,0.1,0.1), omi = c(0,0,0,0))
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


# *****Plot***** ------------------------------------------------------------- #
# diagonal and equal
## get the estimated ZZ
Z_est <- coef(pel_dfa_1_eq, type = "matrix")$Z
# There's no varimax() rotation because one trend
## get the inverse of the rotation matrix
# H_inv <- varimax(Z_est)$rotmat

## ----dfa-rotate-Z-x----------------------------------------------------------------------------
## rotate factor loadings
# Z_rot = Z_est %*% H_inv   
## rotate processes
# proc_rot = solve(H_inv) %*% pel_dfa_1_eq$states
proc_rot = Z_est %*% pel_dfa_1_eq$states

## ----dfa-plot-dfa1, fig.height=9, fig.width=8, eval=TRUE, fig.cap='Estimated states from the DFA model.'----
ylbl <- pelagic_foragers
w_ts <- seq(dim(dat)[2])
layout(matrix(c(1,2), mm, 2), widths = c(2,1))
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
  mtext(paste("State",i,"_pelagic"), side = 3, line = 0.5)
  # axis(1,12*(0:dim(dat_1980)[2])+1,yr_frst+0:dim(dat_1980)[2])
  axis(1,seq(1,42,10),as.character(seq(1982,2023,10)))
}
## plot the loadings
minZ <- 0
ylm <- c(-1,1)*max(abs(Z_rot))
for(i in 1:mm) {
  plot(c(1:N_ts)[abs(Z_est[,i])>minZ], as.vector(Z_est[abs(Z_est[,i])>minZ,i]), type="h",
       lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylm, xlim=c(0.5,N_ts+0.5), 
       col=clr, cex=1.5)
  for(j in 1:N_ts) {
    if(Z_est[j,i] > minZ) {text(j, -0.03, ylbl[j], srt=90, adj=1, cex=1.2, col=clr[j])}
    if(Z_est[j,i] < -minZ) {text(j, 0.03, ylbl[j], srt=90, adj=0, cex=1.2, col=clr[j])}
    abline(h=0, lwd=1.5, col="gray")
  } 
  mtext(paste("Factor loadings on state",i),side=3,line=0.5)
}
abline(h=c(-0.2,0.2), lty=2, col="gray50")

## ----dfa-xy-states12, height=4, width=5, fig.cap='Cross-correlation plot of the two rotations.'----
# par(mai = c(0.9,0.9,0.1,0.1))
# ccf(proc_rot[1,],proc_rot[2,], lag.max = 12, main="")

## ----dfa-plot-dfa-fits, fig.height=9, fig.width=8, fig.cap='Data and fits from the DFA model.'----
## get model fits & CI's
# get_DFA_fits() does not work with one trend
# mod_fit <- get_DFA_fits(pel_dfa_1_eq)
mod_fit <- get_DFA1_fits(pel_dfa_1_eq)

## plot the fits
ylbl <- pelagic_foragers
par(mfrow = c(3,4), mai = c(0.5,0.7,0.1,0.1), omi = c(0,0,0,0))
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


#==============================================================================#
# Plot the two trend models

## number of processes
mm <- 2

# *****Plot***** ------------------------------------------------------------- #
# equalvarcov
## get the estimated ZZ
Z_est <- coef(pel_dfa_2_eq, type = "matrix")$Z
## get the inverse of the rotation matrix
H_inv <- varimax(Z_est)$rotmat

## ----dfa-rotate-Z-x----------------------------------------------------------------------------
## rotate factor loadings
Z_rot = Z_est %*% H_inv   
## rotate processes
proc_rot = solve(H_inv) %*% pel_dfa_2_eq$states

## ----dfa-plot-dfa1, fig.height=9, fig.width=8, eval=TRUE, fig.cap='Estimated states from the DFA model.'----
ylbl <- pelagic_foragers
w_ts <- seq(dim(dat)[2])
layout(matrix(c(1,2,3,4), mm, 2), widths = c(2,1))
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
  mtext(paste("State",i,"_pelagic"), side = 3, line = 0.5)
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
  abline(h=c(-0.2,0.2), lty=2, col="gray50")
  for(j in 1:N_ts) {
    if(Z_rot[j,i] > minZ) {text(j, -0.03, ylbl[j], srt=90, adj=1, cex=1.2, col=clr[j])}
    if(Z_rot[j,i] < -minZ) {text(j, 0.03, ylbl[j], srt=90, adj=0, cex=1.2, col=clr[j])}
    abline(h=0, lwd=1.5, col="gray")
  }
  mtext(paste("Factor loadings on state",i),side=3,line=0.5)
}

## ----dfa-xy-states12, height=4, width=5, fig.cap='Cross-correlation plot of the two rotations.'----
# par(mai = c(0.9,0.9,0.1,0.1))
# ccf(proc_rot[1,],proc_rot[2,], lag.max = 12, main="")

## ----dfa-plot-dfa-fits, fig.height=9, fig.width=8, fig.cap='Data and fits from the DFA model.'----
## get model fits & CI's
mod_fit <- get_DFA_fits(pel_dfa_2_eq)
## plot the fits
ylbl <- pelagic_foragers
par(mfrow = c(3,4), mai = c(0.5,0.7,0.1,0.1), omi = c(0,0,0,0))
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
