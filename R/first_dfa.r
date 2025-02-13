# learning DFA

library(bayesdfa)
library(ggplot2)
library(dplyr)
library(rstan)
chains = 1
iter = 10

# ebs bottom trawl survey data
ebs_dat_file <- "data/ebs_bts_q1.csv"
ebs_dat <- read.csv(ebs_dat_file)
t(ebs_dat[,-1])
pairs(ebs_dat[,2:20])
pairs(ebs_dat[,21:40])

df = data.frame("time" = rep(1:41,46),"y"=c(ebs_dat[,-1]), 
                "ts"=as.factor(sort(rep(1:46,41))))
ggplot(df, aes(time,y,group=ts,col=ts)) + geom_line() + theme_bw() + 
  xlab("Time") + ylab("Observed data")

q100 <- c("sandlance", "atka", "oth_managed_forage", "capelin", "oth_pelagic_smelt",
          "salmon_returning")
bts_na <- c("sandlance", "atka", "oth_managed_forage", "capelin", 
            "oth_pelagic_smelt", "salmon_returning",
            "kamchatka", "dover_sole")
keepers <- c(2:5,7,9,11,13:19,21:27,29,32:36,42:47)
# plot the data z-scored
ebs_dat_z <- scale(ebs_dat[,2:47])
row.names(ebs_dat_z) <- ebs_dat[,1]
plot(1982:2023, runif(42,0,1), type='n', xlab="year", ylab="z-score",
     ylim=c(min(ebs_dat_z, na.rm=TRUE), max(ebs_dat_z, na.rm=TRUE)))
abline(h=0)
for(i in 1:46){
  lines(1982:2019, ebs_dat_z[1:38,i])
}
for(i in 1:46){
  lines(2021:2023, ebs_dat_z[39:41,i])
}

# 1 MCMC chain with 50 iterations
ebs_f1 <- fit_dfa(
  y = t(ebs_dat[,keepers]), num_trends = 1, scale="zscore",
  iter = iter, chains = chains, thin = 1
)
is_converged(ebs_f1, threshold = 1.05)

# more MCMC chains with more iterations
ebs_f1 <- fit_dfa(
  y = t(ebs_dat[,keepers]), num_trends = 1, scale="zscore",
  iter = 5000, chains = 4, thin = 1
)
is_converged(ebs_f1, threshold = 1.05)

# rotate trends
ebs_r <- rotate_trends(ebs_f1)
names(ebs_r)
plot_trends(ebs_r) + theme_bw()
plot_fitted(ebs_f1) + theme_bw()
plot_loadings(ebs_r) + theme_bw()

# now fit 2 trends
ebs_f2 <- fit_dfa(
  y = t(ebs_dat[,keepers]), num_trends = 2, scale="zscore",
  iter = 16000, chains = 4, thin = 1,
  verbose = TRUE,
  estimate_trend_ar = TRUE#,
  # estimate_trend_ma = TRUE,
  # family = "lognormal"
)
is_converged(ebs_f2, threshold = 1.05)
# 3 trends
ebs_f3 <- fit_dfa(
  y = t(ebs_dat[,keepers]), num_trends = 3, scale="zscore",
  iter = 5000, chains = 4, thin = 1
)
is_converged(ebs_f3, threshold = 1.05)

# rotate trends
ebs_r <- rotate_trends(ebs_f2)
names(ebs_r)
plot_trends(ebs_r) + theme_bw()
plot_fitted(ebs_f2) + theme_bw()


# ---------------------------------------------------------------------------- #
# EBS guilds

# ebs bottom trawl survey data
ebs_guilds_file <- "data/ebs_guilds.csv"
ebs_guilds <- read.csv(ebs_guilds_file)
t(ebs_guilds[,-1])
pairs(ebs_guilds[,-1])

df = data.frame("time" = rep(1:41,4),"y"=c(ebs_guilds[,-1]), 
                "ts"=as.factor(sort(rep(1:4,41))))
ggplot(df, aes(time,y,group=ts,col=ts)) + geom_line() + theme_bw() + 
  xlab("Time") + ylab("Observed data")

# 4 MCMC chain with 5000 iterations
guild_f1 <- fit_dfa(
  y = t(ebs_guilds[,-1]), num_trends = 1, scale="zscore",
  iter = iter, chains = chains, thin = 1
)
is_converged(guild_f1, threshold = 1.05)

# more MCMC chains with more iterations
guild_f1 <- fit_dfa(
  y = t(ebs_guilds[,-1]), num_trends = 2, scale="zscore",
  iter = 5000, chains = 4, thin = 1,
  estimate_trend_ar = TRUE
)
is_converged(guild_f1, threshold = 1.05)

