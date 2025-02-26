setwd("//aa.ad.epa.gov/ORD/RTP/USERS/A-D/CLowe/Net MyDocuments/R Data/Bayes_SQ")

library(rjags)
library(R2jags)
library(readxl)

#######
#RF vs IE simple regression
#http://biometry.github.io/APES//LectureNotes/StatsCafe/Linear_models_jags.html

set.seed(42) # Set a random seed for reproducibility of the simulation

# samplesize <- 30 # Number of data points
# b_length <- sort(rnorm(samplesize)) # Body length (explanatory variable)

int_true <- 3000 # True intercept
slope_true <- 1000 # True slope
mu <- int_true + slope_true * b_length # True means of normal distributions
sigma <- 50 # True standard deviation of normal distributions

# b_mass <- rnorm(samplesize, mean = mu, sd = sigma) # Body mass (response variable)
# 
# snakes1 <- data.frame(b_length = b_length, b_mass = b_mass)
# head(snakes1)

# df_ESI_POS <- read_excel("//aa.ad.epa.gov/ORD/RTP/USERS/A-D/CLowe/Net MyDocuments/R Data/Bayes_SQ/Supplemental_Tables_v3.xlsx",
df_ESI_POS <- read_excel(paste0(datadir2,'/Supplemental_Tables_v3.xlsx'), skip=1,
                         sheet = "Table S3", col_types = c("text", 
                                                           "text", "text", "text", "text", "text", 
                                                           "text", "text", "text", "text", "numeric", 
                                                           "numeric", "numeric", "numeric", 
                                                           "numeric", "numeric", "numeric", 
                                                           "numeric", "numeric", "numeric", 
                                                           "numeric", "numeric", "numeric", 
                                                           "numeric", "numeric", "numeric", 
                                                           "numeric", "numeric"))

df1 <- data.frame(df_ESI_POS$`Transformed Response Factor (RF(1/6))`,df_ESI_POS$`log10 Predicted IE`)
colnames(df1) <- c("RF(1/6)","log10Predicted_IE")
RF <- df1$`RF(1/6)`
IE <- df1$log10Predicted_IE

jagsdata_s1 <- with(df1, list(IE = IE, RF = RF, N = length(IE)))

lm1_jags <- function(){
  # Likelihood:
  for (i in 1:N){
    IE[i] ~ dnorm(mu[i], tau) # tau is precision (1 / variance)
    mu[i] <- alpha + beta * RF[i]
  }
  # Priors:
  alpha ~ dnorm(0, 0.01) # intercept
  beta ~ dnorm(0, 0.01) # slope
  sigma ~ dunif(0, 100) # standard deviation
  tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
}

init_values <- function(){
  list(alpha = rnorm(1), beta = rnorm(1), sigma = runif(1))
}

params <- c("alpha", "beta", "sigma")

fit_lm1 <- jags(data = jagsdata_s1, 
                inits = init_values, 
                parameters.to.save = params, 
                model.file = lm1_jags,
                n.chains = 3, 
                n.iter = 12000, 
                n.burnin = 2000, 
                n.thin = 10, 
                DIC = F)

fit_lm1

traceplot(fit_lm1, mfrow = c(2, 2), ask = F)

plot(fit_lm1)

lm1_mcmc <- as.mcmc(fit_lm1)
plot(lm1_mcmc)

nvalues <- 100
RF_new <- seq(min(df1$RF), max(df1$RF), length.out = nvalues)

lm1_mcmc_combi <- as.mcmc(rbind(lm1_mcmc[[1]], lm1_mcmc[[2]], lm1_mcmc[[3]]))

pred_mean_mean <- mean(lm1_mcmc_combi[, "alpha"]) + RF_new * mean(lm1_mcmc_combi[, "beta"])

pred_mean_dist <- matrix(NA, nrow = nrow(lm1_mcmc_combi), ncol = nvalues)
for (i in 1:nrow(pred_mean_dist)){
  pred_mean_dist[i,] <- lm1_mcmc_combi[i,"alpha"] + RF_new * lm1_mcmc_combi[i,"beta"]
}
credible_lower <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.005)
credible_upper <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.995)

lm1_mcmc_combi_rep <- do.call(rbind, rep(list(lm1_mcmc_combi), 50)) # replication

# Draw random values for all parameter combinations (rows) and body length values (columns):
pred_data_dist <- matrix(NA, nrow = nrow(lm1_mcmc_combi_rep), ncol = nvalues)
for (i in 1:nrow(pred_data_dist)){
  pred_data_dist[i,] <- lm1_mcmc_combi_rep[i,"alpha"] + RF_new * lm1_mcmc_combi_rep[i,"beta"] +
    rnorm(nvalues, mean = 0, sd = lm1_mcmc_combi_rep[i, "sigma"])
}

# Calculate quantiles:
uncertain_lower <- apply(pred_data_dist, MARGIN = 2, quantile, prob = 0.005)
uncertain_upper <- apply(pred_data_dist, MARGIN = 2, quantile, prob = 0.995)

plot(IE ~ RF, data = df1,xlab = "Response Factor^(1/6)",ylab = "log10(Predicted IE)")
lines(RF_new, pred_mean_mean)
lines(RF_new, credible_lower, lty = 2)
lines(RF_new, credible_upper, lty = 2)
lines(RF_new, uncertain_lower, lty = 2, col = "red")
lines(RF_new, uncertain_upper, lty = 2, col = "red")

#######
#IE vs RF flipped regression
#http://biometry.github.io/APES//LectureNotes/StatsCafe/Linear_models_jags.html

set.seed(42) # Set a random seed for reproducibility of the simulation

# samplesize <- 30 # Number of data points
# b_length <- sort(rnorm(samplesize)) # Body length (explanatory variable)

int_true <- 3000 # True intercept
slope_true <- 1000 # True slope
mu <- int_true + slope_true * b_length # True means of normal distributions
sigma <- 50 # True standard deviation of normal distributions

# b_mass <- rnorm(samplesize, mean = mu, sd = sigma) # Body mass (response variable)
# 
# snakes1 <- data.frame(b_length = b_length, b_mass = b_mass)
# head(snakes1)

df_ESI_POS <- read_excel(paste0(datadir2,'/Supplemental_Tables_v3.xlsx'), skip=1,
                         sheet = "Table S3", col_types = c("text", 
                                                           "text", "text", "text", "text", "text", 
                                                           "text", "text", "text", "text", "numeric", 
                                                           "numeric", "numeric", "numeric", 
                                                           "numeric", "numeric", "numeric", 
                                                           "numeric", "numeric", "numeric", 
                                                           "numeric", "numeric", "numeric", 
                                                           "numeric", "numeric", "numeric", 
                                                           "numeric", "numeric"))

df2 <- data.frame(df_ESI_POS$`Transformed Response Factor (RF(1/6))`,df_ESI_POS$`log10 Predicted IE`)
colnames(df2) <- c("RF(1/6)","log10Predicted_IE")
df2 <- na.omit(df2)
RF <- df2$`RF(1/6)`
IE <- df2$log10Predicted_IE

jagsdata_s2 <- with(df2, list(RF = RF, IE = IE, N = length(IE)))

lm2_jags <- function(){
  # Likelihood:
  for (i in 1:N){
    RF[i] ~ dnorm(mu[i], tau) # tau is precision (1 / variance)
    mu[i] <- alpha + beta * IE[i]
  }
  # Priors:
  alpha ~ dnorm(0, 0.01) # intercept
  beta ~ dnorm(0, 0.01) # slope
  sigma ~ dunif(0, 100) # standard deviation
  tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
}

init_values <- function(){
  list(alpha = rnorm(1), beta = rnorm(1), sigma = runif(1))
}

params <- c("alpha", "beta", "sigma")

fit_lm2 <- jags(data = jagsdata_s2, 
                inits = init_values, 
                parameters.to.save = params, 
                model.file = lm2_jags,
                n.chains = 3, 
                n.iter = 12000, 
                n.burnin = 2000, 
                n.thin = 10, 
                DIC = F)

fit_lm2

traceplot(fit_lm2, mfrow = c(2, 2), ask = F)

plot(fit_lm2)

lm2_mcmc <- as.mcmc(fit_lm2)
plot(lm2_mcmc)

nvalues <- 100
IE_new <- seq(min(IE), max(IE), length.out = nvalues)

lm2_mcmc_combi <- as.mcmc(rbind(lm2_mcmc[[1]], lm2_mcmc[[2]], lm2_mcmc[[3]]))

pred_mean_mean <- mean(lm2_mcmc_combi[, "alpha"]) + IE_new * mean(lm2_mcmc_combi[, "beta"])

pred_mean_dist <- matrix(NA, nrow = nrow(lm2_mcmc_combi), ncol = nvalues)
for (i in 1:nrow(pred_mean_dist)){
  pred_mean_dist[i,] <- lm2_mcmc_combi[i,"alpha"] + IE_new * lm2_mcmc_combi[i,"beta"]
}
credible_lower <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.005)
credible_upper <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.995)

lm2_mcmc_combi_rep <- do.call(rbind, rep(list(lm2_mcmc_combi), 50)) # replication

# Draw random values for all parameter combinations (rows) and body length values (columns):
pred_data_dist <- matrix(NA, nrow = nrow(lm2_mcmc_combi_rep), ncol = nvalues)
for (i in 1:nrow(pred_data_dist)){
  pred_data_dist[i,] <- lm2_mcmc_combi_rep[i,"alpha"] + IE_new * lm2_mcmc_combi_rep[i,"beta"] +
    rnorm(nvalues, mean = 0, sd = lm2_mcmc_combi_rep[i, "sigma"])
}

# Calculate quantiles:
uncertain_lower <- apply(pred_data_dist, MARGIN = 2, quantile, prob = 0.005)
uncertain_upper <- apply(pred_data_dist, MARGIN = 2, quantile, prob = 0.995)
ggplot()+
  geom_histogram(aes(x=pred_data_dist))
par(mar = c(0, 0, 0, 0))
plot(RF ~ IE, data = df2, 
     ylab = "Response Factor^(1/6)", 
     xlab = "log10(Predicted IE)",
     main = "ESI+")
lines(IE_new, pred_mean_mean)
lines(IE_new, credible_lower, lty = 2)
lines(IE_new, credible_upper, lty = 2)
lines(IE_new, uncertain_lower, lty = 2, col = "red")
lines(IE_new, uncertain_upper, lty = 2, col = "red")
# write.xlsx(cbind(pred_mean_mean,credible_lower,credible_upper,
#                  uncertain_lower,uncertain_upper,IE_new),
#            file='CNL_Bayes_PI_RF6thvsIE.xlsx')


#######################
#Linear mixed model
library(readr)
BritaGC_Inputs_for_CNL_BayesianLR <- read_csv("//aa.ad.epa.gov/ORD/RTP/USERS/A-D/CLowe/Net MyDocuments/R Data/Bayes_SQ/BritaGC_Inputs_for_CNL_BayesianLR.csv")


df2 <- data.frame(log10(BritaGC_Inputs_for_CNL_BayesianLR$RT),(BritaGC_Inputs_for_CNL_BayesianLR$RF)^(1/3),BritaGC_Inputs_for_CNL_BayesianLR$conc_level)
colnames(df2) <- c("RT","RF","CL")
df2 <- na.omit(df2)
RF <- df2$RF
RT <- df2$RT
CL <- df2$CL

set.seed(42)

mu <- c(10,10,10) # True means of normal distributions
sigma <- 5 # True standard deviation of normal distributions

# b_mass <- rnorm(samplesize, mean = mu, sd = sigma) # Body mass (response variable)

# df2 <- data.frame(b_length = b_length, b_mass = b_mass, site = sites)
# head(df2)

plot(RT ~ RF, col = CL, data = df2)

NCLs <- length(levels(as.factor(df2$CL)))
jagsdata_s3 <- with(df2, list(RT = RT, RF = RF, CL = CL,
                                  N = length(RT), NCLs = NCLs))

lm3_jags <- function(){
  # Likelihood:
  for (i in 1:N){
    RT[i] ~ dnorm(mu[i], tau) # tau is precision (1 / variance)
    mu[i] <- alpha + a[CL[i]] + beta * RF[i] # Random intercept for CL
  }
  # Priors:
  alpha ~ dnorm(0, 0.01) # intercept
  sigma_a ~ dunif(0, 100) # standard deviation of random effect (variance between CLs)
  tau_a <- 1 / (sigma_a * sigma_a) # convert to precision
  for (j in 1:NCLs){
    a[j] ~ dnorm(0, tau_a) # random intercept for each CL
  }
  beta ~ dnorm(0, 0.01) # slope
  sigma ~ dunif(0, 100) # standard deviation of fixed effect (variance within CLs)
  tau <- 1 / (sigma * sigma) # convert to precision
}

init_values <- function(){
  list(alpha = rnorm(1), sigma_a = runif(1), beta = rnorm(1), sigma = runif(1))
}

params <- c("alpha", "beta", "sigma", "sigma_a")

fit_lm3 <- jags(data = jagsdata_s3, 
                inits = init_values, 
                parameters.to.save = params, 
                model.file = lm3_jags,
                n.chains = 3, 
                n.iter = 20000, 
                n.burnin = 5000, 
                n.thin = 10, 
                DIC = F)

fit_lm3

traceplot(fit_lm3, mfrow = c(2, 2), ask = F)

plot(fit_lm3)

lm3_mcmc <- as.mcmc(fit_lm3)
par(mar=c(1,1,1,1))
plot(lm3_mcmc)

nvalues <- 100
RF_new <- seq(min(df2$RF), max(df2$RF), length.out = nvalues)

lm3_mcmc_combi <- as.mcmc(rbind(lm3_mcmc[[1]], lm3_mcmc[[2]], lm3_mcmc[[3]]))

pred_mean_mean <- mean(lm3_mcmc_combi[, "alpha"]) + RF_new * mean(lm3_mcmc_combi[, "beta"])

pred_mean_dist <- matrix(NA, nrow = nrow(lm3_mcmc_combi), ncol = nvalues)
for (i in 1:nrow(pred_mean_dist)){
  pred_mean_dist[i,] <- lm3_mcmc_combi[i,"alpha"] + RF_new * lm3_mcmc_combi[i,"beta"]
}
credible_lower <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.005)
credible_upper <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.995)

lm3_mcmc_combi_rep <- do.call(rbind, rep(list(lm3_mcmc_combi), 50)) # replication

# Draw random values for all parameter combinations (rows) and body length values (columns):
pred_data_dist <- matrix(NA, nrow = nrow(lm3_mcmc_combi_rep), ncol = nvalues)
for (i in 1:nrow(pred_data_dist)){
  pred_data_dist[i,] <- lm3_mcmc_combi_rep[i,"alpha"] + RF_new * lm3_mcmc_combi_rep[i,"beta"] +
    rnorm(nvalues, mean = 0, sd = lm3_mcmc_combi_rep[i, "sigma"])
}

# Calculate quantiles:
uncertain_lower <- apply(pred_data_dist, MARGIN = 2, quantile, prob = 0.005)
uncertain_upper <- apply(pred_data_dist, MARGIN = 2, quantile, prob = 0.995)
plot.new()
plot(RT ~ RF, col = CL, data = df2,xlab = "Response Factor^(1/3)",ylab = "log10(Retention Time)",ylim = c(-3,5))
lines(RF_new, pred_mean_mean)
lines(RF_new, credible_lower, lty = 2)
lines(RF_new, credible_upper, lty = 2)
lines(RF_new, uncertain_lower, lty = 2, col = "red")
lines(RF_new, uncertain_upper, lty = 2, col = "red")
