# Mac users: Install the GNU Fortran (gfortran-4.2.3.dmg) library from the CRAN toolsdirectory:http://cran.r-project.org/bin/macosx/tools.3.
# Install JAGS version 3.4.0 from Martyn Plummer’s repository:http://cran.r-project.org/bin/macosx/tools
# Install JAGS version 3.4.0 from Martyn Plummer’s repository:http://sourceforge.net/projects/mcmc-jags/files/JAGS/3.x/
# install.packages("R2jags",dependencies=TRUE,repos="http://cran.us.r-project.org")
# install.packages("runjags",dependencies=TRUE,repos="http://cran.us.r-project.org")
# install.packages("MCMCpack",dependencies=TRUE,repos="http://cran.us.r-project.org")
# install.packages("fitdistrplus", dependencies=TRUE,repos="http://cran.us.r-project.org")
# install.packages("poisson", dependencies=TRUE,repos="http://cran.us.r-project.org")

library("R2jags")
library("readxl")
library("tidyverse")
library("poisson")

set.seed(42)

setwd("/Users/marcelbraasch/RProjects/stochastic_processes/")
#setwd("/Users/marco/dev/stochastic_processes/")
data <- read_excel("data.xlsx")
data <- select(data, -1)  # Drop first

# Create the cumulative data for all five series
N <- nrow(data)
M <- ncol(data)
names <- colnames(data)
SPI_counted <- tibble(a = 1:N)
for (j in 1:M) {
  d <- data[[j]]
  name <- names[j]
  counter <- 0
  cumulative <- c()
  for (i in 1:N) {
    if (d[i] <= -1) { counter <- counter + 1}
    cumulative <- c(cumulative, counter)
  }
  SPI_counted$x <- cumulative
  names(SPI_counted)[ncol(SPI_counted)] <- names[j]  # Rename
}
y <- SPI_counted$"1 mês"

plp.mod.params <- c("alpha1", "sigma1") # Params to be estimated
plp.mod.data <- list("y", "N") # Data to be fed
plp.mod <- function() { # Define model
  for (i in 1:N) { y[i] ~ dpois((i/sigma1)**(alpha1)) } 
  alpha1 ~ dunif(1e-5, 100)
  sigma1 ~ dunif(1e-5, 100)
}

# Fit model
plp.mod.fit <- jags(data = plp.mod.data, 
                    parameters.to.save = plp.mod.params,
                    n.chains = 3, n.iter = 12000,
                    n.burnin = 10000, model.file = plp.mod)
plot(plp.mod.fit)
print(plp.mod.fit)

# Plot process by thinning procedure
# https://stats.stackexchange.com/questions/369288/nonhomogeneous-poisson-process-simulation
sim <- function(){
  sigma_hat <- plp.mod.fit$BUGSoutput[11]$mean$sigma1
  alpha_hat <- plp.mod.fit$BUGSoutput[11]$mean$alpha1
  rate <- function(t) {(t/sigma_hat)**(alpha_hat)}
  intensity <- function(t) {(alpha_hat/sigma_hat)*(t/sigma_hat)**(alpha_hat-1)}
  p <- function(t) {intensity(t)/alpha_hat}
  S <- c()
  t <- 0
  I <- 0
  for (t in 1:N) {
    U_1 <- runif(1)
    t <- t - (log(U_1)/alpha_hat)
    U_2 <- runif(1)
    if (U_2 < p(t)) {I <- I + 1}
    S <- c(S, I)
  }
  S
} 

# Plot the original data + process
plot_data_and_process <- function(y, S) {
  plot(0,0,xlim = c(0,N),ylim = c(0,max(y)), type = "n", main = "")
  lines(stepfun(1:(length(S)-1), S), cex.points = 0.1, lwd=0, col = "#FF0000")
  lines(stepfun(1:(length(y)-1), y), cex.points = 0.1, lwd=0, col = "#000000")
}


###############################################################################
############################# 2 Changepoints ##################################
###############################################################################

# Model for two changepoints
plp.mod <- function() {
  for (i in 1:N) {
    y[i] ~ dpois(ifelse(i<=tau1,
                        (i/sigma1)**alpha1,
                        (tau1/sigma1)**alpha1 + (i/sigma2)**alpha2 - (tau1/sigma2)**alpha2
                        )
                 )
  }
  alpha1 ~ dunif(1e-5, 100)
  sigma1 ~ dunif(1e-5, 100)
  alpha2 ~ dunif(1e-5, 100)
  sigma2 ~ dunif(1e-5, 100)
  tau1 ~ dunif(0,N)
}

plp.mod.params <- c("alpha1", "sigma1", "alpha2", "sigma2", "tau1") # Params to be estimated
plp.mod.data <- list("y", "N")#, "sigma", "alpha") # Data to be fed

# Fit model
plp.mod.fit <- jags(data = plp.mod.data, 
                    parameters.to.save = plp.mod.params,
                    #n.chains = 1, n.iter = 200,
                    #n.burnin = 100, model.file = plp.mod)
                    n.chains = 3, n.iter = 12000,
                    n.burnin = 10000, model.file = plp.mod)
plot(plp.mod.fit)
print(plp.mod.fit)

#simulate <- function(sigmas, alphas, taus) {
#  #sigma_hat <- plp.mod.fit$BUGSoutput[11]$mean$sigma1
#  #alpha_hat <- plp.mod.fit$BUGSoutput[11]$mean$alpha1
#  sigma_hat <- sigmas
#  alpha_hat <- alphas
#  tau_hat <- taus
#  rate <- function(t) {(t/sigma_hat)**(alpha_hat)}
#  intensity <- function(t) {(alpha_hat/sigma_hat)*(t/sigma_hat)**(alpha_hat-1)}
#  p <- function(t) {intensity(t)/alpha_hat}
#  S <- c()
#  t <- 0
#  I <- 0
#  for (t in 1:N) {
#    U_1 <- runif(1)
#    t <- t - (log(U_1)/alpha_hat)
#    U_2 <- runif(1)
#    if (U_2 < p(t)) {I <- I + 1}
#    S <- c(S, I)
#  }
#  S
#}
