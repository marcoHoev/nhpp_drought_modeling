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

# Parameters to be estimated
plp.mod.params <- c("alpha1", "sigma1")

# Provide the data (for `1 mês` series)

y <- SPI_counted$"1 mês"
plp.mod.data <- list("y", "N")

# Define model
plp.mod <- function() {
  # Likelihood
  for (i in 1:N) {
    y[i] ~ dpois((i/sigma1)**(alpha1))
  }
  # Prior
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

# Plot the original data + process
plot(0,0,xlim = c(0,N),ylim = c(0,max(y)), type = "n", main = "")
lines(stepfun(1:(length(S)-1), S), cex.points = 0.1, lwd=0, col = "#FF0000")
lines(stepfun(1:(length(y)-1), y), cex.points = 0.1, lwd=0, col = "#000000")

