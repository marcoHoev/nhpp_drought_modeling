# Mac users: Install the GNU Fortran (gfortran-4.2.3.dmg) library from the CRAN toolsdirectory:http://cran.r-project.org/bin/macosx/tools.3.
# Install JAGS version 3.4.0 from Martyn Plummer’s repository:http://cran.r-project.org/bin/macosx/tools
# Install JAGS version 3.4.0 from Martyn Plummer’s repository:http://sourceforge.net/projects/mcmc-jags/files/JAGS/3.x/
# install.packages("R2jags",dependencies=TRUE,repos="http://cran.us.r-project.org")
# install.packages("runjags",dependencies=TRUE,repos="http://cran.us.r-project.org")
# install.packages("MCMCpack",dependencies=TRUE,repos="http://cran.us.r-project.org")
# install.packages("fitdistrplus", dependencies=TRUE,repos="http://cran.us.r-project.org")

library("R2jags")
library("readxl")
library("fitdistrplus")
library("tidyverse")

set.seed(42)

setwd("/Users/marcelbraasch/RProjects/stochastic_processes/")
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

# Define model
plp.mod <- function() {
  # Likelihood
  for (i in 1:N) {
    y[i] <- lambda[t[i]] * exp(mean1)
    lambda[i] <- (alpha1/sigma1) * ((t[i]/sigma1)**(alpha1-1))
  }
  # Mean until the point T
  mean1 <- (t[N]/sigma1)**alpha1

    # Prior
  alpha1 ~ dunif(0, 100)
  sigma1 ~ dunif(0, 100)
}

# Provide the data (for `1 mês` series)
t <- SPI_counted$`1 mês`
plp.mod.data <- list("t")
plp.mod.params <- c("alpha1", "sigma1")
plp.mod.inits <- function(){
  list("alpha1" = runif(1, min=0, max=100), "sigma1" = runif(1, min=0, max=100))
}

plp.mod.fit <- jags(data = plp.mod.data, 
                    inits = plp.mod.inits,
                    parameters.to.save = plp.mod.params,
                    n.chains = 3, n.iter = 9000,
                    n.burnin = 1000, model.file = plp.mod)
