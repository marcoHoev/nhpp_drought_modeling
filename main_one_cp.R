library("R2jags")
library("readxl")
library("fitdistrplus")
library("tidyverse")

set.seed(42)

setwd("/Users/marco/dev/stochastic_processes/")
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

# Plottinng the count SPI graphs
plot(stepfun(1:(length(SPI_counted$`1 mês`)-1), SPI_counted$`1 mês`), cex.points = 0.1, lwd=0, main="Count SPI for 1 month")
plot(stepfun(1:(length(SPI_counted$`3 meses`)-1), SPI_counted$`3 meses`), cex.points = 0.1, lwd=0, main="Count SPI for 3 months")
plot(stepfun(1:(length(SPI_counted$`6 meses`)-1), SPI_counted$`6 meses`), cex.points = 0.1, lwd=0, main="Count SPI for 6 months")
plot(stepfun(1:(length(SPI_counted$`12 meses`)-1), SPI_counted$`12 meses`), cex.points = 0.1, lwd=0, main="Count SPI for 12 months")

# Parameters to be estimated
plp.mod.params <- c("alpha1", "sigma1", "alpha2", "sigma2", "tau1")

# Provide the data (for `1 mês` series)
y <- SPI_counted$`1 mês`
plp.mod.data <- list("y", "N")

# Define model
plp.mod <- function() {
  # Likelihood
  # for (i in 1:N) {
  #   if (i <= tau1) {
  #     y[i] ~ dpois((i/sigma1)**alpha1)
  #   } else { # mean = m1(tau) + m2(t) - m2(tau)
  #     y[i] ~ dpois((tau/sigma1)**alpha1 + (i/sigma2)**alpha2 - (tau/sigma2)**alpha2)
  #   }
  # }
  for (i in 1:floor(tau1)) {
    y[i] ~ dpois((i/sigma1)**alpha1)
  }
  for (i in ceiling(tau1):N) {
    y[i] ~ dpois((tau/sigma1)**alpha1 + (i/sigma2)**alpha2 - (tau/sigma2)**alpha2)
  }
  
  # Prior
  alpha1 ~ dunif(1e-5, 100)
  sigma1 ~ dunif(1e-5, 100)
  alpha2 ~ dunif(1e-5, 100)
  sigma2 ~ dunif(1e-5, 100)
  tau1 ~ dunif(1, N)
}

# Fit model
plp.mod.fit <- jags(data = plp.mod.data, 
                    parameters.to.save = plp.mod.params,
                    #inits = plp.mod.inits,
                    n.chains = 1, n.iter = 200000,
                    n.burnin = 10000, model.file = plp.mod)
plot(plp.mod.fit)
print(plp.mod.fit)

