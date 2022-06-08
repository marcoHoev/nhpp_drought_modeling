############################# Installations ##################################

# Mac users: Install the GNU Fortran (gfortran-4.2.3.dmg) library from the CRAN toolsdirectory:http://cran.r-project.org/bin/macosx/tools.3.
# Install JAGS version 3.4.0 from Martyn Plummer’s repository:http://cran.r-project.org/bin/macosx/tools
# Install JAGS version 3.4.0 from Martyn Plummer’s repository:http://sourceforge.net/projects/mcmc-jags/files/JAGS/3.x/
# install.packages("R2jags",dependencies=TRUE,repos="http://cran.us.r-project.org")
# install.packages("runjags",dependencies=TRUE,repos="http://cran.us.r-project.org")
# install.packages("MCMCpack",dependencies=TRUE,repos="http://cran.us.r-project.org")
# install.packages("fitdistrplus", dependencies=TRUE,repos="http://cran.us.r-project.org")
# install.packages("poisson", dependencies=TRUE,repos="http://cran.us.r-project.org")

############################# Initialize ####################################

initalize <- function(dir) {
  library("R2jags")
  library("readxl")
  library("tidyverse")
  library("poisson")
  set.seed(42)
  setwd(dir)
  data <- read_excel("./Data/data.xlsx")
  data <- select(data, -1)  # Drop first
}
dir <- "/Users/marcelbraasch/RProjects/stochastic_processes/"
# dir <- "/Users/marco/dev/stochastic_processes/"
data <- initalize(dir)

# Set which series to look at
names <- colnames(data)
current_name <- names[1]

############################# Create data ###################################

create_single_counts <- function(data) {
  N <- nrow(data)
  M <- ncol(data)
  names <- colnames(data)
  SPI_counted <- data.frame(a = 1:N)
  for (name in names) {
    y <- c()
    for (point in data[[name]]) {
      y <- c(y, ifelse(point <= -1, 1, 0))
    }
    SPI_counted[name] <- y
  }
  SPI_counted <- subset(SPI_counted, select = -1)
}

create_cumulative_counts <- function(data) {
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
  SPI_counted
}

single_counts <- create_single_counts(data)
cumulative_counts <- create_cumulative_counts(data)

############################# No change point ###############################

estimate_model_with_no_changepoints <- function(counts, name) {
  y <- counts[[name]]
  # Define parameters
  plp.mod.params <- c("alpha1", "sigma1") 
  plp.mod.data <- list("y", "N") 
  # Define model
  plp.mod <- function() {
    for (i in 1:N) {
      y[i] ~ dpois( (i/sigma1)^alpha1 - ((i-1)/sigma1)^alpha1)
    }
    alpha1 ~ dunif(1e-5, 100)
    sigma1 ~ dunif(1e-5, 100)
  }
  # Fit model
  plp.mod.fit <- jags(data = plp.mod.data, 
                      parameters.to.save = plp.mod.params,
                      n.chains = 3, n.iter = 15000,
                      n.burnin = 10000, model.file = plp.mod)
  print(plp.mod.fit)
  sigma <- plp.mod.fit$BUGSoutput[11]$mean$sigma1
  alpha <- plp.mod.fit$BUGSoutput[11]$mean$alpha1
  c("alpha1" = alpha, "sigma1" = sigma)
}

params <- estimate_model_with_no_changepoints(single_counts, current_name)

# Simulate process by thinning procedure and plot it together with the original data
# https://stats.stackexchange.com/questions/369288/nonhomogeneous-poisson-process-simulation
# simulate_process_with_no_changepoints <- function(alpha1, sigma1){
simulate_process_with_no_changepoints <- function(params) {
  alpha1 <- params["alpha1"]
  sigma1 <- params["sigma1"]
  rate <- function(t) {(t/sigma1)**(alpha_hat)}
  intensity <- function(t) {(alpha1/sigma1)*(t/sigma1)**(alpha1-1)}
  p <- function(t) {intensity(t)/alpha1}
  S <- c()
  t <- 0
  I <- 0
  for (t in 1:N) {
    U_1 <- runif(1)
    t <- t - (log(U_1)/alpha1)
    U_2 <- runif(1)
    if (U_2 < p(t)) {I <- I + 1}
    S <- c(S, I)
  }
  S
}

# Plot the original data + process
plot_data_and_simulation <- function(c, S) {
  plot(0,0,xlim = c(0,N),ylim = c(0,max(c)), type = "n", main = "")
  lines(stepfun(1:(length(S)-1), S), cex.points = 0.1, lwd=0, col = "#FF0000")
  lines(stepfun(1:(length(c)-1), c), cex.points = 0.1, lwd=0, col = "#000000")
}

S <- simulate_process_with_no_changepoints(params)
plot_data_and_simulation(cumulative_counts[[current_name]], S)

############################# 1 Changepoint #################################

# Model for one changepoint
plp.mod <- function() {
  for (i in 1:N) {
    y[i] ~ dpois(m[i])
    m[i] <- ifelse(i<=tau1,
                   (i/sigma1)**alpha1,
                   (tau1/sigma1)**alpha1
                   + (i/sigma2)**alpha2
                   - (tau1/sigma2)**alpha2)
    }
  alpha1 ~ dunif(1e-5, 100)
  sigma1 ~ dunif(1e-5, 100)
  alpha2 ~ dunif(1e-5, 100)
  sigma2 ~ dunif(1e-5, 100)
  tau1 ~ dunif(0,N)
}

plp.mod.params <- c("alpha1", "sigma1", "alpha2", "sigma2", "tau1")
plp.mod.data <- list("y", "N")
plp.mod.fit <- jags(data = plp.mod.data, 
                    parameters.to.save = plp.mod.params,
                    n.chains = 3, n.iter = 30000,
                    n.burnin = 10000, model.file = plp.mod)
print(plp.mod.fit)
