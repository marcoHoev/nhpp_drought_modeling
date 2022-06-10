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
# dir <- "/Users/marcelbraasch/RProjects/stochastic_processes/"
dir <- "/Users/marco/dev/stochastic_processes/"
data <- initalize(dir)

# Set which series to look at
names <- colnames(data)
N <- nrow(data)
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

single_counts <- create_single_counts(data)

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

cumulative_counts <- create_cumulative_counts(data)
current_cumulative <- cumulative_counts[[current_name]]

############################# No change point ###############################

estimate_model_with_no_changepoints <- function(counts, name, burnin, iterations) {
  y <- counts[[name]]
  plp.mod.params <- c("alpha1", "sigma1") 
  plp.mod.data <- list("y", "N") 
  plp.mod <- function() {
    y[1] ~ dpois(m[1])
    m[1] <- (1/sigma1)^alpha1
    for (i in 2:N) {
      y[i] ~ dpois(m[i]-m[i-1])
      m[i] <- (i/sigma1)^alpha1
    }
    alpha1 ~ dunif(1e-5, 100)
    sigma1 ~ dunif(1e-5, 100)
  }
  plp.mod.fit <- jags(data = plp.mod.data, 
                      parameters.to.save = plp.mod.params,
                      n.chains = 3, n.iter = iterations,
                      n.burnin = burnin, model.file = plp.mod)
  print(plp.mod.fit)
  sigma1 <- plp.mod.fit$BUGSoutput[11][["mean"]][["sigma1"]]
  alpha1 <- plp.mod.fit$BUGSoutput[11][["mean"]][["alpha1"]]
  c("alpha1" = alpha1, "sigma1" = sigma1)
}

params <- estimate_model_with_no_changepoints(single_counts, current_name, burnin = 20000, iterations = 25000)

# Simulate process by thinning procedure and plot it together with the original data
# https://stats.stackexchange.com/questions/369288/nonhomogeneous-poisson-process-simulation
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

# Plot the original data + n simulated process
plot_data_and_simulation <- function(cum, n_simulations) {
  plot(0,0,xlim = c(0,N),ylim = c(0,max(cum)), type = "n",
       main = "Original data (black) with simulated processes (red)")
  lines(stepfun(1:(length(cum)-1), cum), cex.points = 0.1, lwd=0, col = "#000000")
  for (i in 1:n_simulations) {
    S <- simulate_process_with_no_changepoints(params)
    lines(stepfun(1:(length(S)-1), S), cex.points = 0.01, lwd=0, col = "#FF0000")
  }
}

plot_data_and_simulation(current_cumulative, 10)

############################# 1 Changepoint #################################

estimate_model_with_one_changepoint <- function(counts, name, burnin, iterations) {
  y <- counts[[name]]
  plp.mod.params <- c("alpha1", "sigma1", "alpha2", "sigma2", "tau1")
  plp.mod.data <- list("y", "N")
  plp.mod <- function() {
    y[1] ~ dpois(m[1])
    m[1] <- (1/sigma1)^alpha1
    for (i in 2:N) {
      y[i] ~ dpois(m[i]-m[i-1])
      m[i] <- ifelse(i<=tau1,
                     (i/sigma1)**alpha1,
                     (tau1/sigma1)**alpha1
                     + (i/sigma2)**alpha2
                     - (tau1/sigma2)**alpha2
                     )
    }
    alpha1 ~ dunif(1e-5, 100)
    sigma1 ~ dunif(1e-5, 100)
    alpha2 ~ dunif(1e-5, 100)
    sigma2 ~ dunif(1e-5, 100)
    tau1 ~ dunif(0,N)
  }
  plp.mod.fit <- jags(data = plp.mod.data, 
                      parameters.to.save = plp.mod.params,
                      n.chains = 3, n.iter = iterations,
                      n.burnin = burnin, model.file = plp.mod)
  print(plp.mod.fit)
  sigma1 <- plp.mod.fit$BUGSoutput[11][["mean"]][["sigma1"]]
  alpha1 <- plp.mod.fit$BUGSoutput[11][["mean"]][["alpha1"]]
  sigma2 <- plp.mod.fit$BUGSoutput[11][["mean"]][["sigma2"]]
  alpha2 <- plp.mod.fit$BUGSoutput[11][["mean"]][["alpha2"]]
  tau1 <- plp.mod.fit$BUGSoutput[11][["mean"]][["tau1"]]
  c("alpha1" = alpha1, "sigma1" = sigma1, "alpha2" = alpha2,
    "sigma2" = sigma2, "tau1" = tau1)
}

params <- estimate_model_with_one_changepoint(single_counts, current_name, burnin = 10000, iterations = 15000)

plot_data_and_mean_1_cp <- function(current_cumulative, params) {
  alpha1 <- params["alpha1"]
  alpha2 <- params["alpha2"]
  sigma1 <- params["sigma1"]
  sigma2 <- params["sigma2"]
  tau1 <- params["tau1"]
  mean_to_tau <- function(t) { (t/sigma1)**alpha1 }
  mean_to_end <- function(t) { ((tau1/sigma1)**alpha1+(t/sigma2)**alpha2-(tau1/sigma2)**alpha2) }
  mean_1_cp <- c()
  for (i in 1:length(current_cumulative)) {
    mean_1_cp <- c(mean_1_cp, ifelse(i <= tau1, mean_to_tau(i), mean_to_end(i)))
  }
  plot(0,0,xlim = c(0,N), ylim = c(0, max(current_cumulative)), type = "n", main = "")
  lines(stepfun(1:(length(current_cumulative)-1), current_cumulative),cex.points = 0.1, lwd=0, col = "#000000")
  lines(stepfun(1:(length(mean_1_cp)-1),mean_1_cp), cex.points = 0.1, lwd=0, col = "#FF0000")
}

plot_data_and_mean_1_cp(current_cumulative, params)

############################# 2 change points ################################

estimate_model_with_two_changepoints <- function(counts, name, burnin, iterations) {
  y <- counts[[name]]
  plp.mod.params <- c("alpha1", "sigma1", "alpha2", "sigma2", "tau1", "alpha3", "sigma3", "tau2")
  plp.mod.data <- list("y", "N")
  plp.mod <- function() {
    y[1] ~ dpois(m[1])
    m[1] <- (1/sigma1)^alpha1
    for (i in 2:N) {
      y[i] ~ dpois(m[i]-m[i-1])
      m[i] <- ifelse(i<=tau1,
                     (i/sigma1)**alpha1,
                     ifelse(i<=tau2,
                            (tau1/sigma1)**alpha1
                            + (i/sigma2)**alpha2
                            - (tau1/sigma2)**alpha2,
                            (tau1/sigma1)**alpha1
                            + (i/sigma3)**alpha3
                            - (tau2/sigma3)**alpha3
                            + (tau2/sigma2)**alpha2
                            - (tau1/sigma2)**alpha2)
                     )
    }
    alpha1 ~ dunif(1e-5, 100)
    sigma1 ~ dunif(1e-5, 100)
    alpha2 ~ dunif(1e-5, 100)
    sigma2 ~ dunif(1e-5, 100)
    alpha3 ~ dunif(1e-5, 100)
    sigma3 ~ dunif(1e-5, 100)
    tau1 ~ dunif(0,N)
    tau2 ~ dunif(0,N)
  }
  plp.mod.fit <- jags(data = plp.mod.data, 
                      parameters.to.save = plp.mod.params,
                      n.chains = 3, n.iter = iterations,
                      n.burnin = burnin, model.file = plp.mod)
  print(plp.mod.fit)
  sigma1 <- plp.mod.fit$BUGSoutput[11][["mean"]][["sigma1"]]
  alpha1 <- plp.mod.fit$BUGSoutput[11][["mean"]][["alpha1"]]
  sigma2 <- plp.mod.fit$BUGSoutput[11][["mean"]][["sigma2"]]
  alpha2 <- plp.mod.fit$BUGSoutput[11][["mean"]][["alpha2"]]
  sigma3 <- plp.mod.fit$BUGSoutput[11][["mean"]][["sigma3"]]
  alpha3 <- plp.mod.fit$BUGSoutput[11][["mean"]][["alpha3"]]
  tau1 <- plp.mod.fit$BUGSoutput[11][["mean"]][["tau1"]]
  tau2 <- plp.mod.fit$BUGSoutput[11][["mean"]][["tau2"]]
  c("alpha1" = alpha1, "sigma1" = sigma1, "alpha2" = alpha2,
    "sigma2" = sigma2, "sigma3" = sigma3, "alpha3" = alpha3,
    "tau1" = tau1, "tau2" = tau2)
}

params <- estimate_model_with_two_changepoints(single_counts, current_name, 10, 20)

plot_data_and_mean_2_cp <- function(current_cumulative, params) {
  alpha1 <- params["alpha1"]
  alpha2 <- params["alpha2"]
  alpha3 <- params["alpha3"]
  sigma1 <- params["sigma1"]
  sigma2 <- params["sigma2"]
  sigma3 <- params["sigma3"]
  tau1 <- params["tau1"]
  tau2 <- params["tau2"]
  mean_to_tau1 <- function(t) { (t/sigma1)**alpha1 }
  mean_to_tau2 <- function(t) { ((tau1/sigma1)**alpha1+(t/sigma2)**alpha2-(tau1/sigma2)**alpha2) }
  mean_to_end <- function(t) { (  (tau1/sigma1)**alpha1
                                + (t/sigma3)**alpha3
                                + (tau2/sigma3)**alpha3
                                + (tau2/sigma2)**alpha2
                                + (tau1/sigma2)**alpha2) }
  mean_2_cp <- c()
  for (i in 1:length(current_cumulative)) {
    if (i <= tau1) {
      current_mean <- mean_to_tau1(i)
    } else if (i <= tau2) {
      current_mean <- mean_to_tau2(i)
    } else {
      current_mean <- mean_to_end(i)
    }
    mean_2_cp <- c(mean_2_cp, current_mean)
  }
  plot(0,0,xlim = c(0,N), ylim = c(0, max(current_cumulative)), type = "n", main = "")
  lines(stepfun(1:(length(current_cumulative)-1), current_cumulative),cex.points = 0.1, lwd=0, col = "#000000")
  lines(stepfun(1:(length(mean_2_cp)-1),mean_2_cp), cex.points = 0.1, lwd=0, col = "#FF0000")
}

plot_data_and_mean_2_cp(current_cumulative, params)