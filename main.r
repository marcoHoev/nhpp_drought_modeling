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
#dir <- "/Users/marco/dev/stochastic_processes/"
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

get_params_0cp_without_recomputing <- function(name) {
  if (name=="1 mês") {
    return <- c("alpha1" = 1.0526, "sigma1" = 7.5730)
  } else if (name=="3 meses") {
    return <- c("alpha1" = 1.0042, "sigma1" = 6.7560)
  } else if (name=="6 meses") {
    return <- c("alpha1" = 1.0498, "sigma1" = 8.8409)
  } else if (name=="12 meses") {
    return <- c("alpha1" = 1.2474, "sigma1" = 18.3014)
  }
  return
}

params <- get_params_0cp_without_recomputing(current_name)

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

get_params_1cp_without_recomputing <- function(name) {
  if (name=="1 mês") {
    return <- c("alpha1" = 1.1984, "sigma1" = 11.6938,
                "alpha2" = 1.2936, "sigma2" = 56.4024, "tau1" = 652.52)
  } else if (name=="3 meses") {
    return <- c("alpha1" = 1.0803, "sigma1" = 8.5861,
                "alpha2" = 1.3638, "sigma2" = 47.7410, "tau1" = 565.89)
  } else if (name=="6 meses") {
    return <- c("alpha1" = 0.7811, "sigma1" = 3.5134,
                "alpha2" = 1.1760, "sigma2" = 13.9759, "tau1" = 558.49)
  } else if (name=="12 meses") {
    return <- c("alpha1" = 0.8579, "sigma1" = 7.0486,
                "alpha2" = 1.0061, "sigma2" = 6.1740, "tau1" = 559.93)
  }
  return
}

params <- get_params_1cp_without_recomputing(current_name)

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
    tau1 ~ dunif(400,600)
    tau2 ~ dunif(600,N)
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

params <- estimate_model_with_two_changepoints(single_counts, current_name, 5000, 15000)

get_params_2cp_without_recomputing <- function(name) {
  if (name=="1 mês") {
    return <- c("alpha1" = 0.9864, "sigma1" = 6.4222, "alpha2" = 2.0600,
                "sigma2" = 56.1728, "sigma3" = 55.2291, "alpha3" = 1.3428,
                "tau1" = 552.0760, "tau2" = 632.6092)
  } else if (name=="3 meses") {
    return <- c("alpha1" = 0.8502, "sigma1" = 3.9576, "alpha2" = 2.1417,
                "sigma2" = 57.1688, "sigma3" = 56.0449, "alpha3" = 1.4342,
                "tau1" = 558.8111, "tau2" = 620.6246)
  } else if (name=="6 meses") {
    return <- c("alpha1" = 0.7765, "sigma1" = 3.4736, "alpha2" = 2.2520,
                "sigma2" = 55.9966, "sigma3" = 57.6732, "alpha3" = 1.3676,
                "tau1" = 562.8702, "tau2" = 617.0276)
  } else if (name=="12 meses") {
    return <- c("alpha1" = 0.8607, "sigma1" = 7.2052, "alpha2" = 2.1979,
                "sigma2" = 49.5556, "sigma3" = 58.9024, "alpha3" = 1.0946,
                "tau1" = 563.0705, "tau2" = 622.9102)
  }
  return
}

params <- get_params_2cp_without_recomputing(current_name)

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
  mean_to_end <- function(t) { ( (tau1/sigma1)**alpha1
                                +(t   /sigma3)**alpha3
                                -(tau2/sigma3)**alpha3
                                +(tau2/sigma2)**alpha2
                                -(tau1/sigma2)**alpha2)
    }
  m <- c()
  for (i in 1:length(current_cumulative)) {
    if (i <= tau1) {
      next_m <- mean_to_tau1(i)
    } else if (i <= tau2) {
      next_m <- mean_to_tau2(i)
    } else {
      next_m <- mean_to_end(i)
    }
    m <- c(m, next_m)
  }
  plot(0,0,xlim = c(0,N), ylim = c(0, max(current_cumulative)), type = "n", main = "")
  lines(stepfun(1:(length(current_cumulative)-1), current_cumulative),cex.points = 0.1, lwd=0, col = "#000000")
  lines(stepfun(1:(length(m)-1),m), cex.points = 0.1, lwd=0, col = "#FF0000")
}

plot_data_and_mean_2_cp(current_cumulative, params)
