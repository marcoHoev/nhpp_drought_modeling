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
#dir <- "/Users/marcelbraasch/RProjects/stochastic_processes/"
dir <- "/Users/marco/dev/stochastic_processes/"
data <- initalize(dir)

# Set which series to look at
names <- colnames(data)
N <- nrow(data)
current_name <- names[5]

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

############################# 0 change point plp ###############################

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

params <- estimate_model_with_no_changepoints(single_counts, current_name, burnin = 10000, iterations = 15000)

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
  for (i in 1:n_simulations) {
    S <- simulate_process_with_no_changepoints(params)
    lines(stepfun(1:(length(S)-1), S), cex.points = 0.01, lwd=0, col = "#0000FF")
  }
  lines(stepfun(1:(length(cum)-1), cum), cex.points = 0.1, lwd=0, col = "#000000")
  alpha1 <- params["alpha1"]
  sigma1 <- params["sigma1"]
  mean_to_end <- function(t) { (t/sigma1)**alpha1 }
  m <- c()
  for (i in 1:length(current_cumulative)) {
    m <- c(m, mean_to_end(i))
  }
  lines(stepfun(1:(length(m)-1),m), cex.points = 0.1, lwd=0, col = "#FF0000")
}

plot_data_and_simulation(current_cumulative, 10)

############################# 0 change point quadratic ###############################

estimate_quadratic_model_with_no_changepoints <- function(counts, name, burnin, iterations) {
  y <- counts[[name]]
  plp.mod.params <- c("a", "b") 
  plp.mod.data <- list("y", "N") 
  plp.mod <- function() {
    y[1] ~ dpois(m[1])
    m[1] <- a+b
    for (i in 2:N) {
      y[i] ~ dpois(m[i]-m[i-1])
      m[i] <- a*i+b*i**2
    }
    a ~ dunif(1e-5, 100)
    b ~ dunif(1e-5, 100)
  }
  plp.mod.fit <- jags(data = plp.mod.data, 
                      parameters.to.save = plp.mod.params,
                      n.chains = 3, n.iter = iterations,
                      n.burnin = burnin, model.file = plp.mod)
  print(plp.mod.fit)
  sigma1 <- plp.mod.fit$BUGSoutput[11][["mean"]][["a"]]
  alpha1 <- plp.mod.fit$BUGSoutput[11][["mean"]][["b"]]
  c("a" = a, "b" = b)
}

params <- estimate_quadratic_model_with_no_changepoints(single_counts, current_name, burnin = 10000, iterations = 15000)

get_params_0cp_without_recomputing <- function(name) {
  if (name=="1 mês") {
    return <- c("a" = 1.0526, "b" = 7.5730)
  } else if (name=="3 meses") {
    return <- c("a" = 1.0042, "b" = 6.7560)
  } else if (name=="6 meses") {
    return <- c("a" = 1.0498, "b" = 8.8409)
  } else if (name=="12 meses") {
    return <- c("a" = 0.124, "b" = 0.000)
  }
  return
}

params <- get_params_0cp_without_recomputing(current_name)

plot_data <- function(cum) {
  plot(0,0,xlim = c(0,N),ylim = c(0,max(cum)), type = "n",
       main = "Original data (black) with simulated processes (red)")
  lines(stepfun(1:(length(cum)-1), cum), cex.points = 0.1, lwd=0, col = "#000000")
  alpha1 <- params["a"]
  sigma1 <- params["b"]
  mean_to_end <- function(t) { (t/sigma1)**alpha1 }
  m <- c()
  for (i in 1:length(current_cumulative)) {
    m <- c(m, mean_to_end(i))
  }
  lines(stepfun(1:(length(m)-1),m), cex.points = 0.1, lwd=0, col = "#FF0000")
}

plot_data(current_cumulative)


############################# 1 changepoint plp #################################

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

params_1cp <- get_params_1cp_without_recomputing(current_name)

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

plot_data_and_mean_1_cp(current_cumulative, params_1cp)

############################# 2 change points plp ################################

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

params_2cp <- get_params_2cp_without_recomputing(current_name)

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

plot_data_and_mean_2_cp(current_cumulative, params_2cp)

############################# 3 change points plp ################################

estimate_model_with_three_changepoints <- function(counts, name, burnin, iterations) {
  y <- counts[[name]]
  plp.mod.params <- c("alpha1", "sigma1", "alpha2", "sigma2", "tau1", "alpha3", "sigma3", "tau2", "alpha4", "sigma4", "tau3")
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
                            ifelse(i<=tau3,
                                    (tau1/sigma1)**alpha1
                                    + (i/sigma3)**alpha3
                                    - (tau2/sigma3)**alpha3
                                    + (tau2/sigma2)**alpha2
                                    - (tau1/sigma2)**alpha2,
                                     (tau1/sigma1)**alpha1
                                     + (i/sigma4)**alpha4
                                     - (tau2/sigma3)**alpha3
                                     + (tau2/sigma2)**alpha2
                                     - (tau1/sigma2)**alpha2
                                     + (tau3/sigma3)**alpha3
                                     - (tau3/sigma4)**alpha4
                                    ))
      )
    }
    alpha1 ~ dunif(1e-5, 100)
    sigma1 ~ dunif(1e-5, 100)
    alpha2 ~ dunif(1e-5, 100)
    sigma2 ~ dunif(1e-5, 100)
    alpha3 ~ dunif(1e-5, 100)
    sigma3 ~ dunif(1e-5, 100)
    alpha4 ~ dunif(1e-5, 100)
    sigma4 ~ dunif(1e-5, 100)
    tau1 ~ dunif(0,400)
    tau2 ~ dunif(400,600)
    tau3 ~ dunif(600,N)
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
  sigma4 <- plp.mod.fit$BUGSoutput[11][["mean"]][["sigma4"]]
  alpha4 <- plp.mod.fit$BUGSoutput[11][["mean"]][["alpha4"]]
  tau1 <- plp.mod.fit$BUGSoutput[11][["mean"]][["tau1"]]
  tau2 <- plp.mod.fit$BUGSoutput[11][["mean"]][["tau2"]]
  tau3 <- plp.mod.fit$BUGSoutput[11][["mean"]][["tau3"]]
  c("alpha1" = alpha1, "sigma1" = sigma1, "alpha2" = alpha2,
    "sigma2" = sigma2, "sigma3" = sigma3, "alpha3" = alpha3,
    "alpha4" = alpha4, "sigma4" = sigma4, 
    "tau1" = tau1, "tau2" = tau2, "tau3" = tau3)
}

params <- estimate_model_with_three_changepoints(single_counts, current_name, 5000, 15000)

plot_data_and_mean_3_cp <- function(current_cumulative, params) {
  alpha1 <- params["alpha1"]
  alpha2 <- params["alpha2"]
  alpha3 <- params["alpha3"]
  alpha4 <- params["alpha4"]
  sigma1 <- params["sigma1"]
  sigma2 <- params["sigma2"]
  sigma3 <- params["sigma3"]
  sigma4 <- params["sigma4"]
  tau1 <- params["tau1"]
  tau2 <- params["tau2"]
  tau3 <- params["tau3"]
  mean_to_tau1 <- function(t) { (t/sigma1)**alpha1 }
  mean_to_tau2 <- function(t) { (tau1/sigma1)**alpha1+(t/sigma2)**alpha2-(tau1/sigma2)**alpha2 }
  mean_to_tau3 <- function(t) { (tau1/sigma1)**alpha1
                                 +(t   /sigma3)**alpha3
                                 -(tau2/sigma3)**alpha3
                                 +(tau2/sigma2)**alpha2
                                 -(tau1/sigma2)**alpha2
  }
  mean_to_end <- function(t) { (tau1/sigma1)**alpha1
                                 + (i/sigma4)**alpha4
                                 - (tau2/sigma3)**alpha3
                                 + (tau2/sigma2)**alpha2
                                 - (tau1/sigma2)**alpha2
                                 + (tau3/sigma3)**alpha3
                                 - (tau3/sigma4)**alpha4
  }
  m <- c()
  for (i in 1:length(current_cumulative)) {
    if (i <= tau1) {
      next_m <- mean_to_tau1(i)
    } else if (i <= tau2) {
      next_m <- mean_to_tau2(i)
    } else if (i <= tau3) {
      next_m <- mean_to_tau3(i)
    } else {
      next_m <- mean_to_end(i)
    }
    m <- c(m, next_m)
  }
  plot(0,0,xlim = c(0,N), ylim = c(0, max(current_cumulative)), type = "n", main = "")
  lines(stepfun(1:(length(current_cumulative)-1), current_cumulative),cex.points = 0.1, lwd=0, col = "#000000")
  lines(stepfun(1:(length(m)-1),m), cex.points = 0.1, lwd=0, col = "#FF0000")
}

plot_data_and_mean_3_cp(current_cumulative, params)

############################# 2 change points homogeneous #######################

estimate_linear_model_with_2cp <- function(counts, name, burnin, iterations) {
  y <- counts[[name]]
  plp.mod.params <- c("b1", "b2", "b3", "tau1", "tau2")
  plp.mod.data <- list("y", "N")
  plp.mod <- function() {
    y[1] ~ dpois(m[1])
    m[1] <- b1
    for (i in 2:N) {
      y[i] ~ dpois(m[i]-m[i-1])
      m[i] <- ifelse(i<=tau1,
                     b1*i,
                     ifelse(i<=tau2,
                            (b1*tau1)+(b2*i)-(b2*tau1),
                            (b1*tau1)+(b2*tau2)-(b2*tau1)+(b3*i)-(b3*tau2)
                     )
      )
    }
    b1 ~ dunif(1e-5, 100)
    b2 ~ dunif(1e-5, 100)
    b3 ~ dunif(1e-5, 100)
    tau1 ~ dunif(400,600)
    tau2 ~ dunif(600,N)
  }
  plp.mod.fit <- jags(data = plp.mod.data, 
                      parameters.to.save = plp.mod.params,
                      n.chains = 3, n.iter = iterations,
                      n.burnin = burnin, model.file = plp.mod)
  print(plp.mod.fit)
  
  b1 <- plp.mod.fit$BUGSoutput[11][["mean"]][["b1"]]
  b2 <- plp.mod.fit$BUGSoutput[11][["mean"]][["b2"]]
  b3 <- plp.mod.fit$BUGSoutput[11][["mean"]][["b3"]]
  tau1 <- plp.mod.fit$BUGSoutput[11][["mean"]][["tau1"]]
  tau2 <- plp.mod.fit$BUGSoutput[11][["mean"]][["tau2"]]
  c("b1" = b1, "b2" = b2, "b3" = b3,
    "tau1" = tau1, "tau2" = tau2)
}

params <- estimate_linear_model_with_2cp(single_counts, current_name, 10000, 15000)


get_params_linear_without_recomputing <- function(name) {
  if (name=="1 mês") {
    return <- c("b1" = 0.1618, "b2" = 0.5388, "b3" = 0.0800,
                "tau1" = 554.378, "tau2" = 629.485)
  } else if (name=="3 meses") {
    return <- c("b1" = 0.13534647, "b2" = 0.62597876, "b3" = 0.09277784,
                "tau1" = 559.881, "tau2" = 620.929)
  } else if (name=="6 meses") {
    return <- c("b1" = 0.1074, "b2" = 0.8814, "b3" = 0.0728,
                "tau1" = 562.561, "tau2" = 616.469)
  } else if (name=="12 meses") {
    return <- c("b1" = 0.08763865, "b2" = 1.00477610, "b3" = 0.03446260,
                "tau1" = 562.785, "tau2" = 622.393)
  }
  return
}

params <- get_params_linear_without_recomputing(current_name)

plot_data_and_mean_2_cp_linear <- function(current_cumulative, params) {
  b1 <- params["b1"]
  b2 <- params["b2"]
  b3 <- params["b3"]
  tau1 <- params["tau1"]
  tau2 <- params["tau2"]
  mean_to_tau1 <- function(t) { b1*t }
  mean_to_tau2 <- function(t) { ((b1*tau1)+(b2*t)-(b2*tau1)) }
  mean_to_end <- function(t) { ((b1*tau1)+(b2*tau2)-(b2*tau1)+(b3*t)-(b3*tau2))}
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

plot_data_and_mean_2_cp_linear(current_cumulative, params)

############################# 2 change points quadratic #######################

estimate_quad_model_with_2cp <- function(counts, name, burnin, iterations) {
  y <- counts[[name]]
  plp.mod.params <- c("b1", "b2", "b3", "c1", "c2", "c3", "tau1", "tau2")
  plp.mod.data <- list("y", "N")
  plp.mod <- function() {
    y[1] ~ dpois(m[1])
    m[1] <- b1 + c1
    for (i in 2:N) {
      y[i] ~ dpois(m[i]-m[i-1])
      m[i] <- ifelse(i<=tau1,
                     b1*i + c1*i^2,
                     ifelse(i<=tau2,
                            (b1*tau1 + c1*tau1^2) + (b2*i + c2*i^2) - (b2*tau1 + c2*tau1^2),
                            (b1*tau1 + c1*tau1^2) + (b2*tau2 + c2*tau2^2) - (b2*tau1 + c2*tau1^2) + (b3*i + c3*i^2) - (b3*tau2 + c3*tau2^2)
                     )
      )
    }
    b1 ~ dunif(1e-5, 1)
    b2 ~ dunif(1e-5, 1)
    b3 ~ dunif(1e-5, 1)
    c1 ~ dunif(-1e-1, 1e-1)
    c2 ~ dunif(-1e-1, 1e-1)
    c3 ~ dunif(-1e-1, 1e-1)
    tau1 ~ dunif(400,600)
    tau2 ~ dunif(600,N)
  }
  plp.mod.fit <- jags(data = plp.mod.data, 
                      parameters.to.save = plp.mod.params,
                      n.chains = 3, n.iter = iterations,
                      n.burnin = burnin, model.file = plp.mod)
  print(plp.mod.fit)
  
  b1 <- plp.mod.fit$BUGSoutput[11][["mean"]][["b1"]]
  b2 <- plp.mod.fit$BUGSoutput[11][["mean"]][["b2"]]
  b3 <- plp.mod.fit$BUGSoutput[11][["mean"]][["b3"]]
  c1 <- plp.mod.fit$BUGSoutput[11][["mean"]][["c1"]]
  c2 <- plp.mod.fit$BUGSoutput[11][["mean"]][["c2"]]
  c3 <- plp.mod.fit$BUGSoutput[11][["mean"]][["c3"]]
  tau1 <- plp.mod.fit$BUGSoutput[11][["mean"]][["tau1"]]
  tau2 <- plp.mod.fit$BUGSoutput[11][["mean"]][["tau2"]]
  c("b1" = b1, "b2" = b2, "b3" = b3,
    "c1" = c1, "c2" = c2, "c3" = c3,
    "tau1" = tau1, "tau2" = tau2)
}

params <- estimate_quad_model_with_2cp(single_counts, current_name, 10000, 15000)

get_params_quadratic_without_recomputing <- function(name) {
  if (name=="1 mês") {
    return <- c("b1" = 0.1785, "b2" = 0.5018, "b3" = 0.6402,
                "c1" = -2.837944e-05, "c2" = 8.174496e-05, "c3" = -3.895594e-04,
                "tau1" = 554.576, "tau2" = 615.436)
  } else if (name=="3 meses") {
    return <- c("b1" = 1.998658e-01, "b2" = 4.809733e-01, "b3" = 3.866826e-01,
                "c1" = -1.131177e-04, "c2" = 1.217071e-04, "c3" = -2.057753e-04,
                "tau1" = 557.997, "tau2" = 619.481)
  } else if (name=="6 meses") {
    return <- c("b1" = 2.001705e-01, "b2" = 5.027049e-01, "b3" = 2.053414e-01,
                "c1" = -1.618228e-04, "c2" = 3.280839e-04, "c3" = -9.065157e-05,
                "tau1" = 562.372, "tau2" = 616.379)
  } else if (name=="12 meses") {
    return <- c("b1" = 1.450162e-01, "b2" = 4.910334e-01, "b3" = 8.017581e-02,
                "c1" = -9.954617e-05, "c2" = 4.314868e-04, "c3" = -2.861638e-05,
                "tau1" = 562.623, "tau2" = 622.345)
  }
  return
}

params <- get_params_quadratic_without_recomputing(current_name)

plot_data_and_mean_2_cp_quad <- function(current_cumulative, params) {
  b1 <- params["b1"]
  b2 <- params["b2"]
  b3 <- params["b3"]
  c1 <- params["c1"]
  c2 <- params["c2"]
  c3 <- params["c3"]
  tau1 <- params["tau1"]
  tau2 <- params["tau2"]
  mean_to_tau1 <- function(t) { b1*t + c1*t^2 }
  mean_to_tau2 <- function(t) { ((b1*tau1 + c1*tau1^2)+(b2*t + c2*t^2)-(b2*tau1 + c2*tau1^2)) }
  mean_to_end <- function(t)  { ((b1*tau1 + c1*tau1^2)+(b2*tau2 + c2*tau2^2)-(b2*tau1+c2*tau1^2)+(b3*t+c3*t^2)-(b3*tau2+c3*tau2^2))}
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

plot_data_and_mean_2_cp_quad(current_cumulative, params)

############################# 2 change points root #######################

estimate_root_model_with_2cp <- function(counts, name, burnin, iterations) {
  y <- counts[[name]]
  plp.mod.params <- c("b1", "b2", "b3", "c1", "c2", "c3", "tau1", "tau2")
  plp.mod.data <- list("y", "N")
  plp.mod <- function() {
    y[1] ~ dpois(m[1])
    m[1] <- b1 + c1
    for (i in 2:N) {
      y[i] ~ dpois(m[i]-m[i-1])
      m[i] <- ifelse(i<=tau1,
                     b1*i + c1*i^(1/2),
                     ifelse(i<=tau2,
                            (b1*tau1 + c1*tau1^(1/2)) + (b2*i + c2*i^(1/2)) - (b2*tau1 + c2*tau1^(1/2)),
                            (b1*tau1 + c1*tau1^(1/2)) + (b2*tau2 + c2*tau2^(1/2)) - (b2*tau1 + c2*tau1^(1/2)) + (b3*i + c3*i^(1/2)) - (b3*tau2 + c3*tau2^(1/2))
                     )
      )
    }
    b1 ~ dunif(1e-5, 1)
    b2 ~ dunif(1e-5, 1)
    b3 ~ dunif(1e-5, 1)
    c1 ~ dunif(-1e-1, 1e-1)
    c2 ~ dunif(-1e-1, 1e-1)
    c3 ~ dunif(-1e-1, 1e-1)
    tau1 ~ dunif(400,600)
    tau2 ~ dunif(600,N)
  }
  plp.mod.fit <- jags(data = plp.mod.data, 
                      parameters.to.save = plp.mod.params,
                      n.chains = 3, n.iter = iterations,
                      n.burnin = burnin, model.file = plp.mod)
  print(plp.mod.fit)
  
  b1 <- plp.mod.fit$BUGSoutput[11][["mean"]][["b1"]]
  b2 <- plp.mod.fit$BUGSoutput[11][["mean"]][["b2"]]
  b3 <- plp.mod.fit$BUGSoutput[11][["mean"]][["b3"]]
  c1 <- plp.mod.fit$BUGSoutput[11][["mean"]][["c1"]]
  c2 <- plp.mod.fit$BUGSoutput[11][["mean"]][["c2"]]
  c3 <- plp.mod.fit$BUGSoutput[11][["mean"]][["c3"]]
  tau1 <- plp.mod.fit$BUGSoutput[11][["mean"]][["tau1"]]
  tau2 <- plp.mod.fit$BUGSoutput[11][["mean"]][["tau2"]]
  c("b1" = b1, "b2" = b2, "b3" = b3,
    "c1" = c1, "c2" = c2, "c3" = c3,
    "tau1" = tau1, "tau2" = tau2)
}

params <- estimate_root_model_with_2cp(single_counts, current_name, 10000, 15000)

plot_data_and_mean_2_cp_root <- function(current_cumulative, params) {
  b1 <- params["b1"]
  b2 <- params["b2"]
  b3 <- params["b3"]
  c1 <- params["c1"]
  c2 <- params["c2"]
  c3 <- params["c3"]
  tau1 <- params["tau1"]
  tau2 <- params["tau2"]
  mean_to_tau1 <- function(t) { b1*t + c1*t^(1/2) }
  mean_to_tau2 <- function(t) { ((b1*tau1 + c1*tau1^(1/2))+(b2*t + c2*t^(1/2))-(b2*tau1 + c2*tau1^(1/2))) }
  mean_to_end <- function(t)  { ((b1*tau1 + c1*tau1^(1/2))+(b2*tau2 + c2*tau2^(1/2))-(b2*tau1+c2*tau1^(1/2))+(b3*t+c3*t^(1/2))-(b3*tau2+c3*tau2^(1/2)))}
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

plot_data_and_mean_2_cp_root(current_cumulative, params)

############################# Plotting for report #######################


plot_all <- function(number) {
  current_name <- names[[number]]
  
  # Data
  cum <- cumulative_counts[[current_name]]
  plot(0,0,xlim = c(0,N),ylim = c(0,max(cum)), type = "l",
       xlab="Month", ylab="Cumulative values for SPI-12 ≤ -1.0",)
  lines(stepfun(1:(length(cum)-1), cum), cex.points = 0.15, lwd=0,
         col = "#000000")
  
  # 0CP
  alpha1 <- params["alpha1"]
  sigma1 <- params["sigma1"]
  mean_to_end <- function(t) { (t/sigma1)**alpha1 }
  m <- c()
  for (i in 1:length(current_cumulative)) {
    m <- c(m, mean_to_end(i))
  }
  lines(stepfun(1:(length(m)-1),m), cex.points = 0.1, lwd=0, col = "#16D881")
  
  # 1CP
  alpha1 <- params_1cp["alpha1"]
  alpha2 <- params_1cp["alpha2"]
  sigma1 <- params_1cp["sigma1"]
  sigma2 <- params_1cp["sigma2"]
  tau1 <- params_1cp["tau1"]
  mean_to_tau <- function(t) { (t/sigma1)**alpha1 }
  mean_to_end <- function(t) { ((tau1/sigma1)**alpha1+(t/sigma2)**alpha2-(tau1/sigma2)**alpha2) }
  mean_1_cp <- c()
  for (i in 1:length(current_cumulative)) {
    mean_1_cp <- c(mean_1_cp, ifelse(i <= tau1, mean_to_tau(i), mean_to_end(i)))
  }
  lines(stepfun(1:(length(mean_1_cp)-1),mean_1_cp), cex.points = 0.1, lwd=0, col = "#EF054A")
  
  # 2CP
  alpha1 <- params_2cp["alpha1"]
  alpha2 <- params_2cp["alpha2"]
  alpha3 <- params_2cp["alpha3"]
  sigma1 <- params_2cp["sigma1"]
  sigma2 <- params_2cp["sigma2"]
  sigma3 <- params_2cp["sigma3"]
  tau1 <- params_2cp["tau1"]
  tau2 <- params_2cp["tau2"]
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
  lines(stepfun(1:(length(m)-1),m), cex.points = 0.1, lwd=0, col = "#3E71F8")
}
number <- 5
current_name <- names[number]
params_2cp <- get_params_2cp_without_recomputing(current_name)
params_1cp <- get_params_1cp_without_recomputing(current_name)
params <- get_params_0cp_without_recomputing(current_name)
plot_all(number)

############################# Plotting extensions #######################


plot_all_ext <- function(number) {
  current_name <- names[[number]]
  
  # Data
  cum <- cumulative_counts[[current_name]]
  plot(0,0,xlim = c(0,N),ylim = c(0,max(cum)), type = "l",
       xlab="Month", ylab="Cumulative values for SPI-12 ≤ -1.0",)
  lines(stepfun(1:(length(cum)-1), cum), cex.points = 0.15, lwd=0,
        col = "#000000")
  
  # 0CP
  alpha1 <- params["alpha1"]
  sigma1 <- params["sigma1"]
  mean_to_end <- function(t) { (t/sigma1)**alpha1 }
  m <- c()
  for (i in 1:length(current_cumulative)) {
    m <- c(m, mean_to_end(i))
  }
  lines(stepfun(1:(length(m)-1),m), cex.points = 0.1, lwd=0, col = "#16D881")
  
  # 1CP
  alpha1 <- params_1cp["alpha1"]
  alpha2 <- params_1cp["alpha2"]
  sigma1 <- params_1cp["sigma1"]
  sigma2 <- params_1cp["sigma2"]
  tau1 <- params_1cp["tau1"]
  mean_to_tau <- function(t) { (t/sigma1)**alpha1 }
  mean_to_end <- function(t) { ((tau1/sigma1)**alpha1+(t/sigma2)**alpha2-(tau1/sigma2)**alpha2) }
  mean_1_cp <- c()
  for (i in 1:length(current_cumulative)) {
    mean_1_cp <- c(mean_1_cp, ifelse(i <= tau1, mean_to_tau(i), mean_to_end(i)))
  }
  lines(stepfun(1:(length(mean_1_cp)-1),mean_1_cp), cex.points = 0.1, lwd=0, col = "#EF054A")
  
  # 2CP
  alpha1 <- params_2cp["alpha1"]
  alpha2 <- params_2cp["alpha2"]
  alpha3 <- params_2cp["alpha3"]
  sigma1 <- params_2cp["sigma1"]
  sigma2 <- params_2cp["sigma2"]
  sigma3 <- params_2cp["sigma3"]
  tau1 <- params_2cp["tau1"]
  tau2 <- params_2cp["tau2"]
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
  lines(stepfun(1:(length(m)-1),m), cex.points = 0.1, lwd=0, col = "#3E71F8")
}
number <- 5
current_name <- names[number]
params_2cp <- get_params_2cp_without_recomputing(current_name)
params_1cp <- get_params_1cp_without_recomputing(current_name)
params <- get_params_0cp_without_recomputing(current_name)
plot_all(number)



