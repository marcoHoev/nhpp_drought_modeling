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

estimate_model_extended <- function(counts, name, burnin, iterations) {
  y <- counts[[name]]
  plp.mod.params <- c("a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8") 
  plp.mod.data <- list("y", "N") 
  plp.mod <- function() {
    y[1] ~ dpois(m[1])
    m[1] <- a1*1 + a2*1^2 + a3*1^3 + a4*1^4 + a5*1^5 + a6*1^6 + a7*1^7 + a8*1^8 
    for (i in 2:N) {
      y[i] ~ dpois(m[i]-m[i-1])
      m[i] <- a1*i + a2*i^2 + a3*i^3 + a4*i^4 + a5*i^5 + a6*i^6 + a7*i^7 + a8*i^8 
    }
    a1 ~ dunif(1e-10, 1)
    a2 ~ dunif(-1000, 1000)
    a3 ~ dunif(-1000, 1000)
    a4 ~ dunif(-1000, 1000)
    a5 ~ dunif(-1000, 1000)
    a6 ~ dunif(-1000, 1000)
    a7 ~ dunif(-1000, 1000)
    a8 ~ dunif(-1000, 1000)
  }
  plp.mod.fit <- jags(data = plp.mod.data, 
                      parameters.to.save = plp.mod.params,
                      n.chains = 3, n.iter = iterations,
                      n.burnin = burnin, model.file = plp.mod)
  print(plp.mod.fit)
  a0 <- plp.mod.fit$BUGSoutput[11][["mean"]][["a0"]]
  a1 <- plp.mod.fit$BUGSoutput[11][["mean"]][["a1"]]
  a2 <- plp.mod.fit$BUGSoutput[11][["mean"]][["a2"]]
  a3 <- plp.mod.fit$BUGSoutput[11][["mean"]][["a3"]]
  a4 <- plp.mod.fit$BUGSoutput[11][["mean"]][["a4"]]
  a5 <- plp.mod.fit$BUGSoutput[11][["mean"]][["a5"]]
  a6 <- plp.mod.fit$BUGSoutput[11][["mean"]][["a6"]]
  a7 <- plp.mod.fit$BUGSoutput[11][["mean"]][["a7"]]
  a8 <- plp.mod.fit$BUGSoutput[11][["mean"]][["a8"]]
  c("a1" = a1, "a2" = a2, "a3" = a3, "a4" = a4, "a5" = a5, "a6" = a6, "a7" = a7, "a8" = a8)
}

params <- estimate_model_extended(single_counts, current_name, burnin = 20000, iterations = 25000)

plot_data_and_mean_extended <- function(current_cumulative, params) {
  a1 <- params["a1"]
  a2 <- params["a2"]
  a3 <- params["a3"]
  a4 <- params["a4"]
  a5 <- params["a5"]
  a6 <- params["a6"]
  a7 <- params["a7"]
  a8 <- params["a8"]
  mean_to_end <- function(t) { a1*t + a2*t^2 + a3*t^3 + a4*t^4 + a5*t^5 + a6*t^6 + a7*t^7 + a8*t^8}
  mean_ext <- c()
  for (i in 1:length(current_cumulative)) {
    mean_ext <- c(mean_ext, mean_to_end(i))
  }
  plot(0,0,xlim = c(0,N), ylim = c(0, max(current_cumulative)), type = "n", main = "")
  lines(stepfun(1:(length(current_cumulative)-1), current_cumulative),cex.points = 0.1, lwd=0, col = "#000000")
  lines(stepfun(1:(length(mean_ext)-1),mean_ext), cex.points = 0.1, lwd=0, col = "#FF0000")
}

plot_data_and_mean_extended(current_cumulative, params)
print(params)
