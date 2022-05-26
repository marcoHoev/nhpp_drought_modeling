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

setwd("/Users/marcelbraasch/RProjects/stochastic_processes/")
data <- read_excel("data.xlsx")
data.1 <- data$`1 mês`
data.3 <- data$`3 meses`
data.6 <- data$`6 meses`
data.12 <- data$`12 meses`
data.24 <- data$`24 meses`

# Create the cumulative data
counter <- 0
cumulative <- c()
for (i in 1:length(data.1)) {
  if (data.1[i] <= -1) { counter <- counter + 1}
  cumulative <- c(cumulative, counter)
}
data.1.cumulative <- cumulative
