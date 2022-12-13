library(BayesFMMM)
library(ggplot2)
library(gridExtra)
library(latex2exp)
library(pracma)

### Set working dir
setwd()

Y <- readRDS(system.file("test-data", "Sim_data.RDS", package = "BayesFMMM"))
time <- readRDS(system.file("test-data", "time.RDS", package = "BayesFMMM"))

## Set Hyperparameters
tot_mcmc_iters <- 150
n_try <- 1
k <- 2
n_funct <- 40
basis_degree <- 3
n_eigen <- 3
boundary_knots <- c(0, 1000)
internal_knots <- c(200, 400, 600, 800)
Y[[1]] <- Y[[1]][seq(1,100, 2)]
time[[1]] <- time[[1]][seq(1,100, 2)]
## Run function
x <- BFMMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                             basis_degree, n_eigen, boundary_knots,
                             internal_knots)
B <- x$B[[1]]