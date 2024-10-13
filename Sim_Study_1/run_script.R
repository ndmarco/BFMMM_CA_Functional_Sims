library(BayesFMMM)
library(MASS)
library(DirichletReg)
library(future.apply)

####################
### 2 Covariates ###
####################

run_sim_mean_adj_2 <- function(iter){
  set.seed(iter)
  n_obs_vec <- c(60, 120, 240)
  for(n in 1:3){
    dir.create("data")
    dir.create(paste0("data/", n_obs_vec[n],"_obs"))
    dir.create(paste0("data/", n_obs_vec[n],"_obs/sim",iter))
    dir.create("CA")
    dir.create(paste0("CA/",n_obs_vec[n],"_obs"))
    dir.create(paste0("CA/",n_obs_vec[n],"_obs/sim",iter))
    ## Generate Data
    n_obs <- n_obs_vec[n]
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
    Y[[1]] <- Y[[1]][seq(1,100, 4)]
    time[[1]] <- time[[1]][seq(1,100, 4)]
    ## Run function
    x <- BFMMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                                 basis_degree, n_eigen, boundary_knots,
                                 internal_knots)
    B <- x$B[[1]]
    if(n == 1){
      nu <- matrix(0,nrow =2, ncol = 8)
      eta <- array(0, c(8,2,2))
      p <- matrix(0, 8, 8)
      for(i in 1:8){
        p[1,1] = 1
        if(i > 1){
          p[i,i] = 2
          p[i-1,i] = -1
          p[i,i-1] = -1
        }
      }

      nu[1,] <- mvrnorm(n=1, mu = seq(6,-8, -2), Sigma = 4*p)
      nu[2,] <- mvrnorm(n=1, mu = seq(-8, 6, 2), Sigma = 4*p)
      for(l in 1:2){
        eta[,1,l] <- mvrnorm(n=1, mu = rep(1, 8), Sigma = p)
        eta[,2,l] <- mvrnorm(n=1, mu = seq(3,-4, -1), Sigma = p)
      }

      Phi_i <- matrix(0, nrow = 2, ncol = 8)
      Phi_i[1,] <- t(rnorm(8, 0, 1.5))
      Phi_i[2,] <- t(rnorm(8, 0, 1.5))
      Phi_1 <- Phi_i
      Phi_i[1,] <- t(rnorm(8, 0, 1))
      Phi_i[2,] <- t(rnorm(8, 0, 1))
      Phi_2 <- Phi_i
      Phi_i[1,] <- t(rnorm(8, 0, 0.7))
      Phi_i[2,] <- t(rnorm(8, 0, 0.7))
      Phi_3 <- Phi_i
      Phi <- array(0, dim = c(2, 8, 3))
      Phi[,,1] <- matrix(Phi_1, nrow = 2)
      Phi[,,2] <- matrix(Phi_2, nrow = 2)
      Phi[,,3] <- matrix(Phi_3, nrow = 2)
    }
    if(n != 1){
      dat <- readRDS(paste("./data/", n_obs_vec[1],"_obs/sim",iter, "/truth.RDS", sep = ""))
      nu <- dat$nu
      eta <- dat$eta
      Phi <- dat$Phi
    }


    chi <- matrix(rnorm(n_obs *3, 0, 1), ncol = 3, nrow=n_obs)

    Z <- matrix(0, nrow = n_obs, ncol = 2)
    alpha <- c(10, 1)
    for(i in 1:(n_obs * 0.3)){
      Z[i,] <- rdirichlet(1, alpha)
    }
    alpha <- c(1, 10)
    for(i in (n_obs * 0.3 + 1):(n_obs * 0.6)){
      Z[i,] <- rdirichlet(1, alpha)
    }
    alpha <- c(1, 1)
    for(i in (n_obs * 0.6 + 1):n_obs){
      Z[i,] <- rdirichlet(1, alpha)
    }

    X <- matrix(rnorm(n = 2 * n_obs), ncol = 2, nrow = n_obs)
    y <- rep(0,25)
    y <- rep(list(y), n_obs)
    time <- time[[1]]
    time <- rep(list(time), n_obs)
    for(i in 1:n_obs){
      mean = rep(0,10)
      for(j in 1:2){
        mean = mean + Z[i,j] * B %*% t(nu[j,] + t(eta[, , j] %*% t(t(X[i,]))))
        for(m in 1:3){
          mean = mean + Z[i,j] * chi[i,m] * B %*% Phi[j, ,m]
        }
      }
      y[[i]] = mvrnorm(n = 1, mean, diag(0.001, 25))
    }

    x <- list("y" = y, "nu" = nu, "Z" = Z, "Phi" = Phi, "Chi" = chi, "eta" = eta, "X" = X)

    saveRDS(x, paste("./data/", n_obs_vec[n],"_obs/sim",iter, "/truth.RDS", sep = ""))


    ## Set Hyperparameters
    tot_mcmc_iters <- 2000
    n_try <- 50
    k <- 2
    n_funct <- n_obs
    basis_degree <- 3
    n_eigen <- 3
    boundary_knots <- c(0, 1000)
    internal_knots <- c(200, 400, 600, 800)
    Y <-y

    ## Get Estimates of Z and nu
    est1 <- BFMMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                                    basis_degree, n_eigen, boundary_knots,
                                    internal_knots, X= x$X)
    tot_mcmc_iters <- 20000
    n_try <- 10
    ## Get estimates of other parameters
    est2 <- BFMMM_Theta_est(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                            basis_degree, n_eigen, boundary_knots,
                            internal_knots, est1, X = x$X)
    dir_i <- paste("./CA/", n_obs_vec[n],"_obs/sim",iter, "/", sep="")
    tot_mcmc_iters <- 500000
    MCMC.chain <-BFMMM_warm_start(tot_mcmc_iters, k, Y, time, n_funct,
                                  basis_degree, n_eigen, boundary_knots,
                                  internal_knots, est1, est2, dir = dir_i, X = x$X,
                                  thinning_num = 100, r_stored_iters = 2500)
  }
}



##### Run Simulation

### Set working dir
setwd("./two_cov")

ncpu <- min(6, availableCores())
#
plan(multisession, workers = ncpu)

already_ran <- dir(paste0(getwd(), "/CA/240_obs"))
to_run <- which(!paste0("sim", 1:50) %in% already_ran)
seeds <- to_run
future_lapply(seeds, function(this_seed) run_sim_mean_adj_2(this_seed))



###################
### 1 Covariate ###
###################

run_sim_mean_adj_1 <- function(iter){
  set.seed(iter)
  n_obs_vec <- c(50, 100, 200)
  for(n in 1:3){
    dir.create("data")
    dir.create(paste0("data/", n_obs_vec[n],"_obs"))
    dir.create(paste0("data/", n_obs_vec[n],"_obs/sim",iter))
    dir.create("CA")
    dir.create(paste0("CA/", n_obs_vec[n],"_obs"))
    dir.create(paste0("CA/", n_obs_vec[n],"_obs/sim",iter))
    ## Generate Data
    n_obs <- n_obs_vec[n]
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
    Y[[1]] <- Y[[1]][seq(1,100, 4)]
    time[[1]] <- time[[1]][seq(1,100, 4)]
    ## Run function
    x <- BFMMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                                 basis_degree, n_eigen, boundary_knots,
                                 internal_knots)
    B <- x$B[[1]]
    if(n == 1){
      nu <- matrix(0,nrow =2, ncol = 8)
      eta <- array(0, c(8,1,2))
      p <- matrix(0, 8, 8)
      for(i in 1:8){
        p[1,1] = 1
        if(i > 1){
          p[i,i] = 2
          p[i-1,i] = -1
          p[i,i-1] = -1
        }
      }

      nu[1,] <- mvrnorm(n=1, mu = seq(6,-8, -2), Sigma = 4*p)
      nu[2,] <- mvrnorm(n=1, mu = seq(-8, 6, 2), Sigma = 4*p)
      eta[,1,1] <- mvrnorm(n=1, mu = rep(2, 8), Sigma = p)
      eta[,1,2] <- mvrnorm(n=1, mu = rep(-2, 8), Sigma = p)

      Phi_i <- matrix(0, nrow = 2, ncol = 8)
      Phi_i[1,] <- t(rnorm(8, 0, 1.5))
      Phi_i[2,] <- t(rnorm(8, 0, 1.5))
      Phi_1 <- Phi_i
      Phi_i[1,] <- t(rnorm(8, 0, 1))
      Phi_i[2,] <- t(rnorm(8, 0, 1))
      Phi_2 <- Phi_i
      Phi_i[1,] <- t(rnorm(8, 0, 0.7))
      Phi_i[2,] <- t(rnorm(8, 0, 0.7))
      Phi_3 <- Phi_i
      Phi <- array(0, dim = c(2, 8, 3))
      Phi[,,1] <- matrix(Phi_1, nrow = 2)
      Phi[,,2] <- matrix(Phi_2, nrow = 2)
      Phi[,,3] <- matrix(Phi_3, nrow = 2)
    }
    if(n != 1){
      dat <- readRDS(paste("./data/", n_obs_vec[1],"_obs/sim",iter, "/truth.RDS", sep = ""))
      nu <- dat$nu
      eta <- dat$eta
      Phi <- dat$Phi
    }

    chi <- matrix(rnorm(n_obs * 3, 0, 1), ncol = 3, nrow=n_obs)

    Z <- matrix(0, nrow = n_obs, ncol = 2)
    alpha <- c(10, 1)
    for(i in 1:(n_obs * 0.3)){
      Z[i,] <- rdirichlet(1, alpha)
    }
    alpha <- c(1, 10)
    for(i in (n_obs * 0.3 + 1):(n_obs * 0.6)){
      Z[i,] <- rdirichlet(1, alpha)
    }
    alpha <- c(1, 1)
    for(i in (n_obs * 0.6 + 1):n_obs){
      Z[i,] <- rdirichlet(1, alpha)
    }

    X <- matrix(rnorm(n = n_obs), ncol = 1, nrow = n_obs)
    y <- rep(0,25)
    y <- rep(list(y), n_obs)
    time <- time[[1]]
    time <- rep(list(time), n_obs)
    for(i in 1:n_obs){
      mean = rep(0,25)
      for(j in 1:2){
        mean = mean + Z[i,j] * B %*% t(nu[j,] + t(eta[, , j] %*% t(t(X[i,]))))
        for(m in 1:3){
          mean = mean + Z[i,j] * chi[i,m] * B %*% Phi[j, ,m]
        }
      }
      y[[i]] = mvrnorm(n = 1, mean, diag(0.001, 25))
    }

    x <- list("y" = y, "nu" = nu, "Z" = Z, "Phi" = Phi, "Chi" = chi, "eta" = eta, "X" = X)

    saveRDS(x, paste("./data/", n_obs_vec[n],"_obs/sim",iter, "/truth.RDS", sep = ""))


    ## Set Hyperparameters
    tot_mcmc_iters <- 2000
    n_try <- 10
    k <- 2
    n_funct <- n_obs
    basis_degree <- 3
    n_eigen <- 3
    boundary_knots <- c(0, 1000)
    internal_knots <- c(200, 400, 600, 800)
    Y <-y

    ## Get Estimates of Z and nu
    est1 <- BFMMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                                    basis_degree, n_eigen, boundary_knots,
                                    internal_knots, X = x$X)
    tot_mcmc_iters <- 20000
    n_try <- 10
    ## Get estimates of other parameters
    est2 <- BFMMM_Theta_est(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                            basis_degree, n_eigen, boundary_knots,
                            internal_knots, est1, X = x$X)
    dir_i <- paste("./CA/", n_obs_vec[n],"_obs/sim",iter, "/", sep="")
    tot_mcmc_iters <- 500000
    MCMC.chain <-BFMMM_warm_start(tot_mcmc_iters, k, Y, time, n_funct,
                                  basis_degree, n_eigen, boundary_knots,
                                  internal_knots, est1, est2, dir = dir_i, X = x$X,
                                  thinning_num = 100, r_stored_iters = 2500)
  }
  for(n in 1:3){
    n_obs <- n_obs_vec[n]
    dir.create("Unadjusted")
    dir.create(paste0("Unadjusted/", n_obs_vec[n],"_obs"))
    dir.create(paste0("Unadjusted/",n_obs_vec[n],"_obs/sim",iter))

    x <- readRDS(paste("./data/", n_obs_vec[n],"_obs/sim",iter, "/truth.RDS", sep = ""))

    ## Set Hyperparameters
    tot_mcmc_iters <- 2000
    n_try <- 10
    k <- 2
    n_funct <- n_obs
    basis_degree <- 3
    n_eigen <- 3
    boundary_knots <- c(0, 1000)
    internal_knots <- c(200, 400, 600, 800)
    Y <-x$y

    ## Get Estimates of Z and nu
    est1 <- BFMMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                                    basis_degree, n_eigen, boundary_knots,
                                    internal_knots)
    tot_mcmc_iters <- 20000
    n_try <- 10
    ## Get estimates of other parameters
    est2 <- BFMMM_Theta_est(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                            basis_degree, n_eigen, boundary_knots,
                            internal_knots, est1)
    dir_i <- paste("./Unadjusted/", n_obs_vec[n],"_obs/sim",iter, "/", sep="")
    tot_mcmc_iters <- 500000
    MCMC.chain <-BFMMM_warm_start(tot_mcmc_iters, k, Y, time, n_funct,
                                  basis_degree, n_eigen, boundary_knots,
                                  internal_knots, est1, est2, dir = dir_i,
                                  thinning_num = 100, r_stored_iters = 2500)
  }
}



##### Run Simulation

### Set working dir
setwd("./one_cov")

ncpu <- min(6, availableCores())
#
plan(multisession, workers = ncpu)

already_ran <- dir(paste0(getwd(), "/Unadjusted/200_obs"))
to_run <- which(!paste0("sim", 1:50) %in% already_ran)
seeds <- to_run
future_lapply(seeds, function(this_seed) run_sim_mean_adj_1(this_seed))



####################
### 0 Covariates ###
####################

run_sim_mean_adj_0 <- function(iter){
  set.seed(iter)
  n_obs_vec <- c(40, 80, 160)
  for(n in 1:3){
    dir.create("data")
    dir.create(paste0("data/", n_obs_vec[n],"_obs"))
    dir.create(paste0("data/", n_obs_vec[n],"_obs/sim",iter))
    dir.create("CA")
    dir.create(paste0("CA/", n_obs_vec[n],"_obs"))
    dir.create(paste0("CA/", n_obs_vec[n],"_obs/sim",iter))
    ## Generate Data
    n_obs <- n_obs_vec[n]
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
    Y[[1]] <- Y[[1]][seq(1,100, 4)]
    time[[1]] <- time[[1]][seq(1,100, 4)]
    ## Run function
    x <- BFMMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                                 basis_degree, n_eigen, boundary_knots,
                                 internal_knots)
    B <- x$B[[1]]

    if(n == 1){
      nu <- matrix(0,nrow =2, ncol = 8)
      p <- matrix(0, 8, 8)
      for(i in 1:8){
        p[1,1] = 1
        if(i > 1){
          p[i,i] = 2
          p[i-1,i] = -1
          p[i,i-1] = -1
        }
      }
      nu[1,] <- mvrnorm(n=1, mu = seq(6,-8, -2), Sigma = 4*p)
      nu[2,] <- mvrnorm(n=1, mu = seq(-8, 6, 2), Sigma = 4*p)

      Phi_i <- matrix(0, nrow = 2, ncol = 8)
      Phi_i[1,] <- t(rnorm(8, 0, 1.5))
      Phi_i[2,] <- t(rnorm(8, 0, 1.5))
      Phi_1 <- Phi_i
      Phi_i[1,] <- t(rnorm(8, 0, 1))
      Phi_i[2,] <- t(rnorm(8, 0, 1))
      Phi_2 <- Phi_i
      Phi_i[1,] <- t(rnorm(8, 0, 0.7))
      Phi_i[2,] <- t(rnorm(8, 0, 0.7))
      Phi_3 <- Phi_i
      Phi <- array(0, dim = c(2, 8, 3))
      Phi[,,1] <- matrix(Phi_1, nrow = 2)
      Phi[,,2] <- matrix(Phi_2, nrow = 2)
      Phi[,,3] <- matrix(Phi_3, nrow = 2)
    }
    if(n != 1){
      dat <- readRDS(paste("./data/", n_obs_vec[1],"_obs/sim",iter, "/truth.RDS", sep = ""))
      nu <- dat$nu
      Phi <- dat$Phi
    }



    chi <- matrix(rnorm(n_obs * 3, 0, 1), ncol = 3, nrow=n_obs)

    Z <- matrix(0, nrow = n_obs, ncol = 2)
    alpha <- c(10, 1)
    for(i in 1:(n_obs * 0.3)){
      Z[i,] <- rdirichlet(1, alpha)
    }
    alpha <- c(1, 10)
    for(i in (n_obs * 0.3 + 1):(n_obs * 0.6)){
      Z[i,] <- rdirichlet(1, alpha)
    }
    alpha <- c(1, 1)
    for(i in (n_obs * 0.6 + 1):n_obs){
      Z[i,] <- rdirichlet(1, alpha)
    }

    X <- matrix(rnorm(n = n_obs, 0, 2), ncol = 1, nrow = n_obs)
    y <- rep(0,25)
    y <- rep(list(y), n_obs)
    time <- time[[1]]
    time <- rep(list(time), n_obs)
    for(i in 1:n_obs){
      mean = rep(0,25)
      for(j in 1:2){
        mean = mean + Z[i,j] * B %*% nu[j,]
        for(m in 1:3){
          mean = mean + Z[i,j] * chi[i,m] * B %*% Phi[j, ,m]
        }
      }
      y[[i]] = mvrnorm(n = 1, mean, diag(0.001, 25))
    }

    x <- list("y" = y, "nu" = nu, "Z" = Z, "Phi" = Phi, "Chi" = chi, "X" = X)

    saveRDS(x, paste("./data/", n_obs_vec[n],"_obs/sim",iter, "/truth.RDS", sep = ""))


    ## Set Hyperparameters
    tot_mcmc_iters <- 2000
    n_try <- 50
    k <- 2
    n_funct <- n_obs
    basis_degree <- 3
    n_eigen <- 3
    boundary_knots <- c(0, 1000)
    internal_knots <- c(200, 400, 600, 800)
    Y <-y

    ## Get Estimates of Z and nu
    est1 <- BFMMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                                            basis_degree, n_eigen, boundary_knots,
                                            internal_knots, X = x$X)
    tot_mcmc_iters <- 20000
    n_try <- 10
    ## Get estimates of other parameters
    est2 <- BFMMM_Theta_est(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                            basis_degree, n_eigen, boundary_knots,
                            internal_knots, est1, X = x$X)
    dir_i <- paste("./CA/", n_obs_vec[n],"_obs/sim",iter, "/", sep="")
    tot_mcmc_iters <- 500000
    MCMC.chain <-BFMMM_warm_start(tot_mcmc_iters, k, Y, time, n_funct,
                                  basis_degree, n_eigen, boundary_knots,
                                  internal_knots, est1, est2, dir = dir_i,
                                  thinning_num = 100, r_stored_iters = 2500, X = x$X)
  }
  for(n in 1:3){
    n_obs <- n_obs_vec[n]
    dir.create("Unadjusted")
    dir.create(paste0("Unadjusted/", n_obs_vec[n],"_obs"))
    dir.create(paste0("Unadjusted/",n_obs_vec[n],"_obs/sim",iter))

    x <- readRDS(paste("./data/", n_obs_vec[n],"_obs/sim",iter, "/truth.RDS", sep = ""))

    ## Set Hyperparameters
    tot_mcmc_iters <- 2000
    n_try <- 50
    k <- 2
    n_funct <- n_obs
    basis_degree <- 3
    n_eigen <- 3
    boundary_knots <- c(0, 1000)
    internal_knots <- c(200, 400, 600, 800)
    Y <-x$y

    ## Get Estimates of Z and nu
    est1 <- BFMMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                                    basis_degree, n_eigen, boundary_knots,
                                    internal_knots)
    tot_mcmc_iters <- 20000
    n_try <- 10
    ## Get estimates of other parameters
    est2 <- BFMMM_Theta_est(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                            basis_degree, n_eigen, boundary_knots,
                            internal_knots, est1)
    dir_i <- paste("./Unadjusted/", n_obs_vec[n],"_obs/sim",iter, "/", sep="")
    tot_mcmc_iters <- 500000
    MCMC.chain <-BFMMM_warm_start(tot_mcmc_iters, k, Y, time, n_funct,
                                  basis_degree, n_eigen, boundary_knots,
                                  internal_knots, est1, est2, dir = dir_i,
                                  thinning_num = 100, r_stored_iters = 2500)
  }
}



##### Run Simulation

### Set working dir
setwd("./zero_cov")

ncpu <- min(6, availableCores())
#
plan(multisession, workers = ncpu)

already_ran <- dir(paste0(getwd(), "/Unadjusted/160_obs"))
to_run <- which(!paste0("sim", 1:50) %in% already_ran)
seeds <- to_run
future_lapply(seeds, function(this_seed) run_sim_mean_adj_0(this_seed))










