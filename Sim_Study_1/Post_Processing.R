library(BayesFMMM)
library(ggplot2)
library(gridExtra)
library(latex2exp)
library(pracma)
library(stringr)


#######################################################
### Adjust the directories based on the simulations ###
#######################################################

################################
### Covariate Adjusted #########
################################

### Set working dir
setwd("")

# Set number of MCMC runs
already_ran <- dir(paste0(getwd(), "/240_obs"))
ran <- which(paste0("sim", 1:50) %in% already_ran)
n_sim <- length(ran)

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

err_Z <- matrix(0, n_sim, 3)
int_err_mean1 <- matrix(0, n_sim, 3)
int_err_mean2 <- matrix(0, n_sim, 3)
int_err_cov12 <- matrix(0, n_sim, 3)
int_err_cov1 <- matrix(0, n_sim, 3)
int_err_cov2 <- matrix(0, n_sim, 3)


X_grid <- matrix(seq(-3, 3, 0.2), nrow = 31, ncol = 1)
########################################
## Uncomment below for 2 covariates ####
########################################

# X_grid <- matrix(0, nrow = 31*31, ncol = 2)
# x_1_seq <- seq(-3,3, 0.2)
# x_2_seq <- seq(-3,3, 0.2)
# index <- 1
# for(i in 1:length(x_1_seq)){
#   for(j in 1:length(x_2_seq)){
#     X_grid[index,1] <- x_1_seq[i]
#     X_grid[index,2] <- x_2_seq[j]
#     index <- index + 1
#   }
# }

norm_mu1 <- matrix(1, nrow = n_sim, ncol = 3)
norm_mu2 <- matrix(1, nrow = n_sim, ncol = 3)
norm_C1 <- matrix(1, nrow = n_sim, ncol = 3)
norm_C2 <- matrix(1, nrow = n_sim, ncol = 3)
norm_C12 <- matrix(1, nrow = n_sim, ncol = 3)

data_dir <- c("60_obs", "120_obs", "240_obs")

for(j in 1:3){
  if(j == 1){
    dir <- "./60_obs/"
  }
  if(j == 2){
    dir <- "./120_obs/"
  }
  if(j == 3){
    dir <- "./240_obs/"
  }
  for(i in 1:n_sim){
    x <- readRDS(paste("./", data_dir[j],"/sim", ran[i], "/truth.RDS", sep=""))
    x$Phi_true <- x$Phi
    x$nu_true <- x$nu
    Z_true <- x$Z
    nu_1_true <-  B %*% t(t(x$nu_true[1,]))
    nu_2_true <-  B %*% t(t(x$nu_true[2,]))
    mean_1_true <- X_grid %*% (t(x$eta[,,1]) %*% t(B))
    mean_2_true <- X_grid %*% t(x$eta[,,2]) %*% t(B)
    for(q in 1:nrow(X_grid)){
      mean_1_true[q,] = mean_1_true[q,] + t(nu_1_true)
      mean_2_true[q,] = mean_2_true[q,] + t(nu_2_true)
    }

    cov1_true <- matrix(0, 25, 25)
    for(q in 1:25){
      for(k in 1:25){
        for(m in 1:3){
          cov1_true[q,k] <- cov1_true[q,k] + x$Phi_true[1, ,m] %*% t(t(B[q,])) %*% B[k,] %*% t(t(x$Phi_true[1, ,m]))
        }
      }
    }

    cov2_true <- matrix(0, 25, 25)
    for(q in 1:25){
      for(k in 1:25){
        for(m in 1:3){
          cov2_true[q,k] <- cov2_true[q,k] + x$Phi_true[2, ,m] %*% t(t(B[q,])) %*% B[k,] %*% t(t(x$Phi_true[2, ,m]))
        }
      }
    }

    cov12_true <- matrix(0, 25, 25)
    for(q in 1:25){
      for(k in 1:25){
        for(m in 1:3){
          cov12_true[q,k] <- cov12_true[q,k] + x$Phi_true[1, ,m] %*% t(t(B[q,])) %*% B[k,] %*% t(t(x$Phi_true[2, ,m]))
        }
      }
    }

    dir_ph <- dir(paste0(dir, "sim", ran[i], "/"))
    n_files <- length(dir_ph[str_detect(dir_ph, "Z")])
    basis_degree <- 3
    boundary_knots <- c(0, 1000)
    internal_knots <- c(200, 400, 600, 800)
    dir_i <- paste(dir, "sim", ran[i], "/", sep = "")
    nu_1 <- FMeanCI_Adj(dir_i, n_files, time[[1]], X_grid, basis_degree, boundary_knots, internal_knots, 1, burnin_prop = 0.5)
    nu_2 <- FMeanCI_Adj(dir_i, n_files, time[[1]], X_grid, basis_degree, boundary_knots, internal_knots, 2, burnin_prop = 0.5)
    cov1 <- FCovCI(dir_i, n_files, 25, time[[1]], time[[1]], basis_degree,
                   boundary_knots, internal_knots, 1, 1, burnin_prop = 0.5)
    cov2 <- FCovCI(dir_i, n_files, 25, time[[1]], time[[1]], basis_degree,
                   boundary_knots, internal_knots, 2, 2, burnin_prop = 0.5)
    cov12 <- FCovCI(dir_i, n_files, 25, time[[1]], time[[1]], basis_degree,
                    boundary_knots, internal_knots, 1, 2, burnin_prop = 0.5)

    Z_est <- ZCI(dir_i, n_files)
    ## Label Switching
    if(sum(abs(nu_1$CI_50 - mean_1_true)) > sum(abs(nu_1$CI_50 - mean_2_true))){
      Z_est_i <- Z_est
      Z_est$CI_50[,1] <- Z_est_i$CI_50[,2]
      Z_est$CI_50[,2] <- Z_est_i$CI_50[,1]
      nu_i <- nu_1
      nu_1 <- nu_2
      nu_2 <- nu_i
      cov_i <- cov1
      cov1 <- cov2
      cov2 <- cov_i
      cov12$CI_50<- t(cov12$CI_50)
      cov12$CI_Upper<- t(cov12$CI_Upper)
      cov12$CI_Lower<- t(cov12$CI_Lower)
    }
    err_Z[i,j] <- sqrt(norm(Z_est$CI_50 - Z_true, "F") / (2 * nrow(Z_est$CI_50)))
    for(l in 1:nrow(X_grid)){
      int_err_mean1[i,j] <- int_err_mean1[i,j] + (sum((nu_1$CI_50[l,] - mean_1_true[l,])^2) * 20 * .2 * .2)
      int_err_mean2[i,j] <- int_err_mean2[i,j] + (sum((nu_2$CI_50[l,] - mean_2_true[l,])^2) * 20 * .2 * .2)
      norm_mu1[i,j] <- norm_mu1[i,j] + (sum((mean_1_true[l,])^2) * 20 * .2 * .2)
      norm_mu2[i,j] <- norm_mu2[i,j] + (sum((mean_2_true[l,])^2) * 20 * .2 * .2)
    }
    int_err_cov1[i,j] <- sum((cov1$CI_50 - cov1_true)^2) * 400
    int_err_cov2[i,j] <- sum((cov2$CI_50 - cov2_true)^2) * 400
    int_err_cov12[i,j] <- sum((cov12$CI_50 - cov12_true)^2) * 400
    norm_C1[i,j] <- sum((cov1_true)^2) * 400
    norm_C2[i,j] <- sum((cov2_true)^2) * 400
    norm_C12[i,j] <- sum((cov12_true)^2) * 400
    print(i)
    print(j)

  }
}

number_ticks <- function(n) {function(limits) 10^(pretty(log(limits,10), n))}
## Visualizations
C1_RMSE <- matrix(0, 3*n_sim, 2)
C1_RMSE[1:n_sim,1] <- (int_err_cov1[,1] / norm_C1[,1])
C1_RMSE[1:n_sim,2] <- 50
C1_RMSE[(n_sim + 1):(2*n_sim),1] <- (int_err_cov1[,2] / norm_C1[,2])
C1_RMSE[(n_sim + 1):(2*n_sim),2] <- 100
C1_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (int_err_cov1[,3] / norm_C1[,3])
C1_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 200
C1_RMSE <- as.data.frame(C1_RMSE)
colnames(C1_RMSE) <- c("R-MISE", "N")
C1_RMSE$N <- as.factor(C1_RMSE$N)
p1 <- ggplot(C1_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = number_ticks(6)) + ggtitle(TeX("$C^{(1,1)}$")) +
  geom_boxplot() +  theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                            legend.position = "none",plot.title = element_text(hjust = 0.5))
C2_RMSE <- matrix(0, 3*n_sim, 2)
C2_RMSE[1:n_sim,1] <- (int_err_cov2[,1] / norm_C2[,1])
C2_RMSE[1:n_sim,2] <- 50
C2_RMSE[(n_sim + 1):(2*n_sim),1] <- (int_err_cov2[,2] / norm_C2[,2])
C2_RMSE[(n_sim + 1):(2*n_sim),2] <- 100
C2_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (int_err_cov2[,3] / norm_C2[,3])
C2_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 200
C2_RMSE <- as.data.frame(C2_RMSE)
colnames(C2_RMSE) <- c("R-MISE", "N")
C2_RMSE$N <- as.factor(C2_RMSE$N)
p2 <- ggplot(C2_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = number_ticks(6)) + ggtitle(TeX("$C^{(2,2)}$")) +
  geom_boxplot() +  theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                            legend.position = "none",plot.title = element_text(hjust = 0.5))
C12_RMSE <- matrix(0, 3*n_sim, 2)
C12_RMSE[1:n_sim,1] <- (int_err_cov12[,1] / norm_C12[,1])
C12_RMSE[1:n_sim,2] <- 50
C12_RMSE[(n_sim + 1):(2*n_sim),1] <- (int_err_cov12[,2] / norm_C12[,2])
C12_RMSE[(n_sim + 1):(2*n_sim),2] <- 100
C12_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (int_err_cov12[,3] / norm_C12[,3])
C12_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 200
C12_RMSE <- as.data.frame(C12_RMSE)
colnames(C12_RMSE) <- c("R-MISE", "N")
C12_RMSE$N <- as.factor(C12_RMSE$N)
p3 <- ggplot(C12_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = number_ticks(6)) + ggtitle(TeX("$C^{(1,2)}$")) +
  geom_boxplot() +  theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                            legend.position = "none",plot.title = element_text(hjust = 0.5))
mu1_RMSE <- matrix(0, 3*n_sim, 2)
mu1_RMSE[1:n_sim,1] <- (int_err_mean1[,1] / norm_mu1[,1])
mu1_RMSE[1:n_sim,2] <- 50
mu1_RMSE[(n_sim + 1):(2*n_sim),1] <- (int_err_mean1[,2] / norm_mu1[,2])
mu1_RMSE[(n_sim + 1):(2*n_sim),2] <- 100
mu1_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (int_err_mean1[,3] / norm_mu1[,3])
mu1_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 200
mu1_RMSE <- as.data.frame(mu1_RMSE)
colnames(mu1_RMSE) <- c("R-MISE", "N")
mu1_RMSE$N <- as.factor(mu1_RMSE$N)
p4 <- ggplot(mu1_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = number_ticks(6))  + ggtitle(TeX("$\\mu_1$")) +
  geom_boxplot() +  theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                            legend.position = "none",plot.title = element_text(hjust = 0.5))
mu2_RMSE <- matrix(0, 3*n_sim, 2)
mu2_RMSE[1:n_sim,1] <- (int_err_mean2[,1] / norm_mu2[,1])
mu2_RMSE[1:n_sim,2] <- 50
mu2_RMSE[(n_sim + 1):(2*n_sim),1] <- (int_err_mean2[,2] / norm_mu2[,2])
mu2_RMSE[(n_sim + 1):(2*n_sim),2] <- 100
mu2_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (int_err_mean2[,3] / norm_mu2[,3])
mu2_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 200
mu2_RMSE <- as.data.frame(mu2_RMSE)
colnames(mu2_RMSE) <- c("R-MISE", "N")
mu2_RMSE$N <- as.factor(mu2_RMSE$N)
p5 <- ggplot(mu2_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = number_ticks(6)) + ggtitle(TeX("$\\mu_2$"))+
  geom_boxplot() +  theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                            legend.position = "none",plot.title = element_text(hjust = 0.5))


Z_RMSE <- matrix(0, 3*n_sim, 2)
Z_RMSE[1:n_sim,1] <- err_Z[,1]
Z_RMSE[1:n_sim,2] <- 50
Z_RMSE[(n_sim + 1):(2*n_sim),1] <- err_Z[,2]
Z_RMSE[(n_sim + 1):(2*n_sim),2] <- 100
Z_RMSE[(2*n_sim + 1):(3*n_sim),1] <- err_Z[,3]
Z_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 200
Z_RMSE <- as.data.frame(Z_RMSE)
colnames(Z_RMSE) <- c("RMSE", "N")
Z_RMSE$N <- as.factor(Z_RMSE$N)
p6 <- ggplot(Z_RMSE, aes(x=N, y=`RMSE`)) + ggtitle("Z") +
  geom_boxplot() +  theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                            legend.position = "none",plot.title = element_text(hjust = 0.5))
grid.arrange(p1, p2, p3, p4, p5, p6, layout_matrix = rbind(c(1, 1, 2, 2, 3, 3),
                                                           c(4, 4, 5, 5, 6, 6)))





################################
### Not Covariate Adjusted #####
################################

setwd("")

# Set number of MCMC runs
already_ran <- dir(paste0(getwd(), "/Unadjusted/200_obs"))
ran <- which(paste0("sim", 1:50) %in% already_ran)
n_sim <- length(ran)

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

err_Z <- matrix(0, n_sim, 3)
int_err_mean1 <- matrix(0, n_sim, 3)
int_err_mean2 <- matrix(0, n_sim, 3)
int_err_cov12 <- matrix(0, n_sim, 3)
int_err_cov1 <- matrix(0, n_sim, 3)
int_err_cov2 <- matrix(0, n_sim, 3)

norm_mu1 <- matrix(1, nrow = n_sim, ncol = 3)
norm_mu2 <- matrix(1, nrow = n_sim, ncol = 3)
norm_C1 <- matrix(1, nrow = n_sim, ncol = 3)
norm_C2 <- matrix(1, nrow = n_sim, ncol = 3)
norm_C12 <- matrix(1, nrow = n_sim, ncol = 3)

data_dir <- c("50_obs", "100_obs", "200_obs")

for(j in 1:3){
  if(j == 1){
    dir <- "./Unadjusted/50_obs/"
  }
  if(j == 2){
    dir <- "./Unadjusted/100_obs/"
  }
  if(j == 3){
    dir <- "./Unadjusted/200_obs/"
  }
  for(i in 1:n_sim){
    x <- readRDS(paste("./data/", data_dir[j],"/sim", ran[i], "/truth.RDS", sep=""))
    x$Phi_true <- x$Phi
    x$nu_true <- x$nu
    Z_true <- x$Z
    nu_1_true <-  B %*% t(t(x$nu_true[1,]))
    nu_2_true <-  B %*% t(t(x$nu_true[2,]))

    cov1_true <- matrix(0, 25, 25)
    for(q in 1:25){
      for(k in 1:25){
        for(m in 1:2){
          cov1_true[q,k] <- cov1_true[q,k] + x$Phi_true[1, ,m] %*% t(t(B[q,])) %*% B[k,] %*% t(t(x$Phi_true[1, ,m]))
        }
      }
    }

    cov2_true <- matrix(0, 25, 25)
    for(q in 1:25){
      for(k in 1:25){
        for(m in 1:2){
          cov2_true[q,k] <- cov2_true[q,k] + x$Phi_true[2, ,m] %*% t(t(B[q,])) %*% B[k,] %*% t(t(x$Phi_true[2, ,m]))
        }
      }
    }

    cov12_true <- matrix(0, 25, 25)
    for(q in 1:25){
      for(k in 1:25){
        for(m in 1:2){
          cov12_true[q,k] <- cov12_true[q,k] + x$Phi_true[1, ,m] %*% t(t(B[q,])) %*% B[k,] %*% t(t(x$Phi_true[2, ,m]))
        }
      }
    }


    n_files <- 120
    basis_degree <- 3
    boundary_knots <- c(0, 1000)
    internal_knots <- c(200, 400, 600, 800)
    dir_i <- paste(dir, "sim", ran[i], "/", sep = "")
    nu_1 <- FMeanCI(dir_i, n_files, time[[1]], basis_degree, boundary_knots, internal_knots, 1, burnin_prop = 0.5)
    nu_2 <- FMeanCI(dir_i, n_files, time[[1]], basis_degree, boundary_knots, internal_knots, 2, burnin_prop = 0.5)
    cov1 <- FCovCI(dir_i, n_files, 25, time[[1]], time[[1]], basis_degree,
                   boundary_knots, internal_knots, 1, 1, burnin_prop = 0.5)
    cov2 <- FCovCI(dir_i, n_files, 25, time[[1]], time[[1]], basis_degree,
                   boundary_knots, internal_knots, 2, 2, burnin_prop = 0.5)
    cov12 <- FCovCI(dir_i, n_files, 25, time[[1]], time[[1]], basis_degree,
                    boundary_knots, internal_knots, 1, 2, burnin_prop = 0.5)

    Z_est <- ZCI(dir_i, n_files)
    ## Label Switching
    if(sum(abs(nu_1$CI_50 - nu_1_true)) > sum(abs(nu_1$CI_50 - nu_2_true))){
      Z_est_i <- Z_est
      Z_est$CI_50[,1] <- Z_est_i$CI_50[,2]
      Z_est$CI_50[,2] <- Z_est_i$CI_50[,1]
      nu_i <- nu_1
      nu_1 <- nu_2
      nu_2 <- nu_i
      cov_i <- cov1
      cov1 <- cov2
      cov2 <- cov_i
      cov12$CI_50<- t(cov12$CI_50)
      cov12$CI_Upper<- t(cov12$CI_Upper)
      cov12$CI_Lower<- t(cov12$CI_Lower)
    }
    err_Z[i,j] <- sqrt(norm(Z_est$CI_50 - Z_true, "F") / (2 * nrow(Z_est$CI_50)))
    int_err_mean1[i,j] <- sum((nu_1$CI_50 - nu_1_true)^2) * 20
    int_err_mean2[i,j] <- sum((nu_2$CI_50 - nu_2_true)^2) * 20
    norm_mu1[i,j] <- sum((nu_1_true)^2) * 20
    norm_mu2[i,j] <- sum((nu_2_true)^2) * 20
    int_err_cov1[i,j] <- sum((cov1$CI_50 - cov1_true)^2) * 400
    int_err_cov2[i,j] <- sum((cov2$CI_50 - cov2_true)^2) * 400
    int_err_cov12[i,j] <- sum((cov12$CI_50 - cov12_true)^2) * 400
    norm_C1[i,j] <- sum((cov1_true)^2) * 400
    norm_C2[i,j] <- sum((cov2_true)^2) * 400
    norm_C12[i,j] <- sum((cov12_true)^2) * 400
    print(i)
    print(j)

  }
}

