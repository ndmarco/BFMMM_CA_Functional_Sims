library(BayesFMMM)
library(ggplot2)
library(gridExtra)
library(latex2exp)
library(pracma)
library(stringr)
library(ggpattern)
library(ggpubr)

#######################################################
### Adjust the directories based on the simulations ###
#######################################################

################################
### Covariate Adjusted #########
################################


################################
####### Two Parameters #########
################################

### Set working dir
setwd("./two_cov")

# Set number of MCMC runs
already_ran <- dir(paste0(getwd(), "/CA/240_obs"))
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

X_grid <- matrix(0, nrow = 31*31, ncol = 2)
x_1_seq <- seq(-3,3, 0.2)
x_2_seq <- seq(-3,3, 0.2)
index <- 1
for(i in 1:length(x_1_seq)){
  for(j in 1:length(x_2_seq)){
    X_grid[index,1] <- x_1_seq[i]
    X_grid[index,2] <- x_2_seq[j]
    index <- index + 1
  }
}

norm_mu1 <- matrix(1, nrow = n_sim, ncol = 3)
norm_mu2 <- matrix(1, nrow = n_sim, ncol = 3)
norm_C1 <- matrix(1, nrow = n_sim, ncol = 3)
norm_C2 <- matrix(1, nrow = n_sim, ncol = 3)
norm_C12 <- matrix(1, nrow = n_sim, ncol = 3)

data_dir <- c("data/60_obs", "data/120_obs", "data/240_obs")

for(j in 1:3){
  if(j == 1){
    dir <- "./CA/60_obs/"
  }
  if(j == 2){
    dir <- "./CA/120_obs/"
  }
  if(j == 3){
    dir <- "./CA/240_obs/"
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
    nu_1 <- FMeanCI(dir_i, n_files, time[[1]], basis_degree, boundary_knots, internal_knots, 1, burnin_prop = 0.5, X = X_grid)
    nu_2 <- FMeanCI(dir_i, n_files, time[[1]], basis_degree, boundary_knots, internal_knots, 2, burnin_prop = 0.5, X = X_grid)
    cov1 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
                   boundary_knots, internal_knots, 1, 1, burnin_prop = 0.5)
    cov2 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
                   boundary_knots, internal_knots, 2, 2, burnin_prop = 0.5)
    cov12 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
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
C1_RMSE[1:n_sim,2] <- 60
C1_RMSE[(n_sim + 1):(2*n_sim),1] <- (int_err_cov1[,2] / norm_C1[,2])
C1_RMSE[(n_sim + 1):(2*n_sim),2] <- 120
C1_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (int_err_cov1[,3] / norm_C1[,3])
C1_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 240
C1_RMSE <- as.data.frame(C1_RMSE)
colnames(C1_RMSE) <- c("R-MISE", "N")
C1_RMSE$N <- as.factor(C1_RMSE$N)
p12CA <- ggplot(C1_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(10, 0.5,  0.1, 0.02), minor_breaks = scales::pretty_breaks(n = 10)) + ggtitle(TeX("$C^{(1,1)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5))
C2_RMSE <- matrix(0, 3*n_sim, 2)
C2_RMSE[1:n_sim,1] <- (int_err_cov2[,1] / norm_C2[,1])
C2_RMSE[1:n_sim,2] <- 60
C2_RMSE[(n_sim + 1):(2*n_sim),1] <- (int_err_cov2[,2] / norm_C2[,2])
C2_RMSE[(n_sim + 1):(2*n_sim),2] <- 120
C2_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (int_err_cov2[,3] / norm_C2[,3])
C2_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 240
C2_RMSE <- as.data.frame(C2_RMSE)
colnames(C2_RMSE) <- c("R-MISE", "N")
C2_RMSE$N <- as.factor(C2_RMSE$N)
p22CA <- ggplot(C2_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(10, 0.5,  0.1, 0.02), minor_breaks = scales::pretty_breaks(n = 10)) + ggtitle(TeX("$C^{(2,2)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5))
C12_RMSE <- matrix(0, 3*n_sim, 2)
C12_RMSE[1:n_sim,1] <- (int_err_cov12[,1] / norm_C12[,1])
C12_RMSE[1:n_sim,2] <- 60
C12_RMSE[(n_sim + 1):(2*n_sim),1] <- (int_err_cov12[,2] / norm_C12[,2])
C12_RMSE[(n_sim + 1):(2*n_sim),2] <- 120
C12_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (int_err_cov12[,3] / norm_C12[,3])
C12_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 240
C12_RMSE <- as.data.frame(C12_RMSE)
colnames(C12_RMSE) <- c("R-MISE", "N")
C12_RMSE$N <- as.factor(C12_RMSE$N)
p32CA <- ggplot(C12_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(10, 0.5,  0.1, 0.02), minor_breaks = scales::pretty_breaks(n = 10)) + ggtitle(TeX("$C^{(1,2)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5))
mu1_RMSE <- matrix(0, 3*n_sim, 2)
mu1_RMSE[1:n_sim,1] <- (int_err_mean1[,1] / norm_mu1[,1])
mu1_RMSE[1:n_sim,2] <- 60
mu1_RMSE[(n_sim + 1):(2*n_sim),1] <- (int_err_mean1[,2] / norm_mu1[,2])
mu1_RMSE[(n_sim + 1):(2*n_sim),2] <- 120
mu1_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (int_err_mean1[,3] / norm_mu1[,3])
mu1_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 240
mu1_RMSE <- as.data.frame(mu1_RMSE)
colnames(mu1_RMSE) <- c("R-MISE", "N")
mu1_RMSE$N <- as.factor(mu1_RMSE$N)
p42CA <- ggplot(mu1_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(0.1, 0.02, 0.005, 0.001), minor_breaks = scales::pretty_breaks(n = 10))  + ggtitle(TeX("$\\mu_1$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5))
mu2_RMSE <- matrix(0, 3*n_sim, 2)
mu2_RMSE[1:n_sim,1] <- (int_err_mean2[,1] / norm_mu2[,1])
mu2_RMSE[1:n_sim,2] <- 60
mu2_RMSE[(n_sim + 1):(2*n_sim),1] <- (int_err_mean2[,2] / norm_mu2[,2])
mu2_RMSE[(n_sim + 1):(2*n_sim),2] <- 120
mu2_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (int_err_mean2[,3] / norm_mu2[,3])
mu2_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 240
mu2_RMSE <- as.data.frame(mu2_RMSE)
colnames(mu2_RMSE) <- c("R-MISE", "N")
mu2_RMSE$N <- as.factor(mu2_RMSE$N)
p52CA <- ggplot(mu2_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(0.1, 0.02, 0.005, 0.001), minor_breaks = scales::pretty_breaks(n = 10))  + ggtitle(TeX("$\\mu_2$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5))


Z_RMSE <- matrix(0, 3*n_sim, 2)
Z_RMSE[1:n_sim,1] <- err_Z[,1]
Z_RMSE[1:n_sim,2] <- 60
Z_RMSE[(n_sim + 1):(2*n_sim),1] <- err_Z[,2]
Z_RMSE[(n_sim + 1):(2*n_sim),2] <- 120
Z_RMSE[(2*n_sim + 1):(3*n_sim),1] <- err_Z[,3]
Z_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 240
Z_RMSE <- as.data.frame(Z_RMSE)
colnames(Z_RMSE) <- c("RMSE", "N")
Z_RMSE$N <- as.factor(Z_RMSE$N)
p62CA <- ggplot(Z_RMSE, aes(x=N, y=`RMSE`)) + ggtitle("Z") + scale_y_continuous( trans = 'log2',  breaks = c(0.08, 0.04, 0.02, 0.01), minor_breaks = scales::pretty_breaks(n = 10)) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5))
grid.arrange(p12CA, p22CA, p32CA, p42CA, p52CA, p62CA, layout_matrix = rbind(c(1, 1, 2, 2, 3, 3),
                                                           c(4, 4, 5, 5, 6, 6)))

CA2_res <- list("int_err_mean1" = int_err_mean1, "int_err_mean2" = int_err_mean2,
                "norm_mu1" = norm_mu1, "norm_mu2" = norm_mu2,
                "int_err_cov1" = int_err_cov1, "int_err_cov2" = int_err_cov2,
                "int_err_cov12" = int_err_cov12, "norm_C1" = norm_C1,
                "norm_C2" = norm_C2, "norm_C12" = norm_C12, "err_Z" = err_Z)


###############################
####### One Parameter #########
###############################

### Set working dir
setwd("./one_cov")

# Set number of MCMC runs
already_ran <- dir(paste0(getwd(), "/CA/200_obs"))
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

data_dir <- c("data/50_obs", "data/100_obs", "data/200_obs")

for(j in 1:3){
  if(j == 1){
    dir <- "./CA/50_obs/"
  }
  if(j == 2){
    dir <- "./CA/100_obs/"
  }
  if(j == 3){
    dir <- "./CA/200_obs/"
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
    nu_1 <- FMeanCI(dir_i, n_files, time[[1]], basis_degree, boundary_knots, internal_knots, 1, burnin_prop = 0.5, X = X_grid)
    nu_2 <- FMeanCI(dir_i, n_files, time[[1]], basis_degree, boundary_knots, internal_knots, 2, burnin_prop = 0.5, X = X_grid)
    cov1 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
                   boundary_knots, internal_knots, 1, 1, burnin_prop = 0.5)
    cov2 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
                   boundary_knots, internal_knots, 2, 2, burnin_prop = 0.5)
    cov12 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
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
p1CA1 <- ggplot(C1_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(10, 0.5,  0.1, 0.02), minor_breaks = scales::pretty_breaks(n = 10)) + ggtitle(TeX("$C^{(1,1)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
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
p2CA1 <- ggplot(C2_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(10, 0.5,  0.1, 0.02), minor_breaks = scales::pretty_breaks(n = 10)) + ggtitle(TeX("$C^{(2,2)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
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
p3CA1 <- ggplot(C12_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(10, 0.5,  0.1, 0.02), minor_breaks = scales::pretty_breaks(n = 10)) + ggtitle(TeX("$C^{(1,2)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
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
p4CA1 <- ggplot(mu1_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(0.1, 0.02, 0.005, 0.001), minor_breaks = scales::pretty_breaks(n = 10))  + ggtitle(TeX("$\\mu_1$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
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
p5CA1 <- ggplot(mu2_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(0.1, 0.02, 0.005, 0.001), minor_breaks = scales::pretty_breaks(n = 10)) + ggtitle(TeX("$\\mu_2$"))+
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
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
p6CA1 <- ggplot(Z_RMSE, aes(x=N, y=`RMSE`)) + ggtitle("Z") + scale_y_continuous( trans = 'log2',  breaks = c(0.08, 0.04, 0.02, 0.01), minor_breaks = scales::pretty_breaks(n = 10)) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                            legend.position = "none",plot.title = element_text(hjust = 0.5))
grid.arrange(p1CA1, p2CA1, p3CA1, p4CA1, p5CA1, p6CA1, layout_matrix = rbind(c(1, 1, 2, 2, 3, 3),
                                                           c(4, 4, 5, 5, 6, 6)))


CA1_res <- list("int_err_mean1" = int_err_mean1, "int_err_mean2" = int_err_mean2,
                "norm_mu1" = norm_mu1, "norm_mu2" = norm_mu2,
                "int_err_cov1" = int_err_cov1, "int_err_cov2" = int_err_cov2,
                "int_err_cov12" = int_err_cov12, "norm_C1" = norm_C1,
                "norm_C2" = norm_C2, "norm_C12" = norm_C12, "err_Z" = err_Z)

################################
### Not Covariate Adjusted #####
################################

setwd("./one_cov")

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

X_grid <- matrix(seq(-3, 3, 0.2), nrow = 31, ncol = 1)

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
    mean_1_true <- X_grid %*% (t(x$eta[,,1]) %*% t(B))
    mean_2_true <- X_grid %*% t(x$eta[,,2]) %*% t(B)
    for(q in 1:nrow(X_grid)){
      mean_1_true[q,] = mean_1_true[q,] + t(nu_1_true)
      mean_2_true[q,] = mean_2_true[q,] + t(nu_2_true)
    }

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

    dir_ph <- dir(paste0(dir, "sim", ran[i], "/"))
    n_files <- length(dir_ph[str_detect(dir_ph, "Z")])
    basis_degree <- 3
    boundary_knots <- c(0, 1000)
    internal_knots <- c(200, 400, 600, 800)
    dir_i <- paste(dir, "sim", ran[i], "/", sep = "")
    nu_1 <- FMeanCI(dir_i, n_files, time[[1]], basis_degree, boundary_knots, internal_knots, 1, burnin_prop = 0.5)
    nu_2 <- FMeanCI(dir_i, n_files, time[[1]], basis_degree, boundary_knots, internal_knots, 2, burnin_prop = 0.5)
    cov1 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
                   boundary_knots, internal_knots, 1, 1, burnin_prop = 0.5)
    cov2 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
                   boundary_knots, internal_knots, 2, 2, burnin_prop = 0.5)
    cov12 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
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
    for(l in 1:nrow(X_grid)){
      int_err_mean1[i,j] <- int_err_mean1[i,j] + (sum((nu_1$CI_50 - mean_1_true[l,])^2) * 20 * .2 * .2)
      int_err_mean2[i,j] <- int_err_mean2[i,j] + (sum((nu_2$CI_50 - mean_2_true[l,])^2) * 20 * .2 * .2)
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
p1CA0 <- ggplot(C1_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(10000, 1000,  100, 10), minor_breaks = scales::pretty_breaks(n = 10)) + ggtitle(TeX("$C^{(1,1)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
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
p2CA0 <- ggplot(C2_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(10000, 1000,  100, 10), minor_breaks = scales::pretty_breaks(n = 10)) + ggtitle(TeX("$C^{(2,2)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
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
p3CA0 <- ggplot(C12_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(10000, 1000,  100, 10), minor_breaks = scales::pretty_breaks(n = 10)) + ggtitle(TeX("$C^{(1,2)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
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
p4CA0 <- ggplot(mu1_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(10, 4, 1, 0.5, 0.25), minor_breaks = scales::pretty_breaks(n = 10))  + ggtitle(TeX("$\\mu_1$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
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
p5CA0 <- ggplot(mu2_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(10, 4, 1, 0.5, 0.25), minor_breaks = scales::pretty_breaks(n = 10)) + ggtitle(TeX("$\\mu_2$"))+
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
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
p6CA0 <- ggplot(Z_RMSE, aes(x=N, y=`RMSE`)) + ggtitle("Z") + scale_y_continuous( trans = 'log2',  breaks = c(0.12, 0.08, 0.04, 0.02), minor_breaks = scales::pretty_breaks(n = 10)) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5))
grid.arrange(p1CA0, p2CA0, p3CA0, p4CA0, p5CA0, p6CA0, layout_matrix = rbind(c(1, 1, 2, 2, 3, 3),
                                                           c(4, 4, 5, 5, 6, 6)))


Unadjusted1_res <- list("int_err_mean1" = int_err_mean1, "int_err_mean2" = int_err_mean2,
                "norm_mu1" = norm_mu1, "norm_mu2" = norm_mu2,
                "int_err_cov1" = int_err_cov1, "int_err_cov2" = int_err_cov2,
                "int_err_cov12" = int_err_cov12, "norm_C1" = norm_C1,
                "norm_C2" = norm_C2, "norm_C12" = norm_C12, "err_Z" = err_Z)

################# Zero True covariates #########################################

###############################
####### One Parameter #########
###############################

### Set working dir
setwd("./zero_cov")

# Set number of MCMC runs
already_ran <- dir(paste0(getwd(), "/CA/160_obs"))
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

norm_mu1 <- matrix(1, nrow = n_sim, ncol = 3)
norm_mu2 <- matrix(1, nrow = n_sim, ncol = 3)
norm_C1 <- matrix(1, nrow = n_sim, ncol = 3)
norm_C2 <- matrix(1, nrow = n_sim, ncol = 3)
norm_C12 <- matrix(1, nrow = n_sim, ncol = 3)

data_dir <- c("data/40_obs", "data/80_obs", "data/160_obs")

for(j in 1:3){
  if(j == 1){
    dir <- "./CA/40_obs/"
  }
  if(j == 2){
    dir <- "./CA/80_obs/"
  }
  if(j == 3){
    dir <- "./CA/160_obs/"
  }
  for(i in 1:n_sim){
    x <- readRDS(paste("./", data_dir[j],"/sim", ran[i], "/truth.RDS", sep=""))
    x$Phi_true <- x$Phi
    x$nu_true <- x$nu
    Z_true <- x$Z
    nu_1_true <-  B %*% t(t(x$nu_true[1,]))
    nu_2_true <-  B %*% t(t(x$nu_true[2,]))
    mean_1_true <- X_grid %*% (t(rep(0,8)) %*% t(B))
    mean_2_true <- X_grid %*% (t(rep(0,8)) %*% t(B))
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
    nu_1 <- FMeanCI(dir_i, n_files, time[[1]], basis_degree, boundary_knots, internal_knots, 1, burnin_prop = 0.5, X = X_grid)
    nu_2 <- FMeanCI(dir_i, n_files, time[[1]], basis_degree, boundary_knots, internal_knots, 2, burnin_prop = 0.5, X = X_grid)
    cov1 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
                   boundary_knots, internal_knots, 1, 1, burnin_prop = 0.5)
    cov2 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
                   boundary_knots, internal_knots, 2, 2, burnin_prop = 0.5)
    cov12 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
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
C1_RMSE[1:n_sim,2] <- 40
C1_RMSE[(n_sim + 1):(2*n_sim),1] <- (int_err_cov1[,2] / norm_C1[,2])
C1_RMSE[(n_sim + 1):(2*n_sim),2] <- 80
C1_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (int_err_cov1[,3] / norm_C1[,3])
C1_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 160
C1_RMSE <- as.data.frame(C1_RMSE)
colnames(C1_RMSE) <- c("R-MISE", "N")
C1_RMSE$N <- as.factor(C1_RMSE$N)
p1U1 <- ggplot(C1_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(10, 0.5,  0.1, 0.02), minor_breaks = scales::pretty_breaks(n = 10)) + ggtitle(TeX("$C^{(1,1)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5))
C2_RMSE <- matrix(0, 3*n_sim, 2)
C2_RMSE[1:n_sim,1] <- (int_err_cov2[,1] / norm_C2[,1])
C2_RMSE[1:n_sim,2] <- 40
C2_RMSE[(n_sim + 1):(2*n_sim),1] <- (int_err_cov2[,2] / norm_C2[,2])
C2_RMSE[(n_sim + 1):(2*n_sim),2] <- 80
C2_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (int_err_cov2[,3] / norm_C2[,3])
C2_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 160
C2_RMSE <- as.data.frame(C2_RMSE)
colnames(C2_RMSE) <- c("R-MISE", "N")
C2_RMSE$N <- as.factor(C2_RMSE$N)
p2U1 <- ggplot(C2_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(10, 0.5,  0.1, 0.02), minor_breaks = scales::pretty_breaks(n = 10)) + ggtitle(TeX("$C^{(2,2)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5))
C12_RMSE <- matrix(0, 3*n_sim, 2)
C12_RMSE[1:n_sim,1] <- (int_err_cov12[,1] / norm_C12[,1])
C12_RMSE[1:n_sim,2] <- 40
C12_RMSE[(n_sim + 1):(2*n_sim),1] <- (int_err_cov12[,2] / norm_C12[,2])
C12_RMSE[(n_sim + 1):(2*n_sim),2] <- 80
C12_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (int_err_cov12[,3] / norm_C12[,3])
C12_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 160
C12_RMSE <- as.data.frame(C12_RMSE)
colnames(C12_RMSE) <- c("R-MISE", "N")
C12_RMSE$N <- as.factor(C12_RMSE$N)
p3U1 <- ggplot(C12_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(10, 0.5,  0.1, 0.02), minor_breaks = scales::pretty_breaks(n = 10)) + ggtitle(TeX("$C^{(1,2)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5))
mu1_RMSE <- matrix(0, 3*n_sim, 2)
mu1_RMSE[1:n_sim,1] <- (int_err_mean1[,1] / norm_mu1[,1])
mu1_RMSE[1:n_sim,2] <- 40
mu1_RMSE[(n_sim + 1):(2*n_sim),1] <- (int_err_mean1[,2] / norm_mu1[,2])
mu1_RMSE[(n_sim + 1):(2*n_sim),2] <- 80
mu1_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (int_err_mean1[,3] / norm_mu1[,3])
mu1_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 160
mu1_RMSE <- as.data.frame(mu1_RMSE)
colnames(mu1_RMSE) <- c("R-MISE", "N")
mu1_RMSE$N <- as.factor(mu1_RMSE$N)
p4U1 <- ggplot(mu1_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(0.1, 0.02, 0.005, 0.001), minor_breaks = scales::pretty_breaks(n = 10))  + ggtitle(TeX("$\\mu_1$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5))
mu2_RMSE <- matrix(0, 3*n_sim, 2)
mu2_RMSE[1:n_sim,1] <- (int_err_mean2[,1] / norm_mu2[,1])
mu2_RMSE[1:n_sim,2] <- 40
mu2_RMSE[(n_sim + 1):(2*n_sim),1] <- (int_err_mean2[,2] / norm_mu2[,2])
mu2_RMSE[(n_sim + 1):(2*n_sim),2] <- 80
mu2_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (int_err_mean2[,3] / norm_mu2[,3])
mu2_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 160
mu2_RMSE <- as.data.frame(mu2_RMSE)
colnames(mu2_RMSE) <- c("R-MISE", "N")
mu2_RMSE$N <- as.factor(mu2_RMSE$N)
p5U1 <- ggplot(mu2_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(0.1, 0.02, 0.005, 0.001), minor_breaks = scales::pretty_breaks(n = 10)) + ggtitle(TeX("$\\mu_2$"))+
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5))


Z_RMSE <- matrix(0, 3*n_sim, 2)
Z_RMSE[1:n_sim,1] <- err_Z[,1]
Z_RMSE[1:n_sim,2] <- 40
Z_RMSE[(n_sim + 1):(2*n_sim),1] <- err_Z[,2]
Z_RMSE[(n_sim + 1):(2*n_sim),2] <- 80
Z_RMSE[(2*n_sim + 1):(3*n_sim),1] <- err_Z[,3]
Z_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 160
Z_RMSE <- as.data.frame(Z_RMSE)
colnames(Z_RMSE) <- c("RMSE", "N")
Z_RMSE$N <- as.factor(Z_RMSE$N)
p6U1<- ggplot(Z_RMSE, aes(x=N, y=`RMSE`)) + ggtitle("Z") + scale_y_continuous( trans = 'log2',  breaks = c(0.08, 0.04, 0.02, 0.01), minor_breaks = scales::pretty_breaks(n = 10)) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5))
grid.arrange(p1U1, p2U1, p3U1, p4U1, p5U1, p6U1, layout_matrix = rbind(c(1, 1, 2, 2, 3, 3),
                                                           c(4, 4, 5, 5, 6, 6)))


CA0_res <- list("int_err_mean1" = int_err_mean1, "int_err_mean2" = int_err_mean2,
                "norm_mu1" = norm_mu1, "norm_mu2" = norm_mu2,
                "int_err_cov1" = int_err_cov1, "int_err_cov2" = int_err_cov2,
                "int_err_cov12" = int_err_cov12, "norm_C1" = norm_C1,
                "norm_C2" = norm_C2, "norm_C12" = norm_C12, "err_Z" = err_Z)

################################
### Not Covariate Adjusted #####
################################

setwd("./zero_cov")

# Set number of MCMC runs
already_ran <- dir(paste0(getwd(), "/Unadjusted/160_obs"))
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

data_dir <- c("40_obs", "80_obs", "160_obs")

for(j in 1:3){
  if(j == 1){
    dir <- "./Unadjusted/40_obs/"
  }
  if(j == 2){
    dir <- "./Unadjusted/80_obs/"
  }
  if(j == 3){
    dir <- "./Unadjusted/160_obs/"
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
    nu_1 <- FMeanCI(dir_i, n_files, time[[1]], basis_degree, boundary_knots, internal_knots, 1, burnin_prop = 0.5)
    nu_2 <- FMeanCI(dir_i, n_files, time[[1]], basis_degree, boundary_knots, internal_knots, 2, burnin_prop = 0.5)
    cov1 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
                   boundary_knots, internal_knots, 1, 1, burnin_prop = 0.5)
    cov2 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
                   boundary_knots, internal_knots, 2, 2, burnin_prop = 0.5)
    cov12 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
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

number_ticks <- function(n) {function(limits) 10^(pretty(log(limits,10), n))}
## Visualizations
C1_RMSE <- matrix(0, 3*n_sim, 2)
C1_RMSE[1:n_sim,1] <- (int_err_cov1[,1] / norm_C1[,1])
C1_RMSE[1:n_sim,2] <- 40
C1_RMSE[(n_sim + 1):(2*n_sim),1] <- (int_err_cov1[,2] / norm_C1[,2])
C1_RMSE[(n_sim + 1):(2*n_sim),2] <- 80
C1_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (int_err_cov1[,3] / norm_C1[,3])
C1_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 160
C1_RMSE <- as.data.frame(C1_RMSE)
colnames(C1_RMSE) <- c("R-MISE", "N")
C1_RMSE$N <- as.factor(C1_RMSE$N)
p1U0 <- ggplot(C1_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(10, 0.5,  0.1, 0.02), minor_breaks = scales::pretty_breaks(n = 10)) + ggtitle(TeX("$C^{(1,1)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5))
C2_RMSE <- matrix(0, 3*n_sim, 2)
C2_RMSE[1:n_sim,1] <- (int_err_cov2[,1] / norm_C2[,1])
C2_RMSE[1:n_sim,2] <- 40
C2_RMSE[(n_sim + 1):(2*n_sim),1] <- (int_err_cov2[,2] / norm_C2[,2])
C2_RMSE[(n_sim + 1):(2*n_sim),2] <- 80
C2_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (int_err_cov2[,3] / norm_C2[,3])
C2_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 160
C2_RMSE <- as.data.frame(C2_RMSE)
colnames(C2_RMSE) <- c("R-MISE", "N")
C2_RMSE$N <- as.factor(C2_RMSE$N)
p2U0 <- ggplot(C2_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(10, 0.5,  0.1, 0.02), minor_breaks = scales::pretty_breaks(n = 10)) + ggtitle(TeX("$C^{(2,2)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5))
C12_RMSE <- matrix(0, 3*n_sim, 2)
C12_RMSE[1:n_sim,1] <- (int_err_cov12[,1] / norm_C12[,1])
C12_RMSE[1:n_sim,2] <- 40
C12_RMSE[(n_sim + 1):(2*n_sim),1] <- (int_err_cov12[,2] / norm_C12[,2])
C12_RMSE[(n_sim + 1):(2*n_sim),2] <- 80
C12_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (int_err_cov12[,3] / norm_C12[,3])
C12_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 160
C12_RMSE <- as.data.frame(C12_RMSE)
colnames(C12_RMSE) <- c("R-MISE", "N")
C12_RMSE$N <- as.factor(C12_RMSE$N)
p3U0 <- ggplot(C12_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(10, 0.5,  0.1, 0.02), minor_breaks = scales::pretty_breaks(n = 10)) + ggtitle(TeX("$C^{(1,2)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5))
mu1_RMSE <- matrix(0, 3*n_sim, 2)
mu1_RMSE[1:n_sim,1] <- (int_err_mean1[,1] / norm_mu1[,1])
mu1_RMSE[1:n_sim,2] <- 40
mu1_RMSE[(n_sim + 1):(2*n_sim),1] <- (int_err_mean1[,2] / norm_mu1[,2])
mu1_RMSE[(n_sim + 1):(2*n_sim),2] <- 80
mu1_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (int_err_mean1[,3] / norm_mu1[,3])
mu1_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 160
mu1_RMSE <- as.data.frame(mu1_RMSE)
colnames(mu1_RMSE) <- c("R-MISE", "N")
mu1_RMSE$N <- as.factor(mu1_RMSE$N)
p4U0 <- ggplot(mu1_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(0.1, 0.02, 0.005, 0.001), minor_breaks = scales::pretty_breaks(n = 10))  + ggtitle(TeX("$\\mu_1$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5))
mu2_RMSE <- matrix(0, 3*n_sim, 2)
mu2_RMSE[1:n_sim,1] <- (int_err_mean2[,1] / norm_mu2[,1])
mu2_RMSE[1:n_sim,2] <- 40
mu2_RMSE[(n_sim + 1):(2*n_sim),1] <- (int_err_mean2[,2] / norm_mu2[,2])
mu2_RMSE[(n_sim + 1):(2*n_sim),2] <- 80
mu2_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (int_err_mean2[,3] / norm_mu2[,3])
mu2_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 160
mu2_RMSE <- as.data.frame(mu2_RMSE)
colnames(mu2_RMSE) <- c("R-MISE", "N")
mu2_RMSE$N <- as.factor(mu2_RMSE$N)
p5U0 <- ggplot(mu2_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(0.1, 0.02, 0.005, 0.001), minor_breaks = scales::pretty_breaks(n = 10)) + ggtitle(TeX("$\\mu_2$"))+
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5))


Z_RMSE <- matrix(0, 3*n_sim, 2)
Z_RMSE[1:n_sim,1] <- err_Z[,1]
Z_RMSE[1:n_sim,2] <- 40
Z_RMSE[(n_sim + 1):(2*n_sim),1] <- err_Z[,2]
Z_RMSE[(n_sim + 1):(2*n_sim),2] <- 80
Z_RMSE[(2*n_sim + 1):(3*n_sim),1] <- err_Z[,3]
Z_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 160
Z_RMSE <- as.data.frame(Z_RMSE)
colnames(Z_RMSE) <- c("RMSE", "N")
Z_RMSE$N <- as.factor(Z_RMSE$N)
p6U0 <- ggplot(Z_RMSE, aes(x=N, y=`RMSE`)) + ggtitle("Z") + scale_y_continuous( trans = 'log2',  breaks = c(0.08, 0.04, 0.02, 0.01), minor_breaks = scales::pretty_breaks(n = 10)) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5))
grid.arrange(p1U0, p2U0, p3U0, p4U0, p5U0, p6U0, layout_matrix = rbind(c(1, 1, 2, 2, 3, 3),
                                                           c(4, 4, 5, 5, 6, 6)))


Unadjusted0_res <- list("int_err_mean1" = int_err_mean1, "int_err_mean2" = int_err_mean2,
                        "norm_mu1" = norm_mu1, "norm_mu2" = norm_mu2,
                        "int_err_cov1" = int_err_cov1, "int_err_cov2" = int_err_cov2,
                        "int_err_cov12" = int_err_cov12, "norm_C1" = norm_C1,
                        "norm_C2" = norm_C2, "norm_C12" = norm_C12, "err_Z" = err_Z)



#### Paper Figures

### Figure for 2 Covariates
mu_RMSE <- matrix(0, 6*n_sim, 3)
mu_RMSE[1:n_sim,1] <- (CA2_res$int_err_mean1[,1] / CA2_res$norm_mu1[,1])
mu_RMSE[1:n_sim,2] <- 60
mu_RMSE[1:n_sim,3] <- 2
mu_RMSE[(n_sim + 1):(2*n_sim),1] <- (CA2_res$int_err_mean1[,2] / CA2_res$norm_mu1[,2])
mu_RMSE[(n_sim + 1):(2*n_sim),2] <- 120
mu_RMSE[(n_sim + 1):(2*n_sim),3] <- 2
mu_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (CA2_res$int_err_mean1[,3] / CA2_res$norm_mu1[,3])
mu_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 240
mu_RMSE[(2*n_sim + 1):(3*n_sim),3] <- 2
mu_RMSE[(3*n_sim + 1):(4*n_sim),1] <- (CA2_res$int_err_mean2[,1] / CA2_res$norm_mu2[,1])
mu_RMSE[(3*n_sim + 1):(4*n_sim),2] <- 60
mu_RMSE[(3*n_sim + 1):(4*n_sim),3] <- 2
mu_RMSE[(4*n_sim + 1):(5*n_sim),1] <- (CA2_res$int_err_mean2[,2] / CA2_res$norm_mu2[,2])
mu_RMSE[(4*n_sim + 1):(5*n_sim),2] <- 120
mu_RMSE[(4*n_sim + 1):(5*n_sim),3] <- 2
mu_RMSE[(5*n_sim + 1):(6*n_sim),1] <- (CA2_res$int_err_mean2[,3] / CA2_res$norm_mu2[,3])
mu_RMSE[(5*n_sim + 1):(6*n_sim),2] <- 240
mu_RMSE[(5*n_sim + 1):(6*n_sim),3] <- 2
mu_RMSE <- as.data.frame(mu_RMSE)
colnames(mu_RMSE) <- c("R-MISE", "N", "# Covariates")
mu_RMSE$N <- as.factor(mu_RMSE$N)
mu_RMSE$`# Covariates` <- as.factor(mu_RMSE$`# Covariates`)
mu_CA2 <- ggplot(mu_RMSE, aes(x=N, y=`R-MISE`, pattern = `# Covariates`)) +
  geom_boxplot_pattern(aes(fill=`# Covariates`),
                       position = position_dodge(preserve = "single"),
                       fill = "lightgreen",
                       color = "black", 
                       pattern_fill = "black",
                       pattern_angle = 45,
                       pattern_density = 0.001,
                       pattern_spacing = 0.1,
                       pattern_key_scale_factor = 0.2) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(10, 1, 0.1, 0.01, 0.001), minor_breaks = scales::pretty_breaks(n = 10))  + ggtitle(TeX("Mean")) +
    theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))

cov_RMSE <- matrix(0, 9*n_sim, 3)
cov_RMSE[1:n_sim,1] <- (CA2_res$int_err_cov1[,1] / CA2_res$norm_C1[,1])
cov_RMSE[1:n_sim,2] <- 60
cov_RMSE[1:n_sim,3] <- 2
cov_RMSE[(n_sim + 1):(2*n_sim),1] <- (CA2_res$int_err_cov1[,2] / CA2_res$norm_C1[,2])
cov_RMSE[(n_sim + 1):(2*n_sim),2] <- 120
cov_RMSE[(n_sim + 1):(2*n_sim),3] <- 2
cov_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (CA2_res$int_err_cov1[,3] / CA2_res$norm_C1[,3])
cov_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 240
cov_RMSE[(2*n_sim + 1):(3*n_sim),3] <- 2
cov_RMSE[(3*n_sim + 1):(4*n_sim),1] <- (CA2_res$`int_err_cov2`[,1] / CA2_res$norm_C2[,1])
cov_RMSE[(3*n_sim + 1):(4*n_sim),2] <- 60
cov_RMSE[(3*n_sim + 1):(4*n_sim),3] <- 2
cov_RMSE[(4*n_sim + 1):(5*n_sim),1] <- (CA2_res$`int_err_cov2`[,2] / CA2_res$norm_C2[,2])
cov_RMSE[(4*n_sim + 1):(5*n_sim),2] <- 120
cov_RMSE[(4*n_sim + 1):(5*n_sim),3] <- 2
cov_RMSE[(5*n_sim + 1):(6*n_sim),1] <- (CA2_res$`int_err_cov2`[,3] / CA2_res$norm_C2[,3])
cov_RMSE[(5*n_sim + 1):(6*n_sim),2] <- 240
cov_RMSE[(5*n_sim + 1):(6*n_sim),3] <- 2
cov_RMSE[(6*n_sim + 1):(7*n_sim),1] <- (CA2_res$int_err_cov12[,1] / CA2_res$norm_C12[,1])
cov_RMSE[(6*n_sim + 1):(7*n_sim),2] <- 60
cov_RMSE[(6*n_sim + 1):(7*n_sim),3] <- 2
cov_RMSE[(7*n_sim + 1):(8*n_sim),1] <- (CA2_res$int_err_cov12[,2] / CA2_res$norm_C12[,2])
cov_RMSE[(7*n_sim + 1):(8*n_sim),2] <- 120
cov_RMSE[(7*n_sim + 1):(8*n_sim),3] <- 2
cov_RMSE[(8*n_sim + 1):(9*n_sim),1] <- (CA2_res$int_err_cov12[,3] / CA2_res$norm_C12[,3])
cov_RMSE[(8*n_sim + 1):(9*n_sim),2] <- 240
cov_RMSE[(8*n_sim + 1):(9*n_sim),3] <- 2
cov_RMSE <- as.data.frame(cov_RMSE)
colnames(cov_RMSE) <- c("R-MISE", "N", "# Covariates")
cov_RMSE$N <- as.factor(cov_RMSE$N)
cov_RMSE$`# Covariates` <- as.factor(cov_RMSE$`# Covariates`)
cov_CA2 <- ggplot(cov_RMSE, aes(x=N, y=`R-MISE`, pattern = `# Covariates`)) +
  geom_boxplot_pattern(aes(fill=`# Covariates`),
                       position = position_dodge(preserve = "single"),
                       fill = "lightgreen",
                       color = "black", 
                       pattern_fill = "black",
                       pattern_angle = 45,
                       pattern_density = 0.001,
                       pattern_spacing = 0.1,
                       pattern_key_scale_factor = 0.2) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(1000, 100, 10, 1, 0.1, 0.01), minor_breaks = scales::pretty_breaks(n = 10))  + ggtitle(TeX("Covariance")) +
  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))

Z_RMSE <- matrix(0, 3*n_sim, 3)
Z_RMSE[1:n_sim,1] <- CA2_res$err_Z[,1]
Z_RMSE[1:n_sim,2] <- 60
Z_RMSE[1:n_sim,3] <- 2
Z_RMSE[(n_sim + 1):(2*n_sim),1] <- CA2_res$err_Z[,2]
Z_RMSE[(n_sim + 1):(2*n_sim),2] <- 120
Z_RMSE[(n_sim + 1):(2*n_sim),3] <- 2
Z_RMSE[(2*n_sim + 1):(3*n_sim),1] <- CA2_res$err_Z[,3]
Z_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 240
Z_RMSE[(2*n_sim + 1):(3*n_sim),3] <- 2
Z_RMSE <- as.data.frame(Z_RMSE)
colnames(Z_RMSE) <- c("RMSE", "N", "# Covariates")
Z_RMSE$N <- as.factor(Z_RMSE$N)
Z_RMSE$`# Covariates` <- as.factor(Z_RMSE$`# Covariates`)
levels(Z_RMSE$`# Covariates`) <- c("2", "0")
Z_CA2 <- ggplot(Z_RMSE, aes(x=N, y=`RMSE`, pattern = `# Covariates`)) +
  geom_boxplot_pattern(aes(fill=`# Covariates`),
                       position = position_dodge(preserve = "single"),
                       fill = "lightgreen",
                       color = "black", 
                       pattern_fill = "black",
                       pattern_angle = 45,
                       pattern_density = 0.001,
                       pattern_spacing = 0.1,
                       pattern_key_scale_factor = 0.2) + ggtitle("Z") + scale_y_continuous( trans = 'log2',  breaks = c(0.16, 0.08, 0.04, 0.02, 0.01), minor_breaks = scales::pretty_breaks(n = 10)) +
  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5)) + scale_fill_discrete(guide = guide_legend()) 


### Figure for 1 Covariate
mu_RMSE <- matrix(0, 12*n_sim, 3)
mu_RMSE[1:n_sim,1] <- (CA1_res$int_err_mean1[,1] / CA1_res$norm_mu1[,1])
mu_RMSE[1:n_sim,2] <- 50
mu_RMSE[1:n_sim,3] <- 1
mu_RMSE[(n_sim + 1):(2*n_sim),1] <- (CA1_res$int_err_mean1[,2] / CA1_res$norm_mu1[,2])
mu_RMSE[(n_sim + 1):(2*n_sim),2] <- 100
mu_RMSE[(n_sim + 1):(2*n_sim),3] <- 1
mu_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (CA1_res$int_err_mean1[,3] / CA1_res$norm_mu1[,3])
mu_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 200
mu_RMSE[(2*n_sim + 1):(3*n_sim),3] <- 1
mu_RMSE[(3*n_sim + 1):(4*n_sim),1] <- (CA1_res$int_err_mean2[,1] / CA1_res$norm_mu2[,1])
mu_RMSE[(3*n_sim + 1):(4*n_sim),2] <- 50
mu_RMSE[(3*n_sim + 1):(4*n_sim),3] <- 1
mu_RMSE[(4*n_sim + 1):(5*n_sim),1] <- (CA1_res$int_err_mean2[,2] / CA1_res$norm_mu2[,2])
mu_RMSE[(4*n_sim + 1):(5*n_sim),2] <- 100
mu_RMSE[(4*n_sim + 1):(5*n_sim),3] <- 1
mu_RMSE[(5*n_sim + 1):(6*n_sim),1] <- (CA1_res$int_err_mean2[,3] / CA1_res$norm_mu2[,3])
mu_RMSE[(5*n_sim + 1):(6*n_sim),2] <- 200
mu_RMSE[(5*n_sim + 1):(6*n_sim),3] <- 1
mu_RMSE[(6*n_sim + 1):(7*n_sim),1] <- (Unadjusted1_res$int_err_mean1[,1] / Unadjusted1_res$norm_mu1[,1])
mu_RMSE[(6*n_sim + 1):(7*n_sim),2] <- 50
mu_RMSE[(6*n_sim + 1):(7*n_sim),3] <- 0
mu_RMSE[(7*n_sim + 1):(8*n_sim),1] <- (Unadjusted1_res$int_err_mean1[,2] / Unadjusted1_res$norm_mu1[,2])
mu_RMSE[(7*n_sim + 1):(8*n_sim),2] <- 100
mu_RMSE[(7*n_sim + 1):(8*n_sim),3] <- 0
mu_RMSE[(8*n_sim + 1):(9*n_sim),1] <- (Unadjusted1_res$int_err_mean1[,3] / Unadjusted1_res$norm_mu1[,3])
mu_RMSE[(8*n_sim + 1):(9*n_sim),2] <- 200
mu_RMSE[(8*n_sim + 1):(9*n_sim),3] <- 0
mu_RMSE[(9*n_sim + 1):(10*n_sim),1] <- (Unadjusted1_res$int_err_mean2[,1] / Unadjusted1_res$norm_mu2[,1])
mu_RMSE[(9*n_sim + 1):(10*n_sim),2] <- 50
mu_RMSE[(9*n_sim + 1):(10*n_sim),3] <- 0
mu_RMSE[(10*n_sim + 1):(11*n_sim),1] <- (Unadjusted1_res$int_err_mean2[,2] / Unadjusted1_res$norm_mu2[,2])
mu_RMSE[(10*n_sim + 1):(11*n_sim),2] <- 100
mu_RMSE[(10*n_sim + 1):(11*n_sim),3] <- 0
mu_RMSE[(11*n_sim + 1):(12*n_sim),1] <- (Unadjusted1_res$int_err_mean2[,3] / Unadjusted1_res$norm_mu2[,3])
mu_RMSE[(11*n_sim + 1):(12*n_sim),2] <- 200
mu_RMSE[(11*n_sim + 1):(12*n_sim),3] <- 0
mu_RMSE <- as.data.frame(mu_RMSE)
colnames(mu_RMSE) <- c("R-MISE", "N", "# Covariates")
mu_RMSE$N <- as.factor(mu_RMSE$N)
mu_RMSE$`# Covariates` <- as.factor(mu_RMSE$`# Covariates`)
mu_CA1 <- ggplot(mu_RMSE, aes(x=N, y=`R-MISE`, pattern = `# Covariates`))  +
  geom_boxplot_pattern(aes(fill=`# Covariates`),
                       position = position_dodge(preserve = "single"),
                       color = "black", 
                       pattern_fill = "black",
                       pattern_angle = 45,
                       pattern_density = 0.1,
                       pattern_spacing = 0.025,
                       pattern_key_scale_factor = 0.6) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(10, 1, 0.1, 0.01, 0.001), minor_breaks = scales::pretty_breaks(n = 10))  + ggtitle(TeX("Mean")) +
  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       plot.title = element_text(hjust = 0.5))

cov_RMSE <- matrix(0, 18*n_sim, 3)
cov_RMSE[1:n_sim,1] <- (CA1_res$int_err_cov1[,1] / CA1_res$norm_C1[,1])
cov_RMSE[1:n_sim,2] <- 50
cov_RMSE[1:n_sim,3] <- 1
cov_RMSE[(n_sim + 1):(2*n_sim),1] <- (CA1_res$int_err_cov1[,2] / CA1_res$norm_C1[,2])
cov_RMSE[(n_sim + 1):(2*n_sim),2] <- 100
cov_RMSE[(n_sim + 1):(2*n_sim),3] <- 1
cov_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (CA1_res$int_err_cov1[,3] / CA1_res$norm_C1[,3])
cov_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 200
cov_RMSE[(2*n_sim + 1):(3*n_sim),3] <- 1
cov_RMSE[(3*n_sim + 1):(4*n_sim),1] <- (CA1_res$`int_err_cov2`[,1] / CA1_res$norm_C2[,1])
cov_RMSE[(3*n_sim + 1):(4*n_sim),2] <- 50
cov_RMSE[(3*n_sim + 1):(4*n_sim),3] <- 1
cov_RMSE[(4*n_sim + 1):(5*n_sim),1] <- (CA1_res$`int_err_cov2`[,2] / CA1_res$norm_C2[,2])
cov_RMSE[(4*n_sim + 1):(5*n_sim),2] <- 100
cov_RMSE[(4*n_sim + 1):(5*n_sim),3] <- 1
cov_RMSE[(5*n_sim + 1):(6*n_sim),1] <- (CA1_res$`int_err_cov2`[,3] / CA1_res$norm_C2[,3])
cov_RMSE[(5*n_sim + 1):(6*n_sim),2] <- 200
cov_RMSE[(5*n_sim + 1):(6*n_sim),3] <- 1
cov_RMSE[(6*n_sim + 1):(7*n_sim),1] <- (CA1_res$int_err_cov12[,1] / CA1_res$norm_C12[,1])
cov_RMSE[(6*n_sim + 1):(7*n_sim),2] <- 50
cov_RMSE[(6*n_sim + 1):(7*n_sim),3] <- 1
cov_RMSE[(7*n_sim + 1):(8*n_sim),1] <- (CA1_res$int_err_cov12[,2] / CA1_res$norm_C12[,2])
cov_RMSE[(7*n_sim + 1):(8*n_sim),2] <- 100
cov_RMSE[(7*n_sim + 1):(8*n_sim),3] <- 1
cov_RMSE[(8*n_sim + 1):(9*n_sim),1] <- (CA1_res$int_err_cov12[,3] / CA1_res$norm_C12[,3])
cov_RMSE[(8*n_sim + 1):(9*n_sim),2] <- 200
cov_RMSE[(8*n_sim + 1):(9*n_sim),3] <- 1
cov_RMSE[(9*n_sim + 1):(10*n_sim),1] <- (Unadjusted1_res$int_err_cov1[,1] / Unadjusted1_res$norm_C1[,1])
cov_RMSE[(9*n_sim + 1):(10*n_sim),2] <- 50
cov_RMSE[(9*n_sim + 1):(10*n_sim),3] <- 0
cov_RMSE[(10*n_sim + 1):(11*n_sim),1] <- (Unadjusted1_res$int_err_cov1[,2] / Unadjusted1_res$norm_C1[,2])
cov_RMSE[(10*n_sim + 1):(11*n_sim),2] <- 100
cov_RMSE[(10*n_sim + 1):(11*n_sim),3] <- 0
cov_RMSE[(11*n_sim + 1):(12*n_sim),1] <- (Unadjusted1_res$int_err_cov1[,3] / Unadjusted1_res$norm_C1[,3])
cov_RMSE[(11*n_sim + 1):(12*n_sim),2] <- 200
cov_RMSE[(11*n_sim + 1):(12*n_sim),3] <- 0
cov_RMSE[(12*n_sim + 1):(13*n_sim),1] <- (Unadjusted1_res$`int_err_cov2`[,1] / Unadjusted1_res$norm_C2[,1])
cov_RMSE[(12*n_sim + 1):(13*n_sim),2] <- 50
cov_RMSE[(12*n_sim + 1):(13*n_sim),3] <- 0
cov_RMSE[(13*n_sim + 1):(14*n_sim),1] <- (Unadjusted1_res$`int_err_cov2`[,2] / Unadjusted1_res$norm_C2[,2])
cov_RMSE[(13*n_sim + 1):(14*n_sim),2] <- 100
cov_RMSE[(13*n_sim + 1):(14*n_sim),3] <- 0
cov_RMSE[(14*n_sim + 1):(15*n_sim),1] <- (Unadjusted1_res$`int_err_cov2`[,3] / Unadjusted1_res$norm_C2[,3])
cov_RMSE[(14*n_sim + 1):(15*n_sim),2] <- 200
cov_RMSE[(14*n_sim + 1):(15*n_sim),3] <- 0
cov_RMSE[(15*n_sim + 1):(16*n_sim),1] <- (Unadjusted1_res$int_err_cov12[,1] / Unadjusted1_res$norm_C12[,1])
cov_RMSE[(15*n_sim + 1):(16*n_sim),2] <- 50
cov_RMSE[(15*n_sim + 1):(16*n_sim),3] <- 0
cov_RMSE[(16*n_sim + 1):(17*n_sim),1] <- (Unadjusted1_res$int_err_cov12[,2] / Unadjusted1_res$norm_C12[,2])
cov_RMSE[(16*n_sim + 1):(17*n_sim),2] <- 100
cov_RMSE[(16*n_sim + 1):(17*n_sim),3] <- 0
cov_RMSE[(17*n_sim + 1):(18*n_sim),1] <- (Unadjusted1_res$int_err_cov12[,3] / Unadjusted1_res$norm_C12[,3])
cov_RMSE[(17*n_sim + 1):(18*n_sim),2] <- 200
cov_RMSE[(17*n_sim + 1):(18*n_sim) ,3] <- 0
cov_RMSE <- as.data.frame(cov_RMSE)
colnames(cov_RMSE) <- c("R-MISE", "N", "# Covariates")
cov_RMSE$N <- as.factor(cov_RMSE$N)
cov_RMSE$`# Covariates` <- as.factor(cov_RMSE$`# Covariates`)
cov_CA1 <- ggplot(cov_RMSE, aes(x=N, y=`R-MISE`, pattern = `# Covariates`)) +
  geom_boxplot_pattern(aes(fill=`# Covariates`),
                       position = position_dodge(preserve = "single"),
                       color = "black", 
                       pattern_fill = "black",
                       pattern_angle = 45,
                       pattern_density = 0.1,
                       pattern_spacing = 0.025,
                       pattern_key_scale_factor = 0.6) +
scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(1000, 100, 10, 1, 0.1, 0.01), minor_breaks = scales::pretty_breaks(n = 10))  + ggtitle(TeX("Covariance")) +
    theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       plot.title = element_text(hjust = 0.5)) 
Z_RMSE <- matrix(0, 6*n_sim, 3)
Z_RMSE[1:n_sim,1] <- CA1_res$err_Z[,1]
Z_RMSE[1:n_sim,2] <- 50
Z_RMSE[1:n_sim,3] <- 1
Z_RMSE[(n_sim + 1):(2*n_sim),1] <- CA1_res$err_Z[,2]
Z_RMSE[(n_sim + 1):(2*n_sim),2] <- 100
Z_RMSE[(n_sim + 1):(2*n_sim),3] <- 1
Z_RMSE[(2*n_sim + 1):(3*n_sim),1] <- CA1_res$err_Z[,3]
Z_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 200
Z_RMSE[(2*n_sim + 1):(3*n_sim),3] <- 1
Z_RMSE[(3*n_sim + 1):(4*n_sim),1] <- Unadjusted1_res$err_Z[,1]
Z_RMSE[(3*n_sim + 1):(4*n_sim),2] <- 50
Z_RMSE[(3*n_sim + 1):(4*n_sim),3] <- 0
Z_RMSE[(4*n_sim + 1):(5*n_sim),1] <- Unadjusted1_res$err_Z[,2]
Z_RMSE[(4*n_sim + 1):(5*n_sim),2] <- 100
Z_RMSE[(4*n_sim + 1):(5*n_sim),3] <- 0
Z_RMSE[(5*n_sim + 1):(6*n_sim),1] <- Unadjusted1_res$err_Z[,3]
Z_RMSE[(5*n_sim + 1):(6*n_sim),2] <- 200
Z_RMSE[(5*n_sim + 1):(6*n_sim),3] <- 0
Z_RMSE <- as.data.frame(Z_RMSE)
colnames(Z_RMSE) <- c("RMSE", "N", "# Covariates")
Z_RMSE$N <- as.factor(Z_RMSE$N)
Z_RMSE$`# Covariates` <- as.factor(Z_RMSE$`# Covariates`)
Z_CA1 <- ggplot(Z_RMSE, aes(x=N, y=`RMSE`, pattern = `# Covariates`)) +
  geom_boxplot_pattern(aes(fill=`# Covariates`),
                       position = position_dodge(preserve = "single"),
                       color = "black", 
                       pattern_fill = "black",
                       pattern_angle = 45,
                       pattern_density = 0.1,
                       pattern_spacing = 0.025,
                       pattern_key_scale_factor = 0.6) + ggtitle("Z") + scale_y_continuous( trans = 'log2',  breaks = c(0.16, 0.08, 0.04, 0.02, 0.01), minor_breaks = scales::pretty_breaks(n = 10)) +
  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       plot.title = element_text(hjust = 0.5))



### Figure for 0 Covariate
mu_RMSE <- matrix(0, 12*n_sim, 3)
mu_RMSE[1:n_sim,1] <- (CA0_res$int_err_mean1[,1] / CA0_res$norm_mu1[,1])
mu_RMSE[1:n_sim,2] <- 40
mu_RMSE[1:n_sim,3] <- 1
mu_RMSE[(n_sim + 1):(2*n_sim),1] <- (CA0_res$int_err_mean1[,2] / CA0_res$norm_mu1[,2])
mu_RMSE[(n_sim + 1):(2*n_sim),2] <- 80
mu_RMSE[(n_sim + 1):(2*n_sim),3] <- 1
mu_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (CA0_res$int_err_mean1[,3] / CA0_res$norm_mu1[,3])
mu_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 160
mu_RMSE[(2*n_sim + 1):(3*n_sim),3] <- 1
mu_RMSE[(3*n_sim + 1):(4*n_sim),1] <- (CA0_res$int_err_mean2[,1] / CA0_res$norm_mu2[,1])
mu_RMSE[(3*n_sim + 1):(4*n_sim),2] <- 40
mu_RMSE[(3*n_sim + 1):(4*n_sim),3] <- 1
mu_RMSE[(4*n_sim + 1):(5*n_sim),1] <- (CA0_res$int_err_mean2[,2] / CA0_res$norm_mu2[,2])
mu_RMSE[(4*n_sim + 1):(5*n_sim),2] <- 80
mu_RMSE[(4*n_sim + 1):(5*n_sim),3] <- 1
mu_RMSE[(5*n_sim + 1):(6*n_sim),1] <- (CA0_res$int_err_mean2[,3] / CA0_res$norm_mu2[,3])
mu_RMSE[(5*n_sim + 1):(6*n_sim),2] <- 160
mu_RMSE[(5*n_sim + 1):(6*n_sim),3] <- 1
mu_RMSE[(6*n_sim + 1):(7*n_sim),1] <- (Unadjusted0_res$int_err_mean1[,1] / Unadjusted0_res$norm_mu1[,1])
mu_RMSE[(6*n_sim + 1):(7*n_sim),2] <- 40
mu_RMSE[(6*n_sim + 1):(7*n_sim),3] <- 0
mu_RMSE[(7*n_sim + 1):(8*n_sim),1] <- (Unadjusted0_res$int_err_mean1[,2] / Unadjusted0_res$norm_mu1[,2])
mu_RMSE[(7*n_sim + 1):(8*n_sim),2] <- 80
mu_RMSE[(7*n_sim + 1):(8*n_sim),3] <- 0
mu_RMSE[(8*n_sim + 1):(9*n_sim),1] <- (Unadjusted0_res$int_err_mean1[,3] / Unadjusted0_res$norm_mu1[,3])
mu_RMSE[(8*n_sim + 1):(9*n_sim),2] <- 160
mu_RMSE[(8*n_sim + 1):(9*n_sim),3] <- 0
mu_RMSE[(9*n_sim + 1):(10*n_sim),1] <- (Unadjusted0_res$int_err_mean2[,1] / Unadjusted0_res$norm_mu2[,1])
mu_RMSE[(9*n_sim + 1):(10*n_sim),2] <- 40
mu_RMSE[(9*n_sim + 1):(10*n_sim),3] <- 0
mu_RMSE[(10*n_sim + 1):(11*n_sim),1] <- (Unadjusted0_res$int_err_mean2[,2] / Unadjusted0_res$norm_mu2[,2])
mu_RMSE[(10*n_sim + 1):(11*n_sim),2] <- 80
mu_RMSE[(10*n_sim + 1):(11*n_sim),3] <- 0
mu_RMSE[(11*n_sim + 1):(12*n_sim),1] <- (Unadjusted0_res$int_err_mean2[,3] / Unadjusted0_res$norm_mu2[,3])
mu_RMSE[(11*n_sim + 1):(12*n_sim),2] <- 160
mu_RMSE[(11*n_sim + 1):(12*n_sim),3] <- 0
mu_RMSE <- as.data.frame(mu_RMSE)
colnames(mu_RMSE) <- c("R-MISE", "N", "# Covariates")
mu_RMSE$N <- as.factor(mu_RMSE$N)
mu_RMSE$`# Covariates` <- as.factor(mu_RMSE$`# Covariates`)
mu_CA0 <- ggplot(mu_RMSE, aes(x=N, y=`R-MISE`, pattern = `# Covariates`))  +
  geom_boxplot_pattern(aes(fill=`# Covariates`),
                       position = position_dodge(preserve = "single"),
                       color = "black", 
                       pattern_fill = "black",
                       pattern_angle = 45,
                       pattern_density = 0.1,
                       pattern_spacing = 0.025,
                       pattern_key_scale_factor = 0.6) + scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(10, 1, 0.1, 0.01, 0.001), minor_breaks = scales::pretty_breaks(n = 10))  + ggtitle(TeX("Mean")) +
  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5))

cov_RMSE <- matrix(0, 18*n_sim, 3)
cov_RMSE[1:n_sim,1] <- (CA0_res$int_err_cov1[,1] / CA0_res$norm_C1[,1])
cov_RMSE[1:n_sim,2] <- 40
cov_RMSE[1:n_sim,3] <- 1
cov_RMSE[(n_sim + 1):(2*n_sim),1] <- (CA0_res$int_err_cov1[,2] / CA0_res$norm_C1[,2])
cov_RMSE[(n_sim + 1):(2*n_sim),2] <- 80
cov_RMSE[(n_sim + 1):(2*n_sim),3] <- 1
cov_RMSE[(2*n_sim + 1):(3*n_sim),1] <- (CA0_res$int_err_cov1[,3] / CA0_res$norm_C1[,3])
cov_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 160
cov_RMSE[(2*n_sim + 1):(3*n_sim),3] <- 1
cov_RMSE[(3*n_sim + 1):(4*n_sim),1] <- (CA0_res$`int_err_cov2`[,1] / CA0_res$norm_C2[,1])
cov_RMSE[(3*n_sim + 1):(4*n_sim),2] <- 40
cov_RMSE[(3*n_sim + 1):(4*n_sim),3] <- 1
cov_RMSE[(4*n_sim + 1):(5*n_sim),1] <- (CA0_res$`int_err_cov2`[,2] / CA0_res$norm_C2[,2])
cov_RMSE[(4*n_sim + 1):(5*n_sim),2] <- 80
cov_RMSE[(4*n_sim + 1):(5*n_sim),3] <- 1
cov_RMSE[(5*n_sim + 1):(6*n_sim),1] <- (CA0_res$`int_err_cov2`[,3] / CA0_res$norm_C2[,3])
cov_RMSE[(5*n_sim + 1):(6*n_sim),2] <- 160
cov_RMSE[(5*n_sim + 1):(6*n_sim),3] <- 1
cov_RMSE[(6*n_sim + 1):(7*n_sim),1] <- (CA0_res$int_err_cov12[,1] / CA0_res$norm_C12[,1])
cov_RMSE[(6*n_sim + 1):(7*n_sim),2] <- 40
cov_RMSE[(6*n_sim + 1):(7*n_sim),3] <- 1
cov_RMSE[(7*n_sim + 1):(8*n_sim),1] <- (CA0_res$int_err_cov12[,2] / CA0_res$norm_C12[,2])
cov_RMSE[(7*n_sim + 1):(8*n_sim),2] <- 80
cov_RMSE[(7*n_sim + 1):(8*n_sim),3] <- 1
cov_RMSE[(8*n_sim + 1):(9*n_sim),1] <- (CA0_res$int_err_cov12[,3] / CA0_res$norm_C12[,3])
cov_RMSE[(8*n_sim + 1):(9*n_sim),2] <- 160
cov_RMSE[(8*n_sim + 1):(9*n_sim),3] <- 1
cov_RMSE[(9*n_sim + 1):(10*n_sim),1] <- (Unadjusted0_res$int_err_cov1[,1] / Unadjusted0_res$norm_C1[,1])
cov_RMSE[(9*n_sim + 1):(10*n_sim),2] <- 40
cov_RMSE[(9*n_sim + 1):(10*n_sim),3] <- 0
cov_RMSE[(10*n_sim + 1):(11*n_sim),1] <- (Unadjusted0_res$int_err_cov1[,2] / Unadjusted0_res$norm_C1[,2])
cov_RMSE[(10*n_sim + 1):(11*n_sim),2] <- 80
cov_RMSE[(10*n_sim + 1):(11*n_sim),3] <- 0
cov_RMSE[(11*n_sim + 1):(12*n_sim),1] <- (Unadjusted0_res$int_err_cov1[,3] / Unadjusted0_res$norm_C1[,3])
cov_RMSE[(11*n_sim + 1):(12*n_sim),2] <- 160
cov_RMSE[(11*n_sim + 1):(12*n_sim),3] <- 0
cov_RMSE[(12*n_sim + 1):(13*n_sim),1] <- (Unadjusted0_res$`int_err_cov2`[,1] / Unadjusted0_res$norm_C2[,1])
cov_RMSE[(12*n_sim + 1):(13*n_sim),2] <- 40
cov_RMSE[(12*n_sim + 1):(13*n_sim),3] <- 0
cov_RMSE[(13*n_sim + 1):(14*n_sim),1] <- (Unadjusted0_res$`int_err_cov2`[,2] / Unadjusted0_res$norm_C2[,2])
cov_RMSE[(13*n_sim + 1):(14*n_sim),2] <- 80
cov_RMSE[(13*n_sim + 1):(14*n_sim),3] <- 0
cov_RMSE[(14*n_sim + 1):(15*n_sim),1] <- (Unadjusted0_res$`int_err_cov2`[,3] / Unadjusted0_res$norm_C2[,3])
cov_RMSE[(14*n_sim + 1):(15*n_sim),2] <- 160
cov_RMSE[(14*n_sim + 1):(15*n_sim),3] <- 0
cov_RMSE[(15*n_sim + 1):(16*n_sim),1] <- (Unadjusted0_res$int_err_cov12[,1] / Unadjusted0_res$norm_C12[,1])
cov_RMSE[(15*n_sim + 1):(16*n_sim),2] <- 40
cov_RMSE[(15*n_sim + 1):(16*n_sim),3] <- 0
cov_RMSE[(16*n_sim + 1):(17*n_sim),1] <- (Unadjusted0_res$int_err_cov12[,2] / Unadjusted0_res$norm_C12[,2])
cov_RMSE[(16*n_sim + 1):(17*n_sim),2] <- 80
cov_RMSE[(16*n_sim + 1):(17*n_sim),3] <- 0
cov_RMSE[(17*n_sim + 1):(18*n_sim),1] <- (Unadjusted0_res$int_err_cov12[,3] / Unadjusted0_res$norm_C12[,3])
cov_RMSE[(17*n_sim + 1):(18*n_sim),2] <- 160
cov_RMSE[(17*n_sim + 1):(18*n_sim) ,3] <- 0
cov_RMSE <- as.data.frame(cov_RMSE)
colnames(cov_RMSE) <- c("R-MISE", "N", "# Covariates")
cov_RMSE$N <- as.factor(cov_RMSE$N)
cov_RMSE$`# Covariates` <- as.factor(cov_RMSE$`# Covariates`)
cov_CA0 <- ggplot(cov_RMSE, aes(x=N, y=`R-MISE`, pattern = `# Covariates`)) +
  geom_boxplot_pattern(aes(fill=`# Covariates`),
                       position = position_dodge(preserve = "single"),
                       color = "black", 
                       pattern_fill = "black",
                       pattern_angle = 45,
                       pattern_density = 0.1,
                       pattern_spacing = 0.025,
                       pattern_key_scale_factor = 0.6) +
  scale_y_continuous(trans='log2', labels = scales::percent, breaks = c(1000, 100, 10, 1, 0.1, 0.01), minor_breaks = scales::pretty_breaks(n = 10))  + ggtitle(TeX("Covariance")) +
  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                     plot.title = element_text(hjust = 0.5)) 
Z_RMSE <- matrix(0, 6*n_sim, 3)
Z_RMSE[1:n_sim,1] <- CA0_res$err_Z[,1]
Z_RMSE[1:n_sim,2] <- 40
Z_RMSE[1:n_sim,3] <- 1
Z_RMSE[(n_sim + 1):(2*n_sim),1] <- CA0_res$err_Z[,2]
Z_RMSE[(n_sim + 1):(2*n_sim),2] <- 80
Z_RMSE[(n_sim + 1):(2*n_sim),3] <- 1
Z_RMSE[(2*n_sim + 1):(3*n_sim),1] <- CA0_res$err_Z[,3]
Z_RMSE[(2*n_sim + 1):(3*n_sim),2] <- 160
Z_RMSE[(2*n_sim + 1):(3*n_sim),3] <- 1
Z_RMSE[(3*n_sim + 1):(4*n_sim),1] <- Unadjusted0_res$err_Z[,1]
Z_RMSE[(3*n_sim + 1):(4*n_sim),2] <- 40
Z_RMSE[(3*n_sim + 1):(4*n_sim),3] <- 0
Z_RMSE[(4*n_sim + 1):(5*n_sim),1] <- Unadjusted0_res$err_Z[,2]
Z_RMSE[(4*n_sim + 1):(5*n_sim),2] <- 80
Z_RMSE[(4*n_sim + 1):(5*n_sim),3] <- 0
Z_RMSE[(5*n_sim + 1):(6*n_sim),1] <- Unadjusted0_res$err_Z[,3]
Z_RMSE[(5*n_sim + 1):(6*n_sim),2] <- 160
Z_RMSE[(5*n_sim + 1):(6*n_sim),3] <- 0
Z_RMSE <- as.data.frame(Z_RMSE)
colnames(Z_RMSE) <- c("RMSE", "N", "# Covariates")
Z_RMSE$N <- as.factor(Z_RMSE$N)
Z_RMSE$`# Covariates` <- as.factor(Z_RMSE$`# Covariates`)
Z_CA0 <- ggplot(Z_RMSE, aes(x=N, y=`RMSE`, pattern = `# Covariates`)) +
  geom_boxplot_pattern(aes(fill=`# Covariates`),
                       position = position_dodge(preserve = "single"),
                       color = "black", 
                       pattern_fill = "black",
                       pattern_angle = 45,
                       pattern_density = 0.1,
                       pattern_spacing = 0.025,
                       pattern_key_scale_factor = 0.6) + ggtitle("Z") + scale_y_continuous( trans = 'log2',  breaks = c(0.16, 0.08, 0.04, 0.02, 0.01), minor_breaks = scales::pretty_breaks(n = 10)) +
  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5))


x2 <- ggarrange(mu_CA2, cov_CA2, Z_CA2, ncol = 3, nrow = 1, common.legend = T, legend = "right")
x1 <- ggarrange(mu_CA1, cov_CA1, Z_CA1, ncol = 3, nrow = 1, common.legend = T, legend = "right")
x0 <- ggarrange(mu_CA0, cov_CA0, Z_CA0, ncol = 3, nrow = 1, common.legend = T, legend = "right")
 

grid.arrange(x2, x1, x0,  nrow = 3)
