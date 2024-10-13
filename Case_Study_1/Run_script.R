### Real Case study
library(BayesFMMM)
setwd(".")

#################################################################
## Change relevant directories and make folders before running ##
#################################################################

######################
#### Age Adjusted ####
######################

### Peak alpha data
library(pracma)
library(gridExtra)
# Subject ID
subj_id <- sort(c(10,	11,	13,	14,	15,	23,	26,	30,	31,	35,	48,	49,	50,
                  53,	54,	55,	161,165,	184,	188,	189,	195,	201,
                  # 202,	excluded due to low counts
                  207,	210,	213,	214,	242,	255,	261,	282,	283,
                  284,	286,	287,	289,	290,	343,	351,	2,	3,	5,	6,
                  7,	8,	9,	12,	18,	19,	22,	24,	25,	27,	33,	34,	37,	38,
                  40,	41,	42,	43,	44,	47,	51,	401,	405,	406,	408,	411,
                  415,	416,	417,	418,	423,	426,	427,	430,
                  #431,	excluded due to low counts
                  433,	436,	438,	439,	440,	442,	444,	445,	446,	447,
                  448,	450,	451,	452,	453,	3019,	3024,	3026,	3029,	3032))
# Channel ID (order of chan_id corresponds to 1:25 labeling of regions)
chan_id <- c('FP1', 'FP2','F9','F7','F3','Fz','F4','F8','F10','T9','T7',
             'C3','CZ','C4','T8','T10','P9','P7','P3','PZ','P4','P8','P10','O1','O2')

# Demographic Data
demDat <- read.csv(file='demographic_data.csv', header = TRUE)
colnames(demDat) <- c("ID", "Gender", "Age", "Group", "VIQ", "NVIQ")
demDat <- demDat[which(demDat$ID %in% subj_id), ]

# Peak Alpha Data
load("pa.dat.Rdata")
# ID: subject ID
# group: TD(1) or ASD (2)
# func: frequency domain
# reg: electrode (order corresponds to chan_id above)
# Age: age in months
# y: alpha spectra density
out1 <- unique(pa.dat$func)
out3 <- unique(pa.dat$reg)
matplot(matrix(pa.dat$y, nrow = length(out1)), type = "l") # data
trapz(out1, pa.dat$y[1:33]) # all functional observations integrate to 1 (normalized across electordes, subjects)

### Convert to wide format
data <- matrix(0,97,33)
freq <- seq(6,14,0.25)
for(i in 1:length(subj_id)){
  dat_i <- pa.dat[pa.dat$ID == subj_id[[i]], ]
  for(j in 1:length(freq)){
    data_ij <- dat_i[dat_i$func == freq[j],]
    data[i,j] = mean(data_ij$y)
  }
}

Y <- data

X <- matrix(log(demDat$Age), ncol = 1)


time <- seq(6, 14, 0.25)
matplot(t(matrix(rep(time,20),nrow =20, byrow = T)),t(Y[1:20,]), type = "l",
        xlab = "Frequency (Hz)", ylab = "Relative Power")

#get rid of ID value
library(reshape2)
library(ggplot2)

Y <- split(Y, seq(nrow(Y)))
time <- seq(6, 14, 0.25)
time <- rep(list(time), 97)


tot_mcmc_iters <- 4000
n_try <- 20
k <- 2
n_funct <- 97
basis_degree <- 3
n_eigen <- 3
boundary_knots <- c(6, 14)
internal_knots <- c(7, 8, 9, 10, 11, 12, 13)

## Get Estimates of Z and nu
est1 <- BFMMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                                basis_degree, n_eigen, boundary_knots,
                                internal_knots, X = X)
tot_mcmc_iters <- 20000
n_try <- 5
## Get estimates of other parameters
est2 <- BFMMM_Theta_est(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                        basis_degree, n_eigen, boundary_knots,
                        internal_knots, est1, X = X)

dir = "./Covariate_Adjusted_ASD/"
tot_mcmc_iters <- 500000
MCMC.chain <-BFMMM_warm_start(tot_mcmc_iters, k, Y, time, n_funct,
                              basis_degree, n_eigen, boundary_knots,
                              internal_knots, est1, est2, dir = dir, X = X,
                              thinning_num = 10, r_stored_iters = 10000)






################################################
### Age and clinical cohort adjusted BFMM ######
################################################

### Real Case study
library(BayesFMMM)


#################################################################
## Change relevant directories and make folders before running ##
#################################################################


### Peak alpha data
library(pracma)
library(gridExtra)
# Subject ID
subj_id <- sort(c(10,	11,	13,	14,	15,	23,	26,	30,	31,	35,	48,	49,	50,
                  53,	54,	55,	161,165,	184,	188,	189,	195,	201,
                  # 202,	excluded due to low counts
                  207,	210,	213,	214,	242,	255,	261,	282,	283,
                  284,	286,	287,	289,	290,	343,	351,	2,	3,	5,	6,
                  7,	8,	9,	12,	18,	19,	22,	24,	25,	27,	33,	34,	37,	38,
                  40,	41,	42,	43,	44,	47,	51,	401,	405,	406,	408,	411,
                  415,	416,	417,	418,	423,	426,	427,	430,
                  #431,	excluded due to low counts
                  433,	436,	438,	439,	440,	442,	444,	445,	446,	447,
                  448,	450,	451,	452,	453,	3019,	3024,	3026,	3029,	3032))
# Channel ID (order of chan_id corresponds to 1:25 labeling of regions)
chan_id <- c('FP1', 'FP2','F9','F7','F3','Fz','F4','F8','F10','T9','T7',
             'C3','CZ','C4','T8','T10','P9','P7','P3','PZ','P4','P8','P10','O1','O2')

# Demographic Data
demDat <- read.csv(file='demographic_data.csv', header = TRUE)
colnames(demDat) <- c("ID", "Gender", "Age", "Group", "VIQ", "NVIQ")
demDat <- demDat[which(demDat$ID %in% subj_id), ]

# Peak Alpha Data
load("pa.dat.Rdata")
# ID: subject ID
# group: TD(1) or ASD (2)
# func: frequency domain
# reg: electrode (order corresponds to chan_id above)
# Age: age in months
# y: alpha spectra density
out1 <- unique(pa.dat$func)
out3 <- unique(pa.dat$reg)
matplot(matrix(pa.dat$y, nrow = length(out1)), type = "l") # data
trapz(out1, pa.dat$y[1:33]) # all functional observations integrate to 1 (normalized across electrodes, subjects)

### Convert to wide format
data <- matrix(0,97,33)
freq <- seq(6,14,0.25)
for(i in 1:length(subj_id)){
  dat_i <- pa.dat[pa.dat$ID == subj_id[[i]], ]
  for(j in 1:length(freq)){
    data_ij <- dat_i[dat_i$func == freq[j],]
    data[i,j] = mean(data_ij$y)
  }
}

Y <- data
X <- cbind(log(demDat$Age), demDat$Group -1, rep(0, length(demDat$Age)))
X[,3][demDat$Group == 2] <- log(demDat$Age[demDat$Group == 2])

#get rid of ID value
library(reshape2)
library(ggplot2)

Y <- split(Y, seq(nrow(Y)))
time <- seq(6, 14, 0.25)
time <- rep(list(time), 97)

tot_mcmc_iters <- 4000
n_try <- 20
k <- 2
n_funct <- 97
basis_degree <- 3
n_eigen <- 3
boundary_knots <- c(6, 14)
internal_knots <- c(7, 8, 9, 10, 11, 12, 13)

## Get Estimates of Z and nu
est1 <- BFMMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                                basis_degree, n_eigen, boundary_knots,
                                internal_knots, X = X)
tot_mcmc_iters <- 20000
n_try <- 5
## Get estimates of other parameters
est2 <- BFMMM_Theta_est(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                        basis_degree, n_eigen, boundary_knots,
                        internal_knots, est1, X = X)

dir = "./Covariate_Adjusted_ASD2/"
tot_mcmc_iters <- 500000
MCMC.chain <-BFMMM_warm_start(tot_mcmc_iters, k, Y, time, n_funct,
                              basis_degree, n_eigen, boundary_knots,
                              internal_knots, est1, est2, dir = dir, X = X,
                              thinning_num = 10, r_stored_iters = 10000)

