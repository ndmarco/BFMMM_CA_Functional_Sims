library(BayesFMMM)
library(latex2exp)
library(pracma)
library(gridExtra)
#################################################
## Change relevant directories  before running ##
#################################################
setwd("/Users/nicholasmarco/Documents")
### Peak alpha data

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
demDat <- read.csv(file='./demographic_data.csv', header = TRUE)
colnames(demDat) <- c("ID", "Gender", "Age", "Group", "VIQ", "NVIQ")
demDat <- demDat[which(demDat$ID %in% subj_id), ]

# Peak Alpha Data
load("./pa.dat.Rdata")
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
Y <- pa.dat
## paper used T8 electrode
Y <- Y[Y$reg == 15,]
Y$ID <- paste(Y$ID, Y$reg, sep = ".")
Y <- reshape(Y[,c(1,3,6)], idvar = "ID", timevar = "func", direction = "wide")
Y <- Y[,-1]
Y <- as.matrix(Y)

#get rid of ID value
library(reshape2)
library(ggplot2)

Y <- split(Y, seq(nrow(Y)))
time <- seq(6, 14, 0.25)
time <- rep(list(time), 97)
basis_degree <- 3
n_eigen <- 3
boundary_knots <- c(6, 14)
internal_knots <- c(7.6, 9.2, 10.8, 12.4)
time <- seq(6, 14, 0.05)
dir <- paste0(getwd(),"/Covariate_Adjusted_ASD/")
X <- matrix(seq(25,145, 10), ncol = 1)
### get credible intervals for mean
mean_1 <- FMeanCI_Adj(dir, 50, time, X, basis_degree, boundary_knots, internal_knots, 2)
mean_1s <- FMeanCI_Adj(dir, 50, time, X, basis_degree, boundary_knots, internal_knots, 2, simultaneous = T)

###############################################
### Plot trajectory of feature 1 for video ####
###############################################

for(i in 1:nrow(mean_1$CI_50)){
  df <- data.frame(freq = time,
                   median=mean_1$CI_50[i,],lwr=mean_1$CI_Lower[i,],
                   upr=mean_1$CI_Upper[i,], lwr_s = mean_1s$CI_Lower[i,],
                   upr_s = mean_1s$CI_Upper[i,])
  p1 <- ggplot(df, aes(freq, median))+
    geom_line(col = "blue")+
    geom_ribbon(data=df,aes(ymin=lwr,ymax=upr),alpha=0.3, fill = "black") + geom_ribbon(data=df,aes(ymin=lwr_s,ymax=upr_s),alpha=0.4, fill = "dark grey")  + ylab("Power") +
    xlab("Frequency (Hz)") + ylim(c(-0.2, 0.8)) + xlim(c(6,14)) + ggtitle(paste0("Feature 1 (", X[i], " months)")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0("feature1_", X[i],".png"),
         p1, device = "png", paste0(getwd(),"/Covariate_Adjusted_ASD/Images/"),
         width = 6, height = 4, units = "in", dpi = "print")
}



#############################################################################################################################

###############################################
### Plot trajectory of feature 1 over time ####
###############################################
feature1_df <- matrix(0, nrow = length(mean_1$CI_50), ncol = 3)
colnames(feature1_df) <- c("Frequency", "Age", "Power")
feature1_df <- as.data.frame(feature1_df)
for(i in 1:nrow(mean_1$CI_50)){
  feature1_df$Frequency[((i-1)*ncol(mean_1$CI_50) + 1):(i*ncol(mean_1$CI_50))] <- time
  feature1_df$Age[((i-1)*ncol(mean_1$CI_50) + 1):(i*ncol(mean_1$CI_50))] <- X[i,1]
  feature1_df$Power[((i-1)*ncol(mean_1$CI_50) + 1):(i*ncol(mean_1$CI_50))] <- mean_1$CI_50[i,]
}
feat1_plot <- ggplot(feature1_df, aes(x = Frequency, y =  Power, color = Age, group = Age)) + geom_line() + theme_bw() +
  ggtitle("Feature 1") + ylab("Relative Power") + xlab("Frequency (Hz)")
feat1_plot

#############################################################################################################################

mean_2 <- FMeanCI_Adj(dir, 50, time, X, basis_degree, boundary_knots, internal_knots, 1)
mean_2s <- FMeanCI_Adj(dir, 50, time, X, basis_degree, boundary_knots, internal_knots, 1, simultaneous = T)

###############################################
### Plot trajectory of feature 2 for video ####
###############################################

for(i in 1:nrow(mean_2$CI_50)){
  df <- data.frame(freq = time,
                   median=mean_2$CI_50[i,],lwr=mean_2$CI_Lower[i,],
                   upr=mean_2$CI_Upper[i,], lwr_s = mean_2s$CI_Lower[i,],
                   upr_s = mean_2s$CI_Upper[i,])
  p1 <- ggplot(df, aes(freq, median))+
    geom_line(col = "blue")+
    geom_ribbon(data=df,aes(ymin=lwr,ymax=upr),alpha=0.3, fill = "black") + geom_ribbon(data=df,aes(ymin=lwr_s,ymax=upr_s),alpha=0.4, fill = "dark grey")  + ylab("Power") +
    xlab("Frequency (Hz)") + ylim(c(-0.2, 0.8)) + xlim(c(6,14)) + ggtitle(paste0("Feature 2 (", X[i], " months)")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0("feature2_", X[i],".png"),
         p1, device = "png", paste0(getwd(),"/Covariate_Adjusted_ASD/Images/Feature_2/"),
         width = 6, height = 4, units = "in", dpi = "print")
}


###############################################
### Plot trajectory of feature 2 over time ####
###############################################
feature2_df <- matrix(0, nrow = length(mean_2$CI_50), ncol = 3)
colnames(feature2_df) <- c("Frequency", "Age", "Power")
feature2_df <- as.data.frame(feature2_df)
for(i in 1:nrow(mean_2$CI_50)){
  feature2_df$Frequency[((i-1)*ncol(mean_2$CI_50) + 1):(i*ncol(mean_2$CI_50))] <- time
  feature2_df$Age[((i-1)*ncol(mean_2$CI_50) + 1):(i*ncol(mean_2$CI_50))] <- X[i,1]
  feature2_df$Power[((i-1)*ncol(mean_2$CI_50) + 1):(i*ncol(mean_2$CI_50))] <- mean_2$CI_50[i,]
}
feat2_plot <- ggplot(feature2_df, aes(x = Frequency, y =  Power, color = Age, group = Age)) + geom_line() +
  theme_bw() + ggtitle("Feature 2") + ylab("Relative Power") + xlab("Frequency (Hz)")
feat2_plot

grid.arrange(feat1_plot, feat2_plot, ncol = 2)

#############################################################################################################################


Z_post <- ZCI(dir, 50)
Z_post$CI_50[,2][Z_post$CI_50[,2] < 0] <- 0
data_Z <- data.frame("Cluster 1" = Z_post$CI_50[,2], "Clinical Diagnosis" = demDat$Group)
data_Z$Clinical.Diagnosis[data_Z$Clinical.Diagnosis == 2] <- "ASD"
data_Z$Clinical.Diagnosis[data_Z$Clinical.Diagnosis == 1] <- "TD"

p3 <- ggplot(data= data_Z, aes(x = `Cluster.1` , y = Clinical.Diagnosis)) + geom_violin(trim = F, xlim = c(0,1)) + geom_point() + xlab("Feature 1") + ylab("Clinical Diagnosis") +
  stat_summary(
    geom = "point",
    fun.x = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "red")+ xlim(c(0,1)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                        panel.background = element_blank(),axis.line = element_line(colour = "black"),
                                        plot.title = element_text(hjust = 0.5))


grid.arrange(feat1_plot, feat2_plot, p3,  layout_matrix = rbind(c(1,2),c(1,2), c(1,2), c(3,3),c(3,3)))


data_VIQ <- data_Z <- data.frame("Feature 1" = Z_post$CI_50[,3], "VIQ" = demDat$VIQ)
p1 <- ggplot(data= data_Z, aes(x = `Cluster.1` , y = VIQ)) + geom_point() + xlab("Cluster 1") + ylab("Verbal IQ") +
  geom_smooth(method='lm', colour = "red") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                   panel.background = element_blank(),axis.line = element_line(colour = "black"),
                                                   plot.title = element_text(hjust = 0.5))
data_NVIQ <- data_Z <- data.frame("Feature 1" = Z_post$CI_50[,3], "NVIQ" = demDat$NVIQ)
p2 <- ggplot(data= data_Z, aes(x = `Cluster.1` , y = NVIQ)) + geom_point() + xlab("Cluster 1") + ylab("Nonverbal IQ") +
  geom_smooth(method='lm', colour = "red") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                   panel.background = element_blank(),axis.line = element_line(colour = "black"),
                                                   plot.title = element_text(hjust = 0.5))

grid.arrange(p1, p2, ncol = 2)


################################################
### Age and clinical cohort adjusted BFMM ######
################################################

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
demDat <- read.csv(file='./demographic_data.csv', header = TRUE)
colnames(demDat) <- c("ID", "Gender", "Age", "Group", "VIQ", "NVIQ")
demDat <- demDat[which(demDat$ID %in% subj_id), ]

# Peak Alpha Data
load("./pa.dat.Rdata")
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
Y <- pa.dat
## paper used T8 electrode
Y <- Y[Y$reg == 15,]
Y$ID <- paste(Y$ID, Y$reg, sep = ".")
Y <- reshape(Y[,c(1,3,6)], idvar = "ID", timevar = "func", direction = "wide")
Y <- Y[,-1]
Y <- as.matrix(Y)

#get rid of ID value
library(reshape2)
library(ggplot2)

Y <- split(Y, seq(nrow(Y)))
time <- seq(6, 14, 0.25)
time <- rep(list(time), 97)
basis_degree <- 3
n_eigen <- 3
boundary_knots <- c(6, 14)
internal_knots <- c(7.6, 9.2, 10.8, 12.4)
time <- seq(6, 14, 0.05)
dir <- paste0(getwd(),"/Covariate_Adjusted_ASD2/")
X <- cbind(seq(25,145, 10), rep(1, 13),seq(25,145, 10))
### get credible intervals for mean
mean_1 <- FMeanCI_Adj(dir, 50, time, X, basis_degree, boundary_knots, internal_knots, 2)
#mean_1s <- FMeanCI_Adj(dir, 50, time, X, basis_degree, boundary_knots, internal_knots, 2, simultaneous = T)


###############################################
### Plot trajectory of feature 1 over time ####
###############################################
feature1_df <- matrix(0, nrow = length(mean_1$CI_50), ncol = 3)
colnames(feature1_df) <- c("Frequency", "Age", "Power")
feature1_df <- as.data.frame(feature1_df)
for(i in 1:nrow(mean_1$CI_50)){
  feature1_df$Frequency[((i-1)*ncol(mean_1$CI_50) + 1):(i*ncol(mean_1$CI_50))] <- time
  feature1_df$Age[((i-1)*ncol(mean_1$CI_50) + 1):(i*ncol(mean_1$CI_50))] <- X[i,1]
  feature1_df$Power[((i-1)*ncol(mean_1$CI_50) + 1):(i*ncol(mean_1$CI_50))] <- mean_1$CI_50[i,]
}
feat1_plot_ASD <- ggplot(feature1_df, aes(x = Frequency, y =  Power, color = Age, group = Age)) + geom_line() +
  theme_bw() + ggtitle("Feature 1 (ASD)") + ylab("Relative Power") + xlab("Frequency (Hz)")
feat1_plot_ASD

#############################################################################################################################

mean_2 <- FMeanCI_Adj(dir, 50, time, X, basis_degree, boundary_knots, internal_knots, 1)
#mean_2s <- FMeanCI_Adj(dir, 50, time, X, basis_degree, boundary_knots, internal_knots, 1, simultaneous = T)

###############################################
### Plot trajectory of feature 2 over time ####
###############################################
feature2_df <- matrix(0, nrow = length(mean_2$CI_50), ncol = 3)
colnames(feature2_df) <- c("Frequency", "Age", "Power")
feature2_df <- as.data.frame(feature2_df)
for(i in 1:nrow(mean_2$CI_50)){
  feature2_df$Frequency[((i-1)*ncol(mean_2$CI_50) + 1):(i*ncol(mean_2$CI_50))] <- time
  feature2_df$Age[((i-1)*ncol(mean_2$CI_50) + 1):(i*ncol(mean_2$CI_50))] <- X[i,1]
  feature2_df$Power[((i-1)*ncol(mean_2$CI_50) + 1):(i*ncol(mean_2$CI_50))] <- mean_2$CI_50[i,]
}
feat2_plot_ASD <- ggplot(feature2_df, aes(x = Frequency, y =  Power, color = Age, group = Age)) +
  geom_line() + theme_bw() + ggtitle("Feature 2 (ASD)") + ylab("Relative Power") + xlab("Frequency (Hz)")
feat2_plot_ASD

###############################################
### Plot trajectory of mean of ASD over time ##
###############################################
feature2_df <- matrix(0, nrow = length(mean_2$CI_50), ncol = 3)
colnames(feature2_df) <- c("Frequency", "Age", "Power")
feature2_df <- as.data.frame(feature2_df)
for(i in 1:nrow(mean_2$CI_50)){
  feature2_df$Frequency[((i-1)*ncol(mean_2$CI_50) + 1):(i*ncol(mean_2$CI_50))] <- time
  feature2_df$Age[((i-1)*ncol(mean_2$CI_50) + 1):(i*ncol(mean_2$CI_50))] <- X[i,1]
  feature2_df$Power[((i-1)*ncol(mean_2$CI_50) + 1):(i*ncol(mean_2$CI_50))] <- ((0.4155785) * mean_1$CI_50[i,]) + ((1- 0.4155785) * mean_2$CI_50[i,])
}
mean_ASD <- ggplot(feature2_df, aes(x = Frequency, y =  Power, color = Age, group = Age)) +
  geom_line() + theme_bw() + ggtitle("Mean Alpha Oscillations for ASD") + ylab("Relative Power") + xlab("Frequency (Hz)") + ylim(c(0,0.3))
mean_ASD

####### TD Children #############

X <- cbind(seq(25,145, 10), rep(0, 13),rep(0, 13))
### get credible intervals for mean
mean_1 <- FMeanCI_Adj(dir, 50, time, X, basis_degree, boundary_knots, internal_knots, 2)
#mean_1s <- FMeanCI_Adj(dir, 50, time, X, basis_degree, boundary_knots, internal_knots, 2, simultaneous = T)


###############################################
### Plot trajectory of feature 1 over time ####
###############################################
feature1_df <- matrix(0, nrow = length(mean_1$CI_50), ncol = 3)
colnames(feature1_df) <- c("Frequency", "Age", "Power")
feature1_df <- as.data.frame(feature1_df)
for(i in 1:nrow(mean_1$CI_50)){
  feature1_df$Frequency[((i-1)*ncol(mean_1$CI_50) + 1):(i*ncol(mean_1$CI_50))] <- time
  feature1_df$Age[((i-1)*ncol(mean_1$CI_50) + 1):(i*ncol(mean_1$CI_50))] <- X[i,1]
  feature1_df$Power[((i-1)*ncol(mean_1$CI_50) + 1):(i*ncol(mean_1$CI_50))] <- mean_1$CI_50[i,]
}
feat1_plot_TD <- ggplot(feature1_df, aes(x = Frequency, y =  Power, color = Age, group = Age)) + geom_line() +
  theme_bw() + ggtitle("Feature 1 (TD)") + ylab("Relative Power") + xlab("Frequency (Hz)")
feat1_plot_TD

#############################################################################################################################

mean_2 <- FMeanCI_Adj(dir, 50, time, X, basis_degree, boundary_knots, internal_knots, 1)
#mean_2s <- FMeanCI_Adj(dir, 50, time, X, basis_degree, boundary_knots, internal_knots, 1, simultaneous = T)

###############################################
### Plot trajectory of feature 2 over time ####
###############################################
feature2_df <- matrix(0, nrow = length(mean_2$CI_50), ncol = 3)
colnames(feature2_df) <- c("Frequency", "Age", "Power")
feature2_df <- as.data.frame(feature2_df)
for(i in 1:nrow(mean_2$CI_50)){
  feature2_df$Frequency[((i-1)*ncol(mean_2$CI_50) + 1):(i*ncol(mean_2$CI_50))] <- time
  feature2_df$Age[((i-1)*ncol(mean_2$CI_50) + 1):(i*ncol(mean_2$CI_50))] <- X[i,1]
  feature2_df$Power[((i-1)*ncol(mean_2$CI_50) + 1):(i*ncol(mean_2$CI_50))] <- mean_2$CI_50[i,]
}
feat2_plot_TD <- ggplot(feature2_df, aes(x = Frequency, y =  Power, color = Age, group = Age)) +
  geom_line() + theme_bw() + ggtitle("Feature 2 (TD)") + ylab("Relative Power") + xlab("Frequency (Hz)")
feat2_plot_TD
#############################################################################################################################

###############################################
### Plot trajectory of mean of TD over time ###
###############################################
feature2_df <- matrix(0, nrow = length(mean_2$CI_50), ncol = 3)
colnames(feature2_df) <- c("Frequency", "Age", "Power")
feature2_df <- as.data.frame(feature2_df)
for(i in 1:nrow(mean_2$CI_50)){
  feature2_df$Frequency[((i-1)*ncol(mean_2$CI_50) + 1):(i*ncol(mean_2$CI_50))] <- time
  feature2_df$Age[((i-1)*ncol(mean_2$CI_50) + 1):(i*ncol(mean_2$CI_50))] <- X[i,1]
  feature2_df$Power[((i-1)*ncol(mean_2$CI_50) + 1):(i*ncol(mean_2$CI_50))] <- ((0.3826806) * mean_1$CI_50[i,]) + ((1- 0.3826806) * mean_2$CI_50[i,])
}
mean_TD <- ggplot(feature2_df, aes(x = Frequency, y =  Power, color = Age, group = Age)) +
  geom_line() + theme_bw() + ggtitle("Mean Alpha Oscillations for TD") + ylab("Relative Power") + xlab("Frequency (Hz)") + ylim(c(0,0.3))
mean_TD


grid.arrange(mean_ASD, mean_TD, nrow = 1)


Z_post <- ZCI(dir, 50)
Z_post$CI_50[,2][Z_post$CI_50[,2] < 0] <- 0
data_Z <- data.frame("Cluster 1" = Z_post$CI_50[,2], "Clinical Diagnosis" = demDat$Group)
data_Z$Clinical.Diagnosis[data_Z$Clinical.Diagnosis == 2] <- "ASD"
data_Z$Clinical.Diagnosis[data_Z$Clinical.Diagnosis == 1] <- "TD"

p3 <- ggplot(data= data_Z, aes(x = `Cluster.1` , y = Clinical.Diagnosis)) + geom_violin(trim = F, xlim = c(0,1)) + geom_point() + xlab("Feature 1") + ylab("Clinical Diagnosis") +
  stat_summary(
    geom = "point",
    fun.x = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "red")+ xlim(c(0,1)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                        panel.background = element_blank(),axis.line = element_line(colour = "black"),
                                        plot.title = element_text(hjust = 0.5))

grid.arrange(feat1_plot_ASD, feat2_plot_ASD, feat1_plot_TD, feat2_plot_TD, p3,
             layout_matrix = rbind(c(1,2), c(1,2), c(1,2),c(3,4), c(3,4), c(3,4), c(5,5),c(5,5)))



####################
## Calculate CPO ###
####################

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
demDat <- read.csv(file='./demographic_data.csv', header = TRUE)
colnames(demDat) <- c("ID", "Gender", "Age", "Group", "VIQ", "NVIQ")
demDat <- demDat[which(demDat$ID %in% subj_id), ]

# Peak Alpha Data
load("./pa.dat.Rdata")
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
Y <- pa.dat
## paper used T8 electrode
Y <- Y[Y$reg == 15,]
Y$ID <- paste(Y$ID, Y$reg, sep = ".")
Y <- reshape(Y[,c(1,3,6)], idvar = "ID", timevar = "func", direction = "wide")
Y <- Y[,-1]
Y <- as.matrix(Y)

Y <- split(Y, seq(nrow(Y)))
time <- seq(6, 14, 0.25)
time <- rep(list(time), 97)

setwd("/Users/nicholasmarco/Documents")
time <- seq(6, 14, 0.25)
time <- rep(list(time), 97)
basis_degree <- 3
n_eigen <- 3
boundary_knots <- c(6, 14)
internal_knots <- c(7.6, 9.2, 10.8, 12.4)
dir = "./trace/"
CPO <- Conditional_Predictive_Ordinates(dir, 50, 1000, basis_degree, boundary_knots,
                                        internal_knots, time, Y)

time <- seq(6, 14, 0.25)
time <- rep(list(time), 97)
basis_degree <- 3
n_eigen <- 3
boundary_knots <- c(6, 14)
internal_knots <- c(7.6, 9.2, 10.8, 12.4)
X <- matrix(demDat$Age, ncol = 1)
dir = "./Covariate_Adjusted_ASD/"
CPO1 <- Conditional_Predictive_Ordinates(dir, 50, 1000, basis_degree, boundary_knots,
                                        internal_knots, time, Y, X = X)

X <- cbind(demDat$Age, demDat$Group -1, rep(0, length(demDat$Age)))
X[,3][demDat$Group == 2] <- demDat$Age[demDat$Group == 2]
dir = "./Covariate_Adjusted_ASD2/"
CPO2 <- Conditional_Predictive_Ordinates(dir, 50, 1000, basis_degree, boundary_knots,
                                         internal_knots, time, Y, X = X)

X <- matrix(demDat$Age, ncol = 1)
dir = "./Covariate_Adjusted_ASD_Covariance/"
CPO3 <- Conditional_Predictive_Ordinates(dir, 50, 1000, basis_degree, boundary_knots,
                                         internal_knots, time, Y, X = X, cov_adj = T)

#### Plots of CPO
df <- matrix(0, nrow = 97, ncol = 3)
colnames(df) <- c("CPO1", "CP02", "Clinical Group")
df <- as.data.frame(df)
df[,1] <- CPO
df[,2] <- CPO1
df[,3] <- "ASD"
df[demDat$Group==2,3] <-"TD"

p1 <- ggplot(data = df, aes(x = `CPO2`, y = `CPO1`, color = `Clinical Group`)) + geom_abline(slope=1,intercept = 0) +
  geom_point() + theme_bw() + xlab(TeX("$M_{1}\\;\\;\\log{(CPO)}$")) + ylab(TeX("$M_{0}\\;\\;\\log{(CPO)}$"))

df <- matrix(0, nrow = 97, ncol = 3)
colnames(df) <- c("CPO1", "CP02", "Clinical Group")
df <- as.data.frame(df)
df[,1] <- CPO1
df[,2] <- CPO2
df[,3] <- "ASD"
df[demDat$Group==2,3] <-"TD"

p2 <- ggplot(data = df, aes(x = `CPO2`, y = `CPO1`, color = `Clinical Group`)) + geom_abline(slope=1,intercept = 0) +
  geom_point() + theme_bw() + xlab(TeX("$M_{2}\\;\\;\\log{(CPO)}$")) + ylab(TeX("$M_{1}\\;\\;\\log{(CPO)}$"))

df <- matrix(0, nrow = 97, ncol = 3)
colnames(df) <- c("CPO2", "CP03", "Clinical Group")
df <- as.data.frame(df)
df[,1] <- CPO2
df[,2] <- CPO3
df[,3] <- "ASD"
df[demDat$Group==2,3] <-"TD"

p3 <- ggplot(data = df, aes(x = `CPO3`, y = `CPO2`, color = `Clinical Group`)) + geom_abline(slope=1,intercept = 0) +
  geom_point() + theme_bw() + xlab(TeX("$M_{3}\\;\\;\\log{(CPO)}$")) + ylab(TeX("$M_{2}\\;\\;\\log{(CPO)}$"))

grid.arrange(p1, p2, ncol =2)

sum(CPO)
sum(CPO1)
sum(CPO2)
sum(CPO3)
