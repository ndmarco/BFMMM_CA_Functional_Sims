### Real Case study
library(BayesFMMM)
library(fda)
library(refund)
setwd(".")

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
trapz(out1, pa.dat$y[1:33]) # all functional observations integrate to 1 (normalized across electordes, subjects)

### Convert to wide format
Y <- pa.dat
## paper used T8 electrode
Y <- Y[Y$reg == 15,]
Y$ID <- paste(Y$ID, Y$reg, sep = ".")
Y <- reshape(Y[,c(1,3,6)], idvar = "ID", timevar = "func", direction = "wide")
Y <- Y[,-1]
Y <- as.matrix(Y)
X <- matrix((demDat$Age), ncol = 1)
time <- seq(6, 14, 0.25)

daybasis20 <- create.bspline.basis(rangeval=c(6,14), nbasis = 20)
Temp.fd <- smooth.basisPar(time, t(Y), daybasis20)$fd
modmat <- cbind(1, X)
olsmod <- fosr(fdobj = Temp.fd, X = modmat)
fd <- olsmod$fd
basis <- getbasismatrix(time, daybasis20)
plot(basis %*% fd$coefs[,1])

df <- as.data.frame(matrix(0, nrow = 33*length(X_seq), ncol = 3))
colnames(df) <- c("Y", "t", "Age (Months)")
X_seq <- (seq(25,145, 10))

for(i in 1:length(X_seq)){
  df[((i-1)*33 + 1):(i*33),1] <- (basis %*% fd$coefs[,1]) + (X_seq[i] * (basis %*% fd$coefs[,2]))
  df[((i-1)*33 + 1):(i*33),2] <- time
  df[((i-1)*33 + 1):(i*33),3] <- X_seq[i]
}

library(ggplot2)
p1 <- ggplot(df, aes(x = t, y =  Y, color = `Age (Months)`, group = `Age (Months)`)) +
  geom_line() + theme_bw() + ggtitle("Function-on-Scalar Regression") + ylab("Relative Power") + xlab("Frequency (Hz)") + ylim(c(0,0.5))


df4 <- as.data.frame(matrix(0, nrow = 33*30, ncol = 4))
colnames(df4) <- c("Y", "t", "Age (Months)", "id")
for(i in 1:30){
  df4[((i-1)*33 + 1):(i*33),1] <- Y[i,]
  df4[((i-1)*33 + 1):(i*33),2] <- time
  df4[((i-1)*33 + 1):(i*33),3] <- X[i]
  df4[((i-1)*33 + 1):(i*33),4] <- i
}
p4 <- ggplot(df4, aes(x = t, y =  Y, color = `Age (Months)`, group = `id`)) +
  geom_line() + theme_bw() + ggtitle("Raw Data") + ylab("Relative Power") + xlab("Frequency (Hz)") + ylim(c(0,0.5))

grid.arrange(p4, p1, nrow = 1)

daybasis20 <- create.bspline.basis(rangeval=c(6,14), nbasis = 20)
Temp.fd <- smooth.basisPar(time, t(Y), daybasis20)$fd
ASD_bin <- rep(0,97)
ASD_bin[demDat$Group == 2] <- 1
modmat <- cbind(1, X, ASD_bin, ASD_bin * X)
olsmod <- fosr(fdobj = Temp.fd, X = modmat)
fd <- olsmod$fd
time <- seq(6,14, 0.05)
basis <- getbasismatrix(time, daybasis20)
plot(basis %*% fd$coefs[,1])

df <- as.data.frame(matrix(0, nrow = 161*length(X_seq), ncol = 3))
colnames(df) <- c("Y", "t", "Age (Months)")
X_seq <- (seq(25,145, 10))
for(i in 1:length(X_seq)){
  df[((i-1)*161 + 1):(i*161),1] <- (basis %*% fd$coefs[,1]) + (X_seq[i] * (basis %*% fd$coefs[,2]))
  df[((i-1)*161 + 1):(i*161),2] <- time
  df[((i-1)*161 + 1):(i*161),3] <- X_seq[i]
}

library(ggplot2)
p2 <- ggplot(df, aes(x = t, y =  Y, color = `Age (Months)`, group = `Age (Months)`)) +
  geom_line() + theme_bw() + ggtitle("Function-on-Scalar Regression (TD)") + ylab("Relative Power") + xlab("Frequency (Hz)") + ylim(c(0,0.3))


df <- as.data.frame(matrix(0, nrow = 161*length(X_seq), ncol = 3))
colnames(df) <- c("Y", "t", "Age (Months)")
X_seq <- (seq(25,145, 10))
for(i in 1:length(X_seq)){
  df[((i-1)*161 + 1):(i*161),1] <- (basis %*% fd$coefs[,1]) + (X_seq[i] * (basis %*% fd$coefs[,2])) + (basis %*% fd$coefs[,3]) + (X_seq[i] * (basis %*% fd$coefs[,4]))
  df[((i-1)*161 + 1):(i*161),2] <- time
  df[((i-1)*161 + 1):(i*161),3] <- X_seq[i]
}

library(ggplot2)
p3 <- ggplot(df, aes(x = t, y =  Y, color = `Age (Months)`, group = `Age (Months)`)) +
  geom_line() + theme_bw() + ggtitle("Function-on-Scalar Regression (ASD)") + ylab("Relative Power") + xlab("Frequency (Hz)") + ylim(c(0,0.3))

grid.arrange(p3, p2, nrow = 1)
