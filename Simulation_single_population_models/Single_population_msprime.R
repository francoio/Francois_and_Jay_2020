library(tfa)

##
# rm(list = ls())

## read age and genotypes
Y <- read.csv2("genotype_msprime_one_population.csv", header = F)
ages <- Y[,1]
Y <- as.matrix(Y[,-(1:2)])

# color palette
ind = round(40*ages/max(ages)) + 1
cc = colorRampPalette(c("grey", "gold", "palegreen", "lightblue", "darkblue"))
col = cc(41)[ind]

# factor analysis
fa <- tfa(sample_ages = ages, Y = Y, k = 2, lambda = 1e-3)
Un <- fa$u

# pca
pc <- prcomp(Y, scale = FALSE) 

### covariance models

# conversion of ages in (0,1)
range_ages <- max(ages) - min(ages)
tn <- 1 - ages/range_ages

# normalisation of times
var_Y <- apply(Y, 1, FUN = var)
tn <- min(var_Y) +  (max(var_Y) - min(var_Y))*tn

# Brownian covariance model
n <- length(tn)
C <- matrix(NA, n, n)

for (i in 1:n){
  for (j in 1:n)  C[i,j] <- min(tn[i], tn[j])}
diag(C) <- NA

# empirical covariance
C_obs <- cov(t(Y))
diag(C_obs) <- NA


#### FIGURE
par(mfrow = c(2,2))

image(x = 1:41, y = -(41:1) , z = C_obs[,41:1], 
      col = rev(col), las = 1, main = "Covariance", 
      xlab = "Samples",
      ylab = "Samples")

plot(pc$x, cex = 1.5, col = rev(col), pch = 19, 
     main = "PCA", las = 1)

image(x = 1:41, y = -(41:1) , 
      z= C[,41:1], col = rev(col), las = 1,
      main = "Drift model", 
      xlab = "Samples", ylab = "Samples")

plot(Un, col = rev(col), pch = 19, cex = 1.5, main = "FA", 
     xlab = "Factor 1", 
     ylab = "Factor 2", las = 1)



