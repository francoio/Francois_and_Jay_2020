
library(tfa)


#####

## Read data
## simu ms_prime DIVERGENCE model (1500 generations) 
Y <- read.csv2("genotype_msprime_divergence_model_divtime_1500g.csv", header = FALSE)
ages <- Y[,1]
pop <- Y[,2]
Y <- as.matrix(Y[,-(1:2)])

#################
cc = colorRampPalette(c("grey80", "gold", "palegreen", "lightblue", "darkblue"))
ind = round(10*ages/max(ages)) + 1
col = cc(12)[ind]

fa <- tfa(sample_ages = ages, k = 2, Y = Y, lambda = 1e-2)
Un <- fa$u

pc <- prcomp(Y) 

par(mfrow = c(1,2))
plot(pc$x, 
     cex = 1.5, 
     col = col, 
     pch = 19, 
     xlim = c(-17,17),
     ylim = c(-22,19), 
     main = "PCA",
     las = 1)

plot(Un, col = col, pch = 19, cex = 1.5, main = "FA", 
     xlab = "Factor 1", 
     ylab = "Factor 2",
     las = 1)

