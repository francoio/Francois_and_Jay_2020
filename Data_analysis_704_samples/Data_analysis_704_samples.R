#####
# Global analysis of ancient DNA samples (704 samples)
#####

library(data.table)
library(LEA)
library(tfa)
library(lfmm)

## Download filtered genotypes and metadata for ancient and present-day samples

file_url <- "http://membres-timc.imag.fr/Olivier.Francois/Francois_and_Jay_2020/Ancient_DNA.zip"
download.file(url = file_url)

## Loading genotypes for ancient samples (filtered data)
geno_filt <- fread("./geno/geno_filter.lfmm", header = FALSE)
geno_filt <- as.matrix(geno_filt)
meta_filt <- read.csv("./geno/meta_filter.csv")
age_filt <- as.matrix(read.table("./geno/age_filter.txt"))[,1]

## Removing non Eurasian samples (YRI) 
meta <- meta_filt[meta_filt$Group.ID != "YRI.SG",]
geno <- geno_filt[meta_filt$Group.ID != "YRI.SG",]
age <- age_filt[meta_filt$Group.ID != "YRI.SG"]


cc = colorRampPalette(c("grey", "yellow","orange", "red", "brown"))
ind = round(30*age/max(age)) + 1
col =  cc(31)[ind]
plot(age, pch = 19, cex = 1.2, col = col)
 


############## FA of ancient samples ############################################

library(tfa)


plot(table(meta_filt[age_filt > 0, ]$Country), 
     las = 3, 
     ylab = "Frequency",
     cex.axis = .3)

# correction for coverage
coverage = as.numeric(as.character(meta_filt[age_filt > 0, ]$Coverage))

geno_ancient_cor <- coverage_adjust(geno_filt[age_filt > 0,], 
                                    coverage, K = 13, log = TRUE)

geno_ancient <- geno_ancient_cor[log10(coverage) > -0.1,]
age_ancient <- age_filt[age_filt > 0][log10(coverage) > -0.1]
meta_ancient <-  meta_filt[age_filt > 0,][log10(coverage) > -0.1,]


# FA analysis
fa_ancient = tfa(age_ancient, 
                geno_ancient, 
                lambda = 1e-1,
                k = 4,
                center = TRUE,
                coverage = coverage[log10(coverage) > -0.1],
                log = TRUE)

# sign change for a representation that agree with geography
fa_ancient$u[,1] <- - fa_ancient$u[,1]


# Build figure
steppe <- c("Russia_Yamnaya_Samara",
            "Russia_Yamnaya_Samara_published",
            "Russia_Yamnaya_Kalmykia.SG",
            "Russia_Poltavka")

plot(fa_ancient$u[,2:1],
     pch = 19, cex = .4, col = "grey80",
     xlab = "Factor 2", ylab = "Factor 1",
     xlim = c(-60, 70),
     ylim = c(-50, 90),
     las = 1)

points(fa_ancient$u[meta_ancient$Group.ID == "Israel_C",2:1], 
       pch = 8, cex = .6, col = "salmon")

m_anatolia_eba <- apply(fa_ancient$u[meta_ancient$Group.ID %in% 
                          c("Anatolia_EBA", 
                            "Anatolia_EBA.SG",
                            "Anatolia_C"),2:1], 2, mean)

points(matrix(m_anatolia_eba, nrow = 1), 
       pch = 21, cex = 1, lwd = 2, col = "salmon")

points(fa_ancient$u[meta_ancient$Country == "Armenia",2:1], 
       pch = 8, cex = .6, col = "pink")

points(fa_ancient$u[meta_ancient$Group.ID == "Iran_Seh_Gabi_C",2:1], 
       pch = 8, cex = .6, col = "pink")

m_iran_n <- apply(fa_ancient$u[meta_ancient$Group.ID %in% 
                                       c("Iran_Tepe_Abdul_Hosein_N.SG", 
                                         "Iran_Wezmeh_N.SG",
                                         "Iran_Ganj_Dareh_N"),2:1], 2, mean)

points(matrix(m_iran_n, nrow = 1), 
       pch = 21, cex = 1, lwd = 2, col = "pink")

points(fa_ancient$u[meta_ancient$Group.ID %in% steppe,2:1],
       pch = 8, cex = .6, col = "darkblue")

points(fa_ancient$u[meta_ancient$Group.ID == "Ukraine_N",2:1],
       pch = 2, cex = .6, col = "green4")

points(fa_ancient$u[meta_ancient$Group.ID == "Latvia_HG",2:1], 
       pch = 2, cex = .6, col = "olivedrab")

m_sweden_hg <- apply(fa_ancient$u[meta_ancient$Group.ID %in% 
                               c("Sweden_Motala_HG"),2:1], 2, mean)

points(matrix(m_sweden_hg, nrow = 1), 
       pch = 21, cex = 1, lwd = 2, col = "green4")

#points(fa_ancient$u[meta_ancient$Group.ID == "Sweden_Motala_HG",2:1], 
#       pch = 21, cex = .7, col = "green4")

points(fa_ancient$u[meta_ancient$Group.ID == "Serbia_Iron_Gates_HG",2:1], 
       pch = 8, cex = .6, col = "olivedrab")

m_gb_n <- apply(fa_ancient$u[meta_ancient$Group.ID %in% 
                                 c("Scotland_N", 
                                   "England_N"),2:1], 2, mean)

points(matrix(m_gb_n, nrow = 1), 
       pch = 21, cex = 1.2, lwd = 2, col = "salmon2")

m_germany_n <- apply(fa_ancient$u[meta_ancient$Group.ID %in% 
                               c("Germany_LBK_EN", 
                                 "Germany_Esperstedt_MN",
                                 "Germany_Blatterhohle_MN",
                                 "Germany_LBK_EN_Stuttgart_published.DG"),2:1], 2, mean)

points(matrix(m_germany_n, nrow = 1), 
       pch = 21, cex = 1.2, lwd = 2, col = "salmon2")

#points(fa_ancient$u[meta_ancient$Group.ID == "Germany_LBK_EN",2:1], 
#       pch = 19, cex = .6, col = "salmon1")

m_hungary_n <- apply(fa_ancient$u[meta_ancient$Group.ID %in% 
                                    c("Hungary_LBK_MN", 
                                      "Hungary_ALPc_MN",
                                      "Hungary_LBK_MN.SG",
                                      "Hungary_ALPc_Tiszadob_MN",
                                      "Hungary_ALPc_Szatmar_MN",
                                      "Hungary_Starcevo_EN_all",
                                      "Hungary_Vinca_MN",
                                      "Hungary_Koros_EN"),2:1], 2, mean)

points(matrix(m_hungary_n, nrow = 1), 
       pch = 0, cex = 1, lwd = 2, col = "salmon")

m_czech_n <- apply(fa_ancient$u[meta_ancient$Group.ID %in% 
                                    c("Czech_MN", 
                                      "Czech_N"),2:1], 2, mean)

points(matrix(m_czech_n, nrow = 1), 
       pch = 0, cex = 1, lwd = 2, col = "salmon")

points(fa_ancient$u[meta_ancient$Group.ID %in% c("Anatolia_N", "Anatolia_N.SG"),2:1], 
       pch = 8, cex = .6, col = "salmon3")


m_gb_bb <- apply(fa_ancient$u[meta_ancient$Group.ID %in% 
                                  c("England_Bell_Beaker"),2:1], 2, mean)

points(matrix(m_gb_bb, nrow = 1), 
       pch = 21, cex = 1.2, lwd = 2, col = "yellow4")

m_germany_bb <- apply(fa_ancient$u[meta_ancient$Group.ID %in% 
                                c("Germany_Bell_Beaker"),2:1], 2, mean)

points(matrix(m_germany_bb, nrow = 1), 
       pch = 21, cex = 1.2, lwd = 2, col = "yellow4")

m_hungary_bb <- apply(fa_ancient$u[meta_ancient$Group.ID %in% 
                                     c("Hungary_Bell_Beaker",
                                       "Hungary_Bell_Beaker_EBA"),2:1], 2, mean)

points(matrix(m_hungary_bb, nrow = 1), 
       pch = 0, cex = 1, lwd = 2, col = "yellow4")

m_hungary_lg <- apply(fa_ancient$u[meta_ancient$Group.ID %in% 
                                     c("Hungary_Langobard.SG",
                                       "Hungary_Langobard"),2:1], 2, mean)
#points(matrix(m_hungary_lg, nrow = 1), 
#       pch = 5, cex = .8, lwd = 2, col = "yellow4")


m_czech_bb <- apply(fa_ancient$u[meta_ancient$Group.ID %in% 
                                     c("Czech_Bell_Beaker"),2:1], 2, mean)

#points(matrix(m_czech_bb, nrow = 1), 
#       pch = 0, cex = 1, lwd = 2, col = "yellow4")

#points(fa_ancient$u[meta_ancient$Group.ID == "England_Bell_Beaker",2:1], 
#     pch = 19, cex = .6, col = "yellow3")


m_sweden_vi <- apply(fa_ancient$u[meta_ancient$Group.ID %in% 
                                   c("Sweden_Viking.SG"),2:1], 2, mean)

points(matrix(m_sweden_vi, nrow = 1), 
       pch = 5, cex = .8, lwd = 2, col = "green4")

#points(fa_ancient$u[meta_ancient$Group.ID == "Sweden_Viking.SG",2:1], 
#       pch = 19, cex = .6, col = "olivedrab")





############ PC projections #################

### Projections of ancient samples on PCs of 
### present-day populations


# computing pca of present-day samples 
pc0 = prcomp(geno[age==0,], scale = FALSE)
c_mean <- apply(geno[age==0,], 2, mean) 
G <- t(t(geno) - c_mean)

# projections of ancient samples on 
who <- rep(TRUE, 1225) ## all samples
Wn <- G[who,][age[who] > 0,] %*% pc0$rotation

# sign change for comparison with FA 
Wn[,1] = - Wn[,1]

# Generating Figure 
plot(Wn[,1:2], pch = 19, cex = .3, 
     col = "grey80",
     xlab = "Axis 1", ylab = "Axis 2",
     xlim = c(-30, 22),
     ylim = c(-35, 43),    
     las = 1)

W_israel <- G[who,][meta[who,]$Group.ID == "Israel_C",] %*% pc0$rotation
W_israel[,1] = - W_israel[,1]
points(W_israel[,], pch = 8, cex = .6, col = "salmon3")

W_anatolia_eba <- G[who,][meta[who,]$Group.ID %in% 
                          c("Anatolia_EBA", 
                            "Anatolia_EBA.SG",
                            "Anatolia_C"),] %*% pc0$rotation

W_anatolia_eba[,1] = - W_anatolia_eba[,1]
W_anatolia_eba <- W_anatolia_eba[,]

m_anatolia_eba <- apply(W_anatolia_eba, 2, mean)
points(matrix(m_anatolia_eba, nrow = 1), 
       pch = 21, cex = 1, lwd = 2, col = "salmon")

W_armenia <- G[who,][meta[who,]$Country == "Armenia",] %*% pc0$rotation
W_armenia[,1] = - W_armenia[,1]
points(W_armenia[,], 
       pch = 8, cex = .6, col = "pink")

W_iran_c <- G[who,][meta[who,]$Group.ID == "Iran_Seh_Gabi_C",] %*% pc0$rotation
W_iran_c[,1] = - W_iran_c[,1]
points(W_iran_c[,], pch = 8, cex = .6, col = "pink")

m_anatolia_eba <- apply(W_anatolia_eba[,], 2, mean)
points(matrix(m_anatolia_eba, nrow = 1), 
       pch = 21, cex = 1, lwd = 2, col = "salmon")

W_iran_n <- G[who,][meta[who,]$Group.ID  %in% 
                      c("Iran_Tepe_Abdul_Hosein_N.SG", 
                        "Iran_Wezmeh_N.SG",
                        "Iran_Ganj_Dareh_N"),] %*% pc0$rotation
W_iran_n[,1] = - W_iran_n[,1]
W_iran_n <- W_iran_n[,]
m_iran_n <- apply(W_iran_n, 2, mean)
points(matrix(m_iran_n, nrow = 1), 
       pch = 21, cex = 1, lwd = 2, col = "pink")

steppe <- c("Russia_Yamnaya_Samara",
            "Russia_Yamnaya_Samara_published",
            "Russia_Yamnaya_Kalmykia.SG",
            "Russia_Poltavka")
W_steppe <- G[who,][meta[who,]$Group.ID %in% steppe,] %*% pc0$rotation
W_steppe[,1] = - W_steppe[,1]
W_steppe = W_steppe[,]
points(W_steppe, pch = 8, cex = .6, col = "darkblue")

W_ukraine_n <- G[who,][meta[who,]$Group.ID == "Ukraine_N",] %*% pc0$rotation
W_ukraine_n[,1] = - W_ukraine_n[,1]
W_ukraine_n = W_ukraine_n[,]
points(W_ukraine_n, pch = 8, cex = .6, col = "green4")

W_latvia_hg <- G[who,][meta[who,]$Group.ID == "Latvia_HG",] %*% pc0$rotation
W_latvia_hg[,1] = - W_latvia_hg[,1]
W_latvia_hg = W_latvia_hg[,]
points(W_latvia_hg, pch = 19, cex = .6, col = "olivedrab")

W_serbia_hg <- G[who,][meta[who,]$Group.ID == "Serbia_Iron_Gates_HG",] %*% pc0$rotation
W_serbia_hg[,1] = - W_serbia_hg[,1]
W_serbia_hg = W_serbia_hg[,]
points(W_serbia_hg, pch = 8, cex = .6, col = "olivedrab")

W_sweden_hg <- G[who,][meta[who,]$Group.ID   %in% 
                         c("Sweden_Motala_HG"),] %*% pc0$rotation
W_sweden_hg[,1] = - W_sweden_hg[,1]
W_sweden_hg <- W_sweden_hg[,]
m_sweden_hg <- apply(W_sweden_hg, 2, mean)
points(matrix(m_sweden_hg, nrow = 1), 
       pch = 21, cex = 1, lwd = 2, col = "green4")

W_gb_n <- G[who,][meta[who,]$Group.ID   %in% 
                         c("Scotland_N", 
                           "England_N"),] %*% pc0$rotation
W_gb_n[,1] = - W_gb_n[,1]
W_gb_n <- W_gb_n[, ]
m_gb_n <- apply(W_gb_n, 2, mean)
points(matrix(m_gb_n, nrow = 1), 
       pch = 21, cex = 1.2, lwd = 2, col = "salmon2")

W_germany_n <- G[who,][meta[who,]$Group.ID   %in% 
                         c("Germany_LBK_EN", 
                           "Germany_Esperstedt_MN",
                           "Germany_Blatterhohle_MN",
                           "Germany_LBK_EN_Stuttgart_published.DG"),] %*% pc0$rotation
W_germany_n[,1] = - W_germany_n[,1]
W_germany_n <- W_germany_n[,]
m_germany_n<- apply(W_germany_n, 2, mean)
points(matrix(m_germany_n, nrow = 1), 
       pch = 21, cex = 1.2, lwd = 2, col = "salmon2")

W_hungary_n <- G[who,][meta[who,]$Group.ID   %in% 
                         c("Hungary_LBK_MN", 
                           "Hungary_ALPc_MN",
                           "Hungary_LBK_MN.SG",
                           "Hungary_ALPc_Tiszadob_MN",
                           "Hungary_ALPc_Szatmar_MN",
                           "Hungary_Starcevo_EN_all",
                           "Hungary_Vinca_MN",
                           "Hungary_Koros_EN"),] %*% pc0$rotation
W_hungary_n[,1] = - W_hungary_n[,1]
W_hungary_n <- W_hungary_n[,]
m_hungary_n <- apply(W_hungary_n, 2, mean)
points(matrix(m_hungary_n, nrow = 1), 
       pch = 0, cex = 1, lwd = 2, col = "salmon")

W_czech_n <- G[who,][meta[who,]$Group.ID   %in% 
                       c("Czech_MN", 
                         "Czech_N"),] %*% pc0$rotation
W_czech_n[,1] = - W_czech_n[,1]
W_czech_n <- W_czech_n[,]
m_czech_n <- apply(W_czech_n, 2, mean)
points(matrix(m_czech_n, nrow = 1), 
       pch = 0, cex = 1, lwd = 2, col = "salmon")

W_anatolia_n <- G[who,][meta[who,]$Group.ID %in% 
                          c("Anatolia_N", "Anatolia_N.SG"),] %*% pc0$rotation
W_anatolia_n[,1] = - W_anatolia_n[,1]
W_anatolia_n = W_anatolia_n[,]
points(W_anatolia_n, 
       pch = 8, cex = .6, col = "salmon3")

W_gb_bb <- G[who,][meta[who,]$Group.ID   %in% 
                     c("England_Bell_Beaker"),] %*% pc0$rotation
W_gb_bb[,1] = - W_gb_bb[,1]
W_gb_bb <- W_gb_bb[,]
m_gb_bb <- apply(W_gb_bb, 2, mean)
points(matrix(m_gb_bb, nrow = 1), 
         pch = 21, cex = 1.2, lwd = 2, col = "yellow4")

W_germany_bb <- G[who,][meta[who,]$Group.ID   %in% 
                     c("Germany_Bell_Beaker"),] %*% pc0$rotation
W_germany_bb[,1] = - W_germany_bb[,1]
W_germany_bb <- W_germany_bb[,]
m_germany_bb <- apply(W_germany_bb, 2, mean)
points(matrix(m_germany_bb, nrow = 1), 
       pch = 21, cex = 1.2, lwd = 2, col = "yellow4")

W_hungary_bb <- G[who,][meta[who,]$Group.ID   %in% 
                          c("Hungary_Bell_Beaker",
                            "Hungary_Bell_Beaker_EBA"),] %*% pc0$rotation
W_hungary_bb[,1] = - W_hungary_bb[,1]
W_hungary_bb <- W_hungary_bb[,]
m_hungary_bb <- apply(W_hungary_bb, 2, mean)
points(matrix(m_hungary_bb, nrow = 1), 
       pch = 0, cex = 1, lwd = 2, col = "yellow4")

W_hungary_lg <- G[who,][meta[who,]$Group.ID   %in% 
                          c("Hungary_Langobard.SG",
                            "Hungary_Langobard"),] %*% pc0$rotation
W_hungary_lg[,1] = - W_hungary_lg[,1]
W_hungary_lg <- W_hungary_lg[,]
#m_hungary_lg <- apply(W_hungary_lg, 2, mean)
#points(matrix(m_hungary_lg, nrow = 1), 
#       pch = 5, cex = .8, lwd = 2, col = "yellow4")

W_czech_bb <- G[who,][meta[who,]$Group.ID   %in% 
                          c("Czech_Bell_Beaker"),] %*% pc0$rotation
W_czech_bb[,1] = - W_czech_bb[,1]
W_czech_bb <- W_czech_bb[,]
m_czech_bb <- apply(W_czech_bb, 2, mean)
points(matrix(m_czech_bb, nrow = 1), 
       pch = 0, cex = 1, lwd = 2, col = "yellow4")

W_sweden_vi <- G[who,][meta[who,]$Group.ID   %in% 
                         c("Sweden_Viking.SG"),] %*% pc0$rotation
W_sweden_vi[,1] = - W_sweden_vi[,1]
W_sweden_vi <- W_sweden_vi[,]
m_sweden_vi <- apply(W_sweden_vi, 2, mean)
points(matrix(m_sweden_vi, nrow = 1), 
       pch = 5, cex = .8, lwd = 2, col = "green4")

 


