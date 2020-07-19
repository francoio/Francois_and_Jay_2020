### 2-way admixture analysis for Great Britain samples
library(tfa)
library(LEA)
library(admixr)

trunc <- function(x){min(max(0,x),1)}

geno_ancient <- geno_filt[age_filt >0,]
meta_ancient <- meta_filt[age_filt >0,]
age_ancient  <- age_filt[age_filt >0]
coverage_ancient = as.numeric(as.character(meta_ancient$Coverage))

lst_keep <- which( (meta_ancient$Group.ID ==  "England_Bell_Beaker" |
                    meta_ancient$Group.ID == "England_MBA" |
                    meta_ancient$Group.ID == "England_CA_EBA" |
                    meta_ancient$Group.ID == "England_N" | 
                    meta_ancient$Group.ID =="Scotland_N" |
                    meta_ancient$Group.ID == "Russia_Yamnaya_Samara" |
                    meta_ancient$Group.ID == "Russia_Yamnaya_Samara_published" |
                    meta_ancient$Group.ID ==  "Russia_Yamnaya_Kalmykia.SG"|
                    meta_ancient$Group.ID ==  "Russia_Poltavka"|
                    meta_ancient$Group.ID ==  "Russia_Srubnaya_published"|
                    meta_ancient$Group.ID ==  "Russia_Sintashta_MLBA.SG"|
                    meta_ancient$Group.ID ==  "Anatolia_N"|
                    meta_ancient$Group.ID ==  "Anatolia_N.SG" )  & coverage > 0.7
)

geno_ancient <- geno_ancient_cor[lst_keep, ]
meta_ancient <- meta_ancient[lst_keep, ]
age_ancient <- age_ancient[lst_keep]
coverage_ancient <- coverage_ancient[lst_keep] 

plot(table(meta_ancient$Country), las = 3, cex.axis = .5)

r_squared <- NULL
llambda = c(0, -.30103, -1, -1.30103, -1.69897, -2, -2.30103,
            -3, -3.30103,-4)
lambda <- 10^llambda

for (l in lambda){
fa_rm = tfa::tfa(age_ancient, 
                    geno_ancient, 
                    lambda = l, 
                    center = TRUE,
                    coverage = coverage_ancient,
                    log = FALSE)

mod <- lm(age_ancient ~ fa_rm$u)
r_squared <- c(r_squared, summary(mod)$adj.r.squared)
print(r_squared)
}

r_sq_gb <- r_squared

plot(llambda, 
     r_squared, 
     xlab = "log10(lambda)",
     ylab = "Variance explained by factors",
     main = "Great Britain", 
     las = 1)

points(llambda, 
       r_squared, 
     col = "lightblue", 
     pch = 19,
     cex = .8)

abline(v = log10(5e-3), col = "orange", lty = 2)

### FA analysis 
fa <-  tfa(age_ancient, 
            geno_ancient, 
            k = 2,
            lambda = 5e-3, #4e-3,                     
            center = TRUE,
            coverage = coverage_ancient,
            log = TRUE)

plot(fa$u, 
     pch = 19, cex = 2, col = "grey90",
     xlab = "Factor 1", ylab = "Factor 2", 
     main = "FA",
     ylim = c(-10, 15), 
     las = 1)


points(fa$u[meta_ancient$Group.ID == "England_Bell_Beaker",], 
       pch = 19, cex = 1, col = "yellow3")

points(fa$u[meta_ancient$Group.ID %in% c("England_N","Scotland_N"),], 
       pch = 19, cex = 1, col = "salmon1")

points(fa$u[meta_ancient$Country == "Russia",], 
       pch = 8, cex = .6, col = "darkblue")

points(fa$u[meta_ancient$Country == "Turkey",], 
       pch = 8, cex = .6, col = "salmon3")

points(fa_rm$u[meta_ancient$Group.ID ==  "England_MBA",], 
       pch = 8, cex = .6, col = "yellow4")
abline(h = 0, lty = 2, col = "orange")


### Evaluation of ancestry 
ne = mean(fa$u[meta_ancient$Country == "Turkey",1])
gba = mean(fa$u[meta_ancient$Group.ID == "England_Bell_Beaker",1])
r = mean(fa$u[meta_ancient$Country == "Russia",1])

cat("Yamnaya ancestry in England Bell Beaker from FA",
(ne - gba)/(ne - r), "\n") 


ancestry_gb <- (ne - fa$u[meta_ancient$Group.ID == "England_Bell_Beaker"
                          ,1])/(ne - r)
ancestry_gb <- sapply(ancestry_gb, trunc)

barplot(rbind(ancestry_gb, 1 - ancestry_gb), 
        space = 0, border = NA,
        col = c("darkblue", "salmon3"),
        ylab = "Proportion",
        xlab = "Individuals",
        main = "Yamnaya ancestry")

abline(h = 0.52, col = "yellow", lty = 2, lwd = 2)

axis(1, 
     at = 1:length(ancestry_gb), 
     labels = as.character(meta_ancient[meta_ancient$Group.ID == "England_Bell_Beaker",]$Instance.ID),
     las = 3, 
     cex.axis = .4)



###### Three population models ###############@
### Admixture GB + HG
geno_ancient <- geno_filt[age_filt >0,]
meta_ancient <- meta_filt[age_filt >0,]
age_ancient  <- age_filt[age_filt >0]
coverage_ancient = as.numeric(as.character(meta_ancient$Coverage))


lst_keep <- which( (meta_ancient$Group.ID ==  "England_Bell_Beaker" |
                      meta_ancient$Group.ID == "England_MBA" |
                      meta_ancient$Group.ID == "England_CA_EBA" |
                      meta_ancient$Group.ID == "England_N" | 
                      meta_ancient$Group.ID =="Scotland_N" |
                      meta_ancient$Group.ID == "Russia_Yamnaya_Samara" |
                      meta_ancient$Group.ID == "Russia_Yamnaya_Samara_published" |
                      meta_ancient$Group.ID ==  "Russia_Yamnaya_Kalmykia.SG"|
                      meta_ancient$Group.ID ==  "Russia_Poltavka"|
                      meta_ancient$Group.ID ==  "Serbia_Iron_Gates_HG"|
                    #meta_ancient$Group.ID ==  "Russia_Sintashta_MLBA.SG"|
                      meta_ancient$Group.ID ==  "Anatolia_N"|
                      meta_ancient$Group.ID ==  "Anatolia_N.SG" )  & coverage > 0.7
)


geno_ancient <- geno_ancient[lst_keep, ]


coverage_ancient <- coverage_ancient[lst_keep] 
geno_ancient <- geno_ancient_cor[lst_keep, ]

meta_ancient <- meta_ancient[lst_keep, ]
age_ancient <- age_ancient[lst_keep]


fa = tfa(age_ancient, 
            geno_ancient, 
            k = 4,
            lambda = 5e-1,#1e-3,                     
            center = TRUE,
            coverage = coverage_ancient,
            log = TRUE)


plot(fa$u, pch = 19, cex = 2, col = "grey90",
     xlab = "Factor 1", ylab = "Factor 2", main = "FA")
     #xlim = c(-70, 100), ylim = c(-80, 55))


m_yamnaya <- apply(fa$u[meta_ancient$Country == "Russia",],
                   2, mean)

points(matrix(m_yamnaya, nrow = 1), 
       pch = 21, cex = 1, lwd = 2, col = "black")

m_anatolia <- apply(fa$u[meta_ancient$Country == "Turkey",],
                    2, mean)

points(matrix(m_anatolia, nrow = 1), 
       pch = 21, cex = 1, lwd = 2, col = "black")

lines(rbind(m_anatolia,m_yamnaya))
m_hg <- apply(fa$u[meta_ancient$Group.ID == "Serbia_Iron_Gates_HG",],
              2, mean)

points(matrix(m_hg, nrow = 1), 
       pch = 21, cex = 1, lwd = 2, col = "black")

lines(rbind(m_hg,m_yamnaya))
lines(rbind(m_anatolia,m_hg))

points(fa$u[meta_ancient$Country == "Great Britain" & age_ancient < 4350,], 
       pch = 19, cex = 1, col = "yellow3")

m_gb_ba <- apply(fa$u[meta_ancient$Country == "Great Britain" & age_ancient < 4350,],
              2, mean)
points(t(as.matrix(m_gb_ba)))

m_gb_bb <- apply(fa$u[meta_ancient$Group.ID ==  "England_Bell_Beaker",],
                 2, mean)
points(t(as.matrix(m_gb_bb)))

points(fa$u[meta_ancient$Country == "Great Britain" & age_ancient > 4350,], 
       pch = 19, cex = 1, col = "salmon1")

points(fa$u[meta_ancient$Country == "Russia",], pch = 8, cex = .6, col = "darkblue")

points(fa$u[meta_ancient$Country == "Turkey",], pch = 8, cex = .6, col = "salmon3")

points(fa$u[meta_ancient$Group.ID ==  "England_MBA",], 
       pch = 8, cex = .6, col = "yellow4")

points(fa$u[meta_ancient$Group.ID ==  "Serbia_Iron_Gates_HG",], 
       pch = 8, cex = .6, col = "olivedrab")

abline(h = 0, lty = 2, col = "orange")







#### STRUCTURE analysis with LEA::snmf

write.geno(R = geno_filt[age_filt > 0,][lst_keep,], 
           output.file = "target.geno")

proj <- snmf("target.geno", 
             K = 3, 
             entropy = TRUE, 
             alpha = 1000, 
             project = "new")

barchart(proj, K = 3, 
         space = 0, border = NA,
         col = c("darkblue","salmon3","green4"), 
         xlab = "", 
         ylab = "Ancestry proportions", 
         main = "Ancestry matrix") -> bp

axis(1, at = 1:length(bp$order), 
     labels = as.character(meta_ancient$Group.ID)[bp$order], 
     las = 3, cex.axis = .4)

Qm <- Q(proj, K = 3)[meta_ancient$Country %in% c("Great Britain"),]
## Check Yamnaya ancestry from the barplot (max or min)
cat("Yamnaya ancestry from STRUCTURE =", 
    round(apply(Qm,2,mean),3),"\n")


## clean environment
remove.snmfProject("target.snmfProject")















### ADMIXTOOLS ANALYSIS
library(admixr)


snps <- eigenstrat("./snps/snps")


result <- qpAdm(
  target = c("England_Bell_Beaker", "England_MBA", "England_CA_EBA", "England_N", "Scotland_N"),
  sources = c("Russia_Yamnaya_Samara",  "Anatolia_N"),
  outgroups = c("YRI.SG", "Russia_Sidelkino_HG.SG", "France_Ranchot88_published"),
  data = snps
)

View(result$proportions)
