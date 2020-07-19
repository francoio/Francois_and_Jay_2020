### 2-way admixture analysis for Germany/Czech republic samples
library(tfa)
library(LEA)
library(admixr)

geno_ancient <- geno_filt[age_filt >0,]
meta_ancient <- meta_filt[age_filt >0,]
age_ancient  <- age_filt[age_filt >0]
coverage_ancient = as.numeric(as.character(meta_ancient$Coverage))

lst_keep <- which( (  meta_ancient$Group.ID ==  "Germany_Bell_Beaker" |
                      meta_ancient$Group.ID ==  "Germany_Corded_Ware" |
                      meta_ancient$Group.ID == "Czech_Bell_Beaker" |
                      meta_ancient$Group.ID == "Czech_Corded_Ware" |
                      meta_ancient$Group.ID == "Czech_EBA" |
                      meta_ancient$Group.ID == "Germany_Unetice_EBA" |
                      meta_ancient$Group.ID == "Russia_Yamnaya_Samara" |
                      meta_ancient$Group.ID == "Russia_Yamnaya_Samara_published" |
                      meta_ancient$Group.ID ==  "Russia_Yamnaya_Kalmykia.SG"|
                      meta_ancient$Group.ID ==  "Russia_Poltavka"|
                      meta_ancient$Group.ID ==  "Anatolia_N"|
                      meta_ancient$Group.ID ==  "Anatolia_N.SG" )  & coverage > 0.7
)

geno_ancient <- geno_ancient_cor[lst_keep, ]
meta_ancient <- meta_ancient[lst_keep, ]
age_ancient <- age_ancient[lst_keep]
coverage_ancient <- coverage_ancient[lst_keep] 
meta_ancient <- meta_ancient[,c("Instance.ID", 
                                "Group.ID", 
                                "Country", 
                                "Coverage")]

plot(table(meta_ancient$Country), las = 3, cex.axis = .5)

r_squared <- NULL
llambda = c(0, -.30103, -1, -1.30103, -1.69897, -2, -2.30103,
            -3, -3.30103,-4, -4.30103, -5, -5.30103 )
lambda <- 10^llambda

for (l in lambda){
  mod = tfa(age_ancient, 
            geno_ancient, 
            lambda = l, 
            center = TRUE,
            coverage = coverage_ancient,
            log = TRUE)
  
  mod_lm <- lm(age_ancient ~ mod$u)
  r_squared <- c(r_squared, summary(mod_lm)$adj.r.squared)
  cat("l = ", l, "\n")
  cat(r_squared, "\n")
}

r_sq_ger <- r_squared

plot(llambda, 
     r_sq_ger, 
     xlab = "log(lambda)",
     ylab = "Variance explained by factors",
     main = "Germany/Czech R.",
     las = 1)

points(llambda, 
       r_sq_ger, 
       col = "lightblue", 
       pch = 19,
       cex = .8)

abline(v = log10(1e-3), col = "orange", lty = 2)


mod  =  tfa(age_ancient, 
                   geno_ancient, 
                   k = 2,
                   lambda = 5e-2,                     
                   center = TRUE,
                   coverage = coverage_ancient,
                   log = TRUE)

plot(mod$u, pch = 19, cex = 2, col = "grey90",
     xlab = "Factor 1", ylab = "Factor 2" ,
     main = "Germany/Czech R.",
     xlim = c(-50,50), ylim = c(-16, 25), las = 1)

points(mod$u[meta_ancient$Country == "Germany",], 
       pch = 19, cex = 1, col = "yellow2")
points(mod$u[meta_ancient$Country == "Czech Republic",], 
       pch = 19, cex = 1, col = "yellow3")

points(mod$u[meta_ancient$Country == "Russia",], pch = 8, cex = .6, col = "darkblue")
points(mod$u[meta_ancient$Country == "Turkey",], pch = 8, cex = .6, col = "salmon3")

points(mod$u[meta_ancient$Group.ID ==  "Germany_Corded_Ware",], pch = 8, cex = .6, col = "yellow4")
points(mod$u[meta_ancient$Group.ID ==  "Czech_Corded_Ware",], pch = 8, cex = .6, col = "yellow4")

points(mod$u[meta_ancient$Group.ID == "Germany_Unetice_EBA",], pch = 2, 
       cex = .6, col = "yellow4")
points(mod$u[meta_ancient$Group.ID == "Czech_EBA",], pch = 2, 
       cex = .6, col = "yellow4")

abline(h = 0, lty = 2, col = "orange")

### Evaluation of ancestries 
ne = mean(mod$u[meta_ancient$Country == "Turkey",1])
ru = mean(mod$u[meta_ancient$Country == "Russia",1])

de = mean(mod$u[meta_ancient$Group.ID ==  "Germany_Bell_Beaker",1])
cat("Yamnaya ancestry in German Bell Beakers", 
    (ne - de)/(ne - ru), "\n")

cz = mean(mod$u[meta_ancient$Group.ID ==  "Czech_Bell_Beaker",1])
cat("Yamnaya ancestry in Czech Bell Beakers", 
    (ne - cz)/(ne - ru), "\n")


cw = mean(mod$u[meta_ancient$Group.ID %in% c("Germany_Corded_Ware", 
                                             "Czech_Corded_Ware"),1])
cat("Yamnaya ancestry in Corded Ware samples", 
    (ne - cw)/(ne - ru), "\n")

eba = mean(mod$u[meta_ancient$Group.ID %in% c("Germany_Unetice_EBA", 
                                             "Czech_EBA"),1])
cat("Yamnaya ancestry in Unetice EBA samples", 
    (ne - eba)/(ne - ru), "\n")

ancestry <- (ne - mod$u[meta_ancient$Country %in% c("Germany", 
                                                    "Czech Republic"),1])/(ne - ru)
ancestry <- sapply(ancestry, trunc)

barplot(rbind(ancestry, 1 - ancestry), 
        space = 0, border = NA,
        col = c("darkblue", "salmon3"),
        ylab = "Ancestry Proportion",
        xlab = "Individuals",
        main = "German/Czech Bronze age")

axis(1, 
     at = 1:length(ancestry), 
     labels = as.character(meta_ancient[meta_ancient$Country %in% c("Germany", 
                                                                   "Czech Republic"),]$Instance.ID),
     las = 3, 
     cex.axis = .4)



#### STRUCTURE analysis with LEA::snmf

write.geno(R = geno_filt[age_filt > 0,][lst_keep,], 
           output.file = "target.geno")

proj <- snmf("target.geno", 
             K = 2, 
             entropy = TRUE, 
             alpha = 1000, 
             project = "new")

barchart(proj, K = 2, 
         space = 0, border = NA,
         col = c("darkblue","salmon3"), 
         xlab = "", 
         ylab = "Ancestry proportions", 
         main = "Ancestry matrix") -> bp

axis(1, at = 1:length(bp$order), 
     labels = as.character(meta_ancient$Group.ID)[bp$order], 
     las = 3, cex.axis = .4)

Qm <- Q(proj, K = 2)[meta_ancient$Country %in% c("Germany", 
                                                  "Czech Republic"),]
## Check Yamnaya ancestry from the barplot (max or min)
cat("Yamnaya ancestry from STRUCTURE =", 
    round(max(mean(Qm[,1]),mean(Qm[,2])),3),"\n")
#round(sd(Qm[,1]),3)

## clean environment
remove.snmfProject("target.snmfProject")


### ADMIXTOOLS ANALYSIS


snps <- eigenstrat("./snps/snps")


result <- qpAdm(
  target = c("Germany_Bell_Beaker", 
             "Czech_Bell_Beaker", 
             "Germany_Corded_Ware", 
             "Czech_Corded_Ware", 
             "Germany_Unetice_EBA",
             "Czech_EBA"),
  sources = c("Russia_Yamnaya_Samara",  "Anatolia_N"),
  outgroups = c("YRI.SG", "Russia_Sidelkino_HG.SG", "France_Ranchot88_published"),
  data = snps
)

View(result$proportions)




