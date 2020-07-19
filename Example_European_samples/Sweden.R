### 2-way admixture analysis for Sweden 
library(tfa)
library(LEA)
library(admixr)

geno_ancient <- geno_filt[age_filt >0,]
meta_ancient <- meta_filt[age_filt >0,]
age_ancient  <- age_filt[age_filt >0]
coverage_ancient = as.numeric(as.character(meta_ancient$Coverage))

lst_keep <- which( (  meta_ancient$Group.ID ==  "Sweden_Viking.SG" |
                      meta_ancient$Group.ID == "Sweden_Motala_HG" |
                      meta_ancient$Group.ID == "Latvia_HG" |                       
                      meta_ancient$Group.ID == "Russia_Yamnaya_Samara" |
                      meta_ancient$Group.ID == "Russia_Yamnaya_Samara_published" |
                      meta_ancient$Group.ID ==  "Russia_Yamnaya_Kalmykia.SG"|
                      meta_ancient$Group.ID ==  "Russia_Poltavka"|
                      meta_ancient$Group.ID ==  "Anatolia_N"|
                      meta_ancient$Group.ID ==  "Anatolia_N.SG" ) # & coverage > 0.7
)

geno_ancient <- geno_ancient_cor[lst_keep, ]
meta_ancient <- meta_ancient[lst_keep, ]
age_ancient <- age_ancient[lst_keep]
coverage_ancient <- coverage_ancient[lst_keep] 
meta_ancient <- meta_ancient[,c("Instance.ID", "Group.ID", "Country", "Coverage")]

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
}

r_sw <- r_squared

plot(llambda, 
     r_sw, 
     xlab = "log(lambda)",
     ylab = "Variance explained by factors",
     main = "Scandinavia")

points(llambda, 
       r_sw, 
       col = "lightblue", 
       pch = 19,
       cex = .8)

abline(v = log10(1e-4), col = "orange", lty = 2)




mod  =  tfa(age_ancient, 
                   geno_ancient, 
                   k = 2,
                   lambda = 3e-3,                     
                   center = TRUE,
                   coverage = coverage_ancient,
                   log = TRUE)


plot(mod$u, pch = 19, cex = 2, col = "grey90",
     xlab = "Factor 1", ylab = "Factor 2",
     main = "Scandinavia",las = 1)
  

points(mod$u[meta_ancient$Country == "Russia",], pch = 8, cex = .6, col = "darkblue")
points(mod$u[meta_ancient$Country == "Turkey",], pch = 8, cex = .6, col = "salmon3")

points(mod$u[meta_ancient$Group.ID ==  "Latvia_HG",], pch = 8, 
       cex = .6, col = "olivedrab")
points(mod$u[meta_ancient$Group.ID ==  "Sweden_Viking.SG",],
       pch = 8, cex = .6, col = "yellow3")

points(mod$u[meta_ancient$Group.ID == "Sweden_Motala_HG",], pch = 2, 
       cex = .6, col = "green4")



abline(h = 0, lty = 2, col = "orange")







### Evaluation of ancestry 
ne = mean(mod$u[meta_ancient$Country == "Turkey",1])
ya = mean(mod$u[meta_ancient$Country == "Russia",1])
hg = mean(mod$u[meta_ancient$Group.ID ==  "Latvia_HG",1])

(ne - ya)/(ne - hg)  #


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
