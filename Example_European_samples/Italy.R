### 2-way admixture analysis for Roman Italy

library(tfa)

geno_ancient <- geno_filt[age_filt >0,]
meta_ancient <- meta_filt[age_filt >0,]
age_ancient  <- age_filt[age_filt >0]

coverage_ancient = as.numeric(as.character(meta_ancient$Coverage))

lst_keep <- which( (  meta_ancient$Group.ID ==  "Italy_Bell_Beaker" |
                      meta_ancient$Group.ID ==  "Italy_Langobard" |
                      meta_ancient$Group.ID == "Italy_Iceman_MN.SG" |
                      meta_ancient$Group.ID == "Russia_Yamnaya_Samara" |
                      meta_ancient$Group.ID == "Russia_Yamnaya_Samara_published" |
                      meta_ancient$Group.ID ==  "Russia_Yamnaya_Kalmykia.SG"|
                      meta_ancient$Group.ID ==  "Russia_Poltavka"|
                      meta_ancient$Group.ID ==  "Anatolia_N"|
                      meta_ancient$Group.ID ==  "Anatolia_N.SG" )  
                   & coverage > 0.7
                   & !(meta_ancient$Instance.ID %in% c("CL146","CL84","CL93") )   
)

geno_ancient <- geno_ancient_cor[lst_keep, ]
meta_ancient <- meta_ancient[lst_keep, ]
age_ancient <- age_ancient[lst_keep]
coverage_ancient <- coverage_ancient[lst_keep] 
meta_ancient <- meta_ancient[,c("Instance.ID", 
                                "Group.ID", 
                                "Country", 
                                "Coverage")]

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

r_sq_it <- r_squared

plot(llambda, 
     r_sq_it, 
     xlab = "log(lambda)",
     ylab = "Variance explained by factors",
     main = "Italy/Langobard",
     las = 1)

points(llambda, 
       r_sq_it, 
       col = "lightblue", 
       pch = 19,
       cex = .8)

abline(v = log10(7e-5), col = "orange", lty = 2)




mod  =  tfa(age_ancient, 
                   geno_ancient, 
                   k = 2,
                   lambda = 7e-5,                     
                   center = TRUE,
                   coverage = coverage_ancient,
                   log = TRUE)


plot(mod$u, pch = 19, cex = 2, col = "grey90", main = "Italy/Langobard",
     xlab = "Factor 1", ylab = "Factor 2" , ylim = c(-20,15), las = 1)


points(mod$u[meta_ancient$Country == "Russia",], pch = 8, cex = .6, col = "darkblue")
points(mod$u[meta_ancient$Country == "Turkey",], pch = 8, cex = .6, col = "salmon3")

points(mod$u[meta_ancient$Group.ID ==  "Italy_Bell_Beaker",], 
       pch = 19, cex = .6, col = "yellow3")

points(mod$u[meta_ancient$Group.ID ==  "Italy_Langobard",], 
       pch = 8, cex = .6, col = "yellow3")

m_iceman <- matrix(mod$u[meta_ancient$Group.ID == "Italy_Iceman_MN.SG",], nrow = 1)
m_iceman[,1] <- -16.23

points(m_iceman, pch = 21, 
       cex = 1.2, col = "salmon")


abline(h = 0, lty = 2, col = "orange")

