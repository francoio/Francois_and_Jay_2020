### Admixture estimates in two-way admixture analyses
### Compared methods: FA, qpAdm, Structure/snmf

library(tfa)
library(LEA)
library(admixr)

target = "England_Bell_Beaker"
meta_ancient <- meta_filt[age_filt >0,]
age_ancient  <- age_filt[age_filt >0]
coverage_ancient <- coverage

## studied population samples
lst_keep <- which( (meta_ancient$Group.ID ==  target |
                        meta_ancient$Group.ID == "Russia_Yamnaya_Samara" |
                        meta_ancient$Group.ID == "Russia_Yamnaya_Samara_published" |
                        meta_ancient$Group.ID ==  "Russia_Yamnaya_Kalmykia.SG"|
                        meta_ancient$Group.ID ==  "Russia_Poltavka"|
                        meta_ancient$Group.ID ==  "Anatolia_N"|
                        meta_ancient$Group.ID ==  "Anatolia_N.SG"#|
#                       meta_ancient$Group.ID ==  "Serbia_Iron_Gates_HG"
                        ) & coverage > 0.7
)

## Outlier individuals in 2 population samples (to be removed)
## Hungary_Longobard: lst_keep <- lst_keep[-c(31,36,44)]
## Italy_Longobard: lst_keep <- lst_keep[-c(32,42, 43)]

geno_ancient <- geno_ancient_cor[lst_keep, ]
meta_ancient <- meta_ancient[lst_keep, ]
age_ancient <- age_ancient[lst_keep]
coverage_ancient <- coverage_ancient[lst_keep] 
meta_ancient <- meta_ancient[,c("Instance.ID", 
                                "Group.ID", 
                                "Country", 
                                "Coverage")]

## Adjust FA model

drift_param = 2e-3

mod  <-  tfa(age_ancient, 
            geno_ancient, 
            k = 2,
            lambda = drift_param,                     
            center = TRUE,
            coverage = coverage_ancient,
            log = TRUE)

## check that temporal dependencies are removed from factor 2
summary(lm(mod$u[,2]~age_ancient))

## results for FA

plot(mod$u, 
     pch = 19, cex = 2, col = "grey90",
     xlab = "Factor 1", ylab = "Factor 2")

points(mod$u[meta_ancient$Group.ID == target,], 
       pch = 19, cex = 1, col = "yellow3")
points(mod$u[meta_ancient$Country == "Russia",], 
       pch = 8, cex = .6, col = "darkblue")
points(mod$u[meta_ancient$Country == "Turkey",], 
       pch = 8, cex = .6, col = "salmon3")
abline(h = 0, lty = 2, col = "orange")

### Yamnaya ancestry evaluation  
ne = mean(mod$u[meta_ancient$Country == "Turkey",1])
ru = mean(mod$u[meta_ancient$Country == "Russia",1])
tg = mean(mod$u[meta_ancient$Group.ID ==  target,1])

cat("Yamnaya ancestry from FA =",(ne - tg)/(ne - ru),"\n")

## Individual ancestry
trunc01 <- function(x){min(1, max(x,0))} 

ancestry <- (ne - mod$u[meta_ancient$Group.ID == target,1])/(ne - ru)
ancestry <- sapply(ancestry, trunc01)

barplot(rbind(ancestry, 1 - ancestry), 
        space = 0,
        col = c("darkblue", "salmon3"),
        ylab = "Ancestry Proportion",
        xlab = "Individuals")

axis(1, 
     at = 1:length(ancestry), 
     labels = as.character(meta_ancient[meta_ancient$Group.ID == target,]$Instance.ID),
     las = 3, 
     cex.axis = .4)

## mean ancestry 

cat("Yamnaya ancestry from FA =",round(mean(ancestry),3),"\n")
#round(sd(ancestry),3)


#### STRUCTURE analyse with LEA::snmf

write.geno(R = geno_filt[age_filt > 0,][lst_keep,], 
           output.file = "target.geno")

proj <- snmf("target.geno", 
             K = 2, 
             entropy = TRUE, 
             alpha = 1000, 
             project = "new")

barchart(proj, K = 2, 
         space = 0, 
         col = c("darkblue","salmon3"), 
         xlab = "Individuals", 
         ylab = "Ancestry proportions", 
         main = "Ancestry matrix") -> bp

axis(1, at = 1:length(bp$order), 
     labels = as.character(meta_ancient$Group.ID)[bp$order], 
     las = 3, cex.axis = .4)

Qm <- Q(proj, K = 2)[meta_ancient$Group.ID == target,]
#Qm <- Q(proj, K = 3)[meta_ancient$Group.ID == target,]


## Check Yamnaya ancestry from the barplot (max or min)
cat("Yamnaya ancestry from STRUCTURE =", 
    round(max(mean(Qm[,1]),mean(Qm[,2])),3),"\n")
#round(sd(Qm[,1]),3)

## clean environment
remove.snmfProject("target.snmfProject")



##########################
### ADMIXTOOLS ANALYSIS
### qpAdm

snps <- eigenstrat("snps/snps")

result <- qpAdm(
  target = c(target),
  sources = c("Russia_Yamnaya_Samara","Anatolia_N"),
  outgroups = c("YRI.SG", 
                "Russia_Sidelkino_HG.SG", 
                #"Serbia_Iron_Gates_HG"
                "France_Ranchot88_published"
                #"Sweden_Motala_HG"
                #"Russia_Samara_HG"
                #"Latvia_HG"
                ),
  data = snps
)

result$proportions

cat("Yamnaya ancestry from qpAdm =", 
    as.numeric(result$proportions[2]),"\n")

#######################################################






###### two-way ancestry analyses for European populations
###### Panels of FIGURE S10

tab <- read.table("table_ancestry.txt", head = TRUE)
name <- as.character(tab$target)

barplot(t(tab[,5:6]), 
        col = c("darkblue","salmon2"), 
        space = 0,
        main = "FA")
axis(1, at = 1:nrow(tab) - 0.5, labels = name, las = 3, cex.axis = .4)

barplot(t(tab[,2:3]), 
        col = c("darkblue","salmon2"), 
        space = 0,
        main = "qpAdm")
axis(1, at = 1:nrow(tab) - 0.5, labels = name, las = 3, cex.axis = .4)

barplot(t(tab[,8:9]), 
        col = c("darkblue","salmon2"), 
        space = 0,
        main = "Structure/snmf")
axis(1, at = 1:nrow(tab) - 0.5, labels = name, las = 3, cex.axis = .4)



######### FIGURE S11

tab <- read.table("table_ancestry.txt", head = TRUE)
name <- as.character(tab$target)

attach(tab)
par(mfrow = c(1,2))
plot(Yamnaya1, Yamnaya2, pch = 19, col = "yellow3",
     xlab = "qpAdm", ylab = "Factor Analysis")
abline(lm(Yamnaya2 ~ Yamnaya1), col = "red")
abline(0,1, col = "grey50", lwd = 2, lty = 2)
points(Yamnaya1, Yamnaya2, pch = 21)

plot(Yamnaya1, Yamnaya3, pch = 19, col = "yellow3",
     xlab = "qpAdm", ylab = "Structure")
abline(lm(Yamnaya3 ~ Yamnaya1), col = "red")
abline(0,1, col = "grey50", lwd = 2, lty = 2)
points(Yamnaya1, Yamnaya3, pch = 21)
detach(tab)

## Figure S12

library(mapplots)
library(maps)

tab <- read.table("table_ancestry.txt", head = TRUE)
name <- as.character(tab$target)


coord <- NULL
for (i in name[3:9]){
  long <- as.numeric(as.character(meta_filt[meta_filt$Group.ID == i,"Long."]))
  lat <- as.numeric(as.character(meta_filt[meta_filt$Group.ID == i,"Lat."]))
  coord <- rbind(coord, c(mean(long),mean(lat)))
}

coord <- rbind(coord, c(-10, 37))
coord <- rbind(coord, c(22, 62))

plot(coord, xlab = "Longitude", ylab = "Latitude", type = "n")
map(add = T, col = "grey90", fill = TRUE)

for (i in 3:9){
  add.pie(z = as.numeric(tab[i,5:6]), 
          x = coord[i-2,1], y = coord[i-2,2], radius = 1.5,
          labels = "",col = c("darkblue","salmon2"))}


