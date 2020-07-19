#####
# Read present-day and ancient DNA samples (704 samples)
#####

library(data.table)
library(tfa)


## Download filtered genotypes and metadata for ancient and present-day samples

file_url <- "http://membres-timc.imag.fr/Olivier.Francois/Francois_and_Jay_2020/Ancient_DNA.zip"
download.file(url = file_url)

## Loading genotypes for ancient samples (filtered data)
geno_filt <- fread("../geno/geno_filter.lfmm", header = FALSE)
geno_filt <- as.matrix(geno_filt)
meta_filt <- read.csv("../geno/meta_filter.csv")
age_filt <- as.matrix(read.table("../geno/age_filter.txt"))[,1]

## Adjust for coverage
coverage = as.numeric(as.character(meta_filt[age_filt > 0,]$Coverage))
geno_ancient_cor <- coverage_adjust(Y = geno_filt[age_filt > 0,], 
                                    coverage = coverage, 
                                    K = 13,
                                    log = TRUE)
