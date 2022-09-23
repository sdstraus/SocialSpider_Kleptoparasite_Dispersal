library(stringr)
library(vcfR)
library(adegenet)
library(SNPRelate)
library(tidyverse)
library(hierfstat)
library(poppr)

library(dartR)

setwd("/Users/sam/PhDThesis/GBS/Rscripts/scripts")

### Eximius ####
# set up genlight #

# done, don't need to redo, just read in RDS
# for eximius
vcf <- read.vcfR("../../filtering/exi/minDP5.maxmiss.60.80.snps.recode.vcf")
pop.data <- read.table("../../pop_maps/eximius_popmap_colony_site_cluster.txt", sep = '\t', header = FALSE)
# 
# 
# #convert to genlight object
gl.exi <- vcfR2genlight(vcf)
ploidy(gl.exi) <- 2
# 
# # make sure popmap matches genlight individuals
pop.new <- as.data.frame(as.factor(indNames(gl.exi)))
colnames(pop.new)[1] <- "Ind"
cleaned <- levels(pop.new$Ind)
pop.data.cleaned <- pop.data %>% filter(V1 %in% cleaned)
# 
# # define population strata
strata_df <- pop.data.cleaned[, c(1, 4, 3, 2)]
strata_df <- strata_df %>%
  dplyr::rename(Site = V2, Cluster = V3, Nest = V4)

# how many clusters?
temp <- strata_df
temp$Nest <- as.factor(temp$Nest)
temp$Cluster <- as.factor(temp$Cluster)
temp$Site <- as.factor(temp$Site)
str(temp)


# # some ways to manipulate strata
strata(gl.exi) <- strata_df
head(strata(gl.exi))
nest.site <- (strata(gl.exi, ~Nest/Site)) #show hierarchically
(strata(gl.exi, ~Nest/Site))
length(indNames(gl.exi))


# have to use genind object
gi.exi <- gl2gi(gl.exi)
strata(gi.exi) <- strata_df


setPop(gi.exi) <- ~Site
wc(gi.exi)


#### F1 ####

# 
vcf.n1 <- read.vcfR("../../filtering/n1/minDP5.maxmiss.60.80.snps.recode.vcf")
pop.data.n1 <- read.table("../../pop_maps/n1_popmap_colony_site_cluster.txt", sep = '\t', header = FALSE)

#convert to genlight object
gl.n1 <- vcfR2genlight(vcf.n1)
ploidy(gl.n1) <- 2

# # make sure popmap matches genlight individuals
pop.new.n1 <- as.data.frame(as.factor(indNames(gl.n1)))
colnames(pop.new.n1)[1] <- "Ind"
cleaned <- levels(pop.new.n1$Ind)
pop.data.cleaned.n1 <- pop.data.n1 %>% filter(V1 %in% cleaned)

strata_df_n1 <- pop.data.cleaned.n1
strata_df_n1 <- strata_df_n1 %>%
  dplyr::rename(Nest = V2, Site = V4, Cluster = V3)

strata(gl.n1) <- strata_df_n1
head(strata(gl.n1))


# have to use genind object
gi.n1 <- gl2gi(gl.n1)
strata(gi.n1) <- strata_df_n1

# Pop as Nest
setPop(gi.n1) <- ~Site
wc(gi.n1)

### F2 ###


vcf.f2 <- read.vcfR("../../filtering/f2/minDP5.maxmiss.60.80.snps.recode.vcf")
pop.data.f2 <- read.table("../../pop_maps/f2_popmap_colony_site.txt", sep = '\t', header = FALSE)
#
# #convert to genlight object
gl.f2 <- vcfR2genlight(vcf.f2)
ploidy(gl.f2) <- 2

# make sure popmap matches genlight individuals
pop.new.f2 <- as.data.frame(as.factor(indNames(gl.f2)))
colnames(pop.new.f2)[1] <- "Ind"
cleaned <- levels(pop.new.f2$Ind)
pop.data.cleaned.f2 <- pop.data.f2 %>% filter(V1 %in% cleaned)
pop.data.cleaned.f2$V3 <- "JatunSacha"


strata_df_f2 <- pop.data.cleaned.f2
strata_df_f2 <- strata_df_f2 %>%
  dplyr::rename(Nest = V2) %>%
  dplyr::rename(Site = V3)

strata(gl.f2) <- strata_df_f2
head(strata(gl.f2))

# have to use genind object
gi.f2 <- gl2gi(gl.f2)
strata(gi.f2) <- strata_df_f2

# Pop as Nest
setPop(gi.f2) <- ~Nest
wc(gi.f2)
