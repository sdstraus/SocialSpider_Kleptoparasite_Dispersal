library(stringr)
library(vcfR)
library(adegenet)
library(SNPRelate)
library(tidyverse)
library(hierfstat)
library(poppr)
library(stringr)
library(dartR)

setwd("/Users/sam/PhDThesis/GBS/github")

# Eximius -------------------------
#### set up genlight ###### 

# done, don't need to redo, just read in RDS
# for eximius
vcf <- read.vcfR("data/filtered_vcfs/exi_minDP3.maxmiss.50.70.snps.recode.vcf")

#convert to genlight object
gl.exi <- vcfR2genlight(vcf)
ploidy(gl.exi) <- 2

pop.data <- read.table("data/pop_maps/eximius_popmap_colony_site_cluster.txt", sep = '\t', header = FALSE)
# gl.exi <- readRDS("rds/gl.exi.RDS")


## make sure popmap matches genlight individuals
pop.new <- as.data.frame(as.factor(indNames(gl.exi)))
colnames(pop.new)[1] <- "Ind"
cleaned <- levels(pop.new$Ind)
pop.data.cleaned <- pop.data %>% filter(V1 %in% cleaned)
 
## define population strata
strata_df <- pop.data.cleaned[, c(1, 4, 3, 2)]
strata_df <- strata_df %>%
  dplyr::rename(Site = V2, Cluster = V3, Nest = V4)

# how many clusters?
temp <- strata_df
temp$Nest <- as.factor(temp$Nest)
temp$Cluster <- as.factor(temp$Cluster)
temp$Site <- as.factor(temp$Site)
str(temp)


saveRDS(gl.exi, file="rds/gl.exi.RDS")
write.csv(pop.data.cleaned, "data/output_files/eximius.pop.data.cleaned.csv")



###### Subset populations ######

setPop(gl.exi) <- ~Site
popNames(gl.exi)
exi.js <- popsub(gl.exi, "JatunSacha")
popNames(exi.js)
setPop(exi.js) <- ~Nest
indNames(exi.js)
# gl.exi <- exi.js
# count of how many left after filtering per nest
indNames(exi.js) %>% str_subset("Ex02") %>% length()

exi.archi <- popsub(gl.exi, "Archidona")
setPop(exi.archi) <- ~Nest
popNames(exi.archi)
indNames(exi.archi)
# remove NAs
toRemove <- is.na(glMean(exi.archi, alleleAsUnit = FALSE))
which(toRemove)
exi.archi <- exi.archi[, !toRemove]
indNames(exi.archi) %>% str_subset("Exi_Archi_11.8.5") %>% length()


exi.vl <- popsub(gl.exi, sublist=c("ViaLoreto"))
indNames(exi.vl)
# remove NAs
toRemove <- is.na(glMean(exi.vl, alleleAsUnit = FALSE))
which(toRemove)
exi.vl <- exi.vl[, !toRemove]
setPop(exi.vl) <- ~Nest
popNames(exi.vl)
indNames(exi.vl) %>% str_subset("Exi_VL_10.5") %>% length()


######### Observed Heterozygosity #########

mean(gl.Ho(gl.exi))
mean(gl.Ho(exi.js))
mean(gl.Ho(exi.archi))
mean(gl.Ho(exi.vl))

mean(gl.He(gl.exi))
mean(gl.He(exi.js))
mean(gl.He(exi.archi))
mean(gl.He(exi.vl))

####### Weir & Cockrham F-stats ########

# have to use genind object
gi.exi <- gl2gi(gl.exi)
strata(gi.exi) <- strata_df

gi.exi <- missingno(gi.exi)
exi.rich <- PopGenReport::allel.rich(gi.exi)$mean.richness

js.exi.gi <- gl2gi(exi.js)
strata(js.exi.gi) <- strata_df %>% filter(Site == "JatunSacha")


archi.exi.gi <- gl2gi(exi.archi)
strata(archi.exi.gi) <- strata_df %>% filter(Site == "Archidona")


vl.exi.gi <- gl2gi(exi.vl)
strata(vl.exi.gi) <- strata_df %>% filter(Site == "ViaLoreto")


saveRDS(gi.exi, file="rds/gi.exi.RDS")
saveRDS(js.exi.gi, file="rds/js.exi.gi.RDS")
saveRDS(archi.exi.gi, file="rds/archi.exi.gi.RDS")
saveRDS(vl.exi.gi, file="rds/vl.exi.gi.RDS")


gi.exi <- readRDS("rds/gi.exi.RDS")
js.exi.gi <- readRDS("rds/js.exi.gi.RDS")
archi.exi.gi <- readRDS("rds/archi.exi.gi.RDS")
vl.exi.gi <- readRDS("rds/vl.exi.gi.RDS")


#gi.exi<- readRDS("Rscripts/gi.exi.RDS")

# Pop as Site (FIR)
setPop(gi.exi) <- ~Site
wc(gi.exi)


# Pop as Nest
setPop(js.exi.gi) <- ~Nest
pop(js.exi.gi)
wc(js.exi.gi)

setPop(archi.exi.gi) <- ~Nest
pop(archi.exi.gi)
wc(archi.exi.gi)

setPop(vl.exi.gi) <- ~Nest
pop(vl.exi.gi)
wc(vl.exi.gi)


#### AMOVA #####

exi.amova <- poppr.amova(gi.exi, ~Site/Cluster/Nest)
amova.test <- randtest(exi.amova)
plot(amova.test)

#### isolation by distance ######

## set up gps coords for eximius genind 
exi.coords <- read.csv("data/eximius_coordinates.csv")
exi.coords$Nest <- exi.coords$Nest %>% str_replace_all("Ex", "Exi")
exi.coords$Nest <- exi.coords$Nest %>% str_replace_all("Exii", "Exi")

strata_df$lat <- 0
strata_df$long <- 0

for(i in 1:length(strata_df$Nest)){
  for(j in 1:length(exi.coords$Nest)){
    if(strata_df$Nest[i] == exi.coords$Nest[j]){
      strata_df$lat[i] <- exi.coords$latitude[j]
      strata_df$long[i] <- exi.coords$longitude[j]
    }}}

## all archidona colonies share one set of coordinates (all within ~30 meters)
row <- which(exi.coords$Nest == "Exi_Archi_11.8")
archi.lat <- exi.coords$latitude[row]
archi.long <- exi.coords$longitude[row]

rows <- which(strata_df$Site == "Archidona")
strata_df$lat[rows] <- archi.lat
strata_df$long[rows] <- archi.long

## VL_9.5.1 didn't get it's coords, do manually
row <- which(exi.coords$Nest == "Exi_VL_9.5")
vl_9.5_lat <- exi.coords$latitude[row]
vl_9.5.long <- exi.coords$longitude[row]

rows <- which(strata_df$Cluster == "Exi_VL_9.5")
strata_df$lat[rows] <- vl_9.5_lat
strata_df$long[rows] <- vl_9.5.long


# add to strata
xy <- strata_df %>% dplyr::select(long, lat)
gi.exi@other$xy <- xy
gi.exi@other$xy


## setting up distance matrices 
Dgen <- poppr::provesti.dist(gi.exi) # genetic distance

Dgeo <- dist(gi.exi@other$xy, method = "euclidean") # geographical distance
library(ade4)
ibd <- mantel.randtest(Dgen, Dgeo)
ibd

plot(ibd, main = "Mantel's test")


## js
js.exi.gi <- gl2gi(exi.js)
strata(js.exi.gi) <- strata_df %>% filter(Site == "JatunSacha")

xy <- strata(js.exi.gi) %>% dplyr::select(long, lat) %>% filter()
js.exi.gi@other$xy <- xy
js.exi.gi@other$xy

Dgen <- poppr::provesti.dist(js.exi.gi) # genetic distance

Dgeo <- dist(js.exi.gi@other$xy, method = "euclidean") # geographical distance
library(ade4)
ibd <- mantel.randtest(Dgen, Dgeo)
ibd

plot(ibd, main = "Mantel's test")


# archidona
## all nests share same coordinate -- too close in space

archi.exi.gi <- gl2gi(exi.archi)
strata(archi.exi.gi) <- strata_df %>% filter(Site == "Archidona")

xy <- strata(archi.exi.gi) %>% dplyr::select(long, lat) %>% filter()
archi.exi.gi@other$xy <- xy
archi.exi.gi@other$xy

Dgen <- poppr::provesti.dist(archi.exi.gi) # genetic distance

Dgeo <- dist(archi.exi.gi@other$xy, method = "euclidean") # geographical distance
library(ade4)
ibd <- mantel.randtest(Dgen, Dgeo)
ibd

plot(ibd, main = "Mantel's test")


# via loreto
vl.exi.gi <- gl2gi(exi.vl)
strata(vl.exi.gi) <- strata_df %>% filter(Site == "ViaLoreto")


xy <- strata(vl.exi.gi) %>% dplyr::select(long, lat)
vl.exi.gi@other$xy <- xy
vl.exi.gi@other$xy

Dgen <- poppr::provesti.dist(vl.exi.gi) # genetic distance

Dgeo <- dist(vl.exi.gi@other$xy, method = "euclidean") # geographical distance
library(ade4)
ibd <- mantel.randtest(Dgen, Dgeo)
ibd

plot(ibd, main = "Mantel's test")


#### admixture plots ####

CVs <- read.table("data/admixture/exi/CV.csv", sep = " ")
CVs <- CVs[, 3:4] ## drop the first two columns
## Remove the formatting around the K values:
CVs[, 1] <- gsub(x = CVs[, 1], pattern = "\\(K=",
                 replacement = "")
CVs[, 1] <- gsub(x = CVs[, 1], pattern = "\\):",
                 replacement = "") 
head(CVs)


plot(CVs, xlab = "K", ylab = "A. eximius CV error")
abline(v = c(6, 15))


## K = 6 ##
tbl=read.table("data/admixture/exi/exi.scaffolds.6.Q")
barplot(t(as.matrix(tbl)), col=getPalette(15), xlab="Individual", ylab="Ancestry", border=NA)

popmap <-  read.table("data/pop_maps/eximius_popmap_colony_site_cluster.txt", sep = '\t', header = FALSE)
popmap <- popmap %>% rename("ind" = "V1", "site" = "V2", "nest" = "V4", "cluster" = "V3")
head(popmap)

names <-gl.exi@ind.names
names2 <- popmap$ind
difference <- setdiff(names2, names)

for(i in 1:length(popmap$ind)){
  for(j in 1:length(difference)){
    if(popmap$ind[i] == difference[j]){
      popmap <- popmap[-i,]
    }
  }
}

exi.df <- cbind(popmap, tbl)
rows <- which(exi.df$ind == "Exi_VL_9.5.01")
exi.df <- exi.df[-rows,]
long.dat <- exi.df %>% pivot_longer(cols = V1:V6, names_to = "Kclust", values_to = "fraction", names_prefix = "V")

# trim "Exi_Archi" off nest names
long.dat <- long.dat %>% mutate(nest = str_remove(nest, "Exi_Archi_"))
# trim "Exi_VL" off nest names
long.dat <- long.dat %>% mutate(nest = str_remove(nest, "Exi_VL_"))
# _JS
long.dat <- long.dat %>% mutate(nest = str_remove(nest, "_JS"))

# get individual #
long.dat$ind.num <- 0
for(i in 1:length(long.dat$site)){
  if(long.dat$site[i] == "Archidona"){
    long.dat$ind.num[i] <- str_sub(long.dat$ind[i], -2)
  }
  else if(long.dat$site[i] == "ViaLoreto"){
    long.dat$ind.num[i] <- str_sub(long.dat$ind[i], -2)
  }
  else if(long.dat$site[i] == "JatunSacha"){
    long.dat$ind.num[i] <- str_sub(long.dat$ind[i], -2)
  }
}

long.dat <- long.dat %>% mutate(ind.num2 = as.numeric(ind.num))
long.dat$ind.num2 <- if_else(long.dat$ind.num2 < 1.0, true = (long.dat$ind.num2*10), false = long.dat$ind.num2)

long.dat <- long.dat %>% mutate(ind.num2 = as.character(ind.num2)) %>% 
  mutate(site = as.factor(site))

# rename levels within a site
levels(long.dat$site)
levels(long.dat$site) <- c("A", "JS", "VL")

## reorder VL nest order
long.dat <- long.dat %>% mutate(nest = as.factor(nest))
levels(long.dat$nest)
## archidona in numerical order
## via loreto by km
## js north to south/ distance from road
long.dat$nest <- factor(long.dat$nest, levels = c("11.8.1", "11.8.2", "11.8.4", 
                                                    "11.8.5", "9.5.1", "10.5", "17.5","24.1", 
                                                    "Bridge.1", "41.2", "Exi06A", "Exi06B", "Exi07", 
                                                    "Exi11", "Exi12", "Exi02", "Exi10", "Exi05", "Exi08", 
                                                    "Exi14"))

library(RColorBrewer)
display.brewer.pal(12, "Paired")
brewer.pal(12, "Paired")
pal <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C")

exi.admix6 <-  ggplot(long.dat, aes(x=ind.num2, y=fraction, fill=Kclust)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(. ~ site + nest, drop=TRUE, space="free", scales="free")+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(fill = "K cluster", y = "Ancestry", x = "Individual")+
  scale_fill_manual(values = pal)+
  theme(panel.grid=element_blank()) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.background=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),)+
  theme(strip.text.x.top = element_text(angle = 90))
exi.admix6

exi.admix.site.level <- ggpubr::ggarrange(exi.admix6, exi.admix15, nrow = 2)
ggsave(exi.admix.site.level, filename = "../../figures/exi.admix.site.level.jpeg", dpi = "retina",
       units = "in", width = 12, height = 10)

####### exi - js ####### 

CVs <- read.table("data/admixture/exi_js/CV.csv", sep = " ")
CVs <- CVs[, 3:4] ## drop the first two columns
## Remove the formatting around the K values:
CVs[, 1] <- gsub(x = CVs[, 1], pattern = "\\(K=",
                 replacement = "")
CVs[, 1] <- gsub(x = CVs[, 1], pattern = "\\):",
                 replacement = "") 
head(CVs)


plot(CVs, xlab = "K", ylab = "A. eximius CV error")
abline(v = c(6, 8))


library(RColorBrewer)
tbl=read.table("data/admixture/exi_js/exi_subset.js.6.Q")
barplot(t(as.matrix(tbl)), col=brewer.pal(6, "Dark2"), xlab="Individual", ylab="Ancestry", border=NA)

popmap <-  read.table("data/pop_maps/eximius_popmap_colony_site_cluster.txt", sep = '\t', header = FALSE)
popmap <- popmap %>% rename("ind" = "V1", "site" = "V2", "nest" = "V3", "cluster" = "V4") %>% 
  filter(site == "JatunSacha")
head(popmap)

names <-exi.js@ind.names
names2 <- popmap$ind
difference <- setdiff(names2, names)

for(i in 1:length(popmap$ind)){
  for(j in 1:length(difference)){
    if(popmap$ind[i] == difference[j]){
      popmap <- popmap[-i,]
    }
  }
}

exi.js.df <- cbind(popmap, tbl)

long.dat.js <- exi.js.df %>% pivot_longer(cols = V1:V6, names_to = "Kclust", values_to = "fraction", names_prefix = "V")

# trim column names
long.dat.js <- long.dat.js %>% 
  mutate(nest = str_sub(nest,0,5)) %>% 
  mutate(ind.num = str_sub(ind, -2)) 

# we have numbers 0-9 as 0.1-0.9 because of naming convention
# tail(long.dat %>% mutate(ind.num = str_replace(ind.num, ".", " ")))  # this doesn't work, turns 14 to 4

long.dat.js <- long.dat.js %>% mutate(ind.num2 = as.numeric(ind.num))

long.dat.js$ind.num2 <- if_else(long.dat.js$ind.num2 < 1.0, true = (long.dat.js$ind.num2*10), false = long.dat.js$ind.num2)

# also, 06A and 06B naming is causing issues
rows <- which(str_detect(long.dat.js$ind, "06A"))
long.dat.js$ind.num2[rows] <- paste(long.dat.js$ind.num2[rows], "A", sep = "")

rows <- which(str_detect(long.dat.js$ind, "06B"))
long.dat.js$ind.num2[rows] <- paste(long.dat.js$ind.num2[rows], "B", sep = "")

exi.js.admix <-  ggplot(long.dat.js, aes(x=ind.num2, y=fraction, fill=Kclust)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(. ~ nest, drop=TRUE, space="free", scales="free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(fill = "K cluster", y = "Ancestry", x = "Individual")
exi.js.admix <- exi.js.admix + scale_fill_brewer(palette="BrBG")
exi.js.admix


### exi - archi ####
CVs <- read.table("data/admixture/exi_archi/CV.csv", sep = " ")
CVs <- CVs[, 3:4] ## drop the first two columns
## Remove the formatting around the K values:
CVs[, 1] <- gsub(x = CVs[, 1], pattern = "\\(K=",
                 replacement = "")
CVs[, 1] <- gsub(x = CVs[, 1], pattern = "\\):",
                 replacement = "") 
head(CVs)


plot(CVs, xlab = "K", ylab = "A. eximius CV error")
abline(v = c(3, 7)) # 7 doesn't make much sense


tbl=read.table("data/admixture/exi_archi/exi_archi_subset2.3.Q")
barplot(t(as.matrix(tbl)), col=brewer.pal(6, "Dark2"), xlab="Individual", ylab="Ancestry", border=NA)

popmap <-  read.table("data/pop_maps/eximius_popmap_colony_site_cluster.txt", sep = '\t', header = FALSE)
popmap <- popmap %>% rename("ind" = "V1", "site" = "V2", "nest" = "V4", "cluster" = "V3") %>% 
  filter(site == "Archidona")
head(popmap)

names <-exi.archi@ind.names
names2 <- popmap$ind
difference <- setdiff(names2, names)

for(i in 1:length(popmap$ind)){
  for(j in 1:length(difference)){
    if(popmap$ind[i] == difference[j]){
      popmap <- popmap[-i,]
    }
  }
}

exi.archi.df <- cbind(popmap, tbl)

long.dat.archi <- exi.archi.df %>% pivot_longer(cols = V1:V3, names_to = "Kclust", values_to = "fraction", names_prefix = "V")

# trim "Exi_Archi" off nest names
long.dat.archi <- long.dat.archi %>% mutate(nest = str_remove(nest, "Exi_Archi_"))

## want just last 2 characters of individudal # (e.g. 01, 02)
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

long.dat.archi <- long.dat.archi %>% mutate(ind.num = substrRight(ind, 2))

exi.archi.admix <-  ggplot(long.dat.archi, aes(x=ind.num, y=fraction, fill=Kclust)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(. ~ nest, drop=TRUE, space="free", scales="free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(fill = "K cluster", y = "Ancestry", x = "Individual")
exi.archi.admix <- exi.archi.admix + scale_fill_brewer(palette="Spectral")
exi.archi.admix

### exi - vl ####

CVs <- read.table("data/admixture/exi_vl/CV.csv", sep = " ")
CVs <- CVs[, 3:4] ## drop the first two columns
## Remove the formatting around the K values:
CVs[, 1] <- gsub(x = CVs[, 1], pattern = "\\(K=",
                 replacement = "")
CVs[, 1] <- gsub(x = CVs[, 1], pattern = "\\):",
                 replacement = "") 
head(CVs)


plot(CVs, xlab = "K", ylab = "A. eximius CV error")
abline(v = c(3, 5)) # 10 doesn't make much sense


tbl=read.table("data/admixture/exi_vl/exi_vl_subset2.5.Q")
barplot(t(as.matrix(tbl)), col=brewer.pal(6, "Dark2"), xlab="Individual", ylab="Ancestry", border=NA)

popmap <-  read.table("data/pop_maps/eximius_popmap_colony_site_cluster.txt", sep = '\t', header = FALSE)
popmap <- popmap %>% rename("ind" = "V1", "site" = "V2", "nest" = "V3", "cluster" = "V4") %>% 
  filter(site == "ViaLoreto")
head(popmap)

# Exi_VL_9.5.01 shouldn't exist? it's diff from Exi_VL_9.5.1.01, removing bc it's notation error
# rows <- which(popmap$ind == "Exi_VL_9.5.01")
# popmap <- popmap[-rows,]
# dartR::gl.drop.ind(exi.vl, ind.list = c("Exi_VL_9.5.01"))

names <-exi.vl@ind.names
names2 <- popmap$ind
difference <- setdiff(names2, names)

for(i in 1:length(popmap$ind)){
  for(j in 1:length(difference)){
    if(popmap$ind[i] == difference[j]){
      popmap <- popmap[-i,]
    }
  }
}

exi.vl.df <- cbind(popmap, tbl)
rows <- which(exi.vl.df$ind == "Exi_VL_9.5.01")
exi.vl.df <- exi.vl.df[-rows,]



long.dat.vl <- exi.vl.df %>% pivot_longer(cols = V1:V5, names_to = "Kclust", values_to = "fraction", names_prefix = "V")


# trim "Exi_VL" off nest names
long.dat.vl <- long.dat.vl %>% mutate(nest = str_remove(nest, "Exi_VL_"))

# ## want just last 2 characters of individudal # (e.g. 01, 02)
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

long.dat.vl <- long.dat.vl %>% mutate(ind.num = substrRight(ind, 1))


# rename levels within a site
long.dat.vl <- long.dat.vl %>% mutate(site = as.factor(site))

levels(long.dat.vl$site)
levels(long.dat.vl$site) <- c("VL")


long.dat.vl$nest <- as.factor(long.dat.vl$nest)
levels(long.dat.vl$nest )
# reorder levels
long.dat.vl$nest <- factor(long.dat.vl$nest, levels = c("9.5", "10.5", "17.5","24.1", 
                                                  "Bridge.1", "41.2"))

RColorBrewer::display.brewer.pal(12, "Paired")
RColorBrewer::brewer.pal(12, "Paired") # last 6, minus pale yellow
pal <- c("#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928")

exi.vl.admix <-  ggplot(long.dat.vl, aes(x=ind.num, y=fraction, fill=Kclust)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(. ~ nest + site, drop=TRUE, space="free", scales="free", switch = 'x')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(fill = "K cluster", y = "Ancestry", x = "Individual")+
  scale_fill_manual(values = pal)+
  theme(panel.grid=element_blank()) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.background=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),)+
  theme(strip.text.x.bottom = element_text(angle = 90))
exi.vl.admix

exi.within.site <- ggpubr::ggarrange(exi.admix6, exi.js.admix, exi.archi.admix, exi.vl.admix, 
                                     nrow=4, labels = "auto")
ggsave(exi.within.site, filename = "../../figures/exi.within.site.jpeg", dpi = "retina",
       units = "in", width = 10, height = 11)

# Faiditius sp.1 ---------------------------------
# don't need to repeat step
# keep code
# 
vcf.n1 <- read.vcfR("data/filtered_vcfs/f1_minDP3.maxmiss.40.70.snps.recode.vcf")
pop.data.n1 <- read.table("data/pop_maps/n1_popmap_colony_site_cluster.txt", sep = '\t', header = FALSE)

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

saveRDS(gl.n1, file="rds/gl.n1.RDS")
write.csv(pop.data.cleaned.n1, "data/output_files/N1.pop.data.cleaned.csv")
 
gl.n1 <- readRDS("rds/gl.n1.RDS")
 

### subset population ########
indNames(gl.n1) %>% length()
setPop(gl.n1) <- ~Nest
popNames(gl.n1) %>% length()

# get full list of names -> for full filenames of individuals kept after filtering
n1_list <- indNames(gl.n1)
write_lines(paste(n1_list, ".*", sep =""), file = "n1_list.txt")



#JS
setPop(gl.n1) <- ~Site
popNames(gl.n1)
n1.js <- popsub(gl.n1, "JatunSacha")
popNames(n1.js)
setPop(n1.js) <- ~Nest
indNames(n1.js)
# remove NAs
toRemove <- is.na(glMean(n1.js, alleleAsUnit = FALSE))
which(toRemove)
n1.js <- n1.js[, !toRemove]

indNames(n1.js) %>% str_subset("Ex08") %>% length()

# Archi
setPop(gl.n1) <- ~Site
popNames(gl.n1)
n1.archi <- popsub(gl.n1, "Archidona")
popNames(n1.archi)
setPop(n1.archi) <- ~Nest
# remove NAs
toRemove <- is.na(glMean(n1.archi, alleleAsUnit = FALSE))
which(toRemove)
n1.archi <- n1.archi[, !toRemove]
indNames(n1.archi) %>% str_subset("Exi_Archi_11.8.1") %>% length()

# VL
setPop(gl.n1) <- ~Site
popNames(gl.n1)
n1.vl <- popsub(gl.n1, "ViaLoreto")
popNames(n1.vl)
setPop(n1.vl) <- ~Nest
# remove NAs
toRemove <- is.na(glMean(n1.vl, alleleAsUnit = FALSE))
which(toRemove)
n1.vl <- n1.vl[, !toRemove]
indNames(n1.vl)
indNames(n1.vl) %>% str_subset("Exi_VL_24.1") %>% length()


######### Observed Heterozygosity #########

mean(gl.Ho(gl.n1))
mean(gl.Ho(n1.js))
mean(gl.Ho(n1.archi))
mean(gl.Ho(n1.vl))

mean(gl.He(gl.n1))
mean(gl.He(n1.js))
mean(gl.He(n1.archi))
mean(gl.He(n1.vl))
####### Weir & Cockrham F-stats ########

library(dartR)
library(hierfstat)

# have to use genind object
gi.n1 <- gl2gi(gl.n1)
strata(gi.n1) <- strata_df_n1

#strata(gi.n1) <- strata_df_n1
setPop(gi.n1) <- ~Nest/Site
n1.rich <- PopGenReport::allel.rich(gi.n1)$mean.richness

js.n1.gi <- gl2gi(n1.js)
strata(js.n1.gi) <- strata_df_n1 %>% filter(Site == "JatunSacha")

archi.n1.gi <- gl2gi(n1.archi)
strata(archi.n1.gi) <- strata_df_n1 %>% filter(Site == "Archidona")

vl.n1.gi <- gl2gi(n1.vl)
strata(vl.n1.gi) <- strata_df_n1 %>% filter(Site == "ViaLoreto")


# saveRDS(gi.n1, file="rds/gi.n1.RDS")
# saveRDS(js.n1.gi, file="rds/js.n1.gi.RDS")
# saveRDS(archi.n1.gi, file="rds/archi.n1.gi.RDS")
# saveRDS(vl.n1.gi, file="rds/vl.n1.gi.RDS")
# 
# gi.n1 <- readRDS("rds/gi.n1.RDS")
# js.n1.gi <- readRDS("rds/js.n1.gi.RDS")
# archi.n1.gi <- readRDS("rds/archi.n1.gi.RDS")
# vl.n1.gi <- readRDS("rds/vl.n1.gi.RDS")

# Pop as Nest
setPop(gi.n1) <- ~Site
pop(gi.n1)
wc(gi.n1)


# Pop as Nest/Cluster
setPop(js.n1.gi) <- ~Nest
pop(js.n1.gi)
wc(js.n1.gi)

setPop(archi.n1.gi) <- ~Nest
pop(archi.n1.gi)
wc(archi.n1.gi)


setPop(vl.n1.gi) <- ~Nest
pop(vl.n1.gi)
wc(vl.n1.gi)



#### AMOVA #####

n1.amova <- poppr.amova(gi.n1, ~Site/Nest)
n1.test.amova <- randtest(n1.amova)
plot(n1.test.amova)


### IBD ####
## set up gps coords for eximius genind 
exi.coords <- read.csv("data/eximius_coordinates.csv")
exi.coords$Nest <- exi.coords$Nest %>% str_replace_all("Ex", "Exi")
exi.coords$Nest <- exi.coords$Nest %>% str_replace_all("Exii", "Exi")

strata_df_n1$lat <- 0
strata_df_n1$long <- 0

for(i in 1:length(strata_df_n1$Nest)){
  for(j in 1:length(exi.coords$Nest)){
    if(strata_df_n1$Nest[i] == exi.coords$Nest[j]){
      strata_df_n1$lat[i] <- exi.coords$latitude[j]
      strata_df_n1$long[i] <- exi.coords$longitude[j]
    }}}

## all archidona colonies share one set of coordinates (all within ~20 meters)
row <- which(exi.coords$Nest == "Exi_Archi_11.8")
archi.lat <- exi.coords$latitude[row]
archi.long <- exi.coords$longitude[row]

rows <- which(strata_df_n1$Site == "Archidona")
strata_df_n1$lat[rows] <- archi.lat
strata_df_n1$long[rows] <- archi.long

## VL_9.5.1 didn't get it's coords, do manually
row <- which(exi.coords$Nest == "Exi_VL_9.5")
vl_9.5_lat <- exi.coords$latitude[row]
vl_9.5.long <- exi.coords$longitude[row]

rows <- which(strata_df_n1$Cluster == "Exi_VL_9.5")
strata_df_n1$lat[rows] <- vl_9.5_lat
strata_df_n1$long[rows] <- vl_9.5.long

xy <- strata_df_n1 %>% dplyr::select(long, lat)
gi.n1@other$xy <- xy
gi.n1@other$xy

Dgen <- poppr::provesti.dist(gi.n1)
# eucl.m <- (dist(gi.exi@other$xy))
Dgeo <- dist(gi.n1@other$xy, method = "euclidean")
library(ade4)
ibd <- mantel.randtest(Dgen, Dgeo)
ibd

plot(ibd, main = "Mantel's test")


# js only
js.n1.gi <- gl2gi(n1.js)
strata(js.n1.gi) <- strata_df_n1 %>% filter(Site == "JatunSacha")


xy <- strata(js.n1.gi) %>% dplyr::select(long, lat)
js.n1.gi@other$xy <- xy
js.n1.gi@other$xy

Dgen <- poppr::provesti.dist(js.n1.gi) # genetic distance

Dgeo <- dist(js.n1.gi@other$xy, method = "euclidean") # geographical distance
library(ade4)
ibd <- mantel.randtest(Dgen, Dgeo)
ibd

plot(ibd, main = "Mantel's test")


# vl only
vl.n1.gi <- gl2gi(n1.vl)
strata(vl.n1.gi) <- strata_df_n1 %>% filter(Site == "ViaLoreto")


xy <- strata(vl.n1.gi) %>% dplyr::select(long, lat)
vl.n1.gi@other$xy <- xy
vl.n1.gi@other$xy

Dgen <- poppr::provesti.dist(vl.n1.gi) # genetic distance

Dgeo <- dist(vl.n1.gi@other$xy, method = "euclidean") # geographical distance
library(ade4)
ibd <- mantel.randtest(Dgen, Dgeo)
ibd

plot(ibd, main = "Mantel's test")

#### admixture plots ####

CVs <- read.table("data/admixture/f1/CV.csv", sep = " ")
CVs <- CVs[, 3:4] ## drop the first two columns
## Remove the formatting around the K values:
CVs[, 1] <- gsub(x = CVs[, 1], pattern = "\\(K=",
                 replacement = "")
CVs[, 1] <- gsub(x = CVs[, 1], pattern = "\\):",
                 replacement = "") 
head(CVs)


plot(CVs, xlab = "K", ylab = "Faiditus sp.1 CV error")
abline(v = c(8))

## using K = 5 based off BIC

tbl=read.table("data/admixture/f1/f1.scaffolds.5.Q")
barplot(t(as.matrix(tbl)), col=brewer.pal(6, "Dark2"), xlab="Individual", ylab="Ancestry", border=NA)

popmap <-  read.table("data/pop_maps/n1_popmap_colony_site_cluster.txt", sep = '\t', header = FALSE)
popmap <- popmap %>% rename("ind" = "V1", "nest" = "V2", "cluster" = "V3", "site" = "V4")
head(popmap)

names <-gl.n1@ind.names
names2 <- popmap$ind
difference <- setdiff(names2, names)

for(i in 1:length(popmap$ind)){
  for(j in 1:length(difference)){
    if(popmap$ind[i] == difference[j]){
      popmap <- popmap[-i,]
    }
  }
}

f1.df <- cbind(popmap, tbl)

long.dat <- f1.df %>% pivot_longer(cols = V1:V5, names_to = "Kclust", values_to = "fraction", names_prefix = "V")

# trim "Exi_Archi" off nest names
long.dat <- long.dat %>% mutate(nest = str_remove(nest, "Exi_Archi_"))
# trim "Exi_VL" off nest names
long.dat <- long.dat %>% mutate(nest = str_remove(nest, "Exi_VL_"))
# _JS
long.dat <- long.dat %>% mutate(nest = str_remove(nest, "_JS"))

# get individual #
long.dat$ind.num <- 0
for(i in 1:length(long.dat$site)){
  if(long.dat$site[i] == "Archidona"){
    long.dat$ind.num[i] <- str_sub(long.dat$ind[i], -2)
  }
  else if(long.dat$site[i] == "ViaLoreto"){
    long.dat$ind.num[i] <- str_sub(long.dat$ind[i], -2)
  }
  else if(long.dat$site[i] == "JatunSacha"){
    long.dat$ind.num[i] <- str_sub(long.dat$ind[i], -2)
  }
}

long.dat <- long.dat %>% mutate(ind.num2 = as.numeric(ind.num))
long.dat$ind.num2 <- if_else(long.dat$ind.num2 < 1.0, true = (long.dat$ind.num2*10), false = long.dat$ind.num2)

long.dat <- long.dat %>% mutate(ind.num2 = as.character(ind.num2)) %>% 
  mutate(site = as.factor(site))

# rename levels within a site
levels(long.dat$site)
levels(long.dat$site) <- c("A", "JS", "VL")

## reorder VL nest order
long.dat <- long.dat %>% mutate(nest = as.factor(nest))
levels(long.dat$nest)
## archidona in numerical order
## via loreto by km
## js north to south/ distance from road
long.dat$nest <- factor(long.dat$nest, levels = c("11.8.1", "11.8.3", 
                                                  "11.8.5", "9.5.2","24.1", "Exi06B","Exi12", 
                                                  "Exi02", "Exi08", "Exi17"))


f1.admix <-  ggplot(long.dat, aes(x=ind.num2, y=fraction, fill=Kclust)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(. ~ site + nest, drop=TRUE, space="free", scales="free", switch = "x")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(fill = "K cluster", y = "Ancestry", x = "Individual")+
  scale_fill_brewer(palette = "Paired")+
  theme(panel.grid=element_blank()) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.background=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),)+
  theme(strip.text.x.bottom = element_text(angle = 90))
f1.admix

bottom  <- ggpubr::ggarrange(exi.vl.admix, f1.admix, labels = c("b", "c"), nrow=1)
top <- ggpubr::ggarrange(exi.admix6, bottom, labels = c("a", ""), nrow=2)

pdf(file = "../../figures/admixture_plots.pdf", height = 8, width = 12)
top
dev.off()

jpeg(file = "../../figures/Figure5_02Sept22.jpeg", height = 8, width = 12, units = "in", res = 1000)
top
dev.off()

####### f1 - js ########
CVs <- read.table("data/admixture/f1_js/CV.csv", sep = " ")
CVs <- CVs[, 3:4] ## drop the first two columns
## Remove the formatting around the K values:
CVs[, 1] <- gsub(x = CVs[, 1], pattern = "\\(K=",
                 replacement = "")
CVs[, 1] <- gsub(x = CVs[, 1], pattern = "\\):",
                 replacement = "") 
head(CVs)


plot(CVs, xlab = "K", ylab = "A. eximius CV error")
abline(v = c(4, 7)) #sample size probs too low for 7


library(RColorBrewer)
tbl=read.table("data/admixture/f1_js/f1_js_subset.4.Q")
barplot(t(as.matrix(tbl)), col=brewer.pal(6, "Dark2"), xlab="Individual", ylab="Ancestry", border=NA)

popmap <-  read.table("data/pop_maps/n1_popmap_colony_site_cluster.txt", sep = '\t', header = FALSE)
popmap <- popmap %>% rename("ind" = "V1", "site" = "V4", "nest" = "V2", "cluster" = "V3") %>% 
  filter(site == "JatunSacha")
head(popmap)

names <-n1.js@ind.names
names2 <- popmap$ind
difference <- setdiff(names2, names)

for(i in 1:length(popmap$ind)){
  for(j in 1:length(difference)){
    if(popmap$ind[i] == difference[j]){
      popmap <- popmap[-i,]
    }
  }
}

f1.js.df <- cbind(popmap, tbl)

long.f1.js <- f1.js.df %>% pivot_longer(cols = V1:V4, names_to = "Kclust", values_to = "fraction", names_prefix = "V")

# trim column names
long.f1.js <- long.f1.js %>% 
  mutate(nest = str_sub(nest,0,5)) %>% 
  mutate(ind.num = str_sub(ind, -2)) 

# we have numbers 0-9 as 0.1-0.9 because of naming convention
# tail(long.dat %>% mutate(ind.num = str_replace(ind.num, ".", " ")))  # this doesn't work, turns 14 to 4

long.f1.js <- long.f1.js %>% mutate(ind.num2 = as.numeric(ind.num))

long.f1.js$ind.num2 <- if_else(long.f1.js$ind.num2 < 1.0, true = (long.f1.js$ind.num2*10), false = long.f1.js$ind.num2)

# also, 06A and 06B naming is causing issues
rows <- which(str_detect(long.f1.js$ind, "06A"))
long.f1.js$ind.num2[rows] <- paste(long.f1.js$ind.num2[rows], "A", sep = "")

rows <- which(str_detect(long.f1.js$ind, "06B"))
long.f1.js$ind.num2[rows] <- paste(long.f1.js$ind.num2[rows], "B", sep = "")

long.f1.js <- long.f1.js %>% mutate(ind.num2 = as.character(ind.num2))

f1.js.admix <-  ggplot(long.f1.js, aes(x=ind.num2, y=fraction, fill=Kclust)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(. ~ nest, drop=TRUE, space="free", scales="free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(fill = "K cluster", y = "Ancestry", x = "Individual")
f1.js.admix


f1.admix.plots <- ggpubr::ggarrange(f1.admix, f1.js.admix, nrow=2, labels = "auto")
ggsave(f1.admix.plots, filename = "../../figures/f1.admix.plots.jpeg", dpi = "retina",
       units = "in", width = 10, height = 8)

# Faiditus sp. 2 ---------------------------

vcf.f2 <- read.vcfR("data/filtered_vcfs/f2_minDP3.maxmiss.40.70.snps.recode.vcf")
pop.data.f2 <- read.table("data/pop_maps/f2_popmap_colony_site.txt", sep = '\t', header = FALSE)
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

# 
# saveRDS(gl.f2, file="rds/gl.f2.RDS")
write.csv(pop.data.cleaned.f2, "data/output_files/F2.pop.data.cleaned.csv")
# 
gl.f2 <- readRDS("rds/gl.f2.RDS")
# basic.stats(gl.f2)
# # remove NAs
toRemove <- is.na(glMean(gl.f2, alleleAsUnit = FALSE))
which(toRemove)
gl.f2 <- gl.f2[, !toRemove]

setPop(gl.f2) <- ~Nest
popNames(gl.f2) %>% length()
indNames(gl.f2) %>% str_subset("Ex08") %>% length()
indNames(gl.f2) %>% length()

####### Weir & Cockrham F-stats ########

mean(gl.Ho(gl.f2))
mean(gl.He(gl.f2))

library(dartR)
library(hierfstat)

# have to use genind object
gi.f2 <- gl2gi(gl.f2)
strata(gi.f2) <- strata_df_f2
#gi.f2 <- missingno(gi.f2)

saveRDS(gi.f2, file="rds/gi.f2.RDS")
gi.f2 <- readRDS("rds/gi.f2.RDS")
strata(gi.f2) <- strata_df_f2

# Pop as Nest
setPop(gi.f2) <- ~Nest
pop(gi.f2)
wc(gi.f2)



bs <- basic.stats(gi.f2)
pop.freq <- bs$pop.freq

f2.rich <- PopGenReport::allel.rich(gi.f2)$mean.richness
freq <- makefreq(gi.f2)

#### AMOVA #####

f2.amova <- poppr.amova(gi.f2, ~Nest)
f2.test.amova <- randtest(f2.amova)
plot(f2.test.amova)


###### IBD #########

## set up gps coords for eximius genind 
exi.coords <- read.csv("data/eximius_coordinates.csv")
exi.coords$Nest <- exi.coords$Nest %>% str_replace_all("Ex", "Exi")
exi.coords$Nest <- exi.coords$Nest %>% str_replace_all("Exii", "Exi")

strata_df_f2$lat <- 0
strata_df_f2$long <- 0

for(i in 1:length(strata_df_f2$Nest)){
  for(j in 1:length(exi.coords$Nest)){
    if(strata_df_f2$Nest[i] == exi.coords$Nest[j]){
      strata_df_f2$lat[i] <- exi.coords$latitude[j]
      strata_df_f2$long[i] <- exi.coords$longitude[j]
    }}}

## all archidona colonies share one set of coordinates (all within ~20 meters)
row <- which(exi.coords$Nest == "Exi_Archi_11.8")
archi.lat <- exi.coords$latitude[row]
archi.long <- exi.coords$longitude[row]

rows <- which(strata_df_f2$Site == "Archidona")
strata_df_f2$lat[rows] <- archi.lat
strata_df_f2$long[rows] <- archi.long

xy <- strata_df_f2 %>% dplyr::select(long, lat)
gi.f2@other$xy <- xy
gi.f2@other$xy



Dgen <- poppr::provesti.dist(gi.f2)
Dgeo <- dist(gi.f2@other$xy, method = "euclidean")

library(ade4)
ibd <- mantel.randtest(Dgen, Dgeo)
ibd

plot(ibd, main = "Mantel's test")


#### admixture plots ####

CVs <- read.table("data/admixture/f2/CV.csv", sep = " ")
CVs <- CVs[, 3:4] ## drop the first two columns
## Remove the formatting around the K values:
CVs[, 1] <- gsub(x = CVs[, 1], pattern = "\\(K=",
                 replacement = "")
CVs[, 1] <- gsub(x = CVs[, 1], pattern = "\\):",
                 replacement = "") 
head(CVs)


plot(CVs, xlab = "K", ylab = "Faiditus sp.1 CV error")
abline(v = c(2,11))



## using K = 2 based off BIC

tbl=read.table("data/admixture/f2/f2.scaffolds.2.Q")
barplot(t(as.matrix(tbl)), col=brewer.pal(6, "Dark2"), xlab="Individual", ylab="Ancestry", border=NA)

popmap <-  read.table("data/pop_maps/F2_popmap_colony_site.txt", sep = '\t', header = FALSE)
popmap <- popmap %>% rename("ind" = "V1", "nest" = "V2") %>% mutate(site = "JatunSacha")
head(popmap)

names <-gl.f2@ind.names
names2 <- popmap$ind
difference <- setdiff(names2, names)

for(i in 1:length(popmap$ind)){
  for(j in 1:length(difference)){
    if(popmap$ind[i] == difference[j]){
      popmap <- popmap[-i,]
    }
  }
}

f2.df <- cbind(popmap, tbl)

long.dat <- f2.df %>% pivot_longer(cols = V1:V2, names_to = "Kclust", values_to = "fraction", names_prefix = "V")


p <-  ggplot(long.dat, aes(x=ind, y=fraction, fill=Kclust)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(. ~ site, drop=TRUE, space="free", scales="free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(fill = "K cluster", y = "Ancestry", x = "Individual")
p
