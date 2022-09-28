# sam's laptop
setwd("/Users/sam/PhDThesis/GBS/Rscripts/scripts")

library(vcfR)
library(tidyr)
library(adegenet)
library(SNPRelate)
library(ggplot2)
library(parallel)
library(poppr)

######## create genlight ########
# # read vcf
# vcf <- read.vcfR("filtering/n1/filtered-n1.vcf")

# # prune down pop.data to match inds in the vcf file
# don't need to redo 
# tidy.vcf <- vcfR2tidy(vcf)
# tidy.vcf.gt <- as.data.frame(tidy.vcf$gt) # extract genotype data

# read genlight
gl.n1 <- readRDS("../rds/gl.n1.RDS")

pop.data <- read.table("../../pop_maps/n1_popmap_colony_site_cluster.txt", sep = '\t', header = FALSE)


strata_df_n1 <- pop.data.cleaned[,-1]
strata_df_n1 <- strata_df %>%
  dplyr::rename(Nest = V2, Site = V4, Cluster = V3)

# strata(gl.n1) <- strata_df
# head(strata(gl.n1))

####### subset population ########


#JS
setPop(gl.n1) <- ~Site
popNames(gl.n1)
n1.js <- popsub(gl.n1, "JatunSacha")
popNames(n1.js)
setPop(n1.js) <- ~Nest
# remove NAs
toRemove <- is.na(glMean(n1.js, alleleAsUnit = FALSE))
which(toRemove)
n1.js <- n1.js[, !toRemove]


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

##### All sites ########

setPop(gl.n1) <- ~Nest
# 
# detectCores()
n1.pca <- glPca(gl.n1, nf = 3, parallel = TRUE, n.core = 4)

#
# saveRDS(n1.pca, "../rds/n1.pca.rds")

n1.pca <- readRDS("../rds/n1.pca.rds")
n1.pca.scores <- as.data.frame(n1.pca$scores)
n1.pca.scores$pop <- pop(gl.n1)

pop.data.cleaned.n1 <- read.csv("../data/N1.pop.data.cleaned.csv")


set.seed(9)
barplot(100*n1.pca$eig/sum(n1.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")

ggplot(n1.pca.scores, aes(x=PC1, y=PC2, colour=pop)) + 
  geom_point(size=2) +
  stat_ellipse(level = 0.95, size = 1) +
  theme_bw()


############## K- means clustering ##############
maxK <- 10
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(gl.n1, n.pca = 40, choose.n.clust = FALSE,  
                       max.n.clust = maxK, glPca = n1.pca, stat = "BIC")
  myMat[i,] <- grp$Kstat
}


library(reshape2)
my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)

p1 <- ggplot(my_df, aes(x = K, y = BIC)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Number of groups (K)")
p1

n1.elbow <- p1
ggsave(n1.elbow, filename = "figures/n1.elbow.jpeg", dpi = "retina",
       units = "in", height = 4.5, width = 6)

# DAPC
my_k <- c(5)

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] <- find.clusters(gl.n1, n.pca = 40, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(gl.n1, pop = grp_l[[i]]$grp, n.pca = 40, n.da = my_k[i], glPca = n1.pca)
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}


my_df <- as.data.frame(dapc_l[[ 1 ]]$ind.coord)
my_df$Group <- dapc_l[[ 1 ]]$grp
head(my_df)

names <- rownames(my_df)
my_df <- cbind(names, my_df)
my_df <- dplyr::full_join(my_df, pop.data.cleaned.n1, by = c("names" = "V1"))

my_pal <- RColorBrewer::brewer.pal(n=6, name = "Dark2")

p2 <- my_df %>% 
  rename(Site = V4) %>% 
  ggplot(aes(x = LD1, y = LD2, color = Site)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values=c(my_pal)) +
  scale_shape_manual(values=c(16, 17, 23, 3, 7))+
  # scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))+
  # geom_text(aes(label=names), nudge_x=0.3, nudge_y=0.5)+
  stat_ellipse(inherit.aes = FALSE, mapping = aes(x = LD1, y = LD2, shape = Group), 
               level = 0.95, size = 0.5)+
  labs(shape = "Group")+
  # guides(shape=guide_legend(ncol=2), color = guide_legend(override.aes = list(shape =  15)))+
  cowplot::theme_cowplot()
p2
n1.all.plot <- p2

# N1.Plot <- p2
# ggsave(N1.Plot, filename="figures/N1.pca.JS.jpeg", dpi = "retina", units = "in", 
#        height = 4.5, width = 6)

tmp <- as.data.frame(dapc_l[[2]]$posterior)
tmp$K <- my_k[1]
tmp$Isolate <- rownames(tmp)
tmp <- melt(tmp, id = c("Isolate", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$Region <- pop.data.cleaned$V2
my_df <- tmp


pop.data.cleaned.n1 <- read.csv("../data/N1.pop.data.cleaned.csv")

for(i in 2:length(dapc_l)){
  tmp <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K <- my_k[i]
  tmp$Isolate <- rownames(tmp)
  tmp <- melt(tmp, id = c("Isolate", "K"))
  names(tmp)[3:4] <- c("Group", "Posterior")
  tmp$Region <- pop.data.cleaned$V2
  
  my_df <- rbind(my_df, tmp)
}

grp.labs <- paste("K =", my_k)
names(grp.labs) <- my_k

p3 <- ggplot(my_df, aes(x = Isolate, y = Posterior, fill = Group)) +
  geom_bar(stat = "identity") +
  facet_grid(K ~ Region, scales = "free_x", space = "free", 
                      labeller = labeller(K = grp.labs)) +
  theme_bw() +
  ylab("Posterior membership probability") +
  theme(legend.position='none')+
#p3 <- p3 + scale_color_brewer(palette="Dark2")
  scale_fill_manual(values=c(my_pal)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p3



#### N1 - js #######
setPop(gl.n1) <- ~Site
popNames(gl.n1)
n1.js <- popsub(gl.n1, "JatunSacha")
popNames(n1.js)
setPop(n1.js) <- ~Nest
# remove NAs
toRemove <- is.na(glMean(n1.js, alleleAsUnit = FALSE))
which(toRemove)
n1.js <- n1.js[, !toRemove]

pop.cleaned.n1.js <- pop.data.cleaned.n1 %>% filter(V4 == "JatunSacha")

# 
# detectCores()
n1.pca.js <- glPca(n1.js, nf = 3, parallel = TRUE, n.core = 4)
n1.pca.js.scores <- as.data.frame(n1.pca.js$scores)
n1.pca.js.scores$pop <- pop(n1.js)

# saveRDS(n1.pca.js, "../rds/n1.pca.js.rds")

n1.pca.js <- readRDS("../rds/n1.pca.js.rds")

############## K- means clustering ##############
maxK <- 6
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(n1.js, n.pca = 40, choose.n.clust = FALSE,  
                       max.n.clust = maxK, glPca = n1.pca.js)
  myMat[i,] <- grp$Kstat
}


library(reshape2)
my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)

p1 <- ggplot(my_df, aes(x = K, y = BIC)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Number of groups (K)")
p1

# n1.elbow <- p1
# ggsave(n1.elbow, filename = "figures/n1.elbow.jpeg", dpi = "retina",
#        units = "in", height = 4.5, width = 6)

# DAPC
my_k <- c(3)

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] <- find.clusters(n1.js, n.pca = 40, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(n1.js, pop = grp_l[[i]]$grp, n.pca = 40, n.da = my_k[i], glPca = n1.pca.js)
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}


my_df <- as.data.frame(dapc_l[[ 1 ]]$ind.coord)
my_df$Group <- dapc_l[[ 1 ]]$grp
head(my_df)

names <- rownames(my_df)
my_df <- cbind(names, my_df)
my_df <- dplyr::full_join(my_df, pop.cleaned.n1.js, by = c("names" = "V1"))

my_pal <- RColorBrewer::brewer.pal(n=6, name = "Dark2")

# drop "_JS"
my_df <- my_df %>%  rename(Nest = V2)
my_df <- my_df %>% mutate(Nest = stringr::str_remove(Nest, "_JS"))

# reorder nest levels
my_df <- my_df %>% mutate(Nest = as.factor(Nest))
levels(my_df$Nest)

my_df$Nest <- factor(my_df$Nest, levels = c("Exi06B","Exi12", "Exi08", "Exi17"))
# reorder nest levels

p2 <- my_df %>% 
  ggplot(aes(x = LD1, y = LD2, color = Nest)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values=c(my_pal)) +
  scale_shape_manual(values=c(16, 17, 23))+
  scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))+
  # geom_text(aes(label=names), nudge_x=0.3, nudge_y=0.5)+
  stat_ellipse(inherit.aes = FALSE, mapping = aes(x = LD1, y = LD2, shape = Group), 
               level = 0.95, size = 0.5)+
  labs(shape = "Group")+
  # guides(shape=guide_legend(ncol=2), color = guide_legend(override.aes = list(shape =  15)))+
  cowplot::theme_cowplot()
p2

n1.js.plot <- p2
N1.Plot <- ggpubr::ggarrange(n1.all.plot, n1.js.plot, labels = "auto", nrow=2)

ggsave(N1.Plot, filename="../../figures/N1.pca.JS.jpeg", dpi = "retina", units = "in", 
       height = 9, width = 4)



## plotting with equal sizing of plots
library(gtable)
library(grid)
library(gridExtra)

# Get the gtables
gE <- ggplotGrob(n1.all.plot)
gF <- ggplotGrob(n1.js.plot)


gE$widths <- gA$widths ## same width as eximius plots
gF$widths <- gA$widths


# Arrange the two charts.
# The legend boxes are centered
grid.newpage()
jpeg(file = "../../figures/Figure4_02Sept22.jpeg", width = 6, height = 5, units = "in", res = 1000)
grid.arrange(gE, gF, nrow = 2, padding = 2)
dev.off()



pop.data.cleaned.n1 <- read.csv("../data/N1.pop.data.cleaned.csv")
pop.data.cleaned.js <- pop.data.cleaned.n1 %>% filter(V4 == "JatunSacha")

tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]
tmp$Isolate <- rownames(tmp)
tmp <- melt(tmp, id = c("Isolate", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$Region <- pop.data.cleaned.js$V2
my_df <- tmp


# pop.data.cleaned <- read.csv("../data/N1.pop.data.cleaned.csv")

for(i in 1:length(dapc_l)){
  tmp <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K <- my_k[i]
  tmp$Isolate <- rownames(tmp)
  tmp <- melt(tmp, id = c("Isolate", "K"))
  names(tmp)[3:4] <- c("Group", "Posterior")
  tmp$Region <- pop.data.cleaned.js$V2
  
  my_df <- rbind(my_df, tmp)
}

grp.labs <- paste("K =", my_k)
names(grp.labs) <- my_k

p3 <- ggplot(my_df, aes(x = Isolate, y = Posterior, fill = Group)) +
  geom_bar(stat = "identity") +
  facet_grid(K ~ Region, scales = "free_x", space = "free", 
             labeller = labeller(K = grp.labs)) +
  theme_bw() +
  ylab("Posterior membership probability") +
  theme(legend.position='none')+
  #p3 <- p3 + scale_color_brewer(palette="Dark2")
  scale_fill_manual(values=c(my_pal)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p3

