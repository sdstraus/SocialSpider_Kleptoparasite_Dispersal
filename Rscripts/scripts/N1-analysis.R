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
# all(colnames(vcf@gt)[-1] == pop.data$AccessID)
# 
# new.list <- as.data.frame(as.factor(tidy.vcf.gt$Indiv)) # list of individuals
# colnames(new.list)[1] <- "Ind"

# cleaned <- levels(new.list[,1])
# dirty <- levels(as.factor(pop.data$V1))
# d.in.c <- dirty %in% cleaned
# rows <- which(d.in.c == TRUE)
# pop.data.cleaned <- pop.data[rows,]


#convert to genlight object
# gl.n1 <- vcfR2genlight(vcf)
# ploidy(gl.n1) <- 2
# pop(gl.n1) <- pop.data.cleaned$V2


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
  ggplot(aes(x = LD1, y = LD2, color = Site, shape = Group)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values=c(my_pal)) +
  scale_shape_manual(values=c(16, 17, 23, 3, 7))+
  # scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))+
  # geom_text(aes(label=names), nudge_x=0.3, nudge_y=0.5)+
  stat_ellipse(inherit.aes = FALSE, mapping = aes(x = LD1, y = LD2, shape = Group), 
               level = 0.95, size = 0.5)+
  labs(shape = "Group")+
  guides(shape=guide_legend(ncol=2), color = guide_legend(override.aes = list(shape =  15)))+
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



### minimum spanning network ####

euc.n1 <- bitwise.dist(gl.n1, euclidean = TRUE)
n1.msn <- poppr.msn(gl.n1, euc.n1, showplot = FALSE, include.ties= TRUE)
node.size <- rep(2, times = nInd(gl.n1))
names(node.size) <- indNames(gl.n1)
igraph::vertex.attributes(n1.msn$graph)$size <- node.size
set.seed(9)
plot_poppr_msn(gl.n1, n1.msn, palette = RColorBrewer::brewer.pal(n = nPop(gl.n1), 
                                                                 name = "Dark2"), gadj = 70, inds = "ALL")



##### compoplot
n1.dapc <- dapc(gl.n1, n.pca = 3, n.da = 2, glPca = n1.pca)

scatter(n1.dapc)
compoplot(n1.dapc)

dapc.results <- as.data.frame(n1.dapc$posterior)
dapc.results$pop <- pop(gl.n1)
dapc.results$indNames <- rownames(dapc.results)


dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))
head(dapc.results, n = 6)
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
write.csv(dapc.results, "dapc.results.csv")

p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop)) +
  geom_bar(stat='identity') + facet_grid(~Original_Pop, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p
n1.compoplot <- p
ggsave(n1.compoplot, filename = "figures/n1.compoplot.jpeg", dpi = "retina",
       height = 6, width = 8, units = 'in')


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

p2 <- my_df %>% 
  rename(Nest = V2) %>% 
  ggplot(aes(x = LD1, y = LD2, color = Nest, shape = Group)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values=c(my_pal)) +
  scale_shape_manual(values=c(16, 17, 23))+
  scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))+
  # geom_text(aes(label=names), nudge_x=0.3, nudge_y=0.5)+
  stat_ellipse(inherit.aes = FALSE, mapping = aes(x = LD1, y = LD2, shape = Group), 
               level = 0.95, size = 0.5)+
  labs(shape = "Group")+
  guides(shape=guide_legend(ncol=2), color = guide_legend(override.aes = list(shape =  15)))+
  cowplot::theme_cowplot()
p2

n1.js.plot <- p2
N1.Plot <- ggpubr::ggarrange(n1.all.plot, n1.js.plot, labels = "auto", nrow=2)

ggsave(N1.Plot, filename="../../figures/N1.pca.JS.jpeg", dpi = "retina", units = "in", 
       height = 9, width = 4)

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



### minimum spanning network ####

euc.n1.js <- bitwise.dist(n1.js, euclidean = TRUE)
n1.msn.js <- poppr.msn(n1.js, euc.n1.js, showplot = FALSE, include.ties= TRUE)
node.size <- rep(2, times = nInd(n1.js))
names(node.size) <- indNames(n1.js)
igraph::vertex.attributes(n1.msn.js$graph)$size <- node.size
set.seed(9)
plot_poppr_msn(n1.js, n1.msn.js, palette = RColorBrewer::brewer.pal(n = nPop(n1.js), name = "Dark2"), gadj = 70)




# N1 - archi ####

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

pop.cleaned.n1.archi <- pop.data.cleaned.n1 %>% filter(V4 == "Archidona")
# 
# 
# detectCores()
# n1.pca.archi <- glPca(n1.archi, nf = 3, parallel = TRUE, n.core = 4)
# n1.pca.archi.scores <- as.data.frame(n1.pca.archi$scores)
# n1.pca.archi.scores$pop <- pop(n1.archi)
# 
# saveRDS(n1.pca.archi, "../rds/n1.pca.archi.rds")
n1.pca.archi <- readRDS("../rds/n1.pca.archi.rds")

############## K- means clustering ##############
maxK <- 5
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(n1.archi, n.pca = 40, choose.n.clust = FALSE,  
                       max.n.clust = maxK, glPca = n1.pca.archi)
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
  grp_l[[i]] <- find.clusters(n1.archi, n.pca = 40, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(n1.archi, pop = grp_l[[i]]$grp, n.pca = 40, n.da = my_k[i], glPca = n1.pca.archi)
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}


my_df <- as.data.frame(dapc_l[[ 1 ]]$ind.coord)
my_df$Group <- dapc_l[[ 1 ]]$grp
head(my_df)

names <- rownames(my_df)
my_df <- cbind(names, my_df)
my_df <- dplyr::full_join(my_df, pop.cleaned.n1.archi, by = c("names" = "V1"))


my_pal <- RColorBrewer::brewer.pal(n=6, name = "Dark2")

p2 <- ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group, shape = V2)) +
  geom_point(size =2) +
  theme_bw() +
  scale_color_manual(values=c(my_pal)) +
  scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))+
  # geom_text(aes(label=names), nudge_x=0.3, nudge_y=0.5)+
  stat_ellipse(inherit.aes = FALSE, mapping = aes(x = LD1, y = LD2, color = Group), 
               level = 0.95, size = 1)+
  labs(shape = "Nest")
p2

N1.Plot <- p2
ggsave(N1.Plot, filename="figures/N1.pca.JS.jpeg", dpi = "retina", units = "in", 
       height = 4.5, width = 6)

# pop.data.cleaned.n1 <- read.csv("../data/N1.pop.data.cleaned.csv")
pop.data.cleaned.js <- pop.data.cleaned.n1 %>% filter(V4 == "Archidona")

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



### minimum spanning network ####

euc.n1.archi <- bitwise.dist(n1.archi, euclidean = TRUE)
n1.msn.archi <- poppr.msn(n1.archi, euc.n1.archi, showplot = FALSE, include.ties= TRUE)
node.size <- rep(2, times = nInd(n1.archi))
names(node.size) <- indNames(n1.archi)
igraph::vertex.attributes(n1.msn.archi$graph)$size <- node.size
set.seed(9)
plot_poppr_msn(n1.archi, n1.msn.archi, palette = RColorBrewer::brewer.pal(n = nPop(n1.archi), name = "Dark2"), gadj = 70)


# VL ####

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

pop.cleaned.n1.vl <- pop.data.cleaned.n1 %>% filter(V4 == "ViaLoreto")
# 
# 
# detectCores()
# n1.pca.vl <- glPca(n1.vl, nf = 3, parallel = TRUE, n.core = 4)
# n1.pca.vl.scores <- as.data.frame(n1.pca.archi$scores)
# n1.pca.vl.scores$pop <- pop(n1.vl)
# 
# saveRDS(n1.pca.vl, "../rds/n1.pca.vl.rds")
n1.pca.vl <- readRDS("../rds/n1.pca.vl.rds")

############## K- means clustering ##############
maxK <- 5
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(n1.vl, n.pca = 40, choose.n.clust = FALSE,  
                       max.n.clust = maxK, glPca = n1.pca.vl)
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
  grp_l[[i]] <- find.clusters(n1.vl, n.pca = 40, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(n1.vl, pop = grp_l[[i]]$grp, n.pca = 40, n.da = my_k[i], glPca = n1.pca.vl)
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}


my_df <- as.data.frame(dapc_l[[ 1 ]]$ind.coord)
my_df$Group <- dapc_l[[ 1 ]]$grp
head(my_df)

names <- rownames(my_df)
my_df <- cbind(names, my_df)
my_df <- full_join(my_df, pop.cleaned.n1.vl, by = c("names" = "V1"))

my_pal <- RColorBrewer::brewer.pal(n=6, name = "Dark2")

p2 <- ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group, shape = V2)) +
  geom_point(size = 2) +
  theme_bw() +
  scale_color_manual(values=c(my_pal)) +
  scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))+
  # geom_text(aes(label=names), nudge_x=0.3, nudge_y=0.5)+
  stat_ellipse(inherit.aes = FALSE, mapping = aes(x = LD1, y = LD2, color = Group), 
               level = 0.95, size = 1)+
  labs(shape = "Nest")
p2

N1.Plot <- p2
ggsave(N1.Plot, filename="figures/N1.pca.JS.jpeg", dpi = "retina", units = "in", 
       height = 4.5, width = 6)

# pop.data.cleaned.n1 <- read.csv("../data/N1.pop.data.cleaned.csv")
pop.data.cleaned.vl <- pop.data.cleaned.n1 %>% filter(V4 == "ViaLoreto")

tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]
tmp$Isolate <- rownames(tmp)
tmp <- melt(tmp, id = c("Isolate", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$Region <- pop.data.cleaned.vl$V2
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



### minimum spanning network ####

euc.n1.vl <- bitwise.dist(n1.vl, euclidean = TRUE)
n1.msn.vl <- poppr.msn(n1.vl, euc.n1.vl, showplot = FALSE, include.ties= TRUE)
node.size <- rep(2, times = nInd(n1.vl))
names(node.size) <- indNames(n1.vl)
igraph::vertex.attributes(n1.msn.vl$graph)$size <- node.size
set.seed(9)
plot_poppr_msn(n1.vl, n1.msn.vl, palette = RColorBrewer::brewer.pal(n = nPop(n1.vl), name = "Dark2"), gadj = 70)



n1.pca.plots <- ggpubr::ggarrange(n1.all.plot, n1.js.plot, labels = "auto")
ggsave(n1.pca.plots, filename = "../../figures/n1.plot.pcas.jpeg", dpi = 'retina', units = "in",
       width = 7, height = 3)
