# sam's laptop
setwd("/Users/samstraus/ownCloud/PhDThesis/GBS/Rscripts/scripts")

library(vcfR)
library(tidyr)
library(adegenet)
library(SNPRelate)
library(ggplot2)
library(parallel)
library(poppr)

# Eximius ########

gl.exi <- readRDS("../rds/gl.exi.RDS")
pop.data.cleaned <- read.csv("../data/eximius.pop.data.cleaned.csv")

##### exi all - pop ~ nest ######
library(parallel)
detectCores()


setPop(gl.exi) <- ~Nest
exi.pca.nest <- glPca(gl.exi, nf = 3, parallel = TRUE, n.core = 4)

saveRDS(exi.pca.nest, "../rds/exi.pca.nest.rds")
exi.pca.nest <- readRDS("../rds/exi.pca.nest.rds")

exi.pca.nest.scores <- as.data.frame(exi.pca.nest$scores)
exi.pca.nest.scores$pop <- unlist(strata(gl.exi, ~Nest))

maxK <- 10
myMat <- matrix(nrow=10, ncol=(maxK))
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(gl.exi, n.pca = 40, choose.n.clust = FALSE,  
                       max.n.clust = maxK, glPca = exi.pca.nest)
  myMat[i,] <- grp$Kstat
}

library(ggplot2)
library(reshape2)
my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)

p1 <- ggplot(my_df, aes(x = K, y = BIC))+
  geom_boxplot()+
  theme_bw()+
  xlab("Number of groups (K)")
 p1
# exi.elbow <- p1
# ggsave(exi.elbow, filename = "figures/exi.elbow.jpeg", dpi = "retina",
#        units = "in", height = 4.5, width = 6)

# DAPC
my_k <- c(5)
grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] <- find.clusters(gl.exi, n.pca = 40, n.clust = my_k[i], glPca = exi.pca.nest)
  dapc_l[[i]] <- dapc(gl.exi, pop = grp_l[[i]]$grp, n.pca = 40, n.da = my_k[i], glPca = exi.pca.nest)
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}


my_df <- as.data.frame(dapc_l[[ 1 ]]$ind.coord)
my_df$Group <- dapc_l[[ 1 ]]$grp
head(my_df)

names <- rownames(my_df)
my_df <- cbind(names, my_df)
my_df <- dplyr::full_join(my_df, pop.data.cleaned, by = c("names"="V1"))

my_pal <- RColorBrewer::brewer.pal(n=6, name = "Dark2")

p2 <- my_df %>% 
  rename(Site = V2) %>% 
  ggplot( aes(x = LD1, y = LD2, color = Site)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values=c(my_pal)) +
  # stat_ellipse(level = 0.95, size = 0.5) +
  scale_shape_manual(values=c(16, 17, 23, 3, 7))+
  scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))+
  # geom_text(aes(label=names), nudge_x=0.3, nudge_y=0.6)+
  stat_ellipse(inherit.aes = FALSE, mapping = aes(x = LD1, y = LD2, shape = Group), 
               level = 0.95, size = 0.5)+
  # guides(shape=guide_legend(ncol=2), color = guide_legend(override.aes = list(shape =  15)))+
  cowplot::theme_cowplot()
  
p2
exi.all.plot <- p2

tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]
tmp$Isolate <- rownames(tmp)
tmp <- melt(tmp, id = c("Isolate", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$Region <- pop.data.cleaned$V2
my_df <- tmp

for(i in 1:length(dapc_l)){
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
  ylab("Posterior membership probability")+
  theme(legend.position='none') +
  #p3 <- p3 + scale_color_brewer(palette="Dark2")
  scale_fill_manual(values=c(my_pal)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p3



#### exi - all - pop ~ site #####
setPop(gl.exi) <- ~Site

euc.exi <- bitwise.dist(gl.exi, euclidean = TRUE)
exi.msn <- poppr.msn(gl.exi, euc.exi, showplot = FALSE, include.ties= TRUE)
node.size <- rep(2, times = nInd(gl.exi))
# setPop(gl.exi) <- ~Nest
# names(node.size) <- pop(gl.exi)
# setPop(gl.exi) <- ~Site
igraph::vertex.attributes(exi.msn$graph)$size <- node.size
set.seed(9)
exi.msn.plot <- plot_poppr_msn(gl.exi, exi.msn, 
               palette = RColorBrewer::brewer.pal(n = nPop(gl.exi), 
                                                  name = "Dark2"), gadj = 70, inds = "ALL")

ggsave(exi.msn.plot, filename = "exi.msn.plot.jpeg")

##### exi - js ######
setPop(gl.exi) <- ~Site
popNames(gl.exi)
exi.js <- popsub(gl.exi, "JatunSacha")
popNames(exi.js)
setPop(exi.js) <- ~Nest
indNames(exi.js)

pop.data.cleaned.js <- strata_df %>% filter(Site == "JatunSacha")

setPop(exi.js) <- ~Nest
js.pca.nest <- glPca(exi.js, nf = 3, parallel = TRUE, n.core = 4)
# 
saveRDS(js.pca.nest, "../rds/js.pca.nest.rds")
js.pca.nest <- readRDS("../rds/js.pca.nest.rds")

js.pca.nest.scores <- as.data.frame(js.pca.nest$scores)
js.pca.nest.scores$pop <- unlist(strata(exi.js, ~Nest))

maxK <- 10
myMat <- matrix(nrow=10, ncol=(maxK))
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(exi.js, n.pca = 40, choose.n.clust = FALSE,  
                       max.n.clust = maxK, glPca = js.pca.nest)
  myMat[i,] <- grp$Kstat
}

library(ggplot2)
library(reshape2)
my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)

p1 <- ggplot(my_df, aes(x = K, y = BIC))+
  geom_boxplot()+
  theme_bw()+
  xlab("Number of groups (K)")
p1
# exi.elbow <- p1
# ggsave(exi.elbow, filename = "figures/exi.elbow.jpeg", dpi = "retina",
#        units = "in", height = 4.5, width = 6)

# DAPC
my_k <- c(6)
grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] <- find.clusters(exi.js, n.pca = 40, n.clust = my_k[i], glPca = js.pca.nest)
  dapc_l[[i]] <- dapc(exi.js, pop = grp_l[[i]]$grp, n.pca = 40, n.da = my_k[i], glPca = js.pca.nest)
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}


my_df <- as.data.frame(dapc_l[[ 1 ]]$ind.coord)
my_df$Group <- dapc_l[[ 1 ]]$grp
head(my_df)

names <- rownames(my_df)
my_df <- cbind(names, my_df)
my_df <- dplyr::full_join(my_df, pop.data.cleaned.js, by = c("names"="V1"))

my_pal <- RColorBrewer::brewer.pal(n=12, name = "Paired")

my_df <- my_df %>% mutate(Nest = str_remove(Nest, "_JS"))

# reorder nest levels
my_df <- my_df %>% mutate(Nest = as.factor(Nest))
levels(my_df$Nest)

my_df$Nest <- factor(my_df$Nest, levels = c("Exi06A", "Exi06B", "Exi07", 
                                                  "Exi11", "Exi12", "Exi02", "Exi10", "Exi05", "Exi08", 
                                                  "Exi14"))


p2 <- my_df %>% 
  # rename(Nest = V2) %>% 
  ggplot(aes(x = LD1, y = LD2, color = Nest)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values=c(my_pal)) +
  scale_shape_manual(values=c(16, 17, 23, 3, 7, 8))+
  # scale_shape_manual(values = c(0,1,2,3,4,5,7,8,15,16,17))+
  # stat_ellipse(level = 0.95, size = 0.5) +
  scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))+
  # geom_text(aes(label=names), nudge_x=0.3, nudge_y=0.6)+
  stat_ellipse(inherit.aes = FALSE, mapping = aes(x = LD1, y = LD2, shape = Group),
               level = 0.95, size = 0.5)+
  # guides(shape=guide_legend(ncol=2), fill = guide_legend(ncol=2), color = guide_legend(ncol=1, 
  #                                                             override.aes = list(shape =  15)))+
  # labs(shape = "Group")+
  cowplot::theme_cowplot()
p2
exi.js.plot <- p2

tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]
tmp$Isolate <- rownames(tmp)
tmp <- melt(tmp, id = c("Isolate", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$Region <- pop.data.cleaned.js$V2
my_df <- tmp


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
  ylab("Posterior membership probability")+
  theme(legend.position='none') +
  #p3 <- p3 + scale_color_brewer(palette="Dark2")
  scale_fill_manual(values=c(my_pal)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p3


## minimum spanning network ###

euc.js <- bitwise.dist(exi.js, euclidean = TRUE)
js.msn <- poppr.msn(exi.js, euc.js, showplot = FALSE, include.ties= TRUE)
node.size <- rep(2, times = nInd(exi.js))
names(node.size) <- indNames(exi.js)
igraph::vertex.attributes(js.msn$graph)$size <- node.size
set.seed(9)
plot_poppr_msn(exi.js, js.msn, palette = RColorBrewer::brewer.pal(n = nPop(exi.js), name = "Dark2"), gadj = 70)



##### exi - archi #######
exi.archi <- popsub(gl.exi, "Archidona")
setPop(exi.archi) <- ~Nest
popNames(exi.archi)
indNames(exi.archi)
# remove NAs
toRemove <- is.na(glMean(exi.archi, alleleAsUnit = FALSE))
which(toRemove)
exi.archi <- exi.archi[, !toRemove]


setPop(exi.archi) <- ~Nest
archi.pca.nest <- glPca(exi.archi, nf = 3, parallel = TRUE, n.core = 4)

pop.data.cleaned.archi <- strata_df %>% filter(Site == "Archidona")

saveRDS(archi.pca.nest, "../rds/archi.pca.nest.rds")
archi.pca.nest <- readRDS("../rds/archi.pca.nest.rds")

archi.pca.nest.scores <- as.data.frame(archi.pca.nest$scores)
archi.pca.nest.scores$pop <- unlist(strata(exi.archi, ~Nest))

maxK <- 10
myMat <- matrix(nrow=10, ncol=(maxK))
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(exi.archi, n.pca = 40, choose.n.clust = FALSE,  
                       max.n.clust = maxK, glPca = archi.pca.nest)
  myMat[i,] <- grp$Kstat
}

library(ggplot2)
library(reshape2)
my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)

p1 <- ggplot(my_df, aes(x = K, y = BIC))+
  geom_boxplot()+
  theme_bw()+
  xlab("Number of groups (K)")
p1
# exi.elbow <- p1
# ggsave(exi.elbow, filename = "figures/exi.elbow.jpeg", dpi = "retina",
#        units = "in", height = 4.5, width = 6)

# DAPC
my_k <- c(3)
grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] <- find.clusters(exi.archi, n.pca = 40, n.clust = my_k[i], glPca = archi.pca.nest)
  dapc_l[[i]] <- dapc(exi.archi, pop = grp_l[[i]]$grp, n.pca = 40, n.da = my_k[i], glPca = archi.pca.nest)
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}


my_df <- as.data.frame(dapc_l[[ 1 ]]$ind.coord)
my_df$Group <- dapc_l[[ 1 ]]$grp
head(my_df)

names <- rownames(my_df)
my_df <- cbind(names, my_df)
my_df <- dplyr::full_join(my_df, pop.data.cleaned.archi, by = c("names"="V1"))

# trim "Exi_Archi" off nest names
my_df <- my_df %>% mutate(Nest = str_remove(Nest, "Exi_Archi_"))

my_pal <- RColorBrewer::brewer.pal(n=7, name = "Dark2")

p2 <- my_df %>% 
  # rename(Nest = V2) %>% 
  ggplot(aes(x = LD1, y = LD2, color = Nest)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values=c(my_pal)) +
  scale_shape_manual(values=c(16, 17, 23))+
  # stat_ellipse(level = 0.95, size = 1) +
  scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))+
  # geom_text(aes(label=names), nudge_x=0.3, nudge_y=0.6)+
  stat_ellipse(inherit.aes = FALSE, mapping = aes(x = LD1, y = LD2, shape = Group), 
               level = 0.95, size = 0.5)+
  # guides(shape=guide_legend(ncol=2, order = 1), color = guide_legend(override.aes = list(shape =  15)))+
  labs(shape = "Group")+
  cowplot::theme_cowplot()
p2
exi.archi.plot <- p2

tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]
tmp$Isolate <- rownames(tmp)
tmp <- melt(tmp, id = c("Isolate", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$Region <- pop.data.cleaned.archi$V2
my_df <- tmp


for(i in 1:length(dapc_l)){
  tmp <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K <- my_k[i]
  tmp$Isolate <- rownames(tmp)
  tmp <- melt(tmp, id = c("Isolate", "K"))
  names(tmp)[3:4] <- c("Group", "Posterior")
  tmp$Region <- pop.data.cleaned.archi$V2
  
  my_df <- rbind(my_df, tmp)
}

grp.labs <- paste("K =", my_k)
names(grp.labs) <- my_k

p3 <- ggplot(my_df, aes(x = Isolate, y = Posterior, fill = Group)) +
  geom_bar(stat = "identity") +
  facet_grid(K ~ Region, scales = "free_x", space = "free", 
             labeller = labeller(K = grp.labs)) +
  theme_bw() +
  ylab("Posterior membership probability")+
  theme(legend.position='none') +
  #p3 <- p3 + scale_color_brewer(palette="Dark2")
  scale_fill_manual(values=c(my_pal)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p3


euc.archi <- bitwise.dist(exi.archi, euclidean = TRUE)
archi.msn <- poppr.msn(exi.archi, euc.archi, showplot = FALSE, include.ties= TRUE)
node.size <- rep(2, times = nInd(exi.archi))
names(node.size) <- indNames(exi.archi)
igraph::vertex.attributes(archi.msn$graph)$size <- node.size
set.seed(9)
plot_poppr_msn(exi.archi, archi.msn, palette = RColorBrewer::brewer.pal(n = nPop(exi.archi), name = "Dark2"), gadj = 70)



#### exi - vl #####

exi.vl <- popsub(gl.exi, sublist=c("ViaLoreto"))
indNames(exi.vl)
# remove NAs
toRemove <- is.na(glMean(exi.vl, alleleAsUnit = FALSE))
which(toRemove)
exi.vl <- exi.vl[, !toRemove]
setPop(exi.vl) <- ~Nest
popNames(exi.vl)

# keep
setPop(exi.vl) <- ~Nest
vl.pca.nest <- glPca(exi.vl, nf = 3, parallel = TRUE, n.core = 4)
# 
# 
saveRDS(vl.pca.nest, "../rds/vl.pca.nest.rds")

vl.pca.nest <- readRDS("../rds/vl.pca.nest.rds")
pop.data.cleaned.vl <- strata_df %>% filter(Site == "ViaLoreto")

vl.pca.nest.scores <- as.data.frame(vl.pca.nest$scores)
vl.pca.nest.scores$pop <- unlist(strata(exi.vl, ~Nest))

maxK <- 10
myMat <- matrix(nrow=10, ncol=(maxK))
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(exi.vl, n.pca = 40, choose.n.clust = FALSE,  
                       max.n.clust = maxK, glPca = vl.pca.nest)
  myMat[i,] <- grp$Kstat
}

library(ggplot2)
library(reshape2)
my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)

p1 <- ggplot(my_df, aes(x = K, y = BIC))+
  geom_boxplot()+
  theme_bw()+
  xlab("Number of groups (K)")
p1
# exi.elbow <- p1
# ggsave(exi.elbow, filename = "figures/exi.elbow.jpeg", dpi = "retina",
#        units = "in", height = 4.5, width = 6)

# DAPC
my_k <- c(4)
grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] <- find.clusters(exi.vl, n.pca = 40, n.clust = my_k[i], glPca = vl.pca.nest)
  dapc_l[[i]] <- dapc(exi.vl, pop = grp_l[[i]]$grp, n.pca = 40, n.da = my_k[i], glPca = vl.pca.nest)
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}


my_df <- as.data.frame(dapc_l[[ 1 ]]$ind.coord)
my_df$Group <- dapc_l[[ 1 ]]$grp
head(my_df)

names <- rownames(my_df)
my_df <- cbind(names, my_df)
my_df <- dplyr::full_join(my_df, pop.data.cleaned.vl, by = c("names" = "V1"))

# trim "Exi_VL" off nest names
my_df <- my_df %>% mutate(Nest = str_remove(Nest, "Exi_VL_"))

# reorder nest levels
my_df$Nest <- as.factor(my_df$Nest)
levels(my_df$Nest )
# reorder levels
my_df$Nest <- factor(my_df$Nest, levels = c("9.5.1", "10.5", "17.5","24.1", 
                                                        "Bridge.1", "41.2"))

my_pal <- RColorBrewer::brewer.pal(n=11, name = "Paired")

p2 <- my_df %>% 
  # rename(Nest = V2) %>% 
  ggplot(aes(x = LD1, y = LD2, color = Nest)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values=c(my_pal)) +
  scale_shape_manual(values=c(16, 17, 23, 3))+
  # scale_shape_manual(values = c(0,1,2,3,4,7,8,15,16))+
  # stat_ellipse(level = 0.95, size = 1) +
  scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))+
  # geom_text(aes(label=names), nudge_x=0.3, nudge_y=0.6)+
  stat_ellipse(inherit.aes = FALSE, aes(x = LD1, y = LD2, shape = Group), 
               level = 0.95, size = 0.5)+
  labs(shape = "Group")+
  # guides(shape=guide_legend(ncol=2, order = 1), color = guide_legend(override.aes = list(shape =  15)))+
  cowplot::theme_cowplot()
p2
exi.vl.plot <- p2


exi.plot <- ggpubr::ggarrange(exi.all.plot, exi.js.plot, exi.archi.plot, exi.vl.plot, labels = "auto",nrow = 4)
ggsave(exi.plot, filename = "../../figures/exi.plot.pcas.jpeg", dpi = 'retina', units = "in",
       width = 7, height = 12)

## trying to make the widths all the same
library(gtable)
library(grid)
library(gridExtra)

# Get the gtables
gA <- ggplotGrob(exi.all.plot)
gB <- ggplotGrob(exi.js.plot)
gC <- ggplotGrob(exi.archi.plot)
gD <- ggplotGrob(exi.vl.plot)

gC$widths <- gA$widths
gB$widths <- gA$widths
gD$widths <- gA$widths
# Arrange the two charts.
# The legend boxes are centered
grid.newpage()
jpeg(file = "../../figures/Figure3_02Sept22.jpeg", width = 6, height = 9, units = "in", res = 1000)
grid.arrange(gA, gB, gC, gD, nrow = 4, padding = 2)
dev.off()
        
tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]
tmp$Isolate <- rownames(tmp)
tmp <- melt(tmp, id = c("Isolate", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$Region <- pop.data.cleaned.vl$V2
my_df <- tmp


for(i in 1:length(dapc_l)){
  tmp <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K <- my_k[i]
  tmp$Isolate <- rownames(tmp)
  tmp <- melt(tmp, id = c("Isolate", "K"))
  names(tmp)[3:4] <- c("Group", "Posterior")
  tmp$Region <- pop.data.cleaned.vl$V2
  
  my_df <- rbind(my_df, tmp)
}

grp.labs <- paste("K =", my_k)
names(grp.labs) <- my_k

p3 <- ggplot(my_df, aes(x = Isolate, y = Posterior, fill = Group)) +
  geom_bar(stat = "identity") +
  facet_grid(K ~ Region, scales = "free_x", space = "free", 
             labeller = labeller(K = grp.labs)) +
  theme_bw() +
  ylab("Posterior membership probability")+
  theme(legend.position='none') +
  #p3 <- p3 + scale_color_brewer(palette="Dark2")
  scale_fill_manual(values=c(my_pal)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p3


euc.vl <- bitwise.dist(exi.vl, euclidean = TRUE)
vl.msn <- poppr.msn(exi.vl, euc.vl, showplot = FALSE, include.ties= TRUE)
node.size <- rep(2, times = nInd(exi.vl))
names(node.size) <- indNames(exi.vl)
igraph::vertex.attributes(vl.msn$graph)$size <- node.size
set.seed(9)
plot_poppr_msn(exi.vl, vl.msn, palette = RColorBrewer::brewer.pal(n = nPop(exi.vl), name = "Dark2"), gadj = 70)


