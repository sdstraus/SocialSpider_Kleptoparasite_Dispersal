setwd("/Users/sam/PhDThesis/GBS")
library(stringr)
# 
###### OLD #######
# #First we install SNPRelate
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# BiocManager::install(c("SNPRelate"))
# 
# #Then load them
# library(SNPRelate)
# library(tidyverse)
# 
# # convert our vcf into a gds as the first step.
# snpgdsVCF2GDS("denovo_exi_col_site/stacks_output/populations.snps.vcf",
#               "denovo_exi_col_site/stacks_output/populations.snps.gds",
#               method="biallelic.only")
# 
# #Load the gds file
# genofile <- snpgdsOpen("denovo_exi_col_site/stacks_output/populations.snps.gds")
# snpgdsSummary("denovo_exi_col_site/stacks_output/populations.snps.gds")
# 
# #Prune for linkage
# # takes ~15 min
# snpset_pruned <- snpgdsLDpruning(genofile, autosome.only=F)
# names(snpset_pruned)
# 
# #Make a list of sites we're keeping.
# snpset.id <- unlist(unname(snpset_pruned))
# 
# #Run the PCA
# pca <- snpgdsPCA(genofile, num.thread = 1, eigen.cnt = 16, snp.id = snpset.id ,autosome.only = F)
# 
# 
# #Here's the percent variance explained for each eigenvector
# pc.percent <- pca$varprop*100
# round(pc.percent, 2)
# 
# 
# #Make a dataframe of your PCA results
# tab <- data.frame(sample = pca$sample.id,
#                   pop = factor(pop_code)[match(pca$sample.id, sample.id)],
#                   PC1 = pca$eigenvect[,1],    # the first eigenvector
#                   PC2 = pca$eigenvect[,2],    # the second eigenvector
#                   stringsAsFactors = FALSE)
# 
# #Plot a PCA image
# tab %>%
#   ggplot(.,aes(x=PC1,y=PC2)) + geom_point()
# 
# 
# tab %>%
#   mutate(group = substr(sample,1,3)) %>%
#   ggplot(.,aes(x=PC1,y=PC2)) + 
#   geom_point(aes(color=group))
# 
# 
# 
# sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
# #pop_code <- read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))
# 
# pop_map <- read_tsv("pop_maps/eximius_popmap_colony_site.txt", col_names = FALSE)
# pop_code <- pop_map[2]$X2
# 
# #install.packages("diveRsity")
# library(diveRsity)
# #diffCalc(infile = "stacks/populations.snps.gds", outfile = 'results', fst = TRUE)

######## prep genlight ########
# 
# # convert vcf to bed file
# # install.packages("adegenet")
# # library(bedr)
# library(vcfR)
# library(adegenet)
# library(SNPRelate)
# library(tidyverse)
# 
# # vcf <- read.vcfR("denovo_exi_col_site/stacks_output/populations.snps.vcf")
# # bed <- vcf2bed(vcf, filename = "denovo_exi_col_site/stacks_output/populations.snps.bed")
# 
# # for eximius
# vcf <- read.vcfR("filtering/exi/filtered.recode.vcf")
# pop.data <- read.table("pop_maps/eximius_popmap_colony_site_cluster.txt", sep = '\t', header = FALSE)
# all(colnames(vcf@gt)[-1] == pop.data$AccessID)
# 
# # filtered list has fewer individuals
# tidy.filtered <- vcfR2tidy(vcf)
# tidy.filtered.gt <- as.data.frame(tidy.filtered$gt)
# 
# # max(tidy.filtered.gt$gt_DP, na.rm = TRUE)
# # min(tidy.filtered.gt$gt_DP, na.rm = TRUE)
# # mean(tidy.filtered.gt$gt_DP, na.rm = TRUE)
# 
# new.list <- as.data.frame(as.factor(tidy.filtered.gt$Indiv))
# colnames(new.list)[1] <- "Ind"
# 
# cleaned <- levels(new.list[,1])
# dirty <- levels(as.factor(pop.data$V1))
# d.in.c <- dirty %in% cleaned
# rows <- which(d.in.c == TRUE)
# pop.data.cleaned <- pop.data[rows,]
# 
# #convert to genlight object
# gl.exi <- vcfR2genlight(vcf)
# ploidy(gl.exi) <- 2
# # pop(gl.exi) <- pop.data.cleaned$V2
# 
# pop.new <- as.data.frame(as.factor(indNames(gl.exi)))
# colnames(pop.new)[1] <- "Ind"
# # pop.new$Site <- ""
# #
# 
# cleaned <- levels(pop.new$Ind)
# dirty <- levels(as.factor(pop.data$V1))
# d.in.c <- dirty %in% cleaned
# 
# pop.data.cleanaed <- pop.data %>% filter(V1 %in% cleaned)
# 
# 
# # pop.new$Site[which(str_detect(pop.new[,1], "Archi")==TRUE)] <- "Archidona"
# # pop.new$Site[which(str_detect(pop.new[,1], "VL")==TRUE)] <- "ViaLoreto"
# # pop.new$Site[which(str_detect(pop.new[,1], "")==TRUE)] <- "JatunSacha"
# 
# strata_df <- pop.data.cleanaed
# strata_df <- strata_df %>%
#   dplyr::rename(Nest = V2, Cluster = V3, Site = V4)
# 
# strata(gl.exi) <- strata_df
# head(strata(gl.exi))
# nest.site <- (strata(gl.exi, ~Nest/Site)) #show hierarchically
# (strata(gl.exi, ~Nest/Site))
# indNames(gl.exi)
# 
# 
# 
# ###### Subset populations ######
# library(poppr)
# 
# setPop(gl.exi) <- ~Site
# popNames(gl.exi)
# exi.js <- popsub(gl.exi, "JatunSacha")
# popNames(exi.js)
# setPop(exi.js) <- ~Nest
# indNames(exi.js)
# # gl.exi <- exi.js
# 
# exi.archi <- popsub(gl.exi, "Archidona")
# # archi <- which(str_detect(pop.data.cleaned$V1, "Archi") == TRUE)
# # archi.list <- as.list(as.character(pop.data.cleaned$V1[archi]))
# # exi.archi <- gl.exi[indNames(gl.exi) == archi.list]
# 
# setPop(exi.archi) <- ~Nest
# popNames(exi.archi)
# indNames(exi.archi)
# # remove NAs
# toRemove <- is.na(glMean(exi.archi, alleleAsUnit = FALSE))
# which(toRemove)
# exi.archi <- exi.archi[, !toRemove]
# # gl.exi <- exi.archi
# 
# exi.vl <- popsub(gl.exi, sublist=c("ViaLoreto"))
# indNames(exi.vl)
# 
# # remove NAs
# toRemove <- is.na(glMean(exi.vl, alleleAsUnit = FALSE))
# which(toRemove)
# exi.vl <- exi.vl[, !toRemove]
# 
# setPop(exi.vl) <- ~Nest
# popNames(exi.vl)
# #gl.exi <- exi.vl2

######### PCA ########
library(parallel)
detectCores()

setPop(gl.exi) <- ~Nest
exi.pca.nest <- glPca(gl.exi, nf = 3, parallel = TRUE, n.core = 4)
exi.pca.nest.scores <- as.data.frame(exi.pca.nest$scores)
exi.pca.nest.scores$pop <- unlist(strata(gl.exi, ~Nest))

ggplot(exi.pca.nest.scores, aes(x=PC1, y=PC2, colour=pop)) + 
  geom_point(size=2) +
  # stat_ellipse(level = 0.95, size = 1) +
  theme_bw()

setPop(gl.exi) <- ~Cluster
popNames(gl.exi)
exi.pca.cluster <- glPca(gl.exi, nf = 3, parallel = TRUE, n.core = 4)
exi.pca.cluster.scores <- as.data.frame(exi.pca.cluster$scores)
exi.pca.cluster.scores$pop <- unlist(strata(gl.exi, ~Cluster))


ggplot(exi.pca.cluster.scores, aes(x=PC1, y=PC2, colour=pop)) + 
  geom_point(size=2) +
  # stat_ellipse(level = 0.95, size = 1) +
  theme_bw()

setPop(gl.exi) <- ~Site
popNames(gl.exi)
exi.pca.site <- glPca(gl.exi, nf = 3, parallel = TRUE, n.core = 4)
exi.pca.site.scores <- as.data.frame(exi.pca.site$scores)
exi.pca.site.scores$pop <- unlist(strata(gl.exi, ~Site))

ggplot(exi.pca.site.scores, aes(x=PC1, y=PC2, colour=pop)) + 
  geom_point(size=2) +
  # stat_ellipse(level = 0.95, size = 1) +
  theme_bw()

library(ggplot2)
set.seed(9)

barplot(100*exi.pca.site$eig/sum(exi.pca.site$eig), col = heat.colors(50), main="PCA Eigenvalues")
barplot(100*exi.pca.cluster$eig/sum(exi.pca.cluster$eig), col = heat.colors(50), main="PCA Eigenvalues")
barplot(100*exi.pca.nest$eig/sum(exi.pca.nest$eig), col = heat.colors(50), main="PCA Eigenvalues")



ggsave(p1, filename = "figures/pca_col.jpg", dpi='retina', units='in',
       width = 8, height = 6)


p2 <- ggplot(exi.pca.scores, aes(x=PC1, y=PC2, colour=site)) + 
  geom_point(size=2) +
  stat_ellipse(level = 0.95, size = 1) +
  theme_bw()

ggsave(p2, filename = "figures/pca_site.jpg", dpi='retina', units='in',
       width = 8, height = 6)

####### DAPC #######

exi.dapc <- dapc(gl.exi, pop = ~Nest, n.pca = 3, n.da = 3, glPca = exi.pca.nest)
scatter(exi.dapc, cex = 2, clab=0, label.inds = list(air = 2, pch = NA))
compoplot(exi.dapc)


dapc.results <- as.data.frame(exi.dapc$posterior)
dapc.results$pop <- pop(gl.exi)
dapc.results$indNames <- rownames(dapc.results)


library(tidyr)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))
head(dapc.results, n = 6)
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
write.csv(dapc.results, "dapc.results.csv")

p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))+
  geom_bar(stat='identity') +
  facet_grid(~Original_Pop, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p

ggsave(p, filename = "../figures/exi_compoplot.jpg", dpi = "retina", height = 8, width = 15, units = "in")



## subset JS ##
library(poppr)
popNames(gl.exi)
js.exi <- popsub(gl.exi, sublist = 1)
setPop(js.exi) <- ~Cluster

exi.js.cluster <- glPca(js.exi, nf = 3, parallel = TRUE, n.core = 4)
exi.js.cluster.scores <- as.data.frame(exi.js.cluster$scores)
exi.js.cluster.scores$pop <- unlist(strata(js.exi, ~Cluster))
ggplot(exi.js.cluster.scores, aes(x=PC1, y=PC2, colour=pop)) + 
  geom_point(size=2) +
  stat_ellipse(level = 0.95, size = 1) +
  theme_bw()

exi.js.dapc <- dapc(js.exi, n.pca = 3, n.da = 2, glPca = exi.js.cluster)
scatter(exi.js.dapc, cex = 2, clab=0, label.inds = list(air = 2, pch = NA))
compoplot(exi.js.dapc)

dapc.results <- as.data.frame(exi.js.dapc$posterior)
dapc.results$pop <- pop(js.exi)
dapc.results$indNames <- rownames(dapc.results)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))
head(dapc.results, n = 6)
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
#write.csv(dapc.results, "dapc.results.csv")

p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))+
  geom_bar(stat='identity') +
  facet_grid(~Original_Pop, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p


setPop(js.exi) <- ~Nest
exi.js.nest <- glPca(js.exi, nf = 3, parallel = TRUE, n.core = 4)
exi.js.nest.scores <- as.data.frame(exi.js.nest$scores)
exi.js.nest.scores$pop <- unlist(strata(js.exi, ~Nest))
ggplot(exi.js.nest.scores, aes(x=PC1, y=PC2, colour=pop)) + 
  geom_point(size=2) +
  # stat_ellipse(level = 0.95, size = 1) +
  theme_bw()

exi.js.dapc <- dapc(js.exi, n.pca = 3, n.da = 2, glPca = exi.js.nest)
scatter(exi.js.dapc, cex = 2, clab=0, label.inds = list(air = 2, pch = NA))
compoplot(exi.js.dapc)

dapc.results <- as.data.frame(exi.js.dapc$posterior)
dapc.results$pop <- pop(js.exi)
dapc.results$indNames <- rownames(dapc.results)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))
head(dapc.results, n = 6)
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
#write.csv(dapc.results, "dapc.results.csv")

p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))+
  geom_bar(stat='identity') +
  facet_grid(~Original_Pop, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p

######## old #######
# 
# library('StAMPP')
# 
# stamppFst(js.exi, nclusters = 4)
# 
# ###### playing with vcfr to get some diagnostic plots ######
# library(vcfR)
# 
# exi.unfiltered <- read.vcfR("../filtering/populations.snps.vcf")
# exi.filtered <- read.vcfR("../filtering/raw.dp3mdp5m40.recode.vcf")
# 
# tidy.unfiltered <- vcfR2tidy(exi.unfiltered)
# tidy.unfiltered.gt <- as.data.frame(tidy.unfiltered$gt)
# hist(tidy.unfiltered.gt$gt_DP)
# table(tidy.unfiltered.gt$gt_DP)
# 
# 
# max(tidy.unfiltered.gt$gt_DP, na.rm = TRUE)
# min(tidy.unfiltered.gt$gt_DP, na.rm = TRUE)
# mean(tidy.unfiltered.gt$gt_DP, na.rm = TRUE)
# hist(tidy.unfiltered.gt$gt_GQ)
# 
# tidy.filtered <- vcfR2tidy(exi.filtered)
# tidy.filtered.gt <- as.data.frame(tidy.filtered$gt)
# hist(tidy.filtered.gt$gt_DP)
# max(tidy.filtered.gt$gt_DP, na.rm = TRUE)
# mean(tidy.filtered.gt$gt_DP, na.rm = TRUE)
# min(tidy.filtered.gt$gt_DP, na.rm = TRUE)
# 
# 
# max(tidy.filtered.gt$gt_GQ, na.rm = TRUE)
# min(tidy.filtered.gt$gt_GQ, na.rm = TRUE)
# hist(tidy.filtered.gt$gt_GQ)
# table(tidy.filtered.gt$gt_GQ)
# 
# # chrom.exi <- create.chromR(vcf)
# 
# dp <- extract.gt(exi.unfiltered, element = "DP", as.numeric = TRUE)
# bp1 <- boxplot(dp, col=2:8, las=3)
# 
# dp2 <- extract.gt(exi.filtered, element = "DP", as.numeric = TRUE)
# bp2 <- boxplot(dp2, col=2:8, las=3)
# 
# 
# gdepth <- read.csv("filtering/out.gdepth", header = TRUE, sep = "\t")
# 
# 



############## K- means clustering tutorial ##############
# use gl.exi

library(adegenet)
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
exi.elbow <- p1
ggsave(exi.elbow, filename = "figures/exi.elbow.jpeg", dpi = "retina",
       units = "in", height = 4.5, width = 6)

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

my_pal <- RColorBrewer::brewer.pal(n=6, name = "Dark2")

p2 <- ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group)) +
  geom_point(size = 4, shape = 21) +
  theme_bw() +
  scale_color_manual(values=c(my_pal)) +
  # stat_ellipse(level = 0.95, size = 1) +
  scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))+
  geom_text(aes(label=names), nudge_x=0.3, nudge_y=0.6)+
  stat_ellipse(level = 0.95, size = 1)
p2

exi.pca.plot <- p2
ggsave(exi.pca.plot, filename = "../../figures/exi.pca.JS.jpeg", dpi = "retina",
       units = "in", height = 6, width = 8)

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


library("ggpubr")
#tiff('dapc__k3_5_dapc.tiff', width=6.5, height=6.5, units='in', compression='lzw', res=300)
ggarrange(ggarrange(p1,
                    p2,
                    ncol = 2, labels = c("A", "B")),
          p3,
          nrow = 2,
          labels = c("", "C"),
          heights = c(1, 2)
)