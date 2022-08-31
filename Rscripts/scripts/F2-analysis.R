# sam's laptop
setwd("/Users/sam/PhDThesis/GBS")

library(vcfR)
library(tidyr)
library(adegenet)
library(SNPRelate)
library(ggplot2)
library(parallel)

######## create genlight ########
# read vcf
vcf <- read.vcfR("../../filtering/f2/minDP3.maxmiss.40.70.snps.recode.vcf")
pop.data <- read.table("../../pop_maps/f2_popmap_colony_site.txt", sep = '\t', header = FALSE)
all(colnames(vcf@gt)[-1] == pop.data$AccessID)

# prune down pop.data to match inds in the vcf file
tidy.vcf <- vcfR2tidy(vcf)
tidy.vcf.gt <- as.data.frame(tidy.vcf$gt) # extract genotype data

new.list <- as.data.frame(as.factor(tidy.vcf.gt$Indiv)) # list of individuals
colnames(new.list)[1] <- "Ind"

cleaned <- levels(new.list[,1])
dirty <- levels(as.factor(pop.data$V1))
d.in.c <- dirty %in% cleaned
rows <- which(d.in.c == TRUE)
pop.data.cleaned <- pop.data[rows,]


#convert to genlight object
gl.f2 <- vcfR2genlight(vcf)
ploidy(gl.f2) <- 2
pop(gl.f2) <- pop.data.cleaned$V2

##### quick PCA ########

detectCores()
f2.pca <- glPca(gl.f2, nf = 3, parallel = TRUE, n.core = 4)
f2.pca.scores <- as.data.frame(f2.pca$scores)
f2.pca.scores$pop <- pop(gl.f2)

set.seed(9)
barplot(100*f2.pca$eig/sum(f2.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")

ggplot(f2.pca.scores, aes(x=PC1, y=PC2, colour=pop)) + 
  geom_point(size=2) +
  stat_ellipse(level = 0.95, size = 1) +
  theme_bw()

############## K- means clustering ##############
gl.f2 <- readRDS("../rds/gl.f2.RDS")
pop.data.cleaned.f2 <- read.csv("../data/F2.pop.data.cleaned.csv")

f2.pca <- glPca(gl.f2, nf = 3, parallel = TRUE, n.core = 4)
f2.pca.scores <- as.data.frame(f2.pca$scores)
f2.pca.scores$pop <- pop(gl.f2)



maxK <- 10
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(gl.f2, n.pca = 40, choose.n.clust = FALSE,  
                       max.n.clust = maxK, glPca = f2.pca)
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

f2.elbow <- p1
ggsave(f2.elbow, filename = "figures/f2.elbow.jpeg", dpi = "retina",
       units = "in", height = 4.5, width = 6)

# DAPC
my_k <- c(3)

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] <- find.clusters(gl.f2, n.pca = 40, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(gl.f2, pop = grp_l[[i]]$grp, n.pca = 40, n.da = my_k[i], glPca = f2.pca)
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}


my_df <- as.data.frame(dapc_l[[ 1 ]]$ind.coord)
my_df$Group <- dapc_l[[ 1 ]]$grp
head(my_df)

names <- rownames(my_df)
my_df <- cbind(names, my_df)
my_df <- full_join(my_df, pop.data.cleaned.f2, by = c("names" = "V1"))

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

f2.Plot <- p2
ggsave(f2.Plot, filename="figures/f2.pca.jpeg", dpi = "retina", units = "in", 
       height = 4.5, width = 6)




##### compoplot
f2.dapc <- dapc(gl.f2, n.pca = 3, n.da = 2, glPca = f2.pca)

scatter(f2.dapc)
compoplot(f2.dapc)

dapc.results <- as.data.frame(f2.dapc$posterior)
dapc.results$pop <- pop(gl.f2)
dapc.results$indNames <- rownames(dapc.results)


dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))
head(dapc.results, n = 6)
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
write.csv(dapc.results, "dapc.results.csv")

p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop)) +
  geom_bar(stat='identity') + 
  facet_grid(~Original_Pop, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p
f2.compoplot <- p
ggsave(f2.compoplot, filename = "figures/f2.compoplot.jpeg", dpi = "retina",
       height = 6, width = 8, units = 'in')




tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]
tmp$Isolate <- rownames(tmp)
tmp <- melt(tmp, id = c("Isolate", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$Region <- pop.data.cleaned.f2$V2
my_df <- tmp


# pop.data.cleaned <- read.csv("../data/N1.pop.data.cleaned.csv")

for(i in 1:length(dapc_l)){
  tmp <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K <- my_k[i]
  tmp$Isolate <- rownames(tmp)
  tmp <- melt(tmp, id = c("Isolate", "K"))
  names(tmp)[3:4] <- c("Group", "Posterior")
  tmp$Region <- pop.data.cleaned.f2$V2
  
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
setPop(gl.f2) <- ~Site
euc.f2 <- bitwise.dist(gl.f2, euclidean = TRUE)
f2.msn <- poppr.msn(gl.f2, euc.f2, showplot = FALSE, include.ties= TRUE)
node.size <- rep(2, times = nInd(gl.f2))
names(node.size) <- indNames(gl.f2)
igraph::vertex.attributes(f2.msn$graph)$size <- node.size
set.seed(9)
f2.msn.plot <- plot_poppr_msn(gl.f2, f2.msn, palette = RColorBrewer::brewer.pal(n = nPop(gl.f2), name = "Dark2"), gadj = 70)

