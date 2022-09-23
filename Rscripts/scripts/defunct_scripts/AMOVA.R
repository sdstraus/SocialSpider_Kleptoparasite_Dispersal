setwd("/Users/sam/PhDThesis/GBS")

library(tidyr)
library(adegenet)
library(poppr)
library(dartR)

# 
# # for eximius
# exi.vcf <- read.vcfR("filtering/exi/filtered.recode.vcf")
# exi.pop.data <- read.table("pop_maps/eximius_popmap_colony_site_cluster.txt", sep = '\t', header = FALSE)
# all(colnames(exi.vcf@gt)[-1] == exi.pop.data$AccessID)
# 
# # filtered list has fewer individuals
# exi.tidy.filtered <- vcfR2tidy(exi.vcf)
# exi.tidy.filtered.gt <- as.data.frame(exi.tidy.filtered$gt)
# 
# # max(exi.tidy.filtered.gt$gt_DP, na.rm = TRUE)
# # min(exi.tidy.filtered.gt$gt_DP, na.rm = TRUE)
# # mean(exi.tidy.filtered.gt$gt_DP, na.rm = TRUE)
# 
# new.list <- as.data.frame(as.factor(exi.tidy.filtered.gt$Indiv))
# colnames(new.list)[1] <- "Ind"
# 
# cleaned <- levels(new.list[,1])
# dirty <- levels(as.factor(exi.pop.data$V1))
# d.in.c <- dirty %in% cleaned
# rows <- which(d.in.c == TRUE)
# exi.pop.data.cleaned <- exi.pop.data[rows,]
# 
# #convert to genlight object
# gl.exi <- vcfR2genlight(exi.vcf)
# ploidy(gl.exi) <- 2
# pop(gl.exi) <- exi.pop.data.cleaned$V2
# 
# strata_df <- exi.pop.data.cleaned[,-1]
# strata_df <- strata_df %>% 
#   rename(Nest = V2, Cluster = V3, Site = V4)
# 
# strata(gl.exi) <- strata_df
# gl.exi@other$loc.metrics.flags
# 


gl.dist.pop(gl.exi, method = "euclidean")


euc.exi <- bitwise.dist(gl.exi, euclidean = TRUE)
exi.msn <- poppr.msn(gl.exi, euc.exi, showplot = FALSE, include.ties= TRUE)
node.size <- rep(2, times = nInd(gl.exi))
names(node.size) <- indNames(gl.exi)
igraph::vertex.attributes(exi.msn$graph)$size <- node.size
set.seed(9)
plot_poppr_msn(gl.exi, exi.msn, palette = RColorBrewer::brewer.pal(n = nPop(gl.exi), name = "Dark2"), gadj = 70)

gl <- gl.filter.callrate(gl.exi)
gl.report.monomorphs(gl.exi)
gl.compliance.check(gl.exi)

# need to remove loci with missing loci to do amova
# have to use genind object
gi.exi <- gl2gi(gl.exi)
strata(gi.exi) <- strata_df
missingno(gi.exi)

exi.amova <- poppr.amova(gi.exi, ~Site/Cluster/Nest)
amova.test <- randtest(exi.amova)
plot(amova.test)

exi.amova2 <- poppr.amova(gi.exi, ~Site/Nest)
amova.test2 <- randtest(exi.amova2)
plot(amova.test2)
