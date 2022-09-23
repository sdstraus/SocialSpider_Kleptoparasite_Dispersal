library(dplyr)
library(stringr)
library(ggplot2)

barcodes <- read.table("sample_lists/barcodes.txt", sep="\t", header=F)

complete <- as.data.frame(as.character(barcodes[,3]))

write.table(complete, "sampleList_complete.txt", row.names = FALSE, col.names = FALSE, quote=FALSE)
# subset just eximius colonies
m <- which(str_detect(barcodes$V3, "Ex") == TRUE)
ex_master <- complete[m]

# subset N1
N1 <- which(str_detect(barcodes$V3, "N1") == TRUE)
N1_list <- as.data.frame(barcodes[N1,]$V3)

# subset F2
F2 <- which(str_detect(barcodes$V3, "F2") == TRUE)
F2_list <- as.data.frame(barcodes[F2,]$V3)

# subset eximius individuals
ex <-ex_master[-c(N1,F2)]

write.table(as.factor(ex), "eximius_prefix.txt", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(as.factor(F2_list), "F2_prefix.txt", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(as.factor(N1_list), "N1_prefix.txt", row.names = FALSE, col.names = FALSE, quote=FALSE)

## create population map
ex <- read.csv("sample_lists/samplelist.txt", header = FALSE)
ex[,1] <- as.character(ex[,1])
ex$site <- ""

## by colony
ex$colony <- ""
ex$colony[which(str_detect(ex[,1], "Ex11")==TRUE)] <- "Exi11_JS"
ex$colony[which(str_detect(ex[,1], "Ex12")==TRUE)] <- "Exi12_JS"
ex$colony[which(str_detect(ex[,1], "Ex14")==TRUE)] <- "Exi14_JS"
ex$colony[which(str_detect(ex[,1], "Ex01")==TRUE)] <- "Exi01_JS"
ex$colony[which(str_detect(ex[,1], "Ex02")==TRUE)] <- "Exi02_JS"
ex$colony[which(str_detect(ex[,1], "Ex05")==TRUE)] <- "Exi05_JS"
ex$colony[which(str_detect(ex[,1], "Ex06A")==TRUE)] <- "Exi06A_JS"
ex$colony[which(str_detect(ex[,1], "Ex06B")==TRUE)] <- "Exi06B_JS"
ex$colony[which(str_detect(ex[,1], "Ex07")==TRUE)] <- "Exi07_JS"
ex$colony[which(str_detect(ex[,1], "Ex08")==TRUE)] <- "Exi08_JS"
ex$colony[which(str_detect(ex[,1], "Ex10")==TRUE)] <- "Exi10_JS"
ex$colony[which(str_detect(ex[,1], "Ex17")==TRUE)] <- "Exi17_JS"
ex$colony[which(str_detect(ex[,1], "Ex19")==TRUE)] <- "Exi19_JS"

ex$colony[which(str_detect(ex[,1], "Exi_VL_29.0.1")==TRUE)] <- "Exi_VL_29.0.1"
ex$colony[which(str_detect(ex[,1], "Exi_VL_29.0.2")==TRUE)] <- "Exi_VL_29.0.2"
ex$colony[which(str_detect(ex[,1], "Exi_VL_47.3.1")==TRUE)] <- "Exi_VL_47.3.1"
ex$colony[which(str_detect(ex[,1], "Exi_VL_9.5.2")==TRUE)] <- "Exi_VL_9.5.2"
ex$colony[which(str_detect(ex[,1], "Exi_VL_9.5.1")==TRUE)] <- "Exi_VL_9.5.1"
ex$colony[which(str_detect(ex[,1], "Exi_VL_43.6.2")==TRUE)] <- "Exi_VL_43.6.2"
ex$colony[which(str_detect(ex[,1], "Exi_VL_24.0")==TRUE)] <- "Exi_VL_24.0"
ex$colony[which(str_detect(ex[,1], "Exi_VL_10.5")==TRUE)] <- "Exi_VL_10.5"
ex$colony[which(str_detect(ex[,1], "Exi_VL_41.2")==TRUE)] <- "Exi_VL_41.2"
ex$colony[which(str_detect(ex[,1], "Exi_VL_17.5")==TRUE)] <- "Exi_VL_17.5"
ex$colony[which(str_detect(ex[,1], "Exi_VL_Bridge.1")==TRUE)] <- "Exi_VL_Bridge.1"
ex$colony[which(str_detect(ex[,1], "Exi_VL_24.1")==TRUE)] <- "Exi_VL_24.1"
ex$colony[114] <- "Exi_VL_24.1"

ex$colony[which(str_detect(ex[,1], "Exi_Archi_11.8.1")==TRUE)] <- "Exi_Archi_11.8.1"
ex$colony[which(str_detect(ex[,1], "Exi_Archi_11.8.2")==TRUE)] <- "Exi_Archi_11.8.2"
ex$colony[which(str_detect(ex[,1], "Exi_Archi_11.8.3")==TRUE)] <- "Exi_Archi_11.8.3"
ex$colony[which(str_detect(ex[,1], "Exi_Archi_11.8.4")==TRUE)] <- "Exi_Archi_11.8.4"
ex$colony[which(str_detect(ex[,1], "Exi_Archi_11.8.5")==TRUE)] <- "Exi_Archi_11.8.5"

#ex$site <- NULL
#write.table(ex, "eximius_popmap_colony.txt", row.names = FALSE, col.names = FALSE, quote=FALSE)


# add sites
archidona <- which(str_detect(ex[,1], "Archi") == TRUE)
ex$site[archidona] <- "Archidona"

vl <- which(str_detect(ex[,1], "VL") == TRUE)
ex$site[vl] <- "ViaLoreto"

ex$site[c(-vl, -archidona)] <- "JatunSacha"
#write.table(ex, "eximius_popmap_site.txt", row.names = FALSE, col.names = FALSE, quote=FALSE)


write.table(ex, "pop_maps/eximius_popmap_colony_site.txt", row.names = FALSE, col.names = FALSE, 
            quote=FALSE, sep = "\t")


########## N1 #########
## repeat above with with N1
## by colony
N1_list$colony <- ""
N1_list$colony[which(str_detect(N1_list[,1], "Ex02")==TRUE)] <- "Exi02_JS"
N1_list$colony[which(str_detect(N1_list[,1], "Ex10")==TRUE)] <- "Exi10_JS"
N1_list$colony[which(str_detect(N1_list[,1], "Ex08")==TRUE)] <- "Exi08_JS"
N1_list$colony[which(str_detect(N1_list[,1], "Ex06A")==TRUE)] <- "Exi06A_JS"
N1_list$colony[which(str_detect(N1_list[,1], "Ex06B")==TRUE)] <- "Exi06B_JS"
N1_list$colony[which(str_detect(N1_list[,1], "Ex12")==TRUE)] <- "Exi12_JS"
N1_list$colony[which(str_detect(N1_list[,1], "Ex17")==TRUE)] <- "Exi17_JS"
N1_list$colony[which(str_detect(N1_list[,1], "Ex19")==TRUE)] <- "Exi19_JS"


N1_list$colony[which(str_detect(N1_list[,1], "Exi_Archi_11.8.1")==TRUE)] <- "Exi_Archi_11.8.1"
N1_list$colony[which(str_detect(N1_list[,1], "Exi_Archi_11.8.3")==TRUE)] <- "Exi_Archi_11.8.3"
N1_list$colony[which(str_detect(N1_list[,1], "Exi_Archi_11.8.4")==TRUE)] <- "Exi_Archi_11.8.4"
N1_list$colony[which(str_detect(N1_list[,1], "Exi_Archi_11.8.5")==TRUE)] <- "Exi_Archi_11.8.5"


N1_list$colony[which(str_detect(N1_list[,1], "Exi_VL_9.5.2")==TRUE)] <- "Exi_VL_9.5.2"
N1_list$colony[which(str_detect(N1_list[,1], "Exi_VL_9.5.1")==TRUE)] <- "Exi_VL_9.5.1"
N1_list$colony[which(str_detect(N1_list[,1], "Exi_VL_43.6")==TRUE)] <- "Exi_VL_43.6"
N1_list$colony[which(str_detect(N1_list[,1], "Exi_VL_24.1")==TRUE)] <- "Exi_VL_24.1"



# add sites
N1_list$site <- ""
archidona <- which(str_detect(N1_list[,1], "Archi") == TRUE)
N1_list$site[archidona] <- "Archidona"

vl <- which(str_detect(N1_list[,1], "VL") == TRUE)
N1_list$site[vl] <- "ViaLoreto"

N1_list$site[c(-vl, -archidona)] <- "JatunSacha"


write.table(N1_list, "pop_maps/N1_popmap_colony_site.txt", row.names = FALSE, col.names = FALSE, 
            quote=FALSE, sep = "\t")




########### F2 #############
## repeat above with with F2
## by colony


F2_list$colony <- ""
F2_list$colony[which(str_detect(F2_list[,1], "Ex02")==TRUE)] <- "Exi02_JS"
F2_list$colony[which(str_detect(F2_list[,1], "Ex17")==TRUE)] <- "Exi17_JS"
F2_list$colony[which(str_detect(F2_list[,1], "Ex16")==TRUE)] <- "Exi16_JS"
F2_list$colony[which(str_detect(F2_list[,1], "Ex08")==TRUE)] <- "Exi08_JS"
F2_list$colony[which(str_detect(F2_list[,1], "Ex12")==TRUE)] <- "Exi12_JS"
F2_list$colony[which(str_detect(F2_list[,1], "Ex14")==TRUE)] <- "Exi14_JS"
F2_list$colony[which(str_detect(F2_list[,1], "Ex07")==TRUE)] <- "Exi07_JS"


write.table(F2_list, "pop_maps/F2_popmap_colony_site.txt", row.names = FALSE, col.names = FALSE, 
            quote=FALSE, sep = "\t")



####### Cluster ########
# Eximius
coord <- read.csv("Rscripts/data/2017_col_locations.csv")

coord$latitude <- str_replace(coord$latitude, "\xa1", "")
coord$longitude <- str_replace(coord$longitude, "\xa1", "")

rows <- which(str_detect(coord$Nest, "ex") == TRUE)
exi.coord <- coord[rows,]
rows <- which(str_detect(coord$Nest, "Ex") == TRUE)
exi.coord <- rbind(exi.coord, coord[rows,])


exi.coord$latitude <- str_replace(exi.coord$latitude, "S0", "")
exi.coord$longitude <- str_replace(exi.coord$longitude, "W0", "")
exi.coord$longitude <- str_replace(exi.coord$longitude, "S0", "")

exi.coord$Nest <- str_replace(exi.coord$Nest, "2017_", "")
exi.coord$Nest <- str_replace(exi.coord$Nest, "ex_", "Ex")

exi.coord$latitude <- as.numeric(exi.coord$latitude) * -1
exi.coord$longitude <- as.numeric(exi.coord$longitude) * -1

write.csv(exi.coord, "../data/eximius_coordiates.csv")

map <- ggplot() + 
  geom_point(data = exi.coord, aes(x = longitude, y = latitude, colour = Nest, shape = locale))+
  geom_text()+
  cowplot::theme_cowplot()

pop.map <- read.table("pop_maps/eximius_popmap_colony_site.txt", sep = "\t")
pop.map$cluster <- ""

#js
# 06A and 06B a cluster
pop.map$cluster <- pop.map[,3]
pop.map$cluster[which(str_detect(pop.map[,3], "Exi06")==TRUE)] <- "Exi06_JS"
pop.map$cluster[which(str_detect(pop.map[,3], "29.0")==TRUE)] <- "Exi_VL_29.0"
pop.map$cluster[which(str_detect(pop.map[,3], "9.5")==TRUE)] <- "Exi_VL_9.5"
pop.map$cluster[which(str_detect(pop.map[,3], "11.8.2")==TRUE | 
                  str_detect(pop.map[,3], "11.8.3")==TRUE)] <- "Exi_Archi_cluster"


pop.map2 <- pop.map[, c(1, 2, 4, 3)]


table(pop.map2$V2, pop.map2$cluster)

write.table(pop.map2, "pop_maps/eximius_popmap_colony_site_cluster.txt", row.names = FALSE, col.names = FALSE, 
            quote=FALSE, sep = "\t")



#### N1
pop.map.n1 <- read.table("pop_maps/n1_popmap_colony_site.txt", sep = '\t', header = FALSE)
pop.map.n1$cluster <- ""

pop.map.n1$cluster <- pop.map.n1[,2]
pop.map.n1$cluster[which(str_detect(pop.map.n1[,2], "Exi06")==TRUE)] <- "Exi06_JS"
pop.map.n1$cluster[which(str_detect(pop.map.n1[,2], "29.0")==TRUE)] <- "Exi_VL_29.0"
pop.map.n1$cluster[which(str_detect(pop.map.n1[,2], "9.5")==TRUE)] <- "Exi_VL_9.5"
pop.map.n1$cluster[which(str_detect(pop.map.n1[,2], "11.8.2")==TRUE | 
                        str_detect(pop.map.n1[,2], "11.8.3")==TRUE)] <- "Exi_Archi_cluster"


pop.map2 <- pop.map.n1[, c(1, 2, 4, 3)]

write.table(pop.map2, "pop_maps/n1_popmap_colony_site_cluster.txt", row.names = FALSE, col.names = FALSE, 
            quote=FALSE, sep = "\t")




#### Map ####
library(ggmap)
library(geom_sf)

map <- get_map('ecuador')
