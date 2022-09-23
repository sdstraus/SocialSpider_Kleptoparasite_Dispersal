setwd("/Users/sam/PhDThesis/GBS")

library(poppr)
library(adegenet)
library(hierfstat)

######## pairwise FST ########
library(SNPRelate)
snpgdsVCF2GDS("filtering/exi/filtered.recode.vcf",
              "filtering/exi/filtered.gds",
              method="biallelic.only")

#Load the gds file
genofile <- snpgdsOpen("filtering/exi/filtered.gds")
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
test.pop <- pop.data.cleaned[,2]
test.pop <- pop.data[,2]

v <- snpgdsFst(genofile, sample.id = sample.id, population = as.factor(test.pop), method = "W&H02")
summary(v$FstSNP)
beta <- v$Beta

pop_code <- pop.data.cleaned[,2]
pop_code <- test.pop
poplevels <- levels(as.factor(test.pop))
# pairwise populations matrix creation
res= outer(X= poplevels , Y= poplevels, 
           FUN = function(h,k){
             paste(h,k,sep = "/")
           }
)
colnames(res)=poplevels
rownames(res)=poplevels

as.data.frame(res)

# pairwise population matrix FST calculation

for(i in poplevels) {
  for(j in poplevels) {
    popelem= unlist(strsplit(res[i,j],"/"))
    
    #takes selection of samples an population to use for each pair
    flag<- pop_code %in% c(popelem[1],popelem[2])
    samp.sel<- sample.id[flag]
    pop.sel<- pop_code[flag]
    
    if (popelem[1]==popelem[2]){
      result="0"}else{
        result = snpgdsFst(genofile, sample.id=samp.sel, population=as.factor(pop.sel), 
                           autosome.only=FALSE, method="W&C84")$MeanFst
      }
    res[i,j]=as.character(result)
  }
}


#https://github.com/rgiannico/RpairwiseFST/blob/master/RpairwiseFST.R
# final edits and prints

######## caculate FST #######

transformNegativeFSTtozero<- TRUE
keepLowerTriangularMatrixOnly<- TRUE


# cell transformation from characters to numerics

res=data.frame(apply(res, 2, function(x) as.numeric(as.character(x))))
rownames(res)=poplevels

# trasforming negative FST into zero
# Negative Fst are technical artifact of the computation (see Roesti el al. 2012) and will be automatically replaced with zero inside this function.
if (transformNegativeFSTtozero==TRUE){
  for(i in poplevels) {
    for(j in poplevels) {
      if (res[i,j]<0){res[i,j]<-0} 
    }
  }
}

# keep only the lower triangular matrix
# (set the upper triangular to zero)
if(keepLowerTriangularMatrixOnly == TRUE){
  res[upper.tri(res)]<-0
}


# PRINT pairwise FST matrix to text file

timewithspaces= Sys.time()
timeAttr= gsub(pattern = " ",replacement = "_",x =  timewithspaces)
outfile <- paste("pairFstMatrix_", timeAttr, ".txt", sep="")
write.table(x = res, file = outfile, sep = "\t", dec = ".", 
            quote = F, row.names = T,col.names = NA)



############## FSTPAIRWISE TABLE PLOT #########################
# taken from Arlequin's pairFstMatrix.r (Author: Heidi Lischer)
# with very few edits


numericMatrix=res

# preliminar functions
#Mirror matrix (left-right)
mirror.matrix <- function(x) {
  xx <- as.data.frame(x);
  xx <- rev(xx);
  xx <- as.matrix(xx);
  xx;
}

#Rotate matrix 270 clockworks
rotate270.matrix <- function(x) {
  mirror.matrix(t(x))
}

Matrix <- rotate270.matrix(numericMatrix)



ColorRamp <- colorRampPalette(c("white", "steelblue1", "blue3"))

timewithspaces= Sys.time()
timeAttr= gsub(pattern = " ",replacement = "_",x =  timewithspaces)
outfileGraphic <- paste("pairFstMatrix_", timeAttr, ".png", sep="")
# outfileGraphic <- paste(outfile, "pairFstMatrix ", timeAttr, ".pdf", sep="")

#save graphic
png(outfileGraphic, width=1300, height=1300, res=144)
# pdf(outfileGraphic, width = 10, height = 10)  

smallplot <- c(0.874, 0.9, 0.18, 0.83)
bigplot <- c(0.13, 0.85, 0.14, 0.87)

old.par <- par(no.readonly = TRUE)

# draw legend
par(plt = smallplot)

# get legend values
Min <- min(Matrix, na.rm=TRUE)
Max <- max(Matrix, na.rm=TRUE)
binwidth <- (Max - Min) / 64
y <- seq(Min + binwidth/2, Max - binwidth/2, by = binwidth)
z <- matrix(y, nrow = 1, ncol = length(y))

image(1, y, z, col = ColorRamp(64),xlab="", ylab="", axes=FALSE)

# adjust axis if only one value exists
if(Min == Max){
  axis(side=4, las = 2, cex.axis=0.8, at=Min, labels=round(Min, 2))
} else {
  axis(side=4, las = 2, cex.axis=0.8)
}

box()
mtext(text=expression(bold(F[ST])), side=4, line=2.5, cex=1.1)


#draw main graphic3
a <- ncol(numericMatrix)
b <- nrow(numericMatrix)

x <- c(1:a)
y <- c(1:b)

par(new = TRUE, plt = bigplot)

image(x,y,as.matrix(Matrix), col=ColorRamp(64),
      main=expression(bold(Matrix~of~pairwise~F[ST])), xlab="",
      ylab="", axes=FALSE)
box()

#add labels
Labels=poplevels
if(is.null(Labels)){
  axis(1, at = c(1:a))
  axis(2, at = c(1:b), labels=c(b:1))
  mtext(side = 1, at =(a/2), line = 2.5, text = "Population", cex=1,
        font=2)
  mtext(side = 2, at =(b/2), line = 2.7, text = "Population", cex=1,
        font=2)
} else{
  axis(1, at = c(1:a), labels=Labels[1:length(Labels)], cex.axis=0.75,
       las=2)
  axis(2, at = c(1:b), labels=Labels[length(Labels):1], cex.axis=0.75,
       las=2)
}


par(old.par)  #reset graphic parameters

dev.off()


## histogram of Fsts
res.list <- unlist(split(res, seq(nrow(res))))
res.list.pos <- res.list[which(res.list>0)]
hist(res.list.pos, breaks = 20, xlab = "Colony pairwise Fst", main = "")
mean(res.list.pos)
sd(res.list.pos) * 1.96
range(res.list.pos)


######### using hierfstat ########
library(dartR)
library(hierfstat)
gi.exi <- gl2gi(gl.exi)
js.exi.gi <- gl2gi(exi.js)
archi.exi.gi <- gl2gi(exi.archi)
exi.vl <- gl2gi(exi.vl)

# 
# x.mat <- as.matrix(exi.archi) # x is a genlight object
# x.mat[x.mat == 0] <- "1/1" # homozygote reference
# x.mat[x.mat == 1] <- "1/2" # heterozygote
# x.mat[x.mat == 2] <- "2/2" # homozygote alternate
# x.gid <- df2genind(exi.archi, sep = "/", ploidy = 2)

gi.exi.2 <- vcfR2genind(vcf)
strata(gi.exi) <- strata_df
head(strata(gi.exi.2))
setPop(gi.exi.2) <- ~Site
exi.archi <- poppr::popsub(gi.exi, "Archidona")
setPop(exi.archi) <- ~Nest
popNames(exi.archi)

saveRDS(gi.exi, file="../rds/gi.exi.RDS")

gi.exi<- readRDS("../rds/gi.exi.RDS")

# Pop as Nest
setPop(gi.exi) <- ~Nest
wc(gi.exi)


# Pop as Nest/Cluster
setPop(js.exi.gi) <- ~Nest/Cluster
wc(js.exi.gi)

setPop(archi.exi.gi) <- ~Nest/Cluster
wc(archi.exi.gi)

setPop(exi.vl) <- ~Nest/Cluster
wc(exi.vl)

setPop(gi.exi) <- ~Nest/Cluster
wc(gi.exi)

setPop(js.exi.gi) <- ~Nest
wc(js.exi.gi)

setPop(archi.exi.gi) <- ~Nest
wc(archi.exi.gi)

setPop(exi.vl) <- ~Nest
wc(exi.vl)




het <- gl.Ho(exi.js)
mean(gl.Ho(gl.exi))
mean(gl.Ho(exi.js))
mean(gl.Ho(exi.vl))
mean(gl.Ho(exi.archi))

het.n1 <- gl.Ho(gl.n1)
mean(het.n1)

het.f2 <- gl.Ho(gl.f2)
mean(het.f2)

bs.nc <- basic.stats(js.exi.gi)
wc.exi <- wc(js.exi.gi)
save(wc.exi, file="Rscripts/Fst.WC.exi.Rdata")


perloc <- as.data.frame(bs.nc$perloc)
perloc.fst <- perloc$Fst
perloc.fst[which(perloc.fst < 0)] <- 0
hist(perloc.fst)
mean(perloc.fst, na.rm = T)


pop.tabl <- poppr::poppr(js.exi.gi)
write.csv(pop.tabl, file="Rscripts/js_poptable.csv")

loc.tabl <- poppr::locus_table(js.exi.gi)

# pairwise
matFst <- genet.dist(gi.exi, method = "Nei87")
matFst.site <- matFst
save(matFst.site, file="Rscripts/matFst.site.Rdata")

# failed
# matFst.WC <- genet.dist(gi.exi, method = "WC84")
# save(matFst.site.WC, file="Rscripts/matFst.site.Rdata")

####### N1 ######

gi.n1 <- gl2gi(gl.n1)
saveRDS(gi.n1, file="Rscripts/gi.n1.RDS")

gi.n1 <- readRDS("Rscripts/gi.n1.RDS")
bs.n1 <- basic.stats(gi.n1)

wc.n1 <- wc(gi.n1)
save(wc.n1, file="Rscripts/Fst.WC.N1.Rdata")

####### F2 ######

gi.f2 <- gl2gi(gl.f2)
saveRDS(gi.f2, file="Rscripts/gi.f2.RDS")

gi.f2 <- readRDS("Rscripts/gi.f2.RDS")
bs.f2 <- basic.stats(gi.f2)


wc.f2 <- wc(gi.f2)
save(wc.f2, file="Rscripts/Fst.WC.F2.Rdata")


