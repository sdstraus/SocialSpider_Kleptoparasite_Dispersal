#### eximius ####
setwd("~/UBC/PhDThesis/GBS/filtering/exi")

## eximius populations run:
# populations -P ./denovo_map/stacks-exi/ --popmap ./pop_maps/eximius_popmap_colony_site.txt -t 14 --vcf -r 70 -p 3 --min-maf 0.05

###### exi filter 1 ########
# minimum percentage of individuals in a population required to process a locus: 70%
# max missing 0.5, minDP 3

## eximius - original submission filtering steps ##
system("vcftools --vcf eximius.snps.vcf --max-missing 0.50 --mac 3 --minDP 3 --min-meanDP 5 --recode --recode-INFO-all --out raw.mm50mac3dp3mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm50mac3dp3mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm50mac3dp3mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out exi.1")

# exi.1.recode.vcf
# 99 out of 138 Individuals
# After filtering, kept 11087 out of a possible 11087 Sites

######## exi filter 2 ########

# minimum percentage of individuals in a population required to process a locus: 70%
# max missing 0.5, minDP 5

## eximius - original submission filtering steps ##
system("vcftools --vcf eximius.snps.vcf  --max-missing 0.50 --mac 3 --recode --recode-INFO-all --out raw.mm50mac3")

# minDP 
system("vcftools --vcf raw.mm50mac3.recode.vcf --minDP 5 --recode --recode-INFO-all --out raw.mm50mac3dp5")

# minimum mean depth
system("vcftools --vcf raw.mm50mac3dp5.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm50mac3dp5mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm50mac3dp5mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm50mac3dp5mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out exi.2")

# exi.2.recode.vcf
# 98 out of 138 Individuals
# After filtering, kept 11087 out of a possible 11087 Sites

###### exi filter 3 ########
# minimum percentage of individuals in a population required to process a locus: 70%
# max missing 0.6, minDP 3

## eximius - original submission filtering steps ##
system("vcftools --vcf eximius.snps.vcf  --max-missing 0.60 --mac 3  --recode --recode-INFO-all --out raw.mm60mac3")

# minDP 
system("vcftools --vcf raw.mm60mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.mm60mac3dp3")

# minimum mean depth
system("vcftools --vcf raw.mm60mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm60mac3dp3mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm60mac3dp3mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm60mac3dp3mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out exi.3")

# exi.3.recode.vcf
# 100 out of 138 Individuals
# After filtering, kept 7595 out of a possible 7595 Sites

###### exi filter 4 ########
# minimum percentage of individuals in a population required to process a locus: 70%
# max missing 0.6, minDP 5

## eximius - original submission filtering steps ##
system("vcftools --vcf eximius.snps.vcf  --max-missing 0.60 --mac 3  --recode --recode-INFO-all --out raw.mm60mac3")

# minDP 
system("vcftools --vcf raw.mm60mac3.recode.vcf --minDP 5 --recode --recode-INFO-all --out raw.mm60mac3dp5")

# minimum mean depth
system("vcftools --vcf raw.mm60mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm60mac3dp5mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm60mac3dp5mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm60mac3dp5mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out exi.4")

# exi.4.recode.vcf
# 100 out of 138 Individuals
# After filtering, kept 7595 out of a possible 7595 Sites


###### exi filter 5 ########
# minimum percentage of individuals in a population required to process a locus: 70%
# max missing 0.4, minDP 3

system("vcftools --vcf eximius.snps.vcf  --max-missing 0.40 --mac 3  --recode --recode-INFO-all --out raw.mm40mac3")

# minDP 
system("vcftools --vcf raw.mm40mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.mm40mac3dp3")

# minimum mean depth
system("vcftools --vcf raw.mm40mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm40mac3dp3mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm40mac3dp3mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm40mac3dp3mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out exi.5")

# exi.5.recode.vcf
# 99 out of 138 Individuals
# After filtering, kept 14213 out of a possible 14213 Sites


###### exi filter 6 ########
# minimum percentage of individuals in a population required to process a locus: 70%
# max missing 0.4, minDP 5

system("vcftools --vcf eximius.snps.vcf  --max-missing 0.40 --mac 3  --recode --recode-INFO-all --out raw.mm40mac3")

# minDP 
system("vcftools --vcf raw.mm40mac3.recode.vcf --minDP 5 --recode --recode-INFO-all --out raw.mm40mac3dp5")

# minimum mean depth
system("vcftools --vcf raw.mm40mac3dp5.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm40mac3dp5mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm40mac3dp5mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm40mac3dp5mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out exi.6")

# exi.6.recode.vcf
# 96 out of 138 Individuals
# After filtering, kept 14213 out of a possible 14213 Sites


############## exi filter 7 ###########
## eximius populations run:
# populations -P ./denovo_map/stacks-exi/ --popmap ./pop_maps/eximius_popmap_colony_site.txt -t 14 --vcf -r 80 -p 3 --min-maf 0.05
#
#minimum percentage of individuals in a population required to process a locus: 80%
# max missing 0.5, minDP 3

## eximius - original submission filtering steps ##
system("vcftools --vcf populations.80.snps.vcf  --max-missing 0.50 --mac 3  --recode --recode-INFO-all --out raw.mm50mac3")

# minDP 
system("vcftools --vcf raw.mm50mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.mm50mac3dp3")

# minimum mean depth
system("vcftools --vcf raw.mm50mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm50mac3dp3mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm50mac3dp3mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm50mac3dp3mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out exi.7")

# exi.7.recode.vcf
# 94 out of 138 Individuals
# After filtering, kept 8866 out of a possible 8866 Sites

######## exi filter 8 ########

# minimum percentage of individuals in a population required to process a locus: 80%
# max missing 0.5, minDP 5

## eximius - original submission filtering steps ##
system("vcftools --vcf populations.80.snps.vcf  --max-missing 0.50 --mac 3  --recode --recode-INFO-all --out raw.mm50mac3")

# minDP 
system("vcftools --vcf raw.mm50mac3.recode.vcf --minDP 5 --recode --recode-INFO-all --out raw.mm50mac3dp5")

# minimum mean depth
system("vcftools --vcf raw.mm50mac3dp5.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm50mac3dp5mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm50mac3dp5mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm50mac3dp5mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out exi.8")

# exi.8.recode.vcf
# 93 out of 138 Individuals
# After filtering, kept 8866 out of a possible 8866 Sites

###### exi filter 9 ########
# minimum percentage of individuals in a population required to process a locus: 80%
# max missing 0.6, minDP 3

## eximius - original submission filtering steps ##
system("vcftools --vcf populations.80.snps.vcf  --max-missing 0.60 --mac 3  --recode --recode-INFO-all --out raw.mm60mac3")

# minDP 
system("vcftools --vcf raw.mm60mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.mm60mac3dp3")

# minimum mean depth
system("vcftools --vcf raw.mm60mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm60mac3dp3mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm60mac3dp3mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm60mac3dp3mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out exi.9")

# exi.9.recode.vcf
# 100 out of 138 Individuals
# After filtering, kept 7595 out of a possible 7595 Sites


###### exi filter 10 ########
# minimum percentage of individuals in a population required to process a locus: 80%
# max missing 0.6, minDP 5

## eximius - original submission filtering steps ##
system("vcftools --vcf populations.80.snps.vcf  --max-missing 0.60 --mac 3 --minDP 5 --min-meanDP 5 --recode --recode-INFO-all --out raw.temp")

# # minDP 
# system("vcftools --vcf raw.mm60mac3.recode.vcf --minDP 5 --recode --recode-INFO-all --out raw.mm60mac3dp5")

# minimum mean depth
# system("vcftools --vcf raw.mm60mac3dp5.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm60mac3dp5mdp5")

# remove missing individuals
system("vcftools --vcf raw.temp.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.temp.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out exi.10")

# exi.10.recode.vcf
# 100 out of 138 Individuals
# After filtering, kept 7595 out of a possible 7595 Sites

###### exi filter 11 ########
# minimum percentage of individuals in a population required to process a locus: 80%
# max missing 0.4, minDP 3

## eximius - original submission filtering steps ##
system("vcftools --vcf populations.80.snps.vcf  --max-missing 0.40 --mac 3 --minDP 3 --min-meanDP 5 --recode --recode-INFO-all --out raw.temp")

# minDP 
# system("vcftools --vcf raw.mm40mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.mm40mac3dp3")
# 
# # minimum mean depth
# system("vcftools --vcf raw.mm40mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm40mac3dp3mdp5")

# remove missing individuals
system("vcftools --vcf raw.temp.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.temp.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out exi.11")

# exi.11.recode.vcf
# 94 out of 138 Individuals
# After filtering, kept 11634 out of a possible 11634 Sites


###### exi filter 12 ########
# minimum percentage of individuals in a population required to process a locus: 80%
# max missing 0.4, minDP 5

## eximius - original submission filtering steps ##
system("vcftools --vcf populations.80.snps.vcf  --max-missing 0.40 --mac 3  --recode --recode-INFO-all --out raw.mm40mac3")

# minDP 
system("vcftools --vcf raw.mm40mac3.recode.vcf --minDP 5 --recode --recode-INFO-all --out raw.mm40mac3dp5")

# minimum mean depth
system("vcftools --vcf raw.mm40mac3dp5.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm40mac3dp5mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm40mac3dp5mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm40mac3dp5mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out exi.12")

# exi.12.recode.vcf
# 91 out of 138 Individuals
# After filtering, kept 11634 out of a possible 11634 Sites

###### exi filter 13 ######
## eximius populations run:
# populations -P ./denovo_map/stacks-exi/ --popmap ./pop_maps/eximius_popmap_colony_site.txt -t 14 --vcf -r 60 -p 3 --min-maf 0.05

#minimum percentage of individuals in a population required to process a locus: 60%
# max missing 0.5, minDP 3

## eximius - original submission filtering steps ##
system("vcftools --vcf populations.60.snps.vcf  --max-missing 0.50 --mac 3  --recode --recode-INFO-all --out raw.mm50mac3")

# minDP 
system("vcftools --vcf raw.mm50mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.mm50mac3dp3")

# minimum mean depth
system("vcftools --vcf raw.mm50mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm50mac3dp3mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm50mac3dp3mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm50mac3dp3mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out exi.13")

# exi.13.recode.vcf
# 102 out of 138 Individuals
# After filtering, kept 13971 out of a possible 13971 Sites

######## exi filter 14 ########

# minimum percentage of individuals in a population required to process a locus: 60%
# max missing 0.5, minDP 5

## eximius - original submission filtering steps ##
system("vcftools --vcf populations.60.snps.vcf  --max-missing 0.50 --mac 3  --recode --recode-INFO-all --out raw.mm50mac3")

# minDP 
system("vcftools --vcf raw.mm50mac3.recode.vcf --minDP 5 --recode --recode-INFO-all --out raw.mm50mac3dp5")

# minimum mean depth
system("vcftools --vcf raw.mm50mac3dp5.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm50mac3dp5mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm50mac3dp5mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm50mac3dp5mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out exi.14")

# exi.14.recode.vcf
# 98 out of 138 Individuals
# After filtering, kept 13971 out of a possible 13971 Sites

###### exi filter 15 ########
# minimum percentage of individuals in a population required to process a locus: 60%
# max missing 0.6, minDP 3

## eximius - original submission filtering steps ##
system("vcftools --vcf populations.60.snps.vcf  --max-missing 0.60 --mac 3  --recode --recode-INFO-all --out raw.mm60mac3")

# minDP 
system("vcftools --vcf raw.mm60mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.mm60mac3dp3")

# minimum mean depth
system("vcftools --vcf raw.mm60mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm60mac3dp3mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm60mac3dp3mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm60mac3dp3mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out exi.15")

# exi.15.recode.vcf
# 106 out of 138 Individuals
# After filtering, kept 9814 out of a possible 9814 Sites


###### exi filter 16 ########
# minimum percentage of individuals in a population required to process a locus: 60%
# max missing 0.6, minDP 5

## eximius - original submission filtering steps ##
system("vcftools --vcf populations.60.snps.vcf  --max-missing 0.60 --mac 3  --recode --recode-INFO-all --out raw.mm60mac3")

# minDP 
system("vcftools --vcf raw.mm60mac3.recode.vcf --minDP 5 --recode --recode-INFO-all --out raw.mm60mac3dp5")

# minimum mean depth
system("vcftools --vcf raw.mm60mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm60mac3dp5mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm60mac3dp5mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm60mac3dp5mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out exi.16")

# exi.16.recode.vcf
# 106 out of 138 Individuals
# After filtering, kept 9814 out of a possible 9814 Sites

###### exi filter 17 ########
# minimum percentage of individuals in a population required to process a locus: 60%
# max missing 0.4, minDP 3

## eximius - original submission filtering steps ##
system("vcftools --vcf populations.60.snps.vcf  --max-missing 0.40 --mac 3  --recode --recode-INFO-all --out raw.mm40mac3")

# minDP 
system("vcftools --vcf raw.mm40mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.mm40mac3dp3")

# minimum mean depth
system("vcftools --vcf raw.mm40mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm40mac3dp3mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm40mac3dp3mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm40mac3dp3mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out exi.17")

# exi.17.recode.vcf
# 99 out of 138 Individuals
# After filtering, kept 17336 out of a possible 17336 Sites


###### exi filter 18 ########
# minimum percentage of individuals in a population required to process a locus: 60%
# max missing 0.4, minDP 5

## eximius - original submission filtering steps ##
system("vcftools --vcf populations.60.snps.vcf  --max-missing 0.40 --mac 3  --recode --recode-INFO-all --out raw.mm40mac3")

# minDP 
system("vcftools --vcf raw.mm40mac3.recode.vcf --minDP 5 --recode --recode-INFO-all --out raw.mm40mac3dp5")

# minimum mean depth
system("vcftools --vcf raw.mm40mac3dp5.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm40mac3dp5mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm40mac3dp5mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm40mac3dp5mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out exi.18")

# exi.18.recode.vcf
# 95 out of 138 Individuals
# After filtering, kept 17336 out of a possible 17336 Sites




##########
setwd("~/UBC/PhDThesis/GBS/admixture")

# use plink to convert to bed format
system("sed '16,$s/^/scaffold_/' exi.1.recode.vcf > exi.scaffolds.vcf")
# awk '{print "chr"$0}' file.in > file.out

# # replace chrom column with zeroes
# system("sed '16,$s/[0-9.]*/0/' filtered.recode.vcf > filtered.recode.scaffolds.zeroed.vcf")

system("./plink2 --vcf exi.scaffolds.vcf --allow-extra-chr 0 --double-id --make-bed --out exi.scaffolds")

## need .map file for admixture bootstrapping
system("./plink2 --bfile exi.scaffolds --recode ped --out exi.scaffolds")

# admixture test
system("./admixture --cv exi.scaffolds.bed 2 > log2.out")
# 
# cat > runadmixture.sh
# 
# for K in 1 2 3 4 5 6 7 8 9 10; \
# do ./admixture \ #location of admixture program.
# --cv=10 \
# exi.scaffolds.bed \
# $K | tee log${K}.out; \
# done
# 
# (ctrlD)

# #to get the CV (Cross validation) errors
system("grep -h CV log*out")

# large drop in error btw K2 and K3, with K=14 being the lowest total 

# CV error (K=1): 0.70637
# CV error (K=10): 0.24010
# CV error (K=11): 0.24220
# CV error (K=12): 0.21849
# CV error (K=13): 0.23627
# CV error (K=14): 0.21753
# CV error (K=15): 0.21257
# CV error (K=16): 0.22658
# CV error (K=17): 0.25221
# CV error (K=18): 0.27262
# CV error (K=19): 0.24057
# CV error (K=2): 0.58312
# CV error (K=20): 0.23612
# CV error (K=3): 0.39135
# CV error (K=4): 0.32527
# CV error (K=5): 0.31247
# CV error (K=6): 0.29089
# CV error (K=7): 0.28422
# CV error (K=8): 0.27574
# CV error (K=9): 0.24693

# run K = for 1000 bootstraps
# cat > runadmixture_boot.sh
# 
# for K in 3; \
# do ./admixture \
# -B1000 \
# --cv=10 \
# exi.scaffolds.bed \
# $K | tee log${K}.out; \
# done
# 
# (ctrlD)

#chmod 755 runadmixture_boot.sh

system("./runadmixture_boot.sh") #### ends after Bootstrap replicate #0? looks like i might need a plink map file?

library(RColorBrewer)
tbl=read.table("../../admixture/exi/exi.scaffolds.6.Q")
barplot(t(as.matrix(tbl)), col=brewer.pal(6, "Dark2"), xlab="Individual", ylab="Ancestry", border=NA)

popmap <-  read.table("../../pop_maps/eximius_popmap_colony_site_cluster.txt", sep = '\t', header = FALSE)
popmap <- popmap %>% rename("ind" = "V1", "site" = "V2", "nest" = "V3", "cluster" = "V4")
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

long.dat <- exi.df %>% pivot_longer(cols = V1:V6, names_to = "Kclust", values_to = "fraction", names_prefix = "V")

p <-  ggplot(long.dat, aes(x=ind, y=fraction, fill=Kclust)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(. ~ site, drop=TRUE, space="free", scales="free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(fill = "K cluster", y = "Ancestry", x = "Individual")
p
