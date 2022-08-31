
####### F1 - filter 1, -r 70, -max-missing 0.5 minDP 3 #########
setwd("~/UBC/PhDThesis/GBS/filtering/n1")

# populations -P ./denovo_map/stacks-n1/ --popmap ./pop_maps/N1_popmap_colony_site.txt -t 14 --vcf -r 70 -p 3 --min-maf 0.05

# minimum percentage of individuals in a population required to process a locus: 70%
# max missing 0.5, minDP 3

system("vcftools --vcf n1.70.snps.vcf  --max-missing 0.50 --mac 3  --recode --recode-INFO-all --out raw.mm50mac3")

# minDP 
system("vcftools --vcf raw.mm50mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.mm50mac3dp3")

# minimum mean depth
system("vcftools --vcf raw.mm50mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm50mac3dp3mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm50mac3dp3mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm50mac3dp3mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out f1.1")

# f.1.recode.vcf
# 19 out of 32 Individuals
# After filtering, kept 3599 sites

###### F1 - filter 2, -r 70, -max-missing 0.5, minDP 5 #########
# minimum percentage of individuals in a population required to process a locus: 70%
# max missing 0.5, minDP 5

system("vcftools --vcf n1.70.snps.vcf  --max-missing 0.50 --mac 3  --recode --recode-INFO-all --out raw.mm50mac3")

# minDP 
system("vcftools --vcf raw.mm50mac3.recode.vcf --minDP 5 --recode --recode-INFO-all --out raw.mm50mac3dp5")

# minimum mean depth
system("vcftools --vcf raw.mm50mac3dp5.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm50mac3dp5mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm50mac3dp5mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm50mac3dp5mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out f1.2")

# f1.2.recode.vcf
# 19 out of 138 Individuals
# After filtering, kept 3599 Sites

###### F1 filter 3, -r 70, -max-missing 0.6, minDP 3 ########

# minimum percentage of individuals in a population required to process a locus: 70%
# max missing 0.6, minDP 3

system("vcftools --vcf n1.70.snps.vcf  --max-missing 0.60 --mac 3  --recode --recode-INFO-all --out raw.mm60mac3")

# minDP 
system("vcftools --vcf raw.mm60mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.mm60mac3dp3")

# minimum mean depth
system("vcftools --vcf raw.mm60mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm60mac3dp3mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm60mac3dp3mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm60mac3dp3mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out f1.3")

# f1.3.recode.vcf
# 22 out of 32 Individuals
# After filtering, kept 491 Sites

##### F1 filter 4, -r 70, -max-misssing 0.6, -minDP 5 ########
# minimum percentage of individuals in a population required to process a locus: 70%
# max missing 0.6, minDP 5

system("vcftools --vcf n1.70.snps.vcf  --max-missing 0.60 --mac 3  --recode --recode-INFO-all --out raw.mm60mac3")

# minDP 
system("vcftools --vcf raw.mm60mac3.recode.vcf --minDP 5 --recode --recode-INFO-all --out raw.mm60mac3dp5")

# minimum mean depth
system("vcftools --vcf raw.mm60mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm60mac3dp5mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm60mac3dp5mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm60mac3dp5mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out f1.4")

# f1.4.recode.vcf
# 22 out of 32 Individuals
# After filtering, kept 491 Sites


###### F1 filter 5, -r 70, -max-missing 40, minDP 3 ########

# minimum percentage of individuals in a population required to process a locus: 70%
# max missing 0.4, minDP 3

system("vcftools --vcf n1.70.snps.vcf  --max-missing 0.40 --mac 3  --recode --recode-INFO-all --out raw.mm40mac3")

# minDP 
system("vcftools --vcf raw.mm40mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.mm40mac3dp3")

# minimum mean depth
system("vcftools --vcf raw.mm40mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm40mac3dp3mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm40mac3dp3mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm40mac3dp3mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out f1.5")

# f1.5.recode.vcf
# 19 out of 31 Individuals
# After filtering, kept 7400 out of a possible 7400 Sites

###### F1 filter 6, -r 70, -max-missing 0.4, -minDP 5 ########
# minimum percentage of individuals in a population required to process a locus: 70%
# max missing 0.4, minDP 5

system("vcftools --vcf n1.70.snps.vcf  --max-missing 0.40 --mac 3  --recode --recode-INFO-all --out raw.mm40mac3")

# minDP 
system("vcftools --vcf raw.mm40mac3.recode.vcf --minDP 5 --recode --recode-INFO-all --out raw.mm40mac3dp5")

# minimum mean depth
system("vcftools --vcf raw.mm40mac3dp5.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm40mac3dp5mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm40mac3dp5mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm40mac3dp5mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out f1.6")

# f1.6.recode.vcf
# 18 out of 31 Individuals
# After filtering, kept 7400 out of a possible 7400 Sites


############## F1 filter 7, -r 80, -max-missing 0.5, -minDP 3 ###########

# populations -P ./denovo_map/stacks-n1/ --popmap ./pop_maps/N1_popmap_colony_site.txt -t 14 --vcf -r 80 -p 3 --min-maf 0.05

#minimum percentage of individuals in a population required to process a locus: 80%
# max missing 0.5, minDP 3

system("vcftools --vcf n1.80.snps.vcf  --max-missing 0.50 --mac 3  --recode --recode-INFO-all --out raw.mm50mac3")

# minDP 
system("vcftools --vcf raw.mm50mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.mm50mac3dp3")

# minimum mean depth
system("vcftools --vcf raw.mm50mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm50mac3dp3mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm50mac3dp3mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm50mac3dp3mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out f1.7")

# f1.7.recode.vcf
# 16 out of 32 Individuals
# After filtering, kept 188 Sites


######## F1 filter 8, -r 80, -max-missing 0.5, -minDP 5 ########

# minimum percentage of individuals in a population required to process a locus: 80%
# max missing 0.5, minDP 5

system("vcftools --vcf n1.80.snps.vcf  --max-missing 0.50 --mac 3  --recode --recode-INFO-all --out raw.mm50mac3")

# minDP 
system("vcftools --vcf raw.mm50mac3.recode.vcf --minDP 5 --recode --recode-INFO-all --out raw.mm50mac3dp5")

# minimum mean depth
system("vcftools --vcf raw.mm50mac3dp5.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm50mac3dp5mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm50mac3dp5mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm50mac3dp5mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out f1.8")

# f1.8.recode.vcf
# 16 out of 138 Individuals
# After filtering, kept 188 Sites



###### F1 filter 9, -r 80, -max-missing 0.6, -minDP 5 ########

# minimum percentage of individuals in a population required to process a locus: 80%
# max missing 0.6, minDP 3

system("vcftools --vcf n1.80.snps.vcf  --max-missing 0.60 --mac 3  --recode --recode-INFO-all --out raw.mm60mac3")

# minDP 
system("vcftools --vcf raw.mm60mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.mm60mac3dp3")

# minimum mean depth
system("vcftools --vcf raw.mm60mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm60mac3dp3mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm60mac3dp3mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm60mac3dp3mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out f11.9")

# f11.9.recode.vcf
# 22 out of 32 Individuals
# After filtering, kept 19 out of a possible 19 Sites




###### F1 filter 10, -r 80, -max-missing 0.6, -minDP 5 ########

# minimum percentage of individuals in a population required to process a locus: 80%
# max missing 0.6, minDP 5

system("vcftools --vcf n1.80.snps.vcf  --max-missing 0.60 --mac 3  --recode --recode-INFO-all --out raw.mm60mac3")

# minDP 
system("vcftools --vcf raw.mm60mac3.recode.vcf --minDP 5 --recode --recode-INFO-all --out raw.mm60mac3dp5")

# minimum mean depth
system("vcftools --vcf raw.mm60mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm60mac3dp5mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm60mac3dp5mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm60mac3dp5mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out f1.10")

# f1.10.recode.vcf
# 22 out of 32 Individuals
# After filtering, kept 19 out of a possible 19 Sites



###### f1 filter 11, -r 80, -max-missing 0.4, -minDP 3 ########

# minimum percentage of individuals in a population required to process a locus: 80%
# max missing 0.4, minDP 3

system("vcftools --vcf n1.80.snps.vcf  --max-missing 0.40 --mac 3  --recode --recode-INFO-all --out raw.mm40mac3")

# minDP 
system("vcftools --vcf raw.mm40mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.mm40mac3dp3")

# minimum mean depth
system("vcftools --vcf raw.mm40mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm40mac3dp3mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm40mac3dp3mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm40mac3dp3mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out f1.11")

# f1.11.recode.vcf
# 13 out of 138 Individuals
# After filtering, kept 1051 Sites



###### f1 filter 12, -r 80, -max-missing 0.4, min-DP 5 ########

# minimum percentage of individuals in a population required to process a locus: 80%
# max missing 0.4, minDP 5

system("vcftools --vcf n1.80.snps.vcf  --max-missing 0.40 --mac 3  --recode --recode-INFO-all --out raw.mm40mac3")

# minDP 
system("vcftools --vcf raw.mm40mac3.recode.vcf --minDP 5 --recode --recode-INFO-all --out raw.mm40mac3dp5")

# minimum mean depth
system("vcftools --vcf raw.mm40mac3dp5.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm40mac3dp5mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm40mac3dp5mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm40mac3dp5mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out f1.12")

# f1.12.recode.vcf
# 13 out of 32 Individuals
# After filtering, kept 1051 Sites


###### F1 filter 13, -r 60, -max-missing 0.5, -minDP 3 ######


# populations -P ./denovo_map/stacks-n1/ --popmap ./pop_maps/N1_popmap_colony_site.txt -t 14 --vcf -r 60 -p 3 --min-maf 0.05

#minimum percentage of individuals in a population required to process a locus: 60%
# max missing 0.5, minDP 3

system("vcftools --vcf n1.60.snps.vcf  --max-missing 0.50 --mac 3  --recode --recode-INFO-all --out raw.mm50mac3")

# minDP 
system("vcftools --vcf raw.mm50mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.mm50mac3dp3")

# minimum mean depth
system("vcftools --vcf raw.mm50mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm50mac3dp3mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm50mac3dp3mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm50mac3dp3mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out f1.13")

# f1.recode.vcf
# 20 out of 32 Individuals
# After filtering, kept 6578 Sites

######### F1 - filter 14, -r 60, -max-missing 50, -minDP 5 ######

# minimum percentage of individuals in a population required to process a locus: 60%
# max missing 0.5, minDP 5

system("vcftools --vcf n1.60.snps.vcf  --max-missing 0.50 --mac 3  --recode --recode-INFO-all --out raw.mm50mac3")

# minDP 
system("vcftools --vcf raw.mm50mac3.recode.vcf --minDP 5 --recode --recode-INFO-all --out raw.mm50mac3dp5")

# minimum mean depth
system("vcftools --vcf raw.mm50mac3dp5.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm50mac3dp5mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm50mac3dp5mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm50mac3dp5mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out f1.14")

# f1.14.recode.vcf
# 20 out of 138 Individuals
# After filtering, kept 6578  Sites


###### F1 filter 15, -r 60, -max-missing 0.6, -minDP 3 ########
# minimum percentage of individuals in a population required to process a locus: 60%
# max missing 0.6, minDP 3

system("vcftools --vcf n1.60.snps.vcf  --max-missing 0.60 --mac 3  --recode --recode-INFO-all --out raw.mm60mac3")

# minDP 
system("vcftools --vcf raw.mm60mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.mm60mac3dp3")

# minimum mean depth
system("vcftools --vcf raw.mm60mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm60mac3dp3mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm60mac3dp3mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm60mac3dp3mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out f1.15")

# f1.15.recode.vcf
# 22 out of 32 Individuals
# After filtering, kept 1510 Sites


###### F1 filter 16, -r 60, -max-missing 0.6, -minDP 5 ########
# minimum percentage of individuals in a population required to process a locus: 60%
# max missing 0.6, minDP 5

## eximius - original submission filtering steps ##
system("vcftools --vcf n1.60.snps.vcf  --max-missing 0.60 --mac 3  --recode --recode-INFO-all --out raw.mm60mac3")

# minDP 
system("vcftools --vcf raw.mm60mac3.recode.vcf --minDP 5 --recode --recode-INFO-all --out raw.mm60mac3dp5")

# minimum mean depth
system("vcftools --vcf raw.mm60mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm60mac3dp5mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm60mac3dp5mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm60mac3dp5mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out f1.16")

# f1.16.recode.vcf
# 22 out of 32 Individuals
# After filtering, kept 1510 Sites



###### f1 filter 17, -r 60, -max-missing 0.4, -minDP 3 ########
# minimum percentage of individuals in a population required to process a locus: 60%
# max missing 0.4, minDP 3

system("vcftools --vcf n1.60.snps.vcf  --max-missing 0.40 --mac 3  --recode --recode-INFO-all --out raw.mm40mac3")

# minDP 
system("vcftools --vcf raw.mm40mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.mm40mac3dp3")

# minimum mean depth
system("vcftools --vcf raw.mm40mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm40mac3dp3mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm40mac3dp3mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm40mac3dp3mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out f1.17")

# f1.17.recode.vcf
# 19 out of 138 Individuals
# After filtering, kept 12098  Sites


###### F1 filter 18, -r 60, -max-missing 0.4, -minDP 5 ########

# minimum percentage of individuals in a population required to process a locus: 60%
# max missing 0.4, minDP 5

system("vcftools --vcf n1.60.snps.vcf  --max-missing 0.40 --mac 3  --recode --recode-INFO-all --out raw.mm40mac3")

# minDP 
system("vcftools --vcf raw.mm40mac3.recode.vcf --minDP 5 --recode --recode-INFO-all --out raw.mm40mac3dp5")

# minimum mean depth
system("vcftools --vcf raw.mm40mac3dp5.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm40mac3dp5mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm40mac3dp5mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm40mac3dp5mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out f1.18")

# f1.18.recode.vcf
# 15 out of 32 Individuals
# After filtering, kept 12098 Sites

##########
# use plink to convert to bed format
system("sed '16,$s/^/scaffold_/' filtered.recode.vcf > filtered.recode.scaffolds.vcf")
# awk '{print "chr"$0}' file.in > file.out

# # replace chrom column with zeroes
# system("sed '16,$s/[0-9.]*/0/' filtered.recode.vcf > filtered.recode.scaffolds.zeroed.vcf")

system("./plink2 --vcf filtered.recode.scaffolds.vcf --allow-extra-chr 0 --double-id --make-bed --out scaffolds")


# admixture test
system("./admixture --cv scaffolds.bed 2 > log2.out")
##############
# 
# #### and repeat for N1 #####
# vcftools --vcf populations.snps.vcf  --max-missing 0.50 --mac 3 --maf 0.05 --recode --recode-INFO-all --out raw.mm50mac3maf05
# # After filtering, kept 3636 out of a possible 80117 Sites, that's too many
# 
# vcftools --vcf populations.snps.vcf  --max-missing 0.40 --mac 3 --maf 0.05 --recode --recode-INFO-all --out raw.mm40mac3maf05
# # After filtering, kept 8417 out of a possible 80117 Sites, better. we'll see
# 
# # minDP 
# vcftools --vcf raw.mm40mac3maf05.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.mm40mac3maf05dp3 
# # kept 8417 out of a possible 8417 Sites
# 
# # min-mean DP
# vcftools --vcf raw.mm40mac3maf05dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm40mac3maf05dp3mdp5
# # kept 7400 out of a possible 8417 Sites
# 
# 
# # remove missing individuals
# vcftools --vcf raw.mm40mac3maf05dp3mdp5.recode.vcf --missing-indv
# # oh yikes, some individuals have a LOT of missing data 
# 
# # filter out any that have more than 0.7 missing data
# awk '$5 > 0.7' out.imiss | cut -f1 > lowDP.indv
# 
# 
# # and remove
# vcftools --vcf raw.mm40mac3maf05dp3mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out filtered
# 
# # After filtering, kept 19 out of 32 Individuals
# # After filtering, kept 7400 out of a possible 7400 Sites
# 
# 
# #### going for F2 ####
# # first, make vcf file
# populations -P ./denovo_map/stacks-f2/ --popmap ./pop_maps/F2_popmap_colony_site.txt -t 14 --vcf
# 
# 
# vcftools --vcf populations.snps.vcf  --max-missing 0.40 --mac 3 --maf 0.05 --recode --recode-INFO-all --out raw.mm40mac3maf05
# # After filtering, kept 17574 out of a possible 175555 Sites
# 
# 
# # minDP 
# vcftools --vcf raw.mm40mac3maf05.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.mm40mac3maf05dp3 
# # After filtering, kept 17574 out of a possible 17574 Sites
# 
# 
# # min-mean DP
# vcftools --vcf raw.mm40mac3maf05dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm40mac3maf05dp3mdp5
# #After filtering, kept 17217 out of a possible 17574 Sites
# 
# 
# # remove missing individuals
# vcftools --vcf raw.mm40mac3maf05dp3mdp5.recode.vcf --missing-indv
# # woah, they're either above 0.75 or below 0.25
# 
# 
# # filter out any that have more than 0.5 missing data
# awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv
# 
# 
# # and remove
# vcftools --vcf raw.mm40mac3maf05dp3mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out filtered
# # After filtering, kept 12 out of 29 Individuals, lol
# # After filtering, kept 17217 out of a possible 17217 Sites
# 
