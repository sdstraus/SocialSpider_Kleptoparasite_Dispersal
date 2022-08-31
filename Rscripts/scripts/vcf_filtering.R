
####### F1 - filter 1 #########
# populations -P ./denovo_map/stacks-n1/ --popmap ./pop_maps/N1_popmap_colony_site.txt -t 14 --vcf -r 70 -p 3 --min-maf 0.05
# minimum percentage of individuals in a population required to process a locus: 70%
# max missing 0.5, minDP 3

## eximius - original submission filtering steps ##
system("vcftools --vcf eximius.snps.vcf  --max-missing 0.50 --mac 3  --recode --recode-INFO-all --out raw.mm50mac3")

# minDP 
system("vcftools --vcf raw.mm50mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.mm50mac3dp3")

# minimum mean depth
system("vcftools --vcf raw.mm50mac3dp3.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out raw.mm50mac3dp3mdp5")

# remove missing individuals
system("vcftools --vcf raw.mm50mac3dp3mdp5.recode.vcf --missing-indv")

# filter out any that have more than 0.5 missing data
system("awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv")

# and remove
system("vcftools --vcf raw.mm50mac3dp3mdp5.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out exi.1")

# f.1.recode.vcf
# 99 out of 138 Individuals
# After filtering, kept 11087 out of a possible 11087 Sites




# populations -P ./denovo_map/stacks-n1/ --popmap ./pop_maps/N1_popmap_colony_site.txt -t 14 --vcf -r 80 -p 3 --min-maf 0.05


# populations -P ./denovo_map/stacks-n1/ --popmap ./pop_maps/N1_popmap_colony_site.txt -t 14 --vcf -r 60 -p 3 --min-maf 0.05


####### F2 - filter 1 #########

# populations -P ./denovo_map/stacks-f2/ --popmap ./pop_maps/F2_popmap_colony_site.txt -t 14 --vcf -r 70 -p 3 --min-maf 0.05

# populations -P ./denovo_map/stacks-f2/ --popmap ./pop_maps/F2_popmap_colony_site.txt -t 14 --vcf -r 80 -p 3 --min-maf 0.05

# populations -P ./denovo_map/stacks-f2/ --popmap ./pop_maps/F2_popmap_colony_site.txt -t 14 --vcf -r 60 -p 3 --min-maf 0.05


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
