## Clear Global Environment
rm(list = ls())

#### Session prep ####
## Install packages and load libraries as required
if(!require(tidyverse)){
  install.packages("tidyverse",dependencies = T)
  library(tidyverse)
}
if(!require(lubridate)){
  install.packages("lubridate",dependencies = T)
  library(lubridate)
}
if(!require(ggplot2)){
  install.packages("ggplot2",dependencies = T)
  library(ggplot2)
}
if(!require(ggrepel)){
  install.packages("ggrepel",dependencies = T)
  library(ggrepel)
}
if(!require(grid)){
  install.packages("grid",dependencies = T)
  library(grid)
}
if(!require(gridExtra)){
  install.packages("gridExtra",dependencies = T)
  library(gridExtra)
}
if(!require(cowplot)){
  install.packages("cowplot")
  library(cowplot)
}
if(!require(reshape2)){
  install.packages("reshape2",dependencies = T)
  library(reshape2)
}
if(!require(ggsignif)){
  install.packages("ggsignif",dependencies = T)
  library(ggsignif)
}
if(!require(vcfR)){
  install.packages("vcfR")
  library(vcfR)
}

#### Functions ####
## data summary
ds <- function(data, varname, groupnames) {
  require(plyr)
  summary_func <- function(x, col) {
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se = sum(sd(x[[col]]) / sqrt(length(x[[1]]))))
  }
  data_sum <- ddply(data, groupnames, .fun=summary_func, varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
# examples:
# x <- as.data.frame(ds(df, 
#                       varname = "varname_of_interest", 
#                       groupnames = c("group_of_interest")))
# x <- as.data.frame(ds(df, 
#                       varname = "varname_of_interest", 
#                       groupnames = c("group1", "group2", "...")))

## not in
'%!in%' <- function(x,y)!('%in%'(x,y))
# useful for removing column values of y from df x
# e.g., x <- filter(x,y %!in% c("unwantedValue1", "unwantedValue2"))


#### Plot standards ####
## factors
#factor_by_gene <- c("PB2", "PB1", "PA", "HA", "HA_antigenic", "HA_nonantigenic", "NP", "Neuraminidase", "M1", "M2", "NS1", "NEP", "PA-X", "PB1-F2")
#factor_by_segment <- c("PB2", "PB1", "PA", "HA", "NP", "Neuraminidase", "M1", "M2", "NS1", "NEP")
#factor_by_gene_cleaned <- c("PB2", "PB1", "PA", "HA", "Anti. HA", "Nonanti. HA", "NP", "NA", "M1", "M2", "NS1", "NEP", "PA-X", "PB1-F2"   )

## theme
axis_formatting <- theme(axis.text.x = element_text(size = 6),
                         axis.text.y = element_text(size = 6),
                         axis.title.x = element_text(size = 6, margin = margin(t = 4)),
                         axis.title.y = element_text(size = 6, margin = margin(r = 4)))

legend_formatting <- theme(legend.text = element_text(size = 6),
                           legend.key.height= unit(0.5, 'cm'),
                           legend.key.width= unit(0.5, 'cm'))

background_formatting <- theme(panel.border = element_rect(color = "grey", fill = NA, size = .5),
                               panel.grid = element_blank(),
                               strip.background = element_blank(),
                               panel.background = element_blank(),
                               legend.background = element_blank())

## plot_grid
Size_adjust = 12
LR_adjust = -0.5 # less = right
UD_adjust = 1.1 # less = down 

## palettes
palette_muts <- c("Nonsynonymous" = "#FF7F20",
                  "Synonymous" = "#4F7899",
                  "Missense" = "#7AD9C2")
# mutation types without stop
palette_muts_NS_S <- c("Nonsynonymous" = "#FF7F20",
                       "Synonymous" = "#4F7899")

#### Import and clean ####
dir_09_consensus_VCF <- paste("/Users/rieshunter/Documents/bioinformatics/Wolbachia/data/reads/data/run/09_consensus_VCF", sep="")
dir_09_reference_VCF <- paste("/Users/rieshunter/Documents/bioinformatics/Wolbachia/data/reads/data/run/09_reference_VCF", sep="")
dir_10_snpdat_TSV <- paste("/Users/rieshunter/Documents/bioinformatics/Wolbachia/data/reads/data/run/10_snpdat_TSV", sep="")
dir_10_snpdat_TXT <- paste("/Users/rieshunter/Documents/bioinformatics/Wolbachia/data/reads/data/run/10_snpdat_TXT", sep="")
dir_12_snpgenie <- paste("/Users/rieshunter/Documents/bioinformatics/Wolbachia/data/reads/data/run/12_snpgenie", sep="")
dir_14_R <- paste("/Users/rieshunter/Documents/bioinformatics/Wolbachia/data/reads/data/run/14_R", sep="")
dir_save <- paste("/Users/rieshunter/Documents/bioinformatics/Wolbachia/data/reads/data/run", sep="")

#### dir_09_consensus_VCF ####
setwd(dir_09_consensus_VCF)
vcf <- dir(pattern="_L001.vcf")
names_trunc <- gsub("_L001.vcf","",vcf)
names_trunc <- gsub("09-consensus_","",names_trunc)
n <- length(vcf)
list <- vector("list",n)
for (i in 1:n) {
  list[[i]] <- read.vcfR(vcf[i], verbose=T)
  list[[i]] <- as.data.frame(list[[i]]@fix)
  ifelse(length(list[[i]]$CHROM)>0, 
         list[[i]]$FILTER <- names_trunc[i],
         print("No variants!"))
  names(list) <- names_trunc}
df_09_consensus_VCF <- Reduce(full_join,list)

#### dir_09_reference_VCF ####
setwd(dir_09_reference_VCF)
vcf <- dir(pattern="_L001.vcf")
names_trunc <- gsub("_L001.vcf","",vcf)
names_trunc <- gsub("09-reference_","",names_trunc)
n <- length(vcf)
list <- vector("list",n)
for (i in 1:n) {
  list[[i]] <- read.vcfR(vcf[i], verbose=T)
  list[[i]] <- as.data.frame(list[[i]]@fix)
  ifelse(length(list[[i]]$CHROM)>0, 
         list[[i]]$FILTER <- names_trunc[i],
         print("No variants!"))
  names(list) <- names_trunc}
df_09_reference_VCF <- Reduce(full_join,list)

#### dir_10_snpdat_TSV ####


#### dir_10_snpdat_TXT ####


#### dir_12_snpgenie ####
#codon_results
setwd(dir_12_snpgenie)
sg_cr <- dir(pattern="_L001_sg_codon_results.txt")
names_trunc <- gsub("_L001_sg_codon_results.txt","",sg_cr)
names_trunc <- gsub("12-","",names_trunc)
n <- length(sg_cr)
list <- vector("list",n)
for (i in 1:n) {
  list[[i]] <- as.data.frame(
    read_tsv(sg_cr[[i]]))
  names(list) <- names_trunc}
df_12_snpgenie_sg_cr <- Reduce(full_join,list)
#product_results
sg_pr <- dir(pattern="_L001_sg_product_results.txt")
names_trunc <- gsub("_L001_sg_product_results.txt","",sg_pr)
names_trunc <- gsub("12-","",names_trunc)
n <- length(sg_pr)
list <- vector("list",n)
for (i in 1:n) {
  list[[i]] <- as.data.frame(
    read_tsv(sg_pr[[i]]))
  names(list) <- names_trunc
  list[[i]]$mean_gdiv_polymorphic <- NULL
  list[[i]]$mean_N_gdiv <- NULL
  list[[i]]$mean_S_gdiv <- NULL}
df_12_snpgenie_sg_pr <- Reduce(full_join,list)
df_12_snpgenie_sg_pr$sample <- df_12_snpgenie_sg_pr$file
df_12_snpgenie_sg_pr$sample <- gsub("./09-consensus_", "", df_12_snpgenie_sg_pr$sample)
df_12_snpgenie_sg_pr$sample <- gsub("_L001.vcf", "", df_12_snpgenie_sg_pr$sample)
df_12_snpgenie_sg_pr <- separate(df_12_snpgenie_sg_pr, "sample", c("sample","S"), sep="_S")

df_12_snpgenie_sg_pr_controls <- df_12_snpgenie_sg_pr[grepl("PC",df_12_snpgenie_sg_pr$sample),]
df_12_snpgenie_sg_pr_ZIKV <- df_12_snpgenie_sg_pr[grepl("ZIKV-seq",df_12_snpgenie_sg_pr$sample),]
df_12_snpgenie_sg_pr_dups <- df_12_snpgenie_sg_pr[grepl("dup",df_12_snpgenie_sg_pr$sample),]
df_12_snpgenie_sg_pr_NC <- df_12_snpgenie_sg_pr[grepl("NC",df_12_snpgenie_sg_pr$sample),]
df_12_snpgenie_sg_pr <- df_12_snpgenie_sg_pr[!grepl("PC",df_12_snpgenie_sg_pr$sample),]
df_12_snpgenie_sg_pr <- df_12_snpgenie_sg_pr[!grepl("ZIKV-seq",df_12_snpgenie_sg_pr$sample),]
df_12_snpgenie_sg_pr <- df_12_snpgenie_sg_pr[!grepl("dup",df_12_snpgenie_sg_pr$sample),]
df_12_snpgenie_sg_pr <- df_12_snpgenie_sg_pr[!grepl("NC",df_12_snpgenie_sg_pr$sample),]

df_12_snpgenie_sg_pr_dups <- separate(df_12_snpgenie_sg_pr_dups, "sample", c("1","2","dup","3","4","5","6"), sep="-")
df_12_snpgenie_sg_pr <- separate(df_12_snpgenie_sg_pr, "sample", c("1","2","3","4","5","6"), sep="-")
df_12_snpgenie_sg_pr$dup <- "not_dup"
df_12_snpgenie_sg_pr <- rbind(df_12_snpgenie_sg_pr, df_12_snpgenie_sg_pr_dups)
df_12_snpgenie_sg_pr$ID <- paste(df_12_snpgenie_sg_pr$`1`,df_12_snpgenie_sg_pr$`2`, sep = "_")
df_12_snpgenie_sg_pr$dpi <- as.integer(df_12_snpgenie_sg_pr$`3`)
df_12_snpgenie_sg_pr$rep <- df_12_snpgenie_sg_pr$`4`
df_12_snpgenie_sg_pr$group <- df_12_snpgenie_sg_pr$`5`
df_12_snpgenie_sg_pr$sample_type <- df_12_snpgenie_sg_pr$`6`
df_12_snpgenie_sg_pr$group_type <- paste(df_12_snpgenie_sg_pr$group, df_12_snpgenie_sg_pr$sample_type, sep = "_")
df_12_snpgenie_sg_pr$`1` <- NULL
df_12_snpgenie_sg_pr$`2` <- NULL
df_12_snpgenie_sg_pr$`3` <- NULL
df_12_snpgenie_sg_pr$`4` <- NULL
df_12_snpgenie_sg_pr$`5` <- NULL
df_12_snpgenie_sg_pr$`6` <- NULL

df_12_snpgenie_sg_pr$pi <- df_12_snpgenie_sg_pr$piN + df_12_snpgenie_sg_pr$piS
df_12_snpgenie_sg_pr$piNpiS <- df_12_snpgenie_sg_pr$piN / df_12_snpgenie_sg_pr$piS
df_12_snpgenie_sg_pr$piNminuspiS <- df_12_snpgenie_sg_pr$piN - df_12_snpgenie_sg_pr$piS

#sliding_window
sg_sw <- dir(pattern="_L001_sg_sliding_window.txt")
names_trunc <- gsub("_L001_sg_sliding_window.txt","",sg_sw)
names_trunc <- gsub("12-","",names_trunc)
n <- length(sg_sw)
list <- vector("list",n)
for (i in 1:n) {
  list[[i]] <- as.data.frame(
    read_tsv(sg_sw[[i]]))
  names(list) <- names_trunc}
sg_sw <- Reduce(full_join,list)
sg_sw$pi <- sg_sw$piN + sg_sw$piS

sg_sw$sample <- sg_sw$file
sg_sw$sample <- gsub("./09-consensus_", "", sg_sw$sample)
sg_sw$sample <- gsub("_L001.vcf", "", sg_sw$sample)
sg_sw <- separate(sg_sw, "sample", c("sample","S"), sep="_S")

sg_sw_controls <- sg_sw[grepl("PC",sg_sw$sample),]
sg_sw_ZIKV <- sg_sw[grepl("ZIKV-seq",sg_sw$sample),]
sg_sw_dups <- sg_sw[grepl("dup",sg_sw$sample),]
sg_sw_NC <- sg_sw[grepl("NC",sg_sw$sample),]
sg_sw <- sg_sw[!grepl("PC",sg_sw$sample),]
sg_sw <- sg_sw[!grepl("ZIKV-seq",sg_sw$sample),]
sg_sw <- sg_sw[!grepl("dup",sg_sw$sample),]
sg_sw <- sg_sw[!grepl("NC",sg_sw$sample),]

sg_sw_dups <- separate(sg_sw_dups, "sample", c("1","2","dup","3","4","5","6"), sep="-")
sg_sw <- separate(sg_sw, "sample", c("1","2","3","4","5","6"), sep="-")
sg_sw$dup <- "not_dup"
sg_sw <- rbind(sg_sw, sg_sw_dups)
sg_sw$ID <- paste(sg_sw$`1`,sg_sw$`2`, sep = "_")
sg_sw$dpi <- as.integer(sg_sw$`3`)
sg_sw$rep <- sg_sw$`4`
sg_sw$group <- sg_sw$`5`
sg_sw$sample_type <- sg_sw$`6`
sg_sw$group_type <- paste(sg_sw$group, sg_sw$sample_type, sep = "_")
sg_sw$`1` <- NULL
sg_sw$`2` <- NULL
sg_sw$`3` <- NULL
sg_sw$`4` <- NULL
sg_sw$`5` <- NULL
sg_sw$`6` <- NULL
levels(as.factor(as.character(sg_sw$product)))
sg_sw$product <- as.factor(sg_sw$product)
sg_sw$file <- as.factor(sg_sw$file)
sg_sw$group_type <- as.factor(sg_sw$group_type)

#### prelim sw plots ####
ds_sg_sw_pi <- ds(sg_sw, varname="pi", groupnames=c("first_site"))
ds_sg_sw_piN <- ds(sg_sw, varname="piN", groupnames=c("first_site"))
ds_sg_sw_piS <- ds(sg_sw, varname="piS", groupnames=c("first_site"))
ds_sg_sw <- merge(ds_sg_sw_piN, ds_sg_sw_piS, by = "first_site")
ds_sg_sw <- merge(ds_sg_sw, ds_sg_sw_pi, by = "first_site")
# calculate sliding window piNpiS
ds_sg_sw <- transform(ds_sg_sw, piNpiS=piN/piS)
# molten
molten_ds_sg_sw <- melt(ds_sg_sw,
                        id.vars = c("first_site"),
                        measure.vars = c("piN", "piS"),
                        variable.name = "Pi")
molten_ds_sg_sw$Pi <- gsub("piN", "Nonsynonymous", molten_ds_sg_sw$Pi)
molten_ds_sg_sw$Pi <- gsub("piS", "Synonymous", molten_ds_sg_sw$Pi)
temp1 <- ggplot(data = molten_ds_sg_sw, aes(x = first_site, y = value, group = Pi, fill = Pi)) + 
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.95) + 
  scale_fill_manual(values = palette_muts_NS_S) +   
  scale_y_continuous(breaks = seq(0, .001, .0002), limits = c(0, .001),
                     labels = function(x) format(x, scientific = TRUE)) +   labs(x="", y="Pi") + 
  labs(x="", y="Nucleotide diversity") + 
  coord_cartesian(xlim = c(0, max(ds_sg_sw$first_site))) + 
  axis_formatting + background_formatting + 
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = c(.875,.9),
        legend.text=element_text(size=6),
        legend.key.size = unit(0.5,"line"))
## PiN / PiS
top_pos <- round(max(ds_sg_sw$piNpiS))+2
right_pos <- max(ds_sg_sw$first_site)
temp2 <- ggplot() + 
  geom_hline(yintercept=1, color = "black", linetype = 'dotted') + 
  ylim(0, top_pos) + 
  xlim(0, round(right_pos)+10) + 
  geom_bar(data = ds_sg_sw, aes(x = first_site, y = piNpiS),
           stat = "identity", position = position_dodge(), fill = "black") + 
  labs(x="", y="PiN / PiS") + 
  axis_formatting + legend_formatting + background_formatting
temp2.1 <- ggplot() + 
  geom_hline(yintercept=1, color = "black", linetype = 'dotted') + 
  ylim(0, top_pos) + 
  xlim(0, round(right_pos)+10) + 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  geom_bar(data = ds_sg_sw, aes(x = first_site, y = piNpiS),
           stat = "identity", position = position_dodge(), fill = "black") + 
  labs(x="", y="PiN / PiS") + 
  axis_formatting + legend_formatting + background_formatting
## PiN - PiS
ds_sg_sw <- transform(ds_sg_sw, PiNminusPiS=piN-piS)
max_value <- round(max(abs(ds_sg_sw$PiNminusPiS)),4)
right_pos2 <- max(ds_sg_sw$first_site)

temp3 <- ggplot() + 
  coord_cartesian(xlim = c(0, right_pos2), ylim = c(-max_value, max_value)) + 
  geom_hline(yintercept=0, color = "black", linetype = 'dotted') + 
  geom_bar(data = ds_sg_sw, aes(x = first_site, y = PiNminusPiS),
           stat = "identity", position = position_dodge(), fill = "black") + 
  labs(x="Nucleotide postion", y="PiN - PiS")          + 
  axis_formatting + legend_formatting + background_formatting
## Pi
top_pos3 <- .001
right_pos3 <- max(ds_sg_sw$first_site)
temp4 <- ggplot() + 
  coord_cartesian(xlim = c(0, right_pos3), ylim = c(0, top_pos3)) + 
  geom_hline(yintercept=0, color = "black", linetype = 'dotted') + 
  geom_bar(data = ds_sg_sw, aes(x = first_site, y = pi),
           stat = "identity", position = position_dodge(), fill = "black") + 
  labs(x="Nucleotide postion", y="Nucleotide diversity") + 
  scale_y_continuous(breaks = seq(0, .001, .0002), limits = c(0, .001),
                     labels = function(x) format(x, scientific = TRUE)) +   labs(x="", y="Pi") + 
  axis_formatting + legend_formatting + background_formatting
snpgenie_plot <- plot_grid(temp1, temp4, temp2, temp2.1, temp3, ncol=1, align="v", axis ="l")

#### dir_14_R ####
#R_cadm
setwd(dir_14_R)
R_cadm <- dir(pattern="_L001_coverage_and_diversity_metrics.csv")
names_trunc <- gsub("_L001_coverage_and_diversity_metrics.csv","",R_cadm)
names_trunc <- gsub("R_","",names_trunc)
n <- length(R_cadm)
list <- vector("list",n)
for (i in 1:n) {
  list[[i]] <- as.data.frame(
    read_csv(R_cadm[[i]]))
  ifelse(length(list[[i]]$RMSD)>0, 
         list[[i]]$sample <- names_trunc[i],
         print("No data!!"))
  ifelse(length(list[[i]]$sample)>0, 
         list[[i]] <- separate(as.data.frame(list[[i]]), "sample", c("sample", "S"), sep="_S"),
         print("No data!!"))
  names(list) <- names_trunc}
df_R_cadm <- Reduce(full_join,list)
df_R_cadm_controls <- df_R_cadm[grepl("PC",df_R_cadm$sample),]
df_R_cadm_ZIKV <- df_R_cadm[grepl("ZIKV-seq",df_R_cadm$sample),]
df_R_cadm_dups <- df_R_cadm[grepl("dup",df_R_cadm$sample),]
df_R_cadm <- df_R_cadm[!grepl("PC",df_R_cadm$sample),]
df_R_cadm <- df_R_cadm[!grepl("ZIKV-seq",df_R_cadm$sample),]
df_R_cadm <- df_R_cadm[!grepl("dup",df_R_cadm$sample),]
df_R_cadm_dups <- separate(df_R_cadm_dups, "sample", c("1","2","dup","3","4","5","6"), sep="-")
df_R_cadm <- separate(df_R_cadm, "sample", c("1","2","3","4","5","6"), sep="-")
df_R_cadm$dup <- "not_dup"
df_R_cadm <- rbind(df_R_cadm, df_R_cadm_dups)
df_R_cadm$ID <- paste(df_R_cadm$`1`,df_R_cadm$`2`, sep = "_")
df_R_cadm$dpi <- as.integer(df_R_cadm$`3`)
df_R_cadm$rep <- df_R_cadm$`4`
df_R_cadm$group <- df_R_cadm$`5`
df_R_cadm$sample_type <- df_R_cadm$`6`
df_R_cadm$group_type <- paste(df_R_cadm$group, df_R_cadm$sample_type, sep = "_")
df_R_cadm$`1` <- NULL
df_R_cadm$`2` <- NULL
df_R_cadm$`3` <- NULL
df_R_cadm$`4` <- NULL
df_R_cadm$`5` <- NULL
df_R_cadm$`6` <- NULL
df_R_cadm <- df_R_cadm[df_R_cadm$`%Covered`>=80,]
df_R_cadm <- df_R_cadm[df_R_cadm$MeanDepth>=2000,]
#R_msf
R_msf <- dir(pattern="_L001_mutational_spectrum_frequencies.csv")
names_trunc <- gsub("_L001_mutational_spectrum_frequencies.csv","",R_msf)
names_trunc <- gsub("R_","",names_trunc)
n <- length(R_msf)
list <- vector("list",n)
for (i in 1:n) {
  list[[i]] <- as.data.frame(
    read_csv(R_msf[[i]]))
  ifelse(length(list[[i]]$`A>U`)>0, 
         list[[i]]$sample <- names_trunc[i],
         print("No data!!"))
  names(list) <- names_trunc}
df_R_msf <- Reduce(full_join,list)
#R_ms
R_ms <- dir(pattern="_L001_mutational_spectrum.csv")
names_trunc <- gsub("_L001_mutational_spectrum.csv","",R_ms)
names_trunc <- gsub("R_","",names_trunc)
n <- length(R_ms)
list <- vector("list",n)
for (i in 1:n) {
  list[[i]] <- as.data.frame(
    read_csv(R_ms[[i]]))
  ifelse(length(list[[i]]$AtoU)>0, 
         list[[i]]$sample <- names_trunc[i],
         print("No data!!"))
  names(list) <- names_trunc}
df_R_ms <- Reduce(full_join,list)

#### Plot ####
ds_df_R_cadm <- as.data.frame(ds(df_R_cadm, 
                                 varname = "Mut.Freq.per.10K", 
                                 groupnames = c("group_type")))
plot_mut <- ggplot(ds_df_R_cadm) + 
  geom_point(aes(x = group_type, y = `Mut.Freq.per.10K`)) + 
  geom_errorbar(aes(x = group_type, y = `Mut.Freq.per.10K`, 
                    ymin = `Mut.Freq.per.10K` - sd, ymax = `Mut.Freq.per.10K` + sd), 
                width = .2, position = position_dodge(.9)) + 
  theme(legend.title = element_blank(), legend.key = element_blank(), legend.position = "right") + 
  axis_formatting + legend_formatting + background_formatting

ds_df_R_cadm <- as.data.frame(ds(df_R_cadm, 
                                 varname = "Mean.Gini.Simpson", 
                                 groupnames = c("group_type")))
plot_mgs <- ggplot(ds_df_R_cadm) + 
  geom_point(aes(x = group_type, y = `Mean.Gini.Simpson`)) + 
  geom_errorbar(aes(x = group_type, y = `Mean.Gini.Simpson`, 
                    ymin = `Mean.Gini.Simpson` - sd, ymax = `Mean.Gini.Simpson` + sd), 
                width = .2, position = position_dodge(.9)) + 
  theme(legend.title = element_blank(), legend.key = element_blank(), legend.position = "right") + 
  axis_formatting + legend_formatting + background_formatting

ds_df_R_cadm <- as.data.frame(ds(df_R_cadm, 
                                 varname = "Mean.Shannon.Entropy", 
                                 groupnames = c("group_type")))
plot_mse <- ggplot(ds_df_R_cadm) + 
  geom_point(aes(x = group_type, y = `Mean.Shannon.Entropy`)) + 
  geom_errorbar(aes(x = group_type, y = `Mean.Shannon.Entropy`, 
                    ymin = `Mean.Shannon.Entropy` - sd, ymax = `Mean.Shannon.Entropy` + sd), 
                width = .2, position = position_dodge(.9)) + 
  theme(legend.title = element_blank(), legend.key = element_blank(), legend.position = "right") + 
  axis_formatting + legend_formatting + background_formatting

## df_12_snpgenie_sg_pr
ds_12_snpgenie_sg_pr_Pi <- as.data.frame(ds(df_12_snpgenie_sg_pr, varname = "pi", groupnames = c("group_type")))
ds_12_snpgenie_sg_pr_PiN <- as.data.frame(ds(df_12_snpgenie_sg_pr, varname = "piN", groupnames = c("group_type")))
ds_12_snpgenie_sg_pr_PiS <- as.data.frame(ds(df_12_snpgenie_sg_pr, varname = "piS", groupnames = c("group_type")))
ds_12_snpgenie_sg_pr_PiNPiS <- as.data.frame(ds(df_12_snpgenie_sg_pr, varname = "piNpiS", groupnames = c("group_type")))
ds_12_snpgenie_sg_pr_PiNminusPiS <- as.data.frame(ds(df_12_snpgenie_sg_pr, varname = "piNminuspiS", groupnames = c("group_type")))

## plot function
new.function <- function (dataframe, x1, y1, 
                          datasummary, x2, y2,
                          x_lab, y_lab, title_lab) {
  ggplot() +
    geom_point(data = dataframe, aes(x = x1, y = y1), color = "black", alpha = .2) + 
    geom_point(data = datasummary, aes(x = x2, y = y2)) + 
    geom_errorbar(data = datasummary, aes(x = x2, y = y2, ymin = y2 - sd, ymax = y2 + sd), 
                  width = .2, position = position_dodge(.9)) + 
    geom_signif(data = dataframe, aes(x = x1, y = y1),
                comparisons = list(c("tet_body", "wmel_body")), 
                map_signif_level = F,
                textsize = 3,
                tip_length = 0,
                vjust = 0,
                margin_top = 0.1,
                test="t.test") + 
    geom_signif(data = dataframe, aes(x = x1, y = y1),
                comparisons = list(c("tet_legs", "wmel_legs")), 
                map_signif_level = F,
                textsize = 3,
                tip_length = 0,
                vjust = 0,
                margin_top = 0.1,
                test="t.test") + 
    labs(x = x_lab, y = y_lab, title = title_lab) + 
    theme(legend.title = element_blank(), legend.key = element_blank(), legend.position = "right") + 
    axis_formatting + legend_formatting + background_formatting
}

plot_1 <- new.function(df_12_snpgenie_sg_pr, df_12_snpgenie_sg_pr$group_type, 
                       df_12_snpgenie_sg_pr$pi, 
                       ds_12_snpgenie_sg_pr_Pi, 
                       ds_12_snpgenie_sg_pr_Pi$group_type,
                       ds_12_snpgenie_sg_pr_Pi$pi,
                       "", "Nucleotide diversity", "")

plot_2 <- new.function(df_12_snpgenie_sg_pr, df_12_snpgenie_sg_pr$group_type, 
                       df_12_snpgenie_sg_pr$piNminuspiS, 
                       ds_12_snpgenie_sg_pr_PiNminusPiS, 
                       ds_12_snpgenie_sg_pr_PiNminusPiS$group_type,
                       ds_12_snpgenie_sg_pr_PiNminusPiS$piNminuspiS,
                       "", "PiN - PiS", "")

plot_3 <- new.function(df_12_snpgenie_sg_pr, df_12_snpgenie_sg_pr$group_type, 
                       df_12_snpgenie_sg_pr$piNpiS, 
                       ds_12_snpgenie_sg_pr_PiNPiS, 
                       ds_12_snpgenie_sg_pr_PiNPiS$group_type,
                       ds_12_snpgenie_sg_pr_PiNPiS$piNpiS,
                       "", "PiNPiS", "")

sg_pr_plot <- plot_grid(plot_1, plot_2, plot_3, ncol=2, align="v", axis ="l")


## diversity
plot_mut #mutation
plot_mgs #mean gini simpson
plot_mse #mean shannon entropy

snpgenie_plot