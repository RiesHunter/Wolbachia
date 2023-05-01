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

if(!require(bayesboot)){
  install.packages("bayesboot")
  library(bayesboot)
}
b <- function(data){
  return(bayesboot(data, mean, R = 10000))
}

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
                  "Stop gained" = "#7AD9C2")
# mutation types without stop
palette_muts_NS_S <- c("Nonsynonymous" = "#FF7F20",
                       "Synonymous" = "#4F7899")

# gl
palette_gl <- c("Tet_saliva" = "#a97cf7", 
                "Tet_body" = "#6b42b3", 
                "Tet_legs" = "#391678", 
                "wmel_saliva" = "#6ef08a", 
                "wmel_body" = "#37a34f", 
                "wmel_legs" = "#0a5e1d")

#### Import and clean ####
dir_09_consensus_VCF <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/Wolbachia/data/reads/data/run/09_consensus_VCF/fn_ann", sep="")
dir_09_reference_VCF <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/Wolbachia/data/reads/data/run/09_reference_VCF/no_BCs", sep="")
dir_10_snpdat_TSV <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/Wolbachia/data/reads/data/run/10_snpdat_TSV", sep="")
dir_10_snpdat_TXT <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/Wolbachia/data/reads/data/run/10_snpdat_TXT", sep="")
dir_14_R <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/Wolbachia/data/reads/data/run/14_R", sep="")
dir_save <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/Wolbachia/figs", sep="")

#### dir_09_consensus_VCF ####
## these are the vcfs from the samples relative to their own consensus
#setwd(dir_09_consensus_VCF); dir()
#vcf <- dir(pattern="_L001_fn_ann.vcf")
#names_trunc <- gsub("_L001_fn_ann.vcf","",vcf)
#names_trunc <- gsub("09-consensus_","",names_trunc)
#n <- length(vcf)
#list <- vector("list",n)
#
### Read all tables in vcf, apply to list, change columns
#for (i in 1:n) {
#  list[[i]] <- read.vcfR(vcf[i], verbose=T)
#  list[[i]] <- as.data.frame(list[[i]]@fix)
#  names(list) <- names_trunc}
#
### remove blank data frames
#list <- Filter(function(x) dim(x)[1] > 0, list)
#
### separate by column
#n <- length(list)
#for (i in 1:n) {
#  list[[i]]$ID <- names_trunc[[i]]
#  
#  # separate the info column into its respective pieces
#  list[[i]] <- separate(list[[i]],"INFO",c("DP","AF","SB","DP4", "snpEff"),
#                        sep=";",convert=FALSE)
#  
#  # separate snpEff column into its respective pieces
#  list[[i]] <- separate(list[[i]],"snpEff",c("Allele","Ann","Ann_Impact","Gene_Name",
#                                             "Gene_ID","Feature_Type","Feature_ID","Transcript_BioType",
#                                             "Rank","HGVS.c","HGVS.p"),
#                        sep="\\|",convert=FALSE)
#  
#  # separate the leftover info column into the label (discarded) and value
#  list[[i]] <- separate(list[[i]],"DP",c("DP_label","DP"),sep="=",convert=FALSE)  
#  list[[i]] <- separate(list[[i]],"AF",c("AF_label","AF"),sep="=",convert=FALSE) 
#  ## remove columns by string in name
#  list[[i]] <- list[[i]] %>% select(-contains(c("_label","Gene_Name","Gene_ID","DP4","Feature_Type","Rank")))
#}
#
#df_09_consensus_VCF <- Reduce(full_join,list)
#df_09_consensus_VCF <- df_09_consensus_VCF[df_09_consensus_VCF$AF > 0.01,]
#df_09_consensus_VCF$POS <- as.integer(df_09_consensus_VCF$POS)
#df_09_consensus_VCF$AF <- as.numeric(df_09_consensus_VCF$AF)
#
#df_09_consensus_VCF$Ann <- as.factor(df_09_consensus_VCF$Ann)
#df_09_consensus_VCF$Ann <- gsub("missense_variant","Nonsynonymous",df_09_consensus_VCF$Ann)
#df_09_consensus_VCF$Ann <- gsub("synonymous_variant","Synonymous",df_09_consensus_VCF$Ann)
#df_09_consensus_VCF$Ann <- gsub("stop_gained","Stop gained",df_09_consensus_VCF$Ann)

#### dir_09_reference_VCF ####
## these are the vcfs from the samples relative to the reference
setwd(dir_09_reference_VCF); dir()
vcf <- dir(pattern="_L001_fn_ann_noBC.vcf")
names_trunc <- gsub("_L001_fn_ann_noBC.vcf","",vcf)
names_trunc <- gsub("09-reference_","",names_trunc)
n <- length(vcf)
list <- vector("list",n)

## Read all tables in vcf, apply to list, change columns
for (i in 1:n) {
  list[[i]] <- read.vcfR(vcf[i], verbose=T)
  list[[i]] <- as.data.frame(list[[i]]@fix)
  names(list) <- names_trunc}

## remove blank data frames
list <- Filter(function(x) dim(x)[1] > 0, list)

## separate by column
n <- length(list)
for (i in 1:n) {
  # separate the info column into its respective pieces
  list[[i]] <- separate(list[[i]],"INFO",c("DP","AF","SB","DP4", "snpEff"),
                        sep=";",convert=FALSE)
  list[[i]]$DP <- gsub(".*=", "", list[[i]]$DP)
  list[[i]]$AF <- gsub(".*=", "", list[[i]]$AF)
  list[[i]]$SB <- gsub(".*=", "", list[[i]]$SB)
  list[[i]]$DP4 <- gsub(".*=", "", list[[i]]$DP4)
  
  # separate snpEff column into its respective pieces
  list[[i]] <- separate(list[[i]],"snpEff",c("Allele","Ann","Ann_Impact","Gene_Name",
                                             "Gene_ID","Feature_Type","Feature_ID","Transcript_BioType",
                                             "Rank","HGVS.c","HGVS.p"),
                        sep="\\|",convert=FALSE)
  
  # remove columns by string in name
  list[[i]] <- list[[i]] %>% select(-contains(c("_label","Gene_Name","Gene_ID","DP4","Feature_Type","Rank")))
}

df_09_reference_VCF <- Reduce(full_join,list)
df_09_reference_VCF <- df_09_reference_VCF[df_09_reference_VCF$AF > 0.01,]
df_09_reference_VCF$POS <- as.integer(df_09_reference_VCF$POS)
df_09_reference_VCF$AF <- as.numeric(df_09_reference_VCF$AF)

df_09_reference_VCF$Ann <- as.factor(df_09_reference_VCF$Ann)
df_09_reference_VCF$Ann <- gsub("missense_variant","Nonsynonymous",df_09_reference_VCF$Ann)
df_09_reference_VCF$Ann <- gsub("synonymous_variant","Synonymous",df_09_reference_VCF$Ann)
df_09_reference_VCF$Ann <- gsub("stop_gained","Stop gained",df_09_reference_VCF$Ann)

#### string matches VCF ####
## consensus
#df_09_consensus_VCF$ID[grepl("7-2-tet-saliva",df_09_consensus_VCF$FILTER)] <- "Tet-saliva-7dpi_2"
#df_09_consensus_VCF$ID[grepl("7-3-tet-saliva",df_09_consensus_VCF$FILTER)] <- "Tet-saliva-7dpi_3"
#df_09_consensus_VCF$ID[grepl("14-1-tet-saliva",df_09_consensus_VCF$FILTER)] <- "Tet-saliva-14dpi_1"
#df_09_consensus_VCF$ID[grepl("14-3-tet-saliva",df_09_consensus_VCF$FILTER)] <- "Tet-saliva-14dpi_3"
#df_09_consensus_VCF$ID[grepl("7-1-tet-legs",df_09_consensus_VCF$FILTER)] <- "Tet-legs-7dpi_1"
#df_09_consensus_VCF$ID[grepl("7-2-tet-legs",df_09_consensus_VCF$FILTER)] <- "Tet-legs-7dpi_2"
#df_09_consensus_VCF$ID[grepl("7-3-tet-legs",df_09_consensus_VCF$FILTER)] <- "Tet-legs-7dpi_3"
#df_09_consensus_VCF$ID[grepl("14-1-tet-legs",df_09_consensus_VCF$FILTER)] <- "Tet-legs-14dpi_1"
#df_09_consensus_VCF$ID[grepl("14-2-tet-legs",df_09_consensus_VCF$FILTER)] <- "Tet-legs-14dpi_2"
#df_09_consensus_VCF$ID[grepl("14-3-tet-legs",df_09_consensus_VCF$FILTER)] <- "Tet-legs-14dpi_3"
#df_09_consensus_VCF$ID[grepl("4-1-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-4dpi_1"
#df_09_consensus_VCF$ID[grepl("4-2-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-4dpi_2"
#df_09_consensus_VCF$ID[grepl("4-3-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-4dpi_3"
#df_09_consensus_VCF$ID[grepl("7-1-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-7dpi_1"
#df_09_consensus_VCF$ID[grepl("7-2-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-7dpi_2"
#df_09_consensus_VCF$ID[grepl("7-3-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-7dpi_3"
#df_09_consensus_VCF$ID[grepl("dup-14-1-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-14dpi_1_dup"
#df_09_consensus_VCF$ID[grepl("14-1-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-14dpi_1"
#df_09_consensus_VCF$ID[grepl("dup-14-2-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-14dpi_2_dup"
#df_09_consensus_VCF$ID[grepl("14-2-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-14dpi_2"
#df_09_consensus_VCF$ID[grepl("dup-14-3-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-14dpi_3_dup"
#df_09_consensus_VCF$ID[grepl("14-3-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-14dpi_3"
#df_09_consensus_VCF$ID[grepl("7-1-wmel-body",df_09_consensus_VCF$FILTER)] <- "wmel-body-7dpi_1"
#df_09_consensus_VCF$ID[grepl("7-2-wmel-body",df_09_consensus_VCF$FILTER)] <- "wmel-body-7dpi_2"
#df_09_consensus_VCF$ID[grepl("7-3-wmel-body",df_09_consensus_VCF$FILTER)] <- "wmel-body-7dpi_3"
#df_09_consensus_VCF$ID[grepl("14-1-wmel-body",df_09_consensus_VCF$FILTER)] <- "wmel-body-14dpi_1"
#df_09_consensus_VCF$ID[grepl("14-2-wmel-body",df_09_consensus_VCF$FILTER)] <- "wmel-body-14dpi_2"
#df_09_consensus_VCF$ID[grepl("14-3-wmel-body",df_09_consensus_VCF$FILTER)] <- "wmel-body-14dpi_3"
#df_09_consensus_VCF$ID[grepl("7-1-wmel-legs",df_09_consensus_VCF$FILTER)] <- "wmel-legs-7dpi_1"
#df_09_consensus_VCF$ID[grepl("14-2-wmel-legs",df_09_consensus_VCF$FILTER)] <- "wmel-legs-14dpi_2"
#df_09_consensus_VCF$ID[grepl("7-2-tet-saliva",df_09_consensus_VCF$FILTER)] <- "Tet-saliva-7dpi_2"
#df_09_consensus_VCF$ID[grepl("7-3-tet-saliva",df_09_consensus_VCF$FILTER)] <- "Tet-saliva-7dpi_3"
#df_09_consensus_VCF$ID[grepl("14-1-tet-saliva",df_09_consensus_VCF$FILTER)] <- "Tet-saliva-14dpi_1"
#df_09_consensus_VCF$ID[grepl("14-3-tet-saliva",df_09_consensus_VCF$FILTER)] <- "Tet-saliva-14dpi_3"
#df_09_consensus_VCF$ID[grepl("7-1-tet-legs",df_09_consensus_VCF$FILTER)] <- "Tet-legs-7dpi_1"
#df_09_consensus_VCF$ID[grepl("7-2-tet-legs",df_09_consensus_VCF$FILTER)] <- "Tet-legs-7dpi_2"
#df_09_consensus_VCF$ID[grepl("7-3-tet-legs",df_09_consensus_VCF$FILTER)] <- "Tet-legs-7dpi_3"
#df_09_consensus_VCF$ID[grepl("14-1-tet-legs",df_09_consensus_VCF$FILTER)] <- "Tet-legs-14dpi_1"
#df_09_consensus_VCF$ID[grepl("14-2-tet-legs",df_09_consensus_VCF$FILTER)] <- "Tet-legs-14dpi_2"
#df_09_consensus_VCF$ID[grepl("14-3-tet-legs",df_09_consensus_VCF$FILTER)] <- "Tet-legs-14dpi_3"
#df_09_consensus_VCF$ID[grepl("4-1-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-4dpi_1"
#df_09_consensus_VCF$ID[grepl("4-2-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-4dpi_2"
#df_09_consensus_VCF$ID[grepl("4-3-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-4dpi_3"
#df_09_consensus_VCF$ID[grepl("7-1-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-7dpi_1"
#df_09_consensus_VCF$ID[grepl("7-2-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-7dpi_2"
#df_09_consensus_VCF$ID[grepl("7-3-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-7dpi_3"
#df_09_consensus_VCF$ID[grepl("dup-14-1-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-14dpi_1_dup"
#df_09_consensus_VCF$ID[grepl("14-1-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-14dpi_1"
#df_09_consensus_VCF$ID[grepl("dup-14-2-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-14dpi_2_dup"
#df_09_consensus_VCF$ID[grepl("14-2-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-14dpi_2"
#df_09_consensus_VCF$ID[grepl("dup-14-3-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-14dpi_3_dup"
#df_09_consensus_VCF$ID[grepl("14-3-tet-body",df_09_consensus_VCF$FILTER)] <- "Tet-body-14dpi_3"
#df_09_consensus_VCF$ID[grepl("7-1-wmel-body",df_09_consensus_VCF$FILTER)] <- "wmel-body-7dpi_1"
#df_09_consensus_VCF$ID[grepl("7-2-wmel-body",df_09_consensus_VCF$FILTER)] <- "wmel-body-7dpi_2"
#df_09_consensus_VCF$ID[grepl("7-3-wmel-body",df_09_consensus_VCF$FILTER)] <- "wmel-body-7dpi_3"
#df_09_consensus_VCF$ID[grepl("14-1-wmel-body",df_09_consensus_VCF$FILTER)] <- "wmel-body-14dpi_1"
#df_09_consensus_VCF$ID[grepl("14-2-wmel-body",df_09_consensus_VCF$FILTER)] <- "wmel-body-14dpi_2"
#df_09_consensus_VCF$ID[grepl("14-3-wmel-body",df_09_consensus_VCF$FILTER)] <- "wmel-body-14dpi_3"
#df_09_consensus_VCF$ID[grepl("7-1-wmel-legs",df_09_consensus_VCF$FILTER)] <- "wmel-legs_7dpi_1"
#df_09_consensus_VCF$ID[grepl("14-2-wmel-legs",df_09_consensus_VCF$FILTER)] <- "wmel-legs_14dpi_2"
#df_09_consensus_VCF_PC <- df_09_consensus_VCF[grepl("PC",df_09_consensus_VCF$FILTER),]
#df_09_consensus_VCF_ZIKV <- df_09_consensus_VCF[grepl("ZIKV",df_09_consensus_VCF$FILTER),]
#df_09_consensus_VCF_NC <- df_09_consensus_VCF[grepl("NC",df_09_consensus_VCF$FILTER),]
#df_09_consensus_VCF <- df_09_consensus_VCF[!grepl("PC",df_09_consensus_VCF$FILTER),]
#df_09_consensus_VCF <- df_09_consensus_VCF[!grepl("ZIKV",df_09_consensus_VCF$FILTER),]
#df_09_consensus_VCF <- df_09_consensus_VCF[!grepl("NC",df_09_consensus_VCF$FILTER),]
#df_09_consensus_VCF <- separate(df_09_consensus_VCF, "ID", c("1","2","3"), sep="-")
#df_09_consensus_VCF <- separate(df_09_consensus_VCF, "3", c("3","4"), sep="_")
#df_09_consensus_VCF$`3` <- as.integer(gsub("dpi", "", df_09_consensus_VCF$`3`))

## reference
df_09_reference_VCF$ID[grepl("7-2-tet-saliva",df_09_reference_VCF$FILTER)] <- "Tet-saliva-7dpi_2"
df_09_reference_VCF$ID[grepl("7-3-tet-saliva",df_09_reference_VCF$FILTER)] <- "Tet-saliva-7dpi_3"
df_09_reference_VCF$ID[grepl("14-1-tet-saliva",df_09_reference_VCF$FILTER)] <- "Tet-saliva-14dpi_1"
df_09_reference_VCF$ID[grepl("14-3-tet-saliva",df_09_reference_VCF$FILTER)] <- "Tet-saliva-14dpi_3"
df_09_reference_VCF$ID[grepl("7-1-tet-legs",df_09_reference_VCF$FILTER)] <- "Tet-legs-7dpi_1"
df_09_reference_VCF$ID[grepl("7-2-tet-legs",df_09_reference_VCF$FILTER)] <- "Tet-legs-7dpi_2"
df_09_reference_VCF$ID[grepl("7-3-tet-legs",df_09_reference_VCF$FILTER)] <- "Tet-legs-7dpi_3"
df_09_reference_VCF$ID[grepl("14-1-tet-legs",df_09_reference_VCF$FILTER)] <- "Tet-legs-14dpi_1"
df_09_reference_VCF$ID[grepl("14-2-tet-legs",df_09_reference_VCF$FILTER)] <- "Tet-legs-14dpi_2"
df_09_reference_VCF$ID[grepl("14-3-tet-legs",df_09_reference_VCF$FILTER)] <- "Tet-legs-14dpi_3"
df_09_reference_VCF$ID[grepl("4-1-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-4dpi_1"
df_09_reference_VCF$ID[grepl("4-2-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-4dpi_2"
df_09_reference_VCF$ID[grepl("4-3-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-4dpi_3"
df_09_reference_VCF$ID[grepl("7-1-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-7dpi_1"
df_09_reference_VCF$ID[grepl("7-2-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-7dpi_2"
df_09_reference_VCF$ID[grepl("7-3-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-7dpi_3"
df_09_reference_VCF$ID[grepl("dup-14-1-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-14dpi_1_dup"
df_09_reference_VCF$ID[grepl("14-1-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-14dpi_1"
df_09_reference_VCF$ID[grepl("dup-14-2-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-14dpi_2_dup"
df_09_reference_VCF$ID[grepl("14-2-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-14dpi_2"
df_09_reference_VCF$ID[grepl("dup-14-3-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-14dpi_3_dup"
df_09_reference_VCF$ID[grepl("14-3-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-14dpi_3"
df_09_reference_VCF$ID[grepl("7-1-wmel-body",df_09_reference_VCF$FILTER)] <- "wmel-body-7dpi_1"
df_09_reference_VCF$ID[grepl("7-2-wmel-body",df_09_reference_VCF$FILTER)] <- "wmel-body-7dpi_2"
df_09_reference_VCF$ID[grepl("7-3-wmel-body",df_09_reference_VCF$FILTER)] <- "wmel-body-7dpi_3"
df_09_reference_VCF$ID[grepl("14-1-wmel-body",df_09_reference_VCF$FILTER)] <- "wmel-body-14dpi_1"
df_09_reference_VCF$ID[grepl("14-2-wmel-body",df_09_reference_VCF$FILTER)] <- "wmel-body-14dpi_2"
df_09_reference_VCF$ID[grepl("14-3-wmel-body",df_09_reference_VCF$FILTER)] <- "wmel-body-14dpi_3"
df_09_reference_VCF$ID[grepl("7-1-wmel-legs",df_09_reference_VCF$FILTER)] <- "wmel-legs-7dpi_1"
df_09_reference_VCF$ID[grepl("14-2-wmel-legs",df_09_reference_VCF$FILTER)] <- "wmel-legs-14dpi_2"
df_09_reference_VCF$ID[grepl("7-2-tet-saliva",df_09_reference_VCF$FILTER)] <- "Tet-saliva-7dpi_2"
df_09_reference_VCF$ID[grepl("7-3-tet-saliva",df_09_reference_VCF$FILTER)] <- "Tet-saliva-7dpi_3"
df_09_reference_VCF$ID[grepl("14-1-tet-saliva",df_09_reference_VCF$FILTER)] <- "Tet-saliva-14dpi_1"
df_09_reference_VCF$ID[grepl("14-3-tet-saliva",df_09_reference_VCF$FILTER)] <- "Tet-saliva-14dpi_3"
df_09_reference_VCF$ID[grepl("7-1-tet-legs",df_09_reference_VCF$FILTER)] <- "Tet-legs-7dpi_1"
df_09_reference_VCF$ID[grepl("7-2-tet-legs",df_09_reference_VCF$FILTER)] <- "Tet-legs-7dpi_2"
df_09_reference_VCF$ID[grepl("7-3-tet-legs",df_09_reference_VCF$FILTER)] <- "Tet-legs-7dpi_3"
df_09_reference_VCF$ID[grepl("14-1-tet-legs",df_09_reference_VCF$FILTER)] <- "Tet-legs-14dpi_1"
df_09_reference_VCF$ID[grepl("14-2-tet-legs",df_09_reference_VCF$FILTER)] <- "Tet-legs-14dpi_2"
df_09_reference_VCF$ID[grepl("14-3-tet-legs",df_09_reference_VCF$FILTER)] <- "Tet-legs-14dpi_3"
df_09_reference_VCF$ID[grepl("4-1-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-4dpi_1"
df_09_reference_VCF$ID[grepl("4-2-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-4dpi_2"
df_09_reference_VCF$ID[grepl("4-3-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-4dpi_3"
df_09_reference_VCF$ID[grepl("7-1-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-7dpi_1"
df_09_reference_VCF$ID[grepl("7-2-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-7dpi_2"
df_09_reference_VCF$ID[grepl("7-3-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-7dpi_3"
df_09_reference_VCF$ID[grepl("dup-14-1-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-14dpi_1_dup"
df_09_reference_VCF$ID[grepl("14-1-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-14dpi_1"
df_09_reference_VCF$ID[grepl("dup-14-2-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-14dpi_2_dup"
df_09_reference_VCF$ID[grepl("14-2-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-14dpi_2"
df_09_reference_VCF$ID[grepl("dup-14-3-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-14dpi_3_dup"
df_09_reference_VCF$ID[grepl("14-3-tet-body",df_09_reference_VCF$FILTER)] <- "Tet-body-14dpi_3"
df_09_reference_VCF$ID[grepl("7-1-wmel-body",df_09_reference_VCF$FILTER)] <- "wmel-body-7dpi_1"
df_09_reference_VCF$ID[grepl("7-2-wmel-body",df_09_reference_VCF$FILTER)] <- "wmel-body-7dpi_2"
df_09_reference_VCF$ID[grepl("7-3-wmel-body",df_09_reference_VCF$FILTER)] <- "wmel-body-7dpi_3"
df_09_reference_VCF$ID[grepl("14-1-wmel-body",df_09_reference_VCF$FILTER)] <- "wmel-body-14dpi_1"
df_09_reference_VCF$ID[grepl("14-2-wmel-body",df_09_reference_VCF$FILTER)] <- "wmel-body-14dpi_2"
df_09_reference_VCF$ID[grepl("14-3-wmel-body",df_09_reference_VCF$FILTER)] <- "wmel-body-14dpi_3"
df_09_reference_VCF$ID[grepl("7-1-wmel-legs",df_09_reference_VCF$FILTER)] <- "wmel-legs-7dpi_1"
df_09_reference_VCF$ID[grepl("14-2-wmel-legs",df_09_reference_VCF$FILTER)] <- "wmel-legs-14dpi_2"
## remove controls
df_09_reference_VCF$FILTER <- gsub("09-reference_", "", df_09_reference_VCF$FILTER)
df_09_reference_VCF_NC <- df_09_reference_VCF[grepl("NC",df_09_reference_VCF$FILTER),]
df_09_reference_VCF_PC <- df_09_reference_VCF[grepl("PC",df_09_reference_VCF$FILTER),]
df_09_reference_VCF_ZIKV <- df_09_reference_VCF[grepl("ZIKV",df_09_reference_VCF$FILTER),]
df_09_reference_VCF_mouse <- df_09_reference_VCF[grepl("mouse",df_09_reference_VCF$FILTER),]
df_09_reference_VCF <- df_09_reference_VCF[!grepl("NC",df_09_reference_VCF$FILTER),]
df_09_reference_VCF <- df_09_reference_VCF[!grepl("PC",df_09_reference_VCF$FILTER),]
df_09_reference_VCF <- df_09_reference_VCF[!grepl("ZIKV",df_09_reference_VCF$FILTER),]
df_09_reference_VCF <- df_09_reference_VCF[!grepl("mouse",df_09_reference_VCF$FILTER),]
## separate cols
df_09_reference_VCF <- separate(df_09_reference_VCF, "ID", c("1","2","3"), sep="-")
df_09_reference_VCF <- separate(df_09_reference_VCF, "3", c("3","4"), sep="_")
df_09_reference_VCF$`3` <- as.integer(gsub("dpi", "", df_09_reference_VCF$`3`))

## bring controls back in
df_09_reference_VCF_PC$`1` <- "Control"
df_09_reference_VCF_mouse$`1` <- "Control"
df_09_reference_VCF_PC$`2` <- "PC"
df_09_reference_VCF_mouse$`2` <- "mouse"
df_09_reference_VCF_PC$`3` <- "Control"
df_09_reference_VCF_mouse$`3` <- "Control"
df_09_reference_VCF_PC$`4` <- "Control"
df_09_reference_VCF_mouse$`4` <- "Control"
df_09_reference_VCF_PC$ID <- NULL
df_09_reference_VCF_mouse$ID <- NULL

df_09_reference_VCF <- rbind(df_09_reference_VCF, 
                             df_09_reference_VCF_PC,
                             df_09_reference_VCF_mouse)

## groups
df_09_reference_VCF$group_location <- paste(df_09_reference_VCF$`1`, 
                                            df_09_reference_VCF$`2`, sep = "_")
df_09_reference_VCF$group_location_dpi <- paste(df_09_reference_VCF$`1`, 
                                                df_09_reference_VCF$`2`, 
                                                df_09_reference_VCF$`3`, 
                                                sep = "_")
df_09_reference_VCF$group_location <- as.factor(df_09_reference_VCF$group_location)
df_09_reference_VCF$group_location_dpi <- as.factor(df_09_reference_VCF$group_location_dpi)

#### enumerate iSNVs in VCF ####
## count
total_iSNVs <- as.data.frame(table(df_09_reference_VCF$FILTER))
df_09_reference_VCF$ID <- paste(df_09_reference_VCF$FILTER, df_09_reference_VCF$Ann, sep = "tk")
total_iSNVs_ann <- as.data.frame(table(df_09_reference_VCF$ID))
total_iSNVs_ann <- separate(total_iSNVs_ann, "Var1", c("sample_name", "Ann"), "tk")
total_iSNVs_ann <- total_iSNVs_ann[total_iSNVs_ann$Ann!="stop_lost&splice_region_variant",]
total_iSNVs_AF <- aggregate(df_09_reference_VCF$AF, by=list(df_09_reference_VCF$FILTER), FUN=sum)
total_iSNVs_AF_ID <- aggregate(df_09_reference_VCF$AF, by=list(df_09_reference_VCF$ID), FUN=sum)
total_iSNVs_AF_ID <- separate(total_iSNVs_AF_ID, "Group.1", c("Group.1", "Ann"), "tk")
total_iSNVs_AF_ID <- total_iSNVs_AF_ID[total_iSNVs_AF_ID$Ann!="stop_lost&splice_region_variant",]

# merge
 # total_iSNVs_AF
total_iSNVs_AF <- merge(total_iSNVs_AF, total_iSNVs_AF_ID, by = "Group.1")
names(total_iSNVs_AF) <- c("sample_name", "Total_AF", "Ann", "AF")
total_iSNVs_AF$ID <- paste(total_iSNVs_AF$sample_name, total_iSNVs_AF$Ann, sep = "tk")
total_iSNVs_AF$Ann <- NULL
total_iSNVs_AF$sample_name <- NULL
total_iSNVs_AF <- total_iSNVs_AF[,c(3,1,2)]

 # df_iSNV_enumeration
df_iSNV_enumeration <- merge(total_iSNVs_ann, total_iSNVs, by.x = "sample_name", by.y = "Var1")
names(df_iSNV_enumeration) <- c("sample_name", "Ann", "Count", "Total")
df_iSNV_enumeration$ID <- paste(df_iSNV_enumeration$sample_name, df_iSNV_enumeration$Ann, sep = "tk")
df_iSNV_enumeration$Ann <- NULL
df_iSNV_enumeration$sample_name <- NULL
df_iSNV_enumeration <- df_iSNV_enumeration[,c(3,1,2)]

 # merge
df_iSNV_enumeration <- merge(df_iSNV_enumeration, total_iSNVs_AF, by = "ID")
df_iSNV_enumeration <- separate(df_iSNV_enumeration, "ID", c("sample_name", "Ann"), "tk")
df_iSNV_enumeration$Prop <- df_iSNV_enumeration$Count/df_iSNV_enumeration$Total

## string match
df_iSNV_enumeration$ID[grepl("7-2-tet-saliva",df_iSNV_enumeration$sample_name)] <- "Tet-saliva-7dpi_2"
df_iSNV_enumeration$ID[grepl("7-3-tet-saliva",df_iSNV_enumeration$sample_name)] <- "Tet-saliva-7dpi_3"
df_iSNV_enumeration$ID[grepl("14-1-tet-saliva",df_iSNV_enumeration$sample_name)] <- "Tet-saliva-14dpi_1"
df_iSNV_enumeration$ID[grepl("14-3-tet-saliva",df_iSNV_enumeration$sample_name)] <- "Tet-saliva-14dpi_3"
df_iSNV_enumeration$ID[grepl("7-1-tet-legs",df_iSNV_enumeration$sample_name)] <- "Tet-legs-7dpi_1"
df_iSNV_enumeration$ID[grepl("7-2-tet-legs",df_iSNV_enumeration$sample_name)] <- "Tet-legs-7dpi_2"
df_iSNV_enumeration$ID[grepl("7-3-tet-legs",df_iSNV_enumeration$sample_name)] <- "Tet-legs-7dpi_3"
df_iSNV_enumeration$ID[grepl("14-1-tet-legs",df_iSNV_enumeration$sample_name)] <- "Tet-legs-14dpi_1"
df_iSNV_enumeration$ID[grepl("14-2-tet-legs",df_iSNV_enumeration$sample_name)] <- "Tet-legs-14dpi_2"
df_iSNV_enumeration$ID[grepl("14-3-tet-legs",df_iSNV_enumeration$sample_name)] <- "Tet-legs-14dpi_3"
df_iSNV_enumeration$ID[grepl("4-1-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-4dpi_1"
df_iSNV_enumeration$ID[grepl("4-2-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-4dpi_2"
df_iSNV_enumeration$ID[grepl("4-3-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-4dpi_3"
df_iSNV_enumeration$ID[grepl("7-1-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-7dpi_1"
df_iSNV_enumeration$ID[grepl("7-2-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-7dpi_2"
df_iSNV_enumeration$ID[grepl("7-3-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-7dpi_3"
df_iSNV_enumeration$ID[grepl("dup-14-1-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-14dpi_1_dup"
df_iSNV_enumeration$ID[grepl("14-1-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-14dpi_1"
df_iSNV_enumeration$ID[grepl("dup-14-2-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-14dpi_2_dup"
df_iSNV_enumeration$ID[grepl("14-2-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-14dpi_2"
df_iSNV_enumeration$ID[grepl("dup-14-3-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-14dpi_3_dup"
df_iSNV_enumeration$ID[grepl("14-3-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-14dpi_3"
df_iSNV_enumeration$ID[grepl("7-1-wmel-body",df_iSNV_enumeration$sample_name)] <- "wmel-body-7dpi_1"
df_iSNV_enumeration$ID[grepl("7-2-wmel-body",df_iSNV_enumeration$sample_name)] <- "wmel-body-7dpi_2"
df_iSNV_enumeration$ID[grepl("7-3-wmel-body",df_iSNV_enumeration$sample_name)] <- "wmel-body-7dpi_3"
df_iSNV_enumeration$ID[grepl("14-1-wmel-body",df_iSNV_enumeration$sample_name)] <- "wmel-body-14dpi_1"
df_iSNV_enumeration$ID[grepl("14-2-wmel-body",df_iSNV_enumeration$sample_name)] <- "wmel-body-14dpi_2"
df_iSNV_enumeration$ID[grepl("14-3-wmel-body",df_iSNV_enumeration$sample_name)] <- "wmel-body-14dpi_3"
df_iSNV_enumeration$ID[grepl("7-1-wmel-legs",df_iSNV_enumeration$sample_name)] <- "wmel-legs-7dpi_1"
df_iSNV_enumeration$ID[grepl("14-2-wmel-legs",df_iSNV_enumeration$sample_name)] <- "wmel-legs-14dpi_2"
df_iSNV_enumeration$ID[grepl("7-2-tet-saliva",df_iSNV_enumeration$sample_name)] <- "Tet-saliva-7dpi_2"
df_iSNV_enumeration$ID[grepl("7-3-tet-saliva",df_iSNV_enumeration$sample_name)] <- "Tet-saliva-7dpi_3"
df_iSNV_enumeration$ID[grepl("14-1-tet-saliva",df_iSNV_enumeration$sample_name)] <- "Tet-saliva-14dpi_1"
df_iSNV_enumeration$ID[grepl("14-3-tet-saliva",df_iSNV_enumeration$sample_name)] <- "Tet-saliva-14dpi_3"
df_iSNV_enumeration$ID[grepl("7-1-tet-legs",df_iSNV_enumeration$sample_name)] <- "Tet-legs-7dpi_1"
df_iSNV_enumeration$ID[grepl("7-2-tet-legs",df_iSNV_enumeration$sample_name)] <- "Tet-legs-7dpi_2"
df_iSNV_enumeration$ID[grepl("7-3-tet-legs",df_iSNV_enumeration$sample_name)] <- "Tet-legs-7dpi_3"
df_iSNV_enumeration$ID[grepl("14-1-tet-legs",df_iSNV_enumeration$sample_name)] <- "Tet-legs-14dpi_1"
df_iSNV_enumeration$ID[grepl("14-2-tet-legs",df_iSNV_enumeration$sample_name)] <- "Tet-legs-14dpi_2"
df_iSNV_enumeration$ID[grepl("14-3-tet-legs",df_iSNV_enumeration$sample_name)] <- "Tet-legs-14dpi_3"
df_iSNV_enumeration$ID[grepl("4-1-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-4dpi_1"
df_iSNV_enumeration$ID[grepl("4-2-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-4dpi_2"
df_iSNV_enumeration$ID[grepl("4-3-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-4dpi_3"
df_iSNV_enumeration$ID[grepl("7-1-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-7dpi_1"
df_iSNV_enumeration$ID[grepl("7-2-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-7dpi_2"
df_iSNV_enumeration$ID[grepl("7-3-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-7dpi_3"
df_iSNV_enumeration$ID[grepl("dup-14-1-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-14dpi_1_dup"
df_iSNV_enumeration$ID[grepl("14-1-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-14dpi_1"
df_iSNV_enumeration$ID[grepl("dup-14-2-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-14dpi_2_dup"
df_iSNV_enumeration$ID[grepl("14-2-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-14dpi_2"
df_iSNV_enumeration$ID[grepl("dup-14-3-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-14dpi_3_dup"
df_iSNV_enumeration$ID[grepl("14-3-tet-body",df_iSNV_enumeration$sample_name)] <- "Tet-body-14dpi_3"
df_iSNV_enumeration$ID[grepl("7-1-wmel-body",df_iSNV_enumeration$sample_name)] <- "wmel-body-7dpi_1"
df_iSNV_enumeration$ID[grepl("7-2-wmel-body",df_iSNV_enumeration$sample_name)] <- "wmel-body-7dpi_2"
df_iSNV_enumeration$ID[grepl("7-3-wmel-body",df_iSNV_enumeration$sample_name)] <- "wmel-body-7dpi_3"
df_iSNV_enumeration$ID[grepl("14-1-wmel-body",df_iSNV_enumeration$sample_name)] <- "wmel-body-14dpi_1"
df_iSNV_enumeration$ID[grepl("14-2-wmel-body",df_iSNV_enumeration$sample_name)] <- "wmel-body-14dpi_2"
df_iSNV_enumeration$ID[grepl("14-3-wmel-body",df_iSNV_enumeration$sample_name)] <- "wmel-body-14dpi_3"
df_iSNV_enumeration$ID[grepl("7-1-wmel-legs",df_iSNV_enumeration$sample_name)] <- "wmel-legs-7dpi_1"
df_iSNV_enumeration$ID[grepl("14-2-wmel-legs",df_iSNV_enumeration$sample_name)] <- "wmel-legs-14dpi_2"

## remove controls
df_iSNV_enumeration_PC <- df_iSNV_enumeration[grepl("PC",df_iSNV_enumeration$sample_name),]
df_iSNV_enumeration_ZIKV <- df_iSNV_enumeration[grepl("ZIKV",df_iSNV_enumeration$sample_name),]
df_iSNV_enumeration_mouse <- df_iSNV_enumeration[grepl("mouse",df_iSNV_enumeration$sample_name),]
df_iSNV_enumeration <- df_iSNV_enumeration[!grepl("PC",df_iSNV_enumeration$sample_name),]
df_iSNV_enumeration <- df_iSNV_enumeration[!grepl("ZIKV",df_iSNV_enumeration$sample_name),]
df_iSNV_enumeration <- df_iSNV_enumeration[!grepl("mouse",df_iSNV_enumeration$sample_name),]

## separate cols
df_iSNV_enumeration <- separate(df_iSNV_enumeration, "ID", c("1","2","3"), sep="-")
df_iSNV_enumeration <- separate(df_iSNV_enumeration, "3", c("3","4"), sep="_")
df_iSNV_enumeration$`3` <- as.integer(gsub("dpi", "", df_iSNV_enumeration$`3`))

## bring controls back in
df_iSNV_enumeration_PC$`1` <- "Control"
df_iSNV_enumeration_mouse$`1` <- "Control"
df_iSNV_enumeration_PC$`2` <- "PC"
df_iSNV_enumeration_mouse$`2` <- "mouse"
df_iSNV_enumeration_PC$`3` <- "Control"
df_iSNV_enumeration_mouse$`3` <- "Control"
df_iSNV_enumeration_PC$`4` <- "Control"
df_iSNV_enumeration_mouse$`4` <- "Control"
df_iSNV_enumeration_PC$ID <- NULL
df_iSNV_enumeration_mouse$ID <- NULL

df_iSNV_enumeration <- rbind(df_iSNV_enumeration, 
                             df_iSNV_enumeration_PC,
                             df_iSNV_enumeration_mouse)

## groups
df_iSNV_enumeration$group_location <- paste(df_iSNV_enumeration$`1`, 
                                            df_iSNV_enumeration$`2`, sep = "_")
df_iSNV_enumeration$group_location_dpi <- paste(df_iSNV_enumeration$`1`, 
                                                df_iSNV_enumeration$`2`, 
                                                df_iSNV_enumeration$`3`, 
                                                sep = "_")
df_iSNV_enumeration$group_location <- as.factor(df_iSNV_enumeration$group_location)
df_iSNV_enumeration$group_location_dpi <- as.factor(df_iSNV_enumeration$group_location_dpi)

#### dir_10_snpdat_TSV #### not_done####


#### dir_10_snpdat_TXT #### not_done####


#### dir_12_snpgenie ####
#codon_results
setwd(paste(dir_09_reference_VCF, "sample_files", sep = "/")); dir()
sg_cr <- dir(pattern="_L001_fn_ann_noBC_sg_codon_results")
names_trunc <- gsub("_L001_fn_ann_noBC_sg_codon_results","",sg_cr)
names_trunc <- gsub("09-reference_","",names_trunc)
n <- length(sg_cr)
list <- vector("list",n)
for (i in 1:n) {
  list[[i]] <- as.data.frame(
    read_tsv(sg_cr[[i]]))
  names(list) <- names_trunc}
df_12_snpgenie_sg_cr <- Reduce(full_join,list)
df_12_snpgenie_sg_cr$sample <- df_12_snpgenie_sg_cr$file
df_12_snpgenie_sg_cr$sample <- gsub("./09-reference_", "", df_12_snpgenie_sg_cr$sample)
df_12_snpgenie_sg_cr$sample <- gsub("_L001_fn_ann_noBC.vcf", "", df_12_snpgenie_sg_cr$sample)
df_12_snpgenie_sg_cr <- separate(df_12_snpgenie_sg_cr, "sample", c("sample","S"), sep="_S")
#df_12_snpgenie_sg_cr$pi <- df_12_snpgenie_sg_cr$piN + df_12_snpgenie_sg_cr$piS
#df_12_snpgenie_sg_cr$piNpiS <- df_12_snpgenie_sg_cr$piN / df_12_snpgenie_sg_cr$piS
#df_12_snpgenie_sg_cr$piNminuspiS <- df_12_snpgenie_sg_cr$piN - df_12_snpgenie_sg_cr$piS

#product_results
sg_pr <- dir(pattern="_L001_fn_ann_noBC_sg_product_results.txt")
names_trunc <- gsub("_L001_fn_ann_noBC_sg_product_results.txt","",sg_pr)
names_trunc <- gsub("09-reference_","",names_trunc)
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
df_12_snpgenie_sg_pr$sample <- gsub("./09-reference_", "", df_12_snpgenie_sg_pr$sample)
df_12_snpgenie_sg_pr$sample <- gsub("_L001_fn_ann_noBC.vcf", "", df_12_snpgenie_sg_pr$sample)
df_12_snpgenie_sg_pr <- separate(df_12_snpgenie_sg_pr, "sample", c("sample","S"), sep="_S")
df_12_snpgenie_sg_pr$pi <- df_12_snpgenie_sg_pr$piN + df_12_snpgenie_sg_pr$piS
df_12_snpgenie_sg_pr$piNpiS <- df_12_snpgenie_sg_pr$piN / df_12_snpgenie_sg_pr$piS
df_12_snpgenie_sg_pr$piNminuspiS <- df_12_snpgenie_sg_pr$piN - df_12_snpgenie_sg_pr$piS

#sliding_window
sg_sw <- dir(pattern="_L001_fn_ann_noBC_sg_sliding_window.txt")
names_trunc <- gsub("_L001_fn_ann_noBC_sg_sliding_window.txt","",sg_sw)
names_trunc <- gsub("09-reference_","",names_trunc)
n <- length(sg_sw)
list <- vector("list",n)
for (i in 1:n) {
  list[[i]] <- as.data.frame(
    read_tsv(sg_sw[[i]]))
  names(list) <- names_trunc}
sg_sw <- Reduce(full_join,list)
sg_sw$pi <- sg_sw$piN + sg_sw$piS

sg_sw$sample <- sg_sw$file
sg_sw$sample <- gsub("./09-reference_", "", sg_sw$sample)
sg_sw$sample <- gsub("_L001_fn_ann_noBC.vcf", "", sg_sw$sample)
sg_sw <- separate(sg_sw, "sample", c("sample","S"), sep="_S")

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

## #### string matches snnpgenie ####
df_12_snpgenie_sg_pr$ID[grepl("7-2-tet-saliva",df_12_snpgenie_sg_pr$file)] <- "Tet-saliva-7dpi_2"
df_12_snpgenie_sg_pr$ID[grepl("7-3-tet-saliva",df_12_snpgenie_sg_pr$file)] <- "Tet-saliva-7dpi_3"
df_12_snpgenie_sg_pr$ID[grepl("14-1-tet-saliva",df_12_snpgenie_sg_pr$file)] <- "Tet-saliva-14dpi_1"
df_12_snpgenie_sg_pr$ID[grepl("14-3-tet-saliva",df_12_snpgenie_sg_pr$file)] <- "Tet-saliva-14dpi_3"
df_12_snpgenie_sg_pr$ID[grepl("7-1-tet-legs",df_12_snpgenie_sg_pr$file)] <- "Tet-legs-7dpi_1"
df_12_snpgenie_sg_pr$ID[grepl("7-2-tet-legs",df_12_snpgenie_sg_pr$file)] <- "Tet-legs-7dpi_2"
df_12_snpgenie_sg_pr$ID[grepl("7-3-tet-legs",df_12_snpgenie_sg_pr$file)] <- "Tet-legs-7dpi_3"
df_12_snpgenie_sg_pr$ID[grepl("14-1-tet-legs",df_12_snpgenie_sg_pr$file)] <- "Tet-legs-14dpi_1"
df_12_snpgenie_sg_pr$ID[grepl("14-2-tet-legs",df_12_snpgenie_sg_pr$file)] <- "Tet-legs-14dpi_2"
df_12_snpgenie_sg_pr$ID[grepl("14-3-tet-legs",df_12_snpgenie_sg_pr$file)] <- "Tet-legs-14dpi_3"
df_12_snpgenie_sg_pr$ID[grepl("4-1-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-4dpi_1"
df_12_snpgenie_sg_pr$ID[grepl("4-2-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-4dpi_2"
df_12_snpgenie_sg_pr$ID[grepl("4-3-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-4dpi_3"
df_12_snpgenie_sg_pr$ID[grepl("7-1-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-7dpi_1"
df_12_snpgenie_sg_pr$ID[grepl("7-2-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-7dpi_2"
df_12_snpgenie_sg_pr$ID[grepl("7-3-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-7dpi_3"
df_12_snpgenie_sg_pr$ID[grepl("dup-14-1-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-14dpi_1_dup"
df_12_snpgenie_sg_pr$ID[grepl("14-1-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-14dpi_1"
df_12_snpgenie_sg_pr$ID[grepl("dup-14-2-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-14dpi_2_dup"
df_12_snpgenie_sg_pr$ID[grepl("14-2-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-14dpi_2"
df_12_snpgenie_sg_pr$ID[grepl("dup-14-3-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-14dpi_3_dup"
df_12_snpgenie_sg_pr$ID[grepl("14-3-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-14dpi_3"
df_12_snpgenie_sg_pr$ID[grepl("7-1-wmel-body",df_12_snpgenie_sg_pr$file)] <- "wmel-body-7dpi_1"
df_12_snpgenie_sg_pr$ID[grepl("7-2-wmel-body",df_12_snpgenie_sg_pr$file)] <- "wmel-body-7dpi_2"
df_12_snpgenie_sg_pr$ID[grepl("7-3-wmel-body",df_12_snpgenie_sg_pr$file)] <- "wmel-body-7dpi_3"
df_12_snpgenie_sg_pr$ID[grepl("14-1-wmel-body",df_12_snpgenie_sg_pr$file)] <- "wmel-body-14dpi_1"
df_12_snpgenie_sg_pr$ID[grepl("14-2-wmel-body",df_12_snpgenie_sg_pr$file)] <- "wmel-body-14dpi_2"
df_12_snpgenie_sg_pr$ID[grepl("14-3-wmel-body",df_12_snpgenie_sg_pr$file)] <- "wmel-body-14dpi_3"
df_12_snpgenie_sg_pr$ID[grepl("7-1-wmel-legs",df_12_snpgenie_sg_pr$file)] <- "wmel-legs-7dpi_1"
df_12_snpgenie_sg_pr$ID[grepl("14-2-wmel-legs",df_12_snpgenie_sg_pr$file)] <- "wmel-legs-14dpi_2"
df_12_snpgenie_sg_pr$ID[grepl("7-2-tet-saliva",df_12_snpgenie_sg_pr$file)] <- "Tet-saliva-7dpi_2"
df_12_snpgenie_sg_pr$ID[grepl("7-3-tet-saliva",df_12_snpgenie_sg_pr$file)] <- "Tet-saliva-7dpi_3"
df_12_snpgenie_sg_pr$ID[grepl("14-1-tet-saliva",df_12_snpgenie_sg_pr$file)] <- "Tet-saliva-14dpi_1"
df_12_snpgenie_sg_pr$ID[grepl("14-3-tet-saliva",df_12_snpgenie_sg_pr$file)] <- "Tet-saliva-14dpi_3"
df_12_snpgenie_sg_pr$ID[grepl("7-1-tet-legs",df_12_snpgenie_sg_pr$file)] <- "Tet-legs-7dpi_1"
df_12_snpgenie_sg_pr$ID[grepl("7-2-tet-legs",df_12_snpgenie_sg_pr$file)] <- "Tet-legs-7dpi_2"
df_12_snpgenie_sg_pr$ID[grepl("7-3-tet-legs",df_12_snpgenie_sg_pr$file)] <- "Tet-legs-7dpi_3"
df_12_snpgenie_sg_pr$ID[grepl("14-1-tet-legs",df_12_snpgenie_sg_pr$file)] <- "Tet-legs-14dpi_1"
df_12_snpgenie_sg_pr$ID[grepl("14-2-tet-legs",df_12_snpgenie_sg_pr$file)] <- "Tet-legs-14dpi_2"
df_12_snpgenie_sg_pr$ID[grepl("14-3-tet-legs",df_12_snpgenie_sg_pr$file)] <- "Tet-legs-14dpi_3"
df_12_snpgenie_sg_pr$ID[grepl("4-1-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-4dpi_1"
df_12_snpgenie_sg_pr$ID[grepl("4-2-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-4dpi_2"
df_12_snpgenie_sg_pr$ID[grepl("4-3-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-4dpi_3"
df_12_snpgenie_sg_pr$ID[grepl("7-1-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-7dpi_1"
df_12_snpgenie_sg_pr$ID[grepl("7-2-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-7dpi_2"
df_12_snpgenie_sg_pr$ID[grepl("7-3-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-7dpi_3"
df_12_snpgenie_sg_pr$ID[grepl("dup-14-1-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-14dpi_1_dup"
df_12_snpgenie_sg_pr$ID[grepl("14-1-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-14dpi_1"
df_12_snpgenie_sg_pr$ID[grepl("dup-14-2-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-14dpi_2_dup"
df_12_snpgenie_sg_pr$ID[grepl("14-2-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-14dpi_2"
df_12_snpgenie_sg_pr$ID[grepl("dup-14-3-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-14dpi_3_dup"
df_12_snpgenie_sg_pr$ID[grepl("14-3-tet-body",df_12_snpgenie_sg_pr$file)] <- "Tet-body-14dpi_3"
df_12_snpgenie_sg_pr$ID[grepl("7-1-wmel-body",df_12_snpgenie_sg_pr$file)] <- "wmel-body-7dpi_1"
df_12_snpgenie_sg_pr$ID[grepl("7-2-wmel-body",df_12_snpgenie_sg_pr$file)] <- "wmel-body-7dpi_2"
df_12_snpgenie_sg_pr$ID[grepl("7-3-wmel-body",df_12_snpgenie_sg_pr$file)] <- "wmel-body-7dpi_3"
df_12_snpgenie_sg_pr$ID[grepl("14-1-wmel-body",df_12_snpgenie_sg_pr$file)] <- "wmel-body-14dpi_1"
df_12_snpgenie_sg_pr$ID[grepl("14-2-wmel-body",df_12_snpgenie_sg_pr$file)] <- "wmel-body-14dpi_2"
df_12_snpgenie_sg_pr$ID[grepl("14-3-wmel-body",df_12_snpgenie_sg_pr$file)] <- "wmel-body-14dpi_3"
df_12_snpgenie_sg_pr$ID[grepl("7-1-wmel-legs",df_12_snpgenie_sg_pr$file)] <- "wmel-legs-7dpi_1"
df_12_snpgenie_sg_pr$ID[grepl("14-2-wmel-legs",df_12_snpgenie_sg_pr$file)] <- "wmel-legs-14dpi_2"
## remove controls
df_12_snpgenie_sg_pr_NC <- df_12_snpgenie_sg_pr[grepl("NC",df_12_snpgenie_sg_pr$file),]
df_12_snpgenie_sg_pr_PC <- df_12_snpgenie_sg_pr[grepl("PC",df_12_snpgenie_sg_pr$file),]
df_12_snpgenie_sg_pr_ZIKV <- df_12_snpgenie_sg_pr[grepl("ZIKV",df_12_snpgenie_sg_pr$file),]
df_12_snpgenie_sg_pr_mouse <- df_12_snpgenie_sg_pr[grepl("mouse",df_12_snpgenie_sg_pr$file),]
df_12_snpgenie_sg_pr <- df_12_snpgenie_sg_pr[!grepl("NC",df_12_snpgenie_sg_pr$file),]
df_12_snpgenie_sg_pr <- df_12_snpgenie_sg_pr[!grepl("PC",df_12_snpgenie_sg_pr$file),]
df_12_snpgenie_sg_pr <- df_12_snpgenie_sg_pr[!grepl("ZIKV",df_12_snpgenie_sg_pr$file),]
df_12_snpgenie_sg_pr <- df_12_snpgenie_sg_pr[!grepl("mouse",df_12_snpgenie_sg_pr$file),]
## separate cols
df_12_snpgenie_sg_pr <- separate(df_12_snpgenie_sg_pr, "ID", c("1","2","3"), sep="-")
df_12_snpgenie_sg_pr <- separate(df_12_snpgenie_sg_pr, "3", c("3","4"), sep="_")
df_12_snpgenie_sg_pr$`3` <- as.integer(gsub("dpi", "", df_12_snpgenie_sg_pr$`3`))

## bring controls back in
df_12_snpgenie_sg_pr_NC$`1` <- "Control"
df_12_snpgenie_sg_pr_PC$`1` <- "Control"
df_12_snpgenie_sg_pr_mouse$`1` <- "Control"
df_12_snpgenie_sg_pr_NC$`2` <- "NC"
df_12_snpgenie_sg_pr_PC$`2` <- "PC"
df_12_snpgenie_sg_pr_mouse$`2` <- "mouse"
df_12_snpgenie_sg_pr_NC$`3` <- "Control"
df_12_snpgenie_sg_pr_PC$`3` <- "Control"
df_12_snpgenie_sg_pr_mouse$`3` <- "Control"
df_12_snpgenie_sg_pr_NC$`4` <- "Control"
df_12_snpgenie_sg_pr_PC$`4` <- "Control"
df_12_snpgenie_sg_pr_mouse$`4` <- "Control"
df_12_snpgenie_sg_pr_NC$ID <- NULL
df_12_snpgenie_sg_pr_PC$ID <- NULL
df_12_snpgenie_sg_pr_mouse$ID <- NULL

df_12_snpgenie_sg_pr <- rbind(df_12_snpgenie_sg_pr,
                              df_12_snpgenie_sg_pr_NC,
                              df_12_snpgenie_sg_pr_PC,
                              df_12_snpgenie_sg_pr_mouse)

## groups
df_12_snpgenie_sg_pr$group_location <- paste(df_12_snpgenie_sg_pr$`1`, 
                                            df_12_snpgenie_sg_pr$`2`, sep = "_")
df_12_snpgenie_sg_pr$group_location_dpi <- paste(df_12_snpgenie_sg_pr$`1`, 
                                                df_12_snpgenie_sg_pr$`2`, 
                                                df_12_snpgenie_sg_pr$`3`, 
                                                sep = "_")
df_12_snpgenie_sg_pr$group_location <- as.factor(df_12_snpgenie_sg_pr$group_location)
df_12_snpgenie_sg_pr$group_location_dpi <- as.factor(df_12_snpgenie_sg_pr$group_location_dpi)


## melt
molten_sg_pr_pi <- melt(df_12_snpgenie_sg_pr,
                        id.vars = c("group_location_dpi"),
                        measure.vars = c("pi"),
                        variable.name = "Pi")
molten_sg_pr_piNpiS <- melt(df_12_snpgenie_sg_pr,
                            id.vars = c("group_location_dpi"),
                            measure.vars = c("piNpiS"),
                            variable.name = "Pi")
molten_sg_pr_piN <- melt(df_12_snpgenie_sg_pr,
                         id.vars = c("group_location_dpi"),
                         measure.vars = c("piN"),
                         variable.name = "Pi")
molten_sg_pr_piS <- melt(df_12_snpgenie_sg_pr,
                         id.vars = c("group_location_dpi"),
                         measure.vars = c("piS"),
                         variable.name = "Pi")
molten_sg_pr_piNminuspiS <- melt(df_12_snpgenie_sg_pr,
                         id.vars = c("group_location_dpi"),
                         measure.vars = c("piNminuspiS"),
                         variable.name = "Pi")
molten_sg_pr <- rbind(molten_sg_pr_pi, molten_sg_pr_piNpiS, molten_sg_pr_piN, molten_sg_pr_piS, molten_sg_pr_piNminuspiS)
molten_sg_pr <- molten_sg_pr[!is.na(molten_sg_pr$value),]

molten_sg_pr <- molten_sg_pr[molten_sg_pr$group_location_dpi!="Control_NC_Control",]
molten_sg_pr$group_location_dpi <- as.factor(as.character(molten_sg_pr$group_location_dpi))
levels <- levels(molten_sg_pr$group_location_dpi)

i=1
temp <- molten_sg_pr[molten_sg_pr$group_location_dpi==levels[i],]
temp_pi <- ds(data.frame("b" = b(temp$value[temp$Pi=="pi"]), "g" = "pi"), varname = "V1", groupnames = "g")
temp_piNpiS <- ds(data.frame("b" = b(temp$value[temp$Pi=="piNpiS"]), "g" = "piNpiS"), varname = "V1", groupnames = "g")
temp_piN <- ds(data.frame("b" = b(temp$value[temp$Pi=="piN"]), "g" = "piN"), varname = "V1", groupnames = "g")
temp_piS <- ds(data.frame("b" = b(temp$value[temp$Pi=="piS"]), "g" = "piS"), varname = "V1", groupnames = "g")
temp_piNminuspiS <- ds(data.frame("b" = b(temp$value[temp$Pi=="piNminuspiS"]), "g" = "piNminuspiS"), varname = "V1", groupnames = "g")
ds_b_sg_pr <- rbind(temp_pi, temp_piNpiS, temp_piN, temp_piS, temp_piNminuspiS)
ds_b_sg_pr$gld <- levels[i]
for (i in 2:length(levels)) {
  print(levels[i])
  temp <- molten_sg_pr[molten_sg_pr$group_location_dpi==levels[i],]
  temp_pi <- ds(data.frame("b" = b(temp$value[temp$Pi=="pi"]), "g" = "pi"), varname = "V1", groupnames = "g")
  temp_piNpiS <- ds(data.frame("b" = b(temp$value[temp$Pi=="piNpiS"]), "g" = "piNpiS"), varname = "V1", groupnames = "g")
  temp_piN <- ds(data.frame("b" = b(temp$value[temp$Pi=="piN"]), "g" = "piN"), varname = "V1", groupnames = "g")
  temp_piS <- ds(data.frame("b" = b(temp$value[temp$Pi=="piS"]), "g" = "piS"), varname = "V1", groupnames = "g")
  temp_piNminuspiS <- ds(data.frame("b" = b(temp$value[temp$Pi=="piNminuspiS"]), "g" = "piNminuspiS"), varname = "V1", groupnames = "g")
  temp2 <- rbind(temp_pi, temp_piNpiS, temp_piN, temp_piS, temp_piNminuspiS)
  temp2$gld <- levels[i]
  ds_b_sg_pr <- rbind(ds_b_sg_pr, temp2)
}
ds_b_sg_pr <- separate(ds_b_sg_pr, "gld", c("group", "location", "dpi"), sep = "_")
ds_b_sg_pr_control <- ds_b_sg_pr[ds_b_sg_pr$group=="Control",]
ds_b_sg_pr_control$gl <- paste(ds_b_sg_pr_control$group, ds_b_sg_pr_control$location)
ds_b_sg_pr <- ds_b_sg_pr[ds_b_sg_pr$group!="Control",]
ds_b_sg_pr$gl <- factor(paste(ds_b_sg_pr$group, ds_b_sg_pr$location, sep = "_"), 
                        levels = c("Tet_saliva", "Tet_body", "wmel_body", "Tet_legs", "wmel_legs"))
ds_b_sg_pr$gld <- as.factor(paste(ds_b_sg_pr$group, ds_b_sg_pr$location, ds_b_sg_pr$dpi, sep = "_"))
ds_b_sg_pr$dpi <- factor(ds_b_sg_pr$dpi, levels = c("4", "7", "14"))
ds_b_sg_pr$g <- factor(ds_b_sg_pr$g, levels = c("pi", "piN", "piS", "piNpiS", "piNminuspiS"))

#### Plot variants ####
## Variant enumeration
df_iSNV_enumeration$gl <- factor(paste(df_iSNV_enumeration$`1`, df_iSNV_enumeration$`2`, sep = "_"), 
                                 levels = c("Tet_body", "Tet_legs", "Tet_saliva", "wmel_body", "wmel_legs"))
table(df_iSNV_enumeration$gl)
df_iSNV_enumeration$dpi <- factor(df_iSNV_enumeration$`3`, 
                                 levels = c(4, 7, 14))
table(df_iSNV_enumeration$dpi)

## ds_iSNV_enumeration_Total
ds_iSNV_enumeration_Total <- as.data.frame(ds(df_iSNV_enumeration, 
                                        varname = "Total", 
                                        groupnames = c("group_location_dpi")))
ds_iSNV_enumeration_Total <- separate(ds_iSNV_enumeration_Total, "group_location_dpi", 
                                    into = c("group", "location", "dpi"), sep = "_")
ds_iSNV_enumeration_Total$group_location <- paste(ds_iSNV_enumeration_Total$group, 
                                                ds_iSNV_enumeration_Total$location, sep = "_")
ds_iSNV_enumeration_Total$group_location_dpi <- paste(ds_iSNV_enumeration_Total$group, 
                                                    ds_iSNV_enumeration_Total$location, 
                                                    ds_iSNV_enumeration_Total$dpi, 
                                                    sep = "_")
ds_iSNV_enumeration_Total$group_location_dpi <- as.factor(ds_iSNV_enumeration_Total$group_location_dpi)
ds_iSNV_enumeration_Total$group_location <- as.factor(ds_iSNV_enumeration_Total$group_location)
ds_iSNV_enumeration_Total$group_location <- factor(ds_iSNV_enumeration_Total$group_location,
                                                   levels = c("Tet_body", "Tet_legs", "Tet_saliva", "wmel_body", "wmel_legs", "Control_mouse", "Control_PC"))
ds_iSNV_enumeration_Total <- separate(ds_iSNV_enumeration_Total, "group_location_dpi", c("g", "l", "dpi"))
ds_iSNV_enumeration_Total$dpi <- factor(ds_iSNV_enumeration_Total$dpi, 
                                      levels = c("4", "7", "14"))


## ds_iSNV_enumeration_TAF
ds_iSNV_enumeration_TAF <- as.data.frame(ds(df_iSNV_enumeration, 
                                              varname = "Total_AF", 
                                              groupnames = c("group_location_dpi")))
ds_iSNV_enumeration_TAF <- separate(ds_iSNV_enumeration_TAF, "group_location_dpi", 
                                      into = c("group", "location", "dpi"), sep = "_")
ds_iSNV_enumeration_TAF$group_location <- paste(ds_iSNV_enumeration_TAF$group, 
                                                  ds_iSNV_enumeration_TAF$location, sep = "_")
ds_iSNV_enumeration_TAF$group_location_dpi <- paste(ds_iSNV_enumeration_TAF$group, 
                                                      ds_iSNV_enumeration_TAF$location, 
                                                      ds_iSNV_enumeration_TAF$dpi, 
                                                      sep = "_")
ds_iSNV_enumeration_TAF$group_location_dpi <- as.factor(ds_iSNV_enumeration_TAF$group_location_dpi)
ds_iSNV_enumeration_TAF$group_location <- factor(ds_iSNV_enumeration_TAF$group_location,
                                                 levels = c("Tet_body", "Tet_legs", "Tet_saliva", "wmel_body", "wmel_legs", "Control_mouse", "Control_PC"))
ds_iSNV_enumeration_TAF <- separate(ds_iSNV_enumeration_TAF, "group_location_dpi", c("g", "l", "dpi"))
ds_iSNV_enumeration_TAF$dpi <- factor(ds_iSNV_enumeration_TAF$dpi, 
                                        levels = c("4", "7", "14"))

## Total 
plot_total_iSNVs <- ggplot(ds_iSNV_enumeration_Total, aes(color = dpi), group = dpi) + 
  geom_point(aes(x = group_location, y = Total),
             position = position_dodge(.75)) + 
  geom_errorbar(aes(x = group_location, y = Total, 
                    ymin = Total - sd, 
                    ymax = Total + sd), 
                width = .15, position = position_dodge(.75)) +
  theme(#legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.position = "right") +
  scale_y_continuous(limits = c(0,30)) + 
  labs(x = "", y = "Number of iSNVs") + 
  axis_formatting + legend_formatting + background_formatting
## TAF
plot_divergence <- ggplot(ds_iSNV_enumeration_TAF, aes(color = dpi), group = dpi) + 
  geom_point(aes(x = group_location, y = Total_AF),
             position = position_dodge(.75)) + 
  geom_errorbar(aes(x = group_location, y = Total_AF, 
                    ymin = Total_AF - sd, 
                    ymax = Total_AF + sd), 
                width = .15, position = position_dodge(.75)) +
  theme(#legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.position = "right") +
  scale_y_continuous(limits = c(0,6.5)) + 
  labs(x = "", y = "Divergence") + 
  axis_formatting + legend_formatting + background_formatting
plot_divergence_counts <- plot_grid(plot_total_iSNVs, plot_divergence, nrow = 2)

# Total iSNVs
plot_enumeration_Total_1 <- ggplot(df_iSNV_enumeration, aes(color = gl), group = dpi) + 
  geom_violin(aes(x = dpi, y = Total),
              position = position_dodge(.75)) + 
  scale_color_manual(values = palette_gl) + 
  theme(legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.position = "right") +
  labs(x = "", y = "Number of iSNVs") + 
  axis_formatting + legend_formatting + background_formatting
plot_enumeration_Total_2 <- ggplot(df_iSNV_enumeration, aes(color = dpi), group = dpi) + 
  geom_violin(aes(x = gl, y = Total),
              position = position_dodge(.75)) + 
  theme(legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.position = "right") +
  labs(x = "", y = "Number of iSNVs") + 
  axis_formatting + legend_formatting + background_formatting
# Total AF
plot_enumeration_TAF_1 <- ggplot(df_iSNV_enumeration, aes(color = gl), group = dpi) + 
  geom_violin(aes(x = dpi, y = Total_AF),
              position = position_dodge(.75)) + 
  scale_color_manual(values = palette_gl) + 
  theme(legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.position = "right") +
  labs(x = "", y = "Divergence") + 
  axis_formatting + legend_formatting + background_formatting
plot_enumeration_TAF_2 <- ggplot(df_iSNV_enumeration, aes(color = dpi), group = dpi) + 
  geom_violin(aes(x = gl, y = Total_AF),
              position = position_dodge(.75)) + 
  theme(legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.position = "right") +
  labs(x = "", y = "Divergence") + 
  axis_formatting + legend_formatting + background_formatting

plots_enumeration_1 <- plot_grid(plot_enumeration_TAF_1, plot_enumeration_Total_1,
                                 nrow = 2, ncol = 1, rel_widths = c(1,1.5))
plots_enumeration_2 <- plot_grid(plot_enumeration_TAF_2, plot_enumeration_Total_2,
                                 nrow = 2, ncol = 1, rel_widths = c(1,1.5))



## Variants across the genome
df_09_reference_VCF_control <- df_09_reference_VCF[df_09_reference_VCF$`3`=="Control",]
df_09_reference_VCF_control$gl <- factor(paste(df_09_reference_VCF_control$`1`, df_09_reference_VCF_control$`2`, sep = "_"), 
                                         levels = c("Control_mouse", "Control_PC"))
df_09_reference_VCF <- df_09_reference_VCF[df_09_reference_VCF$`3`!="Control",]
df_09_reference_VCF$gl <- factor(paste(df_09_reference_VCF$`1`, df_09_reference_VCF$`2`, sep = "_"), 
                                 levels = c("Tet_body", "wmel_body", "Tet_legs", "wmel_legs", "Tet_saliva"))
#table(df_09_reference_VCF$gl)
#n_samples <- as.data.frame(df_09_reference_VCF
#                       [!duplicated(df_09_reference_VCF$FILTER),] 
#                       %>% count(gl, `3`))
#n_samples <- table[order(table[,"3"]),]
#n_samples <- separate(n_samples, "gl", c("g", "l"), sep = "_")
#n_samples$d <- n_samples$`3`
#n_samples$gl <- factor(paste(n_samples$g, n_samples$l, sep = "_"), 
#                       levels = c("Tet_saliva", "Tet_body", "Tet_legs", "wmel_saliva", "wmel_body",  "wmel_legs"))
#
#plot_samples <- ggplot(data = n_samples, aes(x = d, y = n, group = gl, color = gl)) + 
#  geom_point() + 
#  geom_line() + 
#  theme(legend.title = element_blank(), 
#        legend.key = element_blank(), 
#        legend.position = "right") +
#  labs(x = "days post-infection", y = "number of samples") +
#  scale_x_continuous(breaks = c(0, 4, 7, 14), limits = c(0, 14)) + 
#  scale_color_manual(values = palette_gl) + 
#  #facet_grid(rows = vars(gl), cols = vars(`3`)) + 
#  axis_formatting + legend_formatting + background_formatting

df_09_reference_VCF$`3` <- factor(df_09_reference_VCF$`3`, levels = c("4", "7", "14"))
plot_iSNVs <- ggplot(data = df_09_reference_VCF[df_09_reference_VCF$HGVS.p!="p.Pro1683Pro" & 
                                                  df_09_reference_VCF$POS > 4030 | df_09_reference_VCF$POS < 4007,], 
                     aes(x = POS, y = AF, color = Ann)) + 
  geom_point() + 
  labs(x = "", y = "Allele Frequency") +
  scale_y_continuous(n.breaks = 2) + 
  scale_x_continuous(n.breaks = 5) + 
  theme(legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.position = "none") +
  scale_color_manual(values = palette_muts) + 
  facet_grid(rows = vars(gl), cols = vars(`3`)) + 
  axis_formatting + legend_formatting + background_formatting
plot_iSNVs_control <- ggplot(data = df_09_reference_VCF_control[df_09_reference_VCF_control$HGVS.p!="p.Pro1683Pro" & 
                                                                  df_09_reference_VCF_control$POS > 4030 | df_09_reference_VCF_control$POS < 4007,], 
                     aes(x = POS, y = AF, color = Ann)) + 
  geom_point() + 
  labs(x = "Genome position", y = "Allele Frequency") +
  scale_y_continuous(n.breaks = 2) + 
  scale_x_continuous(n.breaks = 5) + 
  theme(legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.position = "bottom") +
  scale_color_manual(values = palette_muts) + 
  facet_grid( cols = vars(`2`)) + 
  axis_formatting + legend_formatting + background_formatting
plot_iSNVs_all <- plot_grid(plot_iSNVs, plot_iSNVs_control, nrow = 2,
                            rel_heights = c(1,.5))

## iSNV frequency spectra
# assign VCF to df_sub50AF
df_sub50AF <- df_09_reference_VCF[df_09_reference_VCF$AF <= .5,]
# add AF bins
df_sub50AF$Bin[df_sub50AF$AF>0.0 & df_sub50AF$AF<=0.1] <- "1-10%"
df_sub50AF$Bin[df_sub50AF$AF>0.1 & df_sub50AF$AF<=0.2] <- "10-20%"
df_sub50AF$Bin[df_sub50AF$AF>0.2 & df_sub50AF$AF<=0.3] <- "20-30%"
df_sub50AF$Bin[df_sub50AF$AF>0.3 & df_sub50AF$AF<=0.4] <- "30-40%"
df_sub50AF$Bin[df_sub50AF$AF>0.4 & df_sub50AF$AF<=0.5] <- "40-50%"
df_sub50AF <- df_sub50AF[df_sub50AF$Ann!="stop_lost&splice_region_variant",]

g <- levels(df_sub50AF$group_location_dpi)
df_sub50AF_backup <- df_sub50AF


## start of 
i_gld = "Tet_body_4"
df_sub50AF <- df_sub50AF_backup[df_sub50AF_backup$group_location_dpi==i_gld,]
df_sub50AF_md_n_SNVs_by_bin_and_mut <- as.data.frame(table(df_sub50AF$Bin,df_sub50AF$Ann))
names(df_sub50AF_md_n_SNVs_by_bin_and_mut)[names(df_sub50AF_md_n_SNVs_by_bin_and_mut) == "Var1"] <- "Bins"
names(df_sub50AF_md_n_SNVs_by_bin_and_mut)[names(df_sub50AF_md_n_SNVs_by_bin_and_mut) == "Var2"] <- "Mutation_type"
names(df_sub50AF_md_n_SNVs_by_bin_and_mut)[names(df_sub50AF_md_n_SNVs_by_bin_and_mut) == "Freq"] <- "Count"
# split df by Mutation_type
list_mut_bins_prop <- split(df_sub50AF_md_n_SNVs_by_bin_and_mut, with(df_sub50AF_md_n_SNVs_by_bin_and_mut, Mutation_type, drop=T))
## if missing a bin, add it in
n <- length(list_mut_bins_prop)
for (i in 1:n) {
  if (length(list_mut_bins_prop[[i]]$Bins)<5) {
    print("Missing a bin!")
    bins <- c("1-10%", "10-20%", "20-30%", "30-40%", "40-50%")
    missing_bin <- bins[which(bins %!in% list_mut_bins_prop[[i]]$Bins[1:length(list_mut_bins_prop[[i]]$Bins)])]
    missing_row <- data.frame("Bins" = missing_bin, "Mutation_type" = list_mut_bins_prop[[i]]$Mutation_type[1], "Count" = 0)
    list_mut_bins_prop[[i]] <- rbind(list_mut_bins_prop[[i]], missing_row)
    list_mut_bins_prop[[i]] <- list_mut_bins_prop[[i]] %>% slice(match(bins, Bins))
  }
}
# calculate proportional column
n <- length(list_mut_bins_prop)
for (i in 1:n) {
  list_mut_bins_prop[[i]]$Prop[1] <- sum(list_mut_bins_prop[[i]][1,3]/sum(list_mut_bins_prop[[i]][,3]))
  list_mut_bins_prop[[i]]$Prop[2] <- sum(list_mut_bins_prop[[i]][2,3]/sum(list_mut_bins_prop[[i]][,3]))
  list_mut_bins_prop[[i]]$Prop[3] <- sum(list_mut_bins_prop[[i]][3,3]/sum(list_mut_bins_prop[[i]][,3]))
  list_mut_bins_prop[[i]]$Prop[4] <- sum(list_mut_bins_prop[[i]][4,3]/sum(list_mut_bins_prop[[i]][,3]))
  list_mut_bins_prop[[i]]$Prop[5] <- sum(list_mut_bins_prop[[i]][5,3]/sum(list_mut_bins_prop[[i]][,3]))
}
# list to df
df_mut_bins_prop_all <- Reduce(full_join,list_mut_bins_prop)
## only include some of the mutation types
df_mut_bins_prop_S <- filter(df_mut_bins_prop_all, Mutation_type=="Synonymous")
df_mut_bins_prop_NS <- filter(df_mut_bins_prop_all, Mutation_type=="Nonsynonymous")
df_mut_bins_prop_Stop_g <- filter(df_mut_bins_prop_all, Mutation_type=="Stop gained")
df_mut_bins_prop_noNeut <- rbind(df_mut_bins_prop_S, df_mut_bins_prop_NS, df_mut_bins_prop_Stop_g)
df_mut_bins_prop_noNeut$gld <- i_gld
temp <- df_mut_bins_prop_noNeut

x <- function(i_gld) {
  df_sub50AF <- df_sub50AF_backup[df_sub50AF_backup$group_location_dpi==i_gld,]
  df_sub50AF_md_n_SNVs_by_bin_and_mut <- as.data.frame(table(df_sub50AF$Bin,df_sub50AF$Ann))
  names(df_sub50AF_md_n_SNVs_by_bin_and_mut)[names(df_sub50AF_md_n_SNVs_by_bin_and_mut) == "Var1"] <- "Bins"
  names(df_sub50AF_md_n_SNVs_by_bin_and_mut)[names(df_sub50AF_md_n_SNVs_by_bin_and_mut) == "Var2"] <- "Mutation_type"
  names(df_sub50AF_md_n_SNVs_by_bin_and_mut)[names(df_sub50AF_md_n_SNVs_by_bin_and_mut) == "Freq"] <- "Count"
  # split df by Mutation_type
  list_mut_bins_prop <- split(df_sub50AF_md_n_SNVs_by_bin_and_mut, with(df_sub50AF_md_n_SNVs_by_bin_and_mut, Mutation_type, drop=T))
  ## if missing a bin, add it in
  n <- length(list_mut_bins_prop)
  for (i in 1:n) {
    if (length(list_mut_bins_prop[[i]]$Bins)<5) {
      print("Missing a bin!")
      bins <- c("1-10%", "10-20%", "20-30%", "30-40%", "40-50%")
      missing_bin <- bins[which(bins %!in% list_mut_bins_prop[[i]]$Bins[1:length(list_mut_bins_prop[[i]]$Bins)])]
      missing_row <- data.frame("Bins" = missing_bin, "Mutation_type" = list_mut_bins_prop[[i]]$Mutation_type[1], "Count" = 0)
      list_mut_bins_prop[[i]] <- rbind(list_mut_bins_prop[[i]], missing_row)
      list_mut_bins_prop[[i]] <- list_mut_bins_prop[[i]] %>% slice(match(bins, Bins))
    }
  }
  # calculate proportional column
  n <- length(list_mut_bins_prop)
  for (i in 1:n) {
    list_mut_bins_prop[[i]]$Prop[1] <- sum(list_mut_bins_prop[[i]][1,3]/sum(list_mut_bins_prop[[i]][,3]))
    list_mut_bins_prop[[i]]$Prop[2] <- sum(list_mut_bins_prop[[i]][2,3]/sum(list_mut_bins_prop[[i]][,3]))
    list_mut_bins_prop[[i]]$Prop[3] <- sum(list_mut_bins_prop[[i]][3,3]/sum(list_mut_bins_prop[[i]][,3]))
    list_mut_bins_prop[[i]]$Prop[4] <- sum(list_mut_bins_prop[[i]][4,3]/sum(list_mut_bins_prop[[i]][,3]))
    list_mut_bins_prop[[i]]$Prop[5] <- sum(list_mut_bins_prop[[i]][5,3]/sum(list_mut_bins_prop[[i]][,3]))
  }
  # list to df
  df_mut_bins_prop_all <- Reduce(full_join,list_mut_bins_prop)
  ## only include some of the mutation types
  df_mut_bins_prop_S <- filter(df_mut_bins_prop_all, Mutation_type=="Synonymous")
  df_mut_bins_prop_NS <- filter(df_mut_bins_prop_all, Mutation_type=="Nonsynonymous")
  df_mut_bins_prop_Stop_g <- filter(df_mut_bins_prop_all, Mutation_type=="Stop gained")
  df_mut_bins_prop_noNeut <- rbind(df_mut_bins_prop_S, df_mut_bins_prop_NS, df_mut_bins_prop_Stop_g)
  df_mut_bins_prop_noNeut$gld <- i_gld
  temp <- rbind(temp, df_mut_bins_prop_noNeut)
  return(temp)
}
temp <- x("Tet_body_7")
temp <- x("Tet_body_14")
temp <- x("Tet_legs_14")
temp <- x("Tet_legs_7")
temp <- x("Tet_saliva_14")
temp <- x("Tet_saliva_7")
temp <- x("wmel_body_14")
temp <- x("wmel_body_7")
temp <- x("wmel_legs_14")
temp <- x("wmel_legs_7")
#temp <- x("Control_PC_Control")
#temp <- x("Control_mouse_Control")
df_mut_bins_prop_noNeut <- temp

## df for neutral expectation
df_mut_bins_neutral <- data.frame("Bins" = c("1-10%", "10-20%", "20-30%", "30-40%", "40-50%"),
                                  "Mutation_type" = c("Neutral expectation","Neutral expectation",
                                                      "Neutral expectation","Neutral expectation",
                                                      "Neutral expectation"),
                                  "Count" = c(NA, NA, NA, NA, NA),
                                  "Prop" = c(0.588592, 0.177184, 0.103646, 0.073538, 0.057040))

## iSNV frequency spectrum
# factor leving
df_mut_bins_prop_noNeut$Mutation_type <- factor(df_mut_bins_prop_noNeut$Mutation_type, levels = c("Nonsynonymous", "Synonymous", "Stop gained"))
df_mut_bins_prop_noNeut <- df_mut_bins_prop_noNeut[df_mut_bins_prop_noNeut$Mutation_type!="Stop gained",]
df_mut_bins_prop_noNeut <- separate(df_mut_bins_prop_noNeut, "gld", c("group", "location", "dpi"), sep = "_")
df_mut_bins_prop_noNeut$gl <- factor(paste(df_mut_bins_prop_noNeut$group, df_mut_bins_prop_noNeut$location, sep = "_"), 
                                     levels = c("Tet_body", "wmel_body", "Tet_legs", "wmel_legs", "Tet_saliva"))
df_mut_bins_prop_noNeut$gld <- paste(df_mut_bins_prop_noNeut$gl, df_mut_bins_prop_noNeut$dpi, sep = "_")
df_mut_bins_prop_noNeut$d <- factor(paste(df_mut_bins_prop_noNeut$dpi), 
                                    levels = c("4", "7", "14"))
df_mut_bins_prop_noNeut$Mutation_type <- factor(paste(df_mut_bins_prop_noNeut$Mutation_type), 
                                    levels = c("Nonsynonymous", "Synonymous"))

plot_spec <- ggplot() + 
  geom_bar(data = df_mut_bins_prop_noNeut, aes(x = Bins, y = Prop, fill = Mutation_type),
           stat = "identity", position = position_dodge2()) +
  geom_text(data = df_mut_bins_prop_noNeut, 
            aes(x = Bins, y = Prop, group = Mutation_type, label = Count),
            position=position_dodge2(0.9), vjust=-0.25, size = 1.5) + 
  geom_point(data = df_mut_bins_neutral, aes(x = Bins, y = Prop)) + 
  geom_line(data = df_mut_bins_neutral, aes(x = Bins, y = Prop, group = 1)) +
  scale_fill_manual(values = palette_muts) + 
  labs(x = "Within-host iSNV frequency bin", 
       y = "Proportion of variants per mutation type") + 
  scale_y_continuous(limits = c(0,1), n.breaks = 2) + 
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "bottom") + 
  facet_grid(rows = vars(gl), cols = vars(d)) + 
  axis_formatting + 
  legend_formatting + 
  background_formatting


#### Plot snpgenie ####
#piN and piS
w <- c(1, .1, .8, .1)
ds_b_sg_pr$gl <- factor(paste(ds_b_sg_pr$group, ds_b_sg_pr$location, sep = "_"), 
                        levels = c("Tet_body", "Tet_legs", "Tet_saliva", "wmel_body", "wmel_legs"))
plot3 <- ggplot(data = ds_b_sg_pr[ds_b_sg_pr$g=="piN" | ds_b_sg_pr$g =="piS",], 
                aes(x = gl, y = V1, group = gld, color = g)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = gl, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  scale_color_manual(values = c("piN" = "#FF7F20",
                                "piS" = "#4F7899")) +
  labs(y = "Nucleotide diversity") + 
  scale_y_continuous(limits = c(0, 0.0015)) +
  facet_grid(cols = vars(dpi)) +
  theme(legend.title = element_blank(), 
        axis.title.x.bottom = element_blank(),
        legend.key = element_blank(), 
        legend.position = "none") +  
  axis_formatting + legend_formatting + background_formatting
plot3.1 <- ggplot(data = ds_b_sg_pr_control[ds_b_sg_pr_control$g=="piN" | ds_b_sg_pr_control$g =="piS",], 
                  aes(x = location, y = V1, group = location, color = g)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = location, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  scale_color_manual(values = c("piN" = "#FF7F20",
                                "piS" = "#4F7899")) +
  labs(y = "Nucleotide diversity") + 
  scale_y_continuous(limits = c(0, 0.0015)) +
  facet_grid(cols = vars(dpi)) +
  theme(legend.title = element_blank(), 
        axis.title.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        legend.key = element_blank(), 
        legend.position = "none") +  
  axis_formatting + legend_formatting + background_formatting
plot4 <- ggplot(data = ds_b_sg_pr[ds_b_sg_pr$g=="piN" | ds_b_sg_pr$g =="piS",], 
                aes(x = dpi, y = V1, group = gld, color = g)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = dpi, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  scale_color_manual(values = c("piN" = "#FF7F20",
                                "piS" = "#4F7899")) +
  labs(y = "Nucleotide diversity") + 
  scale_y_continuous(limits = c(0, 0.0015)) +
  facet_grid(cols = vars(gl)) +
  theme(legend.title = element_blank(), 
        axis.title.x.bottom = element_blank(),
        legend.key = element_blank(), 
        legend.position = "none") +  
  axis_formatting + legend_formatting + background_formatting
plot4.1 <- ggplot(data = ds_b_sg_pr_control[ds_b_sg_pr_control$g=="piN" | ds_b_sg_pr_control$g =="piS",], 
                  aes(x = location, y = V1, group = location, color = g)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = location, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  scale_color_manual(values = c("piN" = "#FF7F20",
                                "piS" = "#4F7899")) +
  labs(y = "Nucleotide diversity") + 
  scale_y_continuous(limits = c(0, 0.0015)) +
  facet_grid(cols = vars(group)) +
  theme(legend.title = element_blank(), 
        axis.title.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        legend.key = element_blank(), 
        legend.position = "none") +  
  axis_formatting + legend_formatting + background_formatting
plot3_4 <- plot_grid(plot3, plot3.1, plot4, plot4.1, rel_widths = w, 
                     hjust = LR_adjust, vjust = UD_adjust, nrow = 1)

# pi
ds_b_sg_pr$gl <- factor(paste(ds_b_sg_pr$group, ds_b_sg_pr$location, sep = "_"), 
                        levels = c("Tet_body", "Tet_legs", "Tet_saliva", "wmel_body", "wmel_legs"))
plot5 <- ggplot(data = ds_b_sg_pr[ds_b_sg_pr$g=="pi",], 
                aes(x = gl, y = V1, group = gld, color = dpi)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = gl, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  facet_grid(cols = vars(dpi)) +
  labs(y = "Nucleotide diversity") + 
  scale_y_continuous(limits = c(0, 0.0015)) +
  theme(legend.title = element_blank(), 
        axis.title.x.bottom = element_blank(),
        legend.key = element_blank(), 
        legend.position = "none") +  
  axis_formatting + legend_formatting + background_formatting
plot5.1 <- ggplot(data = ds_b_sg_pr_control[ds_b_sg_pr_control$g=="pi",], 
                  aes(x = location, y = V1, group = location, color = g)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = location, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  facet_grid(cols = vars(dpi)) +
  labs(y = "Nucleotide diversity") + 
  scale_y_continuous(limits = c(0, 0.0015)) +
  facet_grid(cols = vars(dpi)) +
  theme(legend.title = element_blank(), 
        axis.title.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        legend.key = element_blank(), 
        legend.position = "none") +  
  axis_formatting + legend_formatting + background_formatting
plot6 <- ggplot(data = ds_b_sg_pr[ds_b_sg_pr$g=="pi",], 
                aes(x = dpi, y = V1, group = gld, color = g)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = dpi, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  labs(y = "Nucleotide diversity") + 
  scale_y_continuous(limits = c(0, 0.0015)) +
  facet_grid(cols = vars(gl)) +
  theme(legend.title = element_blank(), 
        axis.title.x.bottom = element_blank(),
        legend.key = element_blank(), 
        legend.position = "none") +  
  axis_formatting + legend_formatting + background_formatting
plot6.1 <- ggplot(data = ds_b_sg_pr_control[ds_b_sg_pr_control$g=="pi",], 
                   aes(x = location, y = V1, group = location, color = g)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = location, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  labs(y = "Nucleotide diversity") + 
  scale_y_continuous(limits = c(0, 0.0015)) +
  facet_grid(cols = vars(group)) +
  theme(legend.title = element_blank(), 
        axis.title.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        legend.key = element_blank(), 
        legend.position = "none") +  
  axis_formatting + legend_formatting + background_formatting
plot5_6 <- plot_grid(plot5, plot5.1, plot6, plot6.1, rel_widths = w, 
                     hjust = LR_adjust, vjust = UD_adjust, nrow = 1)

# piNpiS
ds_b_sg_pr$gl <- factor(paste(ds_b_sg_pr$group, ds_b_sg_pr$location, sep = "_"), 
                        levels = c("Tet_body", "Tet_legs", "Tet_saliva", "wmel_body", "wmel_legs"))
plot7 <- ggplot(data = ds_b_sg_pr[ds_b_sg_pr$g=="piNpiS",], 
                aes(x = gl, y = V1, group = gld, color = dpi)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = gl, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  labs(y = "PiN / PiS") + 
  scale_y_continuous(limits = c(0, 2)) + 
  geom_hline(yintercept = 1) + 
  facet_grid(cols = vars(dpi)) +
  theme(legend.title = element_blank(), 
        axis.title.x.bottom = element_blank(),
        legend.key = element_blank(), 
        legend.position = "none") +  
  axis_formatting + legend_formatting + background_formatting
plot7.1 <- ggplot(data = ds_b_sg_pr_control[ds_b_sg_pr_control$g=="piNpiS",], 
                  aes(x = location, y = V1, group = location, color = g)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = location, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  facet_grid(cols = vars(dpi)) +
  labs(y = "PiN / PiS") + 
  scale_y_continuous(limits = c(0, 2)) + 
  facet_grid(cols = vars(dpi)) +
  theme(legend.title = element_blank(), 
        axis.title.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        legend.key = element_blank(), 
        legend.position = "none") +  
  axis_formatting + legend_formatting + background_formatting
plot8 <- ggplot(data = ds_b_sg_pr[ds_b_sg_pr$g=="piNpiS",], 
                aes(x = dpi, y = V1, group = gld, color = g)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = dpi, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  labs(y = "PiN / PiS") + 
  scale_y_continuous(limits = c(0, 2)) + 
  geom_hline(yintercept = 1) + 
  facet_grid(cols = vars(gl)) +
  theme(legend.title = element_blank(), 
        axis.title.x.bottom = element_blank(),
        legend.key = element_blank(), 
        legend.position = "none") +  
  axis_formatting + legend_formatting + background_formatting
plot8.1 <- ggplot(data = ds_b_sg_pr_control[ds_b_sg_pr_control$g=="piNminuspiS",], 
                   aes(x = location, y = V1, group = gl, color = g)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = location, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  labs(y = "PiN - PiS") + 
  scale_y_continuous(limits = c(-.0015, 0.0015)) + 
  geom_hline(yintercept = 0) + 
  facet_grid(cols = vars(group)) +
  theme(legend.title = element_blank(), 
        axis.title.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        legend.key = element_blank(), 
        legend.position = "none") +  
  axis_formatting + legend_formatting + background_formatting
plot7_8 <- plot_grid(plot7, plot7.1, plot8, plot8.1, rel_widths = w, 
                     hjust = LR_adjust, vjust = UD_adjust, nrow = 1)

# piNminuspiS
ds_b_sg_pr$gl <- factor(paste(ds_b_sg_pr$group, ds_b_sg_pr$location, sep = "_"), 
                        levels = c("Tet_body", "Tet_legs", "Tet_saliva", "wmel_body", "wmel_legs"))
plot9 <- ggplot(data = ds_b_sg_pr[ds_b_sg_pr$g=="piNminuspiS",], 
                aes(x = gl, y = V1, group = gld, color = dpi)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = gl, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  labs(y = "PiN - PiS") + 
  scale_y_continuous(limits = c(-.0015, 0.0015)) + 
  geom_hline(yintercept = 0) + 
  facet_grid(cols = vars(dpi)) +
  theme(legend.title = element_blank(), 
        axis.title.x.bottom = element_blank(),
        legend.key = element_blank(), 
        legend.position = "none") +  
  axis_formatting + legend_formatting + background_formatting
plot9.1 <- ggplot(data = ds_b_sg_pr_control[ds_b_sg_pr_control$g=="piNminuspiS",], 
                aes(x = location, y = V1, group = gl, color = dpi)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = location, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  labs(y = "PiN - PiS") + 
  scale_y_continuous(limits = c(-.0015, 0.0015)) + 
  geom_hline(yintercept = 0) + 
  facet_grid(cols = vars(dpi)) +
  theme(legend.title = element_blank(), 
        axis.title.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        legend.key = element_blank(), 
        legend.position = "none") +  
  axis_formatting + legend_formatting + background_formatting
plot10 <- ggplot(data = ds_b_sg_pr[ds_b_sg_pr$g=="piNminuspiS",], 
                aes(x = dpi, y = V1, group = gld, color = g)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = dpi, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  labs(y = "PiN - PiS") + 
  scale_y_continuous(limits = c(-.0015, 0.0015)) + 
  geom_hline(yintercept = 0) + 
  facet_grid(cols = vars(gl)) +
  theme(legend.title = element_blank(), 
        axis.title.x.bottom = element_blank(),
        legend.key = element_blank(), 
        legend.position = "none") +  
  axis_formatting + legend_formatting + background_formatting
plot10.1 <- ggplot(data = ds_b_sg_pr_control[ds_b_sg_pr_control$g=="piNminuspiS",], 
                 aes(x = location, y = V1, group = gl, color = g)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = location, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  labs(y = "PiN - PiS") + 
  scale_y_continuous(limits = c(-.0015, 0.0015)) + 
  geom_hline(yintercept = 0) + 
  facet_grid(cols = vars(group)) +
  theme(legend.title = element_blank(), 
        axis.title.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        legend.key = element_blank(), 
        legend.position = "none") +  
  axis_formatting + legend_formatting + background_formatting
plot9_10 <- plot_grid(plot9, plot9.1, plot10, plot10.1, rel_widths = w, 
                     hjust = LR_adjust, vjust = UD_adjust, nrow = 1)
plot_snpgenie <- plot_grid(plot3_4, plot9_10, plot7_8, plot5_6, ncol = 1)

#### Plot diversity ####
## R cadm
df_R_cadm$gld <- paste(df_R_cadm$group, df_R_cadm$sample_type, df_R_cadm$dpi, sep = "_")
df_R_cadm$gld <- as.factor(df_R_cadm$gld)
levels <- levels(df_R_cadm$gld)

i=1
temp <- df_R_cadm[df_R_cadm$gld==levels[i] & df_R_cadm$dup=="not_dup",]
temp_mgs <- ds(data.frame("b" = b(temp$Mean.Gini.Simpson), "g" = "mgs"), varname = "V1", groupnames = "g")
temp_mse <- ds(data.frame("b" = b(temp$Mean.Shannon.Entropy), "g" = "mse"), varname = "V1", groupnames = "g")
temp_rmsd <- ds(data.frame("b" = b(temp$RMSD), "g" = "rmsd"), varname = "V1", groupnames = "g")
temp_umf <- ds(data.frame("b" = b(temp$Unique.Mut.Freq), "g" = "umf"), varname = "V1", groupnames = "g")
temp_mfp <- ds(data.frame("b" = b(temp$Mut.Freq.per.10K), "g" = "mfp"), varname = "V1", groupnames = "g")
ds_b_R <- rbind(temp_mgs, temp_mse, temp_rmsd, temp_umf, temp_mfp)
ds_b_R$gld <- levels[i]
for (i in 2:length(levels)) {
  temp <- df_R_cadm[df_R_cadm$gld==levels[i] & df_R_cadm$dup=="not_dup",]
  if (length(temp$MeanDepth)>1) {
    temp_mgs <- ds(data.frame("b" = b(temp$Mean.Gini.Simpson), "g" = "mgs"), varname = "V1", groupnames = "g")
    temp_mse <- ds(data.frame("b" = b(temp$Mean.Shannon.Entropy), "g" = "mse"), varname = "V1", groupnames = "g")
    temp_rmsd <- ds(data.frame("b" = b(temp$RMSD), "g" = "rmsd"), varname = "V1", groupnames = "g")
    temp_umf <- ds(data.frame("b" = b(temp$Unique.Mut.Freq), "g" = "umf"), varname = "V1", groupnames = "g")
    temp_mfp <- ds(data.frame("b" = b(temp$Mut.Freq.per.10K), "g" = "mfp"), varname = "V1", groupnames = "g")
    temp2 <- rbind(temp_mgs, temp_mse, temp_rmsd, temp_umf, temp_mfp)
    temp2$gld <- levels[i]
    ds_b_R <- rbind(ds_b_R, temp2)
  }
}
ds_b_R <- separate(ds_b_R, "gld", c("group", "location", "dpi"), sep = "_")
ds_b_R$gl <- factor(paste(ds_b_R$group, ds_b_R$location, sep = "_"), 
                        levels = c("tet_saliva", "tet_body", "wmel_body", "tet_legs", "wmel_legs"))
ds_b_R$gld <- as.factor(paste(ds_b_R$group, ds_b_R$location, ds_b_R$dpi, sep = "_"))
ds_b_R$dpi <- factor(ds_b_R$dpi, levels = c("4", "7", "14"))
ds_b_R$g <- factor(ds_b_R$g, levels = c("mgs", "mse", "rmsd", "umf", "mfp"))

# mgs
ds_b_R$gl <- factor(paste(ds_b_R$group, ds_b_R$location, sep = "_"), 
                        levels = c("tet_saliva", "tet_body", "wmel_body", "tet_legs", "wmel_legs"))
plot9 <- ggplot(data = ds_b_R[ds_b_R$g=="mgs",], 
                aes(x = gl, y = V1, group = gld, color = g)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = gl, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  facet_grid(cols = vars(dpi)) +
  theme(legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.position = "none") +  
  axis_formatting + legend_formatting + background_formatting

ds_b_R$gl <- factor(paste(ds_b_R$group, ds_b_R$location, sep = "_"), 
                        levels = c("tet_saliva", "tet_body", "tet_legs","wmel_body", "wmel_legs"))
plot10 <- ggplot(data = ds_b_R[ds_b_R$g=="mgs",], 
                aes(x = dpi, y = V1, group = gld, color = g)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = dpi, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  facet_grid(cols = vars(gl)) +
  theme(legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.position = "right") +  
  axis_formatting + legend_formatting + background_formatting
plot_mgs <- plot_grid(plot9, plot10, rel_widths = w, 
                     hjust = LR_adjust, vjust = UD_adjust)

# mse
ds_b_R$gl <- factor(paste(ds_b_R$group, ds_b_R$location, sep = "_"), 
                    levels = c("tet_saliva", "tet_body", "wmel_body", "tet_legs", "wmel_legs"))
plot11 <- ggplot(data = ds_b_R[ds_b_R$g=="mse",], 
                aes(x = gl, y = V1, group = gld, color = g)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = gl, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  facet_grid(cols = vars(dpi)) +
  theme(legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.position = "none") +  
  axis_formatting + legend_formatting + background_formatting

ds_b_R$gl <- factor(paste(ds_b_R$group, ds_b_R$location, sep = "_"), 
                    levels = c("tet_saliva", "tet_body", "tet_legs","wmel_body", "wmel_legs"))
plot12 <- ggplot(data = ds_b_R[ds_b_R$g=="mse",], 
                 aes(x = dpi, y = V1, group = gld, color = g)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = dpi, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  facet_grid(cols = vars(gl)) +
  theme(legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.position = "right") +  
  axis_formatting + legend_formatting + background_formatting
plot_mse <- plot_grid(plot11, plot12, rel_widths = w, 
                      hjust = LR_adjust, vjust = UD_adjust)

# umf
ds_b_R$gl <- factor(paste(ds_b_R$group, ds_b_R$location, sep = "_"), 
                    levels = c("tet_saliva", "tet_body", "wmel_body", "tet_legs", "wmel_legs"))
plot13 <- ggplot(data = ds_b_R[ds_b_R$g=="umf",], 
                 aes(x = gl, y = V1, group = gld, color = g)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = gl, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  facet_grid(cols = vars(dpi)) +
  theme(legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.position = "none") +  
  axis_formatting + legend_formatting + background_formatting

ds_b_R$gl <- factor(paste(ds_b_R$group, ds_b_R$location, sep = "_"), 
                    levels = c("tet_saliva", "tet_body", "tet_legs","wmel_body", "wmel_legs"))
plot14 <- ggplot(data = ds_b_R[ds_b_R$g=="umf",], 
                 aes(x = dpi, y = V1, group = gld, color = g)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = dpi, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  facet_grid(cols = vars(gl)) +
  theme(legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.position = "right") +  
  axis_formatting + legend_formatting + background_formatting
plot_umf <- plot_grid(plot13, plot14, rel_widths = w, 
                      hjust = LR_adjust, vjust = UD_adjust)

# mfp
ds_b_R$gl <- factor(paste(ds_b_R$group, ds_b_R$location, sep = "_"), 
                    levels = c("tet_saliva", "tet_body", "wmel_body", "tet_legs", "wmel_legs"))
plot15 <- ggplot(data = ds_b_R[ds_b_R$g=="mfp",], 
                 aes(x = gl, y = V1, group = gld, color = g)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = gl, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  facet_grid(cols = vars(dpi)) +
  theme(legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.position = "none") +  
  axis_formatting + legend_formatting + background_formatting

ds_b_R$gl <- factor(paste(ds_b_R$group, ds_b_R$location, sep = "_"), 
                    levels = c("tet_saliva", "tet_body", "tet_legs","wmel_body", "wmel_legs"))
plot16 <- ggplot(data = ds_b_R[ds_b_R$g=="mfp",], 
                 aes(x = dpi, y = V1, group = gld, color = g)) + 
  geom_point(position = position_dodge(.9)) + 
  geom_errorbar(aes(x = dpi, y = V1, ymin = V1 - sd, ymax = V1 + sd), 
                width = .2, position = position_dodge(.9)) + 
  facet_grid(cols = vars(gl)) +
  theme(legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.position = "right") +  
  axis_formatting + legend_formatting + background_formatting
plot_mfp <- plot_grid(plot15, plot16, rel_widths = w, 
                      hjust = LR_adjust, vjust = UD_adjust)


plot_R <- plot_grid(plot_mgs, plot_mse, plot_umf, plot_mfp, ncol = 1)


#### save ####
setwd(dir_save)
#Fig1
ggsave("Fig1_snpgenie.pdf", plot_snpgenie,
       width = 15, height = 5, 
       units = "in", dpi = 320)
#Fig2
ggsave("Fig2_diversity.pdf", plot_R,
       width = 15, height = 5, 
       units = "in", dpi = 320)
#Fig3
  ## barcodes on 2-barcode_freqs.R
#Fig4
ggsave("Fig4_iSNVs_noBC.pdf", plot_iSNVs_all,
       width = 8, height = 6, 
       units = "in", dpi = 320)
#Fig5
ggsave("Fig5_spectrum_facet.pdf", plot_spec,
       width = 9, height = 9,
       units = "in", dpi = 320)
#Fig6
ggsave("Fig5_divergence.pdf", plot_divergence_counts,
       width = 5, height = 4,
       units = "in", dpi = 320)
#### Misc ####
#df_09_reference_VCF$group_location_dpi <- factor(df_09_reference_VCF$group_location_dpi,
#                                                 levels = c("Control_mouse_Control",
#                                                            "Control_PC_Control",
#                                                            "Tet_saliva_7", "Tet_saliva_14",
#                                                            "Tet_body_4", "Tet_body_7", "Tet_body_14",
#                                                            "Tet_legs_7", "Tet_legs_14",
#                                                            "wmel_body_7", "wmel_body_14",
#                                                            "wmel_legs_7", "wmel_legs_14"))
#ggplot(df_09_reference_VCF[df_09_reference_VCF$HGVS.p=="p.Pro1683Pro",], aes(x = group_location_dpi, y = AF)) + 
#  geom_point() + theme_bw()