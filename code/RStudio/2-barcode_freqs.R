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
if(!require(ggdendro)){
  install.packages("ggdendro",dependencies = T)
  library(ggdendro)
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
# e.g., x <- filter(x, y %!in% c("unwantedValue1", "unwantedValue2"))


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
dir_15_b <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/Wolbachia/data/reads/data/run/15_barcode_analyses", sep="")
dir_save <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/Wolbachia/figs", sep="")
setwd(dir_15_b); dir()

#### Mean_bc_frequencies.R ####
# R script to combine barcode frequencies from technical duplicate sequencing libraries.
dir_kmercounts <- dir(pattern="counts_BConly_*")
names_trunc <- gsub(".tsv", "", dir_kmercounts)
names_trunc <- gsub("counts_BConly_", "", names_trunc)

# lists
n <- length(dir_kmercounts)
list <- vector("list", n)
list2 <- vector("list", n)

# data frame for totals
sample_name <- rep("foo", n)
sample_total <- rep(0, n)
df_umi_totals <- data.frame(sample_name, sample_total)

# for-loop
for (i in 1:n) {
  list[[i]] <- as.data.frame(
    read_tsv(dir_kmercounts[[i]], show_col_types = F))
  list[[i]] <- list[[i]][c(1,2)]
  colnames(list[[i]]) <- c("count", "barcode")
  list2[[i]]$sample_name <- names_trunc[i]
  ifelse(length(list[[i]]$count)>0, 
         list2[[i]]$sample_total <- sum(list[[i]]$count),
         list2[[i]]$sample_total <- NA)
  list2[[i]]$sample_name <- names_trunc[i]
  names(list) <- names_trunc
  names(list2) <- names_trunc
  df_umi_totals$sample_name[i]  <- list2[[i]]$sample_name
  df_umi_totals$sample_total[i] <- list2[[i]]$sample_total
  list[[i]]$freq <- rep(0, length(list[[i]]$count))
  ifelse(length(list[[i]]$count)>0, 
         list[[i]]$freq <- list[[i]]$count / list2[[i]]$sample_total,
         print("No data!"))
  ifelse(length(list[[i]]$count)>0, 
         list[[i]]$sample_name <- names_trunc[i],
         print("No data!"))
}

## remove empty dfs from list
list <- Filter(function(x) dim(x)[1] > 0, list)

## join
df_barcode_freqs <- Reduce(full_join,list)



#### Euclidean_distance_meanbcfreqs.R ####
 # R script to calculate Euclidean distance between barcode populations of two samples.
all_bcs <- as.data.frame(table(df_barcode_freqs$barcode))

## pairwise comparisons
n <- length(list)
x = 1
df_compiled <- data.frame("foo", "bar", 0, 0)
names(df_compiled) <- c("sample_1", "sample_2", "sumfreqdiffsq", "eucdistance")
df_compiled <- df_compiled[0,]
start = Sys.time()
for (i in 2:n) {
  for (j in 1:x) {
    ## merge and zero-out no-matches
    df_temp <- merge(list[[i]], list[[j]], by = "barcode", all.x = TRUE, all.y = TRUE)
    df_temp[is.na(df_temp)]<-0
    ## squared difference in frequencies
    df_temp$freqdiffsq<-sum(df_temp$freq.x-df_temp$freq.y)^2
    ## Sum the squared differences
    sumfreqdiffsq<-sum(df_temp$freqdiffsq)
    ## Take sq root of total differences squared to calculate Euclidean distance
    eucdistance<-sqrt(sumfreqdiffsq); eucdistance
    ## Create row for final df
    df_temp1 <- data.frame(
      list[[i]]$sample_name[1],
      list[[j]]$sample_name[1],
      sumfreqdiffsq,
      eucdistance
    )
    names(df_temp1) <- c("sample_1", "sample_2", "sumfreqdiffsq", "eucdistance")
    df_compiled <- rbind(df_compiled, df_temp1)
  }
  out <- paste(i, "/", n); print(out)
  x = x + 1
}
finish <- Sys.time()
out1 <- paste("Start: ", start); out2 <- paste("Finish:", finish); print(out1); print(out2)
df_pairwise_euc <- df_compiled


df_pairwise_euc$e <- log10(df_pairwise_euc$eucdistance)
df_pairwise_euc$e[df_pairwise_euc$e=="-Inf"] <- 0

## match string
# sample 1
df_pairwise_euc$sample_1[grepl("7-2-tet-saliva",df_pairwise_euc$sample_1)] <- "Tet-saliva-7dpi_2"
df_pairwise_euc$sample_1[grepl("7-3-tet-saliva",df_pairwise_euc$sample_1)] <- "Tet-saliva-7dpi_3"
df_pairwise_euc$sample_1[grepl("14-1-tet-saliva",df_pairwise_euc$sample_1)] <- "Tet-saliva-14dpi_1"
df_pairwise_euc$sample_1[grepl("14-3-tet-saliva",df_pairwise_euc$sample_1)] <- "Tet-saliva-14dpi_3"
df_pairwise_euc$sample_1[grepl("7-1-tet-legs",df_pairwise_euc$sample_1)] <- "Tet-legs-7dpi_1"
df_pairwise_euc$sample_1[grepl("7-2-tet-legs",df_pairwise_euc$sample_1)] <- "Tet-legs-7dpi_2"
df_pairwise_euc$sample_1[grepl("7-3-tet-legs",df_pairwise_euc$sample_1)] <- "Tet-legs-7dpi_3"
df_pairwise_euc$sample_1[grepl("14-1-tet-legs",df_pairwise_euc$sample_1)] <- "Tet-legs-14dpi_1"
df_pairwise_euc$sample_1[grepl("14-2-tet-legs",df_pairwise_euc$sample_1)] <- "Tet-legs-14dpi_2"
df_pairwise_euc$sample_1[grepl("14-3-tet-legs",df_pairwise_euc$sample_1)] <- "Tet-legs-14dpi_3"
df_pairwise_euc$sample_1[grepl("4-1-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-4dpi_1"
df_pairwise_euc$sample_1[grepl("4-2-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-4dpi_2"
df_pairwise_euc$sample_1[grepl("4-3-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-4dpi_3"
df_pairwise_euc$sample_1[grepl("7-1-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-7dpi_1"
df_pairwise_euc$sample_1[grepl("7-2-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-7dpi_2"
df_pairwise_euc$sample_1[grepl("7-3-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-7dpi_3"
df_pairwise_euc$sample_1[grepl("dup-14-1-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-14dpi_1_dup"
df_pairwise_euc$sample_1[grepl("14-1-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-14dpi_1"
df_pairwise_euc$sample_1[grepl("dup-14-2-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-14dpi_2_dup"
df_pairwise_euc$sample_1[grepl("14-2-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-14dpi_2"
df_pairwise_euc$sample_1[grepl("dup-14-3-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-14dpi_3_dup"
df_pairwise_euc$sample_1[grepl("14-3-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-14dpi_3"
df_pairwise_euc$sample_1[grepl("7-1-wmel-body",df_pairwise_euc$sample_1)] <- "wmel-body-7dpi_1"
df_pairwise_euc$sample_1[grepl("7-2-wmel-body",df_pairwise_euc$sample_1)] <- "wmel-body-7dpi_2"
df_pairwise_euc$sample_1[grepl("7-3-wmel-body",df_pairwise_euc$sample_1)] <- "wmel-body-7dpi_3"
df_pairwise_euc$sample_1[grepl("14-1-wmel-body",df_pairwise_euc$sample_1)] <- "wmel-body-14dpi_1"
df_pairwise_euc$sample_1[grepl("14-2-wmel-body",df_pairwise_euc$sample_1)] <- "wmel-body-14dpi_2"
df_pairwise_euc$sample_1[grepl("14-3-wmel-body",df_pairwise_euc$sample_1)] <- "wmel-body-14dpi_3"
df_pairwise_euc$sample_1[grepl("7-1-wmel-legs",df_pairwise_euc$sample_1)] <- "wmel-legs-7dpi_1"
df_pairwise_euc$sample_1[grepl("14-2-wmel-legs",df_pairwise_euc$sample_1)] <- "wmel-legs-14dpi_2"
df_pairwise_euc$sample_1[grepl("7-2-tet-saliva",df_pairwise_euc$sample_1)] <- "Tet-saliva-7dpi_2"
df_pairwise_euc$sample_1[grepl("7-3-tet-saliva",df_pairwise_euc$sample_1)] <- "Tet-saliva-7dpi_3"
df_pairwise_euc$sample_1[grepl("14-1-tet-saliva",df_pairwise_euc$sample_1)] <- "Tet-saliva-14dpi_1"
df_pairwise_euc$sample_1[grepl("14-3-tet-saliva",df_pairwise_euc$sample_1)] <- "Tet-saliva-14dpi_3"
df_pairwise_euc$sample_1[grepl("7-1-tet-legs",df_pairwise_euc$sample_1)] <- "Tet-legs-7dpi_1"
df_pairwise_euc$sample_1[grepl("7-2-tet-legs",df_pairwise_euc$sample_1)] <- "Tet-legs-7dpi_2"
df_pairwise_euc$sample_1[grepl("7-3-tet-legs",df_pairwise_euc$sample_1)] <- "Tet-legs-7dpi_3"
df_pairwise_euc$sample_1[grepl("14-1-tet-legs",df_pairwise_euc$sample_1)] <- "Tet-legs-14dpi_1"
df_pairwise_euc$sample_1[grepl("14-2-tet-legs",df_pairwise_euc$sample_1)] <- "Tet-legs-14dpi_2"
df_pairwise_euc$sample_1[grepl("14-3-tet-legs",df_pairwise_euc$sample_1)] <- "Tet-legs-14dpi_3"
df_pairwise_euc$sample_1[grepl("4-1-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-4dpi_1"
df_pairwise_euc$sample_1[grepl("4-2-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-4dpi_2"
df_pairwise_euc$sample_1[grepl("4-3-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-4dpi_3"
df_pairwise_euc$sample_1[grepl("7-1-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-7dpi_1"
df_pairwise_euc$sample_1[grepl("7-2-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-7dpi_2"
df_pairwise_euc$sample_1[grepl("7-3-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-7dpi_3"
df_pairwise_euc$sample_1[grepl("dup-14-1-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-14dpi_1_dup"
df_pairwise_euc$sample_1[grepl("14-1-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-14dpi_1"
df_pairwise_euc$sample_1[grepl("dup-14-2-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-14dpi_2_dup"
df_pairwise_euc$sample_1[grepl("14-2-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-14dpi_2"
df_pairwise_euc$sample_1[grepl("dup-14-3-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-14dpi_3_dup"
df_pairwise_euc$sample_1[grepl("14-3-tet-body",df_pairwise_euc$sample_1)] <- "Tet-body-14dpi_3"
df_pairwise_euc$sample_1[grepl("7-1-wmel-body",df_pairwise_euc$sample_1)] <- "wmel-body-7dpi_1"
df_pairwise_euc$sample_1[grepl("7-2-wmel-body",df_pairwise_euc$sample_1)] <- "wmel-body-7dpi_2"
df_pairwise_euc$sample_1[grepl("7-3-wmel-body",df_pairwise_euc$sample_1)] <- "wmel-body-7dpi_3"
df_pairwise_euc$sample_1[grepl("14-1-wmel-body",df_pairwise_euc$sample_1)] <- "wmel-body-14dpi_1"
df_pairwise_euc$sample_1[grepl("14-2-wmel-body",df_pairwise_euc$sample_1)] <- "wmel-body-14dpi_2"
df_pairwise_euc$sample_1[grepl("14-3-wmel-body",df_pairwise_euc$sample_1)] <- "wmel-body-14dpi_3"
df_pairwise_euc$sample_1[grepl("7-1-wmel-legs",df_pairwise_euc$sample_1)] <- "wmel-legs-7dpi_1"
df_pairwise_euc$sample_1[grepl("14-2-wmel-legs",df_pairwise_euc$sample_1)] <- "wmel-legs-14dpi_2"

# sample 2
df_pairwise_euc$sample_2[grepl("7-2-tet-saliva",df_pairwise_euc$sample_2)] <- "Tet-saliva-7dpi_2"
df_pairwise_euc$sample_2[grepl("7-3-tet-saliva",df_pairwise_euc$sample_2)] <- "Tet-saliva-7dpi_3"
df_pairwise_euc$sample_2[grepl("14-1-tet-saliva",df_pairwise_euc$sample_2)] <- "Tet-saliva-14dpi_1"
df_pairwise_euc$sample_2[grepl("14-3-tet-saliva",df_pairwise_euc$sample_2)] <- "Tet-saliva-14dpi_3"
df_pairwise_euc$sample_2[grepl("7-1-tet-legs",df_pairwise_euc$sample_2)] <- "Tet-legs-7dpi_1"
df_pairwise_euc$sample_2[grepl("7-2-tet-legs",df_pairwise_euc$sample_2)] <- "Tet-legs-7dpi_2"
df_pairwise_euc$sample_2[grepl("7-3-tet-legs",df_pairwise_euc$sample_2)] <- "Tet-legs-7dpi_3"
df_pairwise_euc$sample_2[grepl("14-1-tet-legs",df_pairwise_euc$sample_2)] <- "Tet-legs-14dpi_1"
df_pairwise_euc$sample_2[grepl("14-2-tet-legs",df_pairwise_euc$sample_2)] <- "Tet-legs-14dpi_2"
df_pairwise_euc$sample_2[grepl("14-3-tet-legs",df_pairwise_euc$sample_2)] <- "Tet-legs-14dpi_3"
df_pairwise_euc$sample_2[grepl("4-1-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-4dpi_1"
df_pairwise_euc$sample_2[grepl("4-2-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-4dpi_2"
df_pairwise_euc$sample_2[grepl("4-3-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-4dpi_3"
df_pairwise_euc$sample_2[grepl("7-1-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-7dpi_1"
df_pairwise_euc$sample_2[grepl("7-2-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-7dpi_2"
df_pairwise_euc$sample_2[grepl("7-3-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-7dpi_3"
df_pairwise_euc$sample_2[grepl("dup-14-1-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-14dpi_1_dup"
df_pairwise_euc$sample_2[grepl("14-1-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-14dpi_1"
df_pairwise_euc$sample_2[grepl("dup-14-2-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-14dpi_2_dup"
df_pairwise_euc$sample_2[grepl("14-2-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-14dpi_2"
df_pairwise_euc$sample_2[grepl("dup-14-3-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-14dpi_3_dup"
df_pairwise_euc$sample_2[grepl("14-3-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-14dpi_3"
df_pairwise_euc$sample_2[grepl("7-1-wmel-body",df_pairwise_euc$sample_2)] <- "wmel-body-7dpi_1"
df_pairwise_euc$sample_2[grepl("7-2-wmel-body",df_pairwise_euc$sample_2)] <- "wmel-body-7dpi_2"
df_pairwise_euc$sample_2[grepl("7-3-wmel-body",df_pairwise_euc$sample_2)] <- "wmel-body-7dpi_3"
df_pairwise_euc$sample_2[grepl("14-1-wmel-body",df_pairwise_euc$sample_2)] <- "wmel-body-14dpi_1"
df_pairwise_euc$sample_2[grepl("14-2-wmel-body",df_pairwise_euc$sample_2)] <- "wmel-body-14dpi_2"
df_pairwise_euc$sample_2[grepl("14-3-wmel-body",df_pairwise_euc$sample_2)] <- "wmel-body-14dpi_3"
df_pairwise_euc$sample_2[grepl("7-1-wmel-legs",df_pairwise_euc$sample_2)] <- "wmel-legs-7dpi_1"
df_pairwise_euc$sample_2[grepl("14-2-wmel-legs",df_pairwise_euc$sample_2)] <- "wmel-legs-14dpi_2"
df_pairwise_euc$sample_2[grepl("7-2-tet-saliva",df_pairwise_euc$sample_2)] <- "Tet-saliva-7dpi_2"
df_pairwise_euc$sample_2[grepl("7-3-tet-saliva",df_pairwise_euc$sample_2)] <- "Tet-saliva-7dpi_3"
df_pairwise_euc$sample_2[grepl("14-1-tet-saliva",df_pairwise_euc$sample_2)] <- "Tet-saliva-14dpi_1"
df_pairwise_euc$sample_2[grepl("14-3-tet-saliva",df_pairwise_euc$sample_2)] <- "Tet-saliva-14dpi_3"
df_pairwise_euc$sample_2[grepl("7-1-tet-legs",df_pairwise_euc$sample_2)] <- "Tet-legs-7dpi_1"
df_pairwise_euc$sample_2[grepl("7-2-tet-legs",df_pairwise_euc$sample_2)] <- "Tet-legs-7dpi_2"
df_pairwise_euc$sample_2[grepl("7-3-tet-legs",df_pairwise_euc$sample_2)] <- "Tet-legs-7dpi_3"
df_pairwise_euc$sample_2[grepl("14-1-tet-legs",df_pairwise_euc$sample_2)] <- "Tet-legs-14dpi_1"
df_pairwise_euc$sample_2[grepl("14-2-tet-legs",df_pairwise_euc$sample_2)] <- "Tet-legs-14dpi_2"
df_pairwise_euc$sample_2[grepl("14-3-tet-legs",df_pairwise_euc$sample_2)] <- "Tet-legs-14dpi_3"
df_pairwise_euc$sample_2[grepl("4-1-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-4dpi_1"
df_pairwise_euc$sample_2[grepl("4-2-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-4dpi_2"
df_pairwise_euc$sample_2[grepl("4-3-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-4dpi_3"
df_pairwise_euc$sample_2[grepl("7-1-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-7dpi_1"
df_pairwise_euc$sample_2[grepl("7-2-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-7dpi_2"
df_pairwise_euc$sample_2[grepl("7-3-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-7dpi_3"
df_pairwise_euc$sample_2[grepl("dup-14-1-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-14dpi_1_dup"
df_pairwise_euc$sample_2[grepl("14-1-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-14dpi_1"
df_pairwise_euc$sample_2[grepl("dup-14-2-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-14dpi_2_dup"
df_pairwise_euc$sample_2[grepl("14-2-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-14dpi_2"
df_pairwise_euc$sample_2[grepl("dup-14-3-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-14dpi_3_dup"
df_pairwise_euc$sample_2[grepl("14-3-tet-body",df_pairwise_euc$sample_2)] <- "Tet-body-14dpi_3"
df_pairwise_euc$sample_2[grepl("7-1-wmel-body",df_pairwise_euc$sample_2)] <- "wmel-body-7dpi_1"
df_pairwise_euc$sample_2[grepl("7-2-wmel-body",df_pairwise_euc$sample_2)] <- "wmel-body-7dpi_2"
df_pairwise_euc$sample_2[grepl("7-3-wmel-body",df_pairwise_euc$sample_2)] <- "wmel-body-7dpi_3"
df_pairwise_euc$sample_2[grepl("14-1-wmel-body",df_pairwise_euc$sample_2)] <- "wmel-body-14dpi_1"
df_pairwise_euc$sample_2[grepl("14-2-wmel-body",df_pairwise_euc$sample_2)] <- "wmel-body-14dpi_2"
df_pairwise_euc$sample_2[grepl("14-3-wmel-body",df_pairwise_euc$sample_2)] <- "wmel-body-14dpi_3"
df_pairwise_euc$sample_2[grepl("7-1-wmel-legs",df_pairwise_euc$sample_2)] <- "wmel-legs-7dpi_1"
df_pairwise_euc$sample_2[grepl("14-2-wmel-legs",df_pairwise_euc$sample_2)] <- "wmel-legs-14dpi_2"

## remove controls
# 1
df_pairwise_euc_NC_1 <- df_pairwise_euc[grepl("NC",df_pairwise_euc$sample_1),]
df_pairwise_euc_PC_1 <- df_pairwise_euc[grepl("PC",df_pairwise_euc$sample_1),]
df_pairwise_euc_ZIKV_1 <- df_pairwise_euc[grepl("ZIKV",df_pairwise_euc$sample_1),]
df_pairwise_euc_mouse_1 <- df_pairwise_euc[grepl("mouse",df_pairwise_euc$sample_1),]
df_pairwise_euc <- df_pairwise_euc[!grepl("NC",df_pairwise_euc$sample_1),]
df_pairwise_euc <- df_pairwise_euc[!grepl("PC",df_pairwise_euc$sample_1),]
df_pairwise_euc <- df_pairwise_euc[!grepl("ZIKV",df_pairwise_euc$sample_1),]
df_pairwise_euc <- df_pairwise_euc[!grepl("mouse",df_pairwise_euc$sample_1),]
# 2
df_pairwise_euc_NC_2 <- df_pairwise_euc[grepl("NC",df_pairwise_euc$sample_2),]
df_pairwise_euc_PC_2 <- df_pairwise_euc[grepl("PC",df_pairwise_euc$sample_2),]
df_pairwise_euc_ZIKV_2 <- df_pairwise_euc[grepl("ZIKV",df_pairwise_euc$sample_2),]
df_pairwise_euc_mouse_2 <- df_pairwise_euc[grepl("mouse",df_pairwise_euc$sample_2),]
df_pairwise_euc <- df_pairwise_euc[!grepl("NC",df_pairwise_euc$sample_2),]
df_pairwise_euc <- df_pairwise_euc[!grepl("PC",df_pairwise_euc$sample_2),]
df_pairwise_euc <- df_pairwise_euc[!grepl("ZIKV",df_pairwise_euc$sample_2),]
df_pairwise_euc <- df_pairwise_euc[!grepl("mouse",df_pairwise_euc$sample_2),]

## heatmap
#exclude <- c("B1-2-ZIKV-seq-2-2-mouse2-serum_S29", "B1-1-ZIKV-seq-2-2-mouse1-serum_S30", 
#             "PC-ZIKV-PR-SetA_S21", "PC-ZIKV-SetC_S21", "PC-ZIKV-PR-SetD_S21", 
#             "PC-ZIKV-PR-SetA_S25", "PC-ZIKV-PR-SetB_S43", "PC-ZIKV-PR-SetE_S43", "NC-H2O-SetC_S22")
#df_pairwise_euc_filtered <- filter(df_pairwise_euc, sample_1 %!in% exclude)
#df_pairwise_euc_filtered <- filter(df_pairwise_euc_filtered, sample_2 %!in% exclude)
#df_pairwise_euc_filtered$sample_2 <- factor(df_pairwise_euc_filtered$sample_2, 
#                                            levels = c("Tet-saliva-7dpi_3", "Tet-saliva-7dpi_2",
#                                                       "Tet-body-4dpi_1", "Tet-body-4dpi_2", "Tet-body-4dpi_3", 
#                                                       "Tet-body-7dpi_1","Tet-body-7dpi_2", "Tet-body-7dpi_3", 
#                                                       "Tet-legs-7dpi_1", "Tet-legs-7dpi_2", "Tet-legs-7dpi_3", 
#                                                       "Tet-legs-14dpi_1", "Tet-legs-14dpi_3",
#                                                       "Tet-saliva-14dpi_1","Tet-saliva-14dpi_3", 
#                                                       "wmel-body-7dpi_1", "wmel-body-7dpi_2",
#                                                       "wmel-body-14dpi_1", "wmel-body-14dpi_2", "wmel-body-14dpi_3", 
#                                                       "wmel-legs-7dpi_1", 
#                                                       "wmel-legs-14dpi_2"))
#df_pairwise_euc_filtered$sample_1 <- factor(df_pairwise_euc_filtered$sample_1, 
#                                            levels = c("wmel-legs-14dpi_2",
#                                                       "wmel-legs-7dpi_1", 
#                                                       "wmel-body-14dpi_1", "wmel-body-14dpi_2", "wmel-body-14dpi_3", 
#                                                       "wmel-body-7dpi_1", "wmel-body-7dpi_2",
#                                                       "Tet-saliva-14dpi_1","Tet-saliva-14dpi_3", 
#                                                       "Tet-legs-14dpi_1", "Tet-legs-14dpi_3",
#                                                       "Tet-legs-7dpi_1", "Tet-legs-7dpi_2", "Tet-legs-7dpi_3", 
#                                                       "Tet-body-7dpi_1","Tet-body-7dpi_2", "Tet-body-7dpi_3", 
#                                                       "Tet-body-4dpi_1", "Tet-body-4dpi_2", "Tet-body-4dpi_3", 
#                                                       "Tet-saliva-7dpi_2", "Tet-saliva-7dpi_3"))
#
#test2 <- as.matrix(xtabs(
#  eucdistance ~ sample_1 + sample_2, 
#  data=df_pairwise_euc_filtered))
#heatmap(test2, Colv = NA, Rowv = NA)


df_pairwise_euc %>% 
  select(sample_1, sample_2, eucdistance) %>% 
  group_by(sample_1, sample_2) %>% 
  summarize(mean.euc=mean(eucdistance)) -> df_pairwise_euc_mean
#df_pairwise_euc_mean %>% 
#  ggplot(data=., mapping=aes(x=sample_1, y=sample_2, fill=mean.euc)) + 
#  geom_tile() +
#  theme(axis.text.x=element_text(color=rainbow(ncol(df_pairwise_euc_mean))), 
#        axis.text.y=element_text(color=rainbow(ncol(df_pairwise_euc_mean))))


#### Barcodes_overtime_COMPOSITE.R ####
# R script to plot composite barcode trajectories for each passage series.
df_tet_umi_totals <- df_umi_totals[grepl("tet",df_umi_totals$sample_name),]
df_tet_umi_totals <- separate(df_tet_umi_totals, "sample_name", c("sample_name","S"), sep="_S")
df_tet_umi_totals <- separate(df_tet_umi_totals, "sample_name", c("sample_name","location"), sep="-tet-")
df_tet_umi_totals_dups <- df_tet_umi_totals[grepl("dup",df_tet_umi_totals$sample_name),]
df_tet_umi_totals_dups <- separate(df_tet_umi_totals_dups, "sample_name", c("1","2","dup","3","4"), sep="-")
df_tet_umi_totals <- df_tet_umi_totals[!grepl("dup",df_tet_umi_totals$sample_name),]
df_tet_umi_totals <- separate(df_tet_umi_totals, "sample_name", c("1","2","3", "4"), sep="-")
df_tet_umi_totals$dup <- "not_dup"
df_tet_umi_totals <- rbind(df_tet_umi_totals, df_tet_umi_totals_dups)
df_tet_umi_totals$group <- "tet"
head(df_tet_umi_totals)


df_wmel_umi_totals <- df_umi_totals[grepl("wmel",df_umi_totals$sample_name),]
df_wmel_umi_totals <- separate(df_wmel_umi_totals, "sample_name", c("sample_name","S"), sep="_S")
df_wmel_umi_totals <- separate(df_wmel_umi_totals, "sample_name", c("sample_name","location"), sep="-wmel-")
df_wmel_umi_totals_dups <- df_wmel_umi_totals[grepl("dup",df_wmel_umi_totals$sample_name),]
df_wmel_umi_totals_dups <- separate(df_wmel_umi_totals_dups, "sample_name", c("1","2","dup","3","4"), sep="-")
df_wmel_umi_totals <- df_wmel_umi_totals[!grepl("dup",df_wmel_umi_totals$sample_name),]
df_wmel_umi_totals <- separate(df_wmel_umi_totals, "sample_name", c("1","2","3", "4"), sep="-")
df_wmel_umi_totals$dup <- "not_dup"
df_wmel_umi_totals <- rbind(df_wmel_umi_totals, df_tet_umi_totals_dups)
df_wmel_umi_totals$group <- "wmel"
head(df_wmel_umi_totals)

df_controls_umi_totals <- df_umi_totals[grepl("ZIKV",df_umi_totals$sample_name),]
df_controls_umi_totals$sample_name[grepl("mouse",df_controls_umi_totals$sample_name)] <- "Control_mouse"
df_controls_umi_totals$sample_name[grepl("PC",df_controls_umi_totals$sample_name)] <- "Control_PC"
df_controls_umi_totals$S <- ""
df_controls_umi_totals <- separate(df_controls_umi_totals, "sample_name", c("sample_name","location"), sep="_")
df_controls_umi_totals$dup <- "not_dup"
df_controls_umi_totals$group <- "Control"
df_controls_umi_totals$`1` <- "Control"
df_controls_umi_totals$`2` <- df_controls_umi_totals$location
df_controls_umi_totals$`3` <- "Control"
df_controls_umi_totals$`4` <- "Control"
df_controls_umi_totals$sample_name <- NULL
head(df_controls_umi_totals)

df_treatment_umi_totals <- rbind(df_tet_umi_totals, df_wmel_umi_totals, df_controls_umi_totals)
df_treatment_umi_totals$ID <- paste(df_treatment_umi_totals$`1`,df_treatment_umi_totals$`2`, sep = "_")
df_treatment_umi_totals$rep <- df_treatment_umi_totals$`4`
df_treatment_umi_totals$dpi <- df_treatment_umi_totals$`3`
df_treatment_umi_totals$group_location <- paste(df_treatment_umi_totals$group, 
                                                df_treatment_umi_totals$location, sep = "_")
df_treatment_umi_totals$group_location_dpi <- paste(df_treatment_umi_totals$group, 
                                                    df_treatment_umi_totals$location, 
                                                    df_treatment_umi_totals$dpi, 
                                                    sep = "_")
df_treatment_umi_totals$`1` <- NULL
df_treatment_umi_totals$`2` <- NULL
df_treatment_umi_totals$`3` <- NULL
df_treatment_umi_totals$`4` <- NULL
df_treatment_umi_totals$S <- NULL
head(df_treatment_umi_totals)
df_treatment_umi_totals$group_location <- as.factor(df_treatment_umi_totals$group_location)


if(!require(bayesboot)){
  install.packages("bayesboot")
  library(bayesboot)
}
b <- function(data){
  return(bayesboot(data, mean, R = 100))
}
temp_Control_mouse_Control <- ds(data.frame("b" = b(df_treatment_umi_totals$sample_total[df_treatment_umi_totals$group_location_dpi=="Control_mouse_Control"]), "g" = "Control_mouse_Control"), varname = "V1", groupnames = "g")
temp_Control_PC_Control <- ds(data.frame("b" = b(df_treatment_umi_totals$sample_total[df_treatment_umi_totals$group_location_dpi=="Control_PC_Control"]), "g" = "Control_PC_Control"), varname = "V1", groupnames = "g")
temp_tet_body_14 <- ds(data.frame("b" = b(df_treatment_umi_totals$sample_total[df_treatment_umi_totals$group_location_dpi=="tet_body_14"]), "g" = "tet_body_14"), varname = "V1", groupnames = "g")
temp_tet_body_4 <- ds(data.frame("b" = b(df_treatment_umi_totals$sample_total[df_treatment_umi_totals$group_location_dpi=="tet_body_4"]), "g" = "tet_body_4"), varname = "V1", groupnames = "g")
temp_tet_body_7 <- ds(data.frame("b" = b(df_treatment_umi_totals$sample_total[df_treatment_umi_totals$group_location_dpi=="tet_body_7"]), "g" = "tet_body_7"), varname = "V1", groupnames = "g")
temp_tet_legs_14 <- ds(data.frame("b" = b(df_treatment_umi_totals$sample_total[df_treatment_umi_totals$group_location_dpi=="tet_legs_14"]), "g" = "tet_legs_14"), varname = "V1", groupnames = "g")
temp_tet_legs_7 <- ds(data.frame("b" = b(df_treatment_umi_totals$sample_total[df_treatment_umi_totals$group_location_dpi=="tet_legs_7"]), "g" = "tet_legs_7"), varname = "V1", groupnames = "g")
temp_tet_saliva_14 <- ds(data.frame("b" = b(df_treatment_umi_totals$sample_total[df_treatment_umi_totals$group_location_dpi=="tet_saliva_14"]), "g" = "tet_saliva_14"), varname = "V1", groupnames = "g")
temp_tet_saliva_7 <- ds(data.frame("b" = b(df_treatment_umi_totals$sample_total[df_treatment_umi_totals$group_location_dpi=="tet_saliva_7"]), "g" = "tet_saliva_7"), varname = "V1", groupnames = "g")
temp_wmel_body_14 <- ds(data.frame("b" = b(df_treatment_umi_totals$sample_total[df_treatment_umi_totals$group_location_dpi=="wmel_body_14"]), "g" = "wmel_body_14"), varname = "V1", groupnames = "g")
temp_wmel_body_7 <- ds(data.frame("b" = b(df_treatment_umi_totals$sample_total[df_treatment_umi_totals$group_location_dpi=="wmel_body_7"]), "g" = "wmel_body_7"), varname = "V1", groupnames = "g")
#temp_wmel_legs_14 <- ds(data.frame("b" = b(df_treatment_umi_totals$sample_total[df_treatment_umi_totals$group_location_dpi=="wmel_legs_14"]), "g" = "sample_totals"), varname = "V1", groupnames = "g")
#temp_wmel_legs_7 <- ds(data.frame("b" = b(df_treatment_umi_totals$sample_total[df_treatment_umi_totals$group_location_dpi=="wmel_legs_7"]), "g" = "sample_totals"), varname = "V1", groupnames = "g")

ds_treatment_umi_totals <- rbind(temp_Control_mouse_Control,
                                 temp_Control_PC_Control,
                                 temp_tet_body_14,
                                 temp_tet_body_4,
                                 temp_tet_body_7,
                                 temp_tet_legs_14,
                                 temp_tet_legs_7,
                                 temp_tet_saliva_14,
                                 temp_tet_saliva_7,
                                 temp_wmel_body_14,
                                 temp_wmel_body_7)

#ds_treatment_umi_totals <- as.data.frame(ds(df_treatment_umi_totals, 
#                                            varname = "sample_total", 
#                                            groupnames = c("group_location_dpi")))
#
ds_treatment_umi_totals <- separate(ds_treatment_umi_totals, "g", 
                                    into = c("group", "location", "dpi"), sep = "_")
ds_treatment_umi_totals$group_location <- paste(ds_treatment_umi_totals$group, 
                                                ds_treatment_umi_totals$location, sep = "_")
ds_treatment_umi_totals$group_location_dpi <- paste(ds_treatment_umi_totals$group, 
                                                    ds_treatment_umi_totals$location, 
                                                    ds_treatment_umi_totals$dpi, 
                                                    sep = "_")
ds_treatment_umi_totals$group_location_dpi <- as.factor(ds_treatment_umi_totals$group_location_dpi)
ds_treatment_umi_totals$group_location <- as.factor(ds_treatment_umi_totals$group_location)
ds_treatment_umi_totals$dpi <- factor(ds_treatment_umi_totals$dpi, 
                                      levels = c("4", "7", "14", "Control"))
ds_treatment_umi_totals$group_location <- factor(ds_treatment_umi_totals$group_location, 
                                      levels = c("tet_body", "tet_legs", "tet_saliva", "wmel_body", "Control_mouse", "Control_PC"))
df_treatment_umi_totals$dpi <- factor(df_treatment_umi_totals$dpi, 
                                      levels = c("4", "7", "14", "Control"))

ggplot(data = df_barcode_freqs, aes(x = sample_name, y = freq)) + geom_point()

## sample total by group_location; group = dpi
#plot1 <- ggplot(ds_treatment_umi_totals, aes(color = dpi), group = dpi) + 
#  geom_point(aes(x = group_location, y = V1),
#             position = position_dodge(.75)) + 
#  geom_errorbar(aes(x = group_location, y = V1, 
#                    ymin = V1 - sd, 
#                    ymax = V1 + sd), 
#                width = .15, position = position_dodge(.75)) +
#  scale_y_continuous(limits = c(0, 2000)) + 
#  theme_bw()
#
#plot1.1 <- ggplot(df_treatment_umi_totals, aes(color = dpi), group = dpi) + 
#  geom_violin(aes(x = group_location, y = sample_total),
#              position = position_dodge(.75)) + 
#  geom_point(aes(x = group_location, y = sample_total),
#             position = position_dodge(.75)) + 
#  theme_bw()
#
### sample total by dpi; group = group_location
#plot2 <- ggplot(ds_treatment_umi_totals, aes(color = group_location), group = group_location) + 
#  geom_point(aes(x = dpi, y = sample_total),
#             position = position_dodge(.75)) + 
#  geom_errorbar(aes(x = dpi, y = sample_total, 
#                    ymin = sample_total - sd, 
#                    ymax = sample_total + sd), 
#                width = .15, position = position_dodge(.75)) +
#  theme_bw()
#plot <- plot_grid(plot1, plot2, nrow = 2)

#### group df_barcode_freqs ####
df_barcode_freqs$ID[grepl("7-2-tet-saliva",df_barcode_freqs$sample_name)] <- "Tet-saliva-7dpi_2"
df_barcode_freqs$ID[grepl("7-3-tet-saliva",df_barcode_freqs$sample_name)] <- "Tet-saliva-7dpi_3"
df_barcode_freqs$ID[grepl("14-1-tet-saliva",df_barcode_freqs$sample_name)] <- "Tet-saliva-14dpi_1"
df_barcode_freqs$ID[grepl("14-3-tet-saliva",df_barcode_freqs$sample_name)] <- "Tet-saliva-14dpi_3"
df_barcode_freqs$ID[grepl("7-1-tet-legs",df_barcode_freqs$sample_name)] <- "Tet-legs-7dpi_1"
df_barcode_freqs$ID[grepl("7-2-tet-legs",df_barcode_freqs$sample_name)] <- "Tet-legs-7dpi_2"
df_barcode_freqs$ID[grepl("7-3-tet-legs",df_barcode_freqs$sample_name)] <- "Tet-legs-7dpi_3"
df_barcode_freqs$ID[grepl("14-1-tet-legs",df_barcode_freqs$sample_name)] <- "Tet-legs-14dpi_1"
df_barcode_freqs$ID[grepl("14-2-tet-legs",df_barcode_freqs$sample_name)] <- "Tet-legs-14dpi_2"
df_barcode_freqs$ID[grepl("14-3-tet-legs",df_barcode_freqs$sample_name)] <- "Tet-legs-14dpi_3"
df_barcode_freqs$ID[grepl("4-1-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-4dpi_1"
df_barcode_freqs$ID[grepl("4-2-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-4dpi_2"
df_barcode_freqs$ID[grepl("4-3-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-4dpi_3"
df_barcode_freqs$ID[grepl("7-1-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-7dpi_1"
df_barcode_freqs$ID[grepl("7-2-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-7dpi_2"
df_barcode_freqs$ID[grepl("7-3-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-7dpi_3"
df_barcode_freqs$ID[grepl("dup-14-1-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-14dpi_1_dup"
df_barcode_freqs$ID[grepl("14-1-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-14dpi_1"
df_barcode_freqs$ID[grepl("dup-14-2-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-14dpi_2_dup"
df_barcode_freqs$ID[grepl("14-2-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-14dpi_2"
df_barcode_freqs$ID[grepl("dup-14-3-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-14dpi_3_dup"
df_barcode_freqs$ID[grepl("14-3-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-14dpi_3"
df_barcode_freqs$ID[grepl("7-1-wmel-body",df_barcode_freqs$sample_name)] <- "wmel-body-7dpi_1"
df_barcode_freqs$ID[grepl("7-2-wmel-body",df_barcode_freqs$sample_name)] <- "wmel-body-7dpi_2"
df_barcode_freqs$ID[grepl("7-3-wmel-body",df_barcode_freqs$sample_name)] <- "wmel-body-7dpi_3"
df_barcode_freqs$ID[grepl("14-1-wmel-body",df_barcode_freqs$sample_name)] <- "wmel-body-14dpi_1"
df_barcode_freqs$ID[grepl("14-2-wmel-body",df_barcode_freqs$sample_name)] <- "wmel-body-14dpi_2"
df_barcode_freqs$ID[grepl("14-3-wmel-body",df_barcode_freqs$sample_name)] <- "wmel-body-14dpi_3"
df_barcode_freqs$ID[grepl("7-1-wmel-legs",df_barcode_freqs$sample_name)] <- "wmel-legs-7dpi_1"
df_barcode_freqs$ID[grepl("14-2-wmel-legs",df_barcode_freqs$sample_name)] <- "wmel-legs-14dpi_2"
df_barcode_freqs$ID[grepl("7-2-tet-saliva",df_barcode_freqs$sample_name)] <- "Tet-saliva-7dpi_2"
df_barcode_freqs$ID[grepl("7-3-tet-saliva",df_barcode_freqs$sample_name)] <- "Tet-saliva-7dpi_3"
df_barcode_freqs$ID[grepl("14-1-tet-saliva",df_barcode_freqs$sample_name)] <- "Tet-saliva-14dpi_1"
df_barcode_freqs$ID[grepl("14-3-tet-saliva",df_barcode_freqs$sample_name)] <- "Tet-saliva-14dpi_3"
df_barcode_freqs$ID[grepl("7-1-tet-legs",df_barcode_freqs$sample_name)] <- "Tet-legs-7dpi_1"
df_barcode_freqs$ID[grepl("7-2-tet-legs",df_barcode_freqs$sample_name)] <- "Tet-legs-7dpi_2"
df_barcode_freqs$ID[grepl("7-3-tet-legs",df_barcode_freqs$sample_name)] <- "Tet-legs-7dpi_3"
df_barcode_freqs$ID[grepl("14-1-tet-legs",df_barcode_freqs$sample_name)] <- "Tet-legs-14dpi_1"
df_barcode_freqs$ID[grepl("14-2-tet-legs",df_barcode_freqs$sample_name)] <- "Tet-legs-14dpi_2"
df_barcode_freqs$ID[grepl("14-3-tet-legs",df_barcode_freqs$sample_name)] <- "Tet-legs-14dpi_3"
df_barcode_freqs$ID[grepl("4-1-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-4dpi_1"
df_barcode_freqs$ID[grepl("4-2-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-4dpi_2"
df_barcode_freqs$ID[grepl("4-3-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-4dpi_3"
df_barcode_freqs$ID[grepl("7-1-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-7dpi_1"
df_barcode_freqs$ID[grepl("7-2-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-7dpi_2"
df_barcode_freqs$ID[grepl("7-3-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-7dpi_3"
df_barcode_freqs$ID[grepl("dup-14-1-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-14dpi_1_dup"
df_barcode_freqs$ID[grepl("14-1-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-14dpi_1"
df_barcode_freqs$ID[grepl("dup-14-2-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-14dpi_2_dup"
df_barcode_freqs$ID[grepl("14-2-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-14dpi_2"
df_barcode_freqs$ID[grepl("dup-14-3-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-14dpi_3_dup"
df_barcode_freqs$ID[grepl("14-3-tet-body",df_barcode_freqs$sample_name)] <- "Tet-body-14dpi_3"
df_barcode_freqs$ID[grepl("7-1-wmel-body",df_barcode_freqs$sample_name)] <- "wmel-body-7dpi_1"
df_barcode_freqs$ID[grepl("7-2-wmel-body",df_barcode_freqs$sample_name)] <- "wmel-body-7dpi_2"
df_barcode_freqs$ID[grepl("7-3-wmel-body",df_barcode_freqs$sample_name)] <- "wmel-body-7dpi_3"
df_barcode_freqs$ID[grepl("14-1-wmel-body",df_barcode_freqs$sample_name)] <- "wmel-body-14dpi_1"
df_barcode_freqs$ID[grepl("14-2-wmel-body",df_barcode_freqs$sample_name)] <- "wmel-body-14dpi_2"
df_barcode_freqs$ID[grepl("14-3-wmel-body",df_barcode_freqs$sample_name)] <- "wmel-body-14dpi_3"
df_barcode_freqs$ID[grepl("7-1-wmel-legs",df_barcode_freqs$sample_name)] <- "wmel-legs-7dpi_1"
df_barcode_freqs$ID[grepl("14-2-wmel-legs",df_barcode_freqs$sample_name)] <- "wmel-legs-14dpi_2"
## remove controls
df_barcode_freqs_NC <- df_barcode_freqs[grepl("NC",df_barcode_freqs$sample_name),]
df_barcode_freqs_PC <- df_barcode_freqs[grepl("PC",df_barcode_freqs$sample_name),]
df_barcode_freqs_ZIKV <- df_barcode_freqs[grepl("ZIKV",df_barcode_freqs$sample_name),]
df_barcode_freqs_mouse <- df_barcode_freqs[grepl("mouse",df_barcode_freqs$sample_name),]
df_barcode_freqs <- df_barcode_freqs[!grepl("NC",df_barcode_freqs$sample_name),]
df_barcode_freqs <- df_barcode_freqs[!grepl("PC",df_barcode_freqs$sample_name),]
df_barcode_freqs <- df_barcode_freqs[!grepl("ZIKV",df_barcode_freqs$sample_name),]
df_barcode_freqs <- df_barcode_freqs[!grepl("mouse",df_barcode_freqs$sample_name),]
## separate cols
df_barcode_freqs <- separate(df_barcode_freqs, "ID", c("1","2","3"), sep="-")
df_barcode_freqs <- separate(df_barcode_freqs, "3", c("3","4"), sep="_")
df_barcode_freqs$`3` <- as.integer(gsub("dpi", "", df_barcode_freqs$`3`))

## bring controls back in
df_barcode_freqs_PC$`1` <- "Control"
df_barcode_freqs_mouse$`1` <- "Control"
df_barcode_freqs_PC$`2` <- "PC"
df_barcode_freqs_mouse$`2` <- "mouse"
df_barcode_freqs_PC$`3` <- "Control"
df_barcode_freqs_mouse$`3` <- "Control"
df_barcode_freqs_PC$`4` <- "Control"
df_barcode_freqs_mouse$`4` <- "Control"
df_barcode_freqs_PC$ID <- NULL
df_barcode_freqs_mouse$ID <- NULL

df_barcode_freqs <- rbind(df_barcode_freqs, 
                             df_barcode_freqs_PC,
                             df_barcode_freqs_mouse)

## groups
df_barcode_freqs$group_location <- paste(df_barcode_freqs$`1`, 
                                            df_barcode_freqs$`2`, sep = "_")
df_barcode_freqs$group_location_dpi <- paste(df_barcode_freqs$`1`, 
                                                df_barcode_freqs$`2`, 
                                                df_barcode_freqs$`3`, 
                                                sep = "_")
df_barcode_freqs$group_location <- as.factor(df_barcode_freqs$group_location)
df_barcode_freqs$group_location_dpi <- as.factor(df_barcode_freqs$group_location_dpi)

## order
df_barcode_freqs$`3` <- factor(df_barcode_freqs$`3`, levels = c("4", "7", "14", "Control"))
df_barcode_freqs_max <- df_barcode_freqs %>% group_by(sample_name) %>% summarise(Value = max(freq))
df_barcode_freqs <- merge(df_barcode_freqs, df_barcode_freqs_max, by = "sample_name")
df_barcode_freqs$max <- round(df_barcode_freqs$Value, 2)

df_barcode_freqs$sample_name <- reorder(df_barcode_freqs$sample_name, 
                                        df_barcode_freqs$max)
df_barcode_freqs$sample_name <- factor(df_barcode_freqs$sample_name, 
                                       levels=rev(levels(df_barcode_freqs$sample_name)))
#### plots ####
function_bc_prop <- function(n) {
  plot <- ggplot(df_barcode_freqs[df_barcode_freqs$group_location_dpi==n,], 
                 aes(x = sample_name, y = freq, fill = barcode, group = barcode)) + 
    geom_bar(stat="identity") + labs(x = n, y = "Rel. proportion") +
    geom_text(aes(label=max, y = 1), vjust = -0.5, size = 1.5) + 
    scale_fill_grey() + 
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1.2)) + 
    axis_formatting + legend_formatting + background_formatting + 
    theme(legend.position = "none", axis.text.x = element_blank())
  return(plot)
}
n <- "Control_PC_Control";     p1 <- function_bc_prop(n); p1
n <- "Control_mouse_Control";  p2 <- function_bc_prop(n) 
n <- "Tet_saliva_7";           p4 <- function_bc_prop(n) 
n <- "Tet_saliva_14";          p5 <- function_bc_prop(n) 
n <- "Tet_body_4";             p6 <- function_bc_prop(n) 
n <- "Tet_body_7";             p7 <- function_bc_prop(n) 
n <- "Tet_body_14";            p8 <- function_bc_prop(n) 
n <- "Tet_legs_7";             p9 <- function_bc_prop(n) 
n <- "Tet_legs_14";            p10 <- function_bc_prop(n) 
n <- "wmel_body_7";            p11 <- function_bc_prop(n) 
n <- "wmel_body_14";           p12 <- function_bc_prop(n) 
n <- "wmel_legs_7";            p13 <- function_bc_prop(n) 
n <- "wmel_legs_14";           p14 <- function_bc_prop(n) 

plot_rel_prop <- plot_grid(NULL, p1,  p2,
                           NULL, p4,  p5,
                           p6,   p7,  p8,
                           NULL, p9,  p10,
                           NULL, p11, p12,
                           NULL, p13, p14,
                           ncol = 3)

#### save ####
setwd(dir_save)
#Fig1
ggsave("Fig3_barcodes.pdf", plot,
       width = 5, height = 5, 
       units = "in", dpi = 320)
ggsave("Fig3_barcodes_dot.pdf", plot1.1,
       width = 5, height = 2.5, 
       units = "in", dpi = 320)
ggsave("Fig3_barcodes_relprop.pdf", plot_rel_prop,
       width = 10, height = 5, 
       units = "in", dpi = 320)