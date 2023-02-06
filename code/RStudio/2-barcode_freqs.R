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
dir_15_b <- paste("/Users/rieshunter/Documents/bioinformatics/Wolbachia/data/reads/data/run/15_barcode_analyses", sep="")
setwd(dir_15_b)

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
all_bcs <- table[1]

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
## heatmap
exclude <- c("B1-2-ZIKV-seq-2-2-mouse2-serum_S29", "B1-1-ZIKV-seq-2-2-mouse1-serum_S30", 
             "PC-ZIKV-PR-SetA_S21", "PC-ZIKV-SetC_S21", "PC-ZIKV-PR-SetD_S21", 
             "PC-ZIKV-PR-SetA_S25", "PC-ZIKV-PR-SetB_S43", "PC-ZIKV-PR-SetE_S43", "NC-H2O-SetC_S22")
df_pairwise_euc_filtered <- filter(df_pairwise_euc, sample_1 %!in% exclude)
df_pairwise_euc_filtered <- filter(df_pairwise_euc_filtered, sample_2 %!in% exclude)
test2 <- as.matrix(xtabs(
  eucdistance ~ sample_1 + sample_2, 
  data=df_pairwise_euc_filtered))
heatmap(test2)


df_pairwise_euc_filtered$sample_1[grepl("7-2-tet-saliva",df_pairwise_euc_filtered$sample_1)] <- "Tet-saliva-7dpi_2"
df_pairwise_euc_filtered$sample_1[grepl("7-3-tet-saliva",df_pairwise_euc_filtered$sample_1)] <- "Tet-saliva-7dpi_3"
df_pairwise_euc_filtered$sample_1[grepl("14-1-tet-saliva",df_pairwise_euc_filtered$sample_1)] <- "Tet-saliva-14dpi_1"
df_pairwise_euc_filtered$sample_1[grepl("14-3-tet-saliva",df_pairwise_euc_filtered$sample_1)] <- "Tet-saliva-14dpi_3"
df_pairwise_euc_filtered$sample_1[grepl("7-1-tet-legs",df_pairwise_euc_filtered$sample_1)] <- "Tet-legs-7dpi_1"
df_pairwise_euc_filtered$sample_1[grepl("7-2-tet-legs",df_pairwise_euc_filtered$sample_1)] <- "Tet-legs-7dpi_2"
df_pairwise_euc_filtered$sample_1[grepl("7-3-tet-legs",df_pairwise_euc_filtered$sample_1)] <- "Tet-legs-7dpi_3"
df_pairwise_euc_filtered$sample_1[grepl("14-1-tet-legs",df_pairwise_euc_filtered$sample_1)] <- "Tet-legs-14dpi_1"
df_pairwise_euc_filtered$sample_1[grepl("14-2-tet-legs",df_pairwise_euc_filtered$sample_1)] <- "Tet-legs-14dpi_2"
df_pairwise_euc_filtered$sample_1[grepl("14-3-tet-legs",df_pairwise_euc_filtered$sample_1)] <- "Tet-legs-14dpi_3"
df_pairwise_euc_filtered$sample_1[grepl("4-1-tet-body",df_pairwise_euc_filtered$sample_1)] <- "Tet-body-4dpi_1"
df_pairwise_euc_filtered$sample_1[grepl("4-2-tet-body",df_pairwise_euc_filtered$sample_1)] <- "Tet-body-4dpi_2"
df_pairwise_euc_filtered$sample_1[grepl("4-3-tet-body",df_pairwise_euc_filtered$sample_1)] <- "Tet-body-4dpi_3"
df_pairwise_euc_filtered$sample_1[grepl("7-1-tet-body",df_pairwise_euc_filtered$sample_1)] <- "Tet-body-7dpi_1"
df_pairwise_euc_filtered$sample_1[grepl("7-2-tet-body",df_pairwise_euc_filtered$sample_1)] <- "Tet-body-7dpi_2"
df_pairwise_euc_filtered$sample_1[grepl("7-3-tet-body",df_pairwise_euc_filtered$sample_1)] <- "Tet-body-7dpi_3"
df_pairwise_euc_filtered$sample_1[grepl("dup-14-1-tet-body",df_pairwise_euc_filtered$sample_1)] <- "Tet-body-14dpi_1_dup"
df_pairwise_euc_filtered$sample_1[grepl("14-1-tet-body",df_pairwise_euc_filtered$sample_1)] <- "Tet-body-14dpi_1"
df_pairwise_euc_filtered$sample_1[grepl("dup-14-2-tet-body",df_pairwise_euc_filtered$sample_1)] <- "Tet-body-14dpi_2_dup"
df_pairwise_euc_filtered$sample_1[grepl("14-2-tet-body",df_pairwise_euc_filtered$sample_1)] <- "Tet-body-14dpi_2"
df_pairwise_euc_filtered$sample_1[grepl("dup-14-3-tet-body",df_pairwise_euc_filtered$sample_1)] <- "Tet-body-14dpi_3_dup"
df_pairwise_euc_filtered$sample_1[grepl("14-3-tet-body",df_pairwise_euc_filtered$sample_1)] <- "Tet-body-14dpi_3"
df_pairwise_euc_filtered$sample_1[grepl("7-1-wmel-body",df_pairwise_euc_filtered$sample_1)] <- "wmel-body-7dpi_1"
df_pairwise_euc_filtered$sample_1[grepl("7-2-wmel-body",df_pairwise_euc_filtered$sample_1)] <- "wmel-body-7dpi_2"
df_pairwise_euc_filtered$sample_1[grepl("7-3-wmel-body",df_pairwise_euc_filtered$sample_1)] <- "wmel-body-7dpi_3"
df_pairwise_euc_filtered$sample_1[grepl("14-1-wmel-body",df_pairwise_euc_filtered$sample_1)] <- "wmel-body-14dpi_1"
df_pairwise_euc_filtered$sample_1[grepl("14-2-wmel-body",df_pairwise_euc_filtered$sample_1)] <- "wmel-body-14dpi_2"
df_pairwise_euc_filtered$sample_1[grepl("14-3-wmel-body",df_pairwise_euc_filtered$sample_1)] <- "wmel-body-14dpi_3"
df_pairwise_euc_filtered$sample_1[grepl("7-1-wmel-legs",df_pairwise_euc_filtered$sample_1)] <- "wmel-legs_7dpi_1"
df_pairwise_euc_filtered$sample_1[grepl("14-2-wmel-legs",df_pairwise_euc_filtered$sample_1)] <- "wmel-legs_14dpi_2"

df_pairwise_euc_filtered$sample_2[grepl("7-2-tet-saliva",df_pairwise_euc_filtered$sample_2)] <- "Tet-saliva-7dpi_2"
df_pairwise_euc_filtered$sample_2[grepl("7-3-tet-saliva",df_pairwise_euc_filtered$sample_2)] <- "Tet-saliva-7dpi_3"
df_pairwise_euc_filtered$sample_2[grepl("14-1-tet-saliva",df_pairwise_euc_filtered$sample_2)] <- "Tet-saliva-14dpi_1"
df_pairwise_euc_filtered$sample_2[grepl("14-3-tet-saliva",df_pairwise_euc_filtered$sample_2)] <- "Tet-saliva-14dpi_3"
df_pairwise_euc_filtered$sample_2[grepl("7-1-tet-legs",df_pairwise_euc_filtered$sample_2)] <- "Tet-legs-7dpi_1"
df_pairwise_euc_filtered$sample_2[grepl("7-2-tet-legs",df_pairwise_euc_filtered$sample_2)] <- "Tet-legs-7dpi_2"
df_pairwise_euc_filtered$sample_2[grepl("7-3-tet-legs",df_pairwise_euc_filtered$sample_2)] <- "Tet-legs-7dpi_3"
df_pairwise_euc_filtered$sample_2[grepl("14-1-tet-legs",df_pairwise_euc_filtered$sample_2)] <- "Tet-legs-14dpi_1"
df_pairwise_euc_filtered$sample_2[grepl("14-2-tet-legs",df_pairwise_euc_filtered$sample_2)] <- "Tet-legs-14dpi_2"
df_pairwise_euc_filtered$sample_2[grepl("14-3-tet-legs",df_pairwise_euc_filtered$sample_2)] <- "Tet-legs-14dpi_3"
df_pairwise_euc_filtered$sample_2[grepl("4-1-tet-body",df_pairwise_euc_filtered$sample_2)] <- "Tet-body-4dpi_1"
df_pairwise_euc_filtered$sample_2[grepl("4-2-tet-body",df_pairwise_euc_filtered$sample_2)] <- "Tet-body-4dpi_2"
df_pairwise_euc_filtered$sample_2[grepl("4-3-tet-body",df_pairwise_euc_filtered$sample_2)] <- "Tet-body-4dpi_3"
df_pairwise_euc_filtered$sample_2[grepl("7-1-tet-body",df_pairwise_euc_filtered$sample_2)] <- "Tet-body-7dpi_1"
df_pairwise_euc_filtered$sample_2[grepl("7-2-tet-body",df_pairwise_euc_filtered$sample_2)] <- "Tet-body-7dpi_2"
df_pairwise_euc_filtered$sample_2[grepl("7-3-tet-body",df_pairwise_euc_filtered$sample_2)] <- "Tet-body-7dpi_3"
df_pairwise_euc_filtered$sample_2[grepl("dup-14-1-tet-body",df_pairwise_euc_filtered$sample_2)] <- "Tet-body-14dpi_1_dup"
df_pairwise_euc_filtered$sample_2[grepl("14-1-tet-body",df_pairwise_euc_filtered$sample_2)] <- "Tet-body-14dpi_1"
df_pairwise_euc_filtered$sample_2[grepl("dup-14-2-tet-body",df_pairwise_euc_filtered$sample_2)] <- "Tet-body-14dpi_2_dup"
df_pairwise_euc_filtered$sample_2[grepl("14-2-tet-body",df_pairwise_euc_filtered$sample_2)] <- "Tet-body-14dpi_2"
df_pairwise_euc_filtered$sample_2[grepl("dup-14-3-tet-body",df_pairwise_euc_filtered$sample_2)] <- "Tet-body-14dpi_3_dup"
df_pairwise_euc_filtered$sample_2[grepl("14-3-tet-body",df_pairwise_euc_filtered$sample_2)] <- "Tet-body-14dpi_3"
df_pairwise_euc_filtered$sample_2[grepl("7-1-wmel-body",df_pairwise_euc_filtered$sample_2)] <- "wmel-body-7dpi_1"
df_pairwise_euc_filtered$sample_2[grepl("7-2-wmel-body",df_pairwise_euc_filtered$sample_2)] <- "wmel-body-7dpi_2"
df_pairwise_euc_filtered$sample_2[grepl("7-3-wmel-body",df_pairwise_euc_filtered$sample_2)] <- "wmel-body-7dpi_3"
df_pairwise_euc_filtered$sample_2[grepl("14-1-wmel-body",df_pairwise_euc_filtered$sample_2)] <- "wmel-body-14dpi_1"
df_pairwise_euc_filtered$sample_2[grepl("14-2-wmel-body",df_pairwise_euc_filtered$sample_2)] <- "wmel-body-14dpi_2"
df_pairwise_euc_filtered$sample_2[grepl("14-3-wmel-body",df_pairwise_euc_filtered$sample_2)] <- "wmel-body-14dpi_3"
df_pairwise_euc_filtered$sample_2[grepl("7-1-wmel-legs",df_pairwise_euc_filtered$sample_2)] <- "wmel-legs_7dpi_1"
df_pairwise_euc_filtered$sample_2[grepl("14-2-wmel-legs",df_pairwise_euc_filtered$sample_2)] <- "wmel-legs_14dpi_2"

df_pairwise_euc_filtered %>% 
  select(sample_1, sample_2, eucdistance) %>% 
  group_by(sample_1, sample_2) %>% 
  summarize(mean.euc=mean(eucdistance)) -> df_pairwise_euc_filtered
df_pairwise_euc_filtered %>% 
  ggplot(data=., mapping=aes(x=sample_1, y=sample_2, fill=mean.euc)) + 
  geom_tile() +
  theme(axis.text.x=element_text(color=rainbow(ncol(df_pairwise_euc_filtered))), 
        axis.text.y=element_text(color=rainbow(ncol(df_pairwise_euc_filtered))))


#### Barcodes_overtime_COMPOSITE.R ####
# R script to plot composite barcode trajectories for each passage series.
