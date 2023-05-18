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
if(!require(bayesboot)){
  install.packages("bayesboot")
  library(bayesboot)
}
if(!require(BSDA)){
  install.packages("BSDA")
  library(BSDA)
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

## bayesian bootstrap
b <- function(data){
  return(bayesboot(data, mean, R = 10000))
}

## Adds column ID to dataframe x, given x$sample_name
add_ID_from_samplename <- function(x) {
  x$ID[grepl("7-2-tet-saliva",x$sample_name)] <- "Tet-saliva-7dpi_2"
  x$ID[grepl("7-3-tet-saliva",x$sample_name)] <- "Tet-saliva-7dpi_3"
  x$ID[grepl("14-1-tet-saliva",x$sample_name)] <- "Tet-saliva-14dpi_1"
  x$ID[grepl("14-3-tet-saliva",x$sample_name)] <- "Tet-saliva-14dpi_3"
  x$ID[grepl("7-1-tet-legs",x$sample_name)] <- "Tet-legs-7dpi_1"
  x$ID[grepl("7-2-tet-legs",x$sample_name)] <- "Tet-legs-7dpi_2"
  x$ID[grepl("7-3-tet-legs",x$sample_name)] <- "Tet-legs-7dpi_3"
  x$ID[grepl("14-1-tet-legs",x$sample_name)] <- "Tet-legs-14dpi_1"
  x$ID[grepl("14-2-tet-legs",x$sample_name)] <- "Tet-legs-14dpi_2"
  x$ID[grepl("14-3-tet-legs",x$sample_name)] <- "Tet-legs-14dpi_3"
  x$ID[grepl("4-1-tet-body",x$sample_name)] <- "Tet-body-4dpi_1"
  x$ID[grepl("4-2-tet-body",x$sample_name)] <- "Tet-body-4dpi_2"
  x$ID[grepl("4-3-tet-body",x$sample_name)] <- "Tet-body-4dpi_3"
  x$ID[grepl("7-1-tet-body",x$sample_name)] <- "Tet-body-7dpi_1"
  x$ID[grepl("7-2-tet-body",x$sample_name)] <- "Tet-body-7dpi_2"
  x$ID[grepl("7-3-tet-body",x$sample_name)] <- "Tet-body-7dpi_3"
  x$ID[grepl("dup-14-1-tet-body",x$sample_name)] <- "Tet-body-14dpi_1_dup"
  x$ID[grepl("14-1-tet-body",x$sample_name)] <- "Tet-body-14dpi_1"
  x$ID[grepl("dup-14-2-tet-body",x$sample_name)] <- "Tet-body-14dpi_2_dup"
  x$ID[grepl("14-2-tet-body",x$sample_name)] <- "Tet-body-14dpi_2"
  x$ID[grepl("dup-14-3-tet-body",x$sample_name)] <- "Tet-body-14dpi_3_dup"
  x$ID[grepl("14-3-tet-body",x$sample_name)] <- "Tet-body-14dpi_3"
  x$ID[grepl("7-1-wmel-body",x$sample_name)] <- "wmel-body-7dpi_1"
  x$ID[grepl("7-2-wmel-body",x$sample_name)] <- "wmel-body-7dpi_2"
  x$ID[grepl("7-3-wmel-body",x$sample_name)] <- "wmel-body-7dpi_3"
  x$ID[grepl("14-1-wmel-body",x$sample_name)] <- "wmel-body-14dpi_1"
  x$ID[grepl("14-2-wmel-body",x$sample_name)] <- "wmel-body-14dpi_2"
  x$ID[grepl("14-3-wmel-body",x$sample_name)] <- "wmel-body-14dpi_3"
  x$ID[grepl("7-1-wmel-legs",x$sample_name)] <- "wmel-legs-7dpi_1"
  x$ID[grepl("14-2-wmel-legs",x$sample_name)] <- "wmel-legs-14dpi_2"
  x$ID[grepl("7-2-tet-saliva",x$sample_name)] <- "Tet-saliva-7dpi_2"
  x$ID[grepl("7-3-tet-saliva",x$sample_name)] <- "Tet-saliva-7dpi_3"
  x$ID[grepl("14-1-tet-saliva",x$sample_name)] <- "Tet-saliva-14dpi_1"
  x$ID[grepl("14-3-tet-saliva",x$sample_name)] <- "Tet-saliva-14dpi_3"
  x$ID[grepl("7-1-tet-legs",x$sample_name)] <- "Tet-legs-7dpi_1"
  x$ID[grepl("7-2-tet-legs",x$sample_name)] <- "Tet-legs-7dpi_2"
  x$ID[grepl("7-3-tet-legs",x$sample_name)] <- "Tet-legs-7dpi_3"
  x$ID[grepl("14-1-tet-legs",x$sample_name)] <- "Tet-legs-14dpi_1"
  x$ID[grepl("14-2-tet-legs",x$sample_name)] <- "Tet-legs-14dpi_2"
  x$ID[grepl("14-3-tet-legs",x$sample_name)] <- "Tet-legs-14dpi_3"
  x$ID[grepl("4-1-tet-body",x$sample_name)] <- "Tet-body-4dpi_1"
  x$ID[grepl("4-2-tet-body",x$sample_name)] <- "Tet-body-4dpi_2"
  x$ID[grepl("4-3-tet-body",x$sample_name)] <- "Tet-body-4dpi_3"
  x$ID[grepl("7-1-tet-body",x$sample_name)] <- "Tet-body-7dpi_1"
  x$ID[grepl("7-2-tet-body",x$sample_name)] <- "Tet-body-7dpi_2"
  x$ID[grepl("7-3-tet-body",x$sample_name)] <- "Tet-body-7dpi_3"
  x$ID[grepl("dup-14-1-tet-body",x$sample_name)] <- "Tet-body-14dpi_1_dup"
  x$ID[grepl("14-1-tet-body",x$sample_name)] <- "Tet-body-14dpi_1"
  x$ID[grepl("dup-14-2-tet-body",x$sample_name)] <- "Tet-body-14dpi_2_dup"
  x$ID[grepl("14-2-tet-body",x$sample_name)] <- "Tet-body-14dpi_2"
  x$ID[grepl("dup-14-3-tet-body",x$sample_name)] <- "Tet-body-14dpi_3_dup"
  x$ID[grepl("14-3-tet-body",x$sample_name)] <- "Tet-body-14dpi_3"
  x$ID[grepl("7-1-wmel-body",x$sample_name)] <- "wmel-body-7dpi_1"
  x$ID[grepl("7-2-wmel-body",x$sample_name)] <- "wmel-body-7dpi_2"
  x$ID[grepl("7-3-wmel-body",x$sample_name)] <- "wmel-body-7dpi_3"
  x$ID[grepl("14-1-wmel-body",x$sample_name)] <- "wmel-body-14dpi_1"
  x$ID[grepl("14-2-wmel-body",x$sample_name)] <- "wmel-body-14dpi_2"
  x$ID[grepl("14-3-wmel-body",x$sample_name)] <- "wmel-body-14dpi_3"
  x$ID[grepl("7-1-wmel-legs",x$sample_name)] <- "wmel-legs-7dpi_1"
  x$ID[grepl("14-2-wmel-legs",x$sample_name)] <- "wmel-legs-14dpi_2"
  return(x)
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
                  "Missense" = "#7AD9C2")
# mutation types without stop
palette_muts_NS_S <- c("Nonsynonymous" = "#FF7F20",
                       "Synonymous" = "#4F7899")

# group_location
palette_group_location <- c("tet_body" = "#1E88E5",
                            "tet_legs" = "#FFC107",
                            "tet_saliva" = "#7A4128",
                            "wmel_body" = "#004D40",
                            "wmel_legs" = "#5F794E",
                            "Control_mouse" = "#9D7DDC",
                            "Control_PC" = "#9D7DDC")

palette_dpi <- c("4" = "#DCCE42",
                 "7" = "#157050",
                 "14" = "#D06941",
                 "Control" = "#70688C")

#### Import and clean ####
dir_15_b <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/Wolbachia/data/reads/data/run/15_barcode_analyses", sep="")
dir_save <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/Wolbachia/figs", sep="")
setwd(dir_15_b); dir()

#### Derive df_barcode_freqs and df_bc_totals ####
#### Original script: Mean_bc_frequencies.R
### Calculate species richness and barcode frequency
dir_kmercounts <- dir(pattern="counts_BConly_*")
names_trunc <- gsub(".tsv", "", dir_kmercounts)
names_trunc <- gsub("counts_BConly_", "", names_trunc)
## lists
n <- length(dir_kmercounts)
list <- vector("list", n)
list2 <- vector("list", n)
## data frame for totals
sample_name <- rep("foo", n)
sample_total <- rep(0, n)
species_richness <- rep(0, n)
df_bc_totals <- data.frame(sample_name, sample_total, species_richness)
## for-loop
for (i in 1:n) {
  list[[i]] <- as.data.frame(
    read_tsv(dir_kmercounts[[i]], show_col_types = F))
  list[[i]] <- list[[i]][c(1,2)]
  colnames(list[[i]]) <- c("count", "barcode")
  list2[[i]]$sample_name <- names_trunc[i]
  ifelse(length(list[[i]]$count)>0, 
         list2[[i]]$sample_total <- sum(list[[i]]$count),
         list2[[i]]$sample_total <- NA)
  ifelse(length(list[[i]]$count)>0, 
         list2[[i]]$species_richness <- length(list[[i]]$barcode[list[[i]]$count>2]),
         list2[[i]]$species_richness <- NA)
  names(list) <- names_trunc
  names(list2) <- names_trunc
  df_bc_totals$sample_name[i]  <- list2[[i]]$sample_name
  df_bc_totals$sample_total[i] <- list2[[i]]$sample_total
  df_bc_totals$species_richness[i] <- list2[[i]]$species_richness
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
df_barcode_freqs <- Reduce(full_join, list)
head(df_bc_totals)


#### Assign groups to sample_name ####
#### Original script: Barcodes_overtime_COMPOSITE.R
### Tet
df_tet_bc_totals <- df_bc_totals[grepl("tet",df_bc_totals$sample_name),]
df_tet_bc_totals <- separate(df_tet_bc_totals, "sample_name", c("sample_name","S"), sep="_S")
df_tet_bc_totals <- separate(df_tet_bc_totals, "sample_name", c("sample_name","location"), sep="-tet-")
df_tet_bc_totals_dups <- df_tet_bc_totals[grepl("dup",df_tet_bc_totals$sample_name),]
df_tet_bc_totals_dups <- separate(df_tet_bc_totals_dups, "sample_name", c("1","2","dup","3","4"), sep="-")
df_tet_bc_totals <- df_tet_bc_totals[!grepl("dup",df_tet_bc_totals$sample_name),]
df_tet_bc_totals <- separate(df_tet_bc_totals, "sample_name", c("1","2","3", "4"), sep="-")
df_tet_bc_totals$dup <- "not_dup"
df_tet_bc_totals <- rbind(df_tet_bc_totals, df_tet_bc_totals_dups)
df_tet_bc_totals$group <- "tet"
head(df_tet_bc_totals)
### Wmel
df_wmel_bc_totals <- df_bc_totals[grepl("wmel",df_bc_totals$sample_name),]
df_wmel_bc_totals <- separate(df_wmel_bc_totals, "sample_name", c("sample_name","S"), sep="_S")
df_wmel_bc_totals <- separate(df_wmel_bc_totals, "sample_name", c("sample_name","location"), sep="-wmel-")
df_wmel_bc_totals_dups <- df_wmel_bc_totals[grepl("dup",df_wmel_bc_totals$sample_name),]
df_wmel_bc_totals_dups <- separate(df_wmel_bc_totals_dups, "sample_name", c("1","2","dup","3","4"), sep="-")
df_wmel_bc_totals <- df_wmel_bc_totals[!grepl("dup",df_wmel_bc_totals$sample_name),]
df_wmel_bc_totals <- separate(df_wmel_bc_totals, "sample_name", c("1","2","3", "4"), sep="-")
df_wmel_bc_totals$dup <- "not_dup"
df_wmel_bc_totals <- rbind(df_wmel_bc_totals, df_wmel_bc_totals_dups)
df_wmel_bc_totals$group <- "wmel"
head(df_wmel_bc_totals)
### Control
df_controls_bc_totals <- df_bc_totals[grepl("ZIKV",df_bc_totals$sample_name),]
df_controls_bc_totals$sample_name[grepl("mouse",df_controls_bc_totals$sample_name)] <- "Control_mouse"
df_controls_bc_totals$sample_name[grepl("PC",df_controls_bc_totals$sample_name)] <- "Control_PC"
df_controls_bc_totals$S <- ""
df_controls_bc_totals <- separate(df_controls_bc_totals, "sample_name", c("sample_name","location"), sep="_")
df_controls_bc_totals$dup <- "not_dup"
df_controls_bc_totals$group <- "Control"
df_controls_bc_totals$`1` <- "Control"
df_controls_bc_totals$`2` <- df_controls_bc_totals$location
df_controls_bc_totals$`3` <- "Control"
df_controls_bc_totals$`4` <- "Control"
df_controls_bc_totals$sample_name <- NULL
head(df_controls_bc_totals)
### Bind together
df_treatment_bc_totals <- rbind(df_tet_bc_totals, df_wmel_bc_totals, df_controls_bc_totals)
df_tet_bc_totals <- NULL
df_wmel_bc_totals <- NULL
df_controls_bc_totals <- NULL
df_tet_bc_totals_dups <- NULL
df_wmel_bc_totals_dups <- NULL
df_treatment_bc_totals$ID <- paste(df_treatment_bc_totals$`1`,
                                   df_treatment_bc_totals$`2`,
                                   sep = "_")
df_treatment_bc_totals$rep <- df_treatment_bc_totals$`4`
df_treatment_bc_totals$dpi <- df_treatment_bc_totals$`3`
df_treatment_bc_totals$group_location <- paste(df_treatment_bc_totals$group, 
                                                df_treatment_bc_totals$location, sep = "_")
df_treatment_bc_totals$group_location_dpi <- paste(df_treatment_bc_totals$group, 
                                                    df_treatment_bc_totals$location, 
                                                    df_treatment_bc_totals$dpi, 
                                                    sep = "_")
df_treatment_bc_totals$`1` <- NULL
df_treatment_bc_totals$`2` <- NULL
df_treatment_bc_totals$`3` <- NULL
df_treatment_bc_totals$`4` <- NULL
df_treatment_bc_totals$S <- NULL
head(df_treatment_bc_totals)
df_treatment_bc_totals$group_location <- as.factor(df_treatment_bc_totals$group_location)
df_treatment_bc_totals <- df_treatment_bc_totals[!is.na(df_treatment_bc_totals$species_richness),]

## Bayesian bootstrapping of species_richness and mean/sd/se of distribution
temp_Control_mouse_Control <- ds(data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="Control_mouse_Control"]), "g" = "Control_mouse_Control"), varname = "V1", groupnames = "g")
temp_Control_PC_Control <-    ds(data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="Control_PC_Control"]), "g" = "Control_PC_Control"), varname = "V1", groupnames = "g")
temp_tet_body_4 <-            ds(data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="tet_body_4"]), "g" = "tet_body_4"), varname = "V1", groupnames = "g")
temp_tet_body_7 <-            ds(data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="tet_body_7"]), "g" = "tet_body_7"), varname = "V1", groupnames = "g")
temp_tet_body_14 <-           ds(data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="tet_body_14"]), "g" = "tet_body_14"), varname = "V1", groupnames = "g")
temp_tet_legs_7 <-            ds(data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="tet_legs_7"]), "g" = "tet_legs_7"), varname = "V1", groupnames = "g")
temp_tet_legs_14 <-           ds(data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="tet_legs_14"]), "g" = "tet_legs_14"), varname = "V1", groupnames = "g")
temp_tet_saliva_7 <-          ds(data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="tet_saliva_7"]), "g" = "tet_saliva_7"), varname = "V1", groupnames = "g")
temp_tet_saliva_14 <-         ds(data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="tet_saliva_14"]), "g" = "tet_saliva_14"), varname = "V1", groupnames = "g")
temp_wmel_body_7 <-           ds(data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="wmel_body_7"]), "g" = "wmel_body_7"), varname = "V1", groupnames = "g")
temp_wmel_body_14 <-          ds(data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="wmel_body_14"]), "g" = "wmel_body_14"), varname = "V1", groupnames = "g")
## single sample, so no sd/se
temp_wmel_legs_14 <-          data.frame("g" = df_treatment_bc_totals$group_location_dpi[df_treatment_bc_totals$group_location_dpi=="wmel_legs_14"],
                                         "V1"= df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="wmel_legs_14"],
                                         "sd"= 0,
                                         "se"= 0)
temp_wmel_legs_7 <-           data.frame("g" = df_treatment_bc_totals$group_location_dpi[df_treatment_bc_totals$group_location_dpi=="wmel_legs_7"],
                                         "V1"= df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="wmel_legs_7"],
                                         "sd"= 0,
                                         "se"= 0)
#temp_wmel_legs_14 <-          ds(data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="wmel_legs_14"]), "g" = "sample_totals"), varname = "V1", groupnames = "g")
#temp_wmel_legs_7 <-           ds(data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="wmel_legs_7"]), "g" = "sample_totals"), varname = "V1", groupnames = "g")

### Bind to one df
ds_treatment_bc_totals <- rbind(temp_Control_mouse_Control,
                                temp_Control_PC_Control,
                                temp_tet_body_4,
                                temp_tet_body_7,
                                temp_tet_body_14,
                                temp_tet_legs_7,
                                temp_tet_legs_14,
                                temp_tet_saliva_7,
                                temp_tet_saliva_14,
                                temp_wmel_body_7,
                                temp_wmel_body_14,
                                temp_wmel_legs_7,
                                temp_wmel_legs_14)

ds_treatment_bc_totals <- separate(ds_treatment_bc_totals, "g", 
                                    into = c("group", "location", "dpi"), sep = "_")
ds_treatment_bc_totals$group_location <- paste(ds_treatment_bc_totals$group, 
                                                ds_treatment_bc_totals$location, sep = "_")
ds_treatment_bc_totals$group_location_dpi <- paste(ds_treatment_bc_totals$group, 
                                                    ds_treatment_bc_totals$location, 
                                                    ds_treatment_bc_totals$dpi, 
                                                    sep = "_")
ds_treatment_bc_totals$group_location_dpi <- as.factor(ds_treatment_bc_totals$group_location_dpi)
ds_treatment_bc_totals$group_location <- as.factor(ds_treatment_bc_totals$group_location)
ds_treatment_bc_totals$dpi <- factor(ds_treatment_bc_totals$dpi, 
                                      levels = c("4", "7", "14", "Control"))
ds_treatment_bc_totals$group_location <- factor(ds_treatment_bc_totals$group_location, 
                                      levels = c("tet_body", "tet_legs", "tet_saliva", "wmel_body", "wmel_legs", "Control_mouse", "Control_PC"))
df_treatment_bc_totals$dpi <- factor(df_treatment_bc_totals$dpi, 
                                      levels = c("4", "7", "14", "Control"))

#### group df_barcode_freqs ####
df_barcode_freqs <- add_ID_from_samplename(df_barcode_freqs)
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

## calculate max bc freq for each sample
df_barcode_freqs_max <- aggregate(freq ~ sample_name, data = df_barcode_freqs, max)
df_barcode_freqs_max$Max <- df_barcode_freqs_max$freq
df_barcode_freqs_max$freq <- NULL
## calculate avg bc freq for each sample
df_barcode_freqs_avg <- aggregate(freq ~ sample_name, data = df_barcode_freqs, mean)
df_barcode_freqs_avg$Mean <- df_barcode_freqs_avg$freq
df_barcode_freqs_avg$freq <- NULL
## calculate median bc freq for each sample
df_barcode_freqs_med <- aggregate(freq ~ sample_name, data = df_barcode_freqs, median)
df_barcode_freqs_med$Median <- df_barcode_freqs_med$freq
df_barcode_freqs_med$freq <- NULL
## merge all calculations into df_barcode_calcs
df_barcode_calcs <- merge(df_barcode_freqs_max, df_barcode_freqs_avg, by = "sample_name")
df_barcode_calcs <- merge(df_barcode_calcs, df_barcode_freqs_med, by = "sample_name")
## assign groups to sample_names
df_barcode_calcs <- add_ID_from_samplename(df_barcode_calcs)
## remove controls
df_barcode_calcs_NC <- df_barcode_calcs[grepl("NC",df_barcode_calcs$sample_name),]
df_barcode_calcs_PC <- df_barcode_calcs[grepl("PC",df_barcode_calcs$sample_name),]
df_barcode_calcs_ZIKV <- df_barcode_calcs[grepl("ZIKV",df_barcode_calcs$sample_name),]
df_barcode_calcs_mouse <- df_barcode_calcs[grepl("mouse",df_barcode_calcs$sample_name),]
df_barcode_calcs <- df_barcode_calcs[!grepl("NC",df_barcode_calcs$sample_name),]
df_barcode_calcs <- df_barcode_calcs[!grepl("PC",df_barcode_calcs$sample_name),]
df_barcode_calcs <- df_barcode_calcs[!grepl("ZIKV",df_barcode_calcs$sample_name),]
df_barcode_calcs <- df_barcode_calcs[!grepl("mouse",df_barcode_calcs$sample_name),]
## separate cols
df_barcode_calcs <- separate(df_barcode_calcs, "ID", c("1","2","3"), sep="-")
df_barcode_calcs <- separate(df_barcode_calcs, "3", c("3","4"), sep="_")
df_barcode_calcs$`3` <- as.integer(gsub("dpi", "", df_barcode_calcs$`3`))
## bring controls back in
df_barcode_calcs_PC$`1` <- "Control"
df_barcode_calcs_mouse$`1` <- "Control"
df_barcode_calcs_PC$`2` <- "PC"
df_barcode_calcs_mouse$`2` <- "mouse"
df_barcode_calcs_PC$`3` <- "Control"
df_barcode_calcs_mouse$`3` <- "Control"
df_barcode_calcs_PC$`4` <- "Control"
df_barcode_calcs_mouse$`4` <- "Control"
df_barcode_calcs_PC$ID <- NULL
df_barcode_calcs_mouse$ID <- NULL
df_barcode_calcs <- rbind(df_barcode_calcs, 
                          df_barcode_calcs_PC,
                          df_barcode_calcs_mouse)
## groups
df_barcode_calcs$group_location <- paste(df_barcode_calcs$`1`, 
                                         df_barcode_calcs$`2`, sep = "_")
df_barcode_calcs$group_location_dpi <- paste(df_barcode_calcs$`1`, 
                                             df_barcode_calcs$`2`, 
                                             df_barcode_calcs$`3`, 
                                             sep = "_")
df_barcode_calcs$group_location <- as.factor(df_barcode_calcs$group_location)
df_barcode_calcs$group_location_dpi <- as.factor(df_barcode_calcs$group_location_dpi)


## calculate max bc freq for each sample
df_barcode_freqs_max <- aggregate(freq ~ sample_name, data = df_barcode_freqs, max)
df_barcode_freqs_max$Value <- df_barcode_freqs_max$freq
df_barcode_freqs_max$freq <- NULL
# merge to main df
df_barcode_freqs <- merge(df_barcode_freqs, df_barcode_freqs_max, by = "sample_name")
df_barcode_freqs$max <- round(df_barcode_freqs$Value, 2)

df_barcode_freqs$sample_name <- reorder(df_barcode_freqs$sample_name, 
                                        df_barcode_freqs$max)
df_barcode_freqs$sample_name <- factor(df_barcode_freqs$sample_name, 
                                       levels=rev(levels(df_barcode_freqs$sample_name)))
#### Plots - Barcode species richness ####
plot1 <- ggplot(ds_treatment_bc_totals, aes(color = dpi), group = dpi) + 
  geom_point(aes(x = group_location, y = V1),
             position = position_dodge(.75)) + 
  geom_errorbar(data = ds_treatment_bc_totals, 
                aes(x = group_location, 
                    ymin = V1 - sd, 
                    ymax = V1 + sd), 
                width = .15, position = position_dodge(.75)) +
  #geom_text(aes(x = group_location, label = round(V1,0), y = V1), 
  #          vjust = -5, size = 1.5,
  #          position = position_dodge(.75)) +
  scale_color_manual(values = palette_dpi) + 
  scale_y_log10(limits = c(1,1000)) + 
  labs(x = "", y = "Barcode species richness") + 
  annotation_logticks(sides="l") + 
  theme_bw()

plot2 <- ggplot(ds_treatment_bc_totals, aes(color = group_location), group = group_location) + 
  geom_point(aes(x = dpi, y = V1),
             position = position_dodge(.75)) + 
  geom_errorbar(aes(x = dpi, y = V1, 
                    ymin = V1 - sd, 
                    ymax = V1 + sd), 
                width = .15, position = position_dodge(.75)) +
  scale_color_manual(values = palette_group_location) + 
  scale_y_log10(limits = c(1,1000)) + 
  labs(x = "", y = "Barcode species richness") + 
  annotation_logticks(sides="l") + 
  theme_bw()


#### Plots - Barcode frequency ####
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

df_barcode_calcs$gl <- factor(paste(df_barcode_calcs$`1`, df_barcode_calcs$`2`, sep = "_"), 
                        levels = c("Tet_body", "Tet_legs", "Tet_saliva", "wmel_body", "wmel_legs", "Control_mouse", "Control_PC"))
df_barcode_calcs$dpi <- factor(df_barcode_calcs$`3`, levels = c("4", "7", "14", "Control"))

plot_bc_calcs_max <- ggplot(df_barcode_calcs, aes(color = dpi), group = dpi) +
  geom_point(aes(x = gl, y = Max),
             position = position_dodge(.75)) + 
  geom_boxplot(aes(x = gl, y = Max),
               position = position_dodge(.75)) + 
  labs(x = "Group_location") + 
  axis_formatting + legend_formatting + background_formatting + 
  theme(axis.title.x = element_blank(), legend.position = "none")

plot_bc_calcs_mean <- ggplot(df_barcode_calcs, aes(color = dpi), group = dpi) +
  geom_point(aes(x = gl, y = Mean),
             position = position_dodge(.75)) + 
  geom_boxplot(aes(x = gl, y = Mean),
               position = position_dodge(.75)) + 
  labs(x = "Group_location") + 
  axis_formatting + legend_formatting + background_formatting + 
  theme(axis.title.x = element_blank(), legend.position = "none")

plot_bc_calcs_median <- ggplot(df_barcode_calcs, aes(color = dpi), group = dpi) +
  geom_point(aes(x = gl, y = Median),
             position = position_dodge(.75)) + 
  geom_boxplot(aes(x = gl, y = Median),
               position = position_dodge(.75)) + 
  labs(x = "Group_location") + 
  axis_formatting + legend_formatting + background_formatting + 
  theme(axis.title.x = element_blank(),
        legend.position = "bottom", 
        legend.key = element_blank(), 
        legend.title = element_blank())

plot_bc_stats <- plot_grid(plot_bc_calcs_max,
                           plot_bc_calcs_mean,
                           plot_bc_calcs_median,
                           nrow = 3, rel_heights = c(1,1,1.3))


### bc freqs hist
#df_barcode_freqs$gl <- factor(paste(df_barcode_freqs$`1`, df_barcode_freqs$`2`, sep = "_"), 
#                              levels = c("Tet_body", "Tet_legs", "Tet_saliva", "wmel_body", "wmel_legs", "Control_mouse", "Control_PC"))
#df_barcode_freqs$dpi <- factor(df_barcode_freqs$`3`, levels = c("4", "7", "14", "Control"))
#
#i=paste("Tet_body_14")
#plot_bc_spec <- function(i) {
#  plot <- ggplot(df_barcode_freqs[df_barcode_freqs$freq>0 & 
#                                    df_barcode_freqs$group_location_dpi==i,], 
#                 aes(x = freq, group = sample_name)) + 
#    geom_histogram(bins = 50, fill = "blue") + 
#    labs(y = "Number of observations", x = "Frequency of barcode", title = i) + 
#    axis_formatting + legend_formatting + background_formatting + 
#    theme(legend.position = "none") + facet_grid(rows = vars(sample_name))
#  return(plot)
#}
#
#Tet_body <- plot_grid(plot_bc_spec("Tet_body_4"),
#                      plot_bc_spec("Tet_body_7"),
#                      plot_bc_spec("Tet_body_14"),
#                      ncol = 3)
#Tet_legs <- plot_grid(NULL, 
#                      plot_bc_spec("Tet_legs_7"),
#                      plot_bc_spec("Tet_legs_14"),
#                      ncol = 3)
#Tet_saliva <- plot_grid(NULL, plot_bc_spec("Tet_saliva_7"),
#                        plot_bc_spec("Tet_saliva_14"),
#                        ncol = 3)
#wmel_body <- plot_grid(NULL, plot_bc_spec("wmel_body_7"),
#                       plot_bc_spec("wmel_body_14"),
#                       ncol = 3)
#wmel_legs <- plot_grid(NULL, plot_bc_spec("wmel_legs_7"),
#                       plot_bc_spec("wmel_legs_14"),
#                       ncol = 3)
#control <- plot_grid(NULL, NULL, plot_bc_spec("Control_mouse_Control"), ncol = 3)
#
#plot_bc_freq_spec <- plot_grid(Tet_body,
#                               Tet_legs,
#                               Tet_saliva,
#                               wmel_body,
#                               wmel_legs,
#                               control,
#                               ncol = 1,
#                               rel_heights = c(12,8,3,6,2,1))





#### save ####
setwd(dir_save); dir()

temp_Control_mouse_Control <- data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="Control_mouse_Control"]), "sd" = ds_treatment_bc_totals$sd[ds_treatment_bc_totals$group_location_dpi=="Control_mouse_Control"])
temp_Control_PC_Control <-    data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="Control_PC_Control"]), "sd" = ds_treatment_bc_totals$sd[ds_treatment_bc_totals$group_location_dpi=="Control_PC_Control"])
temp_tet_body_4 <-            data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="tet_body_4"]), "sd" = ds_treatment_bc_totals$sd[ds_treatment_bc_totals$group_location_dpi=="tet_body_4"])
temp_tet_body_7 <-            data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="tet_body_7"]), "sd" = ds_treatment_bc_totals$sd[ds_treatment_bc_totals$group_location_dpi=="tet_body_7"])
temp_tet_body_14 <-           data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="tet_body_14"]), "sd" = ds_treatment_bc_totals$sd[ds_treatment_bc_totals$group_location_dpi=="tet_body_14"])
temp_tet_legs_7 <-            data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="tet_legs_7"]), "sd" = ds_treatment_bc_totals$sd[ds_treatment_bc_totals$group_location_dpi=="tet_legs_7"])
temp_tet_legs_14 <-           data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="tet_legs_14"]), "sd" = ds_treatment_bc_totals$sd[ds_treatment_bc_totals$group_location_dpi=="tet_legs_14"])
temp_tet_saliva_7 <-          data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="tet_saliva_7"]), "sd" = ds_treatment_bc_totals$sd[ds_treatment_bc_totals$group_location_dpi=="tet_saliva_7"])
temp_tet_saliva_14 <-         data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="tet_saliva_14"]), "sd" = ds_treatment_bc_totals$sd[ds_treatment_bc_totals$group_location_dpi=="tet_saliva_14"])
temp_wmel_body_7 <-           data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="wmel_body_7"]), "sd" = ds_treatment_bc_totals$sd[ds_treatment_bc_totals$group_location_dpi=="wmel_body_7"])
temp_wmel_body_14 <-          data.frame("b" = b(df_treatment_bc_totals$species_richness[df_treatment_bc_totals$group_location_dpi=="wmel_body_14"]), "sd" = ds_treatment_bc_totals$sd[ds_treatment_bc_totals$group_location_dpi=="wmel_body_14"])

temp <- t.test(x           = temp_tet_legs_7$V1,
               y           = temp_tet_body_14$V1,
               mu          = 0,
               conf.level  = 0.95,
               paired      = F,
               alternative = "two.sided")


#Fig1
ggsave("Fig3_barcodes.pdf", plot,
       width = 5, height = 5, 
       units = "in", dpi = 320)
ggsave("Fig3_barcodes_dot.pdf", plot1,
       width = 7.5, height = 2.5, 
       units = "in", dpi = 320)
ggsave("Fig3_barcodes_relprop.pdf", plot_rel_prop,
       width = 10, height = 5, 
       units = "in", dpi = 320)
ggsave("Fig4_barcodes_calcs.pdf", plot_bc_stats,
       width = 5, height = 5, 
       units = "in", dpi = 320)

#Supp 1


#Supp 2