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