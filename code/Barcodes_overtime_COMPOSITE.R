require(dplyr)
require(ggplot2)

#### Mouse serial passages ####

#### Import and merge Lineage A barcode frequencies (output from Mean_bc_frequencies.R) ####

## File naming structure is as follows: "SCP" = mouse serial passage, ##
## "1_1" = passage 1_lineage 1 (lineage names were interchangeably alphabet or numeric, so 1=A, 2=B, and so on), ##
## rep12_combined = barcode frequencies were averaged across two sequencing replicates ##

## variable naming structure is similar..."p1_L1" = passage 1, lineage 1 ##

p1_L1<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP1_1_rep12_combined.csv", header=TRUE, sep=",")
p1_L1<-p1_L1[c(1,2)]
colnames(p1_L1)<-c("barcode","count_p1")

p2_L1<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP2_1_rep12_combined.csv", header=TRUE, sep=",")
p2_L1<-p2_L1[c(1,2)]
colnames(p2_L1)<-c("barcode","count_p2")

p1to2_L1<-merge(p1_L1, p2_L1, by = "barcode", all.x = TRUE, all.y = TRUE)

p3_L1<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP3_1_rep12_combined.csv", header=TRUE, sep=",")
p3_L1<-p3_L1[c(1,2)]
colnames(p3_L1)<-c("barcode","count_p3")

p1to3_L1<-merge(p1to2_L1, p3_L1, by = "barcode", all.x = TRUE, all.y = TRUE)

p4_L1<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP4_1_rep12_combined.csv", header=TRUE, sep=",")
p4_L1<-p4_L1[c(1,2)]
colnames(p4_L1)<-c("barcode","count_p4")

p1to4_L1<-merge(p1to3_L1, p4_L1, by = "barcode", all.x = TRUE, all.y = TRUE)

p5_L1<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP5_1_rep12_combined.csv", header=TRUE, sep=",")
p5_L1<-p5_L1[c(1,2)]
colnames(p5_L1)<-c("barcode","count_p5")

p1to5_L1<-merge(p1to4_L1, p5_L1, by = "barcode", all.x = TRUE, all.y = TRUE)

p6_L1<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP6_1_rep12_combined.csv", header=TRUE, sep=",")
p6_L1<-p6_L1[c(1,2)]
colnames(p6_L1)<-c("barcode","count_p6")

p1to6_L1<-merge(p1to5_L1, p6_L1, by = "barcode", all.x = TRUE, all.y = TRUE)

p7_L1<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP7_1_rep12_combined.csv", header=TRUE, sep=",")
p7_L1<-p7_L1[c(1,2)]
colnames(p7_L1)<-c("barcode","count_p7")

p1to7_L1<-merge(p1to6_L1, p7_L1, by = "barcode", all.x = TRUE, all.y = TRUE)

p8_L1<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP8_1_rep12_combined.csv", header=TRUE, sep=",")
p8_L1<-p8_L1[c(1,2)]
colnames(p8_L1)<-c("barcode","count_p8")

p1to8_L1<-merge(p1to7_L1, p8_L1, by = "barcode", all.x = TRUE, all.y = TRUE)

p9_L1<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP9_1_rep12_combined.csv", header=TRUE, sep=",")
p9_L1<-p9_L1[c(1,2)]
colnames(p9_L1)<-c("barcode","count_p9")

p1to9_L1<-merge(p1to8_L1, p9_L1, by = "barcode", all.x = TRUE, all.y = TRUE)

p10_L1<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP10_1_rep12_combined.csv", header=TRUE, sep=",")
p10_L1<-p10_L1[c(1,2)]
colnames(p10_L1)<-c("barcode","count_p10")

p1to10_L1<-merge(p1to9_L1, p10_L1, by = "barcode", all.x = TRUE, all.y = TRUE)

p1to10_L1<-subset(p1to10_L1,p1to10_L1$count_p1>0)

p1to10_L1[is.na(p1to10_L1)]<-0
p1total_L1<-sum(p1to10_L1$count_p1)
p2total_L1<-sum(p1to10_L1$count_p2)
p3total_L1<-sum(p1to10_L1$count_p3)
p4total_L1<-sum(p1to10_L1$count_p4)
p5total_L1<-sum(p1to10_L1$count_p5)
p6total_L1<-sum(p1to10_L1$count_p6)
p7total_L1<-sum(p1to10_L1$count_p7)
p8total_L1<-sum(p1to10_L1$count_p8)
p9total_L1<-sum(p1to10_L1$count_p9)
p10total_L1<-sum(p1to10_L1$count_p10)

p1to10_L1$percent_p1<-p1to10_L1$count_p1/p1total_L1*100
p1to10_L1$percent_p2<-p1to10_L1$count_p2/p2total_L1*100
p1to10_L1$percent_p3<-p1to10_L1$count_p3/p3total_L1*100
p1to10_L1$percent_p4<-p1to10_L1$count_p4/p4total_L1*100
p1to10_L1$percent_p5<-p1to10_L1$count_p5/p5total_L1*100
p1to10_L1$percent_p6<-p1to10_L1$count_p6/p6total_L1*100
p1to10_L1$percent_p7<-p1to10_L1$count_p7/p7total_L1*100
p1to10_L1$percent_p8<-p1to10_L1$count_p8/p8total_L1*100
p1to10_L1$percent_p9<-p1to10_L1$count_p9/p9total_L1*100
p1to10_L1$percent_p10<-p1to10_L1$count_p10/p10total_L1*100

#### Import and merge Lineage B barcode frequencies ####

p1_L2<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP1_2_rep12_combined.csv", header=TRUE, sep=",")
p1_L2<-p1_L2[c(1,2)]
colnames(p1_L2)<-c("barcode","count_p1")

p2_L2<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP2_2_rep12_combined.csv", header=TRUE, sep=",")
p2_L2<-p2_L2[c(1,2)]
colnames(p2_L2)<-c("barcode","count_p2")

p1to2_L2<-merge(p1_L2, p2_L2, by = "barcode", all.x = TRUE, all.y = TRUE)

p3_L2<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP3_2_rep12_combined.csv", header=TRUE, sep=",")
p3_L2<-p3_L2[c(1,2)]
colnames(p3_L2)<-c("barcode","count_p3")

p1to3_L2<-merge(p1to2_L2, p3_L2, by = "barcode", all.x = TRUE, all.y = TRUE)

p4_L2<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP4_2_rep12_combined.csv", header=TRUE, sep=",")
p4_L2<-p4_L2[c(1,2)]
colnames(p4_L2)<-c("barcode","count_p4")

p1to4_L2<-merge(p1to3_L2, p4_L2, by = "barcode", all.x = TRUE, all.y = TRUE)

p5_L2<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP5_2_rep12_combined.csv", header=TRUE, sep=",")
p5_L2<-p5_L2[c(1,2)]
colnames(p5_L2)<-c("barcode","count_p5")

p1to5_L2<-merge(p1to4_L2, p5_L2, by = "barcode", all.x = TRUE, all.y = TRUE)

p6_L2<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP6_2_rep12_combined.csv", header=TRUE, sep=",")
p6_L2<-p6_L2[c(1,2)]
colnames(p6_L2)<-c("barcode","count_p6")

p1to6_L2<-merge(p1to5_L2, p6_L2, by = "barcode", all.x = TRUE, all.y = TRUE)

p7_L2<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP7_2_rep12_combined.csv", header=TRUE, sep=",")
p7_L2<-p7_L2[c(1,2)]
colnames(p7_L2)<-c("barcode","count_p7")

p1to7_L2<-merge(p1to6_L2, p7_L2, by = "barcode", all.x = TRUE, all.y = TRUE)

p8_L2<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP8_2_rep12_combined.csv", header=TRUE, sep=",")
p8_L2<-p8_L2[c(1,2)]
colnames(p8_L2)<-c("barcode","count_p8")

p1to8_L2<-merge(p1to7_L2, p8_L2, by = "barcode", all.x = TRUE, all.y = TRUE)

p9_L2<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP9_2_rep12_combined.csv", header=TRUE, sep=",")
p9_L2<-p9_L2[c(1,2)]
colnames(p9_L2)<-c("barcode","count_p9")

p1to9_L2<-merge(p1to8_L2, p9_L2, by = "barcode", all.x = TRUE, all.y = TRUE)

p10_L2<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP10_2_rep12_combined.csv", header=TRUE, sep=",")
p10_L2<-p10_L2[c(1,2)]
colnames(p10_L2)<-c("barcode","count_p10")

p1to10_L2<-merge(p1to9_L2, p10_L2, by = "barcode", all.x = TRUE, all.y = TRUE)

p1to10_L2<-subset(p1to10_L2,p1to10_L2$count_p1>0)

p1to10_L2[is.na(p1to10_L2)]<-0
p1total_L2<-sum(p1to10_L2$count_p1)
p2total_L2<-sum(p1to10_L2$count_p2)
p3total_L2<-sum(p1to10_L2$count_p3)
p4total_L2<-sum(p1to10_L2$count_p4)
p5total_L2<-sum(p1to10_L2$count_p5)
p6total_L2<-sum(p1to10_L2$count_p6)
p7total_L2<-sum(p1to10_L2$count_p7)
p8total_L2<-sum(p1to10_L2$count_p8)
p9total_L2<-sum(p1to10_L2$count_p9)
p10total_L2<-sum(p1to10_L2$count_p10)

p1to10_L2$percent_p1<-p1to10_L2$count_p1/p1total_L2*100
p1to10_L2$percent_p2<-p1to10_L2$count_p2/p2total_L2*100
p1to10_L2$percent_p3<-p1to10_L2$count_p3/p3total_L2*100
p1to10_L2$percent_p4<-p1to10_L2$count_p4/p4total_L2*100
p1to10_L2$percent_p5<-p1to10_L2$count_p5/p5total_L2*100
p1to10_L2$percent_p6<-p1to10_L2$count_p6/p6total_L2*100
p1to10_L2$percent_p7<-p1to10_L2$count_p7/p7total_L2*100
p1to10_L2$percent_p8<-p1to10_L2$count_p8/p8total_L2*100
p1to10_L2$percent_p9<-p1to10_L2$count_p9/p9total_L2*100
p1to10_L2$percent_p10<-p1to10_L2$count_p10/p10total_L2*100

#### Import and merge Lineage C barcode frequencies ####

p1_L3<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP1_3_rep12_combined.csv", header=TRUE, sep=",")
p1_L3<-p1_L3[c(1,2)]
colnames(p1_L3)<-c("barcode","count_p1")

p2_L3<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP2_3_rep12_combined.csv", header=TRUE, sep=",")
p2_L3<-p2_L3[c(1,2)]
colnames(p2_L3)<-c("barcode","count_p2")

p1to2_L3<-merge(p1_L3, p2_L3, by = "barcode", all.x = TRUE, all.y = TRUE)

p3_L3<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP3_3_rep12_combined.csv", header=TRUE, sep=",")
p3_L3<-p3_L3[c(1,2)]
colnames(p3_L3)<-c("barcode","count_p3")

p1to3_L3<-merge(p1to2_L3, p3_L3, by = "barcode", all.x = TRUE, all.y = TRUE)

p4_L3<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP4_3_rep12_combined.csv", header=TRUE, sep=",")
p4_L3<-p4_L3[c(1,2)]
colnames(p4_L3)<-c("barcode","count_p4")

p1to4_L3<-merge(p1to3_L3, p4_L3, by = "barcode", all.x = TRUE, all.y = TRUE)

p5_L3<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP5_3_rep12_combined.csv", header=TRUE, sep=",")
p5_L3<-p5_L3[c(1,2)]
colnames(p5_L3)<-c("barcode","count_p5")

p1to5_L3<-merge(p1to4_L3, p5_L3, by = "barcode", all.x = TRUE, all.y = TRUE)

p6_L3<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP6_3_rep12_combined.csv", header=TRUE, sep=",")
p6_L3<-p6_L3[c(1,2)]
colnames(p6_L3)<-c("barcode","count_p6")

p1to6_L3<-merge(p1to5_L3, p6_L3, by = "barcode", all.x = TRUE, all.y = TRUE)

p7_L3<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP7_3_rep12_combined.csv", header=TRUE, sep=",")
p7_L3<-p7_L3[c(1,2)]
colnames(p7_L3)<-c("barcode","count_p7")

p1to7_L3<-merge(p1to6_L3, p7_L3, by = "barcode", all.x = TRUE, all.y = TRUE)

p8_L3<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP8_3_rep12_combined.csv", header=TRUE, sep=",")
p8_L3<-p8_L3[c(1,2)]
colnames(p8_L3)<-c("barcode","count_p8")

p1to8_L3<-merge(p1to7_L3, p8_L3, by = "barcode", all.x = TRUE, all.y = TRUE)

p9_L3<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP9_3_rep12_combined.csv", header=TRUE, sep=",")
p9_L3<-p9_L3[c(1,2)]
colnames(p9_L3)<-c("barcode","count_p9")

p1to9_L3<-merge(p1to8_L3, p9_L3, by = "barcode", all.x = TRUE, all.y = TRUE)

p10_L3<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP10_3_rep12_combined.csv", header=TRUE, sep=",")
p10_L3<-p10_L3[c(1,2)]
colnames(p10_L3)<-c("barcode","count_p10")

p1to10_L3<-merge(p1to9_L3, p10_L3, by = "barcode", all.x = TRUE, all.y = TRUE)

p1to10_L3<-subset(p1to10_L3,p1to10_L3$count_p1>0)

p1to10_L3[is.na(p1to10_L3)]<-0
p1total_L3<-sum(p1to10_L3$count_p1)
p2total_L3<-sum(p1to10_L3$count_p2)
p3total_L3<-sum(p1to10_L3$count_p3)
p4total_L3<-sum(p1to10_L3$count_p4)
p5total_L3<-sum(p1to10_L3$count_p5)
p6total_L3<-sum(p1to10_L3$count_p6)
p7total_L3<-sum(p1to10_L3$count_p7)
p8total_L3<-sum(p1to10_L3$count_p8)
p9total_L3<-sum(p1to10_L3$count_p9)
p10total_L3<-sum(p1to10_L3$count_p10)

p1to10_L3$percent_p1<-p1to10_L3$count_p1/p1total_L3*100
p1to10_L3$percent_p2<-p1to10_L3$count_p2/p2total_L3*100
p1to10_L3$percent_p3<-p1to10_L3$count_p3/p3total_L3*100
p1to10_L3$percent_p4<-p1to10_L3$count_p4/p4total_L3*100
p1to10_L3$percent_p5<-p1to10_L3$count_p5/p5total_L3*100
p1to10_L3$percent_p6<-p1to10_L3$count_p6/p6total_L3*100
p1to10_L3$percent_p7<-p1to10_L3$count_p7/p7total_L3*100
p1to10_L3$percent_p8<-p1to10_L3$count_p8/p8total_L3*100
p1to10_L3$percent_p9<-p1to10_L3$count_p9/p9total_L3*100
p1to10_L3$percent_p10<-p1to10_L3$count_p10/p10total_L3*100

#### Import and merge Lineage D barcode frequencies ####

p1_L4<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP1_4_rep12_combined.csv", header=TRUE, sep=",")
p1_L4<-p1_L4[c(1,2)]
colnames(p1_L4)<-c("barcode","count_p1")

p2_L4<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP2_4_rep12_combined.csv", header=TRUE, sep=",")
p2_L4<-p2_L4[c(1,2)]
colnames(p2_L4)<-c("barcode","count_p2")

p1to2_L4<-merge(p1_L4, p2_L4, by = "barcode", all.x = TRUE, all.y = TRUE)

p3_L4<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP3_4_rep12_combined.csv", header=TRUE, sep=",")
p3_L4<-p3_L4[c(1,2)]
colnames(p3_L4)<-c("barcode","count_p3")

p1to3_L4<-merge(p1to2_L4, p3_L4, by = "barcode", all.x = TRUE, all.y = TRUE)

p4_L4<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP4_4_rep12_combined.csv", header=TRUE, sep=",")
p4_L4<-p4_L4[c(1,2)]
colnames(p4_L4)<-c("barcode","count_p4")

p1to4_L4<-merge(p1to3_L4, p4_L4, by = "barcode", all.x = TRUE, all.y = TRUE)

p5_L4<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP5_4_rep12_combined.csv", header=TRUE, sep=",")
p5_L4<-p5_L4[c(1,2)]
colnames(p5_L4)<-c("barcode","count_p5")

p1to5_L4<-merge(p1to4_L4, p5_L4, by = "barcode", all.x = TRUE, all.y = TRUE)

p6_L4<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP6_4_rep12_combined.csv", header=TRUE, sep=",")
p6_L4<-p6_L4[c(1,2)]
colnames(p6_L4)<-c("barcode","count_p6")

p1to6_L4<-merge(p1to5_L4, p6_L4, by = "barcode", all.x = TRUE, all.y = TRUE)

p7_L4<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP7_4_rep12_combined.csv", header=TRUE, sep=",")
p7_L4<-p7_L4[c(1,2)]
colnames(p7_L4)<-c("barcode","count_p7")

p1to7_L4<-merge(p1to6_L4, p7_L4, by = "barcode", all.x = TRUE, all.y = TRUE)

p8_L4<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP8_4_rep12_combined.csv", header=TRUE, sep=",")
p8_L4<-p8_L4[c(1,2)]
colnames(p8_L4)<-c("barcode","count_p8")

p1to8_L4<-merge(p1to7_L4, p8_L4, by = "barcode", all.x = TRUE, all.y = TRUE)

p9_L4<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP9_4_rep12_combined.csv", header=TRUE, sep=",")
p9_L4<-p9_L4[c(1,2)]
colnames(p9_L4)<-c("barcode","count_p9")

p1to9_L4<-merge(p1to8_L4, p9_L4, by = "barcode", all.x = TRUE, all.y = TRUE)

p10_L4<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP10_4_rep12_combined.csv", header=TRUE, sep=",")
p10_L4<-p10_L4[c(1,2)]
colnames(p10_L4)<-c("barcode","count_p10")

p1to10_L4<-merge(p1to9_L4, p10_L4, by = "barcode", all.x = TRUE, all.y = TRUE)

p1to10_L4<-subset(p1to10_L4,p1to10_L4$count_p1>0)

p1to10_L4[is.na(p1to10_L4)]<-0
p1total_L4<-sum(p1to10_L4$count_p1)
p2total_L4<-sum(p1to10_L4$count_p2)
p3total_L4<-sum(p1to10_L4$count_p3)
p4total_L4<-sum(p1to10_L4$count_p4)
p5total_L4<-sum(p1to10_L4$count_p5)
p6total_L4<-sum(p1to10_L4$count_p6)
p7total_L4<-sum(p1to10_L4$count_p7)
p8total_L4<-sum(p1to10_L4$count_p8)
p9total_L4<-sum(p1to10_L4$count_p9)
p10total_L4<-sum(p1to10_L4$count_p10)

p1to10_L4$percent_p1<-p1to10_L4$count_p1/p1total_L4*100
p1to10_L4$percent_p2<-p1to10_L4$count_p2/p2total_L4*100
p1to10_L4$percent_p3<-p1to10_L4$count_p3/p3total_L4*100
p1to10_L4$percent_p4<-p1to10_L4$count_p4/p4total_L4*100
p1to10_L4$percent_p5<-p1to10_L4$count_p5/p5total_L4*100
p1to10_L4$percent_p6<-p1to10_L4$count_p6/p6total_L4*100
p1to10_L4$percent_p7<-p1to10_L4$count_p7/p7total_L4*100
p1to10_L4$percent_p8<-p1to10_L4$count_p8/p8total_L4*100
p1to10_L4$percent_p9<-p1to10_L4$count_p9/p9total_L4*100
p1to10_L4$percent_p10<-p1to10_L4$count_p10/p10total_L4*100

#### Import and merge Lineage E barcode frequencies ####

p1_L5<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP1_5_rep12_combined.csv", header=TRUE, sep=",")
p1_L5<-p1_L5[c(1,2)]
colnames(p1_L5)<-c("barcode","count_p1")

p2_L5<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP2_5_rep12_combined.csv", header=TRUE, sep=",")
p2_L5<-p2_L5[c(1,2)]
colnames(p2_L5)<-c("barcode","count_p2")

p1to2_L5<-merge(p1_L5, p2_L5, by = "barcode", all.x = TRUE, all.y = TRUE)

p3_L5<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP3_5_rep12_combined.csv", header=TRUE, sep=",")
p3_L5<-p3_L5[c(1,2)]
colnames(p3_L5)<-c("barcode","count_p3")

p1to3_L5<-merge(p1to2_L5, p3_L5, by = "barcode", all.x = TRUE, all.y = TRUE)

p4_L5<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP4_5_rep12_combined.csv", header=TRUE, sep=",")
p4_L5<-p4_L5[c(1,2)]
colnames(p4_L5)<-c("barcode","count_p4")

p1to4_L5<-merge(p1to3_L5, p4_L5, by = "barcode", all.x = TRUE, all.y = TRUE)

p5_L5<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP5_5_rep12_combined.csv", header=TRUE, sep=",")
p5_L5<-p5_L5[c(1,2)]
colnames(p5_L5)<-c("barcode","count_p5")

p1to5_L5<-merge(p1to4_L5, p5_L5, by = "barcode", all.x = TRUE, all.y = TRUE)

p6_L5<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP6_5_rep12_combined.csv", header=TRUE, sep=",")
p6_L5<-p6_L5[c(1,2)]
colnames(p6_L5)<-c("barcode","count_p6")

p1to6_L5<-merge(p1to5_L5, p6_L5, by = "barcode", all.x = TRUE, all.y = TRUE)

p7_L5<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP7_5_rep12_combined.csv", header=TRUE, sep=",")
p7_L5<-p7_L5[c(1,2)]
colnames(p7_L5)<-c("barcode","count_p7")

p1to7_L5<-merge(p1to6_L5, p7_L5, by = "barcode", all.x = TRUE, all.y = TRUE)

p8_L5<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP8_5_rep12_combined.csv", header=TRUE, sep=",")
p8_L5<-p8_L5[c(1,2)]
colnames(p8_L5)<-c("barcode","count_p8")

p1to8_L5<-merge(p1to7_L5, p8_L5, by = "barcode", all.x = TRUE, all.y = TRUE)

p9_L5<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP9_5_rep12_combined.csv", header=TRUE, sep=",")
p9_L5<-p9_L5[c(1,2)]
colnames(p9_L5)<-c("barcode","count_p9")

p1to9_L5<-merge(p1to8_L5, p9_L5, by = "barcode", all.x = TRUE, all.y = TRUE)

p10_L5<-read.delim("/PATH/TO/Mouse_Passage/Barcodes/SCP10_5_rep12_combined.csv", header=TRUE, sep=",")
p10_L5<-p10_L5[c(1,2)]
colnames(p10_L5)<-c("barcode","count_p10")

p1to10_L5<-merge(p1to9_L5, p10_L5, by = "barcode", all.x = TRUE, all.y = TRUE)

p1to10_L5<-subset(p1to10_L5,p1to10_L5$count_p1>0)

p1to10_L5[is.na(p1to10_L5)]<-0
p1total_L5<-sum(p1to10_L5$count_p1)
p2total_L5<-sum(p1to10_L5$count_p2)
p3total_L5<-sum(p1to10_L5$count_p3)
p4total_L5<-sum(p1to10_L5$count_p4)
p5total_L5<-sum(p1to10_L5$count_p5)
p6total_L5<-sum(p1to10_L5$count_p6)
p7total_L5<-sum(p1to10_L5$count_p7)
p8total_L5<-sum(p1to10_L5$count_p8)
p9total_L5<-sum(p1to10_L5$count_p9)
p10total_L5<-sum(p1to10_L5$count_p10)

p1to10_L5$percent_p1<-p1to10_L5$count_p1/p1total_L5*100
p1to10_L5$percent_p2<-p1to10_L5$count_p2/p2total_L5*100
p1to10_L5$percent_p3<-p1to10_L5$count_p3/p3total_L5*100
p1to10_L5$percent_p4<-p1to10_L5$count_p4/p4total_L5*100
p1to10_L5$percent_p5<-p1to10_L5$count_p5/p5total_L5*100
p1to10_L5$percent_p6<-p1to10_L5$count_p6/p6total_L5*100
p1to10_L5$percent_p7<-p1to10_L5$count_p7/p7total_L5*100
p1to10_L5$percent_p8<-p1to10_L5$count_p8/p8total_L5*100
p1to10_L5$percent_p9<-p1to10_L5$count_p9/p9total_L5*100
p1to10_L5$percent_p10<-p1to10_L5$count_p10/p10total_L5*100


#### Merge all lineages by passage ####

p1to10_L1<-arrange(p1to10_L1, desc(percent_p10))
p1to10_L1$p10_rank<-1:nrow(p1to10_L1)
p1to10_L2<-arrange(p1to10_L2, desc(percent_p10))
p1to10_L2$p10_rank<-1:nrow(p1to10_L2)
p1to10_L3<-arrange(p1to10_L3, desc(percent_p10))
p1to10_L3$p10_rank<-1:nrow(p1to10_L3)
p1to10_L4<-arrange(p1to10_L4, desc(percent_p10))
p1to10_L4$p10_rank<-1:nrow(p1to10_L4)
p1to10_L5<-arrange(p1to10_L5, desc(percent_p10))
p1to10_L5$p10_rank<-1:nrow(p1to10_L5)

p1_ALL<-data.frame("p10_rank" = 1:2321)
p1_ALL<-merge(p1_ALL, p1to10_L1, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p1_ALL<-p1_ALL[,c(1,3)]
colnames(p1_ALL)<-c("p10_rank","p1_L1")
p1_ALL<-merge(p1_ALL, p1to10_L2, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p1_ALL<-p1_ALL[,c(1,2,4)]
colnames(p1_ALL)<-c("p10_rank","p1_L1","p1_L2")
p1_ALL<-merge(p1_ALL, p1to10_L3, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p1_ALL<-p1_ALL[,c(1,2,3,5)]
colnames(p1_ALL)<-c("p10_rank","p1_L1","p1_L2","p1_L3")
p1_ALL<-merge(p1_ALL, p1to10_L4, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p1_ALL<-p1_ALL[,c(1,2,3,4,6)]
colnames(p1_ALL)<-c("p10_rank","p1_L1","p1_L2","p1_L3","p1_L4")
p1_ALL<-merge(p1_ALL, p1to10_L5, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p1_ALL<-p1_ALL[,c(1,2,3,4,5,7)]
colnames(p1_ALL)<-c("p10_rank","p1_L1","p1_L2","p1_L3","p1_L4","p1_L5")
p1_ALL[is.na(p1_ALL)]<-0

p2_ALL<-data.frame("p10_rank" = 1:2321)
p2_ALL<-merge(p2_ALL, p1to10_L1, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p2_ALL<-p2_ALL[,c(1,4)]
colnames(p2_ALL)<-c("p10_rank","p2_L1")
p2_ALL<-merge(p2_ALL, p1to10_L2, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p2_ALL<-p2_ALL[,c(1,2,5)]
colnames(p2_ALL)<-c("p10_rank","p2_L1","p2_L2")
p2_ALL<-merge(p2_ALL, p1to10_L3, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p2_ALL<-p2_ALL[,c(1,2,3,6)]
colnames(p2_ALL)<-c("p10_rank","p2_L1","p2_L2","p2_L3")
p2_ALL<-merge(p2_ALL, p1to10_L4, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p2_ALL<-p2_ALL[,c(1,2,3,4,7)]
colnames(p2_ALL)<-c("p10_rank","p2_L1","p2_L2","p2_L3","p2_L4")
p2_ALL<-merge(p2_ALL, p1to10_L5, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p2_ALL<-p2_ALL[,c(1,2,3,4,5,8)]
colnames(p2_ALL)<-c("p10_rank","p2_L1","p2_L2","p2_L3","p2_L4","p2_L5")
p2_ALL[is.na(p2_ALL)]<-0

p3_ALL<-data.frame("p10_rank" = 1:2321)
p3_ALL<-merge(p3_ALL, p1to10_L1, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p3_ALL<-p3_ALL[,c(1,5)]
colnames(p3_ALL)<-c("p10_rank","p3_L1")
p3_ALL<-merge(p3_ALL, p1to10_L2, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p3_ALL<-p3_ALL[,c(1,2,6)]
colnames(p3_ALL)<-c("p10_rank","p3_L1","p3_L2")
p3_ALL<-merge(p3_ALL, p1to10_L3, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p3_ALL<-p3_ALL[,c(1,2,3,7)]
colnames(p3_ALL)<-c("p10_rank","p3_L1","p3_L2","p3_L3")
p3_ALL<-merge(p3_ALL, p1to10_L4, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p3_ALL<-p3_ALL[,c(1,2,3,4,8)]
colnames(p3_ALL)<-c("p10_rank","p3_L1","p3_L2","p3_L3","p3_L4")
p3_ALL<-merge(p3_ALL, p1to10_L5, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p3_ALL<-p3_ALL[,c(1,2,3,4,5,9)]
colnames(p3_ALL)<-c("p10_rank","p3_L1","p3_L2","p3_L3","p3_L4","p3_L5")
p3_ALL[is.na(p3_ALL)]<-0

p4_ALL<-data.frame("p10_rank" = 1:2321)
p4_ALL<-merge(p4_ALL, p1to10_L1, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p4_ALL<-p4_ALL[,c(1,6)]
colnames(p4_ALL)<-c("p10_rank","p4_L1")
p4_ALL<-merge(p4_ALL, p1to10_L2, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p4_ALL<-p4_ALL[,c(1,2,7)]
colnames(p4_ALL)<-c("p10_rank","p4_L1","p4_L2")
p4_ALL<-merge(p4_ALL, p1to10_L3, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p4_ALL<-p4_ALL[,c(1,2,3,8)]
colnames(p4_ALL)<-c("p10_rank","p4_L1","p4_L2","p4_L3")
p4_ALL<-merge(p4_ALL, p1to10_L4, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p4_ALL<-p4_ALL[,c(1,2,3,4,9)]
colnames(p4_ALL)<-c("p10_rank","p4_L1","p4_L2","p4_L3","p4_L4")
p4_ALL<-merge(p4_ALL, p1to10_L5, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p4_ALL<-p4_ALL[,c(1,2,3,4,5,10)]
colnames(p4_ALL)<-c("p10_rank","p4_L1","p4_L2","p4_L3","p4_L4","p4_L5")
p4_ALL[is.na(p4_ALL)]<-0

p5_ALL<-data.frame("p10_rank" = 1:2321)
p5_ALL<-merge(p5_ALL, p1to10_L1, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p5_ALL<-p5_ALL[,c(1,7)]
colnames(p5_ALL)<-c("p10_rank","p5_L1")
p5_ALL<-merge(p5_ALL, p1to10_L2, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p5_ALL<-p5_ALL[,c(1,2,8)]
colnames(p5_ALL)<-c("p10_rank","p5_L1","p5_L2")
p5_ALL<-merge(p5_ALL, p1to10_L3, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p5_ALL<-p5_ALL[,c(1,2,3,9)]
colnames(p5_ALL)<-c("p10_rank","p5_L1","p5_L2","p5_L3")
p5_ALL<-merge(p5_ALL, p1to10_L4, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p5_ALL<-p5_ALL[,c(1,2,3,4,10)]
colnames(p5_ALL)<-c("p10_rank","p5_L1","p5_L2","p5_L3","p5_L4")
p5_ALL<-merge(p5_ALL, p1to10_L5, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p5_ALL<-p5_ALL[,c(1,2,3,4,5,11)]
colnames(p5_ALL)<-c("p10_rank","p5_L1","p5_L2","p5_L3","p5_L4","p5_L5")
p5_ALL[is.na(p5_ALL)]<-0

p6_ALL<-data.frame("p10_rank" = 1:2321)
p6_ALL<-merge(p6_ALL, p1to10_L1, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p6_ALL<-p6_ALL[,c(1,8)]
colnames(p6_ALL)<-c("p10_rank","p6_L1")
p6_ALL<-merge(p6_ALL, p1to10_L2, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p6_ALL<-p6_ALL[,c(1,2,9)]
colnames(p6_ALL)<-c("p10_rank","p6_L1","p6_L2")
p6_ALL<-merge(p6_ALL, p1to10_L3, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p6_ALL<-p6_ALL[,c(1,2,3,10)]
colnames(p6_ALL)<-c("p10_rank","p6_L1","p6_L2","p6_L3")
p6_ALL<-merge(p6_ALL, p1to10_L4, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p6_ALL<-p6_ALL[,c(1,2,3,4,11)]
colnames(p6_ALL)<-c("p10_rank","p6_L1","p6_L2","p6_L3","p6_L4")
p6_ALL<-merge(p6_ALL, p1to10_L5, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p6_ALL<-p6_ALL[,c(1,2,3,4,5,12)]
colnames(p6_ALL)<-c("p10_rank","p6_L1","p6_L2","p6_L3","p6_L4","p6_L5")
p6_ALL[is.na(p6_ALL)]<-0

p7_ALL<-data.frame("p10_rank" = 1:2321)
p7_ALL<-merge(p7_ALL, p1to10_L1, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p7_ALL<-p7_ALL[,c(1,9)]
colnames(p7_ALL)<-c("p10_rank","p7_L1")
p7_ALL<-merge(p7_ALL, p1to10_L2, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p7_ALL<-p7_ALL[,c(1,2,10)]
colnames(p7_ALL)<-c("p10_rank","p7_L1","p7_L2")
p7_ALL<-merge(p7_ALL, p1to10_L3, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p7_ALL<-p7_ALL[,c(1,2,3,11)]
colnames(p7_ALL)<-c("p10_rank","p7_L1","p7_L2","p7_L3")
p7_ALL<-merge(p7_ALL, p1to10_L4, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p7_ALL<-p7_ALL[,c(1,2,3,4,12)]
colnames(p7_ALL)<-c("p10_rank","p7_L1","p7_L2","p7_L3","p7_L4")
p7_ALL<-merge(p7_ALL, p1to10_L5, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p7_ALL<-p7_ALL[,c(1,2,3,4,5,13)]
colnames(p7_ALL)<-c("p10_rank","p7_L1","p7_L2","p7_L3","p7_L4","p7_L5")
p7_ALL[is.na(p7_ALL)]<-0

p8_ALL<-data.frame("p10_rank" = 1:2321)
p8_ALL<-merge(p8_ALL, p1to10_L1, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p8_ALL<-p8_ALL[,c(1,10)]
colnames(p8_ALL)<-c("p10_rank","p8_L1")
p8_ALL<-merge(p8_ALL, p1to10_L2, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p8_ALL<-p8_ALL[,c(1,2,11)]
colnames(p8_ALL)<-c("p10_rank","p8_L1","p8_L2")
p8_ALL<-merge(p8_ALL, p1to10_L3, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p8_ALL<-p8_ALL[,c(1,2,3,12)]
colnames(p8_ALL)<-c("p10_rank","p8_L1","p8_L2","p8_L3")
p8_ALL<-merge(p8_ALL, p1to10_L4, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p8_ALL<-p8_ALL[,c(1,2,3,4,13)]
colnames(p8_ALL)<-c("p10_rank","p8_L1","p8_L2","p8_L3","p8_L4")
p8_ALL<-merge(p8_ALL, p1to10_L5, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p8_ALL<-p8_ALL[,c(1,2,3,4,5,14)]
colnames(p8_ALL)<-c("p10_rank","p8_L1","p8_L2","p8_L3","p8_L4","p8_L5")
p8_ALL[is.na(p8_ALL)]<-0

p9_ALL<-data.frame("p10_rank" = 1:2321)
p9_ALL<-merge(p9_ALL, p1to10_L1, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p9_ALL<-p9_ALL[,c(1,11)]
colnames(p9_ALL)<-c("p10_rank","p9_L1")
p9_ALL<-merge(p9_ALL, p1to10_L2, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p9_ALL<-p9_ALL[,c(1,2,12)]
colnames(p9_ALL)<-c("p10_rank","p9_L1","p9_L2")
p9_ALL<-merge(p9_ALL, p1to10_L3, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p9_ALL<-p9_ALL[,c(1,2,3,13)]
colnames(p9_ALL)<-c("p10_rank","p9_L1","p9_L2","p9_L3")
p9_ALL<-merge(p9_ALL, p1to10_L4, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p9_ALL<-p9_ALL[,c(1,2,3,4,14)]
colnames(p9_ALL)<-c("p10_rank","p9_L1","p9_L2","p9_L3","p9_L4")
p9_ALL<-merge(p9_ALL, p1to10_L5, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p9_ALL<-p9_ALL[,c(1,2,3,4,5,15)]
colnames(p9_ALL)<-c("p10_rank","p9_L1","p9_L2","p9_L3","p9_L4","p9_L5")
p9_ALL[is.na(p9_ALL)]<-0

p10_ALL<-data.frame("p10_rank" = 1:2321)
p10_ALL<-merge(p10_ALL, p1to10_L1, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p10_ALL<-p10_ALL[,c(1,12)]
colnames(p10_ALL)<-c("p10_rank","p10_L1")
p10_ALL<-merge(p10_ALL, p1to10_L2, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p10_ALL<-p10_ALL[,c(1,2,13)]
colnames(p10_ALL)<-c("p10_rank","p10_L1","p10_L2")
p10_ALL<-merge(p10_ALL, p1to10_L3, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p10_ALL<-p10_ALL[,c(1,2,3,14)]
colnames(p10_ALL)<-c("p10_rank","p10_L1","p10_L2","p10_L3")
p10_ALL<-merge(p10_ALL, p1to10_L4, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p10_ALL<-p10_ALL[,c(1,2,3,4,15)]
colnames(p10_ALL)<-c("p10_rank","p10_L1","p10_L2","p10_L3","p10_L4")
p10_ALL<-merge(p10_ALL, p1to10_L5, by = "p10_rank", all.x = TRUE, all.y = TRUE)
p10_ALL<-p10_ALL[,c(1,2,3,4,5,16)]
colnames(p10_ALL)<-c("p10_rank","p10_L1","p10_L2","p10_L3","p10_L4","p10_L5")
p10_ALL[is.na(p10_ALL)]<-0

p1_ALL$mean<-(p1_ALL$p1_L1+p1_ALL$p1_L2+p1_ALL$p1_L3+p1_ALL$p1_L4+p1_ALL$p1_L5)/5
p2_ALL$mean<-(p2_ALL$p2_L1+p2_ALL$p2_L2+p2_ALL$p2_L3+p2_ALL$p2_L4+p2_ALL$p2_L5)/5
p3_ALL$mean<-(p3_ALL$p3_L1+p3_ALL$p3_L2+p3_ALL$p3_L3+p3_ALL$p3_L4+p3_ALL$p3_L5)/5
p4_ALL$mean<-(p4_ALL$p4_L1+p4_ALL$p4_L2+p4_ALL$p4_L3+p4_ALL$p4_L4+p4_ALL$p4_L5)/5
p5_ALL$mean<-(p5_ALL$p5_L1+p5_ALL$p5_L2+p5_ALL$p5_L3+p5_ALL$p5_L4+p5_ALL$p5_L5)/5
p6_ALL$mean<-(p6_ALL$p6_L1+p6_ALL$p6_L2+p6_ALL$p6_L3+p6_ALL$p6_L4+p6_ALL$p6_L5)/5
p7_ALL$mean<-(p7_ALL$p7_L1+p7_ALL$p7_L2+p7_ALL$p7_L3+p7_ALL$p7_L4+p7_ALL$p7_L5)/5
p8_ALL$mean<-(p8_ALL$p8_L1+p8_ALL$p8_L2+p8_ALL$p8_L3+p8_ALL$p8_L4+p8_ALL$p8_L5)/5
p9_ALL$mean<-(p9_ALL$p9_L1+p9_ALL$p9_L2+p9_ALL$p9_L3+p9_ALL$p9_L4+p9_ALL$p9_L5)/5
p10_ALL$mean<-(p10_ALL$p10_L1+p10_ALL$p10_L2+p10_ALL$p10_L3+p10_ALL$p10_L4+p10_ALL$p10_L5)/5
p1to10_means<-data.frame(p1_ALL$p10_rank,p1_ALL$mean,p2_ALL$mean,p3_ALL$mean,p4_ALL$mean,p5_ALL$mean,p6_ALL$mean,p7_ALL$mean,p8_ALL$mean,p9_ALL$mean,p10_ALL$mean)
colnames(p1to10_means)<-c("p10_rank","p1_mean","p2_mean","p3_mean","p4_mean","p5_mean","p6_mean","p7_mean","p8_mean","p9_mean","p10_mean")

p1to10_means$percent_p1<-(p1to10_means$p1_mean/sum(p1to10_means$p1_mean))*100
p1to10_means$percent_p2<-(p1to10_means$p2_mean/sum(p1to10_means$p2_mean))*100
p1to10_means$percent_p3<-(p1to10_means$p3_mean/sum(p1to10_means$p3_mean))*100
p1to10_means$percent_p4<-(p1to10_means$p4_mean/sum(p1to10_means$p4_mean))*100
p1to10_means$percent_p5<-(p1to10_means$p5_mean/sum(p1to10_means$p5_mean))*100
p1to10_means$percent_p6<-(p1to10_means$p6_mean/sum(p1to10_means$p6_mean))*100
p1to10_means$percent_p7<-(p1to10_means$p7_mean/sum(p1to10_means$p7_mean))*100
p1to10_means$percent_p8<-(p1to10_means$p8_mean/sum(p1to10_means$p8_mean))*100
p1to10_means$percent_p9<-(p1to10_means$p9_mean/sum(p1to10_means$p9_mean))*100
p1to10_means$percent_p10<-(p1to10_means$p10_mean/sum(p1to10_means$p10_mean))*100


##############   BARCODE TRAJECTORY FOR ALL BARCODES PRESENT AT P1   ##############

p1to10_means<-arrange(p1to10_means, desc(percent_p1))

one <- p1to10_means[c("percent_p1")]
one$passage <- 1
rows<-nrow(one)
one$barcode.id <- as.factor(1:rows)
colnames(one)<-c("percent","passage","barcode_id")
two <- p1to10_means[c("percent_p2")]
two$passage <- 2
two$barcode_id <- as.factor(1:rows)
colnames(two)<-c("percent","passage","barcode_id")
three <- p1to10_means[c("percent_p3")]
three$passage <- 3
three$barcode_id <- as.factor(1:rows)
colnames(three)<-c("percent","passage","barcode_id")
four <- p1to10_means[c("percent_p4")]
four$passage <- 4
four$barcode_id <- as.factor(1:rows)
colnames(four)<-c("percent","passage","barcode_id")
five <- p1to10_means[c("percent_p5")]
five$passage <- 5
five$barcode_id <- as.factor(1:rows)
colnames(five)<-c("percent","passage","barcode_id")
six <- p1to10_means[c("percent_p6")]
six$passage <- 6
six$barcode_id <- as.factor(1:rows)
colnames(six)<-c("percent","passage","barcode_id")
seven <- p1to10_means[c("percent_p7")]
seven$passage <- 7
seven$barcode_id <- as.factor(1:rows)
colnames(seven)<-c("percent","passage","barcode_id")
eight <- p1to10_means[c("percent_p8")]
eight$passage <- 8
eight$barcode_id <- as.factor(1:rows)
colnames(eight)<-c("percent","passage","barcode_id")
nine <- p1to10_means[c("percent_p9")]
nine$passage <- 9
nine$barcode_id <- as.factor(1:rows)
colnames(nine)<-c("percent","passage","barcode_id")
ten <- p1to10_means[c("percent_p10")]
ten$passage <- 10
ten$barcode_id <- as.factor(1:rows)
colnames(ten)<-c("percent","passage","barcode_id")

percent<-rbind(one, two, three, four, five, six, seven, eight, nine, ten)

percent<-arrange(percent, desc(barcode_id))

x<-(rep(c("red","green","blue","yellow","purple","grey","orange","brown","turquoise","pink","cyan","darkgrey"),1000))

ggplot(percent, aes(x=passage, y=percent, fill=barcode_id)) + 
  geom_area(alpha=0.6, show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_line(colour="light grey"), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15)) +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  ggtitle("ZIKV Mouse Composite") +
  labs(x="Passage", y="Percent of Population") +
  scale_x_continuous(breaks=seq(1,10,by=1), expand = c(0, 0)) +
  scale_y_continuous(breaks=seq(0,100,by=10),expand = c(0, 0)) +
  scale_fill_manual(values=x)

setwd("PATH/TO/WD/")
ggsave("FILENAME.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)
