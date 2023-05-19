## Merge barcode frequency tables
rep1<-read.delim("/PATH/TO/REP1_KMERCOUNTS.tsv", header=TRUE, sep="")
rep1<-rep1[c(1,2)]
colnames(rep1)<-c("count_rep1","barcode")
totalumis_rep1<-sum(rep1$count_rep1)
rep1$freq_rep1<-(rep1$count_rep1)/totalumis_rep1



rep2<-read.delim("/PATH/TO/REP2_KMERCOUNTS.tsv", header=TRUE, sep="")
rep2<-rep2[c(1,2)]
colnames(rep2)<-c("count_rep2","barcode")
totalumis_rep2<-sum(rep2$count_rep2)
rep2$freq_rep2<-(rep2$count_rep2)/totalumis_rep2

rep1and2<-merge(rep1, rep2, by = "barcode", all.x = TRUE, all.y = TRUE)
rep1and2[is.na(rep1and2)]<-0
rep1and2$cumulativefreq<-(rep1and2$freq_rep1+rep1and2$freq_rep2)/2
rep1and2<-rep1and2[c(1,6)]

write.csv(rep1and2,"/PATH/TO/OUTPUT/COMBINED_KMER_COUNTS.tsv",row.names = FALSE)
