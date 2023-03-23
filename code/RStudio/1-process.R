## load required packages
list.of.packages <- c("diverse","plyr","readr","dplyr","splitstackshape","RcppRoll")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## upload package libraries
library("diverse")
library("plyr")
library("readr")
library("dplyr")
library("RcppRoll")
library("splitstackshape")

#homedir<-paste("/Users/rieshunter/Documents/bioinformatics/Wolbachia/data/reads/data/test")
#setwd(homedir)
#name<-paste("B1-26-dup-7-2-tet-body_S39")
#### Enter the sample name and output directory for the sample
args<-commandArgs()
homedir<-args[6]
name<-args[7]

sampledir<-paste(name)
sampledir2<-paste(sampledir)
outputdir<-paste(homedir,sampledir2,sep="/")
samplename<-paste(name)
outputsample<-paste(outputdir,samplename,sep="/")

######   DIVERSITY METRICS     #######
## import Perbase text file and lofreq file
table<-paste("13-ntcounts_consensus",samplename,sep="_")
table<-paste(table,".csv",sep="")
table<-paste(homedir,table,sep="/")
ntcounts_c<-read.csv(table, header=TRUE)

vcf<-paste("09-consensus_", samplename, sep="")
vcf<-paste(vcf, ".vcf", sep="")
vcf<-paste(homedir, vcf, sep="/")
var<-read.table(vcf, sep="\t", header=FALSE)

var<-cSplit(var,"V8",sep=";", type.convert=FALSE)
var<-cSplit(var, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
var<-var[,c(2,4,5,11,13)]
colnames(var)<-c("position","ref","alt","depth","freq")

## rename column headers
columns<-c("pos","ref","reads_all","mismatches","deletions","insertions","A","C","T","G")
ntcounts_c<-ntcounts_c[columns]
colnames(ntcounts_c)<-c("position","reference","coverage","mismatches","deletions","insertions","A","C","T","G")

## generate absolute SNV columns
ntcounts_c$A.SNV<-as.numeric(ifelse(ntcounts_c$reference=="A",0,
                                    ifelse(ntcounts_c$reference=="C",ntcounts_c$A,
                                           ifelse(ntcounts_c$reference=="G",ntcounts_c$A,
                                                  ifelse(ntcounts_c$reference=="T",ntcounts_c$A,
                                                         ifelse(ntcounts_c$reference=="",0,"na"))))))

ntcounts_c$C.SNV<-as.numeric(ifelse(ntcounts_c$reference=="C",0,
                                    ifelse(ntcounts_c$reference=="A",ntcounts_c$C,
                                           ifelse(ntcounts_c$reference=="G",ntcounts_c$C,
                                                  ifelse(ntcounts_c$reference=="T",ntcounts_c$C,
                                                         ifelse(ntcounts_c$reference=="",0,"na"))))))

ntcounts_c$G.SNV<-as.numeric(ifelse(ntcounts_c$reference=="G",0,
                                    ifelse(ntcounts_c$reference=="A",ntcounts_c$G,
                                           ifelse(ntcounts_c$reference=="C",ntcounts_c$G,
                                                  ifelse(ntcounts_c$reference=="T",ntcounts_c$G,
                                                         ifelse(ntcounts_c$reference=="",0,"na"))))))

ntcounts_c$T.SNV<-as.numeric(ifelse(ntcounts_c$reference=="T",0,
                                    ifelse(ntcounts_c$reference=="A",ntcounts_c$T,
                                           ifelse(ntcounts_c$reference=="C",ntcounts_c$T,
                                                  ifelse(ntcounts_c$reference=="G",ntcounts_c$T,
                                                         ifelse(ntcounts_c$reference=="",0,"na"))))))

## generate SNV frequency columns
ntcounts_c$snv.freq<-(ntcounts_c$mismatches/ntcounts_c$coverage)
ntcounts_c$snvfreqA<-(ntcounts_c$A.SNV/ntcounts_c$coverage)
ntcounts_c$snvfreqT<-(ntcounts_c$T.SNV/ntcounts_c$coverage)
ntcounts_c$snvfreqC<-(ntcounts_c$C.SNV/ntcounts_c$coverage)
ntcounts_c$snvfreqG<-(ntcounts_c$G.SNV/ntcounts_c$coverage)

## generate Squared Deviation column
ntcounts_c$sq.dev<-((ntcounts_c$snv.freq)^2)

## generate unique SNV column
ntcounts_c$unique.A.SNV<-ifelse(ntcounts_c$A.SNV>0,1,0)
ntcounts_c$unique.C.SNV<-ifelse(ntcounts_c$C.SNV>0,1,0)
ntcounts_c$unique.G.SNV<-ifelse(ntcounts_c$G.SNV>0,1,0)
ntcounts_c$unique.T.SNV<-ifelse(ntcounts_c$T.SNV>0,1,0)

## subset data for coverage >= 300
highcov<-subset(ntcounts_c,ntcounts_c$coverage>=300)
highcov$snv.freq<-((highcov$A.SNV + highcov$C.SNV + highcov$G.SNV + highcov$T.SNV)/highcov$coverage)
highcov$sq.dev<-((highcov$snv.freq)^2)

## subset highcov for positions called by lofreq
variant<-subset(highcov, highcov$position %in% var$position)

## calculate root mean squared deviation
rmsd<-sqrt((sum(variant$sq.dev)/nrow(highcov)))

## calculate shannon entropy
nucvar<-subset(variant,select=c("A","C","G","T"))
nucvar<-as.matrix(nucvar)

entropy.nucvar<-diversity(nucvar, type="entropy")
entropy.mean<-sum(entropy.nucvar$entropy)/nrow(highcov)
entropy.sd<-sd(entropy.nucvar$entropy)

## calculate Gini-Simpson index
ginisimpson<-diversity(nucvar, type="gini-simpson")
gs.mean<-sum(ginisimpson$gini.simpson.C)/nrow(highcov)
ginisimpson$mean.gs.C<-gs.mean
gs.sd<-sd(ginisimpson$gini.simpson.C)

## calculate mutation frequency (SNV/10,000 nt sequenced)
mut.freq<-(sum(variant$A.SNV,variant$C.SNV,variant$G.SNV,variant$T.SNV)/sum(highcov$coverage))*10000
unique.mut.freq<-(sum(variant$unique.A.SNV,variant$unique.C.SNV,variant$unique.G.SNV,variant$unique.T.SNV)/sum(highcov$coverage))*10000

## Calculate coverage metrics
percentcovered<-(nrow(highcov)/10675*100)
meandepth<-mean(highcov$coverage)

## merge diversity indices and write file
diversity<-c(percentcovered,meandepth,mut.freq,unique.mut.freq,rmsd,entropy.mean,gs.mean)
diversity<-data.frame(t(diversity))
colnames(diversity)<-c("%Covered","MeanDepth","Mut.Freq.per.10K","Unique.Mut.Freq","RMSD","Mean.Shannon.Entropy", "Mean.Gini.Simpson")
outputdiversity<-paste("R_", samplename,"_coverage_and_diversity_metrics.csv",sep="")
write.csv(diversity, file = outputdiversity, row.names = FALSE)

######     Mutational Spectrum      #######
## ntcounts_r column headers
ntcounts_r<-paste("13-ntcounts_reference",samplename,sep="_")
ntcounts_r<-paste(table,".csv",sep="")
ntcounts_r<-paste(homedir,table,sep="/")
ntcounts_r<-read.csv(table, header=TRUE)

columns<-c("pos","ref","reads_all","mismatches","deletions","insertions","A","C","T","G")
ntcounts_r<-ntcounts_r[columns]
colnames(ntcounts_r)<-c("position","reference","coverage","mismatches","deletions","insertions","A","C","T","G")


## generate absolute SNV columns
ntcounts_r$A.SNV<-as.numeric(ifelse(ntcounts_r$reference=="A",0,
                                    ifelse(ntcounts_r$reference=="C",ntcounts_r$A,
                                           ifelse(ntcounts_r$reference=="G",ntcounts_r$A,
                                                  ifelse(ntcounts_r$reference=="T",ntcounts_r$A,
                                                         ifelse(ntcounts_r$reference=="",0,"na"))))))

ntcounts_r$C.SNV<-as.numeric(ifelse(ntcounts_r$reference=="C",0,
                                    ifelse(ntcounts_r$reference=="A",ntcounts_r$C,
                                           ifelse(ntcounts_r$reference=="G",ntcounts_r$C,
                                                  ifelse(ntcounts_r$reference=="T",ntcounts_r$C,
                                                         ifelse(ntcounts_r$reference=="",0,"na"))))))

ntcounts_r$G.SNV<-as.numeric(ifelse(ntcounts_r$reference=="G",0,
                                    ifelse(ntcounts_r$reference=="A",ntcounts_r$G,
                                           ifelse(ntcounts_r$reference=="C",ntcounts_r$G,
                                                  ifelse(ntcounts_r$reference=="T",ntcounts_r$G,
                                                         ifelse(ntcounts_r$reference=="",0,"na"))))))

ntcounts_r$T.SNV<-as.numeric(ifelse(ntcounts_r$reference=="T",0,
                                    ifelse(ntcounts_r$reference=="A",ntcounts_r$T,
                                           ifelse(ntcounts_r$reference=="C",ntcounts_r$T,
                                                  ifelse(ntcounts_r$reference=="G",ntcounts_r$T,
                                                         ifelse(ntcounts_r$reference=="",0,"na"))))))

highcov_r<-subset(ntcounts_r,ntcounts_r$coverage>299)

## Calculate frequency of specific point mutations
highcov_r$AtoC<-as.numeric(ifelse(highcov_r$reference=="A",highcov_r$C.SNV,0))
highcov_r$AtoG<-as.numeric(ifelse(highcov_r$reference=="A",highcov_r$G.SNV,0))
highcov_r$AtoT<-as.numeric(ifelse(highcov_r$reference=="A",highcov_r$T.SNV,0))
highcov_r$CtoA<-as.numeric(ifelse(highcov_r$reference=="C",highcov_r$A.SNV,0))
highcov_r$CtoG<-as.numeric(ifelse(highcov_r$reference=="C",highcov_r$G.SNV,0))
highcov_r$CtoT<-as.numeric(ifelse(highcov_r$reference=="C",highcov_r$T.SNV,0))
highcov_r$GtoA<-as.numeric(ifelse(highcov_r$reference=="G",highcov_r$A.SNV,0))
highcov_r$GtoC<-as.numeric(ifelse(highcov_r$reference=="G",highcov_r$C.SNV,0))
highcov_r$GtoT<-as.numeric(ifelse(highcov_r$reference=="G",highcov_r$T.SNV,0))
highcov_r$TtoA<-as.numeric(ifelse(highcov_r$reference=="T",highcov_r$A.SNV,0))
highcov_r$TtoC<-as.numeric(ifelse(highcov_r$reference=="T",highcov_r$C.SNV,0))
highcov_r$TtoG<-as.numeric(ifelse(highcov_r$reference=="T",highcov_r$G.SNV,0))

AtoCfreq<-sum(highcov_r$AtoC)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
AtoGfreq<-sum(highcov_r$AtoG)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
AtoTfreq<-sum(highcov_r$AtoT)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
CtoAfreq<-sum(highcov_r$CtoA)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
CtoGfreq<-sum(highcov_r$CtoG)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
CtoTfreq<-sum(highcov_r$CtoT)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
GtoAfreq<-sum(highcov_r$GtoA)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
GtoCfreq<-sum(highcov_r$GtoC)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
GtoTfreq<-sum(highcov_r$GtoT)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
TtoAfreq<-sum(highcov_r$TtoA)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
TtoCfreq<-sum(highcov_r$TtoC)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
TtoGfreq<-sum(highcov_r$TtoG)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))

## Calculate mutation frequencies for specific subsititutions
highcov_r$Acoverage<-as.numeric(ifelse(highcov_r$reference=="A",highcov_r$coverage,0))
highcov_r$Ccoverage<-as.numeric(ifelse(highcov_r$reference=="C",highcov_r$coverage,0))
highcov_r$Gcoverage<-as.numeric(ifelse(highcov_r$reference=="G",highcov_r$coverage,0))
highcov_r$Tcoverage<-as.numeric(ifelse(highcov_r$reference=="T",highcov_r$coverage,0))

AtoCmutfreq<-sum(highcov_r$AtoC)/(sum(highcov_r$Acoverage))
AtoGmutfreq<-sum(highcov_r$AtoG)/(sum(highcov_r$Acoverage))
AtoTmutfreq<-sum(highcov_r$AtoT)/(sum(highcov_r$Acoverage))
CtoAmutfreq<-sum(highcov_r$CtoA)/(sum(highcov_r$Ccoverage))
CtoGmutfreq<-sum(highcov_r$CtoG)/(sum(highcov_r$Ccoverage))
CtoTmutfreq<-sum(highcov_r$CtoT)/(sum(highcov_r$Ccoverage))
GtoAmutfreq<-sum(highcov_r$GtoA)/(sum(highcov_r$Gcoverage))
GtoCmutfreq<-sum(highcov_r$GtoC)/(sum(highcov_r$Gcoverage))
GtoTmutfreq<-sum(highcov_r$GtoT)/(sum(highcov_r$Gcoverage))
TtoAmutfreq<-sum(highcov_r$TtoA)/(sum(highcov_r$Tcoverage))
TtoCmutfreq<-sum(highcov_r$TtoC)/(sum(highcov_r$Tcoverage))
TtoGmutfreq<-sum(highcov_r$TtoG)/(sum(highcov_r$Tcoverage))

## merge mutation spectrum and write file
spectrum<-c(AtoTfreq,AtoCfreq,AtoGfreq,CtoAfreq,CtoGfreq,CtoTfreq,GtoAfreq,GtoCfreq,GtoTfreq,TtoAfreq,TtoCfreq,TtoGfreq)
spectrum<-t(spectrum)
colnames(spectrum)<-c("AtoU","AtoC","AtoG","CtoA","CtoG","CtoU","GtoA","GtoC","GtoU","UtoA","UtoC","UtoG")
spectrum<-as.data.frame(spectrum)
csvfile<-paste("R_", samplename,"_mutational_spectrum.csv",sep="")
write.csv(spectrum, file = csvfile, row.names = FALSE)

## merge mutation spectrum frequencies and write file
spectrum_freq<-c(AtoTmutfreq,AtoCmutfreq,AtoGmutfreq,CtoAmutfreq,CtoGmutfreq,CtoTmutfreq,GtoAmutfreq,GtoCmutfreq,GtoTmutfreq,TtoAmutfreq,TtoCmutfreq,TtoGmutfreq)
spectrum_freq<-t(spectrum_freq)
colnames(spectrum_freq)<-c("A>U","A>C","A>G","C>A","C>G","C>U","G>A","G>C","G>U","U>A","U>C","U>G")
spectrum_freq<-as.data.frame(spectrum_freq)
csvfile<-paste("R_", samplename,"_mutational_spectrum_frequencies.csv",sep="")
write.csv(spectrum_freq, file = csvfile, row.names = FALSE)