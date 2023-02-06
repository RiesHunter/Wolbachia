### Euclidean distance ###
## Import barcode counts for sample 1
sample1<-read.csv("/PATH/TO/BARCODE_COUNTS_SAMPLE1.csv", header=TRUE, sep=",")

## Import barcode counts for sample 2
sample2<-read.csv("/PATH/TO/BARCODE_COUNTS_SAMPLE2.csv", header=TRUE, sep=",")

## Merge barcode counts for the two samples
sample1and2<-merge(sample1, sample2, by = "barcode", all.x = TRUE, all.y = TRUE)
sample1and2[is.na(sample1and2)]<-0

## Generate column of the squared difference in frequencies
sample1and2$freqdiffsq<-(sample1and2$cumulativefreq.x-sample1and2$cumulativefreq.y)^2

## Sum the squared differences
sumfreqdiffsq<-sum(sample1and2$freqdiffsq)

## Take sq root of total differences squared to calculate Euclidean distance
eucdistance<-sqrt(sumfreqdiffsq)
eucdistance
