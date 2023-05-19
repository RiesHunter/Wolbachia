#!/bin/sh
## Start time
SECONDS=0
echo ""; echo "-------------------- Starting barcode analyses... --------------------"

### ---------------------------------------------------------------------------------------------------- ###
## Parameters
workingdir="/Users/rieshunter/Google Drive/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/Wolbachia/data/reads/data/run/test_r2"
reference="/Users/rieshunter/Google Drive/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/Wolbachia/data/refseqs/PRVABC59.fasta"

## Index
touch _stderr.txt
bwa index ${reference} 2>>_stderr.txt
samtools faidx ${reference} 2>>_stderr.txt

## Working directory
cd "$workingdir"
pwd
ls -lh

### ---------------------------------------------------------------------------------------------------- ###
## Calculate number of pairs
declare -i x=0
for pairs in "$workingdir"/00_raw/*_L001_R2_001.fastq.gz
do
x=$(( x + 1 ))
done
declare -i n=0

### ---------------------------------------------------------------------------------------------------- ###
## For-loop
cd "$workingdir"/00_raw
for pairs in *_L001_R2_001.fastq.gz
do
cd "$workingdir"

# Names
sample=${pairs%%_L001_R2_001.fastq.gz}
n=$(( n + 1 ))
echo "[$n/$x]: ${sample}"

## Align processed paired-end reads to reference with default settings
bwa mem -t 4 "$reference" ./06_norm/06-norm_${sample}_L001_r2.fastq > ./06_norm/relaxed_norm_${sample}.bam  2>>_stderr.txt

### sort and index alignment
samtools sort ./06_norm/relaxed_norm_${sample}.bam > ./06_norm/sorted_relaxed_norm_${sample}.bam 2>>_stderr.txt
samtools index ./06_norm/sorted_relaxed_norm_${sample}.bam #2>>_stderr.txt

### extract aligned reads overlapping the barcode region
samtools view -h ./06_norm/sorted_relaxed_norm_${sample}.bam "KU501215.1:4007-4030" > ./06_norm/BCregion_sorted_relaxed_norm_${sample}.bam 2>>_stderr.txt
 # for some reason, I had to change "KU501215_1:4007-4030" to "KU501215.1:4007-4030" 

## convert aligned reads to fastq
bamToFastq -i ./06_norm/BCregion_sorted_relaxed_norm_${sample}.bam \
  -fq ./06_norm/BCregion_sorted_relaxed_${sample}.fastq 2>>_stderr.txt

## convert fastq to fasta
seqtk seq -a ./06_norm/BCregion_sorted_relaxed_${sample}.fastq > ./06_norm/BCregion_sorted_relaxed_${sample}.fasta 2>>_stderr.txt

## extract just the barcode sequence from all the reads
sed -n '/^>/p; /CT.GC.GC.CT.AC.CC.CT.GC./p' ./06_norm/BCregion_sorted_relaxed_${sample}.fasta > ./06_norm/sed_BCregion_sorted_relaxed_${sample}.fasta 2>>_stderr.txt
sed -e '$!N;/^>.*\n>/D' -e 'P;D' ./06_norm/sed_BCregion_sorted_relaxed_${sample}.fasta > ./06_norm/sedfilter_BCregion_sorted_relaxed_${sample}.fasta 2>>_stderr.txt
sed 's/^.*\(CT.GC.GC.CT.AC.CC.CT.GC.\).*$/\1/' ./06_norm/sedfilter_BCregion_sorted_relaxed_${sample}.fasta > ./06_norm/BConly_${sample}.fasta 2>>_stderr.txt

## Count kmers
kmercountexact.sh \
  in=./06_norm/BConly_${sample}.fasta \
  out=./06_norm/kmer_counts_${sample}.fasta \
  k=24 khist=./06_norm/kmer_freq_hist_${sample}.tsv \
  fastadump=t rcomp=f ow=t 2>>_stderr.txt

## Generate tsv for barcode counts
seqkit fx2tab ./06_norm/kmer_counts_${sample}.fasta -H > ./06_norm/counts_BConly_${sample}.tsv 2>>_stderr.txt
done

### ---------------------------------------------------------------------------------------------------- ###
## Clean up
mkdir ./15_barcode_analyses
mv ./06_norm/*relaxed*.fastq    ./15_barcode_analyses
mv ./06_norm/*.fasta            ./15_barcode_analyses
mv ./06_norm/*relaxed*.bam      ./15_barcode_analyses
mv ./06_norm/*relaxed*.bam.bai  ./15_barcode_analyses
mv ./06_norm/*.tsv              ./15_barcode_analyses

## Time end
echo ""; duration=$(($SECONDS / 60))
echo "$(($duration / 60)) hours and $(($duration)) minutes elapsed."; echo ""
echo "-------------------- Barcode analyses completed! --------------------"; echo ""