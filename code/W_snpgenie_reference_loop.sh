#!/bin/sh
## Start time
SECONDS=0
echo ""; echo "-------------------- Starting... --------------------"

### ---------------------------------------------------------------------------------------------------- ###
## Parameters
workingdir="/Users/rieshunter/Google Drive/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/Wolbachia/data/reads/data/run"
reference="/Users/rieshunter/Google Drive/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/Wolbachia/data/refseqs/PRVABC59.fasta"
gtf="/Users/rieshunter/Google Drive/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/Wolbachia/data/refseqs/PRVABC59.gtf"
bed="/Users/rieshunter/Google Drive/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/Wolbachia/data/refseqs/PRVABC59.bed"

#ls "${workingdir}"
#cat "${reference}"
#cat "${gtf}"
#cat "${bed}"

### ---------------------------------------------------------------------------------------------------- ###
## Calculate number of pairs
declare -i x=0
for pairs in *.vcf
do
x=$(( x + 1 ))
done
declare -i n=0

touch _stderr.txt

# Use SNPgenie to characterize variants and calculate pi
for reads in *.vcf
do
sample=${reads%%.vcf}
echo ${sample}
#cat ${sample}.vcf

snpgenie.pl \
  --vcfformat=2 \
  --snpreport=./${sample}.vcf \
  --fastafile="${reference}" \
  --gtffile="${gtf}" \
  --slidingwindow=30 2>>_stderr.txt  

chmod 755 SNPGenie_Results/product_results.txt
mv SNPGenie_Results/product_results.txt ./${sample}_sg_product_results.txt
mv SNPGenie_Results/codon_results.txt ./${sample}_sg_codon_results.txt
mv SNPGenie_Results/sliding_window_length* ./${sample}_sg_sliding_window.txt
rm -r SNPGenie_Results

n=n+1

echo""; echo "-------------------- [$n/$x]: $sample --------------------"
done


## organize
mkdir output
mv *_product_results.txt output
mv *_codon_results.txt output
mv *_sliding_window.txt output
cd output
mkdir sample_files



## compile product_results.txt
awk /file/ ${sample}*_product_results.txt > header
for reads in *_product_results.txt
do
f=${reads%%_product_results.txt}
# add each file's data to the header
awk /09-reference/ ${f}_product_results.txt >> header
mv ${f}_product_results.txt ./sample_files
done
mv header ./snpgenie_compiled.tsv

## compile codon_results.txt
awk /file/ *_codon_results.txt > header
head -n 1 header > cols
for reads in *_codon_results.txt
do
f=${reads%%_codon_results.txt}
# add each file's data to the header
awk /09-reference/ ${f}_codon_results.txt >> cols
mv ${f}_codon_results.txt ./sample_files
done
mv cols ./snpgenie_compiled_sw.tsv
rm header

## compile sliding_window.txt
awk /file/ *_sliding_window.txt > header
head -n 1 header > cols
for reads in *_sliding_window.txt 
do
f=${reads%%_sliding_window.txt} 
# add each file's data to the header
awk /09-reference/ ${f}_sliding_window.txt >> cols   
mv ${f}_sliding_window.txt ./sample_files 
done
mv cols ./snpgenie_compiled_sliding_window.tsv
rm header




## time end          
echo ""
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
echo ""
