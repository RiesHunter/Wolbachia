#!/bin/sh
## Start time
SECONDS=0
echo ""; echo "-------------------- Starting... --------------------"


### ---------------------------------------------------------------------------------------------------- ###
## Parameters
workingdir="/Users/rieshunter/Documents/bioinformatics/Wolbachia/data/reads/data/run"
reference="/Users/rieshunter/Documents/bioinformatics/Wolbachia/data/refseqs/PRVABC59.fasta"
gtf="/Users/rieshunter/Documents/bioinformatics/Wolbachia/data/refseqs/PRVABC59.gtf"
bed="/Users/rieshunter/Documents/bioinformatics/Wolbachia/data/refseqs/PRVABC59.bed"
## Versions
bwa &> out_bwa; head -3 out_bwa | tail -2 > _versions.txt; rm out_bwa; echo "" >> _versions.txt
samtools --version &> out_sam; head -2 out_sam >> _versions.txt; rm out_sam; echo "" >> _versions.txt
cutadapt &> out_cut; head -1 out_cut >> _versions.txt; rm out_cut; echo "" >> _versions.txt
bbmerge.sh --version &> out_bbm; head -3 out_bbm | tail -1>> _versions.txt; rm out_bbm; echo "" >> _versions.txt
printf "lofreq " > out_lof; lofreq version &> out_lof2; cat out_lof2 >> out_lof; head -1 out_lof >> _versions.txt; rm out_lo*; echo "" >> _versions.txt
perl --version &> out_perl; head -2 out_perl | tail -1 >> _versions.txt; rm out_perl; echo "" >> _versions.txt
perl $HOME/Documents/GitHub/snpdat/SNPdat_v1.0.5.pl -v &> out_snpdat; head -2 out_snpdat | tail -1 >> _versions.txt; rm out_snpdat; echo "" >> _versions.txt
bcftools version &> out_bcf; head -1 out_bcf >> _versions.txt; rm out_bcf; echo "" >> _versions.txt
snpgenie.pl --version &> out_sg; head -1 out_sg >> _versions.txt; rm out_sg; echo "" >> _versions.txt
pysamstats --help &> out_pys; tail -2 out_pys | head -1 >> _versions.txt; rm out_pys; echo "" >> _versions.txt
Rscript --version &> out_Rs; head -1 out_Rs >> _versions.txt; rm out_Rs; echo "" >> _versions.txt
source ~/miniconda3/bin/activate; conda activate pysam; conda info >> _versions.txt; echo "" >> _versions.txt

### ---------------------------------------------------------------------------------------------------- ###
## Calculate number of pairs
declare -i x=0
for pairs in *_R1_001.fastq.gz
do
x=$(( x + 1 ))
done
declare -i n=0

### ---------------------------------------------------------------------------------------------------- ###
## _stderr.txt
touch _stderr.txt
## For-loop
for pairs in *_R1_001.fastq.gz
do
# Index
bwa index ${reference} 2>>_stderr.txt
samtools faidx ${reference} 2>>_stderr.txt

# Names
sample=${pairs%%_R1_001.fastq.gz}
n=$(( n + 1 ))
echo "[$n/$x]: ${sample}"

# Quality trim and adapter trim
cutadapt --quiet -j 4 \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -o 01-at_${sample}_r1_L1.fastq \
  -p 01-at_${sample}_r2_L1.fastq \
  ${sample}_R1_001.fastq.gz ${sample}_R2_001.fastq.gz 2>>_stderr.txt

# Trim non-overlapping regions, filter reads with mismatches, and merge reads
bbmerge.sh \
  in1=01-at_${sample}_r1_L1.fastq \
  in2=01-at_${sample}_r2_L1.fastq \
  out=02-merged_${sample}_r1_L1.fastq \
  out2=02-merged_${sample}_r2_L1.fastq \
  outu=02-unmerged_${sample}_r1_L1.fastq \
  outu2=02-unmerged_${sample}_r2_L1.fastq \
  minoverlap=50 \
  ecco=t \
  mix=f 2>>_stderr.txt

# Concatenate reads from lane 1 and lane 2
cat 02-merged_${sample}_r1_L1.fastq > 03-cat_${sample}_r1.fastq
cat 02-merged_${sample}_r2_L1.fastq > 03-cat_${sample}_r2.fastq

# Quality trim and adapter trim 
cutadapt -q 35,35 -m 50 --max-n 0 -j 4 \
 -o 04-qtatlf_${sample}_r1.fastq \
 -p 04-qtatlf_${sample}_r2.fastq \
 ./03-cat_${sample}_r1.fastq \
 ./03-cat_${sample}_r2.fastq 2>>_stderr.txt

# Primer trim reads
cutadapt --quiet -g XGACAGTTCGAGTTTGAAGCGAAAG -g XAAGAAAGATCTGGCTGCCATGC -g XAGATGACGTCGATTGTTGGTGC \
	-g XTCAGGTGCATAGGAGTCAGCAA -g XAGAACGTTAGTGGACAGAGGCT -g XTTGATTGTGAACCGAGGACAGG -g XTGAAGGGCGTGTCATACTCCTT \
	-g XGGGAGAAGAAGATCACCCACCA -g XGCCTTAGGGGGAGTGTTGATCT -g XACGGTCGTTGTGGGATCTGTAA -g XCAGCCGTTATTGGAACAGCTGT \
	-g XCACTAAGGTCCACGTGGAGGAA -g XTGGCAGTGCTGGTAGCTATGAT -g XCAATGGTTTTGCTTTGGCCTGG -g XCCCTAGCGAAGTACTCACAGCT \
	-g XGTGGCATGAACCCAATAGCCAT -g XGTGGTCCATGGAAGCTAGATGC -g XCTGTTGAGTGCTTCGAGCCTTC -g XTATGGATGAGGCCCACTTCACA \
	-g XGGCTGGAAAACGGGTCATACAG -g XAGAGACTGACGAAGACCATGCA -g XTGGACCAGACACGGAGAGAAAA -g XCGTCTTGATGAGGAACAAGGGC \
	-g XTAATGGGAAGGAGAGAGGAGGG -g XCCCTGACCCTAATAGTGGCCAT -g XACTGGAACTCCTCTACAGCCAC -g XAGTGCAAAGCTGAGATGGTTGG \
	-g XGGTGGGGGATTGGCTTGAAAAA -g XAGGATGTGAATCTCGGCTCTGG -g XAAAAGTGGACACTAGGGTGCCA -g XACAAGGGGAATTTGGAAAGGCC \
	-g XAAATGGAAAAAGGGCACAGGGC -g XCAAACGAATGGCAGTCAGTGGA -g XATTTCCACAGAAGGGACCTCCG -g XACCACCTGGGCTGAGAACATTA \
	-G XAGTATGCACTCCCACGTCTAGT -G XTGATTCCAACCAGGTTTGCGAC -G XTACGGTGACACAACCTCCATGT -G XGGAGCCATGAACTGACAGCATT \
	-G XTGTGCGTCCTTGAACTCTACCA -G XCCATCTGTCCCTGCGTACTGTA -G XCGCCTCCAACTGATCCAAAGTC -G XTTGACTGCTGCTGCCAATCTAC \
	-G XGAGTGGGCATTCCTTCAGTGTG -G XGTGGGACTTTGGCCATTCACAT -G XCCTGGGCCTTATCTCCATTCCA -G XTATCAGCGCCAGATGAGCTACA \
	-G XAGAGAGAGGAGCATAAACCCCC -G XTTTCCCATGTGATGTCACCTGC -G XTACACTCCATCTGTGGTCTCCC -G XGCTCCAATGTCCCCATCCTTTG \
	-G XCCTCTAAGGGCCTCCTCCATTT -G XTGGTGAGTTGGAGTCCGGAAAT -G XGCCATCAAGTATGACCGGCTTT -G XCCTTTGCTCCGTCCTAAGCTTG \
	-G XCTCCAAAAGCCGCTCCTCTTTT -G XATTCTGGCTGGCTCAATTTCCG -G XAAGTGGTCACTGCATGTTGGAC -G XTCTCCACTTGGGGGTCAATTGT \
	-G XCCTTCCATTTCTCTCCCAGGGT -G XACCAGGGCCTCCTTTTGTGTAT -G XATGTGTAGAGTTGCGGGAGAGT -G XGGGCCTCATAGCTTCCATGGTA \
	-G XATGCTGCATTGCTACGAACCTT -G XTAATCCCAGCCCTTCAACACCA -G XCGTAAGTGACAACTTGTCCGCT -G XTGTCCCATCCAGTTGAGGGTTT \
	-G XATCCACACTCTGTTCCACACCA -G XTGACTAGCAGGCCTGACAACAT -G XACCACTAGTCCCTCTTCTGGAG \
	-a CTTTCGCTTCAAACTCGAACTGTCX -a GCATGGCAGCCAGATCTTTCTTX -a GCACCAACAATCGACGTCATCTX -a TTGCTGACTCCTATGCACCTGAX \
	-a AGCCTCTGTCCACTAACGTTCTX -a CCTGTCCTCGGTTCACAATCAAX -a AAGGAGTATGACACGCCCTTCAX -a TGGTGGGTGATCTTCTTCTCCCX \
	-a AGATCAACACTCCCCCTAAGGCX -a TTACAGATCCCACAACGACCGTX -a ACAGCTGTTCCAATAACGGCTGX -a TTCCTCCACGTGGACCTTAGTGX \
	-a ATCATAGCTACCAGCACTGCCAX -a CCAGGCCAAAGCAAAACCATTGX -a AGCTGTGAGTACTTCGCTAGGGX -a ATGGCTATTGGGTTCATGCCACX \
	-a GCATCTAGCTTCCATGGACCACX -a GAAGGCTCGAAGCACTCAACAGX -a TGTGAAGTGGGCCTCATCCATAX -a CTGTATGACCCGTTTTCCAGCCX \
	-a TGCATGGTCTTCGTCAGTCTCTX -a TTTTCTCTCCGTGTCTGGTCCAX -a GCCCTTGTTCCTCATCAAGACGX -a CCCTCCTCTCTCCTTCCCATTAX \
	-a ATGGCCACTATTAGGGTCAGGGX -a GTGGCTGTAGAGGAGTTCCAGTX -a CCAACCATCTCAGCTTTGCACTX -a TTTTTCAAGCCAATCCCCCACCX \
	-a CCAGAGCCGAGATTCACATCCTX -a TGGCACCCTAGTGTCCACTTTTX -a GGCCTTTCCAAATTCCCCTTGTX -a GCCCTGTGCCCTTTTTCCATTTX \
	-a TCCACTGACTGCCATTCGTTTGX -a CGGAGGTCCCTTCTGTGGAAATX -a TAATGTTCTCAGCCCAGGTGGTX -A ACTAGACGTGGGAGTGCATACTX \
	-A GTCGCAAACCTGGTTGGAATCAX -A ACATGGAGGTTGTGTCACCGTAX -A AATGCTGTCAGTTCATGGCTCCX -A TGGTAGAGTTCAAGGACGCACAX \
	-A TACAGTACGCAGGGACAGATGGX -A GACTTTGGATCAGTTGGAGGCGX -A GTAGATTGGCAGCAGCAGTCAAX -A CACACTGAAGGAATGCCCACTCX \
	-A ATGTGAATGGCCAAAGTCCCACX -A TGGAATGGAGATAAGGCCCAGGX -A TGTAGCTCATCTGGCGCTGATAX -A GGGGGTTTATGCTCCTCTCTCTX \
	-A GCAGGTGACATCACATGGGAAAX -A GGGAGACCACAGATGGAGTGTAX -A CAAAGGATGGGGACATTGGAGCX -A AAATGGAGGAGGCCCTTAGAGGX \
	-A ATTTCCGGACTCCAACTCACCAX -A AAAGCCGGTCATACTTGATGGCX -A CAAGCTTAGGACGGAGCAAAGGX -A AAAAGAGGAGCGGCTTTTGGAGX \
	-A CGGAAATTGAGCCAGCCAGAATX -A GTCCAACATGCAGTGACCACTTX -A ACAATTGACCCCCAAGTGGAGAX -A ACCCTGGGAGAGAAATGGAAGGX \
	-A ATACACAAAAGGAGGCCCTGGTX -A ACTCTCCCGCAACTCTACACATX -A TACCATGGAAGCTATGAGGCCCX -A AAGGTTCGTAGCAATGCAGCATX \
	-A TGGTGTTGAAGGGCTGGGATTAX -A AGCGGACAAGTTGTCACTTACGX -A AAACCCTCAACTGGATGGGACAX -A TGGTGTGGAACAGAGTGTGGATX \
	-A ATGTTGTCAGGCCTGCTAGTCAX -A CTCCAGAAGAGGGACTAGTGGTX  \
	   -e 0.1 -O 6 -j 4 \
       -o 05-ptqtatlf_${sample}_r1.fastq \
       -p 05-ptqtatlf_${sample}_r2.fastq \
       ./04-qtatlf_${sample}_r1.fastq \
       ./04-qtatlf_${sample}_r2.fastq 2>>_stderr.txt

# Normalize depth with BBnorm
echo "Grab a coffee, some tea, or go for a brief jog..."
bbnorm.sh \
  in1=05-ptqtatlf_${sample}_r1.fastq \
  in2=05-ptqtatlf_${sample}_r2.fastq \
  out1=06-norm_${sample}_r1.fastq \
  out2=06-norm_${sample}_r2.fastq \
  target=2000 2>>_stderr.txt
echo "Normalization complete!"

# Align reads with BWA MEM
bwa mem -t 8 -B 10 ${reference} ./06-norm_${sample}_r1.fastq ./06-norm_${sample}_r2.fastq > 07-merged_norm_${sample}.bam 2>>_stderr.txt

# Sort sam by coordinate
samtools sort \
  -o ./08-sorted_norm_${sample}.bam \
  ./07-merged_norm_${sample}.bam 2>>_stderr.txt

# Call variants with lofreq
lofreq call \
  -l ${bed} \
  -q 30 -Q 30 -C 300 \
  -f ${reference} \
  -o 09-reference_${sample}.vcf \
  ./08-sorted_norm_${sample}.bam 2>>_stderr.txt

# Use SNPdat to annotate variants
perl $HOME/Documents/GitHub/snpdat/SNPdat_v1.0.5.pl \
  -i 09-reference_${sample}.vcf \
  -g ${gtf} \
  -f ${reference} \
  -s 10-snpdat_${sample}.txt \
  -o 10-snpdat_${sample}.tsv 2>>_stderr.txt

# Bgzip and index lofreq vcf
bgzip ./09-reference_${sample}.vcf 2>>_stderr.txt
tabix -p vcf ./09-reference_${sample}.vcf.gz 2>>_stderr.txt

# Extract consensus sequence for population
cat ${reference} | bcftools consensus -e "INFO/AF<0.50" ./09-reference_${sample}.vcf.gz > 11-consensus_${sample}.fasta 2>>_stderr.txt

# Call variants relative to consensus with lofreq
lofreq call \
  -l ${bed} \
  -q 30 -Q 30 -C 300 \
  -f ./11-consensus_${sample}.fasta \
  -o 09-consensus_${sample}.vcf \
  ./08-sorted_norm_${sample}.bam 2>>_stderr.txt

# Use SNPgenie to characterize variants and calculate pi
snpgenie.pl \
  --vcfformat=2 \
  --snpreport=./09-consensus_${sample}.vcf \
  --fastafile=./11-consensus_${sample}.fasta \
  --gtffile=${gtf} \
  --slidingwindow=30 2>>_stderr.txt
chmod 755 SNPGenie_Results/product_results.txt
mv SNPGenie_Results/product_results.txt ./12-${sample}_sg_product_results.txt
mv SNPGenie_Results/codon_results.txt ./12-${sample}_sg_codon_results.txt
mv SNPGenie_Results/sliding_window_length* ./12-${sample}_sg_sliding_window.txt
rm -r SNPGenie_Results

# Generate nucleotide counts per position
samtools index -b ./08-sorted_norm_${sample}.bam 2>>_stderr.txt

pysamstats -f 11-consensus_${sample}.fasta \
  -t variation \
  -D 1000000 \
  --format csv \
  --output 13-ntcounts_consensus_${sample}.csv \
  ./08-sorted_norm_${sample}.bam 2>>_stderr.txt
pysamstats -f ${reference} \
  -t variation \
  -D 1000000 \
  --format csv \
  --output 13-ntcounts_reference_${sample}.csv \
  ./08-sorted_norm_${sample}.bam 2>>_stderr.txt

# Calculate diversity metrics with R
Rscript /Users/rieshunter/Documents/bioinformatics/Wolbachia/code/1-process.R \
  ${workingdir} \
  ${sample} 2>>_stderr.txt

# Clean
rm *-ptqtatlf_*.fastq *-qtatlf_*.fastq *-cat*.fastq *norm*.bam *-unmerged*.fastq *-merged*.fastq *-at*.fastq *.bai *.fai *.tbi

echo""; echo "-------------------- Done with ${sample} --------------------"
done

### ---------------------------------------------------------------------------------------------------- ###
## Directories
printf "00_raw files\t\t"; ls *.fastq.gz | wc -l
printf "06_norm files\t\t"; ls 06-*.fastq | wc -l
printf "09_consensus_VCF files\t"; ls 09-con* | wc -l
printf "09_reference_VCF files\t"; ls 09-ref* | wc -l
printf "10_snpdat_TSV files\t"; ls 10-*.tsv | wc -l
printf "10_snpdat_TXT files\t"; ls 10-*.txt | wc -l
printf "11_consensus_FASTA files"; ls 11-*.fasta | wc -l
printf "12_snpgenie files\t"; ls 12-* | wc -l
printf "13_ntcounts files\t"; ls 13-* | wc -l
printf "14_R files\t\t"; ls R_* | wc -l

mkdir 00_raw; mv *.fastq.gz 00_raw
mkdir 06_norm; mv 06-* 06_norm
mkdir 09_consensus_VCF; mv 09-con* 09_consensus_VCF
mkdir 09_reference_VCF; mv 09-ref* 09_reference_VCF
mkdir 10_snpdat_TSV; mv 10-*tsv 10_snpdat_TSV
mkdir 10_snpdat_TXT; mv 10-*txt 10_snpdat_TXT
mkdir 11_consensus_FASTA; mv 11-* 11_consensus_FASTA
mkdir 12_snpgenie; mv 12-* 12_snpgenie
mkdir 13_ntcounts; mv 13-* 13_ntcounts
mkdir 14_R; mv R_* 14_R

### ---------------------------------------------------------------------------------------------------- ###
## Time end          
echo ""; duration=$(($SECONDS / 60)); echo $duration
echo "$(($duration / 60)) hours and $(($duration)) minutes elapsed."; echo ""
exit 1
