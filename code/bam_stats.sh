## TIME START
SECONDS=0


## CALCULATE TOTAL N BAMS
declare -i x=0
for reads in *.bam
do
x=$(( x + 1 ))
done
declare -i n=0


## Sample stats for average coverage on segments
printf "FileGene\tAvgCov\tSitesg0\tSites\tPerCov\n" > header

for reads in *.bam
do
f=${reads%%.bam}
samtools sort ${f}.bam > ${f}_sort_temp.bam
samtools index ${f}_sort_temp.bam
samtools depth -a -d 0 ${f}_sort_temp.bam > temp

# NP
grep A.*NP ./temp > ./temp_gene
fgene=$(printf ${f}""_NP)
avgcov=$(awk '{ total += $3; count++ } END { print total/count }' temp_gene)
sitestotal=$(awk '{ total += $3; count++ } END { print count }' temp_gene)
sitesgzero=$(awk '$3>0{c++} END {print c+0}' temp_gene)
percov=$(awk -v a=$sitesgzero -v b=$sitestotal 'BEGIN { print ( a / b ) }')
(printf $fgene"\t"$avgcov"\t"$sitesgzero"\t"$sitestotal"\t"$percov"\n") >> header
echo -ne "[#-------]\r"

# MP
grep A.*MP ./temp > ./temp_gene
fgene=$(printf ${f}""_MP)
avgcov=$(awk '{ total += $3; count++ } END { print total/count }' temp_gene)
sitestotal=$(awk '{ total += $3; count++ } END { print count }' temp_gene)
sitesgzero=$(awk '$3>0{c++} END {print c+0}' temp_gene)
percov=$(awk -v a=$sitesgzero -v b=$sitestotal 'BEGIN { print ( a / b ) }')
(printf $fgene"\t"$avgcov"\t"$sitesgzero"\t"$sitestotal"\t"$percov"\n") >> header
echo -ne "[##------]\r"

# PA
grep A.*PA ./temp > ./temp_gene
fgene=$(printf ${f}""_PA)
avgcov=$(awk '{ total += $3; count++ } END { print total/count }' temp_gene)
sitestotal=$(awk '{ total += $3; count++ } END { print count }' temp_gene)
sitesgzero=$(awk '$3>0{c++} END {print c+0}' temp_gene)
percov=$(awk -v a=$sitesgzero -v b=$sitestotal 'BEGIN { print ( a / b ) }')
(printf $fgene"\t"$avgcov"\t"$sitesgzero"\t"$sitestotal"\t"$percov"\n") >> header
echo -ne "[###-----]\r"

# HA
grep A.*HA ./temp > ./temp_gene
fgene=$(printf ${f}""_HA)
avgcov=$(awk '{ total += $3; count++ } END { print total/count }' temp_gene)
sitestotal=$(awk '{ total += $3; count++ } END { print count }' temp_gene)
sitesgzero=$(awk '$3>0{c++} END {print c+0}' temp_gene)
percov=$(awk -v a=$sitesgzero -v b=$sitestotal 'BEGIN { print ( a / b ) }')
(printf $fgene"\t"$avgcov"\t"$sitesgzero"\t"$sitestotal"\t"$percov"\n") >> header
echo -ne "[####----]\r"

# NA
grep A.*NA ./temp > ./temp_gene
fgene=$(printf ${f}""_NA)
avgcov=$(awk '{ total += $3; count++ } END { print total/count }' temp_gene)
sitestotal=$(awk '{ total += $3; count++ } END { print count }' temp_gene)
sitesgzero=$(awk '$3>0{c++} END {print c+0}' temp_gene)
percov=$(awk -v a=$sitesgzero -v b=$sitestotal 'BEGIN { print ( a / b ) }')
(printf $fgene"\t"$avgcov"\t"$sitesgzero"\t"$sitestotal"\t"$percov"\n") >> header
echo -ne "[#####---]\r"

# NS
grep A.*NS ./temp > ./temp_gene
fgene=$(printf ${f}""_NS)
avgcov=$(awk '{ total += $3; count++ } END { print total/count }' temp_gene)
sitestotal=$(awk '{ total += $3; count++ } END { print count }' temp_gene)
sitesgzero=$(awk '$3>0{c++} END {print c+0}' temp_gene)
percov=$(awk -v a=$sitesgzero -v b=$sitestotal 'BEGIN { print ( a / b ) }')
(printf $fgene"\t"$avgcov"\t"$sitesgzero"\t"$sitestotal"\t"$percov"\n") >> header
echo -ne "[######--]\r"

# PB1
grep A.*PB1 ./temp > ./temp_gene
fgene=$(printf ${f}""_PB1)
avgcov=$(awk '{ total += $3; count++ } END { print total/count }' temp_gene)
sitestotal=$(awk '{ total += $3; count++ } END { print count }' temp_gene)
sitesgzero=$(awk '$3>0{c++} END {print c+0}' temp_gene)
percov=$(awk -v a=$sitesgzero -v b=$sitestotal 'BEGIN { print ( a / b ) }')
(printf $fgene"\t"$avgcov"\t"$sitesgzero"\t"$sitestotal"\t"$percov"\n") >> header
echo -ne "[#######-]\r"

# PB2
grep A.*PB2 ./temp > ./temp_gene
fgene=$(printf ${f}""_PB2)
avgcov=$(awk '{ total += $3; count++ } END { print total/count }' temp_gene)
sitestotal=$(awk '{ total += $3; count++ } END { print count }' temp_gene)
sitesgzero=$(awk '$3>0{c++} END {print c+0}' temp_gene)
percov=$(awk -v a=$sitesgzero -v b=$sitestotal 'BEGIN { print ( a / b ) }')
(printf $fgene"\t"$avgcov"\t"$sitesgzero"\t"$sitestotal"\t"$percov"\n") >> header
echo -ne "[########]\r"

rm ${f}_sort_temp.bam; rm *.bai; rm temp; rm temp_gene

n=$(( n + 1 ))
echo "[$n/$x]: $f"
done
mv header bam_stats.tsv


## REPORT TIME             
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
echo "----- Done -----"
echo ""


# Reads that didn't map at 1 or more segment
grep "\t0\t\t" bam_stats.tsv | cut -f1 -d "_" | sort -u > to_remove

# Reads with less than 200x coverage at 1 or more segment
awk 'FNR==1 || ($2<200)' bam_stats.tsv | cut -f1 -d "_" | sort -u >> to_remove

# Reads with coverage of less than 99% of full genome
cat bam_stats.tsv > bam_stats.temp
sed -i '' '/\t0\t\t/d' bam_stats.temp
awk '{sub(/_.*/,"",$1)} 1' bam_stats.temp | awk '{ a[$1]+=$5 }END{ for(i in a) print a[i],i }' | awk '{$1=$1/8} 1' | awk '{ print $2 " " $1 }' > percent_coverage.tsv
awk -F" " '$2<0.99' percent_coverage.tsv | sed 's/\|/ /' | awk '{print $1}' >> to_remove

# HA coverage
cat bam_stats.tsv | grep HA > HA_coverage.tsv
awk 'FNR==1 || ($2>200)' HA_coverage.tsv | cut -f1 -d "_" | sort -u > HA_CoverageG200
awk 'FNR==1 || ($2<200)' HA_coverage.tsv | cut -f1 -d "_" | sort -u > HA_CoverageL200
awk 'FNR==1 || ($2>100)' HA_coverage.tsv | cut -f1 -d "_" | sort -u > HA_CoverageG100 
awk 'FNR==1 || ($2<100)' HA_coverage.tsv | cut -f1 -d "_" | sort -u > HA_CoverageL100

# Clean up
sort -u to_remove > to_remove.txt
rm to_remove bam_stats.temp
