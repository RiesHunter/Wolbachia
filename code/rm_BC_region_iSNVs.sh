#!/bin/sh
## Start time
SECONDS=0
echo ""; echo "-------------------- Starting barcode removal... --------------------"

### ---------------------------------------------------------------------------------------------------- ###
## Calculate number of pairs
declare -i x=0
for pairs in *.vcf
do
x=$(( x + 1 ))
done
declare -i n=0

## Remove BC iSNVs from reads
for reads in *.vcf
do
# sample name
sample=${reads%%.vcf}
# sample exclusion
bcftools filter \
 -o ${sample}_noBC.vcf \
 -e 'POS>=4007 && POS<=4030' \
 ${sample}.vcf

bcftools filter \
 -o ${sample}_noBC_no5155.vcf \
 -e 'POS=5155' \
 ${sample}_noBC.vcf

bcftools filter \
 -o ${sample}_noBC_no5155_no1964.vcf \
 -e 'POS=1964' \
 ${sample}_noBC_no5155.vcf

 bcftools filter \
 -o ${sample}_cleaned.vcf \
 -e 'POS=8488' \
 ${sample}_noBC_no5155_no1964.vcf

# add one more to the pile!
n=n+1
echo""; echo "-------------------- [$n/$x]: $sample --------------------"
done


### ---------------------------------------------------------------------------------------------------- ###
## Time end
echo ""; duration=$(($SECONDS / 60))
echo "$(($duration / 60)) hours and $(($duration)) minutes elapsed."; echo ""
echo "-------------------- Barcode removal completed! --------------------"; echo ""
