## TIME START
SECONDS=0

## CALCULATE TOTAL N VCFS
declare -i x=0
for reads in *.vcf
do
x=$(( x + 1 ))
done
declare -i n=0

## RENAME VCFS
for reads in *.vcf
do
f=${reads%%.vcf}
sed  "s/PASS/${f}/" ${f}.vcf > ${f}_temp.vcf
sed '/\t\t/d' ${f}_temp.vcf > ${f}_fn.vcf
rm ${f}_temp.vcf
done
echo ""
echo "---------- Renaming complete ----------"

## REPORT TIME
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

## TIME START
SECONDS=0

## ECHO DATABASE OF CHOICE
echo ""
echo $1
echo ""

## SNPEFF WITH DATABASE OF CHOICE
for reads in *_fn.vcf
do
f=${reads%%.vcf}
java -jar snpEff.jar ann $1 ${f}.vcf > ${f}_ann.vcf
rm snpEff_summary.html
mv snpEff_genes.txt ${f}_snpEff_genes.tsv
sed -i '' '/# The following table is formatted as tab separated values./d' ./${f}_snpEff_genes.tsv
sed -i '' 's/#//g' ${f}_snpEff_genes.tsv
n=$(( n + 1 ))
echo "[$n/$x]: $f"
done

## ORGANIZE FILES
mkdir output
cd output
mkdir raw
mkdir fn
mkdir fn_ann
mkdir snpEff_genes
mv ../*_fn.vcf fn
mv ../*_fn_ann.vcf fn_ann
mv ../*_snpEff_genes.tsv snpEff_genes
mv ../*.vcf raw

echo ""
echo "---------- Completed -----------"

## REPORT TIME
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
