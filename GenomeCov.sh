### Create a file with mean reads Depth Coverage and % of bp covering the bin
# Eduardo Gonzalez: eduardo_gonzalez@eva.mpg.de
#Require bedtools v2.30.0

### Inputs:
	#1 Bed file indictaing chr size ex:
	#		1	249250621
	#		2	243199373

	#2 Sample Bam file
	#3 Bin Size (BP)

#How to run:
#	sh geneCov.sh <Hg19.chrSize.bed> <Sample_genome.bam> <10000>

out=${2%.bam}w$3.genCov
### Counting DP
bedtools genomecov -ibam $2 -bg > $2.0.tmp
	# Filtering chr
	awk '($1+0)>0 && ($1+0)<23' $2.0.tmp > $2.0.1.tmp
	awk '$1 = "X"' $2.0.tmp > $2.0.2.tmp
 	awk '$1 = "Y"' $2.0.tmp > $2.0.3.tmp

	cat $2.0.1.tmp $2.0.2.tmp $2.0.3.tmp |sed 's/\s/\t/g' | sort -V -k1,1 -k2,2n > $2.1.tmp 

## Generate Bins
bedtools makewindows -b $1 -w $3 > $2.2.tmp

## Intersecting bins to DP cov files
bedtools intersect -a $2.2.tmp -b $2.1.tmp  -wao | cut -f1,2,3,8 | gawk '{print$1"_"$2"_"$3"\t"$4}' |sort -V -k1,1 > $2.3.tmp


## Counting total bp covered in bin
awk 'NR==1{print;next} {for (i=2;i<=NF;i++) {a[$1][i]+=$i}} END{ for (j in a) {s=j; for (i=2;i<=NF;i++) {s=s" "a[j][i]}; print s}}' $2.3.tmp|sort -V -k1,1 > $2.4.tmp

## Estimating mean of DP by Bins
bedtools map -a $2.2.tmp -b $2.1.tmp -c 4 -o mean -null "0" | gawk '{print$1"_"$2"_"$3"\t"$4}' > $2.5.tmp

### Joining y generando output Final
join -1 1 -2 1 $2.4.tmp $2.5.tmp | sed 's/_/,/g'|sed 's/\s/,/g' > $out

rm $2*tmp



