#!/bin/bash

mkdir -p beds
for BAM in *.bam;
do
  bedtools bamtobed -i $BAM > beds/"${BAM%.bam}.bed"
done

cd beds
for BED in *.bed;
do
  awk 'BEGIN{OFS="\t"}{if($6=="+") print($1,($2-1),$2,$4,1,$6); else print($1,($3),($3+1),$4,1,$6)}' $BED | sort -T . -k1,1 -k2,2n > ${BED%.bed}.xlsite.bed
done

OUTBED=01_bedgraphs_strands
mkdir -p $OUTBED

for BED in *.xlsite.bed
do
  for STRND in + -
  do
    bedtools sort -i $BED | awk -v strand=$STRND '$6==strand{print}' | bedtools merge -i - -s -d -1 -c 5 -o sum > $OUTBED/"${BED%.xlsite.bed}".$STRND.bedgraph
  done
done

OUTBW=01_bigwigs
mkdir -p $OUTBW

cd $OUTBED

TAIR=../../../../reference/TAIR10_sizes.dat
BGTOBW=../../../../sources/bedGraphToBigWig

for BG in *.bedgraph
do
  $BGTOBW $BG $TAIR ../$OUTBW/${BG%.bedgraph}.bw
done

