#!/bin/bash

prefix=02_QC_passed/01_correct_IDs/
ANN=annotation/atRTD3_TS_21Feb22_transfix.gtf
OUTDIR=04_mapped
mkdir -p $OUTDIR

for fastq in 02_QC_passed/01_correct_IDs/*.gz; 
do
  outname="${fastq#$prefix}"
  outname="${outname%.fastq.gz}"
  mkdir -p $OUTDIR/$outname

  STAR \
  --genomeDir STARindex \
  --readFilesIn $fastq \
  --readFilesCommand zcat \
  --runThreadN 12 \
  --outReadsUnmapped None \
  --sjdbOverhang 91 \
  --sjdbGTFfile $ANN \
  --alignIntronMin 11 \
  --alignIntronMax 28000 \
  --outFilterMismatchNoverReadLmax 0.04 \
  --outFilterMismatchNmax 999 \
  --outFilterMultimapNmax 1 \
  --outSJfilterReads Unique \
  --outFileNamePrefix $OUTDIR/$outname/ \
  --outSAMtype BAM SortedByCoordinate \
  --outStd BAM_SortedByCoordinate \
  --alignEndsType Extend5pOfRead1 > $OUTDIR/$outname".bam"
done

