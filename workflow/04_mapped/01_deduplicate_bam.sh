#!/bin/bash

LOGFOLDER=01_dedup_logs
OUTFOLDER=01_deduped_BAMs
mkdir -p $LOGFOLDER
mkdir -p $OUTFOLDER
for BAM in *.bam;
do
  OUTFILE=$BAM
  umi_tools dedup -I $BAM \
  -L $LOGFOLDER/$OUTFILE.log \
  -S $OUTFOLDER/$OUTFILE \
  --extract-umi-method read_id \
  --method unique \
  --umi-separator="#"
done

