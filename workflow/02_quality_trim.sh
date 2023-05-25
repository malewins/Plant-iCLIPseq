#!/bin/bash
OUTDIR=02_QC_passed
mkdir -p $OUTDIR
for FQ in 01_adapterless/*.gz;
do
  OUTPREFIX="${FQ#01_adapterless/}"
  OUTPREFIX="${OUTPREFIX%.fastq.gz}"
  flexbar -r $FQ  --zip-output GZ -t $OUTDIR/$OUTPREFIX  -q WIN -qf sanger  -qt 24 --min-read-length 15 -n 6
done

