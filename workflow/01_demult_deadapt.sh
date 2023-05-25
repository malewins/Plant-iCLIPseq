#!/bin/bash

INFILE=00_raw_fastq/AtGRP7-GFP.fastq.gz
OUTDIR=01_adapterless/
mkdir -p $OUTDIR
flexbar -r $INFILE -t $OUTDIR -n 6 -O $OUTDIR/flexbar.log --zip-output GZ --barcodes 00_raw_fastq/barcodes-rc.fa --barcode-unassigned --barcode-trim-end LTAIL --barcode-error-rate 0 --adapter-seq AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG --adapter-trim-end RIGHT --adapter-error-rate 0.1 --adapter-min-overlap 1 --min-read-length 15 --umi-tags

