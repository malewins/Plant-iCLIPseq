#!/bin/bash

OUTDIR=01_barcode_QC
FILE_IN=AtGRP7-GFP.fastq.gz
BC_OUT=$OUTDIR/iclip.barcodes.qc
mkdir -p $OUTDIR
zcat $FILE_IN | awk -v umi1_len=3 -v exp_bc_len=4 '{if (FNR%4==2) print substr($1,(umi1_len+1),exp_bc_len)} ' | sort | uniq -c | sort -k1,1rn > $BC_OUT

# iCLIP2
#zcat $FILE_IN | awk -v umi1_len=5 -v exp_bc_len=6 '{if (FNR%4==2) print substr($1,(umi1_len+1),exp_bc_len)} ' | sort | uniq -c | sort -k1,1rn > $BC_OUT
