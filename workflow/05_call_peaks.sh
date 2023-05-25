#!/bin/bash

OUTDIR=05_called_peaks
mkdir -p $OUTDIR

BAM=04_mapped/01_deduped_BAMs/AtGRP7-GFP.bam
BAI=04_mapped/01_deduped_BAMs/AtGRP7-GFP.bam.bai
REF=reference/TAIR10_chr_all_canonical.fa
PEAKS=$OUTDIR/AtGRP7-GFP.bed

# AtGRP7-GFP samples
pureclip -iv 'Chr1;Chr2' -ld -i $BAM -bai $BAI -g $REF -o $PEAKS -nt 24 -vv > 05_pureclip.log
