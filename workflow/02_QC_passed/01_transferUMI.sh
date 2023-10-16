#!/bin/bash

OUTDIR=01_correct_IDs
mkdir -p $OUTDIR
for FQ in *.gz;
do
  bioawk -c fastx '{split($comment, UMI,"_"); print "@"$name"#"UMI[2] "\n"$seq"\n+\n"$qual}' $FQ | gzip >$OUTDIR/$FQ
done

