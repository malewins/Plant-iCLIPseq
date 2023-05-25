#!/bin/bash

OUTDIR=01_correct_IDs
mkdir -p $OUTDIR
for FQ in *.gz;
do
  bioawk -c fastx '{print "@"$name"#"substr($comment,8) "\n"$seq"\n+\n"$qual}' $FQ | gzip >$OUTDIR/$FQ
done

