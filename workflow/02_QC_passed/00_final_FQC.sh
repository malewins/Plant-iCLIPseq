#!/bin/bash

OUTDIR=00_final_FQC

mkdir -p $OUTDIR

fastqc -o $OUTDIR -t 3 *.gz
