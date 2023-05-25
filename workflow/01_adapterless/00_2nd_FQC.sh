#!/bin/bash

OUTDIR=00_2nd_FQC

mkdir -p $OUTDIR

fastqc -o $OUTDIR -t 3 *.gz
