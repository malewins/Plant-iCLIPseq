#!/bin/bash

OUTDIR=00_initial_FQC

mkdir -p $OUTDIR

fastqc -o $OUTDIR *.gz
