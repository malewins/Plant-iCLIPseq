#!/bin/bash

for bam in *.bam
do
  samtools index $bam
done
