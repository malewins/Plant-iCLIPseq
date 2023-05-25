#!/bin/bash

samtools merge 01_deduped_BAMs/AtGRP7-GFP.bam 01_deduped_BAMs/AtGRP7-GFP_rep?.bam
samtools index 01_deduped_BAMs/AtGRP7-GFP.bam
