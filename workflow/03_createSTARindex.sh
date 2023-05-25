#!/bin/bash

REF=reference/TAIR10_chr_all.fa
ANN=annotation/atRTD3_TS_21Feb22_transfix.gtf


mkdir -p STARindex
STAR --runMode genomeGenerate --runThreadN 2 --genomeDir STARindex --genomeFastaFiles $REF --sjdbGTFfile $ANN --sjdbOverhang 91 --genomeSAindexNbases 12 --outFileNamePrefix 03_

