#!/bin/bash

# download reference and annotation
mkdir -p reference
mkdir -p annotation

wget -P reference/ https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz
wget -P annotation/ https://ics.hutton.ac.uk/atRTD/RTD3/atRTD3_TS_21Feb22_transfix.gtf
wget -P annotation/ https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-56/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.56.gff3.gz

# unpack GFF3 file
gunzip annotation/*

# prepare reference files
gunzip reference/*
mv reference/TAIR10_chr_all.fas reference/TAIR10_chr_all.fa

samtools faidx reference/TAIR10_chr_all.fa
cut -f1,2 reference/TAIR10_chr_all.fa.fai > reference/TAIR10_sizes.dat

bioawk -c fastx '{gsub(/[YWMKSRD]/,N,$seq);print ">"$name"\n"$seq}' reference/TAIR10_chr_all.fa > reference/TAIR10_chr_all_canonical.fa

# download sources
mkdir -p sources

wget -P sources/ https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
chmod ug+x sources/bedGraphToBigWig
wget -P sources/ https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.2/sratoolkit.3.0.2-ubuntu64.tar.gz
tar -zxvf sources/sratoolkit.3.0.2-ubuntu64.tar.gz

# download sample fastq via SRA toolkit with accession SRR24391474 
mkdir -p 00_raw_fastq
sources/sratoolkit.3.0.2-ubuntu64/bin/fastq-dump SRR24391474 --gzip

mv SRR24391474.fastq.gz 00_raw_fastq/AtGRP7-GFP.fastq.gz




