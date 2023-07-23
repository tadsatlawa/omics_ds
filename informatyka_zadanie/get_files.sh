#!/bin/bash
wget "ftp://ftp.ncbi.nih.gov/refseq/D_rerio/mRNA_Prot/*faa.gz"
wget "ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/*faa.gz"
gzip -d *.gz
