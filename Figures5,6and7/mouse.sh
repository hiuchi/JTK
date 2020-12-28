#!/bin/sh

#mouse
#downloading the genome and sequence files
wget ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna//Mus_musculus.GRCm38.dna.toplevel.fa.gz
gunzip Mus_musculus.GRCm38.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus//Mus_musculus.GRCm38.102.gtf.gz
gunzip Mus_musculus.GRCm38.102.gtf.gz

qsub mouse.sh
qsub mouse_htseq.sh
