#!/bin/sh

#rat
#downloading the genome and sequence files
wget ftp://ftp.ensembl.org/pub/release-102/fasta/rattus_norvegicus//dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
gunzip Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-102/gtf/rattus_norvegicus//Rattus_norvegicus.Rnor_6.0.102.gtf.gz
gunzip Rattus_norvegicus.Rnor_6.0.102.gtf.gz

qsub rat.sh
qsub rat_htseq.sh
