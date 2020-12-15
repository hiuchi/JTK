#!/bin/sh

#rat
#downloading the genome and sequence files
wget ftp://ftp.ensembl.org/pub/release-102/fasta/rattus_norvegicus//dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
gunzip Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-102/gtf/rattus_norvegicus//Rattus_norvegicus.Rnor_6.0.102.gtf.gz
gunzip Rattus_norvegicus.Rnor_6.0.102.gtf.gz

#generating genome indexes
STAR --runMode genomeGenerate \
	--runThreadN 32 \
	--genomeDir STAR_index \
	--genomeFastaFiles Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
	--sjdbGTFfile Rattus_norvegicus.Rnor_6.0.102.gtf \
	--limitGenomeGenerateRAM 53524373088

rsem-prepare-reference --gtf Rattus_norvegicus.Rnor_6.0.102.gtf \
	-p 32 \
	Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
	./RSEM_index/RSEM_index

#mapping and calculating TPM
files=`find ./files -name \*.gz | cut -c 9-27`

for file in ${files[@]};
do
STAR --runThreadN 32 \
	--genomeDir STAR_index \
	--readFilesCommand gunzip -c \
	--readFilesIn ./files/${file} \
	--quantMode TranscriptomeSAM \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix ./STAR_result/${file}.

rsem-calculate-expression -p 32 \
	--alignments \
	--append-names \
	--no-bam-output \
	--bam \
	./STAR_result/${file}.Aligned.toTranscriptome.out.bam \
	./RSEM_index/RSEM_index \
	./results/${file}
done
