#!/bin/sh

#mouse
#downloading the genome and sequence files
#wget ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna//Mus_musculus.GRCm38.dna.toplevel.fa.gz
#gunzip Mus_musculus.GRCm38.dna.toplevel.fa.gz
#wget ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus//Mus_musculus.GRCm38.102.gtf.gz
#gunzip Mus_musculus.GRCm38.102.gtf.gz
#aria2c -i ../urls.txt -x5

#generating genome indexes
STAR --runMode genomeGenerate \
	--runThreadN 64 \
	--genomeDir STAR_index \
	--genomeFastaFiles Mus_musculus.GRCm38.dna.toplevel.fa \
	--sjdbGTFfile Mus_musculus.GRCm38.102.gtf \
	--limitGenomeGenerateRAM 53524373088

rsem-prepare-reference --gtf Mus_musculus.GRCm38.102.gtf \
	-p 64 \
	Mus_musculus.GRCm38.dna.toplevel.fa \
	./RSEM_index/RSEM_index

#mapping and calculating TPM
files=`find ./files -name \*.gz | cut -c 9-27`

for file in ${files[@]};
do
STAR --runThreadN 64 \
	--genomeDir STAR_index \
	--readFilesCommand gunzip -c \
	--readFilesIn ./files/${file} \
	--quantMode TranscriptomeSAM \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix ./STAR_result/${file}.

rsem-calculate-expression -p 64 \
	--alignments \
	--append-names \
	--no-bam-output \
	--bam \
	./STAR_result/${file}.Aligned.toTranscriptome.out.bam \
	./RSEM_index/RSEM_index \
	./results/${file}
done

nohup find . -name 'ERR2588370.*.bam' | parallel -j 32 'htseq-count -s reverse -f bam {} ../Mus_musculus.GRCm38.102.gtf > ../counts/{.}.txt'

