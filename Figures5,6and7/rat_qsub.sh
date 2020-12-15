#!/bin/bash
#$ -S /bin/bash
#$ -t 1-295:1
#$ -l s_vmem=8G -l mem_req=8G
#$ -pe smp 32
#$ -N qsub_rat
#$ -cwd
#$ -o out.log
#$ -e error.log

files=(`find ./files -name \*.gz | cut -c 9-27`)
file=${files[$SGE_TASK_ID]}

STAR --runThreadN 32 \
	--genomeDir STAR_index \
	--readFilesCommand gunzip -c \
	--readFilesIn ./files/${file} \
	--quantMode TranscriptomeSAM GeneCounts \
	--outSAMtype BAM Unsorted \
	--outFileNamePrefix ./STAR_results/${file}.

rsem-calculate-expression -p 32 \
	--alignments \
	--append-names \
	--no-bam-output \
	--bam \
	./STAR_results/${file}.Aligned.toTranscriptome.out.bam \
	./RSEM_index/RSEM_index \
	./results/${file}
