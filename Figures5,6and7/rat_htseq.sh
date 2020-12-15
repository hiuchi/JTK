#!/bin/bash
#$ -S /bin/bash
#$ -t 1-294:1
#$ -l s_vmem=4G -l mem_req=4G
#$ -pe smp 1
#$ -N qsub_rat_htseq
#$ -cwd
#$ -o out.log
#$ -e error.log

files=(`find ~/TEG/201209/STAR_results/ -name \*.fastq.gz.Aligned.out.bam | cut -c 38-47`)
file=${files[$SGE_TASK_ID]}

htseq-count -s reverse -f bam ~/TEG/201209/STAR_results/$file.fastq.gz.Aligned.out.bam ./Rattus_norvegicus.Rnor_6.0.102.gtf > ~/TEG/201211/results/rat/$file.txt
