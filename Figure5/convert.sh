#This is the source code to reproduce Figure 5.

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.ensGene.gtf.gz
gunzip hg19.ensGene.gtf.gz
nohup find . -name './files/*.Human.*.bam' | parallel -j 64 'htseq-count -s reverse -f bam {} ./hg19.ensGene.gtf > ./results/{.}.txt'

wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ensGene.gtf.gz
gunzip mm10.ensGene.gtf.gz
nohup find . -name './files/*.Mouse.*.bam' | parallel -j 64 'htseq-count -s reverse -f bam {} ./mm10.ensGene.gtf > ./results/{.}.txt'
