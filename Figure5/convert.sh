#This is the source code to reproduce Figure 5.

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.ensGene.gtf.gz
gunzip hg19.ensGene.gtf.gz
nohup find . -name './files/*.Human.*.bam' | parallel -j 64 'htseq-count -s reverse -f bam {} ./hg19.ensGene.gtf > ./results/{.}.txt'

wget ftp://ftp.ensembl.org/pub/release-101/gtf/mus_musculus/Mus_musculus.GRCm38.101.gtf.gz
gunzip Mus_musculus.GRCm38.101.gtf.gz
nohup find . -name './files/*.Mouse.*.bam' | parallel -j 64 'htseq-count -s reverse -f bam {} ./Mus_musculus.GRCm38.101.gtf > ./results/{.}.txt'
