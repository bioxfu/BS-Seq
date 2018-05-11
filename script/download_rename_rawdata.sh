scp -r xfu@10.41.25.251:$1 ./
dir=`echo $1|sed -r 's/.*\///g'`
ls $dir/*/* |./script/rush -k 'mv {} {/%}_{%@.+_(R\d).fq.gz}.fastq.gz' --dry-run
ls $dir/*/* |./script/rush -k 'mv {} {/%}_{%@.+_(R\d).fq.gz}.fastq.gz'
mkdir fastq
mv *.fastq.gz fastq

