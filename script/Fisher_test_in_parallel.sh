mkdir -p fisher_test/methyC_CpG
mkdir -p fisher_test/methyC_CHG
mkdir -p fisher_test/methyC_CHH

parallel -a tables/methyC_CpG_counts_filter.tsv --block 1M --header : --pipepart "cat |grep -v '^$' > fisher_test/methyC_CpG/split_{#}.tsv"
parallel -a tables/methyC_CHG_counts_filter.tsv --block 1M --header : --pipepart "cat |grep -v '^$' > fisher_test/methyC_CHG/split_{#}.tsv"
parallel -a tables/methyC_CHH_counts_filter.tsv --block 1M --header : --pipepart "cat |grep -v '^$' > fisher_test/methyC_CHH/split_{#}.tsv"

find fisher_test/methyC_*/split_*.tsv | parallel --gnu 'Rscript script/fisher_test_DMC.R {} {}.RData'

Rscript script/combine_fisher_parallel.R

Rscript script/find_DMR.R

mkdir tables/anno
rm tables/*.tmp.*
find tables/*DMR*|./script/rush -k 'Rscript script/DMR_anno.R {} tables/anno/{%}'
