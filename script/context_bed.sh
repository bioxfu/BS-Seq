zcat $1 |awk '{print $1"\t"$2-1"\t"$2"\t"$4"\t"$5}' > $2
