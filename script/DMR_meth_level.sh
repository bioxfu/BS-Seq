bedtools intersect -a $1 -b $2 -wa -wb|sort -k1,1 -k2,2n|groupBy -g 1,2,3 -c 5,6,4,7,8 -o min,max,count,sum,sum > $3
