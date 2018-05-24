#! /usr/bin/env python

import sys
import gzip

cov_cutoff = int(sys.argv[1])

with gzip.open(sys.argv[2], 'rb') as f:
	for line in f:
		line = line.decode('utf-8')
		lst = line.strip().split('\t')
		chrom = lst[0]
		pos = int(lst[1]) # 1-based
		strand = lst[2]
		methy_count = float(lst[3])
		unmethy_count = float(lst[4])
		total_count = methy_count + unmethy_count
		if total_count >= cov_cutoff:
			ratio = methy_count / total_count
			if strand == '-':
				ratio = -1 * ratio
			print('%s\t%s\t%s\t%.3f' % (chrom, pos-1, pos, ratio))
