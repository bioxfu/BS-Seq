#! /usr/bin/env python

import sys
import gzip

output_file = open(sys.argv[2], 'w')
output_file.write('track type=bedGraph\n')

with gzip.open(sys.argv[1], 'rb') as f:
	for line in f:
		line = line.decode('utf-8')
		lst = line.strip().split('\t')
		chrom = lst[0]
		pos = int(lst[1])
		strand = lst[2]
		methy_count = float(lst[3])
		unmethy_count = float(lst[4])
		total_count = methy_count + unmethy_count
		if total_count > 0 and methy_count > 0:
			ratio = methy_count / total_count
			if strand == '-':
				ratio = -1 * ratio
			output_file.write('%s\t%s\t%s\t%.3f\n' % (chrom, pos, pos+1, ratio))

output_file.close()
