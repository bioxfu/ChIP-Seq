import sys

tot = 0
patch = {}

input1 = sys.argv[1]
input2 = sys.argv[2]
output_file = sys.argv[3]

with open(input1) as f:
	for line in f:
		seqname, start, cnt = line.strip().split('\t')
		if seqname == 'TYLCV-Alm':
			tot += int(cnt)
		elif  seqname == 'TYLCV-Alm_3prime_140bp_5prime_140bp':
			tot += int(cnt)
			if int(start) >= 141:
				patch[int(start)-140] = cnt
			else:
				patch[int(start)-140+2781] = cnt

with open(input2) as f:
	for line in f:
		seqname, start, cnt = line.strip().split('\t')
		if seqname == 'TYLCV-Alm':
			tot += int(cnt)
		elif  seqname == 'TYLCV-Alm_3prime_140bp_5prime_140bp':
			tot += int(cnt)
			if int(start) >= 141:
				patch[int(start)-140] = cnt
			else:
				patch[int(start)-140+2781] = cnt

tot = float(tot) / 1000000

output = open(output_file, 'w')
output.write('seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\n')

with open(input1) as f:
	for line in f:
		seqname, start, cnt = line.strip().split('\t')
		if seqname == 'TYLCV-Alm':
			start = int(start)
			if start in patch:
				score = (float(cnt) + float(patch[start])) 
				score = score / tot
			else:
				score = float(cnt)/tot
			output.write('.\t.\t.\t%s\t%s\t%s\t+\t.\n' % (start, start, score))

with open(input2) as f:
	for line in f:
		seqname, start, cnt = line.strip().split('\t')
		if seqname == 'TYLCV-Alm':
			start = int(start)
			if start in patch:
				score = (float(cnt) + float(patch[start])) 
				score = score / tot
			else:
				score = float(cnt)/tot
			output.write('.\t.\t.\t%s\t%s\t-%s\t-\t.\n' % (start, start, score))

output.close()
