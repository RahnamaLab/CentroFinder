# Example usage: python3 trf2bed.py --dat hg19.fa.2.7.7.80.10.24.6.dat --bed hg19.fa.2.7.7.80.10.24.6_gangstr.bed --tool gangstr

from argparse import (ArgumentParser, FileType)

def parse_args():
	parser = ArgumentParser(description='Convert Tandem Repeat Finder (TRF) dat file to bed format with repeat units for genotyping STRs.')
	parser.add_argument('--dat', type=str, required=True, help='Input dat file produced by Tandem Repeat Finder (TRF) using the -d option.')
	parser.add_argument('--bed', type=str, required=True, help='Output bed file based on Tandem Repeat Finder (TRF) data.')
	parser.add_argument('--tool', type=str, required=True, choices=['lobstr', 'gangstr', 'hipstr', 'repeatseq', 'gatk'], help='Name of the tool for what the file bed file will be generated.')

	return parser.parse_args()

def main():
	args = parse_args()
	datfile = args.dat
	bedfile = args.bed
	toolname = args.tool


	with open(bedfile, 'w') as bed:
		print("Writing bed file.")
		chrom = ""
		with open(datfile, 'r') as dat:
			print("Writing dat file.")
			for line in dat:
				splitline = line.split()
				if line.startswith("Sequence:"):
					chrom = line.split()[1]
				else:
					try:
						try:
							int(splitline[0])
						except ValueError:
							continue
						start = splitline[0]
						end = splitline[1]
						motif_length = splitline[2]
						reference_length = splitline[3]
						motif = splitline[13]
						alignment_score = splitline[7]
						none = '.'

						if toolname == 'lobstr':
							bed.write('\t'.join([chrom, start, end, motif_length, reference_length, none, none, none, alignment_score, none, none, none, none, none, motif]) + '\n') # for LobSTR

						elif toolname == 'gangstr':
							bed.write('\t'.join([chrom, start, end, motif_length, motif]) + '\n') # for GangSTR

						elif toolname == 'hipstr':
							bed.write('\t'.join([chrom, start, end, motif_length, reference_length]) + '\n') # for HipSTR

						elif toolname == 'repeatseq':
							bed.write(chrom + ':' + start + '-' + end + '\t' + motif_length + '_'.join([reference_length, motif_length, start, end, alignment_score, splitline[8], splitline[9], splitline[10], splitline[11], splitline[12], motif]) + '\n') # for RepeatSeq

						elif toolname == 'gatk':
							bed.write(chrom + '\t' + start + '\t' + end + '\n') # for GATK

						else:
							print("Tool not specified")
							print("Unable to write files.")

					except IndexError:
						pass
			print("Dat file creation successful.")
		print("Bed file creation successful.")

if __name__ == '__main__':
	main()
