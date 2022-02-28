import pysam
import sys
import random
import argparse

def GetArgs():

	def ParseArgs(parser):
		class Parser(argparse.ArgumentParser):
			def error(self, message):
				sys.stderr.write('error: %s\n' % message)
				self.print_help()
				sys.exit(2)

		parser.add_argument('-b', '--bam',
							type=str,
							required=True,
							help='Path to bam output from Cell Ranger. Look for .../outs/possorted_genome_bam.bam')
		parser.add_argument('-1', '--read1_range',
							required=False,
							default="L1HS:1-20",
							type=str,
							help='Range in which to look for read 1s. Default is for human. See readme for mouse instructions [L1HS:1-20]')
		parser.add_argument('-2', '--read2_range',
							required=False,
							default="L1HS",
							type=str,
							help='Range in which to look for read 2s. Default is for human. See readme for mouse instructions [L1HS]')
		parser.add_argument('-e', '--max_edit_distance',
							required=False,
							default=2,
							type=int,
							help='Maximum allowed edit distance between read sequence and LINE-1 consensus [2]')
		parser.add_argument('-c', '--clip_max',
							required=False,
							default=20,
							type=int,
							help="Maximum allowable clipping from 5' end of read 1 (to account for TSO). No other clipping is allowed. [20]")
		parser.add_argument('-d', '--downsample_to',
							required=False,
							default=1.0,
							type=float,
							help="Downsample reads to this fraction of the total. If you are planning to look for differential expression, we recommend that you use the aggr program in Cell Ranger and then input the amoung of downsampling down by aggr here. Grep for frac_reads_kept in .../outs/summary.json. [1.0]")
		return parser.parse_args()

	parser = argparse.ArgumentParser()
	args = ParseArgs(parser)

	return args.bam, args.read1_range, args.read2_range, args.max_edit_distance, args.clip_max, args.downsample_to

def main():
	bamfilename, read1_range, read2_range, max_edit_distance, clip_max, downsample_to = GetArgs()
	
	bam = pysam.AlignmentFile(bamfilename,'rb')

	r1_UMI_dict = dict()
	for read in bam.fetch(region=read1_range):
		if random.random() < downsample_to and read.is_read1 and read.has_tag("nM") and read.get_tag("nM") <= max_edit_distance and read.has_tag("CB") and read.has_tag("UB") and not read.is_reverse and read.cigarstring[-1]!='S' and 'N' not in read.cigarstring:
			if read.cigartuples[0][0] != 4 or (read.cigartuples[0][1] <= clip_max):
				CB = read.get_tag("CB")
				UB = read.get_tag("UB")
				if CB not in r1_UMI_dict:
					r1_UMI_dict[CB] = set()
				r1_UMI_dict[CB].add(UB)

	r2_UMI_dict = dict()
	for read in bam.fetch(region=read2_range):
		if read.is_read2 and read.has_tag("nM") and read.get_tag("nM") <= max_edit_distance and read.has_tag("CB") and read.has_tag("UB") and read.is_reverse and 'S' not in read.cigarstring and 'N' not in read.cigarstring:
			CB = read.get_tag("CB")
			UB = read.get_tag("UB")
			if CB not in r2_UMI_dict:
				r2_UMI_dict[CB] = set()
			r2_UMI_dict[CB].add(UB)
			
	UMI_dict = dict()
	for CB in r1_UMI_dict:
		if CB in r2_UMI_dict:
			UMI_dict[CB] = r1_UMI_dict[CB].intersection(r2_UMI_dict[CB])

	for CB in UMI_dict:
		print (CB+'\t'+str(len(UMI_dict[CB])))

if __name__ == '__main__':
	main()
