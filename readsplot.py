import sys
import os
import pysam
import re
from optparse import OptionParser


class ReadPlot(object):
	def __init__(self, bam, outdir=os.getcwd(), readsgroup=None, qual_filter=1, rmdup=False,  **kwargs):
		self.outdir = outdir
		self.bam = bam
		self._bam = pysam.AlignmentFile(bam, "rb")
		self.Rmdup = rmdup
		self.qua_ctr = int(qual_filter)
		self.readsgroup = readsgroup
		self.leftpadding = int(kwargs["leftpadding"]) if "leftpadding" in kwargs else 3
		self.bgcolor = kwargs["bgcolor"] if "bgcolor" in kwargs else "#FFF"

	def __del__(self):
		self._bam.close()

	def get_png(self, region, prefix):
		chrom, tmp_r = region.split(":")
		start, stop = map(int, tmp_r.split("-"))
		reads = self.get_reads(chrom, start, stop)
		for read in reads:
                    print (read)


	@staticmethod
	def get_color(number):
		code = hex(max(60 - int(number), 0))[2:]
		return "#" + code * 3

	def get_reads(self, chrom, start, stop):

                # keep left most smaller than the first read start
                leftMostReadPos = 0
                rightMostReadPos = 0
                for aln_line in self._bam.fetch(chrom, start-1, start):
                    leftMostReadPos = aln_line.reference_start
                    alignment_start = aln_line.query_alignment_start
                    leftMostReadPos -= alignment_start 
                    break

                for aln_line in self._bam.fetch(chrom, stop - 1, stop):
                    rightMostReadPos = aln_line.reference_start + aln_line.query_length # actually,  it excat == match pos end + insert_seq, just cheat it 

                begin = leftMostReadPos - 1
		end = rightMostReadPos + 1 # actually no reads can more than stop, by fetch
		if stop == start:
			start -= 1
		refer = "N" * (end - begin + 1)
		total_reads = list()
		indel = self.check_indel(chrom, start, stop)
		refer_bases = ""
		for i in range(begin, end):
			refer_bases += refer[i - begin]
			if i in indel and indel[i] > 0:
				refer_bases += "-" * indel[i]
		total_reads.append(refer_bases)
		for aln_line in self._bam.fetch(chrom, start, stop):
			if self._filter_reads(aln_line) is True:
				continue
			cigar, pos_start = aln_line.cigartuples, aln_line.reference_start # ref pos
			stand = "-" if aln_line.is_reverse is True else "+"
			sequece = aln_line.query_sequence
			alignment_start = aln_line.query_alignment_start # this is the reads start pos, if clip 4 , start from 4, defalut 0
			pos_start -= alignment_start  
                        refIdx = 0
                        reads = []
                        refMatchLen = 0
                        gaplen = 0 # bug fixed, original version forget to add "-" to ref corrdinate, just add match an insert to POS 
                        while refMatchLen <  pos_start - begin:
                            if refer_bases[refIdx] == "-":
                                reads.append(" ")
                                gaplen += 1
                            else:
                                reads.append(" ")
                                refMatchLen += 1
                            refIdx += 1
                        
			pos = pos_start # ref pos
                        dele = 0
                        readLenIdx = 0
                        # For realign bam will meet bam begin with INS when origianl alignment is clip or match
			for cigar_terms in cigar:
				if cigar_terms[0] == 0:
					base, number = 0, dele +gaplen  # number is the pos of ref
					while base < cigar_terms[1]:
                                                if refer_bases[pos - begin + number ] == "-":
                                                    reads.append("-")
                                                    gaplen += 1
                                                else:
   						    reads.append(sequece[pos + base - pos_start])
						    base += 1
                                                    readLenIdx +=  1
						number += 1
					pos += cigar_terms[1]
				elif cigar_terms[0] == 1:
					for base in range(cigar_terms[1]):
						reads.append(sequece[pos - pos_start + base])
					pos += cigar_terms[1] # add insert shift
                                        readLenIdx +=  cigar_terms[1]
				elif cigar_terms[0] == 2:
					for base in range(cigar_terms[1]):
						reads.append("-")
 					dele += cigar_terms[1]
                                        readLenIdx +=  cigar_terms[1]
				else:
#					for base in range(cigar_terms[1]):
#						reads += sequece[pos - pos_start + base]
                                        reads.append("N" * cigar_terms[1])
					pos += cigar_terms[1]
                                        readLenIdx +=  cigar_terms[1]
			total_reads.append("".join(reads))
		return total_reads

	def _filter_reads(self, reads):
		if self.Rmdup and reads.is_duplicate:
			return True
		if reads.mapping_quality < self.qua_ctr or reads.is_qcfail:
			return True
		if self.readsgroup is not None and reads.get_tag('RG') != self.readsgroup:
			return True
		return False

	def check_indel(self, chrom, start, stop):
		indel = dict()
		for aln_line in self._bam.pileup(chrom, start, stop):
			indel[aln_line.reference_pos] = 0
			for pileup_read in aln_line.pileups:
				if self._filter_reads(pileup_read.alignment) is True:
					continue
				if pileup_read.indel > 0:
					indel[aln_line.reference_pos] = max(indel[aln_line.reference_pos], pileup_read.indel)
				else:
					continue
		return indel


def main():
	usage = "Usage: %prog --bam [bamfile] --region [region] [options]"
	description = "A program for plot the reads covered target position for BAM sequence alignment files."
	author = "Author: joeyhwong@hotmail.com"
	parser = OptionParser(usage=usage, version="%prog 0.1", description=description,
	                      epilog=author)
	parser.add_option('-b', '--bam', help='The sorted bam file (required)',
	                  dest='bamfile')
	parser.add_option('-r', '--region', help='The target region (1 base), example: '
	                                         'chr1:123456-123457', dest='region')
	parser.add_option('--prefix', help='String to add to the start of output graph',
	                  default=os.path.join(os.getcwd(), 'Test'), dest='prefix')
	parser.add_option('--rg', help='The read group you want to draw', default=None,
	                  dest='readsgroup')
	parser.add_option('-f', help='Only show reads which mapping quality > NUM',
	                  default=1, type='int', dest='quality')
	parser.add_option('-d', help='Remove duplication reads, cigar level dups',
	                  action="store_true", dest="Rmdup", default=False)
	(options, args) = parser.parse_args()

	if not options.bamfile or not options.region:
		parser.print_help()
		return 'Expected arguments lost !!!'
	bamfile = os.path.abspath(options.bamfile)
	regions = str(options.region)
	assert os.path.isfile(bamfile), "%s not exit !" % bamfile
	prefix = os.path.abspath(options.prefix).strip('_')
	outdir, prefix = os.path.split(prefix)
	rg = options.readsgroup
	qua_ctl = int(options.quality)
	rmdup = options.Rmdup
	test = ReadPlot(bamfile, outdir, rg, qua_ctl, rmdup)
	test.get_png(region=regions, prefix=prefix)


if __name__ == '__main__':
	sys.exit(main())
