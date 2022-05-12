
from collections import defaultdict
from modules import REGIONS
import logging
from pandas import read_csv
from os import path
# from line_profiler_pycharm import profile


dir_root = path.dirname(path.dirname(__file__))
logging.config.fileConfig(path.join(dir_root, 'logging_config.ini'), disable_existing_loggers=False)
logger = logging.getLogger(__name__)


class SEGMENT(REGIONS.REGION):

	_list_samples = []

	def __init__(self, sample, contig, start, end, numprobes, mean, name="segment", comment=""):
		"""
		Define how the SEGMENTS in a SEG file is. Seg file such as a segment outputted by GATK CNV tool.
		SEGMENTS is an extension of the class REGION because a SEGMENT share the same type of coordinates as a REGION
		:param sample: sample name usually present in column 1 of seg file
		:type sample: string
		:param contig: contig name (column 2 of SEG file)
		:type contig: string
		:param start:  start position of segment (column 3 of SEG file)
		:type start: int
		:param end: stop position of segment (column 4 of SEG file)
		:type end: int
		:param numprobes: number of probes in the segment (column 5 of SEG file)
		:type numprobes: integer
		:param mean: mean log2 ratio of segment (column 6 of SEG file)
		:type mean: float
		:param name: as we deal with copy number regions we call them segments as the extension of SEG file means ".seg" in GATK's header file; segment is a key word used from
		REGION class; Do not change this name unless necessary
		:type name: string
		:param comment: any comment to be added to the instance
		:type comment: string
		"""
		super().__init__(contig, start, end, name, comment)
		self.sample = sample
		self.numprobes = int(numprobes)
		self.mean = float(mean)
		self.flag = "notset"  # Not in use yet;
		self.index = "notset"  # Not in use yet;
		self.ratioNR = self.numprobes / self.slen  # Ratio that also could be use later if researchers are interested in

		self.list_the_samples()  # allows to check if the SEG file contains only and only one sample

	def print_obj(self):
		logger.info('\t'.join([self.contig, str(self.start), str(self.end), str(self.numprobes), str(self.mean), str(self.slen), self.flag, self.index, str(self.ratioNR)]))

	def print_segment(self):
		logger.info('\t'.join([str(x) for x in [self.contig, self.start, self.end, self.numprobes, self.mean, self.slen, self.flag, self.index, self.ratioNR]]))

	def __str__(self):
		return '\t'.join([str(x) for x in [self.sample, self.contig, self.start, self.end, self.numprobes, self.mean, self.slen]])

	# def __repr__(self):
	# 	# If you change the representation of the print for the object, it will not indicate the type anymore; default is for instance <modules.REGION.REGION object at 0x108fb4190>
	#  	return '\t'.join([str(x) for x in [self.contig, str(self.start), str(self.end), str(self.numprobes), str(self.mean), str(self.slen), self.flag, self.index, self.ratioNR]])

	def list_the_samples(self):
		"""
		make a unique list of samples found in the SEG file
		:return: a self list of samples
		:rtype: None
		"""
		if self.sample not in self._list_samples:
			self._list_samples.append(self.sample)

	def get_list_samples(self):
		return self._list_samples

	def copy(self):
		return SEGMENT(self.sample, self.contig, self.start, self.end, self.numprobes, self.mean)

	def get_distance_btw_segments(self, seg2):
		return super().get_distance_btw_regions(seg2)
	
	@staticmethod
	def mean_log2_to_copy_number(log2_value):
		return round(2 * 2 ** log2_value, 2)

	def get_copy_number_diff_btw_segments(self, prev_seg, th_min_diff_copy_number=0.5):
		if isinstance(prev_seg, SEGMENT):
			if self.contig == prev_seg.contig:
				if abs(self.mean_log2_to_copy_number(prev_seg.mean) - self.mean_log2_to_copy_number(self.mean)) >= th_min_diff_copy_number:
					logger.info("abs({} - {}) = DIFF Copy Number == {} >=? {}".format(str(self.mean_log2_to_copy_number(prev_seg.mean)),
					                                                                  str(self.mean_log2_to_copy_number(self.mean)),
					                                                                  str(abs(self.mean_log2_to_copy_number(prev_seg.mean)-self.mean_log2_to_copy_number(self.mean))),
					                                                                  str(th_min_diff_copy_number)
					                                                                  )
					            )
					return True
				return False
			else:
				logger.warning("Regions provided are on different contig; Regions are {} versus {}; Returning distance between the two regions is impossible; returning None value instead".format(
					self.__str__, prev_seg.__str__))
				return False
		return NotImplemented


# @profile
def read_and_sort_segment_file(fseg, delim=" ", list_exclude_contigs=['chrX', 'chrY', 'chrM'], list_keep_contigs=[]):
	"""
	read the GATK CNV's SEG file and sort the segment by chromosome, start, end (but not necessarily in natural sort) using pandas datafarme
	:param fseg: SEG file
	:type fseg: file
	:param delim: change the delimiter if different from space
	:type delim: basestring
	:param list_exclude_contigs: list of contigs present in the SEG file but need to be excluded from the process ; default is 3 contigs found in humans
	:type list_exclude_contigs: list
	:param list_keep_contigs: if a subset of all the contigs present in the SEG file need to be used for processing; default is ALL contig except the excluded ones
	:type list_keep_contigs: list
	:return: dictionary of contigs as keys and list of SEGMENTS as values
	:rtype: defaultdict
	"""
	dd = defaultdict(list)
	df = read_csv(fseg, sep=delim, skipinitialspace=True, header=0, comment='#')
	df.sort_values(["Chromosome", "Start", "End"], inplace=True)
	counter = 0
	logger.debug("list of contigs to keep given by User: " + str(list_keep_contigs) )
	for line in df.values:
		sample, contig, start, stop, numprobes, mean, *extra = line
		if contig in list_exclude_contigs:
			continue
		if list_keep_contigs != [] and contig not in list_keep_contigs:
			continue
		sego = SEGMENT(sample, contig, start, stop, numprobes, mean)
		counter += 1
		dd[contig].append(sego)  # dictionary of contig with a list of ordered segments
	logger.info("total number of segments to process after reading seg file: "+str(counter))
	return dd


def check_count_samples_in_seg_file(dico_of_segments):
	"""
	We want to capture the list_of_samples variable common to ALL the SEGMENTS instances.
	We need the dico if segments to get access to at least one SEGMENT instance and then get access to the list_of_segments
	## NOTE this is not the pythonic way but it works :-(
	:param dico_of_segments: dictionary of list of segments dict[contig] = list_of_SEGMENTS
	:type dico_of_segments: dict
	:param list_of_samples: list of sample names
	:type list_of_samples: list
	:return: the samplename is it is unique, otherwise ERROR
	:rtype: string
	"""
	for lsegs in dico_of_segments.values():
		for seg in lsegs:
			if len(seg.get_list_samples()) == 1:
				return seg.get_list_samples()[0]
			else:
				logger.error(seg.get_list_samples())
				raise ValueError("ERROR: more than one sample has been found in SEG file; Aborting. Check your input")


def are_segments_overlapping(obj_segment1, obj_segment2):
	"""
	Look whether given segments overlap or not and if so by how much
	example of returned value: (True, 76693100, 100.0, 80.39) ;
	- True, means the segments (or regions) overlap by at least Xbp
	- the number of common bases between segment 1 and 2 is 76693100 (Xbp is 76693100)
	- 100.0 means that 100% of common base pairs overlap segment 1, aka, segment1 overlaps @100% segment2
	- 80.39 means that 76693100 base from common bases represent 80.39% of the segment2.
	:param obj_segment1: segment represented in the seg file we want to know if it overlaps a known segment such as centromere or telomere
	:type obj_segment1: object of class SEGMENT or REGION
	:param obj_segment2: should be the segment we check if we overlap or not such as a centromere or telomere or q or q arm
	:type obj_segment2: object of class REGION or SEGMENT
	:return: tuple with 4 information: yes or no for overlapping segment2, number of overlapping bases in common, \
	and if yes overlapping by how much percentage of segment1, and % of segment 2 overlapped
	:rtype: tuple
	"""
	# ovlp_bp = max(0, min(max1, max2) - max(min1, min2)) # generic statement
	if obj_segment1.contig != obj_segment2.contig:
		raise Exception("contigs DO NOT match; found {} and {}; Check your input".format(obj_segment1.contig, obj_segment2.contig))
	ovlp_bp = max(0, min(obj_segment1.end, obj_segment2.end) - max(obj_segment1.start, obj_segment2.start))
	if obj_segment1.slen == 0 and obj_segment2.slen == 0:
		return False if ovlp_bp == 0 else True, ovlp_bp, 0, 0
	elif obj_segment1.slen == 0 and obj_segment2.slen != 0:
		return False if ovlp_bp == 0 else True, ovlp_bp, 0, ovlp_bp / obj_segment2.slen * 100
	elif obj_segment1.slen != 0 and obj_segment2.slen == 0:
		return False if ovlp_bp == 0 else True, ovlp_bp, ovlp_bp / obj_segment1.slen * 100, 0
	else:
		return False if ovlp_bp == 0 else True, ovlp_bp, ovlp_bp / obj_segment1.slen * 100, ovlp_bp / obj_segment2.slen * 100


def get_min_max_segments_per_contig(dsegs):
	"""
	capture the tuples of min and max position of the segment in their respective contig
	The returned dictionary will be used to re-make the "genoem territory"
	:param dsegs: dictionary of all the segments found in the SEG file
	:type dsegs: defaultdict
	:return: dictionary of tuples per contig; one tpl per contig is expected
	:rtype: dict
	"""
	dico_min_max_segs_per_contig = dict()

	for k, lsegs in dsegs.items():
		min = 10 ** 10
		max = 0
		for seg in lsegs:
			if seg.start < min:
				min = seg.start
			if seg.end > max:
				max = seg.end
		dico_min_max_segs_per_contig[k] = (min, max)

	return dico_min_max_segs_per_contig
