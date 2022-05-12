
import logging

from natsort import natsorted
from collections import defaultdict
from os import path


dir_root = path.dirname(path.dirname(__file__))
logging.config.fileConfig(path.join(dir_root, 'logging_config.ini'), disable_existing_loggers=False)
logger = logging.getLogger(__name__)


class REGION:
	"""
	# NOTE: the assumption is there is only one p arm, one q arm and one centromere per contig
	# NOTE: the assumption is there two telomeres per contig
	"""

	_ACCEPTED_REGION_NAMES = ("p", "q", "telomere", "centromere", "segment", "contig", "chromosome",  ".")

	def __init__(self, contig, start, end, name, comment=""):
		self.contig = contig
		self.start = int(start)
		self.end = int(end)
		self.name = name
		self.comment = comment
		self.slen = int(self.end - self.start)  # BASE HARDCODED

		# if we validate inputs
		# we need to validate the inputs when reading the line of the file; if len(line.split())<4 --> ERROR
		if name not in self._ACCEPTED_REGION_NAMES:
			raise Exception(
				"Unknown name of the region; Region name for Chromosome information MUST be one of the expected names as in {} [CASE SENSITIVE]".format(";".join(self._ACCEPTED_REGION_NAMES)))

	def print_region(self):
		logger.info('\t'.join([self.contig, str(self.start), str(self.end), str(self.name), str(self.slen)]))

	def __str__(self):
		return '\t'.join([str(x) for x in [self.contig, self.start, self.end, self.name, self.slen]])

	def get_distance_btw_regions(self, prev_seg):
		if isinstance(prev_seg, REGION):
			if self.contig == prev_seg.contig:
				if prev_seg.end <= self.start and prev_seg.end <= self.end:
					return abs(self.start - prev_seg.end)
				elif prev_seg.start >= self.end and prev_seg.start >= self.start:
					return abs(self.end - prev_seg.start)
				else:
					logger.warning("Regions provided overlaps; This is not Expected within a SEG file; Regions are {} versus {}; Returning distance between the two regions is possible but make no "
					               "sense in this current context of counting breaks ; returning None value instead".format(self.__str__, prev_seg.__str__))
					return 10**10
			else:
				logger.warning("Regions provided are on different contig; Regions are {} versus {}; Returning distance between the two regions is impossible; returning None value instead".format(
					self.__str__, prev_seg.__str__))
				return 10**10
		return NotImplemented



def read_file_genomic_regions(file_reg):
	"""
	read 4-columns or more tabulated file with genomic regions
	Name of the region can be one of these: p, q, centromere, telomere or a dot (.) if no name can be associated; The regions with dot will not be processed as for now
	:param file_reg: file with genomic regions; at least 4 columns such as col1=CHR or CONTIG, col2=START, col3=END, col4=NAME_of_the_REGION
	:type file_reg: file
	:return: dictionary
	:rtype: defaultdict
	"""
	dd = defaultdict(list)
	counter = 0
	try:
		with open(file_reg, 'r') as f:
			for line in f:
				counter += 1
				if line.startswith("#") or line == "":
					continue
				contig, start, stop, name = line.strip().split("\t")
				rego = REGION(contig, start, stop, name)
				dd[contig].append(rego)  # dictionary of contigs as Key with a list of REGION objects as Values
	except IOError as ioe:
		logger.error("ERROR: Your File may have some spaces; It must be all tabulated; check your input")
		logger.error(ioe)
	logger.info("number of regions captured from input file <{}>: {}".format(file_reg, str(counter)))
	return dd


def subset_dico_by_region(d, region_name):
	"""
	using a name of a region to extract the appropriate dictionary from a namedTuple
	:param d: dictionary of regions
	:type d: dict
	:param region_name: name of a genomic region such as telomer, centromere, q, p, chromosome, etc.
	:type region_name: string
	:return: dictionary of REGION representing region-name's REGIONS
	:rtype: defaultdict
	"""
	subd = defaultdict(list)
	if not isinstance(d, dict):
		raise Exception("Expected Dictionary; found {} ".format(str(type(d))))
	if not isinstance(d[0], list):
		raise Exception("Expected List in dictionary; found {} ".format(str(type(d[0]))))

	for k, v in d.items():
		for obj in v:
			if obj.name == region_name:
				subd[k].append(obj)
	return subd


def list_unique_regions_found_in_genomic_region_file(dico_of_regions):
	"""
	list the NAMES of the GENOMIC REGIONS found in genomic region file in order to check if the names match the names from the expected list of names of region
	See variabel '_ACCEPTED_REGION_NAMES' from REGION class for the list
	:param dico_of_regions: dictionary of regions
	:type dico_of_regions: defaultdict
	:return: tuple of unique names of regions
	:rtype: tuple
	"""
	list_all_regions = []
	if len(list(dico_of_regions.keys())) == 0:
		raise Exception("Dico CANNOT be Empty; No Genomic Regions Provided; Check your input")
	if not isinstance(dico_of_regions[list(dico_of_regions.keys())[0]], list):
		raise TypeError("List expected in dico_of_regions[X]")
	for lreg in dico_of_regions.values():
		for r in lreg:
			if r.name in list_all_regions:
				continue
			list_all_regions.append(r.name)
	return tuple(list_all_regions)


def set_namedtuple_regions(d_rname, dico_of_region_x):
	"""
	create a named tuple to organized data
	:param d_rname: name of the tuple (e.g: telomere, centromere, p, q, etc.)
	:type d_rname: string
	:param dico_of_region_x: dictionary we need to name
	:type dico_of_region_x: defaultdict
	:return: a namedtuple RNAME for RegionNAME and the dictionary associated it
	:rtype: namedtuple
	"""
	from collections import namedtuple
	RNAME = namedtuple('RNAME', 'rname, dict')
	return RNAME(d_rname, dico_of_region_x)


def make_obj_list_of_subset_regions(dico_of_regions):
	"""
	make the object needed for HRD_SCORE_CONTIG
	:param dico_of_regions:
	:type dico_of_regions:
	:return:
	:rtype:
	"""
	lrgs = []
	for region_name in list_unique_regions_found_in_genomic_region_file(dico_of_regions):
		lrgs.append(set_namedtuple_regions(region_name, subset_dico_by_region(dico_of_regions, region_name)))
		logger.debug(lrgs)
	return lrgs

def get_dico_for_given_region(region_name, lrgs):
	"""
	retrieve the dictionary that has all the intervals and contigs associated with the given region_name
	:param region_name: full name of the region such as p, q, telomere or centromere
	:type region_name: string
	:param lrgs: list of namedtupled tuples which have the requested dictionary
	:type lrgs: list
	:return: dictionary with all the contig associated to the given region_name
	:rtype: defaultdict
	"""
	for tpl in lrgs:
		if tpl.rname == region_name:
			return tpl[1]
	logger.warning("We will not check if segments are overlapping Telomere for now\ntuple with region named << {} >> not found; check your region file if the expected name of the region is present "
	               "in the file;\nif you want to deal with {}, modify the genomic region file accordingly;".format(region_name, region_name.upper()))
	return None


def check_if_telomeres_exist_for_contig(dico_genomic_regions):
	"""
	Check if we have both telomeres defined in genomic regions file
	:param dico_genomic_regions:
	:type dico_genomic_regions:
	:return:
	:rtype:
	"""
	d_rnames = defaultdict(list)
	for kreg in dico_genomic_regions.keys():
		_c = 0
		for reg in dico_genomic_regions[kreg]:  # dico_genomic_regions[kreg] is a list of REGIONs
			logger.debug("{} has {}".format(kreg, reg.name))
			d_rnames[kreg].append(reg.name)
			if 'telomere' == reg.name:
				logger.debug("{} has telomere".format(kreg))
				_c += 1
		if _c == 0:
			return False
		if _c > 2:
			raise ValueError(" Wrong Number of Telomere for Contig {}; Expected 2 Telomeres, found: {}; Check your Input".format(kreg, str(_c)))
		logger.debug("OK Found {} Telomeres in contig {}".format(str(_c), kreg))
	return True  # if we wanted to implement the fact that we do not have telomere and make the p_q arm accordingly


def check_if_centromere_exist_for_contig(dico_genomic_regions):
	"""
	we make sure that the centromere exists in the genomic region file for EACH chromosome or contig found in segment file
	:param dico_genomic_regions: dictionary of genomic regions in the form dico[contig] for which the value is a list of regions
	:type dico_genomic_regions: dict
	:return: If only one centromere is missing, we abort the computation here
	:rtype: None
	"""
	l_rnames = []
	for kreg in dico_genomic_regions.keys():
		for reg in dico_genomic_regions[kreg]:  # dico_genomic_regions[kreg] is a list of REGIONs
			logger.debug("{} has {}".format(kreg, reg.name))
			l_rnames.append(reg.name)
		if 'centromere' not in l_rnames:
			raise ValueError("Centromere for Contig {} is NOT defined in Genomic Regions file; Check your Input; "
			                 "Each contig found in SEG file MUST have a centromere defined in the genomic region file".format(kreg))

def check_if_all_contigs_in_seg_file_are_represented_in_genomic_regions(dico_genomic_regions, dico_seg):
	"""
	Check if all the contigs define
	:param dico_genomic_regions:
	:type dico_genomic_regions:
	:param dico_seg:
	:type dico_seg:
	:return:
	:rtype:
	"""
	sub_d_gen_reg = defaultdict(list)
	if all(contig in dico_genomic_regions.keys() for contig in dico_seg.keys()):
		# subset genomic region with contigs in seg file
		for contig in natsorted(dico_seg.keys()):
			sub_d_gen_reg[contig] = dico_genomic_regions[contig]
			check_if_centromere_exist_for_contig(sub_d_gen_reg)
		logger.debug("All the Contigs defined in SEG file have a Centromere defined in Genomic Region file")
		return sub_d_gen_reg
	else:
		logger.error("Contigs found in SEG file are NOT found in GENOMIC REGION file; Check your input")
		return False


def get_data_to_make_p_arm(reg_centro, reg_telo):
	"""
	Re-adjusting a contig's p-arm genomic territory using centromeric region and One telomeric region
	:param reg_centro: list of one instance of a REGION from class REGION defining a centromeric region
	:type reg_centro:  REGION
	:param reg_telo: list of one instance of a REGION from class REGION defining a teloomeric region
	:type reg_telo: REGION
	:return: an instance of class REGION for q arm
	:rtype: REGION
	"""
	# region example: chrXXYY	58467900	62522800	centromere	4054900
	# reg_centro is a list of one region
	# reg_telo is a list of two regions
	logger.debug("Captured Region Centromere : {}".format(reg_centro[0]))
	logger.debug("and captured region Telomere   : {}".format(reg_telo))
	end_pos = reg_centro[0].start
	start_pos = min(reg_telo[0].start, reg_telo[1].start)
	if end_pos == 0 or end_pos < start_pos:
		# this means we deal with acrocentric chromosome, i.e. no p-arm defined
		start_pos = end_pos
	return REGION(reg_centro[0].contig, start_pos, end_pos, "p")


def get_data_to_make_q_arm(reg_centro, reg_telo):
	"""
	Re-adjusting a contig's q-arm genomic territory using centromeric region and One telomeric region
	:param reg_centro: list of one instance of a REGION from class REGION defining a centromeric region
	:type reg_centro:  list
	:param reg_telo: list of one instance of a REGION from class REGION defining a teloomeric region
	:type reg_telo: list
	:return: an instance of class REGION for q arm
	:rtype: REGION
	"""
	start_pos = reg_centro[0].end
	end_pos = max(reg_telo[0].start, reg_telo[1].start)
	return REGION(reg_centro[0].contig, start_pos, end_pos, "q")


def get_data_to_make_p_arm_chrom(reg_centro, reg_telo, reg_chrom):
	"""
	# region example: chrXXYY	58467900	62522800	centromere	4054900
	# reg_centro is a list of one region
	# reg_telo is a list of two regions
	# reg_chrom is a list of one REGION object
	:param reg_centro: REGION representing a centromere
	:type reg_centro: REGION
	:param reg_telo: REGION representing a telomere
	:type reg_telo: REGION
	:param reg_chrom: REGION representing a chromosome
	:type reg_chrom: REGION
	:return: recreated REGION using values from centromere, telomere and/or chromosome
	:rtype: REGION
	"""
	logger.debug("Captured Region Centromere : {}".format(reg_centro[0]))
	logger.debug("and captured region Telomere   : {}".format(reg_telo))
	start_pos = min(reg_telo[0].start, reg_telo[1].start, reg_chrom[0].start)
	end_pos = reg_centro[0].start
	return REGION(reg_centro[0].contig, start_pos, end_pos, "p")


def get_data_to_make_q_arm_chrom(reg_centro, reg_telo, reg_chrom):
	"""
	# same comment as in functions get_data_to_make_p_arm_chrom but for chromosome's q-arm
	:param reg_centro: REGION representing a centromere
	:type reg_centro: REGION
	:param reg_telo: REGION representing a telomere
	:type reg_telo: REGION
	:param reg_chrom: REGION representing a chromosome
	:type reg_chrom: REGION
	:return: recreated REGION using values from centromere, telomere and/or chromosome
	:rtype: REGION
	"""
	logger.info(reg_centro)
	start_pos = reg_centro[0].end
	end_pos = max(reg_telo[0].start, reg_telo[1].start, reg_chrom[0].end)
	return REGION(reg_centro[0].contig, start_pos, end_pos, "q")


def get_data_to_make_p_arm_from_centromere_only(reg_centro, tpl_min_max):
	"""
	# region example: chrX	58467900	62522800	centromere	4054900
	# reg_centro is a list of one region
	# reg_telo is a list of two regions
	:param reg_centro: REGION representing a centromere
	:type reg_centro: REGION
	:param tpl_min_max: tuple representing the min and the max of the segments found in a contig from the SEG file
	:type tpl_min_max: tuple
	:return: recreated REGION using values from centromere and min and max segments in a contig
	:rtype: REGION
	"""
	logger.debug("Captured Region Centromere : {}".format(reg_centro[0]))
	end_pos = reg_centro[0].start
	start_pos = tpl_min_max[0]
	if end_pos == 0 or end_pos < start_pos:
		# this means we deal with acrocentric chromosome, i.e. no p-arm defined
		# but if defined it should lead to a size of zero when acrocentriq
		start_pos = end_pos
	return REGION(reg_centro[0].contig, start_pos, end_pos, "p")


def get_data_to_make_q_arm_from_centromere_only(reg_centro, tpl_min_max):
	"""
		# same comment as in functions get_data_to_make_p_arm_from_centromere_only but for chromosome's q-arm
	:param reg_centro: REGION representing a centromere
	:type reg_centro: REGION
	:param tpl_min_max: tuple representing the min and the max of the segments found in a contig from the SEG file
	:type tpl_min_max: tuple
	:return: recreated REGION using values from centromere and min and max segments in a contig
	:rtype: REGION
	"""

	start_pos = reg_centro[0].end
	end_pos = tpl_min_max[1]
	return REGION(reg_centro[0].contig, start_pos, end_pos, "q")


def make_pq_arms_from_centromeres(lrgsnt, telomeres_exist=False, chromosomes_exist=False, dico_min_max_segs_per_contig={}):
	"""
	1. Before running this function we had to make sure that each contig has a centromere region defined
	2. We can have Telomeres defined or not
	3. if we do not have telomeres defined, we make p and q arm like this:
	P-arm == 0 or from min_chromosome region up to min_centromere region
	Q-arm == max_centromere region to max contig length
	4. if we do have telomeres defined, we make p and q arm like this:
	P-arm == max telomere region from 1st telomere up to min_centromere region
	Q-arm == max_centromere region up to min of 2nd telomere

	:param lrgsnt: list of namedtupled region where tuple is composed of (region_name, REGION)
	:type lrgsnt: list
	:param telomeres_exist: result from the test if we have telomere or not defined
	:type telomeres_exist: boolean
	:return: dictionary of contig where items are dictionaries of arm (key) and REGION (value)
	:rtype: list
	"""
	d_p = defaultdict(list)
	d_q = defaultdict(list)
	if telomeres_exist and chromosomes_exist:
		for tpl in lrgsnt:
			if tpl.rname == "centromere":
				centromere_dico = tpl.dict
			if tpl.rname == "telomere":
				telomere_dico = tpl.dict
			if tpl.rname == "chromosome" or tpl.rname == "contig":
				chromosome_dico = tpl.dict
		for k in centromere_dico:
			d_p[k].append(get_data_to_make_p_arm_chrom(centromere_dico[k], telomere_dico[k], chromosome_dico[k]))
			d_q[k].append(get_data_to_make_q_arm_chrom(centromere_dico[k], telomere_dico[k], chromosome_dico[k]))
		lrgsnt.append(set_namedtuple_regions("p", d_p))
		lrgsnt.append(set_namedtuple_regions("q", d_q))
	elif telomeres_exist and not chromosomes_exist:
		centromere_dico = get_dico_for_given_region("centromere", lrgsnt)
		telomere_dico = get_dico_for_given_region("telomere", lrgsnt)
		for k in centromere_dico:
			d_p[k].append(get_data_to_make_p_arm(centromere_dico[k], telomere_dico[k]) )
			d_q[k].append(get_data_to_make_q_arm(centromere_dico[k], telomere_dico[k]) )
		lrgsnt.append(set_namedtuple_regions("p", d_p))
		lrgsnt.append(set_namedtuple_regions("q", d_q))
	elif not telomeres_exist and not chromosomes_exist:
		centromere_dico = get_dico_for_given_region("centromere", lrgsnt)
		for k in centromere_dico:
			if k not in dico_min_max_segs_per_contig.keys():
				continue
			d_p[k].append(get_data_to_make_p_arm_from_centromere_only(centromere_dico[k], dico_min_max_segs_per_contig[k]))
			d_q[k].append(get_data_to_make_q_arm_from_centromere_only(centromere_dico[k], dico_min_max_segs_per_contig[k]))
		# by appending data to the lrgsnt we update the current object and need to return it
		lrgsnt.append(set_namedtuple_regions("p", d_p))
		lrgsnt.append(set_namedtuple_regions("q", d_q))
	logger.debug("updated object 'lrgsnt': ")
	logger.debug(lrgsnt)
	return lrgsnt
