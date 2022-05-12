import logging.config
from collections import defaultdict
from operator import attrgetter
from os import path
from math import log10, log2, log
from modules import REGIONS
from modules import SEGMENTS
# from line_profiler_pycharm import profile

# Profiler cam from:
# https://plugins.jetbrains.com/plugin/16536-line-profiler

dir_root = path.dirname(path.dirname(__file__))
logging.config.fileConfig(path.join(dir_root, 'logging_config.ini'), disable_existing_loggers=False)
logger = logging.getLogger(__name__)


def make_new_arm_territory_by_contig(dseg, d_centro, d_arm, arm_name, contig):
	"""
	Capture the Genome Territory as defined by the design.
	from first_segment in encountered in p-arm to min_centromere, and from max_centromere to max_segment in q-arm
	:param dseg: dictionary of segments to capture min and max coordinates on each contig
	:type dseg: defaultdict
	:param d_centro: dictionary of the genomic region corresponding to what is called centromere in chromosomes
	:type d_centro: dict
	:param d_arm: dictionary of the genomic region corresponding to what is called arms in chromosomes
	:type d_arm: dict
	:param arm_name: arm name of the contig
	:type arm_name: basestring
	:param contig: contig name
	:type contig: basestring
	:return: (list_of_segments, region_of_an_arm, count_of_processed_segments_in_list_of_segments)
	:rtype:  (list_of_SEGMENTS, REGIONS, integer, integer)
	"""
	
	dseg_arm = defaultdict(dict)
	dseg_arm[contig][arm_name] = list()
	# dnt_arm = defaultdict(dict)
	# dnt_arm[contig][arm_name] = list()
	c = 0
	lsegs = list()
	size_new_territory = 0
	
	for seg in dseg[contig]:
		res = SEGMENTS.are_segments_overlapping(seg, d_arm[contig][0])
		if res[0]:
			c += 1  # counting overlapping segments
			dseg_arm[contig][arm_name].append(seg)
			# list segments in_current_contig_and_current_arm
			lsegs.append(seg)
	# we need to check if data exist for the two following objects otherwise ERROR will be raised when calculating min and max below
	if not d_centro[contig] or not d_arm[contig][0]:
		raise ValueError("ERROR: missing data in either d_centro[contig] == {}  or d_arm[contig][0] == {}\
		for contig <{}>; Aborting; Check your genomic regions definition inputs".format(d_centro[contig], d_arm[contig],
		                                                                                contig))
	
	# We now have all the segments belonging to an arm;  ## WE ALSO CAN SPLIT THAT INTO ANOTHER FUNCTION, even though it will be redundant to parse the dico segments twice
	# let's make the new territory of the arm based on both centromere and the segments found in that arm
	if arm_name == "p":
		if lsegs:
			min_num = min(lsegs, key=attrgetter('start')).start
			max_num = min(d_centro[contig], key=attrgetter('start')).start
			if max_num <= min_num:
				min_num = 0  # HARDCODED
		else:
			max_num = min(d_centro[contig], key=attrgetter('start')).start
			min_num = max_num  # HARDCODED
	else:
		min_num = min(d_centro[contig], key=attrgetter('end')).end
		if lsegs:
			max_num = max(lsegs, key=attrgetter('end')).end
			if max_num <= min_num:
				max_num = min_num  # HARDCODED
		else:
			max_num = min_num  # HARDCODED
	size_new_territory += (max_num - min_num) + 0  # BASE HARDCODED
	# dnt_arm[contig][arm_name].append(REGION.REGION(contig=contig, start=min_num, end=max_num, name=arm_name))
	print_logger_debug(lsegs, c, size_new_territory)
	return lsegs, REGIONS.REGION(contig=contig, start=min_num, end=max_num, name=arm_name), c, size_new_territory


def is_overlapping_telomeres(telomeres, segment, th_pct_overlapping=0.90):
	"""
	check if the segment overlap a telomere
	:param telomeres: list of REGION or SEGMENT
	:type telomeres: list
	:param segment: a SEGMENT
	:type segment: SEGMENT
	:param th_pct_overlapping: if the % of the segment overlaps the telomere by more of X % of its size, we reject it.
	:type th_pct_overlapping: float
	:return: True or False
	:rtype: Boolean
	"""
	for telomere in telomeres:
		res = SEGMENTS.are_segments_overlapping(segment, telomere)
		if res[0] and res[2] >= th_pct_overlapping:
			logger.info("OVERLAP TELOMERE by 90%" * 50)
			return True
	return False


def count_breaks(dico_segs_by_contig_and_by_arm, th_loh_min=-0.1613, th_min_diff_copy_number=0.5, th_min_seg_length=1 * 10 ** 6, th_max_dist_btw_segs=3 * 10 ** 6):
	"""
	Count the number of breaks between segments if adjacent segments are at least of th_min_seg_length and
	 if their distance is less than th_min_dist_btw_seg and
	 if the log2 ratio between the two adjacent segments is more than the th_min_diff_copy_number


	:param th_max_dist_btw_segs: MAXimal distance between two segments to consider a break
	:type th_max_dist_btw_segs: int
	:param th_min_diff_copy_number: Expected min difference in Copy Number domain to count it as a break between two segments; default half a copy, i.e. 0.5
	:type th_min_diff_copy_number: float
	:param dico_segs_by_contig_and_by_arm: dico with all the segments from each arm and contig
	:type dico_segs_by_contig_and_by_arm: defaultdict
	:param th_loh_min: threshold for the log2r we start considering segments as a deletion or a gain
	:type th_loh_min: float
	:param th_min_seg_length: minimum size of a segment to be considered valid to process; default 1 Million bp
	:type th_min_seg_length: int
	:return: Count of Breaks according to the rules
	:rtype: integer
	"""
	cbks = 0  # count Breaks
	for contig, dico_arms in  dico_segs_by_contig_and_by_arm.items():
		for arm, lsegs in dico_arms.items():
			logger.debug("{} -- {}__{} -- {}".format("&" * 35, contig, arm, "&" * 45))
			
			prev_seg = None
			for segment in lsegs:
				if prev_seg is not None and segment.contig == prev_seg.contig:
					logger.debug(
						"{} || {} || {} || {}".format(segment.mean <= abs(th_loh_min) * -1, prev_seg.mean <= abs(th_loh_min) * -1, segment.mean >= abs(th_loh_min), prev_seg.mean >= abs(th_loh_min)))
					if segment.get_distance_btw_segments(prev_seg) <= th_max_dist_btw_segs \
							and segment.slen >= th_min_seg_length and prev_seg.slen >= th_min_seg_length \
							and (segment.mean <= abs(th_loh_min)*-1 or prev_seg.mean <= abs(th_loh_min)*-1 or segment.mean >= abs(th_loh_min) or prev_seg.mean >= abs(th_loh_min)) \
							and segment.get_copy_number_diff_btw_segments(prev_seg, th_min_diff_copy_number=th_min_diff_copy_number):
						logger.debug("#### {}  ######\nDIFF_OK:  {} ----- VS---- {}".format(contig, str(prev_seg), str(segment)))
						cbks += 1
						prev_seg = segment
					else:
						prev_seg = None
				else:
					prev_seg = segment
	# if contig == "chr1": break
	return cbks


# @profile
def filter_out_segments(dico_segs_by_contig_and_by_arm,
                        dico_new_terr,
                        contig, arm,
                        dico_telomeres=None,
                        th_loh_min=-0.1613,
                        th_min_seg_length=10 ** 6,
                        th_pct_overlapping_arm=0.90,
                        th_pct_overlapping_telomere=0.90):
	"""
	filter_out_segments_by different features:
	1) filter by DEL thresholds (th_loh_min)
	2) filter if segments is overlapping_telomere by at least th_pct_overlapping
	3) filter if the sum of the segments in arm is X% of the arm (th_pct_overlapping)
	:param dico_segs_by_contig_and_by_arm:
	:type dico_segs_by_contig_and_by_arm: defaultdict
	:param dico_new_terr:
	:type dico_new_terr: defaultdict
	:param dico_telomeres: if exist, dictionary of telomere regions is provided here
	:type dico_telomeres: defaultdict
	:param contig: name of the contig to process
	:type contig: string
	:param arm: name of the arm
	:type arm: string
	:param th_loh_min: threshold for the log2r we start considering segments as a deletion
	:type th_loh_min: float
	:param th_min_seg_length: minimum size of a segment to be considered valid to process
	:type th_min_seg_length: int
	:param th_pct_overlapping_arm: threshold for excluding segments if their sum is over th_pct_value represented as a ratio (0.90 == 90%)
	:type th_pct_overlapping_arm: float
	:param th_pct_overlapping_telomere:  threshold for excluding segments if 90% of the segment overlaps the telomere (0.90 == 90%) ; NOT IMPLEMENTED YET
	:type th_pct_overlapping_telomere: float
	:return: dico of filtered SEGMENTS, dico of LENGTH of filtered SEGMENTS, SUM_LENSGTH_of_filtered_SEGMENTS, DICO_of_SEGMENTS_excluded_b/c_Overlapping_more_than_pct_the_arm
	:rtype: dict, dict, int, dict
	"""
	# INIT VARIABLES
	dico_excluded_90 = defaultdict(dict)
	dsegs_flt = defaultdict(dict)
	dsegs_flt_len = defaultdict(dict)
	dico_sum_size_segments_loh_per_arm = defaultdict(dict)
	dico_sum_size_segments_loh_per_arm[contig][arm] = 0
	dsegs_flt[contig][arm] = list()
	dsegs_flt_len[contig][arm] = list()
	dico_excluded_90[contig][arm] = list()
	logger.debug("filtering segments on {}_{}".format(contig, arm))
	# ------------------------------------------------------------
	for segment in dico_segs_by_contig_and_by_arm[contig][arm]:
		print_logger_debug("contig_-_arm", contig + " _-_ " + arm, "dico_new_terr[contig][arm]", dico_new_terr[contig][arm], "segment", segment)
		## CHECK IF OVERLAPPING p or q arm
		res = SEGMENTS.are_segments_overlapping(segment, dico_new_terr[contig][arm])  # res[0] represents True or False values
		## WE TEST HERE IF WE DEAL WITH LOH (aka DELETION) and if MIN LENGTH SEGMENT OK
		if res[0] and segment.mean <= th_loh_min and segment.slen > th_min_seg_length:
			logger.debug("Segment << {} >>  on {} PASSED for TESTs 1) overlap arm + 2) less than {} and 3) segment_length > {} ".format(str(segment), contig+"_"+arm, th_loh_min, th_min_seg_length))
			if dico_telomeres is None:
				dsegs_flt[contig][arm].append(segment)
				dsegs_flt_len[contig][arm].append(segment.slen)
				dico_sum_size_segments_loh_per_arm[contig][arm] += segment.slen
			elif not is_overlapping_telomeres(dico_telomeres[contig], segment, th_pct_overlapping=th_pct_overlapping_telomere):
				logger.debug("YEAH-YEAH  --> IT does NOT overlap telomere .. adding segment" + str(segment))
				dsegs_flt[contig][arm].append(segment)
				dsegs_flt_len[contig][arm].append(segment.slen)
				dico_sum_size_segments_loh_per_arm[contig][arm] += segment.slen
		# Is the SUM of the segments filtered for current contig and arm overlaps by th_pct_overlapping
		len_current_arm = dico_new_terr[contig][arm].slen
		logger.debug("contig + arm  == {} -- {} ; length_current_arm == {} ... and ... sum_size_segments_loh_per_arm == {} ".format(str(contig), str(arm), str(len_current_arm), str(dico_sum_size_segments_loh_per_arm[contig][arm]) ))
		if len_current_arm < dico_sum_size_segments_loh_per_arm[contig][arm]:
			raise Exception("ERROR: WHY is that? len_current_arm < sum_size_segments_loh_per_arm? How Come?")
		if len_current_arm != 0 and dico_sum_size_segments_loh_per_arm[contig][arm] / len_current_arm >= th_pct_overlapping_arm:
			dico_excluded_90[contig][arm] = list()
			# we check if the sum of the segments overlaps more then `th_pct_overlapping_arm` in %
			logger.debug("sum_size_segments_loh_per_arm {} / {} len_current_arm >= {} th_pct_overlapping".format(
				str(dico_sum_size_segments_loh_per_arm[contig][arm]), str(len_current_arm), str(th_pct_overlapping_arm)))
			logger.warning("EXCLUSION: Segments on {}_{} overlap more than {}% the arm".format(str(contig), str(arm), str(th_pct_overlapping_arm*100)))
			dico_excluded_90[contig][arm] = dsegs_flt[contig][arm]
			logger.debug("WHY: {}_{}  WHY: {}".format(contig, arm, '\n'.join([ str(x) for x in dico_excluded_90[contig][arm] ] )))
			dsegs_flt[contig][arm] = []
			dsegs_flt_len[contig][arm] = []
			dico_sum_size_segments_loh_per_arm[contig][arm] = 0
	
	print_logger_debug(dsegs_flt, dsegs_flt, dico_new_terr)
	return dsegs_flt, dsegs_flt_len, dico_sum_size_segments_loh_per_arm, dico_excluded_90


def calculate_hrd_score(tpl_res_values, withBK=True):
	"""
	calculate HRD sccore
	:param tpl_res_values: 3 values:
	:type tpl_res_values:
	:param withBK:
	:type withBK:
	:return:
	:rtype:
	"""
	BKperMb = 0
	try:
		HRDscore = tpl_res_values[0] / tpl_res_values[1]
		if withBK:
			# BKperMb = round((tpl_res_values[2] + 2 * log(tpl_res_values[2])) * 10 ** 6 / tpl_res_values[1], 3)
			BKperMb = round(( tpl_res_values[2] / tpl_res_values[1] ) * 10 ** 6, 6)
	except ZeroDivisionError as zde:
		logger.error("Division by Zero captured; How come the Genome Territory is Zero in space?; check your input ; setting HRDscore to Zero to keep going to next segment; " + zde)
		return 0, 0
	except ArithmeticError as ae:
		logger.error("Zero Breaks ; so log(0) raises an Arithmetic Error to calculate the Breaks per Million bases value; " + ae)
		return 0, HRDscore
	return BKperMb, HRDscore


def get_hrd_scores(dico_scores):
	s_by_p = 0
	s_by_q = 0
	s = 0
	for k in dico_scores.keys():
		if "_p" in k:
			s_by_p += dico_scores[k]
		elif "_q" in k:
			s_by_q += dico_scores[k]
		else:
			raise ValueError("ERROR arm NOT FOUND with Key {} in dico of scores".format(k))
		s += dico_scores[k]
	return s, s_by_p, s_by_q


def print_logger_debug(*args):
	for item in args:
		item_var_name = [k for k, v in locals().items() if v == item][0]
		logger.debug("{} :".format(item_var_name))
		logger.debug(str(item))