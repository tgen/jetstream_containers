
import logging
from collections import defaultdict
from natsort import natsorted
from modules import SEGMENTS
from os import path

dir_root = path.dirname(path.dirname(__file__))
logging.config.fileConfig(path.join(dir_root, 'logging_config.ini'), disable_existing_loggers=False)
logger = logging.getLogger(__name__)


def cut_centromere_out_of_segment(dico_segs, dico_centro):
	"""
	chop out the segments that overlap the centromeric part of the contig.
	:param dico_segs: dictionary of contig of list of segments
	:type dico_segs: dict
	:param dico_centro:  dictionary of contig of list of REGIONs, specifically the REGIONS are only CENTROMERIC ones
	:type dico_centro: dict
	:return: an new dictionary of contigs with chopped segments if any of them were overlapping centromeres
	:rtype: defaultdict
	"""
	l_seg_contigs = dico_segs.keys()
	l_centro_contigs = dico_centro.keys()
	logger.debug(l_seg_contigs)
	logger.debug(l_centro_contigs)
	if not all(contig in l_centro_contigs  for contig in l_seg_contigs):
		raise ValueError("NOT all the contigs present in SEG file can be found in list of Centromere")

	new_dico_segs_cut = defaultdict(list)
	for contig in natsorted(l_seg_contigs):
		new_dico_segs_cut[contig] = []
		for seg in dico_segs[contig]:
			reg_centromere = dico_centro[contig][0]
			ovp = SEGMENTS.are_segments_overlapping(seg, reg_centromere)
			# (True, 4242000, 5.106905580117807, 100.0) means 100% of the centromeric bases are overlapping the segment, so segment overlaps the whole centromere
			if not ovp[0]:
				new_dico_segs_cut[contig].append(seg)
			else:
				if reg_centromere.start < seg.start and reg_centromere.end < seg.end:
					logger.debug("case1")
					logger.debug("contig: {}  and segment: {}".format(contig, str(seg)))
					seg.start = reg_centromere.end
					seg.slen = seg.end - seg.start
					new_dico_segs_cut[contig].append(seg)
					logger.debug("contig: {}  and segment: {}".format(contig, str(seg)))
				elif seg.start < reg_centromere.start and reg_centromere.end < seg.end:
					logger.debug("case2")
					seg1 = seg.copy()
					logger.debug("contig: {}  and segment: {}".format(contig, str(seg1)))
					seg1.end = reg_centromere.start
					seg1.slen = seg1.end - seg1.start
					logger.debug("contig: {}  and segment1: {}".format(contig, str(seg1)))
					new_dico_segs_cut[contig].append(seg1)
					del seg1
					seg.start = reg_centromere.end
					seg.slen = seg.end - seg.start
					new_dico_segs_cut[contig].append(seg)
					logger.debug("contig: {}  and segment2: {}".format(contig, str(seg)))
				elif seg.start < reg_centromere.start and seg.end < reg_centromere.end:
					logger.debug("case3")
					logger.debug("contig: {}  and segment: {}".format(contig, str(seg)))
					seg.end = reg_centromere.start
					seg.slen = seg.end - seg.start
					new_dico_segs_cut[contig].append(seg)
					logger.debug("contig: {}  and segment: {}".format(contig, str(seg)))
				elif reg_centromere.start < seg.start < reg_centromere.end and reg_centromere.start < seg.end < reg_centromere.end:
					logger.warning("case4 : the segment is smaller than centromere is within the centromere: " + str(seg))
				else:
					raise Exception("UnKnown Case: <<centromere|segments>> overlapping unknown case with segment: {};\nCheck you input or report an issue;".format(str(seg)))
		logger.debug([str(seg) for seg in dico_segs[contig]])
		logger.debug([str(seg) for seg in new_dico_segs_cut[contig]])
	return new_dico_segs_cut
