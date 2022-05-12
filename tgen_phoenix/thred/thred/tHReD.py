#!/usr/bin/env python3

import logging
import os
from collections import defaultdict
from logging.config import fileConfig
from os import path
from sys import exit
from time import time

from natsort import natsorted

from modules import functions as funcs
from modules import inputs
from modules import outputs
from modules import REGIONS
from modules import SEGMENTS
from scripts import cut_centromere_out

dir_root = path.dirname(__file__)
logging.config.fileConfig(path.join(dir_root, 'logging_config.ini'), disable_existing_loggers=False)
logger = logging.getLogger(__name__)


def main():
	parser = inputs.make_parser_args()
	logger.info("parsing args ...")
	args = vars(parser.parse_args())
	logger.debug(str(args))
	tpl_args = inputs.capture_and_return_args(args)
	inputs.print_recap_inputs(tpl_args)

	# Assign values of tpl_args to variables
	make_plots = tpl_args[13]
	karyo_file = tpl_args[12]
	list_exclude_contigs = [str(x).strip() for x in tpl_args[11].split(",")]
	any_id = tpl_args[10]
	samplename = tpl_args[9]
	outfile = tpl_args[8]
	th_pct_overlapping = tpl_args[7]
	list_contigs = tpl_args[6]
	minsize_segment = tpl_args[5]
	th_log2r_deletion = tpl_args[4]
	# ncpus = tpl_args[3]  # not used yet
	dirout = tpl_args[2]
	gen_reg_file = tpl_args[1]
	segfile = tpl_args[0]

	# ------------------
	# CUR DIR @ START
	# ------------------
	logger.info("curdir = {}".format(os.path.abspath(os.curdir)))
	# ---------------------------
	# read Genomic Regions File
	# --------------------------
	d_genomic_regions = REGIONS.read_file_genomic_regions(gen_reg_file)
	list_all_genomic_named_regions = REGIONS.list_unique_regions_found_in_genomic_region_file(d_genomic_regions)

	for region_name in ['centromere']:  # Mandatory centromere must be defined in the genomic region file ; We can add more mandatory later if needed
		if region_name not in list_all_genomic_named_regions:
			raise NameError("GENOMIC REGION MISSING: {}. That region MUST be defined in the genomic regions file provided".format(region_name))

	# list of namedtuples of dico_of_keys==contigs with_list_of_regions as values
	lrgs = REGIONS.make_obj_list_of_subset_regions(d_genomic_regions)

	# ---------------------------------------------------
	# read Segments File aka GATK's CNV generated file
	# ---------------------------------------------------

	dsegs = SEGMENTS.read_and_sort_segment_file(segfile, "\t", list_exclude_contigs=list_exclude_contigs, list_keep_contigs=list_contigs)
	list_contigs = natsorted(dsegs.keys())
	logger.info("list_contigs that will be processed ==> {}".format(list_contigs))

	# ---------------------------------------------------
	# check if the SEG file contains only one sample
	# ---------------------------------------------------
	unique_sample = SEGMENTS.check_count_samples_in_seg_file(dsegs)
	if samplename == ".":
		samplename = unique_sample
		
	# ------------------------------------------------------------------------
	# check if centromeres are defined for all the contigs found in seg file
	# ------------------------------------------------------------------------
	d_genomic_regions = REGIONS.check_if_all_contigs_in_seg_file_are_represented_in_genomic_regions(d_genomic_regions, dsegs)

	# ----------------------------------------------------------------------
	# Check if segments overlap centromere, and split segments accordingly
	# ----------------------------------------------------------------------
	dsegs = cut_centromere_out.cut_centromere_out_of_segment(dsegs, REGIONS.get_dico_for_given_region("centromere", lrgs))
	
	# ----------------------------------------------------------------------
	# capture min and max from each contig for making future arm territory
	# ----------------------------------------------------------------------
	dico_min_max_segs_per_contig = SEGMENTS.get_min_max_segments_per_contig(dsegs)
	logger.info("dico_min_max_segs_per_contig : {}".format(str(dico_min_max_segs_per_contig)))

	# ----------------------------------------------------
	# Check if the genomic region file includes the telomere definition
	# ---------------------------------------------------
	are_telomeres_defined = REGIONS.check_if_telomeres_exist_for_contig(d_genomic_regions)
	lrgs = REGIONS.make_pq_arms_from_centromeres(lrgs, telomeres_exist=are_telomeres_defined, chromosomes_exist=False, dico_min_max_segs_per_contig=dico_min_max_segs_per_contig)

	# -------------------------------
	# loops over arms and contigs
	# -------------------------------
	# init variables
	# ------------------------------
	ldics = list()
	dico_excluded_90_all_contigs = defaultdict(dict)
	dico_new_terr = defaultdict(dict)
	dico_segs_by_contig_and_by_arm = defaultdict(dict)
	ct = 0
	size_territory_all_contigs = 0
	d_len_seg_filtered = dict()
	d_filt_segs = dict()
	sum_size_all_segments_all_contigs = 0
	dico_telomeres = REGIONS.get_dico_for_given_region("telomere", lrgs)
	dico_centromeres = REGIONS.get_dico_for_given_region("centromere", lrgs)
	# ------------------------------
	for contig in natsorted(list_contigs):
		d_filt_segs[contig] = dict()
		d_len_seg_filtered[contig] = dict()
		dico_excluded_90_all_contigs[contig] = dict()

		for arm in ["p", "q"]:
			dico_arm = REGIONS.get_dico_for_given_region(arm, lrgs)
			# ---------------------------------------------------------------------
			# RECAPTURE GENOME TERRITORY AND SORT OUT SEGMENTS by Chromosome ARM
			# ---------------------------------------------------------------------
			newdic_segs_in_arm, new_arm_territory, c, size_territory = funcs.make_new_arm_territory_by_contig(dsegs, dico_centromeres, dico_arm, arm, contig)
			ldics.append(newdic_segs_in_arm)
			dico_segs_by_contig_and_by_arm[contig][arm] = newdic_segs_in_arm
			dico_new_terr[contig][arm] = new_arm_territory
			ct += c
			# prints >-)
			logger.debug("{}_{} ; count: {}, size_territory: {}".format(contig, arm, str(c), str(size_territory)))
			funcs.print_logger_debug(
					"contig + arm",
					contig + "_" + arm,
					"newdic_segs_in_arm",
					newdic_segs_in_arm,
					"dico_segs_by_contig_and_by_arm",
					dico_segs_by_contig_and_by_arm,
					"new_arm_territory dico:",
					new_arm_territory)

			size_territory_all_contigs += size_territory
			# ---------------------------------------------------
			# FILTERING of SEGMENTS to COUNT HRD
			# ---------------------------------------------------
			try:
				# -------------------------
				# INIT the Dictionaries
				# -------------------------
				d_filt_segs[contig][arm] = []
				d_len_seg_filtered[contig][arm] = []
				dico_excluded_90_all_contigs[contig][arm] = []

				d_filt_segs_per_contig, d_len_seg_filtered_per_contig, dico_sum_size_segments_loh_per_contig_per_arm, dico_excluded_90 = \
					funcs.filter_out_segments(dico_segs_by_contig_and_by_arm, dico_new_terr, contig, arm, dico_telomeres=dico_telomeres, th_loh_min=th_log2r_deletion,
					                          th_min_seg_length=minsize_segment, th_pct_overlapping_arm=th_pct_overlapping, th_pct_overlapping_telomere=th_pct_overlapping)

				funcs.print_logger_debug(
					"d_filt_segs_per_contig",
					d_filt_segs_per_contig,
					"d_len_seg_filtered_per_contig",
					d_len_seg_filtered_per_contig,
					"dico_excluded_90",
					dico_excluded_90)

				d_filt_segs[contig][arm] = d_filt_segs_per_contig[contig][arm]
				d_len_seg_filtered[contig][arm] = d_len_seg_filtered_per_contig[contig][arm]
				dico_excluded_90_all_contigs[contig][arm] = dico_excluded_90[contig][arm]
				logger.debug(dico_excluded_90_all_contigs)
				sum_size_all_segments_all_contigs += sum(
						int(dico_sum_size_segments_loh_per_contig_per_arm[contig][arm]) for arm in dico_sum_size_segments_loh_per_contig_per_arm[contig])

			except ValueError as VE:
				logger.error("Error in Try filter_out_segments: " + str(VE))
			except TypeError as TE:
				logger.error("Error in Try filter_out_segments: " + str(TE))
			except Exception as E:
				logger.error("Unkown Error while filtering Segments" + str(E))

			logger.debug("current sum_size_all_segments_all_contigs {}".format(str(sum_size_all_segments_all_contigs)))

	# ------------------------------------------------------------------------
	# COUNT BREAKS
	# ------------------------------------------------------------------------
	BK = funcs.count_breaks(dico_segs_by_contig_and_by_arm, th_loh_min=th_log2r_deletion, th_min_diff_copy_number=0.5, th_min_seg_length=1 * 10 ** 6, th_max_dist_btw_segs=3 * 10 ** 6)

	# ---------------------
	# OUTPUT RESULTS
	# ---------------------
	tpl_res_values = (sum_size_all_segments_all_contigs, size_territory_all_contigs, BK)
	hrd_score = funcs.calculate_hrd_score(tpl_res_values, withBK=True)
	#hrd_score = sum_size_all_segments_all_contigs / size_territory_all_contigs

	tpl_res_values = tpl_res_values + hrd_score
	logger.debug("HRD_SCORE = {}".format(str(hrd_score)))
	logger.info("HRD_SCORE_tpl_res_values = {}".format(str(tpl_res_values)))
	logger.debug("BK == {}  _____  and BKPB == {} , and newHRD score == {} ".format(str(BK), str(hrd_score[0]), str(hrd_score[0] + hrd_score[1])))

	# if samplename == "." and len()
	outputs.write_output_scores(tpl_res_values, dirout=dirout, outfile=outfile, sample=samplename, id=any_id)

	# write out segments kept
	outputs.write_output_filtered_segments(d_filt_segs, dirout=dirout, outfile=outfile)
	outputs.write_output_original_segments(dsegs, dirout=dirout, outfile=outfile)

	# write out segment excluded because they overlapped more than XX% the p or q arm
	outputs.write_output_filtered_segments(dico_excluded_90_all_contigs, dirout=dirout, outfile=outfile+"_excluded"+str(int(100*th_pct_overlapping)), default_outfilename="hrd_excluded_segments.txt")

	# write the new genome territory calculated out of min//max segments and min//max centromere
	outputs.write_new_territory(dico_new_terr, dirout=dirout, outfile=outfile)

	# ------------------------------------------
	# MAKE simple Karyotype plots (optional)
	# ------------------------------------------
	if make_plots:
		# Using original segments in SEG file
		plot_filename = samplename + "_original_segments" if samplename != "." else "segments_original"
		outputs.make_karyoplot(karyo_file, plot_filename, th_log2r_deletion, dico_segs_by_contig_and_by_arm, part=1, dirout=dirout)
		outputs.make_karyoplot(karyo_file, plot_filename, th_log2r_deletion, dico_segs_by_contig_and_by_arm, part=2, dirout=dirout)

		# plots after filtering segments
		plot_filename = samplename + "_segments_filtered" if samplename != "." else "karyoplot_flt"
		outputs.make_karyoplot(karyo_file, plot_filename, th_log2r_deletion, d_filt_segs, part=1, dirout=dirout)
		outputs.make_karyoplot(karyo_file, plot_filename, th_log2r_deletion, d_filt_segs, part=2, dirout=dirout)

		# plots with excluded segments
		plot_filename = samplename + "_segments_excluded" if samplename!="." else "karyoplot_excl"
		outputs.make_karyoplot(karyo_file, plot_filename, th_log2r_deletion, dico_excluded_90_all_contigs, part=1, dirout=dirout)
		outputs.make_karyoplot(karyo_file, plot_filename, th_log2r_deletion, dico_excluded_90_all_contigs, part=2, dirout=dirout)


if __name__ == '__main__':
	try:
		ts = time()
		main()
		logger.info("total runtime is: {} secs".format(time() - ts))
	except Exception as e:
		logger.error(e, exc_info=True)
		logger.info("total runtime is: {} secs".format(time() - ts))
		exit(2)
	exit(0)
