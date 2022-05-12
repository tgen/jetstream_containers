#!/usr/bin/env python3

import logging
import clusterings
from time import time
from collections import defaultdict
from makeFastq import process_regions_parallel, make_parser_args, capture_and_return_args, make_dico_from_igs_regions
from functools import partial
from multiprocessing import Pool


def main():
	ts = time()

	# IgH   chr14   ; #IgK      chr2  ; #IgL	chr22
	# NSD2	chr4    ; #CCND3    chr6  ; #MYC	chr8  ; #MAFA	chr8
	# CCND1	chr11   ; #CCND2	chr12 ; #MAF	chr16 ; #MAFB	chr20
	# regions = ["chr4:1798273-1998273", "chr6:41632262-42332262", "chr8:124987758-129487754", "chr8:142918584-143925832",
	# "chr11:68732532-69685232", "chr12:3690834-4690834", "chr16:78096103-79866103", "chr20:39671358-40971360"]
	# igs_regions = ["chr2:88750483-90261139", "chr14:105533663-106881350", "chr22:21995603-23057822"]
	# run command line example:
	# python3 dle.py \
	# -a "${USER}/PycharmProjects/filterMantaSVs/alignment/bwa/MMRF_1157/MMRF_1157_2_BM_CD138pos_T1_KBS5U.bwa.cram" \
	# -g "${USER}/PycharmProjects/filterMantaSVs/reference_genome/GRCh38tgen_decoy_alts_hla.fa" \
	# -r "chr11:5000000-5100000" "chr11:10000000-10100000" \
	# -i "chr2:88750483-90261139"	"chr14:105533663-106881350"	"chr22:21995603-23057822" \
	# -d / ${USER} / Documents / TESTS_RUNS	-p 	TEST0_small__  -q 0 -t 2


	parser = make_parser_args()
	logger.info("parse args")
	args = vars(parser.parse_args())
	logger.info(str(args))
	logger.debug("recap all inputs parameters:\n>>>> " + "\n>>>> ".join("{}\t{}".format(k, v) for k, v in args.items()))
	tbam_file, ref_genome_file, var_min_mapq, dir_out_fq, prefix_out_fq, ncpus, regions, igs_regions, igs_contig_names, \
	min_tuples, left_position_clustering_only, min_reads_per_cluster, min_dist_left_right, clustering_distance = capture_and_return_args(args)

	logger.info("")
	logger.info("Inputs parameters:")
	logger.info(">>> alignment file = " + tbam_file)
	logger.info(">>> ref_genome_file = " + ref_genome_file)
	logger.info(">>> regions = " + str(regions))
	logger.info(">>> igs regions = " + str(igs_regions))
	logger.info(">>> var_min_mapq = " + str(var_min_mapq))
	logger.info(">>> dir_output = " + str(dir_out_fq))
	logger.info(">>> prefix_out_fq = " + str(prefix_out_fq))
	logger.info(">>> threads = " + str(ncpus))
	logger.info(">>> min_tuples = " + str(min_tuples))
	logger.info(">>> left_position_clustering_only = " + str(left_position_clustering_only) )
	logger.info(">>> min_reads_per_cluster = " + str(min_reads_per_cluster))
	logger.info(">>> min_dist_left_right = " + str(min_dist_left_right))
	logger.info(">>> clustering_distance = " + str(clustering_distance))
	logger.info("")

	# ------------------------
	# init clustering object
	# ------------------------
	logger.info("make instance of clustering ...")
	oclust = clusterings.ReadClusteringManagement(None, None, selected_clustering_type="hclust", cluster_left_pos_only=left_position_clustering_only,
		min_tpls=min_tuples, min_dist_leftmost_to_rightmost_reads=min_dist_left_right, min_read_count_per_cluster_to_Pass=min_reads_per_cluster,
		threshold_distance_for_clustering=clustering_distance)

	# make dictionary from igs_regions
	dico_igs_contig_regions = make_dico_from_igs_regions(igs_regions)

	# -------------------------
	# pool regions
	# -------------------------
	logger.info("loop over regions by calling parallel processing ...")
	pool = Pool(ncpus)
	func = partial(process_regions_parallel, var_min_mapq,  tbam_file, ref_genome_file, dir_out_fq, prefix_out_fq, igs_contig_names, dico_igs_contig_regions, oclust)
	pool.map(func, iterable=regions)
	pool.close()
	pool.join()

	logger.info("total runtime is: {} secs".format(time()-ts))


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)s %(asctime)-15s %(module)s %(lineno)d\t %(message)s')
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)

if __name__ == '__main__':
	try:
		main()
	except Exception as e:
		logger.error(e, exc_info=True)
		exit(2)
	exit(0)
