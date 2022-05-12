#!/usr/bin/env python3

import argparse
import gzip
import logging
import os
import pysam
from Bio import Seq
from collections import defaultdict
from sys import exit
from functools import wraps

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
#logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(levelname)s %(asctime)-15s %(module)s %(lineno)d\t %(message)s')
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)


class UniqueStore(argparse.Action):
	"""
	class grabbed from stackOverflow (2018-11-08)
	https://stackoverflow.com/questions/23032514/argparse-disable-same-argument-occurences
	Thanks To the Community
	We override the function __call__ from argparse to check if one option is given more than once
	# WARNING -- WARNING #
	# Does not work with the option_arguments ; if we pass an optioanl args to the command line, \
	# it is considered as a second copy of the default optional arg; --> So does not work for optional args but we kept it for required \
	# to make sure that the user does not pass twice teh same arg
	# WARNING -- WARNING #
	"""

	def __call__(self, parser, namespace, values, option_string):
		if getattr(namespace, self.dest, self.default) is not None:
			parser.error(option_string + " appears several times.  Please modify your options.")
		setattr(namespace, self.dest, values)


def check_file(f):
	if f is None or not os.path.isfile(f) or not os.path.exists(f):
		logger.error("FNF: " + str(f))
		raise FileNotFoundError("FNF: " + str(f))


def make_parser_args():
	parser = argparse.ArgumentParser(add_help=True, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser._action_groups.pop()
	required = parser.add_argument_group('###-- required arguments --###')
	optional = parser.add_argument_group('###-- optional arguments --###')
	required_mutually_exclusive = required.add_mutually_exclusive_group(required=True)
	igs_mutually_exclusive = optional.add_mutually_exclusive_group(required=False)

	required.add_argument('-a', '--alignment-file', required=True,
						action=UniqueStore, help='Alignment file (BAM/CRAM); Index must be present')
	required.add_argument('-g', '--reference-genome', required=True,
						action=UniqueStore, help='Reference genome fasta file. Index must be present in same directory')
	required.add_argument('-d', '--dir-out-fq', required=True, default=".",
						help='Directory where the fastq.gz files are going to be written; default is current directory; ')
	required.add_argument('-p', '--prefix-fq', required=True, default="region_",
						help='prefix for the compressed fastq files; default is : region_  ; \
						NOTE: by default the region will be added to the prefix; so no need to add the region name or value in the prefix ; Will overwrite file if exists')

	required_mutually_exclusive.add_argument('--regions-bed-file', required=False, metavar='file',
						help='list of regions in a BED file with 1 region per row; region format: contig\tstart\tstop\tname ')
	required_mutually_exclusive.add_argument('-r', '--regions', required=False, nargs="+", action=UniqueStore, help='space-separated lists of region with each region of format: contig:start-stop ')

	igs_mutually_exclusive.add_argument('--igs-regions-bed-file', required=False, metavar='file', default=os.path.join(os.path.abspath(os.path.curdir), "data", "igs_regions.bed"),
						help='list of Igs regions (or region-like Igs considered) in a BED file with 1 region per row; region format: contig\tstart\tstop\tname ')
	igs_mutually_exclusive.add_argument('-i', '--igs-regions', required=False, nargs="+",
						help='space-separated list of Immunoglobulins regions same format as with given << --regions >>; ')

	optional.add_argument('-q', '--min-mapq', required=False, default='0', help='minimum mapping quality values for the reads to be selected; ')
	optional.add_argument('-t', '--threads', '--cpus', required=False, default=2,
						help='Number of threads or process to run in parallel; This allows processing regions in parallel to speed up read capture;')
	optional.add_argument('--min-tuples', required=False, default=3, help='Minimum Number of tuples of positions in a region; Allows filtering regions with less than a number of tuples')
	optional.add_argument('--left-position-clustering-only', required=False, default=True,
						help='For specific usage a la pairoscope, clustering of the left position only')
	optional.add_argument('--min-reads-per-cluster', required=False, default=5, help='minimum number of reads per cluster')
	optional.add_argument('--min-dist-left-right', required=False, default='0', help='filter out cluster if distance btw leftmost and rightmost reads is less than given value; \
						to enable that filtering set a value;  default is disabled (aka set to 0).')
	optional.add_argument('--clustering-distance', required=False, default=600,
						help='Define the clustering distance between two read to be grouped together in same cluster; \
						the smaller the value the more numerous amount of small islands will be capture if exist ; the bigger the value, \
						the wider and more spread across a bigger region an island will be.')

	return parser


def capture_and_return_args(args):
	"""
	Capture and returns the arguments from the parser
	:param args: args from the parser
	:type args: argParse
	:return: arguments
	:rtype: tuple
	"""
	try:
		# ------------------------------------------------------------------
		if args["alignment_file"] is not None:
			aln_file = args["alignment_file"]
			check_file(aln_file)
		if args['reference_genome']:
			ref_genome_file = args['reference_genome']
			check_file(ref_genome_file)
		# ------------------------------------------------------------------
		if args["regions"]:
			regions = args["regions"]
		elif args['regions_bed_file']:
			genes_bed_file = args['regions_bed_file']
			check_file(genes_bed_file)
			regions = get_regions_from_bed_file(genes_bed_file)
		if args['igs_regions']:
			igs_regions = args["igs_regions"]
			igs_contig_names = [str(x).split(":")[0] for x in igs_regions]
		elif args['igs_regions_bed_file']:
			igs_bed_file = args['igs_regions_bed_file']
			check_file(igs_bed_file)
			igs_regions = get_regions_from_bed_file(igs_bed_file)
			igs_contig_names = [str(x).split(":")[0] for x in igs_regions]
		# ------------------------------------------------------------------
		if args['min_mapq']:
			min_mapq = int(args["min_mapq"])
		if args["dir_out_fq"]:
			dir_out_fq = args["dir_out_fq"]
			if not os.path.exists(dir_out_fq):
				raise IsADirectoryError("Directory NOT Found; please check your input or create the Directory")
		if args["threads"] or args['cpus']:
			ncpus = int(args["threads"])
		if args["prefix_fq"]:
			prefix_out_fq = args["prefix_fq"]
		if args["min_tuples"]:
			min_tuples = int(args["min_tuples"])
		if args["left_position_clustering_only"]:
			left_position_clustering_only = bool(args["left_position_clustering_only"])
		if args["min_reads_per_cluster"]:
			min_reads_per_cluster = int(args["min_reads_per_cluster"])
		if args["min_dist_left_right"]:
			min_dist_left_right = int(args["min_dist_left_right"])
		if args["clustering_distance"]:
			clustering_distance = int(args["clustering_distance"])
		# ------------------------------------------------------------------
	except TypeError as vte:
		logger.error(vte, exc_info=True)
		exit(2)
	except FileNotFoundError as fnf:
		logger.error(fnf, exc_info=True)
		exit(2)
	except Exception as e:
		logger.error(e, exc_info=True)
		exit(2)

	return aln_file, ref_genome_file, min_mapq, dir_out_fq, prefix_out_fq, ncpus, regions, igs_regions, igs_contig_names,  \
	       min_tuples, left_position_clustering_only, min_reads_per_cluster, min_dist_left_right, clustering_distance


def get_regions_from_bed_file(regions_bed_file):
	try:
		logger.info("parsing BED: {}".format(regions_bed_file))
		list_igs_regions = []
		with open(regions_bed_file, 'r') as bf:
			for line in bf.readlines():
				if line.startswith("#") or line.startswith("@") or line == "\n":
					continue
				contig, start, stop, *_ = line.strip().split('\t')  # So far only the first Three columns are used; we discard the remaining columns; update if needed
				list_igs_regions.append("".join([contig, ":", str(start), "-", str(stop)]))
		return list_igs_regions
	except IOError as ioe:
		logger.error("{}".format(ioe), exc_info=True)
	except ValueError as ve:
		logger.error("{}".format(ve), exc_info=True)


def make_dico_from_igs_regions(igs_regions):
	dico_igs_contig_regions = defaultdict(dict)
	for region in igs_regions:
		chrom, positions = str(region).split(":")
		dico_igs_contig_regions[chrom] = {}
		start, stop = str(positions).split("-")
		dico_igs_contig_regions[chrom]['start'] = start
		dico_igs_contig_regions[chrom]['stop'] = stop
	return dico_igs_contig_regions

def get_runtime_function(my_func):
	from time import time
	@wraps(my_func)
	def wrapper(*args, **kwargs):
		te = time()
		logger.debug("calling function: {}".format(str(my_func.__name__)))
		res = my_func(*args, **kwargs)
		logger.debug("function  << {} >>  ran for {} secs with args: {}".format(my_func.__name__, str(round(time() - te, 2)), str(args)))
		return res
	return wrapper


def merge_dicos(*args):
	"""
	Merge n dictionaries of list of tuples
	:param args: dictionaries
	:type args: dictionaries
	:return: merged dictionary
	:rtype: dictionary
	"""
	dd = defaultdict(list)
	for d in args:  # we can list as many input dicos as we want here
		# for key, value in d.items():
		# 	for my_tuple in value:
		# 		dd[key].append(my_tuple)
		for key in d.keys():
			for my_tuple in d[key]:
				dd[key].append(my_tuple)
	return dict(dd)


def is_read_within_igs_regions(dico_igs_regions, chrom, start, stop, slop=500):
	"""
	Check if the read or the mate is aligned within Igs region
	:param dico_igs_regions
	:rtype: dictionary
	:param chrom: contig
	:type chrom: string
	:param start: start position
	:type start: integer
	:param stop: stop position
	:type stop: integer
	:param slop: add slop to given Igs region in case read near boundaries of Igs regions
	:type slop: integer
	:return: True of False
	:rtype: boolean
	"""

	if chrom in dico_igs_regions.keys():
		if int(dico_igs_regions[chrom]['start'])-slop <= start and stop <= int(dico_igs_regions[chrom]['stop'])+slop:
			logger.debug("IN_IGS_REGIONS")
			return True
	else:
		raise ValueError("Contig {} for Immunoglobulin not found in dico_igs_regions" % chrom)
	return False


def process_discordant_reads_with_two_mapped(var_min_mapq, obj_aln, igs_contig_names, dico_igs_contig_regions, region):
	ddr = defaultdict(list)  # ddr stands for Dictionary Discordant Reads ; The Values iare a list of tuples ddr[chrX#+#chrY#-] = [(100,20),(101,21),(102,22)]
	darp = defaultdict(list)  # darp stands for Dico_Alns_Records_Pairs: contains list of tuples of Alignment Records ; needs a triple for loop to access a aln_rec
	tot_disc_read_count = 0

	logger.debug("calling discordant reads for given Alignment file: " + obj_aln.aln_file)
	logger.debug("calling objAln.get_discordant_reads_from_region ...")
	discordant_reads = obj_aln.get_discordant_reads_from_region(region=region, slop=0, minmapq=var_min_mapq)
	if logger.level == 10:
		logger.debug(len([1 for _ in discordant_reads]))
		discordant_reads = obj_aln.get_discordant_reads_from_region(region=region, slop=0, minmapq=var_min_mapq)

	for aln_rec in discordant_reads:
		logger.debug("in for loop process_discordant_reads_with_two_mapped:  " + str(aln_rec))
		aln_rec_mate = obj_aln.alnh.mate(aln_rec)
		tot_disc_read_count += 1
		chrom, pos, ro, chrom2, pos2, ro2 = obj_aln.extract_read_info_for_dico(aln_rec)
		# ro : stands for read orientation ; chrom : chromosome (aka contig) ; pos = Position;
		logger.debug("{}:{}{}__{}:{}{}".format(chrom, pos, ro, chrom2, pos2, ro2))
		if chrom == chrom2:
			# For Now, we only deal with Translocations
			logger.debug("##### read SKIPPED in TEST chrom == chrom2 : "+str(chrom)+"=="+str(chrom2))
			continue
		if chrom2 not in igs_contig_names:
			# for the Amgen and for the current type of analysis we want to do ...
			# we want only the mate going to an Immunoglobulin region
			# so the contig name must be one of the contigs found in igs_regions list
			logger.debug("##### read SKIPPED in TEST chrom2 not in igs_contig_names -- read SKIPPED ; chrom=={}___chrom2=={}".format(str(chrom), str(chrom2)))
			continue
		elif not is_read_within_igs_regions(dico_igs_contig_regions, chrom2, pos2, pos2, slop=500):
			logger.debug("##### read SKIPPED in 'not is_read_within_igs_regions' ; chrom=={}___chrom2=={}:{}".format(str(chrom), str(chrom2), str(pos2)))
			continue
		# ------------------------------------------------------------------------------------------------
		# Enforcing the clustering to be performed only with the read positions on the gene_regions side
		# ------------------------------------------------------------------------------------------------
		enforce_reads_orientation_both_side = False  # HARDCODED for now
		if not enforce_reads_orientation_both_side:
			k = "#".join([chrom, ro])
		else:
			k = "#".join([chrom, ro, chrom2, ro2])
		v = (pos, pos2)
		logger.debug("######  {}--{} read pair PASSED; pos-pos2::{}".format(str(chrom), str(chrom2), str(v)) )
		pairs_aln_recs = (aln_rec, aln_rec_mate)

		if k in darp.keys():
			ddr[k].append(v)
			darp[k].append(pairs_aln_recs)
		else:
			ddr[k] = [v]
			darp[k] = [pairs_aln_recs]

	logger.debug("count of total discordant reads found: [ {} ] and kept after filtering: [ {} ] in region {}".format(
		str(tot_disc_read_count), str(len(darp.values())), region))

	return ddr, darp


def process_discordant_reads_with_one_unmapped(obj_aln, discordant_reads, region):
	ddr = defaultdict(list)
	darp = defaultdict(list)  # list of pairs of reads that were found in region of interest where the MATE is unmapped
	tot_disc_read_count = 0
	for aln_rec in discordant_reads:
		logger.debug("in for loop process_discordant_reads_with_one_unmapped - UnMaPPeD_CaPTuReD: " + str(aln_rec))
		aln_rec_mate = obj_aln.alnh.mate(aln_rec)  # we capture here, in this case, the counterpart mapped read from the current processed unmapped one
		tot_disc_read_count += 1  # count of unmapped reads or pairs with one read unmapped
		chrom2, pos2, ro2, chrom, pos, ro = obj_aln.extract_read_info_for_dico(aln_rec)
		k = "#".join([chrom, ro])
		# As we deal with unmapped, we set pos2 to ZERO  # HARDCODED
		v = (pos, pos2)
		logger.debug(str(k) + " <---> " + str(v))  # k == chr11#-  ; v == (10075372, 0)  ;  0 means unmapped here
		pairs_aln_recs = (aln_rec_mate, aln_rec)
		if k in darp.keys():
			ddr[k].append(v)
			darp[k].append(pairs_aln_recs)
		else:
			ddr[k] = [v]
			darp[k] = [pairs_aln_recs]
	logger.debug("count of total discordant reads found with mate unmapped : [ {} ] in region {}".format(str(tot_disc_read_count), region))

	return ddr, darp


def process_regions_parallel(var_min_mapq, tbam_file, ref_genome_file, dir_out_fq, prefix_out_fq, igs_contig_names, dico_igs_contig_regions, oclust, region):
	logger.info("processing region: " + str(region))

	obj_aln = Alignments(tbam_file, reference_genome_fasta=ref_genome_file, threads=4)

	# prepare variables for processing
	if not str(prefix_out_fq).endswith("_"):
		prefix_out_fq = prefix_out_fq + "_" + str(region).replace(":", "_").replace("-", "_").replace(" ", "_")
	else:
		prefix_out_fq = prefix_out_fq + str(region).replace(":", "_").replace("-", "_").replace(" ", "_")
	# -----------------------------------------------------------------------------------------------------
	# opening here the files to list all the fastq files generated for current region
	fbycol = open(os.path.join(dir_out_fq, prefix_out_fq + "_list_fastq_files.txt"), 'w')
	fbyitl = open(os.path.join(dir_out_fq, prefix_out_fq + "_list_fastq_files_interleaved.txt"), 'w')
	fhs = (fbycol, fbyitl)
	# -----------------------------------------------------------------------------------------------------

	ddr, darp = process_discordant_reads_with_two_mapped(var_min_mapq, obj_aln, igs_contig_names, dico_igs_contig_regions, region)

	if logger.level == 10:
		logger.debug("DARP")
		logger.debug(str(darp))
		logger.debug("DDR")
		logger.debug(str(ddr))
		for k in darp.keys():
			logger.debug("key == " + k)
			for v in darp[k]:
				logger.debug("value == " + str(len(v)) + " --> " + str(v))
				for aln in v:
					logger.debug(str(aln))
		for k in ddr.keys():
			logger.debug("key == " + k)
			for v in ddr[k]:
				logger.debug("value == " + str(len(v)) + " --> " + str(v))
				for aln in v:
					logger.debug(str(aln))

	# -----------------------------
	# Capture read pairs where the read is unmapped on the Igs region side
	# -----------------------------
	discordant_reads = obj_aln.get_discordant_reads_from_region_with_mate_unmapped(region=region, slop=0, minmapq=0, eof=True)
	ddr_unmapped, darp_unmapped = process_discordant_reads_with_one_unmapped(obj_aln, discordant_reads, region)
	logger.debug("DARP_UNMAPPED")
	logger.debug(str(darp_unmapped))
	logger.debug("DDR_UNMAPPED")
	logger.debug(str(ddr_unmapped))

	if logger.level == 10:
		for raln in darp_unmapped.values():
			for x in raln:
				for item in x:
					logger.debug(str(item))
		for elem in ddr_unmapped.values():
			for tpl in elem:
				logger.debug(str(tpl))

	# --------------------------------------------------
	# Merging dictionaries mapped and unmapped on keys
	# --------------------------------------------------
	dddr = merge_dicos(ddr, ddr_unmapped)
	ddarp = merge_dicos(darp, darp_unmapped)
	logger.debug("DDARP")
	logger.debug(str(ddarp))
	if logger.level == 10:
		for k in ddarp.keys():
			logger.debug("key == " + k)
			logger.debug("value == " + str(ddarp[k]))
			logger.debug("vtype == " + str(type(ddarp[k])))
			logger.debug("ltype == " + str(len(ddarp[k])))
			for v in ddarp[k]:
				logger.debug("value == " + str(len(v)) + " --> " + str(v))
				for aln in v:
					logger.debug(str(aln))

		for k in dddr.keys():
			logger.debug("key == " + k)
			logger.debug("value == " + str(dddr[k]))
			logger.debug("vtype == " + str(type(dddr[k])))
			logger.debug("ltype == " + str(len(dddr[k])))
			for v in dddr[k]:
				logger.debug("value == " + str(len(v)) + " --> " + str(v))
				for aln in v:
					logger.debug(str(aln))
	# ------------------------------
	# CLUSTERING SECTION
	# ------------------------------
	# we can cluster on tuple (pos,pos2) or just pos value which for Amgen makes more sense
	clustering_enable = True
	if clustering_enable:
		try:
			oclust.dico_pairs_of_start = dddr
			oclust.dico_reads_for_pairs_of_start = ddarp
			if oclust.cluster_left_pos_only:
				oclust.reformat_dico_pairs_of_start()
			d_res_clusters, darp_res_clusters = oclust.proceed_clustering()
			if logger.level == 10:
				logger.debug("oclust == {}".format(str(oclust.__dict__)))
				logger.debug("oclust == {}".format(str(oclust)))
				logger.debug(str(oclust.dico_reads_for_pairs_of_start))
				logger.debug(str(oclust.dico_pairs_of_start))
				logger.debug(str(oclust.print_dico_pos()))
				logger.debug("##"*25)
				oclust.filter_keys_with_number_tuples_pos_less_than_th()
				# oclust.print_dico_pos()
				logger.debug("@@" * 25)
				# oclust.print_dico_reads_for_pairs_of_start()
				# oclust.print_dico_clusters_by_pos()

				# oclust.print_dico_clusters_by_pos()
				# print(str(clusterings.get_list_indexes_with_read_count_ge_th(oclust.dico_clusters_by_pos, 10)))
				logger.debug("%"*30)
				logger.debug(str(d_res_clusters))
				logger.debug("%" * 30)
				logger.debug(str(darp_res_clusters))
				logger.debug("%" * 50)
				logger.debug(str(len(d_res_clusters)) + " /// " + str(len(darp_res_clusters)))
			for k in d_res_clusters.keys():
				if k not in darp_res_clusters.keys():
					logger.warning("key " + k + " NOT FOUND in darp; Something is WRONG; Check Implementation!")
				else:
					logger.debug(str(len(d_res_clusters[k])) + " /// " + str(len(darp_res_clusters[k])))
			for k in darp_res_clusters.keys():
				if k not in d_res_clusters.keys():
					logger.warning("key " + k + " NOT FOUND in d_res_clusters; Something is WRONG; Check Implementation!")
				else:
					logger.debug(str(len(d_res_clusters[k])) + " /// " + str(len(darp_res_clusters[k])))

			logger.debug(str(darp_res_clusters.keys()))
			logger.debug(str(darp_res_clusters.values()))
			logger.debug(str(darp_res_clusters.items()))

			# make fastq after clustering:
			if len(darp_res_clusters) != 0:
				logger.debug("YEAH Clusters for region")
				logger.debug(str(darp_res_clusters))
				make_fastq_files_from_cluster(dir_out_fq, prefix_out_fq, darp_res_clusters, region, fhs)
			else:
				logger.debug(str(darp_res_clusters))
				logger.debug("darp_res_clusters dico has NO Clusters for region " + region)
				make_fastq_files(dir_out_fq, prefix_out_fq + "_C00", darp_res_clusters, region, fhs)
		except Exception as e:
			logger.error(e, exc_info=True)

	else:
		# Now we have the list of read pair reads selected, we copy them into a file
		logger.info("Clustering Disabled")
		make_fastq_files(dir_out_fq, prefix_out_fq, darp, region, fhs)

	fbycol.close()
	fbyitl.close()


def get_min_max_pos(list_tpls_aln_segments):
	lpos1 = []
	for t in list_tpls_aln_segments:
		lpos1.append(t[0].pos)
	return min(lpos1), max(lpos1)


def make_fastq_files_from_cluster(dir_out_fq, prefix_out_fq, darp, region, fhs):
	"""

	:param dir_out_fq:
	:type dir_out_fq:
	:param prefix_out_fq:
	:type prefix_out_fq:
	:param darp:
	:type darp:
	:param region:
	:type region:
	:param fhs: filehandles
	:type fhs: typle of opened filehandle
	:return:
	:rtype:
	"""
	for key, dval in darp.items():
		reformat_dico = {}
		for cluster, v in dval.items():
			reformat_dico[key] = v
			# we reformatted the dictionary to get it compatible with the method `make_fastq_files(dir_out_fq, prefix_out_fq, darp, region, fhs, cluster=0)`
			# which was designed before adding clustering as a process
			min_pos, max_pos = get_min_max_pos(v)
			prefix_out_fq_nclust = "{}{}{:02n}_{:n}-{:n}".format(prefix_out_fq, "_C", cluster, min_pos, max_pos)
			logger.debug("reformat dico is: " + str(reformat_dico))
			make_fastq_files(dir_out_fq, prefix_out_fq_nclust, reformat_dico, region, fhs, cluster=cluster)


def make_fastq_files(dir_out_fq, prefix_out_fq, darp, region, fhs, cluster=0):
	"""
	From the list of reads, we make the pair-end fastq files R1 and R2
	:param dir_out_fq: output directory (full or rel path)
	:type dir_out_fq: string
	:param prefix_out_fq: prefix for the future fastq file
	:type prefix_out_fq: string
	:param darp: dictionary of tuples of pairs of reads based on keys (e.g chr4r or chr4f with r for reverse and f for forward);
	NOTE: Actually we could divide more the process here by having R1 and R2 reads as input instead of dico of reads
	:type darp: dict
	:param region: region from where the discrodant reads got extracted
	:type region: string
	:param fhs: open file handle for writing the data in the fastq file
	:type fhs: file handle
	:param cluster: optional cluster value ; default is 0 incase clustering is not performed
	:type cluster: string or int
	:return: None
	:rtype: None
	"""
	# ----------------
	# MAKE FASTQ
	# ----------------
	logger.debug("entered in make_fastq_files_function")

	def make_empty_fastq(fqr1, fqr2):
		open(fqr1, 'w').close()
		open(fqr2, 'w').close()

	if len(darp) == 0:
		logger.debug("no discordant reads found for region {}; making empty files anyway".format(region))
		# ----------------------------------------------------------------------
		# NOTE: if later we decide to make empty files when no discordant reads
		# can be found in a region, please uncomment the next 5 lines
		# ---------------------------------------------------------------------
		# expected_fq_r1_filename = prefix_out_fq + "_R1_001.fastq.gz"
		# expected_fq_r2_filename = prefix_out_fq + "_R2_001.fastq.gz"
		# file_r1 = os.path.join(dir_out_fq, expected_fq_r1_filename)
		# file_r2 = os.path.join(dir_out_fq, expected_fq_r2_filename)
		# make_empty_fastq(file_r1, file_r2)  # uncomment this lien if you want to make empty files for any region processed whether discordant reads were found or not
		logger.warning("No fastq made for region << {} >> ".format(region))
		return None
	else:
		logger.debug(region + "  darp.keys = " + str(len(darp.keys())))
		logger.debug(region + "  darp.values = " + str(len(darp.values())))
	# ----------------------------------------------------
	# WARNING: we assume this is ILLUMINA pair end reads
	# so first read must be forward strand and second read in pair must be reverse strand in the FASTQ file
	# ----------------------------------------------------
	# https://www.biostars.org/p/145590/
	# Illumina paired-end sequencing is based on the idea that you have initial DNA fragments (longer than your actual read length)\
	# and you sequence both its ends. On the Illumina chip, both ends of each sequence are amplified prior to actual sequencing using bridging. \
	# This approach results in two reads per fragment, with the first read in forward orientation and the second read in reverse-complement orientation. \
	# Depending on the initial fragment size and read length, these fragment can either overlap or not

	# -------------------------------------------
	# looping over the reformatted darp dico
	# -------------------------------------------
	lreadname = []  # list to keep record of read names that have been added to the fastq file already ; avoid having duplicates pairs in fastq
	for keyRegions in darp.keys():
		k = str(keyRegions).replace("#", "")  # example: chr4#r to chr4r
		__list_read_pairs = []
		for tuple_rp in darp[keyRegions]:  # darp[keyRegions] contain a list of Tuples with read_left of region and read_right normally in one of the Igs regions
			expected_fq_r1_filename = prefix_out_fq + "_" + k + "_R1_001.fastq.gz"
			expected_fq_r2_filename = prefix_out_fq + "_" + k + "_R2_001.fastq.gz"
			file_r1 = os.path.join(dir_out_fq, expected_fq_r1_filename)
			file_r2 = os.path.join(dir_out_fq, expected_fq_r2_filename)
			a1 = tuple_rp[0]
			a2 = tuple_rp[1]
			logger.debug(region + " <---> " + str(a1))
			logger.debug(region + " <---> " + str(a2))
			fqbk1 = None
			fqbk2 = None
			if a1.is_read1:
				if a1.is_reverse:
					fqbk1 = FastqBlock(1, a1.query_name + "/1", Seq.complement(a1.query_sequence), "".join(reversed(a1.qual)))
				else:
					fqbk1 = FastqBlock(1, a1.query_name + "/1", a1.query_sequence, a1.qual)
				if not a2.is_reverse:
					fqbk2 = FastqBlock(2, a2.query_name + "/2", Seq.complement(a2.query_sequence), "".join(reversed(a2.qual)))
				else:
					fqbk2 = FastqBlock(2, a2.query_name + "/2", a2.query_sequence, a2.qual)
			elif a2.is_read1:
				if a2.is_reverse:
					fqbk1 = FastqBlock(1, a2.query_name + "/1", Seq.complement(a2.query_sequence), "".join(reversed(a2.qual)))
				else:
					fqbk1 = FastqBlock(1, a2.query_name + "/1", a2.query_sequence, a2.qual)
				if not a1.is_reverse:
					fqbk2 = FastqBlock(2, a1.query_name + "/2", Seq.complement(a1.query_sequence), "".join(reversed(a1.qual)))
				else:
					fqbk2 = FastqBlock(2, a1.query_name + "/2", a1.query_sequence, a1.qual)
			# if we wanted to use the method a1.query_alignment_qualities as is we would need to convert to phred scale value
			# quality = ''.join(map(lambda x: chr( x+33 ), read.query_qualities)) ## extracted from pysam github
			__list_read_pairs.append(ReadPair(fqbk1, fqbk2))

			if os.path.exists(file_r1):
				logger.info("removing existing r1 file: " + expected_fq_r1_filename)
				os.remove(file_r1)
			if os.path.exists(file_r2):
				logger.debug("removing existing r2 file: " + expected_fq_r2_filename)
				os.remove(file_r2)

		logger.debug("calling making fastq.gz files ...")
		logger.debug("in " + __name__ + " fullpath file r1: " + str(file_r1) + "  <--> file r2: " + str(file_r2))

		try:
			if len(__list_read_pairs) == 0:
				make_empty_fastq(file_r1, file_r2)
			else:
				with gzip.open(file_r2, "w") as f2, gzip.open(file_r1, "w") as f1:
					for o_paired_reads in __list_read_pairs:
						if o_paired_reads.r1.readname in lreadname:
							logger.debug("R1 already seen once so we skip its addition to the fastq: " + str(o_paired_reads.r1.readname))
							continue
						lreadname.append(o_paired_reads.r1.readname)
						str1_out = "@{} {}\n{}\n+\n{}\n".format(str(o_paired_reads.r1.readname),
						                                        str(o_paired_reads.r1.comment),
						                                        str(o_paired_reads.r1.sequence),
						                                        str(o_paired_reads.r1.rquals))
						str2_out = "@{} {}\n{}\n+\n{}\n".format(str(o_paired_reads.r2.readname),
						                                        str(o_paired_reads.r2.comment),
						                                        str(o_paired_reads.r2.sequence),
						                                        str(o_paired_reads.r2.rquals))
						logger.debug("####"*50)
						logger.debug("str1_out : "+str1_out)
						logger.debug("str2_out : " + str2_out)
						logger.debug("####" * 50)

						# writing to file here
						f2.write(str2_out.encode('utf-8'))
						f1.write(str1_out.encode('utf-8'))
			fhs[0].write("{}\t{}\n".format(file_r1, file_r2))
			fhs[1].write("{}\n{}\n".format(file_r1, file_r2))
		except IOError as ioe:
			logger.error(ioe, exc_info=True)
	logger.info("END making fastq file for region " + str(region) + ", key " + k + ", and cluster " + str(cluster))


class Alignments:
	def __init__(self, alnfile, reference_genome_fasta=None, threads=2, slop=0):
		"""

		:param alnfile: BAM, CRAM or SAM alignment file
		:type alnfile: string or pathname
		:param reference_genome_fasta: filename for the reference genome fasta file
		:type reference_genome_fasta: string or pathname
		:param threads: use this much threads when reading CRAM file; not use when SAM or BAM
		:type threads: integer
		"""
		self.aln_file = alnfile
		self.reference_genome_fasta = reference_genome_fasta
		self.threads = threads
		self.slop = slop
		self.alnh = self.read_aln_file()  # bam handle or cram handle
		self.dict_contig_length = self.make_dico_contig_lengths_from_bam_handler()

	def read_aln_file(self):
		"""
		read the alignment file whether it is a SAM, BAM or CRAM file and returns the bam file handle
		:return: aln read file handle (bamh or alnh)
		"""

		extension = os.path.splitext(self.aln_file)[1]
		try:
			if extension == ".cram":
				if self.reference_genome_fasta is None:
					raise FileNotFoundError("ERROR: reading CRAM file requires a Reference Genome Fasta File To be Provided with its FAI index.")
				return pysam.AlignmentFile(self.aln_file, mode='rc', reference_filename=self.reference_genome_fasta, threads=self.threads)
			elif extension == ".bam":
				return pysam.AlignmentFile(self.aln_file, mode='rb', threads=self.threads)
			elif extension == ".sam":
				return pysam.AlignmentFile(self.aln_file, mode='r')
			else:
				logger.debug("extension found: " + extension)
				raise Exception("EXPECTED EXTENSION for ALIGNMENT FILE NOT FOUND; must be either .cram, .bam or .sam")
		except FileNotFoundError as fnf:
			logger.exception(fnf)
			exit(2)
		except Exception as e:
			logger.exception(e)
			exit(2)

	def make_dico_contig_lengths_from_bam_handler(self):
		"""
		create a dictionary of contig with their length as values
		:param self: aln handle from pysam bam/cram reader
		:return: dictionary of contigs:lengths
		"""
		return dict(zip(self.alnh.references, self.alnh.lengths))

	def get_contig_length(self, contig_name):
		"""
		capture the contig and their length from a dictionary
		:param self.dict_contig_length: dictionary of contig with their length
		:param contig_name: contig name for which we want the length
		:return: length of the contig
		"""
		if contig_name in self.dict_contig_length.keys():
			return int(self.dict_contig_length[contig_name])
		else:
			raise KeyError("ERROR: Contig name << {} >> NOT FOUND ".format(contig_name))

	def is_contig_exist(self, contig):
		"""
		Check if a contig exists in the bam file ; if not raises error and exits ; if yes, do nothing.
		:param contig: contig name such as 'chr4' or 'X'
		:type contig: string
		:return: None
		:rtype: None
		"""
		if contig not in self.dict_contig_length.keys():
			raise ValueError("Given Contig Name <<{}>> not Found In Alignment file; Check your input".format(contig))

	def is_given_region_ok(self, region):
		if ":" not in region:  # we assume that the user wants read for the whole contig
			contig = region
			self.is_contig_exist(contig)
			return region, 0, int(self.get_contig_length(contig))
		contig = region.split(":")[0]
		self.is_contig_exist(contig)
		start, stop = region.split(":")[1].split("-")
		assert isinstance(int(start), int), "start value for region MUST be an INTEGER"
		assert isinstance(int(stop), int), "stop value for region MUST be an INTEGER"
		return contig, int(start), int(stop)

	@get_runtime_function
	def fetch_reads_in_region(self, region, slop=0):
		"""
		Fetch All the reads in the given region
		Region MUST be of the format : contigName:Start-Stop ;
		Stop is Mandatory ; If you only want reads overlapping 1 base, set Start and Stop to the same value
		:param region: region of an existing contig ; region format can be region=contig  OR region=contig:start-stop ; No other format allowed
		:type region: basestring
		:param slop: extra bases given to a region of interest if needed, default is 0 ; often use when capturing variant position automatically
		:type slop: integer
		:return: Pysam Fetch Object (generator)
		:rtype: generator
		"""
		contig, start, stop = self.is_given_region_ok(region)
		if slop != 0 and not isinstance(slop, int):
			raise ValueError("slop value in {} MUST be in Integer; Found {}".format(self.fetch_reads_in_region.__name__, str(slop)))
		if stop+slop > self.get_contig_length(contig):
			stop = self.get_contig_length(contig)
			logger.warning("'stop+slop'value went higher than contig length; stop value got adjusted to {}".format(str(stop)))
		return self.alnh.fetch(contig=contig, start=max(start-slop, 1), end=stop+slop)

	@get_runtime_function
	def get_discordant_reads_from_region(self, region, slop=0, minmapq=0):
		"""
		Fetch All the reads in the given region; We do not filter the reads here ; this is done when parsing the generator
		Region MUST be of the format : contigName:Start-Stop ;
		Stop is Mandatory ; If you only want reads overlapping 1 base, set Start and Stop to the same value
		:param region: region of an existing contig
		:type region: basestring
		:param slop: extra bases given to a region of interest if needed, default is 0 ; often use when capturing variant position automatically
		:type slop: integer
		:param minmapq: minimum mapping quality value for read to be kept
		:type minmapq: integer
		:return: Pysam Fetch Object (generator)
		:rtype: generator
		"""

		contig, start, stop = self.is_given_region_ok(region)
		logger.debug("{} : {} -- {}".format(str(contig), str(start), str(stop)))
		if slop != 0 and not isinstance(slop, int):
			raise ValueError("slop value in {} MUST be in Integer; Found {}".format(self.fetch_reads_in_region.__name__, str(slop)))
		if stop+slop > self.get_contig_length(contig):
			stop = self.get_contig_length(contig)-slop
			logger.warning("'stop+slop' value went higher than contig length; stop value got adjusted to {}".format(str(stop)))
		logger.debug("{} : {} -- {}  (after_slop_adjustment if any)".format(str(contig), str(start), str(stop)))
		fetched_aln_record = self.alnh.fetch(contig=contig, start=max(start - slop, 1), end=stop + slop)

		for x in fetched_aln_record:
			logger.debug("testing if we should yield the following PAIR: " + str(x))
			mq_mate = self.get_mate_mq_value(x)
			mq = x.mapping_quality
			logger.debug("Mapping Qualities for both read and it mate: "+str(mq) +" and " + str(mq_mate))
			# and not x.is_supplementary \
			# and not x.is_secondary \ ## is_secondary means true if not primary alignment in Pysam language
			if not x.is_proper_pair or (x.is_proper_pair and x.is_supplementary):
				logger.debug("PASSED PROPER PAIR test --- testing if we should yield the following PAIR: " + str(x))
				if not x.is_unmapped \
						and not x.is_qcfail \
						and not x.is_duplicate \
						and not x.is_secondary \
						and not x.mate_is_unmapped:
					if mq >= minmapq and mq_mate >= minmapq:
						logger.debug("in YIELD ----  testing if we should yield the following PAIR: " + str(x))
						yield x

	def get_discordant_reads_from_region_with_mate_unmapped(self, region, slop=0, minmapq=0, eof=False, pct_n=0.20):
		"""
		Fetch All the reads in the given region; We do not filter the reads here ; this is done when parsing the generator
		Region MUST be of the format : contigName:Start-Stop ;
		Stop is Mandatory ; If you only want reads overlapping 1 base, set Start and Stop to the same value
		:param region: region of an existing contig
		:type region: basestring
		:param slop: extra bases given to a region of interest if needed, default is 0 ; often use when capturing variant position automatically
		:type slop: integer
		:param minmapq: minimum mapping quality value for read to be kept
		:type minmapq: integer
		:param eof: In order to capture teh UNMAPPED reads we need this option
		:type eof: boolean
		:param pct_n: percentage of maximum N character in the query sequence we allow to pass filter
		:type pct_n: float
		:return: Pysam Fetch Object (generator)
		:rtype: generator
		"""

		contig, start, stop = self.is_given_region_ok(region)
		logger.debug("{} : {} -- {}".format(str(contig), str(start), str(stop)))
		if slop != 0 and not isinstance(slop, int):
			raise ValueError("slop value in {} MUST be in Integer; Found {}".format(self.fetch_reads_in_region.__name__, str(slop)))
		if stop+slop > self.get_contig_length(contig):
			stop = self.get_contig_length(contig)-slop
			logger.warning("'stop+slop' value went higher than contig length; stop value got adjusted to {}".format(str(stop)))
		logger.debug("{} : {} -- {}  (after_slop_adjustment if any)".format(str(contig), str(start), str(stop)))
		fetched_aln_record = self.alnh.fetch(contig=contig, start=max(start - slop, 1), end=stop + slop, until_eof=eof)

		# we capture only the unmapped read for which the mate is mapped within the current region
		for x in fetched_aln_record:
			# we do not need the mapping quality for unmapped reads so we set it at zero, but better off
			if x.is_unmapped and not x.is_qcfail and not x.mate_is_unmapped and x.mapping_quality >= minmapq:
				# we check if the unmapped read has less than a certain percentage of N value, otherwise if only Ns what is the point?
				if x.query_sequence.count("N")/len(x.query_sequence) <= pct_n:
					yield x

	def is_paired(self, alignments=1000):
		"""
		check if an *alignment file* contains paired end data.
		The method reads at most the first *alignments* and returns
		True (n!=alignments) if any of the alignments are paired.
		# from source: https://www.programcreek.com/python/example/90614/pysam.AlignmentFile for isPAired and estimateSizeDistribution
		"""
		n = 0
		for read in self.alnh:
			if read.is_paired:
				break
			n += 1
			if n == alignments:
				break
		return n != alignments


	@staticmethod
	@get_runtime_function
	def count_reads_in_fetched_region(pysam_fetch_obj):
		"""
		from the pysam's fetch obj we count the number of total unfiltered reads found in given region
		warning: the generator fetch is consumed here; user must re-create it with
		:return: total_number_of_reads
		:rtype: integer
		"""
		return sum(1 for _ in pysam_fetch_obj)

	@staticmethod
	def get_mate_mq_value(aln_rec):
		"""
		check if MQ exists in aln_rec or assign default passing value to MQ flag and return value
		:param aln_rec:
		:return: mq_mate value
		"""
		if aln_rec.has_tag('MQ'):  # some alignment records do not have the MQ value in the mate when pair aligned; BWA-known-issue
			mq_mate = aln_rec.get_tag('MQ')
		else:
			mq_mate = 60  # we assume that if the mate quality is not present, it is still a good one.
		return mq_mate

	@staticmethod
	def extract_read_info_for_dico(aln_rec):
		"""
		Manage data extraction from pair of reads whether both reads are mapped or only one of the two is mapped
		:param aln_rec: AlignSegment or Alignment record or Read record
		:type aln_rec: AlignSegment
		:return: chromosome, position, read_orientation, chromosome_of_mate, position_of_mate, read_orientation_of_mate
		:rtype: tuple of the return values above
		"""
		r = aln_rec
		if r.cigarstring is not None:
			chrom = r.reference_name
			pos = r.reference_start
			# ro = "-" if r.is_reverse else "+"
			ro = "r" if r.is_reverse else "f"
			chrom2 = r.next_reference_name
			pos2 = r.next_reference_start
			# ro2 = "-" if r.mate_is_reverse else "+"
			ro2 = "r" if r.mate_is_reverse else "f"
		else:
			chrom = "none"
			pos = 0
			ro = "none"
			chrom2 = r.next_reference_name
			pos2 = r.next_reference_start
			# ro2 = "-" if r.mate_is_reverse else "+"
			ro2 = "r" if r.mate_is_reverse else "f"

		return chrom, pos, ro, chrom2, pos2, ro2


class ReadPair:
	"""
	Object that gather the information about a pair of reads R1 and R2
	"""
	def __init__(self, r1, r2):
		"""

		:param r1: aln_rec
		:type r1: Pysam Aligned Segment
		:param r2: Mate record of r1 (should be unique)
		:type r2: Pysam Aligned Segment
		"""
		self.r1 = r1
		self.r2 = r2


class FastqBlock:
	"""
	### READ_1
	@A00674:63:HLCMKDMXX:1:1101:1371:1016 1:N:0:GCCAAT
	GNCGTAATGGGGCTCTGAAGGTAGCTCAGCGCACTGTTCGTGGAGTGGACAGGGATGGAGAC
	+
	F#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

	### READ_2 (or mate)
	@A00674:63:HLCMKDMXX:1:1101:1371:1016 2:N:0:GCCAAT
	CTTCAGCTTCATGGGCCAGCTGCTGCAGTTTGAGTCCCAGGTGCTGGCTCCGCACTGTTCGG
	+
	FFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

	"""
	def __init__(self, read_number, readname, sequence, rquals, comment=""):
		self.read_number = read_number  # 1 or 2
		self.readname = readname
		self.sequence = sequence
		self.rquals = rquals
		self.comment = comment
