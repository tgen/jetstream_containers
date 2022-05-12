import argparse
import logging.config
from os import path

dir_root = path.dirname(path.dirname(__file__))
logging.config.fileConfig(path.join(dir_root, 'logging_config.ini'), disable_existing_loggers=False)
logger = logging.getLogger(__name__)


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


def make_parser_args():
	
	dir_root = path.dirname(path.dirname(__file__))
	
	parser = argparse.ArgumentParser(add_help=True, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser._action_groups.pop()
	required = parser.add_argument_group('###-- required arguments --###')
	optional = parser.add_argument_group('###-- optional arguments --###')

	required.add_argument('-s', '--seg', required=True, action=UniqueStore,
						  help='SEG file with non-overlapping segments')
	required.add_argument('-g', '--genomic-regions', required=True, action=UniqueStore,
						  help='Genomic Regions File defining Centromere (required) and optionally telomeres, p and q arms. NOTE: defined regions MUST not overlap')

	optional.add_argument('-t', '--th-log2r', required=False, default=-0.1613,
						  help=' log2R Threshold value for considering a Deletion ')
	optional.add_argument('-m', '--minsize', required=False, default=10 ** 6,
						  help=' Minimum size region for a Segment to be considered into HRD score')
	optional.add_argument('-j', '--th-pct-overlapping', required=False, default=0.90,
						  help=' Percentage of overlapping between arm and sum of segements with deletions [allow to exclude deletion being an entire arm for instance ]')
	optional.add_argument('-d', '--dir-out', required=False, default=".",
						  help=' Output directory where the temp and results file will be written (default is current directory)')
	optional.add_argument('-o', '--outfile', required=False, default="hrd_scores.txt",
						  help=' output file where the score HRD will be saved into; Can be rlative or full path')
	optional.add_argument('--threads', required=False, default=1,
						  help=' Number of cpus or threads')
	optional.add_argument('-w', '--contigs', required=False, default=".",
						  help=' comma-separated list of contigs to use with the HRD score; Contigs Must Exist in Seg file and Genomic Region file; If None, list will be captured from SEG file')
	optional.add_argument('-S', '--sample', required=False, default=".",
						  help=' Sample Name to add to the results output table ')
	optional.add_argument('-i', '--id', required=False, default=".",
						  help=' Kit Code or Any ID you want to add to the results output table ')
	optional.add_argument('-e', '--exclude-contigs', required=False, default="chrX,chrY,chrM",
						  help=' comma-separated list of contigs to exclude from HRD scores ; default is "chrX,chrY,chrM" ')
	optional.add_argument('-k', '--karyo-file', required=False, default=path.join(dir_root,"examples/inputs/grch38_ucsc_cytobands.bed"),
	                      help=' Cytoband Information of the Genome [UCSC --> table_browser --> Sequence & Mapping --> Chromosome Band --> allFieldsFromSelectedtable --> '
	                           'getOuput button ; default is karyo_file = '
	                           '"examples/inputs/grch38_ucsc_cytobands.bed" ')
	optional.add_argument('-p', '--plots', required=False, default=False, action='store_true', help=' by default no plot will be created ; use --plots to enabled making plots')

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
		if args["seg"] is not None:
			seg_file = args["seg"]
			check_file(seg_file)
		if args['genomic_regions']:
			ref_genome_file = args['genomic_regions']
			check_file(ref_genome_file)
		# ------------------------------------------------------------------
		# ------------------------------------------------------------------
		if args["dir_out"]:
			dir_out = args["dir_out"]
			if not path.exists(dir_out):
				raise IsADirectoryError(
					"Directory NOT Found; please check your input or create the Directory {}".format(dir_out))
		if args["threads"] or args['cpus']:
			ncpus = int(args["threads"])
			check_cpus(ncpus)
		if args['th_log2r']:
			th_log2r = float(args["th_log2r"])
			check_th_log2r(th_log2r)
		if args["minsize"]:
			minsize = int(args["minsize"])
			check_minsize(minsize)
		if args["contigs"]:
			list_contigs = args["contigs"]
			if list_contigs != ".":
				list_contigs = [str(x).strip() for x in list_contigs.split(",")]
			else:
				list_contigs = []
			if not isinstance(list_contigs, list):
				raise TypeError("list of contigs provided is MALFORMED or of the wrong type; Check your input")
		if args["th_pct_overlapping"]:
			th_pct_overlapping = float(args["th_pct_overlapping"])
		if args["outfile"]:
			outfile = args["outfile"]
		if args["sample"]:
			samplename = args["sample"]
		if args["id"]:
			id = args["id"]
		if args["exclude_contigs"]:
			list_exclude_contigs = args["exclude_contigs"]
		if args["karyo_file"]:
			check_file(path.abspath(args["karyo_file"]))
			karyo_file = path.abspath(args["karyo_file"])

		# Arguments with FALSE or TRUE as Default value
		make_plots = args['plots']

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

	return seg_file, ref_genome_file, dir_out, ncpus, th_log2r, \
		   minsize, list_contigs, th_pct_overlapping, outfile, \
		   samplename, id, list_exclude_contigs, karyo_file, make_plots



def print_recap_inputs(tpl_args):
	logger.info("")
	logger.info("Inputs captured:")
	logger.info(">>> SEG file = " + tpl_args[0])
	logger.info(">>> Genomic Regions File = " + tpl_args[1])
	logger.info(">>> Dir_Out = " + str(tpl_args[2]))
	logger.info(">>> Threads = " + str(tpl_args[3]))
	logger.info(">>> Threshold Log2R = " + str(tpl_args[4]))
	logger.info(">>> Threshold minSize Segment = " + str(tpl_args[5]))
	logger.info(">>> Comma-separated list_contigs provided by User = " + str(tpl_args[6]))
	logger.info(">>> Threshold percentage overlapping = " + str(tpl_args[7]))
	logger.info(">>> OutFile = " + str(tpl_args[8]))
	logger.info(">>> samplename = " + str(tpl_args[9]))
	logger.info(">>> id = " + str(tpl_args[10]))
	logger.info(">>> list_exclude_contigs = " + str(tpl_args[11]))
	logger.info(">>> Karyo file = " + str(tpl_args[12]))
	logger.info(">>> plots = " + str(tpl_args[13]))
	logger.info("")


def check_file(f):
	if f is None or not path.isfile(f) or not path.exists(f):
		logger.error("FNF: " + str(f))
		raise FileNotFoundError("FNF: " + str(f))


def check_minsize(minsize):
	if minsize <=0:
		raise ValueError("ERROR: minimum size for segment length must be positive and greater than zero")
	elif minsize <=10000:
		logger.warning("WARNING: min size is less than ten thousands; you might grab false positive HRD loci with that value")


def check_th_log2r(th_log2r):
	if th_log2r >= 0:
		raise ValueError("ERROR: The threshold entered is positive and therefore does not suit the purpose "
		                 "of that tool; The threshold MUST be a Negative Float Value; Recommended values is equal to or less than -0.1613")


def check_cpus(ncpus):
	if ncpus < 1:
		raise ValueError("ERROR: Number of CPUs or Threads MUST be greater than 0; Aborting")
