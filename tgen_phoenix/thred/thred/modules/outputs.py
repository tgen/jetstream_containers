import logging
from os import path, chdir
from modules import karyoplot
from collections import defaultdict

dir_root = path.dirname(path.dirname(__file__))
logging.config.fileConfig(path.join(dir_root, 'logging_config.ini'), disable_existing_loggers=False)
logger = logging.getLogger(__name__)


def change_dir(curdir_path, newdir_path):
	"""
	modify the working directory from the current to new given directory
	:param curdir_path: current directory
	:type curdir_path: basestring
	:param newdir_path: new target directory
	:type newdir_path: basestring
	:return: Nothing to return; the try statement will return error if needed and nothing is all OK
	:rtype: None
	"""
	logger.debug("old_dir: " + curdir_path)
	if newdir_path.startswith("/"):
		chdir(path=newdir_path)
	elif newdir_path.startswith("./") or newdir_path.startswith("../"):
		try:
			chdir(path=path.join(curdir_path, newdir_path))
		except IOError:
			from os import mkdir
			mkdir(newdir_path)
			chdir(curdir_path, newdir_path)
	else:
		try:
			chdir(path=newdir_path)
		except IOError:
			from os import mkdir
			mkdir(newdir_path)
			chdir(newdir_path)
	logger.debug("new_dir: " + str(path.abspath(path.curdir)))


def write_output_scores(scores, dirout=".", outfile=None, sample=".", id=".", default_outfilename = "hrd_scores.txt"
	):
	"""
	write HRD score and additional information (customizable upon request)
	:param scores: HRD score and also the sum of the segments isze involved in the calculus of the score
	:type scores: tuple
	:param dirout:  directory where the outfile will be written if the given outfile represents only a filename without any path associated to it
	:type dirout: basestring
	:param outfile: outfile where scores and miscs infos are written; can be relative or full path
	:type outfile: basestring
	:param sample: name of the processed sample
	:type sample: basestring
	:param id: represents an ID value; any ID can be used ; User's Choice
	:type id: basestring
	:return: write data into a file
	:rtype: None
	"""
	if outfile is None:
		outfile = default_outfilename
	else:
		outfile = outfile + "_" + default_outfilename
	ori_curdir = path.abspath(path.curdir)
	change_dir(ori_curdir, dirout)
	logger.info("HRD score(s) File can be found at {}".format(path.abspath(outfile)))
	with open(path.join(outfile), 'w') as of:
		header = "\t".join(["#sample", "id", "sizeFilteredSegments", "sizeGenomeTerritory", "BreaksCount", "BreaksPerMb", "HRDscore"])
		line_score = "\t".join([sample, id, "\t".join([str(x) for x in scores])])
		of.write(header + "\n")
		of.write(line_score + "\n")
	change_dir(path.abspath(path.curdir), ori_curdir)


def write_new_territory(dico_new_terr, dirout=".", outfile=None, default_outfilename = "hrd_captured_genome_territory.txt"):
	"""
	write genome territory used for calculated to the given outfile
	:param dico_new_terr: dictionary of contigs with list of REGIONS for each contig representing the genome territory used to calculate teh HRDscore
	:type dico_new_terr: dict
	:param dirout: directory where the outfile will be written if the given outfile represents only a filename without any path associated to it
	:type dirout: basestring
	:param outfile: outfile where filtered segments used to calculate HRDscore are written; can be relative or full path
	:type outfile: basestring
	:return: write the genome territory used into the outfile
	:rtype: None
	"""
	if outfile is None:
		outfile = default_outfilename
	else:
		outfile = outfile + "_" + default_outfilename
	ori_curdir = path.abspath(path.curdir)
	change_dir(ori_curdir, dirout)
	logger.info("File with genome Territory used for HRD score(s) calculi can be found at: {}".format(path.abspath(outfile)))
	with open(outfile, 'w') as of:
		header = "\t".join(["#contig", "start", "end", "region_name", "length"])
		of.write(header + "\n")
		for contig, dico_arm in dico_new_terr.items():
			for reg in dico_arm.values():  # dico_arm.values() contains only ONE region per arm since we should have only one region representing the p or q arm NEW territory
				logger.debug(str(type(reg)))
				of.write(str(reg) + "\n")
	# back to where we were before writing output
	change_dir(path.abspath(path.curdir), ori_curdir)


def write_output_filtered_segments(dico_flt_segments, dirout=".", outfile=None, default_outfilename = "hrd_flt_segments.txt"):
	"""
	write filtered segments to the given outfile
	:param dico_flt_segments: dictionary of contigs with list of segments for each contig
	:type dico_flt_segments: dict
	:param dirout: directory where the outfile will be written if the given outfile represents only a filename without any path associated to it
	:type dirout: basestring
	:param outfile: outfile where filtered segments used to calculate HRDscore are written; can be relative or full path
	:type outfile: basestring
	:return: write filtered segments into the outfile
	:rtype: None
	"""
	if outfile is None:
		outfile = default_outfilename
	else:
		outfile = outfile + "_" + default_outfilename
	ori_curdir = path.abspath(path.curdir)
	change_dir(ori_curdir, dirout)
	logger.info("File with Filtered Segments used to calculate HRD score(s) can be found at {}".format(path.abspath(outfile)))
	with open(outfile, 'w') as of:
		header = "\t".join(["Sample", "Chromosome",	"Start", "End", "Num_Probes", "Segment_Mean", "Segment_Length", "arm"])
		of.write(header + "\n")
		for contig, dico_arm in dico_flt_segments.items():
			for arm, lsegs in dico_arm.items():
				for seg in lsegs:
					logger.debug(str(type(seg)))
					of.write(str(seg) + "\t" + arm + "\n")
	# back to where we were before writing output
	change_dir(path.abspath(path.curdir), ori_curdir)


def write_output_original_segments(dico_ori_segments, dirout=".", outfile=None, default_outfilename = "hrd_ori_segments.txt"):
	"""
	write filtered segments to the given outfile
	:param dico_flt_segments: dictionary of contigs with list of segments for each contig
	:type dico_flt_segments: dict
	:param dirout: directory where the outfile will be written if the given outfile represents only a filename without any path associated to it
	:type dirout: basestring
	:param outfile: outfile where filtered segments used to calculate HRDscore are written; can be relative or full path
	:type outfile: basestring
	:return: write filtered segments into the outfile
	:rtype: None
	"""
	if outfile is None:
		outfile = default_outfilename
	else:
		outfile = outfile + "_" + default_outfilename
	ori_curdir = path.abspath(path.curdir)
	change_dir(ori_curdir, dirout)
	logger.info("File with Filtered Segments used to calculate HRD score(s) can be found at {}".format(path.abspath(outfile)))
	with open(outfile, 'w') as of:
		header = "\t".join(["Sample", "Chromosome",	"Start", "End", "Num_Probes", "Segment_Mean", "Segment_Length"])
		of.write(header + "\n")
		for contig, lsegs in dico_ori_segments.items():
			for seg in lsegs:
				logger.debug(str(type(seg)))
				of.write(str(seg) + "\n")
	# back to where we were before writing output
	change_dir(path.abspath(path.curdir), ori_curdir)

def make_karyoplot(kf, plotname, th_log2r_deletion, dico_segments, part=1, dirout="."):
	"""
	Make Simple Karyoplots and add to the plots the location of the SEGMENTS provided by the given dictionary
	:param kf: Need the Karyotype file representing the Cytobands; See Module header for further details about the file
	:type kf: basestring
	:param plotname: name of the png image file that wil be created
	:type plotname: basestring
	:param th_log2r_deletion:
	:type th_log2r_deletion: float
	:param dico_segments:
	:type dico_segments: dict
	:param part:
	:type part: integer
	:param dirout:
	:type dirout: basestring
	:return: make the png file but nothing is returned
	:rtype: None
	"""
	ori_curdir = path.abspath(path.curdir)
	change_dir(ori_curdir, dirout)
	
	logger.debug("karyo_file is: "+ kf)
	logger.debug(dico_segments)
	dico_metadata = defaultdict(dict)
	for k, d_arms in dico_segments.items():
		dico_metadata[k] = dict()
		for arm, lsegs in d_arms.items():
			if not lsegs:
				continue
			dico_metadata[k][arm] = zip([ int(seg.start) for seg in lsegs ], [ float(seg.mean) for seg in lsegs ], [ int(seg.end) for seg in lsegs ] )
	logger.debug(dico_metadata)
	karyoplot.karyoplot(kf, plotname, th_log2r_deletion, metadata=dico_metadata, part=part)
	# back to where we were before writing output
	change_dir(path.abspath(path.curdir), ori_curdir)