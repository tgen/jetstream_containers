
# Karyoplot with Python...
# SOURCE: https://gist.github.com/kantale/e390cf7a47c4afdff9e4
# Modified by Christophe Legendre for tHReD tool

from os import linesep, environ, path
import logging
from matplotlib.patches import Wedge, Rectangle
import matplotlib.pyplot as plt

dir_root = path.dirname(path.dirname(__file__))
logging.config.fileConfig(path.join(dir_root, 'logging_config.ini'), disable_existing_loggers=False)
logger = logging.getLogger(__name__)


def karyoplot(karyo_filename, samplename, th_log2r_deletion=-0.1613, metadata={}, part=1):
	"""
	To create a karyo_filename go to: http://genome.ucsc.edu/cgi-bin/hgTables
	group: Mapping and Sequencing
	track: Chromosome Band
	[ UCSC --> table_browser --> Sequence & Mapping --> Chromosome Band --> allFieldsFromSelectedtable --> getOuput button ]

	An example of an output (hg19, Human) is here: http://pastebin.com/6nBX6sdE

	The script will plot dots next to loci defined in metadata as:
	metadata = {
		'1' : [2300000, 125000000, 249250621],
	}
	"""

	# -------------------------------------------------------------------------------------
	# check if chromosome in metadata are prefixed with chr; if so; we remove the prefix
	# -------------------------------------------------------------------------------------
	new_metadata = {}
	for k in metadata.keys():
		if "chr" in k:
			new_metadata[k.replace("chr", "")] = metadata[k]
	del metadata
	# metadata.update(new_metadata)
	metadata = new_metadata

	karyo_dict = {}
	
	with open(karyo_filename, 'r') as karyo_f:
		lines = [x.replace(linesep, '').split() for x in karyo_f.readlines()]
		# TODO: Currently the number of Chromosome is tight to the number of chromosome in Human; Make it Generic by adding parameter to function to take a list of contig
		# TODO: that list of contig has to be converted to be compatible in the line below : 1) remove chr, 2) make it a variable
		for chromosome in [str(x) for x in range(1, 23)] + ['X', 'Y']:
			karyo_dict[chromosome] = [[y[0], int(y[1]), int(y[2]), y[3], y[4]] for y in [x for x in lines if x[0] == 'chr' + chromosome]]

	fig, ax = plt.subplots()

	# original value: DIM = 1.0 () not working with our implementation
	DIM = 1.0

	ax.set_xlim([0.0, DIM * (1.3)])
	ax.set_ylim([0.0, DIM])

	def get_chromosome_length(chromosome):
		chromosome_start = float(min([x[1] for x in karyo_dict[chromosome]]))
		chromosome_end = float(max(x[2] for x in karyo_dict[chromosome]))
		chromosome_length = chromosome_end - chromosome_start
		return chromosome_length

	def plot_chromosome(chromosome, order):
		chromosome_length = get_chromosome_length(chromosome)
		chromosome_length_1 = get_chromosome_length('1')

		x_start = order * DIM * 0.1
		x_end = x_start + (DIM * 0.04)
		y_start = DIM * 0.8 * (chromosome_length/chromosome_length_1)
		y_end = DIM * 0.1

		# We use the same colors as: http://circos.ca/tutorials/lessons/2d_tracks/connectors/configuration
		colors = {
			'gpos100' : (0/255.0,0/255.0,0/255.0),
			'gpos'    : (0/255.0,0/255.0,0/255.0),
			'gpos75'  : (130/255.0,130/255.0,130/255.0),
			'gpos66'  : (160/255.0,160/255.0,160/255.0),
			'gpos50'  : (200/255.0,200/255.0,200/255.0),
			'gpos33'  : (210/255.0,210/255.0,210/255.0),
			'gpos25'  : (200/255.0,200/255.0,200/255.0),
			'gvar'    : (220/255.0,220/255.0,220/255.0),
			'gneg'    : (255/255.0,255/255.0,255/255.0),
			'acen'    : (217/255.0,47/255.0,39/255.0),
			'stalk'   : (100/255.0,127/255.0,164/255.0),
		}

		for index, piece in enumerate(karyo_dict[chromosome]):

			current_height = piece[2] - piece[1]
			current_height_sc = ((y_end - y_start) / chromosome_length) * current_height
			
			if index == 0:
				y_previous = y_start

			y_next = y_previous + current_height_sc

			color = colors[piece[4]]
			# ------------------------------
			# plot the Karyotypes
			# ------------------------------
			r = Rectangle((x_start, y_previous), x_end-x_start, current_height_sc, color=color)
			ax.add_patch(r)

			y_previous = y_next
		
		# ---------------------------------------------------------------
		# Plot semicircles at the beginning and end of the chromosomes
		# --------------------------------------------------------------
		center_x = x_start + (x_end-x_start)/2.0
		radius = (x_end-x_start)/2.0
		theta1 = 0.0
		theta2 = 180.0
		w1 = Wedge((center_x, y_start), radius, theta1, theta2, width=0.00001, facecolor='white', edgecolor='black')
		w2 = Wedge((center_x, y_end), radius, theta2, theta1, width=0.00001, facecolor='white', edgecolor='black')
		ax.add_patch(w1)
		ax.add_patch(w2)
		ax.plot([x_start, x_start], [y_start, y_end], ls='-', color='black')
		ax.plot([x_end, x_end], [y_start, y_end], ls='-', color='black')

		# ------------------------------
		# Plot metadata
		# -----------------------------
		th_log2r_deletion_pos = abs(th_log2r_deletion)
		th_log2r_deletion_neg = -th_log2r_deletion_pos
		display_char = "|"
		shift_by_color = 0.005
		count = 0
		if chromosome in metadata:
			for arm in metadata[chromosome]:
				for tpl in metadata[chromosome][arm]:
					if count % 2 == 0:
						shift = 0.00095
						count += 1
					else:
						shift = 0
						count += 1
					# ax.plot([x_end + (DIM*0.015)], [y_start + (y_end-y_start) * (md/chromosome_length)], '.', color='black')  # original line
					if tpl[1] < th_log2r_deletion_neg:
						logger.debug("seg.start={} and seg.end={}".format(str(tpl[0]), str(tpl[2]) ))
						logger.debug(("x_end", "DIM", "y_start", "y_end", "chromosome_length", "chromosome"))
						logger.debug("{}".format(' '.join([ str(x) for x in [x_end, DIM, y_start, y_end, chromosome_length, chromosome]])) )

						ax.plot(
							[( (x_end+shift_by_color+shift) + (DIM * 0.015)), ( (x_end+shift_by_color+shift) + (DIM * 0.015)) ],
							[y_start + (y_end-y_start) * (tpl[0]/chromosome_length), y_start + (y_end-y_start) * (tpl[2]/chromosome_length)],
							color='green'
						)
					elif tpl[1] > th_log2r_deletion_pos:
						ax.plot(
							[( (x_end-shift_by_color-shift) + (DIM * 0.015)), ( (x_end-shift_by_color-shift) + (DIM * 0.015)) ],
							[y_start + (y_end - y_start) * (tpl[0] / chromosome_length), y_start + (y_end - y_start) * (tpl[2] / chromosome_length)],
							color='red'
						)
					else:
						ax.plot(
							[(x_end + (DIM * 0.015)), (x_end + (DIM * 0.015))],
							[y_start + (y_end - y_start) * (tpl[0] / chromosome_length), y_start + (y_end - y_start) * (tpl[2] / chromosome_length)],
							color='black'
						)

		
		# center_x : position of the chromosome name from the vertical center of the ideogram (+/- left or right)
		# ((DIM * 0.07) : offset the chromosome name from the ideogram
		# ax.text(center_x, y_end - (DIM * 0.07), chromosome)   # original line
		ax.text(center_x-0.018, y_end - (DIM * 0.07), chromosome)
		
		ax.annotate('threshold\nlog2r >= {} (red)\nlog2r <= {} (green)\nblack in-btw'.format(str(th_log2r_deletion_pos), str(th_log2r_deletion_neg)),
		            xy=(0.8, 1), xycoords='data',
		            xytext=(0, 0), textcoords='offset points',
		            horizontalalignment='right', verticalalignment='top')

	if part == 1:
		plot_chromosome('1', 1)
		plot_chromosome('2', 2)
		plot_chromosome('3', 3)
		plot_chromosome('4', 4)
		plot_chromosome('5', 5)
		plot_chromosome('6', 6)
		plot_chromosome('7', 7)
		plot_chromosome('8', 8)
		plot_chromosome('9', 9)
		plot_chromosome('10', 10)
		plot_chromosome('11', 11)
		plot_chromosome('12', 12)
	elif part == 2:
		plot_chromosome('13', 1)
		plot_chromosome('14', 2)
		plot_chromosome('15', 3)
		plot_chromosome('16', 4)
		plot_chromosome('17', 5)
		plot_chromosome('18', 6)
		plot_chromosome('19', 7)
		plot_chromosome('20', 8)
		plot_chromosome('21', 9)
		plot_chromosome('22', 10)
		plot_chromosome('X', 11)
		plot_chromosome('Y', 12)
	else:
		raise Exception('plot argument should be either "1" or "2"')

	plt.axis('off')
	# plt.show()
	plot_filename = "{}_karyoplot_{}.png".format(samplename, part)
	plt.savefig(plot_filename, dpi=300)
	logger.info("plot karyo {}".format(plot_filename))




