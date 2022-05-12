
import logging
from makeFastq import get_runtime_function
import scipy.cluster.hierarchy as hcluster
from collections import Counter, defaultdict


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)s %(asctime)-15s %(module)s %(lineno)d\t %(message)s')
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)


class ReadClusteringManagement:
	"""
	Manage the clusterings of the reads for each Key_Contig
	"""
	def __init__(self,
	             dico_pairs_of_start, dico_reads_for_pairs_of_start,
	             selected_clustering_type="hclust",
	             cluster_left_pos_only=True,
	             min_tpls=3,
	             min_dist_leftmost_to_rightmost_reads=150,
	             min_read_count_per_cluster_to_Pass=3,
	             threshold_distance_for_clustering=1000
	             ):

		self.dico_pairs_of_start = dico_pairs_of_start  # this is a list of tuples with two integer values associated to the key_contig; e.g.: [(1000,200),(800,500),(1050,250)]
		self.dico_reads_for_pairs_of_start = dico_reads_for_pairs_of_start  # read names mapping the list_of_pairs_of_tuples_of_start_positions ; need to keep track to later microassmble
		self.selected_clustering_type = selected_clustering_type
		self.dico_clusters_by_pos = {}
		self.min_tpls = min_tpls
		self.min_read_count_per_cluster_to_Pass = min_read_count_per_cluster_to_Pass
		self.min_dist_leftmost_to_rightmost_reads = min_dist_leftmost_to_rightmost_reads
		self.cluster_left_pos_only = cluster_left_pos_only
		self.threshold_distance_for_clustering = threshold_distance_for_clustering
		if self.cluster_left_pos_only and self.dico_pairs_of_start is not None and isinstance(self.dico_pairs_of_start, dict):
			self.reformat_dico_pairs_of_start()

	def filter_clusters(self):
		pass

	def filter_keys_with_number_pairs_less_than(self):
		pass

	def cluster_start_pos(self):
		pass

	def print_dico_clusters_by_pos(self, n=5):

		for k in self.dico_clusters_by_pos.keys():
			my_counter = 0
			for v in self.dico_clusters_by_pos[k]:
				if my_counter > n:
					break
				logger.debug("{} :::::: {}".format(k, str(v)))
				my_counter += 1

	def print_dico_pos(self, n=5):

		for k in self.dico_pairs_of_start.keys():
			my_counter = 0
			for v in self.dico_pairs_of_start[k]:
				if my_counter > n:
					break
				logger.debug("{} :::::: {}".format(k, str(v)))
				my_counter += 1

	def print_dico_reads_for_pairs_of_start(self, n=5):
		for k in self.dico_reads_for_pairs_of_start.keys():
			my_counter = 0
			for v in self.dico_reads_for_pairs_of_start[k]:  # assert v is a tuple!
				if my_counter > n:
					break
				logger.debug("{} :::::: {}".format(k, str(v)))
				my_counter += 1

	def filter_keys_with_number_tuples_pos_less_than_th(self):
		"""
		Filter the keys where the number of tuple of positions is less than threshold
		:param min_tpls is the the number minimum of tuples of positions we want to be clustered ;
		:rtype integer
		:return: None
		:rtype: None
		"""
		dpos = self.dico_pairs_of_start.copy()
		dreads = self.dico_reads_for_pairs_of_start.copy()
		for k, v in self.dico_pairs_of_start.items():
			if len(v) < self.min_tpls:
				dpos.pop(k)
				dreads.pop(k)
		self.dico_pairs_of_start = dpos
		self.dico_reads_for_pairs_of_start = dreads

	def reformat_dico_pairs_of_start(self):
		dpos = defaultdict(list)
		for k, v in self.dico_pairs_of_start.items():
			ltpls_pos = []
			for tpl in v:
				ltpls_pos.append((tpl[0],))
			dpos[k] = ltpls_pos
		self.dico_pairs_of_start = dpos

	@staticmethod
	def get_clusters(LIDXCLUS, lidxrcgeth):
		"""
		Subset original cluster list using another list;
		Here we use the list of indexes of the list of clusters that got filtered by minimum readcount (get_list_indexes_with_read_count_ge_th(c, 10))

		:param LIDXCLUS:
		:type LIDXCLUS:
		:param lidxrcgeth: list of indexes
		:type lidxrcgeth: list
		:return:
		:rtype:
		"""
		# LIDXCLUS stands for List Index of the CLUSters
		for cluster_idx in lidxrcgeth:
			yield [i for i, v in enumerate(LIDXCLUS) if v == cluster_idx]

	def is_read_positions_distance_distribution_ok(self, cluster_values, min_distance_expected=150):

		lpos1 = []
		lpos2 = []
		sum_lpos2 = 0  # We need to take into account the fact taht unmapped read will have pos2 == 0
		if not self.cluster_left_pos_only:
			for tpl in cluster_values:
				lpos1.append(tpl[0])
				lpos2.append(tpl[1])
				sum_lpos2 += tpl[1]
			dist1 = max(lpos1) - min(lpos1)
			dist2 = max(lpos2) - min(lpos2)
			if sum_lpos2 != 0:
				if dist1 >= min_distance_expected and dist2 >= min_distance_expected:
					return True
			elif dist1 >= min_distance_expected:
				return True
			return False
		else:
			for tpl in cluster_values:
				lpos1.append(tpl[0])
			dist1 = max(lpos1) - min(lpos1)
			if dist1 >= min_distance_expected:
				return True
			return False

	# #keep indexes where the number of reads clustered together is more than X
	@staticmethod
	def get_list_indexes_with_read_count_ge_th(counter, th):
		# print("minimum number of reads (including multiple starts (so not only the unique start)) in cluster: " + str(th))
		c = counter
		lidx = []
		for k, v in zip(c.keys(), c.values()):
			logger.debug(str(k) + "  and v=  " + str(v))
			if v >= th:
				lidx.append(k)
		return lidx

	@get_runtime_function
	def proceed_clustering(self, criterion="distance"):
		"""
		Run the process of clustering and return the appropriate dictionaries
		:param criterion:
		:type criterion:
		:return:
		:rtype:
		"""
		ddr_clusters = {}
		ddr_counter_clusters_per_key = {}
		for k, list_of_tuple_positions in self.dico_pairs_of_start.items():
			clustering = hcluster.fclusterdata(list_of_tuple_positions, self.threshold_distance_for_clustering, criterion=criterion)
			logger.debug(str(k) + " <----> " + str(clustering.tolist()))
			ddr_clusters[k] = clustering.tolist()
			ddr_counter_clusters_per_key[k] = Counter(clustering.tolist())
			logger.debug(ddr_counter_clusters_per_key[k])
			most_commons = ddr_counter_clusters_per_key[k].most_common(10)
			logger.debug("most_commons" + str(most_commons))
			# logger.debug(self.get_list_indexes_with_read_count_ge_th(ddr_counter_clusters_per_key[k], 5))
		self.dico_clusters_by_pos = ddr_clusters
		logger.debug(self.dico_clusters_by_pos)

		d_res = defaultdict(dict)
		drap_res = defaultdict(dict)

		for processed_key in self.dico_pairs_of_start.keys():
			LTI = self.dico_pairs_of_start[processed_key]  # LTI == List Tuple IndexKey
			logger.debug(str(len(LTI)) + "-->> " + str(LTI))
			LTALNREC = self.dico_reads_for_pairs_of_start[processed_key]
			logger.debug(str(len(LTALNREC)) + "-->> " + str(LTALNREC))
			LIDXCLUS = self.dico_clusters_by_pos[processed_key]  # LIDXCLUS == List InDeX CLUSter
			logger.debug(str(len(LIDXCLUS)) + "-->> " + str(LIDXCLUS))
			c = Counter(LIDXCLUS)

			if logger.level == 10:
				logger.debug(str(c.most_common(1)))
				# #[(3, 11)]
				logger.debug("cluster with the most number of items: " + str(c.most_common(1)[0][0]))
				logger.debug("number of items in that cluster: " + str(c.most_common(1)[0][1]))
				# # 3
				idxclust = c.most_common(1)[0][0]
				logger.debug(c.items())
				# #dict_items([(3, 11), (2, 2), (1, 3), (4, 1), (5, 1)])
				logger.debug(c.keys())
				# #dict_keys([3, 2, 1, 4, 5])
				logger.debug(c.values())
				# #dict_values([11, 2, 3, 1, 1])

			# lidxrcgeth == List Index Read Count Greater than THreshold
			lidxrcgeth = self.get_list_indexes_with_read_count_ge_th(c, self.min_read_count_per_cluster_to_Pass)
			# [i for i, v in enumerate(LIDXCLUS) if v in lidxrcgeth]
			# #[0, 1, 2, 4, 6, 8, 9, 10, 12, 14, 16]

			# we replaced the data in the ddr_filtered with the clustered data
			# but this is wrong because some keys might be skipped because the number of reads within clusters is low and therefore the number of keys decreases
			# So : create a new dict for keeping the results
			dclust = defaultdict(list)
			dclust_alnrecs = defaultdict(list)

			for ci, e in enumerate(self.get_clusters(LIDXCLUS, lidxrcgeth)):
				for idx in e:
					dclust[ci].append(LTI[idx])
					dclust_alnrecs[ci].append(LTALNREC[idx])
				if not self.is_read_positions_distance_distribution_ok(dclust[ci], min_distance_expected=self.min_dist_leftmost_to_rightmost_reads):
					if self.min_dist_leftmost_to_rightmost_reads != 0:  # 0 disable the check; implement better to disable #TODO
						logger.warning("Cluster of reads DID NOT PASS FILTER of Distance Distribution ... could be a pile of reads only :-(  ..." + str(processed_key) + ":::" + str(dclust[ci]))
						continue
				self.dico_pairs_of_start[processed_key] = dclust
				if processed_key in d_res.keys():
					d_res[processed_key][ci] = dclust[ci]
					drap_res[processed_key][ci] = dclust_alnrecs[ci]
				else:
					d_res[processed_key] = {}
					d_res[processed_key][ci] = dclust[ci]
					drap_res[processed_key] = {}
					drap_res[processed_key][ci] = dclust_alnrecs[ci]

			if logger.level == 10:
				logger.debug("##@@@@@@@@@@@@@@@@@@@@@@@##")
				logger.debug(d_res)
				logger.debug("##@@@@@@@@@@@@@@@@@@@@@@@##")
			logger.debug("End process KEY: " + processed_key)
		logger.debug("*"*15 + "END LOOP to process clusters " + "*"*15)
		logger.debug("@" * 50)
		if logger.level == 10:
			for k, dv in d_res.items():  # d_res.items() is a dictionary of dictionary
				for cluster_id, v in dv.items():
					logger.debug(k + " ==> cluster #" + str(cluster_id) + " ==> len:" + str(len(v)) + " list: " + str(v))

			for k, dv in drap_res.items():  # d_res.items() is a dictionary of dictionary
				for cluster_id, v in dv.items():
					logger.debug(k + " ==> cluster #" + str(cluster_id) + " ==> len:" + str(len(v)) + " list: " + str(v))
					for talns in v:
						logger.debug(k + " ==> cluster #" + str(cluster_id) + " ==> len:" + str(len(v)) + " tuple_pair_aln: " + str(talns[0]) + "%%%%%" + str(talns[1]))

		return dict(d_res), dict(drap_res)

