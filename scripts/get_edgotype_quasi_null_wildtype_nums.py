'''
Print number of edgetic, quasi-null and quasi-wildtype mutations
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

# prints number of edgetic, quasi-wildtype, quasi-null mutations in *_mutation_edgotypes_quasi_null_wildtype.tsv
import os.path as osp
import os
import math
import numpy as np
from simple_tools import get_mutation_info, pickle_dump

class Nums:
	def __init__(self, data_dir, names):
		self.data_dir = data_dir
		self.names = names

	def get_edgotype_quasi_null_wildtype_nums_all(self):
		for name in self.names:
			# print('Finding number of edgetic, quasi-wildtype, and quasi-null mutations for:', name)
			header_dict, mutation_info = get_mutation_info(osp.join(self.data_dir, name + '_mutation_edgotypes_quasi_null_wildtype.tsv'))
			self.get_edgotype_quasi_null_wildtype_nums(name, mutation_info, header_dict)

	def get_avg_blosum_score(self, blosum_scores_list, num):
		if num == 0 and blosum_scores_list == []:
			return 0
		else:
			return sum(blosum_scores_list)/num

	def get_blosum_score_error_bars(self, blosum_scores_list, sqrt_total_num):
		if blosum_scores_list == []:
			return 0
		else:
			return np.std(blosum_scores_list)/sqrt_total_num

	def get_edgotype_quasi_null_wildtype_nums(self, name, mutation_info, header_dict):
		print('-----' + name + '-----')
		num_edgetic, num_quasi_wildtype, num_quasi_null = 0, 0, 0
		total_num_mutations = len(mutation_info) - 1
		edgetic_blosum_scores, quasi_null_blosum_scores, quasi_wildtype_blosum_scores = [], [], []
		edgotypes_blosum_scores_pickle_file = osp.join(self.data_dir, name + '_edgotypes_blosum_scores.pickle')

		for missense in mutation_info[1:]:
			# print(header_dict)
			edgotype = missense[header_dict['edgotype']]
			blosum_score = float(missense[header_dict['blosum62_score']]) # new BLOSUM62 matrix gives scores in float format (instead of int)
			if edgotype == 'non-edgetic':
				quasi_null_wildtype = missense[header_dict['quasi_null_wildtype']]
				if quasi_null_wildtype == 'quasi-wildtype':
					num_quasi_wildtype += 1
					quasi_wildtype_blosum_scores.append(blosum_score)
				else:
					if quasi_null_wildtype != 'quasi-null': # sanity check
						print('Something went wrong with:', missense)
					else:
						num_quasi_null += 1
						quasi_null_blosum_scores.append(blosum_score)
			else: # edgetic mutation
				num_edgetic += 1
				edgetic_blosum_scores.append(blosum_score)
		# avg_blosum_scores = [sum(quasi_wildtype_blosum_scores)/num_quasi_wildtype, sum(edgetic_blosum_scores)/num_edgetic, sum(quasi_null_blosum_scores)/num_quasi_null]
		avg_blosum_scores = [self.get_avg_blosum_score(quasi_wildtype_blosum_scores, num_quasi_wildtype), self.get_avg_blosum_score(edgetic_blosum_scores, num_edgetic), self.get_avg_blosum_score(quasi_null_blosum_scores, num_quasi_null)]
		sqrt_total_num = math.sqrt(total_num_mutations)
		# blosum_scores_error_bars = [np.std(quasi_wildtype_blosum_scores)/sqrt_total_num, np.std(edgetic_blosum_scores)/sqrt_total_num, np.std(quasi_null_blosum_scores)/sqrt_total_num]
		blosum_scores_error_bars = [self.get_blosum_score_error_bars(quasi_wildtype_blosum_scores, sqrt_total_num), self.get_blosum_score_error_bars(edgetic_blosum_scores, sqrt_total_num), self.get_blosum_score_error_bars(quasi_null_blosum_scores, sqrt_total_num)]

		nums_list = [total_num_mutations, num_quasi_wildtype, num_edgetic, num_quasi_null]

		pickle_dump([nums_list, avg_blosum_scores, blosum_scores_error_bars], edgotypes_blosum_scores_pickle_file)
		print('Total number of mutations:', total_num_mutations)
		print('Number of quasi-wildtype mutations:', num_quasi_wildtype, "{:.3f}".format(num_quasi_wildtype/total_num_mutations))
		print('Number of edgetic mutations:', num_edgetic, "{:.3f}".format(num_edgetic/total_num_mutations))
		print('Number of quasi-null mutations:', num_quasi_null, "{:.3f}".format(num_quasi_null/total_num_mutations))
		print('Avg BLOSUM62 scores [quasi-wildtype, edgetic, quasi-null]:', avg_blosum_scores)

def main():
	data_dir = osp.join(osp.dirname(__file__), '..', 'data', 'processed', 'edgotypes')

	interactomes = ['hiunion', 'intact']
	mutation_data = ['dbsnp', 'clinvar', 'cosmic']
	cosmic_datasets = ['genome_screen', 'cgc']
	ir_types = ['structural_templates', 'homology_models']
	mutation_files_list = []


	# ONLY FOR MISSENSE MUTATIONS

	for interactome in interactomes:
		for ir_type in ir_types:
			for data in mutation_data:
				if data != 'cosmic':
					mutation_files_list.append('_'.join([interactome, 'mapped', data, 'missense', 'ir', ir_type]))
				else:
					for cosmic_data in cosmic_datasets:
						mutation_files_list.append('_'.join([interactome, 'mapped', data, cosmic_data, 'missense', 'ir', ir_type]))

	# EXTRA (remove later?) ClinVar mutations causing degenerative diseases
	for interactome in interactomes:
		mutation_files_list.append('_'.join(['degenerative', interactome, 'mapped', 'clinvar', 'missense', 'ir', 'homology_models']))
	# for mutation_file in mutation_files_list:
	# 	print(mutation_file)

	n = Nums(data_dir, mutation_files_list)
	n.get_edgotype_quasi_null_wildtype_nums_all()

if __name__=='__main__':
	main()
