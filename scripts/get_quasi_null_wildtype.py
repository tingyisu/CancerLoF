'''
Find quasi-null and quasi-wildtype mutations using folding DDG
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os.path as osp
import os
from simple_tools import pickle_load, pickle_dump, check_create_dir, get_mutation_info, write_mutation_info_to_file, check_exist
import argparse

class QuasiNullWildtype:
	def __init__(self, names, script_dir):
		self.names = names
		self.data_dir = osp.join(script_dir, 'files')
		self.folding_ddg_dict = pickle_load(osp.join(self.data_dir, 'foldx_buildmodel_folding_ddg_dict.pickle')) # key = FoldX BuildModel mutation, value = folding DDG 

	def get_quasi_null_wildtype_all(self):
		for name in self.names:
			print('Finding quasi null/wildtype mutations for:', name)
			# read in edgotypes file
			header_dict, mutation_info = get_mutation_info(osp.join(self.data_dir, name + '_mutation_edgotypes.tsv'))

			# get quasi-null/wildtype mutations and write to file
			mutation_info_edgotype_quasi_null_wildtype = osp.join(self.data_dir, name + '_mutation_edgotypes_quasi_null_wildtype.tsv')
			self.get_quasi_null_wildtype(name, mutation_info, header_dict, mutation_info_edgotype_quasi_null_wildtype)
		
	def get_quasi_null_wildtype(self, name, mutation_info, header_dict, mutation_info_edgotype_quasi_null_wildtype):
		mutation_info[0].append('foldx_buildmodel_folding_ddg')
		header_dict['foldx_buildmodel_folding_ddg'] = len(header_dict)
		mutation_info[0].append('quasi_null_wildtype')
		header_dict['quasi_null_wildtype'] = len(header_dict)

		for missense in mutation_info[1:]:
			# print(header_dict)
			edgotype = missense[header_dict['edgotype']] # switch to edgotype?
			if edgotype == 'non-edgetic':
				foldx_buildmodel_mutation = missense[header_dict['foldx_buildmodel_mutations']]
				# print('Looking at:', foldx_buildmodel_mutation)
				folding_ddg = self.folding_ddg_dict[foldx_buildmodel_mutation]
				missense.append(folding_ddg)

				# quasi-wildtype (exposed residue, mutation doesn't make protein fall apart)
				# quasi-null (buried residue, mutation makes protein fall apart)
				# add quasi-wildtype/quasi-null to missense
				if folding_ddg >= 2.0: # mutation de-stabilizes protein
					missense.append('quasi-null')
				else: # folding_ddg < 2.0 (mutation isn't enough to de-stabilize protein)
					missense.append('quasi-wildtype')
			else: # edgetic mutation
				missense.append('NA')
				missense.append('NA')
		write_mutation_info_to_file(mutation_info, mutation_info_edgotype_quasi_null_wildtype)


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--script_dir')  #/home/tingyisu/projects/def-yxia/tingyisu/foldx
	args = parser.parse_args()

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


	print('-----Finding quasi null, quasi wildtype based on folding DDG-----')
	q = QuasiNullWildtype(mutation_files_list, args.script_dir)
	q.get_quasi_null_wildtype_all()

	

if __name__=='__main__':
	main()


