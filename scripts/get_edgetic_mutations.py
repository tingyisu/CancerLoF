'''
Finding edgetic mutations (using binding DDG)
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os.path as osp
import os
from simple_tools import pickle_load, check_create_dir, get_mutation_info, write_mutation_info_to_file, pickle_dump, check_exist
import argparse

class BindingDDG:
	def __init__(self, names, script_dir, scratch_dir):
		self.names = names
		self.data_dir = osp.join(script_dir, 'files')
		self.blast_best_alignments =  pickle_load(osp.join(self.data_dir, 'all_blast_best_alignments.pickle'))
		self.all_ppi_templates = pickle_load(osp.join(script_dir, 'interactome', 'all_combined_structural_templates_modeled.pickle')) # key = (uniprot1, uniprot2), value = (pdb_chain1, pdb_chain2)
		self.all_foldx_buildmodel_mutations = pickle_load(osp.join(self.data_dir, 'all_foldx_buildmodel_mutations.pickle')) # update with additional FoldX BuildModel mutations
		self.foldx_pssm_dict = pickle_load(osp.join(self.data_dir, 'foldx_pssm_binding_ddg_dict.pickle'))
		# self.foldx_buildmodel_dir = osp.join(scratch_dir, 'foldx_buildmodel_all')
		self.additional_foldx_buildmodel_mutations = []
		self.additional_foldx_buildmodel_mutations_pickle_file = osp.join(self.data_dir, 'additional_foldx_buildmodel_mutations.pickle') # additional 
		self.all_foldx_buildmodel_mutations_combined_pickle_file = osp.join(self.data_dir, 'all_foldx_buildmodel_mutations_combined.pickle') # plus additional (self.all_foldx_buildmodel_mutations + self.additional_foldx_buildmodel_mutations)

	# takes list of possible templates for uniprot protein with non-edgetic mutation
	# returns best template with the lowest evalue (if more than one, picks the first one)
	def select_non_edgetic_template(self, non_ir_templates):
		evalue_list = []
		for (uniprot_protein, pdb_chain) in non_ir_templates:
			alignment = self.blast_best_alignments[(uniprot_protein, pdb_chain)]
			evalue = float(alignment[8])
			evalue_list.append(evalue)
		min_evalue = min(evalue_list)
		min_index = evalue_list.index(min_evalue)
		# print(evalue_list, min_evalue, min_index, non_ir_templates)
		return min_index, non_ir_templates[min_index][1] # return pdb_chain

	def get_edgotype_all(self):
		for name in self.names:
			header_dict, mutation_info = get_mutation_info(osp.join(self.data_dir, name + '_with_modeller_foldx_mutations.tsv'))
			mutation_info_edgotype_file = osp.join(self.data_dir, name + '_mutation_edgotypes.tsv')
			self.get_binding_ddg_edgotype(name, mutation_info, header_dict, mutation_info_edgotype_file)
		print('------Additional items to run-----')
		# combine self.all_foldx_buildmodel_mutations and self.additional_foldx_buildmodel_mutations and pickle to file
		all_combined_foldx_buildmodel_mutations = list(set(self.additional_foldx_buildmodel_mutations) | set(self.all_foldx_buildmodel_mutations))
		
		print('Number of additional FoldX BuildModel mutations:', len(set(self.additional_foldx_buildmodel_mutations)))
		print('Total number of combined FoldX BuildModel mutations:', len(all_combined_foldx_buildmodel_mutations))
		# pickle both additional and combined 
		pickle_dump(list(set(self.additional_foldx_buildmodel_mutations)), self.additional_foldx_buildmodel_mutations_pickle_file)
		pickle_dump(all_combined_foldx_buildmodel_mutations, self.all_foldx_buildmodel_mutations_combined_pickle_file)

	def get_binding_ddg_edgotype(self, name, mutation_info, header_dict, mutation_info_edgotype_file):
		edgetic, non_edgetic = 0, 0
		ir_non_edgetic = 0
		mutation_info[0].append('foldx_pssm_binding_ddg')
		header_dict['foldx_pssm_binding_ddg'] = len(header_dict)
		mutation_info[0].append('perturbations')
		header_dict['perturbations'] = len(header_dict)
		mutation_info[0].append('edgotype')
		header_dict['edgotype'] = len(header_dict)
		mutation_info[0].append('foldx_buildmodel_mutations')
		header_dict['foldx_buildmodel_mutations'] = len(header_dict)
		for missense in mutation_info[1:]:
			ir_mutation = missense[header_dict['ir_mutation']]
			on_interfacial_res = missense[header_dict['on_interfacial_res']]
			if ir_mutation == '0': # not on interfacial residue
				modeller_foldx_mutations = missense[header_dict['modeller_foldx_mutations']]
				missense.append('NA')
				missense.append(on_interfacial_res)
				missense.append('non-edgetic')
				missense.append(modeller_foldx_mutations)
				non_edgetic += 1
			else: # on interfacial residue
				modeller_foldx_mutations = missense[header_dict['modeller_foldx_mutations']].split(',')
				binding_ddgs = []
				perturbations = []
				foldx_mutation = ''

				for modeller_foldx_mutation in modeller_foldx_mutations:
					if modeller_foldx_mutation != '-1':
						# get binding DDG and perturbation
						binding_ddg = self.foldx_pssm_dict[modeller_foldx_mutation]
						binding_ddgs.append(str(binding_ddg))
						if binding_ddg > 0.8: # updated according to FoldX's estimated DDG error
							perturbations.append('1')
						else:
							perturbations.append('0')
						# get FoldX mutation
						if foldx_mutation == '':
							foldx_mutation = modeller_foldx_mutation.split('-')[-1]
					else:
						binding_ddgs.append('NA')
						perturbations.append('0')

				missense.append(','.join(binding_ddgs))
				missense.append(','.join(perturbations))
				if '1' in perturbations:
					missense.append('edgetic')
					missense.append('NA')
					edgetic += 1
				else:
					missense.append('non-edgetic')
					non_edgetic += 1
					ir_non_edgetic += 1
					non_edgetic_templates = []
					chain_with_mutation_list = [] # store chain with mutation
					dimers = []
					uniprot_protein = missense[header_dict['uniprot_protein']]
					protein_partners = missense[header_dict['protein_partners']].split(',')

					# FIX ACCORDING TO WHAT'S IN get_foldx_mutations.py
	
					for j in range(len(protein_partners)):

						pdb_chain, pdb_chain2 = '', ''
						chain_with_mutation = '' # mutation can either be on chain A or chain B
						dimer = ''
						if (uniprot_protein, protein_partners[j]) in self.all_ppi_templates:
							pdb_chain, pdb_chain2 = self.all_ppi_templates[(uniprot_protein, protein_partners[j])]
							chain_with_mutation = 'A'
							dimer = '_'.join([uniprot_protein, pdb_chain, protein_partners[j], pdb_chain2])
						elif (protein_partners[j], uniprot_protein) in self.all_ppi_templates:
							pdb_chain2, pdb_chain = self.all_ppi_templates[(protein_partners[j], uniprot_protein)]
							chain_with_mutation = 'B'
							dimer = '_'.join([protein_partners[j], pdb_chain2, uniprot_protein, pdb_chain])
						else:
							print('ERROR!!! interaction pair is not in self.selected_templates_dict!', uniprot_protein, protein_partners[j])

						# add non_edgetic templates
						if (uniprot_protein, pdb_chain) not in non_edgetic_templates: # pick first occurrence if more than one interaction has the same chain
							non_edgetic_templates.append((uniprot_protein, pdb_chain))
							chain_with_mutation_list.append(chain_with_mutation)
							dimers.append(dimer)

					# get unique templates
					best_dimer, best_chain_with_mutation = '', ''
					if len(non_edgetic_templates) == 1:
						best_dimer = dimers[0]
						best_chain_with_mutation = chain_with_mutation_list[0]
					else:
						# find the template with the smallest evalue alignment (if more than one, takes first occurrence)
						min_index, _ = self.select_non_edgetic_template(non_edgetic_templates)
						best_dimer = dimers[min_index]
						best_chain_with_mutation = chain_with_mutation_list[min_index]

					foldx_mutation_updated = foldx_mutation[0] + best_chain_with_mutation + foldx_mutation[2:] # update chain (can either be 'A' or 'B')
					foldx_buildmodel_mutation = '-'.join([best_dimer, foldx_mutation_updated])
					# add foldx_build_model_mutation to self.mutation_info
					missense.append(foldx_buildmodel_mutation)
					if foldx_buildmodel_mutation not in self.all_foldx_buildmodel_mutations:
						self.additional_foldx_buildmodel_mutations.append(foldx_buildmodel_mutation) # add to list, to pickle later
		
		print('-----' + name + '-----')
		print('Number of edgetic mutations:', edgetic)
		print('Number of non-edgetic mutations:', non_edgetic)
		print('Number of IR mutations that are non-edgetic:', ir_non_edgetic)
		write_mutation_info_to_file(mutation_info, mutation_info_edgotype_file)


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--script_dir')  #/home/username/projects/def-*/username/foldx
	parser.add_argument('-c', '--scratch_dir') #/home/username/scratch
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


	print('-----Finding edgotypes (edgetic/non-edgetic) using FoldX binding DDG-----')
	b = BindingDDG(mutation_files_list, args.script_dir, args.scratch_dir)
	b.get_edgotype_all()

	

if __name__=='__main__':
	main()


