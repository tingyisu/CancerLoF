'''
Prepare files needed to run FoldX PSSM (IR mutations) and BuildModel (non-IR mutations)
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

# get templates that need to be modelled with MODELLER
# these templates are either templates for disrupted interactions (in edgetic mutations) or from a specifically selected interaction 
# on a protein with a non-edgetic mutation

# do for both IR mutations found using interfacial residues between structural templates & between MODELLER-constructed homology models
# combine & take union of FoldX PSSM mutations and FoldX BuildModel mutations

import pickle
import os.path as osp
import os
from simple_tools import pickle_load, pickle_dump, get_mutation_info, write_mutation_info_to_file

class ModellerTemplates:
	def __init__(self, data_dir, mutation_files_list):
		self.data_dir = data_dir
		self.mutation_files_list = mutation_files_list
		self.blast_best_alignments =  pickle_load(osp.join(self.data_dir, 'interactome', 'all_blast_best_alignments.pickle'))
		self.all_ppi_templates = pickle_load(osp.join(self.data_dir, 'interactome', 'all_combined_structural_templates_modeled.pickle'))

		self.foldx_pssm_mutations = [] # all mutations on which to run FoldX PSSM
		self.foldx_buildmodel_mutations = [] # all mutations on which to run FoldX BuildModel
		self.pdb_structures = []

		# for mutation positions that are over 10000 (cannot be accomodated by .pdb files, so need to take %10000)
		self.pssm_over_9999 = [] # list of tuples (foldx pssm mutation, original mutation position prior to taking %)
		self.buildmodel_over_9999 = [] # list of tuples (foldx buildmodel mutation, original mutation position prior to taking %)

		self.foldx_pssm_mutations_pickle_file = osp.join(self.data_dir, 'edgotypes', 'all_foldx_pssm_mutations.pickle')
		self.foldx_buildmodel_mutations_pickle_file = osp.join(self.data_dir, 'edgotypes', 'all_foldx_buildmodel_mutations.pickle')


	# takes list of possible templates for uniprot protein with non-ir mutation
	# returns best template with the lowest evalue (if more than one, picks the first one)
	def select_non_ir_template(self, non_ir_templates):
		evalue_list = []
		for (uniprot_protein, pdb_chain) in non_ir_templates:
			alignment = self.blast_best_alignments[(uniprot_protein, pdb_chain)]
			evalue = float(alignment[8])
			evalue_list.append(evalue)
		min_evalue = min(evalue_list)
		min_index = evalue_list.index(min_evalue)
		# print(evalue_list, min_evalue, min_index, non_ir_templates)
		return min_index, non_ir_templates[min_index][1] # return pdb_chain (sanity check purposes)

	# mutation is ir if it maps onto an interfacial residue of at least one interaction
	def get_foldx_mutations(self, header_dict, mutation_info):
		# print(header_dict)
		# print('-----Compiling templates to be modelled using MODELLERs-----')
		mutation_info[0].append('modeller_foldx_mutations')
		header_dict['modeller_foldx_mutations'] = len(header_dict)
		for missense in mutation_info[1:]:
			# print(len(missense), len(header_dict))
			over_9999 = False
			uniprot_protein = missense[header_dict['uniprot_protein']]
			# if uniprot_protein == 'Q8WZ42':
			# 	print('Yes! Over 9999')
			ir_mutation = missense[header_dict['ir_mutation']]
			protein_partners = missense[header_dict['protein_partners']].split(',')
			on_ir = missense[header_dict['on_interfacial_res']].split(',')
			ref_res = missense[header_dict['ref_res']]
			alt_res = missense[header_dict['alt_res']]
			prot_res_pos_in_uniprot_protein = missense[header_dict['prot_res_pos_in_uniprot_protein']]
			# get foldx_mutation, take % 10000 if mutation position > 9999
			new_prot_res_pos_in_uniprot_protein = -1
			if int(prot_res_pos_in_uniprot_protein) > 9999:
				new_prot_res_pos_in_uniprot_protein = str(int(prot_res_pos_in_uniprot_protein) % 10000)
				over_9999 = True
			else:
				new_prot_res_pos_in_uniprot_protein = prot_res_pos_in_uniprot_protein

			if ir_mutation == '1': # is an IR mutation
				foldx_pssm_mutations_list = [] # for each mutation (line in mutation file)
				for j in range(len(protein_partners)):
					
					if on_ir[j] == '1':
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
						# .pdb formats cannot handle cases where prot_res_pos_in_uniprot_protein > 9999, so do the following
						foldx_mutation = ref_res + chain_with_mutation + new_prot_res_pos_in_uniprot_protein + alt_res 
						foldx_pssm_mutation = '-'.join([dimer, foldx_mutation])
						if over_9999:
							self.pssm_over_9999.append((foldx_pssm_mutation, prot_res_pos_in_uniprot_protein))
						foldx_pssm_mutations_list.append(foldx_pssm_mutation)
						self.foldx_pssm_mutations.append(foldx_pssm_mutation) # add to list, to pickle later
					
					else:
						foldx_pssm_mutations_list.append('-1')
				
				# add foldx_pssm_mutations_list to self.mutation_info
				missense.append(','.join(foldx_pssm_mutations_list))

			else: # non-ir
				# can calculate folding DDG on dimer, so store partner protein and chain
				non_ir_templates = []
				chain_with_mutation_list = [] # store chain with mutation
				dimers = []
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

					# add non_ir templates
					if (uniprot_protein, pdb_chain) not in non_ir_templates: # pick first occurrence if more than one interaction has the same chain
						non_ir_templates.append((uniprot_protein, pdb_chain))
						chain_with_mutation_list.append(chain_with_mutation)
						dimers.append(dimer)

				# get unique templates
				best_dimer, best_chain_with_mutation = '', ''
				if len(non_ir_templates) == 1:
					best_dimer = dimers[0]
					best_chain_with_mutation = chain_with_mutation_list[0]
				else:
					# find the template with the smallest evalue alignment (if more than one, takes first occurrence)
					min_index, _ = self.select_non_ir_template(non_ir_templates)
					best_dimer = dimers[min_index]
					best_chain_with_mutation = chain_with_mutation_list[min_index]
				
				foldx_mutation = ref_res + best_chain_with_mutation + new_prot_res_pos_in_uniprot_protein + alt_res 
				foldx_buildmodel_mutation = '-'.join([best_dimer, foldx_mutation])
				if over_9999:
					self.buildmodel_over_9999.append((foldx_buildmodel_mutation, prot_res_pos_in_uniprot_protein))
				# add foldx_build_model_mutation to self.mutation_info
				missense.append(foldx_buildmodel_mutation)
				self.foldx_buildmodel_mutations.append(foldx_buildmodel_mutation) # add to list, to pickle later

		# sanity check
		# print('Curr number of FoldX PSSM and BuildModel mutations:', len(set(self.foldx_pssm_mutations)), len(set(self.foldx_buildmodel_mutations)))
		
	def get_all_templates_mutations(self):
		for ir_mutations_file in self.mutation_files_list:
			print('Getting foldx mutations for:', ir_mutations_file)

			# get necessary info
			interactome = ir_mutations_file.split('_')[0]
			ir_mutations_file_path = osp.join(self.data_dir, 'mutations_final', ir_mutations_file)
			fname = osp.basename(ir_mutations_file_path)

			header_dict, mutation_info = get_mutation_info(ir_mutations_file_path)
			if len(header_dict) != len(mutation_info[1]): # sanity check
				print('Header & entries have different lengths!', fname, len(header_dict), len(mutation_info[1]))

			# file to write MODELLER & FoldX templates/mutations to 
			ir_modeller_foldx_mutations_file = osp.join(self.data_dir, 'edgotypes', ir_mutations_file.split('.')[0] + '_with_modeller_foldx_mutations.tsv')

			# updates self.foldx_pssm_mutations and self.foldx_buildmodel_mutations for each ir_mutations_file
			self.get_foldx_mutations(header_dict, mutation_info)
			# write mutation info with included modeller foldx mutations to file
			write_mutation_info_to_file(mutation_info, ir_modeller_foldx_mutations_file)


		# take only unique entries in each
		self.foldx_pssm_mutations = list(set(self.foldx_pssm_mutations))
		self.foldx_buildmodel_mutations = list(set(self.foldx_buildmodel_mutations))
		self.pssm_over_9999 = list(set(self.pssm_over_9999))
		self.buildmodel_over_9999 = list(set(self.buildmodel_over_9999))
		print('Total number of FoldX PSSM mutations:', len(self.foldx_pssm_mutations))
		print('Total number of foldx_buildmodel_mutations:', len(self.foldx_buildmodel_mutations))
		print('PSSM mutations with positions over 9999:', self.pssm_over_9999)
		print('BuildModel mutations with positions over 9999:', self.buildmodel_over_9999)

		# pickle into corresponding files
		pickle_dump(self.foldx_pssm_mutations, self.foldx_pssm_mutations_pickle_file)
		pickle_dump(self.foldx_buildmodel_mutations, self.foldx_buildmodel_mutations_pickle_file)


def main():
	script_dir = osp.dirname(__file__)
	data_dir = osp.join(script_dir, '..', 'data', 'processed')

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
					mutation_files_list.append('_'.join([interactome, 'mapped', data, 'missense', 'ir', ir_type + '.tsv']))
				else:
					for cosmic_data in cosmic_datasets:
						mutation_files_list.append('_'.join([interactome, 'mapped', data, cosmic_data, 'missense', 'ir', ir_type + '.tsv']))


	templates = ModellerTemplates(data_dir, mutation_files_list)
	templates.get_all_templates_mutations()

if __name__=='__main__':
	main()