'''
Find mutations that coincide with interfacial residues (IR mutations)
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import pickle
import os.path as osp
import os
from simple_tools import pickle_load, pickle_dump, get_mutation_info, write_mutation_info_to_file, convert_str_of_items_to_list_of_ints

# removes mutations that do not have protein partners or that don't have any leftoever interactions in which the mutation position lies within the structural template
# THIS VERSION MAPS TO THE FOLLOWING: 
# (1) interfacial residues between structural templates
# (2) interfacial residues between MODELLER constructed homology model (of a heterodimer) that has been repaired using FoldX's RepairPDB

class IRMutation:
	def __init__(self, data_dir, final_data_dir, mutation_files_list):
		self.data_dir = data_dir
		self.final_data_dir = final_data_dir
		self.mutation_files_list = mutation_files_list

		interactome_data_dir = osp.join(data_dir, 'interactome')
		self.all_ppi_templates = pickle_load(osp.join(interactome_data_dir, 'all_combined_structural_templates_modeled.pickle')) # key = (uniprot1, uniprot2), value = (pdb_chain1, pdb_chain2)
		self.all_structural_template_interfacial_residues = pickle_load(osp.join(interactome_data_dir, 'all_combined_interfacial_residues_modeled.pickle')) # key = (uniprot1, uniprot2), value = [residue1_str, residue2_str]]
		self.all_homology_model_interfacial_residues = pickle_load(osp.join(interactome_data_dir, 'modeller_heterodimer_interfacial_residues.pickle')) # key = (uniprot1, uniprot2), value = [residue1_str, residue2_str]]
		self.all_blast_best_alignment = pickle_load(osp.join(interactome_data_dir, 'all_blast_best_alignments.pickle')) # key = (protein, pdb_chain), value = [pident, length, midmatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, protein alignment, template_chain alignment]
		self.hiunion_selected_pdb_structural_template = pickle_load(osp.join(interactome_data_dir, 'hiunion_selected_pdb_structural_template.pickle'))
		self.intact_selected_pdb_structural_template = pickle_load(osp.join(interactome_data_dir, 'intact_selected_pdb_structural_template.pickle'))
		
		self.hiunion_uniprot_pairs_updated = []
		self.intact_uniprot_pairs_updated = []

		self.hiunion_uniprot_pairs_updated_pickle_file = osp.join(interactome_data_dir, 'hiunion_uniprot_pairs_updated.pickle')
		self.intact_uniprot_pairs_updated_pickle_file = osp.join(interactome_data_dir, 'intact_uniprot_pairs_updated.pickle')
	

	# reverse interaction pairs to match the pair-ordering in self.all_ppi_templates
	def get_updated_uniprot_pairs(self, selected_structural_template):
		uniprot_pairs_updated = []
		for (p1, p2) in selected_structural_template:
			if (p1, p2) in self.all_ppi_templates:
				uniprot_pairs_updated.append((p1, p2))
			elif (p2, p1) in self.all_ppi_templates:
				uniprot_pairs_updated.append((p2, p1))
			else:
				print('UniProt pair was unable to be modeled?', (p1, p2))
		return uniprot_pairs_updated
	
	# update HI-union and IntAct selected structural templates from all_combined_structural_templates.pickle
	def get_updated_uniprot_pairs_all(self):
		self.hiunion_uniprot_pairs_updated = self.get_updated_uniprot_pairs(self.hiunion_selected_pdb_structural_template)
		self.intact_uniprot_pairs_updated = self.get_updated_uniprot_pairs(self.intact_selected_pdb_structural_template)

		print('Number of modeled HI-union PPIs:', len(self.hiunion_uniprot_pairs_updated))
		print('Number of modeled IntAct PPIs:', len(self.intact_uniprot_pairs_updated))

		pickle_dump(self.hiunion_uniprot_pairs_updated, self.hiunion_uniprot_pairs_updated_pickle_file)
		pickle_dump(self.intact_uniprot_pairs_updated, self.intact_uniprot_pairs_updated_pickle_file)

	# gets interactors and initializes mutation_info_with_protein_partners (are the mutations on proteins that exist in the list of interactions)
	# only keeps mutations with at least one existing protein partner (interaction)
	def get_interactors(self, header_dict, mutation_info, uniprot_pairs):
		interactors = {}
		mutation_info_with_protein_partners = []

		mutation_info_with_protein_partners.append(mutation_info[0])
		for i in range(1, len(mutation_info)):
			uniprot_protein = mutation_info[i][header_dict['uniprot_protein']]
			if uniprot_protein not in interactors:
				# get interacting partners
				partners = []
				for (p1, p2) in uniprot_pairs:
					if uniprot_protein == p1:
						partners.append(p2)
					elif uniprot_protein == p2:
						partners.append(p1)
				interactors[uniprot_protein] = partners
				# add if has interacting partners
				if partners != []:
					mutation_info_with_protein_partners.append(mutation_info[i])
			else:
				# add if has interacting partners
				if interactors[uniprot_protein] != []:
					mutation_info_with_protein_partners.append(mutation_info[i])

		return interactors, mutation_info_with_protein_partners

		# print(self.interactors)
		# print(len(self.final_mutation_info), self.final_mutation_info[20])
		# print(len(self.header_dict), len(self.final_mutation_info[10]))

	def find_ir_mutations(self, header_dict, interactors, mutation_info_with_protein_partners, uniprot_pairs, interfacial_residues, ir_mutations_file):
		
		mutation_info_position_within_blast_alignment = []

		total_num_mutations = 0
		num_removed = 0
		num_kept = 0
		mutation_info_position_within_blast_alignment.append(mutation_info_with_protein_partners[0] + ['protein_partners', 'on_interfacial_res', 'ir_mutation'])
		print('-----' + osp.basename(ir_mutations_file) + '-----')
		# mutation_info_position_within_blast_alignment[0].append('protein_partners')
		# mutation_info_position_within_blast_alignment[0].append('on_interfacial_res') # on interfacial residue
		# mutation_info_position_within_blast_alignment[0].append('ir_mutation') # 0/1 (1 = on interfacial residue, 0 = not on interfacial residue)

		num_ir, num_non_ir = 0, 0
		for i in range(1, len(mutation_info_with_protein_partners)):
			total_num_mutations += 1
			uniprot_protein = mutation_info_with_protein_partners[i][header_dict['uniprot_protein']] # protein that mutation lies on
			mut_pos = int(mutation_info_with_protein_partners[i][header_dict['prot_res_pos_in_uniprot_protein']]) # mutation position in uniprot numbering
			partners = interactors[uniprot_protein]
			# check that mutation is within blast alignment of for each interaction, remove interaction if is not
			# if no interactions left, then remove mutation
			confirmed_partners = [] # keep only interactions where mutation lies within blast alignment
			perturbations = []
			for p in partners:
				residues = []
				pdb, chain1, chain2 = "", "", ""
				# get residues and perturbations
				if (uniprot_protein, p) in uniprot_pairs:
					residues = convert_str_of_items_to_list_of_ints(interfacial_residues[(uniprot_protein, p)][0])
					# print(residues)
					template, partner_template = self.all_ppi_templates[(uniprot_protein, p)]
					pdb, chain1 = template.split('_')
					_, chain2 = partner_template.split('_')
				elif (p, uniprot_protein) in uniprot_pairs:
					residues = convert_str_of_items_to_list_of_ints(interfacial_residues[(p, uniprot_protein)][1])
					# print(residues)
					partner_template, template = self.all_ppi_templates[(p, uniprot_protein)]
					pdb, chain1 = template.split('_')
					_, chain2 = partner_template.split('_')
				else:
					print('Something went wrong! Cannot find interfacial residues for:', uniprot_protein, p, interactors[uniprot_protein])
					return
				# check that mutation is within blast alignment for template1
				# do not need to check for template2
				template1 = '_'.join([pdb, chain1])
				alignment = self.all_blast_best_alignment[(uniprot_protein, template1)]
				qstart, qend = int(alignment[4]), int(alignment[5])
				if mut_pos >= qstart and mut_pos <= qend:
					confirmed_partners.append(p)
					if mut_pos in residues: # see if mutation position coincides with a residue
						perturbations.append('1')
					else:
						perturbations.append('0')
				# else:
				# 	print('not within blast alignment:', uniprot_protein, template1, mut_pos, qstart, qend)
			if len(confirmed_partners) != len(perturbations):
				print('Something went wrong!')
			else:
				if len(confirmed_partners) == 0:
					# print('no interactions left for mutation:', uniprot_protein, partners)
					num_removed += 1
				else:
					num_kept += 1
					# add protein partners, perturbations, and edgotype to final_mutation_info
					updated_mutation_info = mutation_info_with_protein_partners[i] + [','.join(confirmed_partners), ','.join(perturbations)]
					# mutation_info_with_protein_partners[i].append(','.join(confirmed_partners))
					# mutation_info_with_protein_partners[i].append(','.join(perturbations))
					if '1' in perturbations: 
						# mutation_info_with_protein_partners[i].append('1')
						updated_mutation_info.append('1')
						num_ir += 1
					else:
						# mutation_info_with_protein_partners[i].append('0')
						updated_mutation_info.append('0')
						num_non_ir += 1

					mutation_info_position_within_blast_alignment.append(updated_mutation_info)
					# mutation_info_position_within_blast_alignment.append(mutation_info_with_protein_partners[i])

		print('Total number of mutations with protein partners:', total_num_mutations)
		print('Number of mutations removed:', num_removed)
		print('Number of mutations kept (with mutation positions within blast alignments):', num_kept)
		print('Number of IR mutations:', num_ir, round(num_ir/num_kept*100, 2))
		print('Number of NIR (non-IR) mutations:', num_non_ir, round(num_non_ir/num_kept*100, 2))

		write_mutation_info_to_file(mutation_info_position_within_blast_alignment, ir_mutations_file)

	# mutation is IR if it maps onto an interfacial residue of at least one interaction
	# maps to:
	# (1) interfacial residues between structural templates
	# (2) interfacial residues between MODELLER constructed homology model (of a heterodimer) that has been repaired using FoldX's RepairPDB
	def find_ir_mutations_all(self, mutation_file, ir_mutations_file):

		header_dict, mutation_info = get_mutation_info(osp.join(self.final_data_dir, mutation_file))
		if len(header_dict) != len(mutation_info[1]): # sanity check
				print('Header & entries have different lengths!', fname, len(header_dict), len(mutation_info[1]))
		
		# get necessary info
		if 'hiunion' in mutation_file:
			# info to update
			interactors, mutation_info_with_protein_partners = self.get_interactors(header_dict, mutation_info, self.hiunion_uniprot_pairs_updated) # key = uniprot protein, value = list of protein partners
			# interfacial residues between structural template
			self.find_ir_mutations(header_dict, interactors, mutation_info_with_protein_partners, self.hiunion_uniprot_pairs_updated, self.all_structural_template_interfacial_residues, ir_mutations_file + '_ir_structural_templates.tsv')
			# print(mutation_info_with_protein_partners[0])
			# print(mutation_info_with_protein_partners[1])
			# interfacial residues between homology model
			self.find_ir_mutations(header_dict, interactors, mutation_info_with_protein_partners, self.hiunion_uniprot_pairs_updated, self.all_homology_model_interfacial_residues, ir_mutations_file + '_ir_homology_models.tsv')
		else: # intact
			# info to update
			interactors, mutation_info_with_protein_partners = self.get_interactors(header_dict, mutation_info, self.intact_uniprot_pairs_updated) # key = uniprot protein, value = list of protein partners
			# interfacial residues between structural template
			self.find_ir_mutations(header_dict, interactors, mutation_info_with_protein_partners, self.intact_uniprot_pairs_updated, self.all_structural_template_interfacial_residues, ir_mutations_file + '_ir_structural_templates.tsv')
			# print(mutation_info_with_protein_partners[0])
			# print(mutation_info_with_protein_partners[1])
			# interfacial residues between homology model
			self.find_ir_mutations(header_dict, interactors, mutation_info_with_protein_partners, self.intact_uniprot_pairs_updated, self.all_homology_model_interfacial_residues, ir_mutations_file + '_ir_homology_models.tsv')

	def get_ir_mutations_all(self):
		
		# get updated uniprot pairs for HI-union and IntAct PPIs
		self.get_updated_uniprot_pairs_all()

		for mutation_file in self.mutation_files_list:
			self.find_ir_mutations_all(osp.join(self.final_data_dir, mutation_file), osp.join(self.final_data_dir, mutation_file[:-17])) # -16 to remove 'nonredundant.tsv'


def main():
	script_dir = osp.dirname(__file__)
	data_dir = osp.join(script_dir, '..', 'data', 'processed')
	final_data_dir = osp.join(script_dir, '..', 'data', 'processed', 'mutations_final')

	interactomes = ['hiunion', 'intact']
	mutation_data = ['dbsnp', 'clinvar', 'cosmic']
	cosmic_datasets = ['genome_screen', 'cgc']
	mutation_files_list = []

	# ONLY FOR MISSENSE MUTATIONS

	for interactome in interactomes:
		for data in mutation_data:
			if data != 'cosmic':
				mutation_files_list.append('_'.join([interactome, 'mapped', data, 'missense', 'nonredundant.tsv']))
			else:
				for cosmic_data in cosmic_datasets:
					mutation_files_list.append('_'.join([interactome, 'mapped', data, cosmic_data, 'missense', 'nonredundant.tsv']))

	i = IRMutation(data_dir, final_data_dir, mutation_files_list)
	i.get_ir_mutations_all()

if __name__=='__main__':
	main()