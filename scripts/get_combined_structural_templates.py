'''
Combines HI-union and IntAct structural templates & picks one as representative for building 3D MODELLER structural model
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

from simple_tools import pickle_load, pickle_dump, write_mutation_info_to_file, convert_list_of_ints_to_strs
import os.path as osp

class HeterodimerAlignments:
	def __init__(self, interactome_data_dir):
		self.hiunion_selected_templates = pickle_load(osp.join(interactome_data_dir, 'hiunion_selected_pdb_structural_template.pickle')) # format: {(Protein1, Protein2):'4awl_A_C'}
		self.intact_selected_templates = pickle_load(osp.join(interactome_data_dir, 'intact_selected_pdb_structural_template.pickle'))
		self.hiunion_selected_interfacial_residues = pickle_load(osp.join(interactome_data_dir, 'hiunion_interfacial_residues_for_selected_template.pickle')) # format {(Protein1, Protein2): [residue1_list, residue2_list]}}
		self.intact_selected_interfacial_residues = pickle_load(osp.join(interactome_data_dir, 'intact_interfacial_residues_for_selected_template.pickle'))
		
		self.all_ppi_templates = {} # key = (uniprot1, uniprot2), value = (pdb_chain1, pdb_chain2)
		self.all_interfacial_residues = {} # key = (uniprot1, uniprot2), value = [residue1_str, residue2_str]]
		self.all_ppi_num_interfacial_residues = {} # (uniprot1, uniprot2), value = total number of interfacial residues

		# pickled files
		self.all_ppi_templates_pickle_file = osp.join(interactome_data_dir, 'all_combined_structural_templates.pickle')
		self.all_interfacial_residues_pickle_file = osp.join(interactome_data_dir, 'all_combined_interfacial_residues.pickle')
		self.all_combined_structural_interactome_file = osp.join(interactome_data_dir, 'all_combined_structural_interactome.tsv')

	def combine_selected_templates(self, selected_templates, selected_interfacial_residues):
		for (p1, p2) in selected_templates:
			structural_templates = selected_templates[(p1, p2)] # format: selected_templates[(p1, p2)]
			pdb, chain1, chain2 = structural_templates.split('_')
			curr_pdb_chain1, curr_pdb_chain2 = '_'.join([pdb, chain1]), '_'.join([pdb, chain2])
			residue1_list, residue2_list = selected_interfacial_residues[(p1, p2)]
			if (p1, p2) not in self.all_ppi_templates and (p2, p1) not in self.all_ppi_templates: # add only one (undirected interaction)
				self.all_ppi_templates[(p1, p2)] = (curr_pdb_chain1, curr_pdb_chain2)
				self.all_ppi_num_interfacial_residues[(p1, p2)] = len(residue1_list) + len(residue2_list)
				self.all_interfacial_residues[(p1, p2)] = [','.join(convert_list_of_ints_to_strs(residue1_list)), ','.join(convert_list_of_ints_to_strs(residue2_list))]
			else: # see if have different pdb structure/chain combo for a given interaction
				# yes, does happen, but seems like the two chains are identical, (i.e. both p1 & p2 have the best alignments with the exact same chains)
				prev_pdb_chain1, prev_pdb_chain2 = '', ''
				reverse = False
				if (p1, p2) in self.all_ppi_templates:
					prev_pdb_chain1, prev_pdb_chain2 = self.all_ppi_templates[(p1, p2)]
				elif (p2, p1) in self.all_ppi_templates:
					prev_pdb_chain2, prev_pdb_chain1 = self.all_ppi_templates[(p2, p1)]
					reverse = True
				else:
					print('something went wrong!')

				# sanity check
				# checked that most are similar, so can just choose one representative over both HI-union and IntAct
				if prev_pdb_chain1 != curr_pdb_chain1 or prev_pdb_chain2 != curr_pdb_chain2:
					if {prev_pdb_chain1, curr_pdb_chain1} != {prev_pdb_chain2, curr_pdb_chain2}:
						# pick the templates with the most # of interfacial residues
						# get prev # of interfacial residues
						prev_num = 0
						if reverse:
							prev_num = self.all_ppi_num_interfacial_residues[(p2, p1)]
						else:
							prev_num = self.all_ppi_num_interfacial_residues[(p1, p2)]

						if prev_num < len(residue1_list) + len(residue2_list):
							# update
							self.all_ppi_templates[(p1, p2)] = (curr_pdb_chain1, curr_pdb_chain2)
							self.all_ppi_num_interfacial_residues[(p1, p2)] = len(residue1_list) + len(residue2_list)
							self.all_interfacial_residues[(p1, p2)] = [','.join(residue1_list), ','.join(residue2_list)]
							print('not the same structural templates, different # of interfacial residues; update:', (p1, p2), (prev_pdb_chain1, prev_pdb_chain2), (curr_pdb_chain1, curr_pdb_chain2))
						else:
							print('not the same structural templates, but same # of interfacial residues:', (p1, p2), (prev_pdb_chain1, prev_pdb_chain2), (curr_pdb_chain1, curr_pdb_chain2))
		pickle_dump(self.all_ppi_templates, self.all_ppi_templates_pickle_file)
		pickle_dump(self.all_interfacial_residues, self.all_interfacial_residues_pickle_file)

	def write_combined_structural_interactome_to_file(self):
		info = [['uniprot1', 'uniprot2', 'pdb_structure', 'chain1', 'chain2', 'interfacial_residues1', 'interfacial_residues2']]
		for (p1, p2) in self.all_ppi_templates:
			pdb_chain1, pdb_chain2 = self.all_ppi_templates[(p1, p2)]
			pdb = pdb_chain1.split('_')[0]
			chain1, chain2 = pdb_chain1.split('_')[1], pdb_chain2.split('_')[1]
			residue1_list, residue2_list = self.all_interfacial_residues[(p1, p2)]
			info.append([p1, p2, pdb, chain1, chain2, residue1_list, residue2_list])
		write_mutation_info_to_file(info, self.all_combined_structural_interactome_file)


	def get_modeller_heterodimers(self):
		# retrieves list of heterodimer models (representing PPIs & their selected templates) to build using MODELLER
		self.combine_selected_templates(self.hiunion_selected_templates, self.hiunion_selected_interfacial_residues)
		self.combine_selected_templates(self.intact_selected_templates, self.intact_selected_interfacial_residues)
		self.write_combined_structural_interactome_to_file()
		print('Number of combined PPIs from HI-union and IntAct:', len(self.all_ppi_templates))

def main():
	script_dir = osp.dirname(__file__)
	interactome_data_dir = osp.join(script_dir, '..', 'data', 'processed', 'interactome')

	a = HeterodimerAlignments(interactome_data_dir)
	a.get_modeller_heterodimers()

if __name__ == '__main__':
	main()