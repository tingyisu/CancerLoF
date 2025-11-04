'''
Removed unmodeled PPIs which MODELLER was unable to build structures for
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os.path as osp
from simple_tools import pickle_load, pickle_dump, get_mutation_info, write_mutation_info_to_file

class RemoveUnmodeledPPIs:
	def __init__(self, interactome_data_dir):
		self.uncompleted_models = pickle_load(osp.join(interactome_data_dir, 'uncompleted_models.pickle'))
		self.all_ppi_templates = pickle_load(osp.join(interactome_data_dir, 'all_combined_structural_templates.pickle')) # key = (uniprot1, uniprot2), value = (pdb_chain1, pdb_chain2)
		self.all_interfacial_residues = pickle_load(osp.join(interactome_data_dir, 'all_combined_interfacial_residues.pickle'))
		self.all_combined_structural_interactome_file = osp.join(interactome_data_dir, 'all_combined_structural_interactome.tsv')

		self.all_ppi_templates_modeled_pickle_file = osp.join(interactome_data_dir, 'all_combined_structural_templates_modeled.pickle')
		self.all_interfacial_residues_modeled_pickle_file = osp.join(interactome_data_dir, 'all_combined_interfacial_residues_modeled.pickle')
		self.all_combined_structural_interactome_file_modeled = osp.join(interactome_data_dir, 'all_combined_structural_interactome_modeled.tsv')
	
	def remove_ppis(self):
		uncompleted_ppis = []

		for uncompleted_model in self.uncompleted_models:
			p1, pdb, p1_chain, p2, _, p2_chain = uncompleted_model.split('_')
			uncompleted_ppis.append((p1, p2))

		print('Initial number of combined PPIs:', len(self.all_ppi_templates))

		for uncompleted_ppi in uncompleted_ppis:
			del self.all_ppi_templates[uncompleted_ppi]
			del self.all_interfacial_residues[uncompleted_ppi]

		pickle_dump(self.all_ppi_templates, self.all_ppi_templates_modeled_pickle_file)
		pickle_dump(self.all_interfacial_residues, self.all_interfacial_residues_modeled_pickle_file)

		header_dict, ppi_info = get_mutation_info(self.all_combined_structural_interactome_file)
		new_ppi_info = [ppi_info[0]]
		for ppi in ppi_info[1:]:
			p1, p2 = ppi[header_dict['uniprot1']], ppi[header_dict['uniprot2']]
			if (p1, p2) not in uncompleted_ppis:
				new_ppi_info.append(ppi)

		write_mutation_info_to_file(new_ppi_info, self.all_combined_structural_interactome_file_modeled)

		print('Number of combined PPIs with homology models:', len(self.all_ppi_templates))


def main():
	script_dir = osp.dirname(__file__)
	interactome_data_dir = osp.join(script_dir, '..', 'data', 'processed', 'interactome')

	r = RemoveUnmodeledPPIs(interactome_data_dir)
	r.remove_ppis()

if __name__ == '__main__':
	main()