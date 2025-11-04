'''
Finding interfacial residues in MODELLER heterodimers
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

from Bio.PDB.PDBParser import PDBParser 
import timeit
import argparse
import os.path as osp
from simple_tools import pickle_load, pickle_dump, convert_list_of_ints_to_strs, write_mutation_info_to_file

# NOTE: MODELLER typically converts HETATOM entries to ATOM entries when building homology models!! 
# so don't have any HETATOM entries in newly constructed .pdb files
# (also did not check distances between heteroatoms in structural templates, so also not done here)

class InterfacialResidues:
	def __init__(self, script_dir, split_num):
		self.interactome_data_dir = osp.join(script_dir, 'interactome')
		self.repaired_pdb_dir = osp.join(script_dir, 'repaired_pdb')
		self.foldx_pdb_dir = osp.join(script_dir, 'foldx_pdb')
		self.all_ppi_templates = pickle_load(osp.join(self.interactome_data_dir, 'all_combined_structural_templates_modeled.pickle')) # key = (uniprot1, uniprot2), value = (pdb_chain1, pdb_chain2)
		self.ppi_pairs_in_split = pickle_load(osp.join(self.foldx_pdb_dir, 'split_' + split_num, 'ppi_pairs_for_repair_pdb_' + split_num + '.pickle')) # list of (uniprot1, uniprot2), value = (pdb_chain1, pdb_chain2)
		
		self.homology_model_interfacial_residues = {} # key = (uniprot1, uniprot2), value = [residue1_str, residue2_str]]
		
		# pickle and write to file
		self.homology_model_interfacial_residues_pickle_file = osp.join(self.interactome_data_dir, 'modeller_heterodimer_interfacial_residues_' + split_num + '.pickle')
		self.homology_model_interfacial_residues_file = osp.join(self.interactome_data_dir, 'modeller_heterodimer_structural_interactome_' + split_num + '.tsv')

	
	def get_interfacial_residues(self, struct_name):
		
		# create parser
		parser = PDBParser(QUIET=True) # suppress warnings
		fname = osp.join(self.repaired_pdb_dir, struct_name + '_Repair.pdb')
		# print(fname)
		# print(fname)
		structure = ""
		try:
			structure = parser.get_structure(struct_name, fname)
		except Exception as e: # redownload .cif file if need be
			print('Unable to load structure from provided PDB file')
		# if not osp.exists(osp.join(pdb_download_path, fname)):

		model = structure[0]
		chain1_name, chain2_name = 'A', 'B'
		num_atoms = 0
		start = timeit.default_timer() 
		# get interfacial residues (residues that are within 5 angstroms of each other)
		chain1, chain2 = model[chain1_name], model[chain2_name]
		residue1_list, residue2_list = [], []
		for residue1 in chain1:    
			if residue1.get_full_id()[3][0] == " ": #ignore the heteroatoms 
				res1_num = str(residue1.get_full_id()[3][1])
				res1_id = str(residue1.get_full_id()[3][2]) # add in the auth_seq_id identifier (if exists), otherwise could have multiple label_seq_id for the same auth_seq_id (e.g. for auth_seq_id 156 on chain A of 1kmc --> has 156 and 156A)
				res1 = res1_num
				if res1_id != ' ':
					res1 = res1_num + res1_id
				for residue2 in chain2:
					if residue2.get_full_id()[3][0] == " ": #ignore the heteroatoms 
						res2_num = str(residue2.get_full_id()[3][1])
						res2_id = str(residue2.get_full_id()[3][2]) # add in the auth_seq_id identifier, otherwise could have multiple label_seq_id for the same auth_seq_id (e.g. for auth_seq_id 156 on chain A of 1kmc --> has 156 and 156A)
						res2 = res2_num
						if res2_id != ' ':
							res2 = res2_num + res2_id
						# compute distance between all atoms
						try:
							break_loop = False
							# print(f'Looking at {res1_num}, {res2_num}')
							for atom1 in residue1: 
								for atom2 in residue2:
									# ways to calculate euclidean distance
									# 1. np.linalg.norm(residue1[atom1.name].get_coord() - residue2[atom2.name].get_coord())
									# 2. residue1[atom1.name] - residue2[atom2.name] 
									# atom1 - atom2 < 5 --> might be faster
									# check distances between x, y, z coordinates
									# if atom1[0] - atom2[0] > 5 
									num_atoms += 1
									atom1_coords, atom2_coords = atom1.get_coord(), atom2.get_coord()
									# don't keep searching if distance between each coordinate (x, y, z) > 5
									if abs(atom1_coords[0] - atom2_coords[0]) > 5 or abs(atom1_coords[1] - atom2_coords[1]) > 5 or abs(atom1_coords[2] - atom2_coords[2]) > 5:
										continue
									else:
										if atom1 - atom2 <= 5: # if the euclidean distance between any atoms in residues 1 & 2 <= 5, then they are interacting
											if res1 not in residue1_list: residue1_list.append(res1)
											if res2 not in residue2_list: residue2_list.append(res2)
											break_loop = True
											break
								if break_loop: break
								# distance = residue1['CA'] - residue2['CA']
						except KeyError:
							## no CA atom, e.g. for H_NAG
							continue
								# print(residue1, residue2, distance)
		stop = timeit.default_timer()
		print('\t'.join([struct_name, chain1_name, chain2_name]) + '; time taken: ' + str(stop-start) + '; num atoms:' + str(num_atoms))
		# residue list in int format
		residue1_list_int = [int(res) for res in residue1_list]
		residue2_list_int = [int(res) for res in residue2_list]

		return sorted(residue1_list_int), sorted(residue2_list_int)


	def get_interfacial_residues_for_all_heterodimers(self):
		info = [['uniprot1', 'uniprot2', 'pdb_structure', 'chain1', 'chain2', 'interfacial_residues1', 'interfacial_residues2']]
		for (p1, p2) in self.ppi_pairs_in_split:
			pdb_chain1, pdb_chain2 = self.all_ppi_templates[(p1, p2)]
			pdb, chain1 = pdb_chain1.split('_')
			_, chain2 = pdb_chain2.split('_')
			struct_name = '_'.join([p1, pdb_chain1, p2, pdb_chain2])
			residue1_list_int, residue2_list_int = self.get_interfacial_residues(struct_name)
			residue1_list_str, residue2_list_str = ','.join(convert_list_of_ints_to_strs(residue1_list_int)), ','.join(convert_list_of_ints_to_strs(residue2_list_int))
			self.homology_model_interfacial_residues[(p1, p2)] = [residue1_list_str, residue2_list_str]
			
			# write to file
			info.append([p1, p2, pdb, chain1, chain2, residue1_list_str, residue2_list_str])
		
		write_mutation_info_to_file(info, self.homology_model_interfacial_residues_file)
		pickle_dump(self.homology_model_interfacial_residues, self.homology_model_interfacial_residues_pickle_file)


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--script_dir') # /home/projects/def-*/username/foldx
	parser.add_argument('-i', '--split_num') # 20
	args = parser.parse_args()

	i = InterfacialResidues(args.script_dir, args.split_num)
	i.get_interfacial_residues_for_all_heterodimers()

if __name__ == '__main__':
	main()