'''
Preparing structural models for FoldX energy calculations
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os.path as osp
import argparse
from simple_tools import pickle_load, pickle_dump, check_create_dir, check_exist

# repairing MODELLER contructed PDB files with FoldX RepairPDB
# structures need to be repaired before performing mutagenesis calculations (PSSM and BuildModel)

class Repair:
	def __init__(self, script_dir, foldx_executable_name, num_splits, account, email_address):
		self.interactome_data_dir = osp.join(script_dir, 'interactome')
		self.modeller_dir = osp.join(script_dir, 'modeller')
		self.foldx_pdb_dir = osp.join(script_dir, 'foldx_pdb')
		self.foldx_executable_name = foldx_executable_name
		self.num_splits = int(num_splits)
		self.account = account
		self.email_address = email_address
		self.all_ppi_templates = pickle_load(osp.join(self.interactome_data_dir, 'all_combined_structural_templates_modeled.pickle')) # key = (uniprot1, uniprot2), value = (pdb_chain1, pdb_chain2)
		self.all_ppi_pairs = list(self.all_ppi_templates.keys())

	# splits self.all_ppi_pairs into self.num_splits parts
	# and pickles each part to file
	def create_splits(self):
		total_num = 0
		total = []
		num_ppis_per_part = int(len(self.all_ppi_pairs) / self.num_splits)
		# print('Number PPIs per part:', num_ppis_per_part)
		for i in range(self.num_splits):
			check_create_dir(osp.join(self.foldx_pdb_dir, 'split_' + str(i)))
			ppi_pairs_for_part = []
			if i != self.num_splits-1:
				ppi_pairs_for_part = self.all_ppi_pairs[i*num_ppis_per_part:(i+1)*num_ppis_per_part]
			else:
				ppi_pairs_for_part = self.all_ppi_pairs[i*num_ppis_per_part:len(self.all_ppi_pairs)]

			# pickle
			pickle_dump(ppi_pairs_for_part, osp.join(self.foldx_pdb_dir, 'split_' + str(i), 'ppi_pairs_for_repair_pdb_' + str(i) + '.pickle'))

			print('Part', i, 'contains', len(ppi_pairs_for_part), 'PPI pairs')
			total_num += len(ppi_pairs_for_part)
			total += ppi_pairs_for_part

	def create_repair_pdb_bash_files(self):
		for i in range(self.num_splits):
			with open(osp.join(self.foldx_pdb_dir, 'split_' + str(i), 'repair_pdb_' + str(i) + '.bash'), 'w') as f:
				f.write('#!/bin/bash' + '\n')
				all_ppi_pairs_for_part = pickle_load(osp.join(self.foldx_pdb_dir, 'split_' + str(i), 'ppi_pairs_for_repair_pdb_' + str(i) + '.pickle'))
				for (p1, p2) in all_ppi_pairs_for_part:
					pdb_chain1, pdb_chain2 = self.all_ppi_templates[(p1, p2)]
					repaired_pdb_fname = '_'.join([p1, pdb_chain1, p2, pdb_chain2])
					modeller_fname = '_'.join([p1, pdb_chain1, p2, pdb_chain2]) + '.B99990001.pdb'
					if not check_exist(osp.join(self.foldx_pdb_dir, 'split_' + str(i), repaired_pdb_fname + '_Repair.pdb')):
						f.write('./' + self.foldx_executable_name + ' --command=RepairPDB --pdb=' + modeller_fname + ' --pdb-dir=' + self.modeller_dir + ' --output-dir=' + osp.join(self.foldx_pdb_dir, 'split_' + str(i)) + '\n')
						# f.write('rm ' + osp.join(self.foldx_pdb_dir, 'split_' + str(i), repaired_pdb_fname + '_Repair.fxout' + '\n')) # remove intermediate file

	def create_repair_pdb_slurm_jobs(self):
		for i in range(self.num_splits):
			with open(osp.join(self.foldx_pdb_dir, 'split_' + str(i), 'repair_pdb_' + str(i) + '.slurm'), 'w') as f:
				f.write('#!/bin/bash' + '\n')
				# f.write('#SBATCH -n 5' + '\n')
				f.write('#SBATCH --mem-per-cpu 5G' + '\n')
				f.write('#SBATCH -t 2-00:00' + '\n') # 2 days
				f.write('#SBATCH -o repair_pdb_' + str(i) + '.%N.%j.log' + '\n')
				f.write('#SBATCH -e repair_pdb_' + str(i) + '.%N.%j.log' + '\n')
				f.write('#SBATCH --mail-type=END,FAIL' + '\n') 
				f.write('#SBATCH --mail-user=' + self.email_address + '\n') # ting-yi.su@mail.mcgill.ca
				f.write('#SBATCH --account=' + self.account + '\n') # def-yxia
				f.write('bash repair_pdb_' + str(i) + '.bash' '\n')

	def split_create_repair_pdb_jobs(self):
		self.create_splits()
		self.create_repair_pdb_bash_files()
		self.create_repair_pdb_slurm_jobs()

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--script_dir') # /home/projects/def-*/username/foldx
	parser.add_argument('-f', '--foldx_executable_name') # foldx_20251231
	parser.add_argument('-n', '--num_splits') # 20
	parser.add_argument('-a', '--account') # def-* (compute canada project allocation)
	parser.add_argument('-e', '--email_address') # your email address
	args = parser.parse_args()

	r = Repair(args.script_dir, args.foldx_executable_name, args.num_splits, args.account, args.email_address)
	r.split_create_repair_pdb_jobs()

if __name__ == '__main__':
	main()