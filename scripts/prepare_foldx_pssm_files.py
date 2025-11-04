'''
Split FoldX PSSM jobs and create configuration files
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os
import os.path as osp
import argparse
import pickle
from simple_tools import pickle_load, check_create_dir, check_exist, pickle_dump

class PSSM:
	def __init__(self, script_dir, scratch_dir, account, foldx_executable_name):
		self.account = account
		self.foldx_executable_name = foldx_executable_name
		self.repaired_pdb_dir = osp.join(script_dir, 'repaired_pdb')
		self.data_dir = osp.join(script_dir, 'files')
		self.foldx_pssm_dir = osp.join(scratch_dir, 'foldx_pssm_all')
		check_create_dir(self.foldx_pssm_dir)
		self.all_foldx_pssm_mutations = pickle_load(osp.join(self.data_dir, 'all_foldx_pssm_mutations.pickle'))
		self.num_mutations_per_job = 2000
		self.foldx_pssm_num_splits = int(len(self.all_foldx_pssm_mutations)/2000) + 1

	def create_parts(self, num_splits, dirname, mutations, name):
		total_num = 0
		total = []
		for i in range(num_splits):
			check_create_dir(osp.join(dirname, 'split_' + str(i)))
			mutations_for_part = []
			if i != num_splits:
				mutations_for_part = mutations[i*self.num_mutations_per_job:(i+1)*self.num_mutations_per_job]
			else:
				mutations_for_part = mutations[i*self.num_mutations_per_job:len(mutations)]

			# pickle
			pickle_dump(mutations_for_part, osp.join(dirname, name + '_' + str(i) + '.pickle'))

			print('Part', i, 'contains', len(mutations_for_part), name)
			total_num += len(mutations_for_part)
			total += mutations_for_part

		print('Total number of mutations:', total_num)
		if sorted(total) != sorted(mutations):
			print('ERROR! Something went wrong when splitting the mutations!')

	def split_foldx_mutations(self):
		# first check whether or not a split exists already
		if not check_exist(osp.join(self.foldx_pssm_dir, 'foldx_pssm_mutations_0.pickle')):
			self.create_parts(self.foldx_pssm_num_splits, self.foldx_pssm_dir, self.all_foldx_pssm_mutations, 'foldx_pssm_mutations')

	def create_foldx_pssm_config_files(self):
		print('Creating FoldX PSSM configuration files...')
		for i in range(self.foldx_pssm_num_splits):
			print('Creating FoldX PSSM configuration files for part:', i)
			foldx_pssm_mutations_for_part = pickle_load(osp.join(self.foldx_pssm_dir, 'foldx_pssm_mutations_' + str(i) + '.pickle'))
			for foldx_pssm_mutation in foldx_pssm_mutations_for_part:
				dirname = osp.join(self.foldx_pssm_dir, 'split_' + str(i), foldx_pssm_mutation)
				check_create_dir(dirname)
				dimer, foldx_mutation = foldx_pssm_mutation.split('-')
				mut_res = foldx_mutation[-1]
				position = foldx_mutation[:-1] + 'a'
				# write Pssm config file
				with open(osp.join(dirname, foldx_pssm_mutation + '_pssm_config.cfg'), 'w') as f:
					f.write('command=Pssm' + '\n')
					f.write('analyseComplexChains=A,B' + '\n')
					f.write('aminoacids=' + mut_res + '\n')
					f.write('positions=' + position + '\n')
					f.write('pdb-dir=' + self.repaired_pdb_dir + '\n')
					f.write('output-dir=' + dirname + '\n')
					f.write('pdb=' + dimer + '_Repair.pdb' + '\n')
					# the rest of the parameters are all default

	def create_foldx_pssm_slurm_job(self, num_days):
		for i in range(self.foldx_pssm_num_splits):
			print('Creating FoldX PSSM SLURM job for part:', i)
			with open(osp.join(self.foldx_pssm_dir, 'split_' + str(i), 'foldx_pssm_' + str(i) + '.slurm'), 'w') as f:
				f.write('#!/bin/bash' + '\n')
				f.write('#SBATCH -n 1' + '\n')
				f.write('#SBATCH --mem 20G' + '\n')  # actually does need 20G maybe...1G gave an out of memory error and killed the execution of BuildModel
				f.write('#SBATCH --account=' + self.account + '\n')
				f.write('#SBATCH -t ' + str(num_days) + '-00:00' + '\n')
				f.write('#SBATCH -o ' + 'foldx_pssm' + '.%N.%j.log' + '\n')
				f.write('#SBATCH -e ' + 'foldx_pssm' + '.%N.%j.log' + '\n')
				foldx_pssm_mutations_for_part = pickle_load(osp.join(self.foldx_pssm_dir, 'foldx_pssm_mutations_' + str(i) + '.pickle'))
				for foldx_pssm_mutation in foldx_pssm_mutations_for_part:
					dirname = osp.join(self.foldx_pssm_dir, 'split_' + str(i), foldx_pssm_mutation)
					dimer, _ = foldx_pssm_mutation.split('-')
					f.write('./' + self.foldx_executable_name + ' -f ' + osp.join(dirname, foldx_pssm_mutation + '_pssm_config.cfg') + '\n')
					# remove all intermediate/uneeded output files to conserve space (# of files) in scratch directory (maximum # of files = 1000k)
					# ** output files have the same names as the dimer (.pdb) file, need to rename to include the foldx_pssm_mutation in move_binding_ddg_files.bash
					f.write('rm ' + osp.join(dirname, '*.fxout') + '\n')
					f.write('rm ' + osp.join(dirname, 'PSSM_Clash_'+ dimer + '_Repair.txt') + '\n') 
					f.write('rm ' + osp.join(dirname, dimer + '_Repair_1.pdb') + '\n')
					f.write('rm ' + osp.join(dirname, 'WT_' + dimer + '_Repair_1.pdb') + '\n')
					f.write('rm ' + osp.join(dirname, 'individual_list_0_PSSM.txt') + '\n')

	# FoldX PSSM
	def prepare_modeller_foldx_pssm_pdb_files(self):
		self.create_foldx_pssm_config_files()
		self.create_foldx_pssm_slurm_job(1)
	

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--script_dir')  #/home/username/projects/def-*/username/foldx
	parser.add_argument('-c', '--scratch_dir') #/home/username/scratch
	parser.add_argument('-a', '--account') # def-*
	parser.add_argument('-f', '--foldx_executable_name') # foldx_20251231
	args = parser.parse_args()

	print('Creating FoldX configuration files for calculating binding DDG (FoldX PSSM command) upon mutation')
	p = PSSM(args.script_dir, args.scratch_dir, args.account, args.foldx_executable_name)
	p.split_foldx_mutations()
	p.prepare_modeller_foldx_pssm_pdb_files()


if __name__=='__main__':
	main()