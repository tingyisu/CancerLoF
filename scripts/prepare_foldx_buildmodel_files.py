'''
Split FoldX BuildModel jobs and create configuration files
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os
import os.path as osp
import argparse
import pickle
from simple_tools import pickle_load, check_create_dir, check_exist, pickle_dump

# to calculate folding DDG (with FoldX BuildModel), used dimers and just specified monomer chain (either 'A' or 'B')

class BuildModel:
	def __init__(self, script_dir, scratch_dir, account, foldx_executable_name):
		self.account = account
		self.foldx_executable_name = foldx_executable_name
		self.data_dir = osp.join(script_dir, 'files')
		self.repaired_pdb_dir = osp.join(script_dir, 'repaired_pdb')
		self.foldx_buildmodel_dir = osp.join(scratch_dir, 'foldx_buildmodel_all')
		check_create_dir(self.foldx_buildmodel_dir)
		self.all_foldx_buildmodel_mutations = pickle_load(osp.join(self.data_dir, 'all_foldx_buildmodel_mutations_combined.pickle')) # plus additional (self.all_foldx_buildmodel_mutations + self.additional_foldx_buildmodel_mutations)
		self.num_mutations_per_job = 2000
		self.foldx_buildmodel_num_splits = int(len(self.all_foldx_buildmodel_mutations)/2000) + 1

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
		if not check_exist(osp.join(self.foldx_buildmodel_dir, 'foldx_buildmodel_mutations_0.pickle')):
			self.create_parts(self.foldx_buildmodel_num_splits, self.foldx_buildmodel_dir, self.all_foldx_buildmodel_mutations, 'foldx_buildmodel_mutations')

	def create_foldx_buildmodel_config_files(self):
		for i in range(self.foldx_buildmodel_num_splits):
			print('Creating FoldX BuildModel configuration files for part:', i)
			foldx_buildmodel_mutations_for_part = pickle_load(osp.join(self.foldx_buildmodel_dir, 'foldx_buildmodel_mutations_' + str(i) + '.pickle'))
			for foldx_buildmodel_mutation in foldx_buildmodel_mutations_for_part:
				dirname = osp.join(self.foldx_buildmodel_dir, 'split_' + str(i), foldx_buildmodel_mutation)
				check_create_dir(dirname)
				dimer, foldx_mutation = foldx_buildmodel_mutation.split('-')
				# create mutant file
				with open(osp.join(dirname, 'individual_list.txt'), 'w') as f_mut:
					f_mut.write(foldx_mutation + ';')
				# write BuildModel config file
				with open(osp.join(dirname, foldx_buildmodel_mutation + '_buildmodel_config.cfg'), 'w') as f:
					f.write('command=BuildModel' + '\n')
					f.write('pdb-dir=' + self.repaired_pdb_dir + '\n')
					f.write('output-dir=' + dirname + '\n')
					f.write('pdb=' + dimer + '_Repair.pdb' + '\n')
					f.write('mutant-file=' + osp.join(dirname, 'individual_list.txt') + '\n')
					# the rest of the parameters are all default

	def create_foldx_buildmodel_slurm_job(self, num_days):
		for i in range(self.foldx_buildmodel_num_splits):
			print('Creating FoldX BuildModel SLURM job for part:', i)
			with open(osp.join(self.foldx_buildmodel_dir, 'split_' + str(i), 'foldx_buildmodel_' + str(i) + '.slurm'), 'w') as f:
				f.write('#!/bin/bash' + '\n')
				f.write('#SBATCH -n 1' + '\n')
				f.write('#SBATCH --mem 20G' + '\n') # actually does need 20G maybe...1G gave an out of memory error and killed the execution of BuildModel
				f.write('#SBATCH --account=' + self.account + '\n')
				f.write('#SBATCH -t ' + str(num_days) + '-00:00' + '\n')
				f.write('#SBATCH -o ' + 'foldx_buildmodel' + '.%N.%j.log' + '\n')
				f.write('#SBATCH -e ' + 'foldx_buildmodel' + '.%N.%j.log' + '\n')
				foldx_buildmodel_mutations_for_part = pickle_load(osp.join(self.foldx_buildmodel_dir, 'foldx_buildmodel_mutations_' + str(i) + '.pickle'))
				for foldx_buildmodel_mutation in foldx_buildmodel_mutations_for_part:
					dirname = osp.join(self.foldx_buildmodel_dir, 'split_' + str(i), foldx_buildmodel_mutation)
					dimer, _ = foldx_buildmodel_mutation.split('-')
					f.write('./' + self.foldx_executable_name + ' -f ' + osp.join(dirname, foldx_buildmodel_mutation + '_buildmodel_config.cfg') + '\n')
					# remove all intermediate/uneeded output files to conserve space (# of files) in scratch directory (maximum # of files = 1000k)
					f.write('rm ' + osp.join(dirname, 'Average_' + dimer + '_Repair.fxout') + '\n')
					f.write('rm ' + osp.join(dirname, 'PdbList_' + dimer + '_Repair.fxout') + '\n')
					f.write('rm ' + osp.join(dirname, dimer + '_Repair_1.pdb') + '\n') # this is the PDB file with the mutation residue?
					f.write('rm ' + osp.join(dirname, 'Raw_' + dimer + '_Repair.fxout') + '\n')
					f.write('rm ' + osp.join(dirname, 'WT_' + dimer + '_Repair_1.pdb') + '\n') # this is the same file as the original .pdb file?

	# FoldX BuildModel
	def prepare_modeller_foldx_buildmodel_pdb_files(self):
		self.create_foldx_buildmodel_config_files()
		self.create_foldx_buildmodel_slurm_job(1)

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--script_dir')  #/home/tingyisu/projects/def-yxia/tingyisu/modeller
	parser.add_argument('-c', '--scratch_dir') #/home/tingyisu/scratch
	parser.add_argument('-a', '--account') # def-yxia
	parser.add_argument('-f', '--foldx_executable_name') # foldx_20241231
	args = parser.parse_args()

	print('Preparing modeller PDB files for calculating folding DDG (FoldX BuildModel) upon mutation')
	b = BuildModel(args.script_dir, args.scratch_dir, args.account, args.foldx_executable_name)
	b.split_foldx_mutations()
	b.prepare_modeller_foldx_buildmodel_pdb_files()


if __name__=='__main__':
	main()