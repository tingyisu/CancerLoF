'''
Creating multiple SLURM jobs for filding intefacial residues in MODELLER heterodimers
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import argparse
import os.path as osp

def create_interfacial_residue_slurm_jobs(script_dir, num_splits, account, email_address):
	for i in range(num_splits):
		with open(osp.join(script_dir, 'interfacial_residues_' + str(i) + '.slurm'), 'w') as f:
			f.write('#!/bin/bash' + '\n')
			f.write('#SBATCH --mem-per-cpu 1G' + '\n')
			f.write('#SBATCH -t 10:00:00' + '\n') # 10 hours
			f.write('#SBATCH -o interfacial_residues_' + str(i) + '.%N.%j.log' + '\n')
			f.write('#SBATCH -e interfacial_residues_' + str(i) + '.%N.%j.log' + '\n')
			f.write('#SBATCH --mail-type=END,FAIL' + '\n') 
			f.write('#SBATCH --mail-user=' + email_address + '\n') # ting-yi.su@mail.mcgill.ca
			f.write('#SBATCH --account=' + account + '\n') # def-yxia
			f.write('python -u find_interfacial_residues_in_modeller_heterodimer.py -s ' + script_dir + ' -i ' + str(i) + '\n')

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--script_dir') # /home/projects/def-*/username/foldx
	parser.add_argument('-n', '--num_splits') # 20
	parser.add_argument('-a', '--account') # def-* (compute canada project allocation)
	parser.add_argument('-e', '--email_address') # your email address
	args = parser.parse_args()

	create_interfacial_residue_slurm_jobs(args.script_dir, int(args.num_splits), args.account, args.email_address)

if __name__ == '__main__':
	main()