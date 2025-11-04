'''
Combine pickle files across all SLURM jobs into a single pickle file
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os.path as osp
import argparse
from simple_tools import pickle_load, pickle_dump, combine_two_dicts

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--script_dir') # /home/projects/def-*/username/foldx
	parser.add_argument('-n', '--num_splits') # 20
	args = parser.parse_args()

	interactome_dir = osp.join(args.script_dir, 'interactome')
	combined_dict = {}
	combined_dict_pickle_file = osp.join(interactome_dir, 'modeller_heterodimer_interfacial_residues.pickle')

	for i in range(int(args.num_splits)):
		print('Combining split', str(i), 'dictionary...')
		curr_dict = pickle_load(osp.join(interactome_dir, 'modeller_heterodimer_interfacial_residues_' + str(i) + '.pickle'))
		combined_dict = combine_two_dicts(curr_dict, combined_dict)

	pickle_dump(combined_dict, combined_dict_pickle_file)

if __name__ == '__main__':
	main()