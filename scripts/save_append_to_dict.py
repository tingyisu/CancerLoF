'''
Pickles DDG calculations into a single pickle file (to save on storage space and make edgotyping easier and faster)
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os.path as osp
import os
from simple_tools import pickle_load, pickle_dump, check_exist
import argparse
import shutil
import math

class SaveDict:
	def __init__(self, folder_dir, pickled_dict_to_append_to, datatype):
		self.folder_dir = folder_dir
		self.pickled_dict_to_append_to = pickled_dict_to_append_to
		self.dict_to_append_to = {}
		if check_exist(self.pickled_dict_to_append_to): # if pickled dict file exists, load it; otherwise, create an empty dict
			self.dict_to_append_to = pickle_load(self.pickled_dict_to_append_to)
		self.datatype = datatype # datatype = 'folding' or 'binding' ('folding'/'binding' = DDG energy calculations)

	# get binding DDG
	def append_to_pssm_dict(self):
		print('Appending info from', self.folder_dir, 'to PSSM dict', self.pickled_dict_to_append_to, '...')
		i = 0
		for file in os.listdir(self.folder_dir):
			i += 1
			print('Getting binding DDG from FoldX PSSM file #:', i, file[5:-4])
			pssm_file = osp.join(self.folder_dir, file)
			# print('Getting binding DDG from FoldX PSSM file #:', i, pssm_file)
			with open(pssm_file, 'r') as f:
				lines = f.readlines()
				# mut_res = lines[0].strip()
				_, binding_ddg = lines[1].strip().split('\t')
				foldx_pssm_mutation = file[5:-4]
				print(foldx_pssm_mutation, binding_ddg)
				self.dict_to_append_to[foldx_pssm_mutation] = float(binding_ddg)

		pickle_dump(self.dict_to_append_to, self.pickled_dict_to_append_to)

	# get folding DDG
	def append_to_buildmodel_dict(self): 

		print('Appending info from', self.folder_dir, 'to BuildModel dict', self.pickled_dict_to_append_to, '...')
		i = 0
		for file in os.listdir(self.folder_dir):
			i += 1
			print('Getting folding DDG from FoldX BuildModel file #:', i, file[4:-6])
			buildmodel_file = osp.join(self.folder_dir, file)
			with open(buildmodel_file, 'r') as f:
				folding_ddg = f.readlines()[-1].split('\t')[1]
				foldx_buildmodel_mutation = file[4:-6]
				print(foldx_buildmodel_mutation, folding_ddg)
				self.dict_to_append_to[foldx_buildmodel_mutation] = float(folding_ddg)

		pickle_dump(self.dict_to_append_to, self.pickled_dict_to_append_to)


	def append(self): 
		if self.datatype == 'folding':
			self.append_to_buildmodel_dict()
		elif self.datatype == 'binding':
			self.append_to_pssm_dict()
		else:
			print('ERROR! Datatype provided is not an accepted type! Should be either folding or binding!')

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--folder_dir')  #e.g. /home/username/scratch/folding_ddg
	parser.add_argument('-p', '--pickled_dict_to_append_to') # e.g. /home/username/projects/def-*/username/foldx/files/foldx_buildmodel_folding_ddg_dict.pickle
	parser.add_argument('-d', '--datatype') #  datatype = 'folding', 'binding'
	args = parser.parse_args()

	s = SaveDict(args.folder_dir, args.pickled_dict_to_append_to, args.datatype)
	s.append()

if __name__ == '__main__':
	main()