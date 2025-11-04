'''
Simple functions that are used repeatedly
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

# simple functions
import pickle
import os.path as osp
import os
from Bio.PDB import PDBList

def pickle_load(pickle_file):
	with open(pickle_file, 'rb') as f:
		contents = pickle.load(f)
		return contents

def pickle_dump(contents, pickle_file):
	with open(pickle_file, 'wb') as f:
		pickle.dump(contents, f)

# check whether a file/directory exists or not
def check_exist(dirname):
	return osp.exists(dirname)

# creates a directory if it doesn't exist already
def check_create_dir(dirname):
	if not osp.exists(dirname):
		os.mkdir(dirname)

# remove directory named dirname if it exists
def check_remove(dirname):
	if osp.exists(dirname):
		os.remove(dirname)

def write_list_of_lists_to_tsv(contents, file):
	with open(file, 'w') as f:
		for line in contents:
			line = [str(item) for item in line]
			f.write('\t'.join(line) + '\n')

# reads in items from newline delimited file
# and returns a list of the items
def get_list_from_file(newline_delimited_file):
	items = []
	with open(newline_delimited_file, 'r') as f:
		items = [line.strip() for line in f]
	return items

# read in list with header from newline delimited file
def read_list_with_header_from_file(file):
	items = []
	with open(file, 'r') as f:
		next(f)
		items = [line.strip() for line in f]
	return items
		
def write_list_to_file(list_to_write, file):
	with open(file, 'w') as f:
		for item in list_to_write:
			f.write(item + '\n')

def write_list_to_file_convert_to_str(list_to_write, file):
	with open(file, 'w') as f:
		for item in list_to_write:
			f.write(str(item) + '\n')

# update self.mutation_info and self.header_dict
def get_mutation_info(mutation_file):
	'''
	returns:
	1) header_dict = {'column1_name': column1_index, 'column2_name': column2_index, ...}
	2) mutation_info = [[line1_column1_entry, line1_column2_entry,...], [line2_column1_entry, line2_column2_entry,...]]
	'''
	header_dict = {}
	mutation_info = []
	with open(mutation_file, 'r') as f:
		column_names = f.readline().strip().split('\t')
		mutation_info.append(column_names) # mutation_info also contains a list of the header column names
		i = 0
		for name in column_names:
			header_dict[name] = i
			i += 1
		for line in f:
			mutation_info.append(line.strip().split('\t'))
	return header_dict, mutation_info

def get_header_dict(mutation_file):
	'''
	returns:
	1) header_dict = {'column1_name': column1_index, 'column2_name': column2_index, ...}
	'''
	header_dict = {}
	f = open(mutation_file, 'r')
	column_names = f.readline().strip().split('\t')
	f.close()
	i = 0
	for name in column_names:
		header_dict[name] = i
		i += 1
	return header_dict

def get_mutation_info_no_header(mutation_file):
	'''
	returns:
	1) mutation_info = [[line1_column1_entry, line1_column2_entry,...], [line2_column1_entry, line2_column2_entry,...]]
	'''
	mutation_info = []
	with open(mutation_file, 'r') as f:
		for line in f:
			mutation_info.append(line.strip().split('\t'))
	return mutation_info

# write to file
# assumes that the first item in mutation info is the header
def write_mutation_info_to_file(mutation_info, mutation_file):
	with open(mutation_file, 'w') as f:
		for line in mutation_info:
			line = [str(item) for item in line]
			f.write('\t'.join(line) + '\n')

# download PDB structure in mmcif format
# gets PDB structure and downloads it
def download_pdb_structure(pdb_structure, pdb_download_path):
	print('-----Downloading PDB structures in mmCIF format-----')
	pdb = PDBList()
	# create pdb directory if it doesn't exist already
	check_create_dir(pdb_download_path)
	# download pdb structure if it doesn't exist in pdb directory already
	if not check_exist(osp.join(pdb_download_path, pdb_structure + '.cif')):
		pdb.retrieve_pdb_file(pdb_structure, pdir=pdb_download_path)

# download PDB structures in mmcif format
# gets PDB structures and downloads them
def download_pdb_structures(pdb_structures_list, pdb_download_path):
	print('-----Downloading PDB structures in mmCIF format-----')
	pdb = PDBList()
	# create pdb directory if it doesn't exist already
	check_create_dir(pdb_download_path)
	# download pdb structure if it doesn't exist in pdb directory already
	for struct in pdb_structures_list:
		if not check_exist(osp.join(pdb_download_path, struct + '.cif')):
			pdb.retrieve_pdb_file(struct, pdir=pdb_download_path)

# create dictionary with keys & lists as values
def initialize_dict_of_lists(keys):
	initialized_dict = {}
	for key in keys:
		initialized_dict[key] = []
	return initialized_dict

# combines multiple lists into one list (removes overlapping items in lists)
def combine_lists(list_of_lists):
	all_items = []
	for l in list_of_lists:
		all_items += l
	return list(set(all_items))

# combines two dictionaries into a new 3rd dict and returns it
def combine_two_dicts(dict1, dict2):
	dict3 = dict1.copy()
	for key in dict2:
		if key not in dict3:
			dict3[key] = dict2[key]
		else:
			print('Overlapping key:', key)
	return dict3

# combined 3 dictionaries into a new 4th dict and returns it
def combine_three_dicts(dict1, dict2, dict3):
	dict4 = dict1.copy()
	for key in dict2:
		if key not in dict4:
			dict4[key] = dict2[key]
		else:
			print('Overlapping key:', key)
	for key in dict3:
		if key not in dict4:
			dict4[key] = dict3[key]
		else:
			print('Overlapping key:', key)
	return dict4

# combines two dict of dicts (dictionary within dictionary) into a new 3rd dict and returns it
# for overlapping keys in dict1 & dict2, combines their values and adds them to dict3
def combine_two_dict_of_dicts_same_keys(dict1, dict2):
	dict3 = dict1.copy()
	for key in dict2:
		if key not in dict3:
			dict3[key] = dict2[key]
		else: # overlapping key, exists in both dict1 and dict2
			# add sub_keys (from inner dict of dict2)
			for sub_key in dict2[key]:
				dict3[key][sub_key] = dict2[key][sub_key]
	return dict3

# convert list of ints to list of strs
def convert_list_of_ints_to_strs(int_list):
	return [str(item) for item in int_list]

# convert list of ints to list of strs in 2 decimal format
def convert_list_of_ints_to_strs_two_decimals(int_list):
	return ["{:.2f}".format(item) for item in int_list]

# convert str of items separated by ',' to list of ints
# e.g. for interfacial residues '16,18,28,293' --> [16, 18, 28, 293]
def convert_str_of_items_to_list_of_ints(str_of_items):
	if str_of_items == '': # empty str (can happen when there are no interfacial residues)
		return []
	else:
		list_of_items_str = str_of_items.split(',')
		return [int(item) for item in list_of_items_str]
	
# to count # of columns in a file
# awk -F'\t' '{print NF; exit}' file

# to count average of a column (e.g. column 52 here) 
# awk -F'\t' '{ total += $52; count++ } END { print total/count }' file