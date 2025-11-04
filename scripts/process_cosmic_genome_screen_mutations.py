'''
Processes COSMIC genome screen mutations (in batches)
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os.path as osp
import argparse
from simple_tools import get_mutation_info_no_header, get_header_dict, pickle_load, write_mutation_info_to_file

class COSMIC:
	def __init__(self, original_data_dir, processed_data_dir, split_number):
		self.processed_data_dir = processed_data_dir
		self.original_data_dir = original_data_dir
		self.ensembl_prot_gene_name_uniprot_dict = pickle_load(osp.join(self.processed_data_dir, 'ensembl_prot_gene_name_uniprot_dict.pickle'))
		self.filename_split_number = self.get_filename_split_number(split_number)
		self.mutation_file = osp.join(original_data_dir, 'cosmic_genome_screen_mutants_split_' + self.filename_split_number)
		self.header_dict = get_header_dict(osp.join(original_data_dir, 'header.tsv')) # get header_dict from header.tsv
		self.processed_missense_mutation_file = osp.join(original_data_dir, 'processed_cosmic_genome_screen_mutants_missense_split_' + self.filename_split_number + '.tsv')
		self.processed_synonymous_mutation_file = osp.join(original_data_dir, 'processed_cosmic_genome_screen_mutants_synonymous_split_' + self.filename_split_number + '.tsv')
		self.processed_nonsense_mutation_file = osp.join(original_data_dir, 'processed_cosmic_genome_screen_mutants_nonsense_split_' + self.filename_split_number + '.tsv')

	def get_filename_split_number(self, split_number): # split_number is a str
		filename_split_number = split_number
		if int(split_number)/10 < 1: # if number is smaller than 10
			filename_split_number = '0' + str(split_number)
		return filename_split_number


	def get_cdna_info(self, mutation_cds):
		'''
		example of mutation_cds: c.194G>C
		'''
		# print(mutation_cds)
		to_continue_to_split, cdna_alt_allele = mutation_cds.split('>') # split by '>'
		cdna_allele_pos, cdna_ref_allele = to_continue_to_split[2:-1], to_continue_to_split[-1]
		# print((cdna_ref_allele, cdna_allele_pos, cdna_allele))
		return (cdna_ref_allele, cdna_allele_pos, cdna_alt_allele)	


	def get_prot_info(self, mutation_aa):
		'''
		example of mutation_aa: p.S65T
		'''
		ref_res = mutation_aa[2]
		prot_res_pos = mutation_aa[3:len(mutation_aa)-1]
		alt_res = mutation_aa[-1]
		return (ref_res, prot_res_pos, alt_res)


	def get_info_for_new_mutation_info_columns(self, mutation_aa, mutation_cds, ensembl_protein, gene_name):
		if '>' in mutation_cds and 'fs' not in mutation_aa:
			cdna_ref_allele, cdna_allele_pos, cdna_alt_allele = self.get_mrna_info(mutation_cds)
			ref_res, prot_res_pos, alt_res = self.get_prot_info(mutation_aa)
			if prot_res_pos.isdigit(): # could maybe get an invalid prot_res_pos (e.g. from breaking down p.Sec48L), prot_res_pos becomes ec48
				uniprot_protein = ''
				if (ensembl_protein, gene_name) in self.ensembl_prot_gene_name_uniprot_dict:
					uniprot_protein = self.ensembl_prot_gene_name_uniprot_dict[(ensembl_protein, gene_name)]

				if uniprot_protein != '':
					if alt_res != '=':
						return [cdna_ref_allele, cdna_allele_pos, cdna_alt_allele, ref_res, prot_res_pos, alt_res, uniprot_protein]
					else: # if alt_res == '=', this means that it is a synonymous mutation, so alt_res == ref_res
						return [cdna_ref_allele, cdna_allele_pos, cdna_alt_allele, ref_res, prot_res_pos, ref_res, uniprot_protein]
				else:
					return [] # no corresponding uniprot protein, so discard mutation
			else:
				print('Strange HGVSP format:', mutation_aa) # typically is b/c of Sec...which stands for Selenocysteine and is not one of the 20 naturally occurring amino acids
				return []
		else:
			# print('Strange mutation format:', mutation_aa, mutation_cds)
			return []


	def process_genome_screen_mutations(self):
		mutation_info = get_mutation_info_no_header(self.mutation_file)
		column_names_list = list(self.header_dict.keys())
		column_names_list.extend(['ensembl_mrna', 'cdna_ref_allele', 'cdna_allele_pos', 'cdna_alt_allele', 'ensembl_protein', 'ref_res', 'prot_res_pos', 'alt_res', 'uniprot_protein'])

		missense_info = [column_names_list]
		synonymous_info = [column_names_list]
		nonsense_info = [column_names_list]

		for mutation in mutation_info:
			# first remove mutations from the 'MT' (mitochondrial) chromosome
			chromosome = mutation[self.header_dict['CHROMOSOME']]
			if chromosome != 'MT': 
				# don't want mitochondrial chromosome
				'''
				from wikipedia page on mitochondrial DNA:
				Unlike nuclear DNA, which is inherited from both parents and in which genes are rearranged 
				in the process of recombination, there is usually no change in mtDNA from parent to offspring. 
				'''
				mutation_types = mutation[self.header_dict['MUTATION_DESCRIPTION']].split(',')
				# to save on runtime, first check that mutation is either missense, synonymous, nonsense
				mutation_types_to_keep = ['missense_variant', 'synonymous_variant', 'stop_gained']
				if list(set(mutation_types) & set(mutation_types_to_keep)) != []:
					gene_name = mutation[self.header_dict['GENE_SYMBOL']]
					ensembl_protein = mutation[self.header_dict['HGVSP']].split(':')[0]
					ensembl_mrna = mutation[self.header_dict['HGVSC']].split(':')[0]
					mutation_aa = mutation[self.header_dict['MUTATION_AA']]
					mutation_cds = mutation[self.header_dict['MUTATION_CDS']]
					# print(gene_name, ensembl_protein, ensembl_mrna, mutation_aa, mutation_cds)
					new_mutation_columns = self.get_info_for_new_mutation_info_columns(mutation_aa, mutation_cds, ensembl_protein, gene_name)
					if new_mutation_columns != []:  
						# print(new_mutation_columns)
						cdna_ref_allele, cdna_allele_pos, cdna_alt_allele, ref_res, prot_res_pos, alt_res, uniprot_protein = new_mutation_columns
						mutation.extend([ensembl_mrna, cdna_ref_allele, cdna_allele_pos, cdna_alt_allele, ensembl_protein, ref_res, prot_res_pos, alt_res, uniprot_protein])
						# check num columns
						if len(mutation) != len(column_names_list):
							print('Oh no something went wrong...the number of columns do not match')
						else:
							if 'missense_variant' in mutation_types: 
								missense_info.append(mutation)
							elif 'synonymous_variant' in mutation_types:
								synonymous_info.append(mutation)
							elif 'stop_gained' in mutation_types: # nonsense
								nonsense_info.append(mutation)
				else: # none of the above
					continue

		write_mutation_info_to_file(missense_info, self.processed_missense_mutation_file)
		write_mutation_info_to_file(synonymous_info, self.processed_synonymous_mutation_file)
		write_mutation_info_to_file(nonsense_info, self.processed_nonsense_mutation_file)

def main():
	script_dir = osp.dirname(__file__)
	
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--split_number') # e.g. a number (0-62)
	args = parser.parse_args()

	original_data_dir = osp.join(script_dir, '..', 'data', 'original', 'cosmic_genome_screen_mutants_split_files')
	processed_data_dir = osp.join(script_dir, '..', 'data', 'processed')
	c = COSMIC(original_data_dir, processed_data_dir, args.split_number)
	c.process_genome_screen_mutations()


if __name__=='__main__':
	main()