'''
Processes COSMIC CGC (Cancer Gene Census) mutations
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os.path as osp
from simple_tools import get_mutation_info, pickle_load, write_mutation_info_to_file


class COSMIC:
	def __init__(self, processed_data_dir, mutation_file, mutation_name):
		self.processed_data_dir = processed_data_dir
		self.ensembl_prot_gene_name_uniprot_dict = pickle_load(osp.join(self.processed_data_dir, 'ensembl_prot_gene_name_uniprot_dict.pickle'))
		self.mutation_file = mutation_file
		self.processed_missense_mutation_file = osp.join(self.processed_data_dir, 'mutations', 'cosmic_' + mutation_name + '_missense.tsv')
		self.processed_synonymous_mutation_file = osp.join(self.processed_data_dir, 'mutations', 'cosmic_' + mutation_name + '_synonymous.tsv')
		self.processed_nonsense_mutation_file = osp.join(self.processed_data_dir, 'mutations', 'cosmic_' + mutation_name + '_nonsense.tsv')

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
		if mutation_aa != 'p.?': # remove if having missing mutation_aa info
			if '>' in mutation_cds and 'fs' not in mutation_aa:
				cdna_ref_allele, cdna_allele_pos, cdna_alt_allele = self.get_cdna_info(mutation_cds)
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
		else:
			return []

	def determine_mutation_type(self, ref_res, alt_res): # only for the drug resistance mutations (which don't have a 'MUTATION_DESCRIPTION' column)
		if ref_res != alt_res: # missense or nonsense
			if alt_res == '*': # nonsense
				return 'nonsense'
			else: # missense
				return 'missense'
		else: # synonymous
			return 'synonymous'


	def process_census_gene_mutations(self):
		header_dict, mutation_info = get_mutation_info(self.mutation_file)
		mutation_info[0].extend(['ensembl_mrna', 'cdna_ref_allele', 'cdna_allele_pos', 'cdna_alt_allele', 'ensembl_protein', 'ref_res', 'prot_res_pos', 'alt_res', 'uniprot_protein'])

		missense_info = [mutation_info[0]]
		synonymous_info = [mutation_info[0]]
		nonsense_info = [mutation_info[0]]

		for mutation in mutation_info[1:]:
			# first remove mutations from the 'MT' (mitochondrial) chromosome
			chromosome = mutation[header_dict['CHROMOSOME']]
			if chromosome != 'MT': # quick check shows that both files don't actually have any mutations on the 'MT' chromosome
				# keep missense variants only
				mutation_types = mutation[header_dict['MUTATION_DESCRIPTION']].split(',')
				# to save on runtime, first check that mutation is either missense, synonymous, nonsense
				mutation_types_to_keep = ['missense_variant', 'synonymous_variant', 'stop_gained']

				if list(set(mutation_types) & set(mutation_types_to_keep)) != []:
					gene_name = mutation[header_dict['GENE_SYMBOL']]
					ensembl_protein = mutation[header_dict['HGVSP']].split(':')[0]
					ensembl_mrna = mutation[header_dict['HGVSC']].split(':')[0]
					mutation_aa = mutation[header_dict['MUTATION_AA']]
					mutation_cds = mutation[header_dict['MUTATION_CDS']]
					# print(gene_name, ensembl_protein, ensembl_mrna, mutation_aa, mutation_cds)
					new_mutation_columns = self.get_info_for_new_mutation_info_columns(mutation_aa, mutation_cds, ensembl_protein, gene_name)
					if new_mutation_columns != []:  
						# print(new_mutation_columns)
						cdna_ref_allele, cdna_allele_pos, cdna_alt_allele, ref_res, prot_res_pos, alt_res, uniprot_protein = new_mutation_columns
						mutation.extend([ensembl_mrna, cdna_ref_allele, cdna_allele_pos, cdna_alt_allele, ensembl_protein, ref_res, prot_res_pos, alt_res, uniprot_protein])
						# check num columns
						if len(mutation) != len(mutation_info[0]):
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

	def process_drug_resistance_mutations(self): 
		header_dict, mutation_info = get_mutation_info(self.mutation_file)
		mutation_info[0].extend(['mutation_type', 'ensembl_mrna', 'cdna_ref_allele', 'cdna_allele_pos', 'cdna_alt_allele', 'ensembl_protein', 'ref_res', 'prot_res_pos', 'alt_res', 'uniprot_protein'])

		missense_info = [mutation_info[0]]
		synonymous_info = [mutation_info[0]]
		nonsense_info = [mutation_info[0]]

		for mutation in mutation_info[1:]:
			# first remove mutations from the 'MT' (mitochondrial) chromosome
			chromosome = mutation[header_dict['CHROMOSOME']]
			if chromosome != 'MT':
				gene_name = mutation[header_dict['GENE_SYMBOL']]
				ensembl_protein = mutation[header_dict['HGVSP']].split(':')[0]
				ensembl_mrna = mutation[header_dict['HGVSC']].split(':')[0]
				mutation_aa = mutation[header_dict['MUTATION_AA']]
				mutation_cds = mutation[header_dict['MUTATION_CDS']]
				# print(gene_name, ensembl_protein, ensembl_mrna, mutation_aa, mutation_cds)
				new_mutation_columns = self.get_info_for_new_mutation_info_columns(mutation_aa, mutation_cds, ensembl_protein, gene_name)
				if new_mutation_columns != []:  
					# print(new_mutation_columns)
					cdna_ref_allele, cdna_allele_pos, cdna_alt_allele, ref_res, prot_res_pos, alt_res, uniprot_protein = new_mutation_columns
					mutation_type = self.determine_mutation_type(ref_res, alt_res)
					mutation.extend([mutation_type, ensembl_mrna, cdna_ref_allele, cdna_allele_pos, cdna_alt_allele, ensembl_protein, ref_res, prot_res_pos, alt_res, uniprot_protein])
					# check num columns
					if len(mutation) != len(mutation_info[0]):
						print('Oh no something went wrong...the number of columns do not match')
					else:
						if mutation_type == 'missense': 
							missense_info.append(mutation)
						elif mutation_type == 'synonymous':
							synonymous_info.append(mutation)
						elif mutation_type == 'nonsense':
							nonsense_info.append(mutation)

		write_mutation_info_to_file(missense_info, self.processed_missense_mutation_file)
		write_mutation_info_to_file(synonymous_info, self.processed_synonymous_mutation_file)
		write_mutation_info_to_file(nonsense_info, self.processed_nonsense_mutation_file)

def main():
	script_dir = osp.dirname(__file__)
	original_data_dir = osp.join(script_dir, '..', 'data', 'original')
	processed_data_dir = osp.join(script_dir, '..', 'data', 'processed')

	cancer_gene_census_mutation_file = osp.join(original_data_dir, 'Cosmic_MutantCensus_v100_GRCh38.tsv')
	
	c = COSMIC(processed_data_dir, cancer_gene_census_mutation_file, 'cgc')
	c.process_census_gene_mutations()

if __name__=='__main__':
	main()