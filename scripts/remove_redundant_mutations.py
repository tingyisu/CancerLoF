'''
Removes redundant mutations
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

# keeps only one mutation with the same amino acid change (even if diff nucleotide change) at a given position on a UniProt protein 
# and for dbSNP mutations specifically, which includes the same mutation on diff mRNA/protein transcripts, keeps only one (with the longest protein transcript)

from simple_tools import get_mutation_info, write_mutation_info_to_file, pickle_load
import os.path as osp
import os

class RemoveRedundant:
	def __init__(self, data_dir, final_data_dir, refseq_prot_seq_dict_pickle_file, ensembl_prot_seq_dict_pickle_file, mutation_files_list):
		self.data_dir = data_dir
		self.final_data_dir = final_data_dir
		self.refseq_prot_seq_dict = pickle_load(refseq_prot_seq_dict_pickle_file)
		self.ensembl_prot_seq_dict = pickle_load(ensembl_prot_seq_dict_pickle_file)
		if not osp.exists(self.final_data_dir):
			os.mkdir(self.final_data_dir)
		self.mutation_files_list = mutation_files_list

	# for dbSNP and COSMIC mutations only
	def find_longest_protein_transcript(self, prot_seq_dict, proteins_list):
		protein_lengths = [len(prot_seq_dict[protein]) for protein in proteins_list]
		# print(protein_lengths)
		return protein_lengths.index(max(protein_lengths)) # if more than one max, returns the first occurrence
	
	# for dbSNP and COSMIC mutations only
	def select_dbsnp_cosmic_mutation(self, type, mutation_file):
		# type = 'refseq' or 'ensembl'
		header_dict, mutations = get_mutation_info(mutation_file)
		print('# of mapped mutations:', len(mutations)-1)

		# get only one mutation per (mutation_id, chrom_ref_allele, chrom_alt_allele)
		# for dbsnp, mutation_id = rs number ('ID' in header field)
		# for cosmic, mutation_id = GENOMIC_MUTATION_ID (starts with 'COSV')
		mutations_dict = {} # key = (mutation_id, chrom_ref_allele, chrom_alt_allele), value = list of mutations
		mutation_identifier = ''
		ref_allele_identifier, alt_allele_identifier = '', ''
		if 'dbsnp' in mutation_file:
			mutation_identifier = 'ID'
			ref_allele_identifier = 'chrom_ref_allele'
			alt_allele_identifier = 'chrom_alt_allele'
		elif 'cosmic' in mutation_file:
			mutation_identifier = 'GENOMIC_MUTATION_ID'
			ref_allele_identifier = 'GENOMIC_WT_ALLELE'
			alt_allele_identifier = 'GENOMIC_MUT_ALLELE'
		else:
			print('ERROR! Should only be selecting one mutation per (mutation_id, chrom_ref_allele, chrom_alt_allele) for the dbSNP and COSMIC mutations!')
		for mutation in mutations[1:]:
			mutation_id = mutation[header_dict[mutation_identifier]]
			chrom_ref_allele = mutation[header_dict[ref_allele_identifier]]
			chrom_alt_allele = mutation[header_dict[alt_allele_identifier]]
			if (mutation_id, chrom_ref_allele, chrom_alt_allele) not in mutations_dict:
				mutations_dict[(mutation_id, chrom_ref_allele, chrom_alt_allele)] = [mutation]
			else:
				mutations_dict[(mutation_id, chrom_ref_allele, chrom_alt_allele)].append(mutation)

		# if 'cosmic' in mutation_file:
		# 	print(mutations_dict.keys())

		# sanity check
		# print(len(mutations_dict.keys()))

		protein_type = type + '_protein'
		prot_seq_dict = {}
		if type == 'refseq':
			prot_seq_dict = self.refseq_prot_seq_dict
		elif type == 'ensembl':
			prot_seq_dict = self.ensembl_prot_seq_dict
		else:
			print('ERROR! Proteins can only come from RefSeq or Ensembl sequence databases!')
		# for the same mutation on diff mRNA/protein transcripts, keeps only one (with the longest protein transcript)
		selected_mutations = [mutations[0]]
		for key in mutations_dict:
			proteins_list = [unambiguous_mutation[header_dict[protein_type]] for unambiguous_mutation in mutations_dict[key]]
			unique_proteins_list = list(set(proteins_list))
			selected_mutation = mutations_dict[key][self.find_longest_protein_transcript(prot_seq_dict, unique_proteins_list)]
			selected_mutations.append(selected_mutation)

		return header_dict, selected_mutations

	def remove_redundant_mutations(self, header_dict, mutations):
		# only keep one mutation per (uniprot_protein, prot_res_pos_in_uniprot_protein, ref_res, alt_res)
		mutations_dict = {} # key = (uniprot_protein, prot_res_pos_in_uniprot_protein, ref_res, alt_res), value = list of mutations
		for mutation in mutations[1:]:
			uniprot_protein = mutation[header_dict['uniprot_protein']]
			prot_res_pos_in_uniprot_protein = mutation[header_dict['prot_res_pos_in_uniprot_protein']]
			ref_res = mutation[header_dict['ref_res']]
			alt_res = mutation[header_dict['alt_res']]
			if (uniprot_protein, prot_res_pos_in_uniprot_protein, ref_res, alt_res) not in mutations_dict:
				mutations_dict[(uniprot_protein, prot_res_pos_in_uniprot_protein, ref_res, alt_res)] = [mutation]
			else:
				mutations_dict[(uniprot_protein, prot_res_pos_in_uniprot_protein, ref_res, alt_res)].append(mutation)

		# take first instance 
		selected_mutations = [mutations[0]]
		for key in mutations_dict:
			selected_mutations.append(mutations_dict[key][0]) # take first one

		return selected_mutations


	def remove_redundant_mutations_all(self):
		for mutation_file in self.mutation_files_list:
			mutation_file_path = osp.join(self.data_dir, mutation_file)
			nonredundant_mutations = []
			if 'dbsnp' in mutation_file:
				print('-----', mutation_file, '-----')
				header_dict, dbsnp_selected_mutations = self.select_dbsnp_cosmic_mutation('refseq', mutation_file_path)
				nonredundant_mutations = self.remove_redundant_mutations(header_dict, dbsnp_selected_mutations)
				print('# of selected mutations (one per rs, chrom_ref_allele, chrom_alt_allele):', len(dbsnp_selected_mutations)-1)
				print('# of nonredundant mutations:', len(nonredundant_mutations)-1)

				write_mutation_info_to_file(nonredundant_mutations, osp.join(self.final_data_dir, mutation_file[:-4] + '_nonredundant.tsv'))
			elif 'cosmic' in mutation_file:
				print('-----', mutation_file, '-----')
				header_dict, cosmic_selected_mutations = self.select_dbsnp_cosmic_mutation('ensembl', mutation_file_path)
				nonredundant_mutations = self.remove_redundant_mutations(header_dict, cosmic_selected_mutations)
				print('# of selected mutations (one per GENOMIC_MUTATION_ID, GENOMIC_WT_ALLELE, GENOMIC_MUT_ALLELE):', len(cosmic_selected_mutations)-1)
				print('# of nonredundant mutations:', len(nonredundant_mutations)-1)

				write_mutation_info_to_file(nonredundant_mutations, osp.join(self.final_data_dir, mutation_file[:-4] + '_nonredundant.tsv'))
			else:
				print('-----', mutation_file, '-----')
				header_dict, mutations = get_mutation_info(mutation_file_path)
				nonredundant_mutations = self.remove_redundant_mutations(header_dict, mutations)
				print('# of mapped mutations:', len(mutations)-1)
				print('# of nonredundant mutations:', len(nonredundant_mutations)-1)

				write_mutation_info_to_file(nonredundant_mutations, osp.join(self.final_data_dir, mutation_file[:-4] + '_nonredundant.tsv'))

def main():
	script_dir = osp.dirname(__file__)
	processed_data_dir = osp.join(script_dir, '..', 'data', 'processed')
	data_dir = osp.join(processed_data_dir, 'mutations')
	final_data_dir = osp.join(script_dir, '..', 'data', 'processed', 'mutations_final')
	refseq_prot_seq_dict_pickle_file = osp.join(processed_data_dir, 'refseq_prot_seq_dict.pickle')
	ensembl_prot_seq_dict_pickle_file = osp.join(processed_data_dir, 'ensembl_prot_seq_dict.pickle')

	interactomes = ['hiunion', 'intact']
	mutation_data = ['dbsnp', 'clinvar', 'cosmic']
	cosmic_datasets = ['genome_screen', 'cgc',]
	cosmic_mutation_types = ['missense', 'nonsense'] # no nonstop...
	mutation_types = ['missense', 'nonstop', 'nonsense']
	mutation_files_list = []

	for interactome in interactomes:
		for data in mutation_data:
			if data == 'clinvar' or data == 'dbsnp':
				for mutation_type in mutation_types:
					mutation_files_list.append('_'.join([interactome, 'mapped', data, mutation_type + '.tsv']))
			elif data == 'cosmic':
				for cosmic_dataset in cosmic_datasets:
					for cosmic_mutation_type in cosmic_mutation_types:
						mutation_files_list.append('_'.join([interactome, 'mapped', data, cosmic_dataset, cosmic_mutation_type + '.tsv']))
			else: 
				print('Something went wrong!')


	r = RemoveRedundant(data_dir, final_data_dir, refseq_prot_seq_dict_pickle_file, ensembl_prot_seq_dict_pickle_file, mutation_files_list)
	r.remove_redundant_mutations_all()


if __name__=='__main__':
	main()
