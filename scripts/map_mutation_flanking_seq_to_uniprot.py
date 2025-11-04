'''
Maps mutations onto UniProt sequences in structural interactomes
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os.path as osp
import pickle
from Bio.Align import substitution_matrices 
import math
from Bio import SeqIO
from simple_tools import get_mutation_info, pickle_load
import re

'''
This class includes functions to find mutation flanking sequences (the first 10 and last 10 amino acids flanking a mutation or shorter if mutation is near front or end of protein) 
to map the mutations to UniProt proteins in that could exist in either the HI-union or IntAct structural interactomes
'''
class MapMutation:
	def __init__(self, mutation_name, mutation_nomenclature, mutation_data_dir, mutation_file, prot_seq_dict_file, blast_hiunion_fasta_file, blast_intact_fasta_file, hiunion_uniprot_proteins_file, intact_uniprot_proteins_file):
		self.mutation_name = mutation_name # clinvar, dbsnp, cosmic
		self.mutation_nomenclature = mutation_nomenclature # either 'refseq' or 'ensembl'
		self.mutation_data_dir = mutation_data_dir
		self.header_dict, self.mutation_info = get_mutation_info(mutation_file) 
		'''
		self.header_dict --> key = column name, value = index in list of column names
		self.mutation_info --> is a list of lines from the initial mutation file
		'''
		self.prot_seq_dict = pickle_load(prot_seq_dict_file) # key = Ensembl/RefSeq protein ID, value = protein sequence
		self.blast_hiunion_fasta_file = blast_hiunion_fasta_file
		self.blast_intact_fasta_file = blast_intact_fasta_file
		self.hiunion_uniprot_proteins = pickle_load(hiunion_uniprot_proteins_file) # list of UniProt proteins that may be in final HI-union structural interactome
		self.intact_uniprot_proteins = pickle_load(intact_uniprot_proteins_file) # list of UniProt proteins that may be in final IntAct structural interactome
		self.swissprot_seq_hiunion = {}
		self.swissprot_seq_intact = {}
		# self.hiunion_mapped_mutations = []
		# self.intact_mapped_mutations = []
		# self.mismatched = [] # list of mutations that are mismatched
	
		
	# update self.mutation_info with BLOSUM62 scores
	def get_blosum_score(self):
		print('Adding BLOSUM62 scores...')
		# b62 = MatrixInfo.blosum62
		b62 = substitution_matrices.load('BLOSUM62') # the values in b62 are floats, not integers (e.g. -2.0, instead of -2)
		self.header_dict['blosum62_score'] = len(self.header_dict)
		self.mutation_info[0].append('blosum62_score')
		for i in range(1, len(self.mutation_info)):
			ref_res = self.mutation_info[i][self.header_dict['ref_res']]
			alt_res = self.mutation_info[i][self.header_dict['alt_res']]
			score = math.inf
			if (ref_res, alt_res) in b62:
				score = b62[(ref_res, alt_res)]
			else:
				score = b62[(alt_res, ref_res)]
			if score == math.inf:
				print('Oh no...cannot get Blosum62 score for:', mutation_info[i])
			else:
				self.mutation_info[i].append(score)


	# update self.mutation_info with mutation flanking sequences
	def get_flanking_sequence(self, side_length): # mutation nomenclature can either be 'refseq' or 'ensembl'
		'''
		side_length = number of residues to include on each side of the mutation residue
		'''
		print('Adding mutation flanking sequences...')
		# add to self.header_dict
		self.header_dict['flanking_seq'] = len(self.header_dict)
		self.header_dict['prot_res_pos_in_flanking_seq'] = len(self.header_dict)

		# add flanking_seq and prot_res_pos_in_flanking_seq
		self.mutation_info[0].append('flanking_seq')
		self.mutation_info[0].append('prot_res_pos_in_flanking_seq')

		num_diff_ref_res = 0

		num_no_flanking_seq = 0

		for i in range(1, len(self.mutation_info)):
			prot = self.mutation_info[i][self.header_dict[self.mutation_nomenclature + '_protein']]
			prot_res_pos = int(self.mutation_info[i][self.header_dict['prot_res_pos']])
			ref_res = self.mutation_info[i][self.header_dict['ref_res']]
			# check that ref_res at prot_res_pos matches with corresponding res in self.prot_seq_dict[prot]
			# print(prot, prot_res_pos, len(self.prot_seq_dict[prot]), self.mutation_info[i][self.header_dict['#AlleleID']])
			prot_seq = self.prot_seq_dict[prot]
			# print(prot_res_pos, len(prot_seq))
			# print(prot_res_pos, len(prot_seq))
			# print(self.mutation_info[i])
			ref_res_in_prot_seq = ''
			diff_ref_res = False
			if 'nonstop' not in self.mutation_name:
				# some mutation positions are not within the length of the protein for some reason...(maybe due to insertions/deletions)?
				if prot_res_pos > len(prot_seq):
					num_diff_ref_res += 1
					diff_ref_res = True
				else:
					ref_res_in_prot_seq = prot_seq[prot_res_pos-1]
					if ref_res != ref_res_in_prot_seq: # sanity check
						# print('Different wildtype residue!', prot, ref_res, ref_res_in_prot_seq)
						num_diff_ref_res += 1
						diff_ref_res = True
			else:
				ref_res_in_prot_seq = '*'
			# make sure ref res is the same in order to make flanking sequence
			if not diff_ref_res:
				# find flanking sequence
				front_flank, back_flank = '', ''
				# initialize indices for getting flanking sequence
				front_index, back_index = prot_res_pos-1, prot_res_pos-1 # python index
				# iterate to get front flanking sequence
				while front_index >= 1 and front_index >= (prot_res_pos - side_length):
					front_index -= 1
					front_flank = prot_seq[front_index] + front_flank
				# iterate to get back flanking sequence
				while back_index < (len(prot_seq) - 1) and back_index < (prot_res_pos + side_length - 1):
					back_index += 1
					back_flank += prot_seq[back_index]
				flanking_seq = front_flank + ref_res_in_prot_seq + back_flank
				prot_res_pos_in_flanking_seq = len(front_flank) + 1 # index starting from 1 (to be consistent with prot_res_pos)
				self.mutation_info[i].append(flanking_seq)
				self.mutation_info[i].append(prot_res_pos_in_flanking_seq)
			else:
				self.mutation_info[i].append('') # no flanking seq b/c ref res doesn't even match to what's on the RefSeq/Ensembl protein seq
				self.mutation_info[i].append(-1) # no position in flanking seq b/c ref res doesn't even match
				num_no_flanking_seq += 1


		print('Number of mutations with different reference residues:', num_diff_ref_res, num_no_flanking_seq)


	# get only that swissprot proteins & sequences that exist in inital file for blast
	# use to update self.swissprot_seq_hiunion & self.swissprot_seq_intact
	def get_swissprot_seqs_for_blast(self, blast_fasta, name):
		print('Getting SwissProt-reviewed UniProt sequences in', name + '...')
		swissprot_seq = {}
		with open(blast_fasta, 'r') as f:
			fasta_seqs = SeqIO.parse(f,'fasta')
			for fasta_seq in fasta_seqs:
				swissprot_seq[fasta_seq.id] = str(fasta_seq.seq) # convert to str just in case
		return swissprot_seq

	# finds number of mismatching amino acids for mutations that do not map 100%
	def find_num_mismatching_aas(self, flanking_seq, uniprot_seq, prot_res_pos):
		max_num_mismatches = len(flanking_seq)
		# find number of mismatches
		for i in range(1, max_num_mismatches+1):
			substring = '(' + flanking_seq + ')' + '{e<=' + str(i) + '}'
			# print(substring)
			m = regex.findall(substring, uniprot_seq)
			print(substring, m)
			if m:
				print(flanking_seq, uniprot_seq, prot_res_pos, m)
				break

	# match mutation flanking sequences to corresponding area on swissprot protein sequence
	# check that the wildtype residue has the same index on the swissprot protein
	def match_flanking_sequence(self, mapped_mutation_info_file, swissprot_seq, uniprot_proteins_list):
		print('Matching mutation flanking sequences with corresponding uniprot (swissprot reviewed) protein...')
		num_mismatch_pos, num_mismatch = 0, 0
		mismatch_pos_list = []
		num_mapped = 0
		num_diff_ref_res = 0
		with open(mapped_mutation_info_file, 'w') as f:
			# write header line
			# add in prot_res_pos_in_uniprot_protein
			header = self.mutation_info[0] + ['prot_res_pos_in_uniprot_protein']
			f.write('\t'.join(header) + '\n')
			for i in range(1, len(self.mutation_info)):

				uniprot_protein = self.mutation_info[i][self.header_dict['uniprot_protein']]

				# only keep if uniprot_protein is in swissprot_seq (from blast fasta file) and in uniprot_proteins_list (in interactome)
				if uniprot_protein in swissprot_seq and uniprot_protein in uniprot_proteins_list:
					seq = swissprot_seq[uniprot_protein]
					flanking_seq = self.mutation_info[i][self.header_dict['flanking_seq']]
					prot_res_pos_in_flanking_seq = int(self.mutation_info[i][self.header_dict['prot_res_pos_in_flanking_seq']])
					prot_res_pos = int(self.mutation_info[i][self.header_dict['prot_res_pos']]) # in regards to RefSeq

					if flanking_seq == '' and prot_res_pos_in_flanking_seq == -1: # remove mutations with no mismatching ref res in RefSeq/Ensembl protein seq
						num_diff_ref_res += 1
					else:	
						match_indices = [m.start() for m in re.finditer(flanking_seq, seq)]
						if len(match_indices) > 1:
							for index in match_indices:
								if index == prot_res_pos:
									print('FOUND INDEX in multiple matches',  match_indices, flanking_seq, uniprot_protein)
						elif len(match_indices) == 0:
							num_mismatch += 1
						else:
							if len(match_indices) != 1: # sanity check
								print('Uh oh......')
							else:
								to_write = [str(item) for item in self.mutation_info[i]]
								prot_res_pos_in_uniprot_protein = match_indices[0] + prot_res_pos_in_flanking_seq

								if prot_res_pos != prot_res_pos_in_uniprot_protein: # in Pythonic indexing (?????)
									num_mismatch_pos += 1

									mismatch_pos_list.append(abs(prot_res_pos_in_uniprot_protein - prot_res_pos)) 

								# add to self.mapped_mutation_info
								to_write += [str(prot_res_pos_in_uniprot_protein)]
								f.write('\t'.join(to_write) + '\n')
								num_mapped += 1


		print('Number of mutations with mismatching ref res on RefSeq/Ensembl protein:', num_diff_ref_res)
		print('Number of mismatching mutations:', num_mismatch)
		print('Number of matching mutations with different positions on corresponding uniprot proteins (kept):', num_mismatch_pos)
		if len(mismatch_pos_list) > 0:
			print('Average difference in mutation position on RefSeq vs. uniprot protein:', sum(mismatch_pos_list)/len(mismatch_pos_list))
		print('Number of mapped mutations:', num_mapped)


	def map_mutation_flanking_seq_to_uniprot_all(self):
		print('.....' + self.mutation_name + '.....')
		# update self.missense_info first with BLOSUM62 scores and mutation flanking sequences
		if 'missense' in self.mutation_name: # only get BLOSUM62 scores for missense mutations
			self.get_blosum_score()
		self.get_flanking_sequence(10)

		# map mutation flanking sequences to uniprot proteins that are potentially in HI-union SI
		print('-----HI-union-----')
		self.swissprot_seq_hiunion = self.get_swissprot_seqs_for_blast(self.blast_hiunion_fasta_file, 'HI-union')
		hiunion_mapped_missense_info_file = osp.join(self.mutation_data_dir, '_'.join(['hiunion', 'mapped', self.mutation_name + '.tsv']))
		self.match_flanking_sequence(hiunion_mapped_missense_info_file, self.swissprot_seq_hiunion, self.hiunion_uniprot_proteins)

		print('-----IntAct-----')
		# map mutation flanking sequences to HI-union & IntAct structural interactome separately and write to seperate files
		self.swissprot_seq_intact = self.get_swissprot_seqs_for_blast(self.blast_intact_fasta_file, 'IntAct')
		intact_mapped_missense_info_file = osp.join(self.mutation_data_dir, '_'.join(['intact', 'mapped', self.mutation_name + '.tsv']))
		self.match_flanking_sequence(intact_mapped_missense_info_file, self.swissprot_seq_intact, self.intact_uniprot_proteins)



def main():
	script_dir = osp.dirname(__file__)
	processed_data_dir = osp.join(script_dir, '..', 'data', 'processed')
	mutation_data_dir = osp.join(processed_data_dir, 'mutations')
	data_interactome_dir = osp.join(processed_data_dir, 'interactome')

	# prot seq dicts
	refseq_prot_seq_dict_file = osp.join(processed_data_dir, 'refseq_prot_seq_dict.pickle')
	ensembl_prot_seq_dict_file = osp.join(processed_data_dir, 'ensembl_prot_seq_dict.pickle')

	# HI-union/IntAct UniProt fasta sequences
	blast_hiunion_fasta_file = osp.join(data_interactome_dir, 'hiunion_uniprot_sequences_blast.fasta')
	blast_intact_fasta_file = osp.join(data_interactome_dir, 'intact_uniprot_sequences_blast.fasta')

	# HI-union/IntAct UniProt proteins
	hiunion_uniprot_proteins_file = osp.join(mutation_data_dir, 'hiunion_uniprot_proteins_to_map_mutations_to.pickle')
	intact_uniprot_proteins_file = osp.join(mutation_data_dir, 'intact_uniprot_proteins_to_map_mutations_to.pickle')

	# lists to loop through
	# no synonymous mutation for now
	# no 'pathogenic' dbSNP mutations for now...only dealing with 'non-pathogenic' ones
	mutation_types = ['missense', 'nonstop', 'nonsense']
	mutation_datasets = ['dbsnp', 'clinvar']
	cosmic_mutation_types = ['missense', 'nonsense'] # no nonstop...
	cosmic_datasets = ['genome_screen', 'cgc']

	print('**********ClinVar & dbSNP**********')
	for mutation_dataset in mutation_datasets:
		for mutation_type in mutation_types:
			mutation_name = '_'.join([mutation_dataset, mutation_type])
			mutation_file = osp.join(mutation_data_dir, mutation_name + '.tsv')
			clinvar = MapMutation(mutation_name, 'refseq', mutation_data_dir, mutation_file, refseq_prot_seq_dict_file, blast_hiunion_fasta_file, blast_intact_fasta_file, hiunion_uniprot_proteins_file, intact_uniprot_proteins_file)
			clinvar.map_mutation_flanking_seq_to_uniprot_all()

	print('**********COSMIC**********')
	for cosmic_dataset in cosmic_datasets:
		for cosmic_mutation_type in cosmic_mutation_types:
			mutation_name = '_'.join(['cosmic', cosmic_dataset, cosmic_mutation_type])
			mutation_file = osp.join(mutation_data_dir, mutation_name + '.tsv')
			cosmic = MapMutation(mutation_name, 'ensembl', mutation_data_dir, mutation_file, ensembl_prot_seq_dict_file, blast_hiunion_fasta_file, blast_intact_fasta_file, hiunion_uniprot_proteins_file, intact_uniprot_proteins_file)
			cosmic.map_mutation_flanking_seq_to_uniprot_all()

	
if __name__=='__main__':
	main()