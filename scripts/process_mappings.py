'''
Processes mappings necessary for processing & mapping mutations to the structural interactomes
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os.path as osp
from Bio import SeqIO
from simple_tools import pickle_dump, get_missense_info, pickle_load, get_list_from_file

'''
pickles the following mappings:
(1) Ensembl protein -> uniprot
(2) RefSeq protein -> uniprot
(3) RefSeq mRNA -> RefSeq protein
(4) Ensembl mRNA -> (Ensembl prot, uniprot) 
(5) Ensembl gene -> uniprot (for HI-Union)
(6) RefSeq protein -> seq
(7) Ensembl prot -> seq
(8) (RefSeq prot, gene_name) -> UniProt
(9) (Ensembl prot, gene_name) -> UniProt


these mappings are necessary for processing & mapping mutations to the structural interactomes
will be finding UniProt proteins based on gene name and RefSeq/Ensembl protein ID

NOTE: only takes canonical UniProt IDs (those without '-*' in their IDs) that are SwissProt reviewed
'''

class Mappings:
	def __init__(self, data_dir, orig_data_dir):
		self.uniprot_mappings_file = osp.join(orig_data_dir, 'HUMAN_9606_idmapping.dat')
		self.lrg_ref_seq_gene_file = osp.join(orig_data_dir, 'LRG_RefSeqGene')
		self.ref_seq_protein_seq_file = osp.join(data_dir, 'merged_ref_seq_protein.fa')
		self.ensembl_protein_seq_file = osp.join(orig_data_dir, 'Homo_sapiens.GRCh38.pep.all.fa')
		self.ensembl_uniprot_file = osp.join(orig_data_dir, 'Homo_sapiens.GRCh38.109.uniprot.tsv')
		self.swissprot_fasta_file = osp.join(orig_data_dir, 'uniprot_reviewed_human_proteome.fasta')
		self.data_dir = data_dir
		# attributes to get
		self.swissprot_ids = get_list_from_file(osp.join(orig_data_dir, 'uniprot_reviewed_human_proteome.list')) # list of reviewed Uniprotkb proteins
		self.swissprot_seq_dict = {} # key = swissprot id, value = sequence
		self.ensembl_prot_uniprot_dict = {} # key = Ensembl protein ID, value = list of UniProt protein ID(s)
		self.ensembl_prot_seq_dict = {} # key = Ensembl protein ID, value = protein sequence
		self.refseq_prot_uniprot_dict = {} # key = RefSeq protein ID, value = list of UniProt protein ID(s)
		self.refseq_prot_seq_dict = {} # key = RefSeq protein ID, value = protein sequence
		self.refseq_mrna_refseq_prot_dict = {} # key = RefSeq mRNA ID, value = RefSeq protein ID
		self.ensembl_mrna_prot_uniprot_dict = {} # key = Ensembl mRNA ID, value = (Ensembl Protein, UniProt protein ID)
		self.ensembl_gene_uniprot_dict = {} # key = Ensembl gene, value = UniProt protein
		self.uniprot_gene_name_dict = {} # key = UniProt protein ID, value = list of Gene names
		self.refseq_prot_gene_name_uniprot_dict = {} # key = (RefSeq protein ID, gene name), value = UniProt protein ID
		self.ensembl_prot_gene_name_uniprot_dict = {} # key = (Ensembl protein ID, gene name), value = UniProt protein ID
		# files to pickle to
		self.swissprot_ids_list_pickle_file = osp.join(self.data_dir, 'swissprot_ids_list.pickle')
		self.swissprot_seq_dict_pickle_file = osp.join(self.data_dir, 'swissprot_seq_dict.pickle')
		self.ensembl_prot_uniprot_dict_pickle_file = osp.join(self.data_dir, 'ensembl_prot_uniprot_dict.pickle')
		self.ensembl_prot_seq_dict_pickle_file = osp.join(self.data_dir, 'ensembl_prot_seq_dict.pickle')
		self.refseq_prot_uniprot_dict_pickle_file = osp.join(self.data_dir, 'refseq_prot_uniprot_dict.pickle')
		self.refseq_prot_seq_dict_pickle_file = osp.join(self.data_dir, 'refseq_prot_seq_dict.pickle')
		self.refseq_mrna_refseq_prot_dict_pickle_file = osp.join(self.data_dir, 'refseq_mrna_refseq_prot_dict.pickle')
		# these two are from Ensembl mappings file
		self.ensembl_mrna_prot_uniprot_dict_pickle_file = osp.join(self.data_dir, 'ensembl_mrna_prot_uniprot_dict.pickle')
		self.ensembl_gene_uniprot_dict_pickle_file = osp.join(self.data_dir, 'ensembl_gene_uniprot_dict.pickle')
		# ---
		self.uniprot_gene_name_dict_pickle_file = osp.join(self.data_dir, 'uniprot_gene_name_dict.pickle')
		self.refseq_prot_gene_name_uniprot_dict_pickle_file = osp.join(self.data_dir, 'refseq_prot_gene_name_uniprot_dict.pickle')
		self.ensembl_prot_gene_name_uniprot_dict_pickle_file = osp.join(self.data_dir, 'ensembl_prot_gene_name_uniprot_dict.pickle')

	# updates self.swissprot_seq_dict
	# checks that self.swissprot_ids are all keys in self.swissprot_seq_dict
	# then pickles self.swissprot_seq_dict and  self.swissprot_ids to files (so that don't need to process again for IntAct) (no need for querying now)
	def get_swissprot_id_seq(self):
		print('------Getting swissprot sequences and pickling the ids list and seqs dict-----')
		with open(self.swissprot_fasta_file, 'r') as f:
			fasta_seqs = SeqIO.parse(f,'fasta')
			for fasta_seq in fasta_seqs:
				# print(fasta_seq.id.split('|')[1])
				swissprot = fasta_seq.id.split('|')[1]
				self.swissprot_seq_dict[swissprot] = str(fasta_seq.seq)

		# check keys in self.swssiprot_seq_dict
		# swissprot_ids = list(self.swissprot_seq_dict.keys())
		# print(sorted(swissprot_ids) == sorted(self.swissprot_ids))

		# dump to files
		pickle_dump(self.swissprot_ids, self.swissprot_ids_list_pickle_file)
		pickle_dump(self.swissprot_seq_dict, self.swissprot_seq_dict_pickle_file)

	# use to update self.ensembl_prot_seq_dict & self.refseq_prot_seq_dict
	# get seq_dicts first, so that can retrieve protein ID -> uniprot mappings for only those protein IDs w/ known protein sequences
	def get_protein_seq(self, name, file):
		prot_id_seq_dict = {}
		print('-----Getting protein sequences for', name,  'protein IDs-----')
		with open(file, 'r') as f:
			fasta_seqs = SeqIO.parse(f,'fasta')
			for fasta_seq in fasta_seqs:
				if fasta_seq.id in prot_id_seq_dict:
					if str(fasta_seq.seq) != prot_id_seq_dict[fasta_seq.id]:
						print('Different protein sequence!!!!', fasta_seq.id)
				else:
					prot_id_seq_dict[fasta_seq.id] = str(fasta_seq.seq)
		return prot_id_seq_dict

	# only works for PanCancer driver mutations b/c ENST transcript does not have a version number
	# while the ones in COSMIC do...
	def get_ensembl_specific_mappings(self):
		header_dict, mappings = get_missense_info(self.ensembl_uniprot_file)
		for mapping in mappings[1:]:
			if mapping[header_dict['db_name']] == 'Uniprot/SWISSPROT': # and mapping[header_dict['info_type']] == 'DIRECT'
				gene = mapping[header_dict['gene_stable_id']]
				mrna = mapping[header_dict['transcript_stable_id']]
				protein = mapping[header_dict['protein_stable_id']]
				uniprot = mapping[header_dict['xref']]
				if uniprot in self.swissprot_ids: # sanity check

					# update self.ensembl_mrna_prot_uniprot_dict
					if mrna not in self.ensembl_mrna_prot_uniprot_dict:
						self.ensembl_mrna_prot_uniprot_dict[mrna] = (protein, uniprot)
					else: # if multiple entries, change if find a longer UniProt protein
						prev_entry = self.ensembl_mrna_prot_uniprot_dict[mrna]
						curr_entry = (protein, uniprot)
						if prev_entry != curr_entry:
							if len(self.swissprot_seq_dict[uniprot]) > len(self.swissprot_seq_dict[prev_entry[1]]):
								self.ensembl_mrna_prot_uniprot_dict[mrna] = curr_entry
								# print('updated:', mrna, 'curr_entry:', curr_entry, 'prev_entry:', prev_entry, len(self.swissprot_seq_dict[uniprot]), len(self.swissprot_seq_dict[prev_entry[1]]))
							# print(mrna, 'already present!', 'prev_entry:', prev_entry, 'curr_entry:', curr_entry)

					# update self.ensembl_gene_uniprot_dict
					if gene not in self.ensembl_gene_uniprot_dict:
						self.ensembl_gene_uniprot_dict[gene] = uniprot
					else: # if multiple entries, change if find a longer UniProt protein
						prev_uniprot = self.ensembl_gene_uniprot_dict[gene]
						if prev_uniprot != uniprot:
							if len(self.swissprot_seq_dict[uniprot]) > len(self.swissprot_seq_dict[prev_uniprot]):
								self.ensembl_gene_uniprot_dict[gene] = uniprot
								# print('updated:', gene, 'curr_uniprot:', uniprot, 'prev_uniprot:', prev_uniprot, len(self.swissprot_seq_dict[uniprot]), len(self.swissprot_seq_dict[prev_uniprot]))
							# print(gene, 'already present!', 'prev_uniprot:', prev_uniprot, 'curr_uniprot:', uniprot)

				# else:
				# 	print(uniprot, 'not SWISSPROT reviewed!')

		print('Number of Ensembl mRNA to Ensembl Protein/UniProt mappings:', len(self.ensembl_mrna_prot_uniprot_dict))
		pickle_dump(self.ensembl_mrna_prot_uniprot_dict, self.ensembl_mrna_prot_uniprot_dict_pickle_file)

		print('Number of Ensembl gene/UniProt mappings:', len(self.ensembl_gene_uniprot_dict))
		pickle_dump(self.ensembl_gene_uniprot_dict, self.ensembl_gene_uniprot_dict_pickle_file)


	# update self.refseq_mrna_refseq_prot_dict
	def get_refseq_mrna_refseq_prot(self):
		print('-----Getting RefSeq mRNA to protein ID mappings-----')
		header_dict, refseq_info = get_missense_info(self.lrg_ref_seq_gene_file)
		# print(len(refseq_info))

		for i in range(1, len(refseq_info)):
			mrna_acc = refseq_info[i][header_dict['RNA']]
			prot_acc = refseq_info[i][header_dict['Protein']] 
			# gene_name = items[header_dict['Symbol']]
			if prot_acc in self.refseq_prot_seq_dict: # take only the entries that have protein sequences
				# check that duplicates contain the same info (same prot_acc and gene_name)
				if mrna_acc in self.refseq_mrna_refseq_prot_dict:
					# print('Have already seen this mRNA accession!', line)
					if prot_acc != self.refseq_mrna_refseq_prot_dict[mrna_acc]:
						print('Different protein accession:', mrna_acc, prot_acc, self.refseq_mrna_refseq_prot_dict[mrna_acc])
				else:
					self.refseq_mrna_refseq_prot_dict[mrna_acc] = prot_acc

	# update self.refseq_prot_uniprot_dict and self.ensembl_prot_uniprot_dict
	def get_uniprot_mappings(self):
		print('-----Getting UniProt mappings for Ensembl & RefSeq protein IDs-----')
		# i = 0
		with open(self.uniprot_mappings_file, 'r') as f:
			for line in f:
				items = line.strip().split('\t')
				uniprot_protein = items[0]
				category = items[1]
				other_id = items[2]

				# i += 1
				# print('Reading line:', i)

				if '-' not in uniprot_protein and uniprot_protein in self.swissprot_ids: # take only canonical uniprotkb ID that is swissprot reviewed
					if category == 'RefSeq':
						if other_id in self.refseq_prot_seq_dict: # make sure that refseq_prot has an existing protein sequence (i.e. exists in self.refseq_prot_seq_dict)
							if other_id in self.refseq_prot_uniprot_dict: # already exists, duplicate entry
								self.refseq_prot_uniprot_dict[other_id].append(uniprot_protein)
								# if uniprot_protein != self.refseq_prot_uniprot_dict: # make sure uniprot proteins are the same in duplicate entries
								# 	print('Different uniprot ID for RefSeq prot:', other_id, self.refseq_prot_uniprot_dict[other_id], uniprot_protein)
							else:
								self.refseq_prot_uniprot_dict[other_id] = [uniprot_protein]
					elif category == 'Ensembl_PRO': # make sure that ensembl_prot has an existing protein sequence (i.e. exists in self.ensembl_prot_seq_dict)
						if other_id in self.ensembl_prot_seq_dict:
							if other_id in self.ensembl_prot_uniprot_dict: # already exists, duplicate entry
								self.ensembl_prot_uniprot_dict[other_id].append(uniprot_protein)
								# if uniprot_protein != self.ensembl_prot_uniprot_dict: # make sure uniprot proteins are the same in duplicate entries
								# 		print('Different uniprot ID for Ensembl prot:', other_id, self.ensembl_prot_uniprot_dict[other_id], uniprot_protein)
							else:
								self.ensembl_prot_uniprot_dict[other_id] = [uniprot_protein]
					elif category == 'Gene_Name':
						if uniprot_protein in self.uniprot_gene_name_dict:
							self.uniprot_gene_name_dict[uniprot_protein].append(other_id)
							# if prev_gene_name != other_id:
							# 	print('Different gene name for uniprot protein:', uniprot_protein, prev_gene_name, other_id)
						else:
							self.uniprot_gene_name_dict[uniprot_protein] = [other_id]
					else:
						continue
				else:
					continue

	# for updating self.refseq_prot_gene_name_uniprot_dict & self.ensembl_prot_gene_name_uniprot_dict
	def create_tuple_dict(self, name, prot_uniprot_dict, uniprot_gene_name_dict):
		print('-----Getting tuple dict for', name + '-----')
		tuple_dict = {}
		for prot_id in prot_uniprot_dict:
			for uniprot_protein in prot_uniprot_dict[prot_id]:
				if uniprot_protein in uniprot_gene_name_dict:
					for gene_name in uniprot_gene_name_dict[uniprot_protein]:
						tuple_dict[(prot_id, gene_name)] = uniprot_protein
				else:
					continue
		return tuple_dict


	def test_print_key_with_multiple_values(self, dict_to_print):
		i = 0
		for key in dict_to_print:
			if len(dict_to_print[key]) > 1:
				i += 1
				if i > 5:
					return
				else:
					print(key, dict_to_print[key])



	def get_and_pickle_all_mappings(self):
		# get seq_dicts first, so that can retrieve protein ID -> uniprot mappings for only those protein IDs w/ known protein sequences
		self.ensembl_prot_seq_dict = self.get_protein_seq('Ensembl', self.ensembl_protein_seq_file)
		self.refseq_prot_seq_dict = self.get_protein_seq('RefSeq', self.ref_seq_protein_seq_file)

		print('Number of Ensembl protein sequences:', len(self.ensembl_prot_seq_dict))
		print('Number of RefSeq protein sequences:', len(self.refseq_prot_seq_dict))

		pickle_dump(self.ensembl_prot_seq_dict, self.ensembl_prot_seq_dict_pickle_file)
		pickle_dump(self.refseq_prot_seq_dict, self.refseq_prot_seq_dict_pickle_file)

		# print(list(self.ensembl_prot_seq_dict.keys())[:10])
		# print(self.ensembl_prot_seq_dict['ENSP00000398030.1'])
		self.get_uniprot_mappings()

		# check
		print('Number of RefSeq/UniProt proteins in self.refseq_prot_uniprot_dict:', len(self.refseq_prot_uniprot_dict))
		print('Number of Ensembl/UniProt proteins in self.ensembl_prot_uniprot_dict:', len(self.ensembl_prot_uniprot_dict))
		print('Number of UniProt protein IDs with gene names:', len(self.uniprot_gene_name_dict))

		pickle_dump(self.ensembl_prot_uniprot_dict, self.ensembl_prot_uniprot_dict_pickle_file)
		pickle_dump(self.refseq_prot_uniprot_dict, self.refseq_prot_uniprot_dict_pickle_file)
		pickle_dump(self.uniprot_gene_name_dict, self.uniprot_gene_name_dict_pickle_file)

		# self.ensembl_prot_uniprot_dict = pickle_load(self.ensembl_prot_uniprot_dict_pickle_file)
		# self.refseq_prot_uniprot_dict = pickle_load(self.refseq_prot_uniprot_dict_pickle_file)
		# self.uniprot_gene_name_dict = pickle_load(self.uniprot_gene_name_dict_pickle_file)

		# self.test_print_key_with_multiple_values(self.ensembl_prot_uniprot_dict)
		# self.test_print_key_with_multiple_values(self.refseq_prot_uniprot_dict)
		# self.test_print_key_with_multiple_values(self.uniprot_gene_name_dict)


		# print(list(self.ensembl_prot_uniprot_dict.keys())[:10])
		# print(self.ensembl_prot_uniprot_dict['ENSP00000368109.1'])
		# print(list(self.refseq_prot_uniprot_dict.keys())[:10])
		# print(self.refseq_prot_uniprot_dict['XP_005267884.1'])
		# print(list(self.uniprot_gene_name_dict.keys())[:10])
		# print(self.uniprot_gene_name_dict['A0JNW5'])


		self.refseq_prot_gene_name_uniprot_dict = self.create_tuple_dict('RefSeq', self.refseq_prot_uniprot_dict, self.uniprot_gene_name_dict)
		self.ensembl_prot_gene_name_uniprot_dict = self.create_tuple_dict('Ensembl', self.ensembl_prot_uniprot_dict, self.uniprot_gene_name_dict)

		# print(self.refseq_prot_gene_name_uniprot_dict[('NP_001008222.1', 'AMY1A')])

		# check
		print('Number of RefSeq + gene_name tuples in self.refseq_prot_gene_name_uniprot_dict:', len(self.refseq_prot_gene_name_uniprot_dict))
		print('Number of Ensembl + gene_name tuples in self.ensembl_prot_gene_name_uniprot_dict:', len(self.ensembl_prot_gene_name_uniprot_dict))

		pickle_dump(self.ensembl_prot_gene_name_uniprot_dict, self.ensembl_prot_gene_name_uniprot_dict_pickle_file)
		pickle_dump(self.refseq_prot_gene_name_uniprot_dict, self.refseq_prot_gene_name_uniprot_dict_pickle_file)
		# return

		# self.refseq_prot_seq_dict = pickle_load(self.refseq_prot_seq_dict_pickle_file)
		# get self.refseq_mrna_refseq_prot_dict
		self.get_refseq_mrna_refseq_prot()

		print('Number of RefSeq mRNA to protein mappings:', len(self.refseq_mrna_refseq_prot_dict))
		pickle_dump(self.refseq_mrna_refseq_prot_dict, self.refseq_mrna_refseq_prot_dict_pickle_file)

		self.get_swissprot_id_seq()

		# get self.ensembl_specific_mappings
		self.get_ensembl_specific_mappings()



def main():
	script_dir = osp.dirname(__file__)
	data_dir = osp.join(script_dir, '..', 'data', 'processed')
	orig_data_dir = osp.join(script_dir, '..', 'data', 'original')

	m = Mappings(data_dir, orig_data_dir)
	m.get_and_pickle_all_mappings()

if __name__=='__main__':
	main()
