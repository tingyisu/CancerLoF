'''
Create .ali files (containing BLASTP alignments of heterodimers) for running MODELLER
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import gemmi
import os.path as osp
import requests
from simple_tools import download_pdb_structure, pickle_load, pickle_dump, check_exist, check_create_dir

class HeterodimerFiles:
	def __init__(self, pdb_download_dir, interactome_data_dir, alignments_dir):
		self.pdb_download_dir = pdb_download_dir
		self.pdb_updated_dir = osp.join(pdb_download_dir, 'modified_cif')
		self.alignments_dir = alignments_dir
		check_create_dir(self.pdb_download_dir)
		check_create_dir(self.pdb_updated_dir)
		check_create_dir(self.alignments_dir)
		self.all_blast_best_alignments = pickle_load(osp.join(interactome_data_dir, 'all_blast_best_alignments_reduced_info.pickle')) # key = (protein, pdb_chain), value = [qstart, qend, sstart, send, protein alignment, template_chain alignment]
		self.all_ppi_templates = pickle_load(osp.join(interactome_data_dir, 'all_combined_structural_templates.pickle')) # key = (uniprot1, uniprot2), value = (pdb_chain1, pdb_chain2)
		self.label_auth_dict = {} # {('P23511', '4awl_A'):{label_seq_id: auth_seq_id}}
		self.corrected_blast_best_alignments = {} # corrected label_sstart and label_send for alignments with 'CRO' heteroatom(s) ('TYG' in PDB SeqRes)
		self.converted_blast_alignments = {}
		self.hetatm_records = {} # {('P31323', '4zp3_D'):[auth_seq_id1, auth_seq_id2, etc]} only includes templates with HETATM entries 
		self.CRO_hetatm = {} # {'6wvg_A':[auth_seq_id1, auth_seq_id2, et6c]} only includes templates with CRO HETATM entries; not needed in the end, but storing just in case

		# should probably pickle
		self.label_auth_dict_pickle_file = osp.join(interactome_data_dir, 'label_auth_dict.pickle') # only stores for residues within BLAST local alignments
		self.corrected_blast_best_alignments_pickle_file = osp.join(interactome_data_dir, 'corrected_blast_best_alignments_reduced_info.pickle') # key = (protein, pdb_chain), value = [qstart, qend, sstart, send, protein alignment, template_chain alignment]
		self.converted_blast_alignments_pickle_file = osp.join(interactome_data_dir, 'converted_blast_alignments.pickle')
		self.hetatm_records_pickle_file = osp.join(interactome_data_dir, 'hetatm_records.pickle') # only stores for residues within BLAST local alignments
		self.CRO_hetatm_pickle_file = osp.join(interactome_data_dir, 'CRO_hetatm.pickle') # stores all CRO residues in the corresponding PDB chain
		
		# self.label_auth_dict = pickle_load(self.label_auth_dict_pickle_file) # {(P23511, 4awl_A):{label_seq_id: auth_seq_id}}
		# self.converted_blast_alignments = pickle_load(self.converted_blast_alignments_pickle_file)
	
	# creates new PDB chain using Gemmi
	# also updates self.label_auth_dict
	def create_new_chain_with_specified_label_residue_range(self, uniprot, pdb_chain, structure_chain, label_sstart, label_send):
		
		# create new chain with same identifier
		new_chain = gemmi.Chain(structure_chain.name)

		in_label_auth_dict = False
		if (uniprot, pdb_chain) not in self.label_auth_dict:
			self.label_auth_dict[(uniprot, pdb_chain)] = {}
		else:
			in_label_auth_dict = True

		# get CRO residues, need to know how many there are to shift label_sstart and label_send
		corrected_label_sstart, corrected_label_send = label_sstart, label_send
		CRO_residues = [res for res in structure_chain if res.name == 'CRO' and res.het_flag == 'H']
		for CRO_residue in CRO_residues:
			if CRO_residue.label_seq and CRO_residue.seqid.num: # Ensure the residue has both label_seq_id and auth_seq_id
				label_seq_id = CRO_residue.label_seq  # Gemmi's label_seq_id (integer)
				if label_seq_id + 2 <= label_sstart: # occurs before label_sstart
					corrected_label_sstart -= 2
					corrected_label_send -= 2
				elif label_seq_id + 2 <= label_ssend: # occurs in-between label_sstart and label_send
					corrected_label_send -= 2
					print('Occurs in-between label_sstart and label_send:', uniprot, pdb_chain)

		for residue in structure_chain:
			# Ensure the residue has both label_seq_id and auth_seq_id
			if residue.label_seq and residue.seqid.num:
				label_seq_id = residue.label_seq  # Gemmi's label_seq_id (integer)
				auth_seq_id = residue.seqid.num  # Gemmi's auth_seq_id (integer)

				if label_seq_id >= corrected_label_sstart and label_seq_id <= corrected_label_send: # keep residues in specified range

					# add residue to self.label_auth_dict and new_chain
					if not in_label_auth_dict: # add if (uniprot, pdb_chain) doesn't exist yet
						self.label_auth_dict[(uniprot, pdb_chain)][label_seq_id] = auth_seq_id
						if residue.het_flag == 'H': # is '' for ATOM entries
							# print('HETATM:', uniprot, pdb_chain, auth_seq_id)
							self.hetatm_records.setdefault((uniprot, pdb_chain), []).append((auth_seq_id, residue.name)) # get residue name, so that can check later

					new_chain.add_residue(residue)

				# store CRO entries (even ones outside of specified range b/c they shift the residue numbers)
				if not in_label_auth_dict:
					if residue.name == 'CRO': # CRO hetatm, is 'TYG' in SeqRes
						self.CRO_hetatm.setdefault(pdb_chain, []).append(auth_seq_id)

		return new_chain, corrected_label_sstart, corrected_label_send				
			

	# create new modified .cif file with only the 2 specified chains in chains_to_keep
	# for each chain, only keeps residues within BLASTP local alignment 
	def create_cif_with_two_heterodimer_chains(self, p1, p2, pdb, new_cif_fname, chains_to_keep):
		'''
		p1 = UniProt protein 1
		p2 = UniProt protein 2
		pdb = name of PDB structure
		new_cif_fname = name of new .cif to create with only the two structural template chains (residues within each chain are constrained to a specified range)
		chains_to_keep = list of the two structural template chains
		'''

		# Load the structure from the CIF file
		if not check_exist(osp.join(self.pdb_updated_dir, new_cif_fname + '.cif')):

			# print('Creating .cif for PDB chain-pairs:', pdb, chains_to_keep)
			input_file = osp.join(self.pdb_download_dir, pdb + '.cif')
			output_file = osp.join(self.pdb_updated_dir, new_cif_fname + '.cif') 

			# Read the structure
			structure = gemmi.read_structure(input_file)

			# Create a new structure to hold the selected chains
			new_structure = gemmi.Structure()
			new_structure.name = new_cif_fname

			# Create a new model in the new structure 
			new_model = gemmi.Model('1')

			# necessary variables
			pdb_chain1 = '_'.join([pdb, chains_to_keep[0]])
			pdb_chain2 = '_'.join([pdb, chains_to_keep[1]])
			qstart1, qend1, label_sstart1, label_send1, qalignment1, salignment1 = self.all_blast_best_alignments[(p1, pdb_chain1)]
			qstart2, qend2, label_sstart2, label_send2, qalignment2, salignment2 = self.all_blast_best_alignments[(p2, pdb_chain2)]

			# get chains from the first (and only) model (structure[0])
			structure_chain1 = structure[0][chains_to_keep[0]]
			structure_chain2 = structure[0][chains_to_keep[1]]

			# create new chain with specified label_sstart1 and label_sstart2
			new_chain1, corrected_label_sstart1, corrected_label_send1 = self.create_new_chain_with_specified_label_residue_range(p1, pdb_chain1, structure_chain1, int(label_sstart1), int(label_send1))
			new_chain2, corrected_label_sstart2, corrected_label_send2 = self.create_new_chain_with_specified_label_residue_range(p2, pdb_chain2, structure_chain2, int(label_sstart2), int(label_send2))

			# update self.corrected_blast_best_alignments
			if (p1, pdb_chain1) not in self.corrected_blast_best_alignments:
				self.corrected_blast_best_alignments[(p1, pdb_chain1)] = [qstart1, qend1, corrected_label_sstart1, corrected_label_send1, qalignment1, salignment1]
			if (p2, pdb_chain2) not in self.corrected_blast_best_alignments:
				self.corrected_blast_best_alignments[(p2, pdb_chain2)] = [qstart2, qend2, corrected_label_sstart2, corrected_label_send2, qalignment2, salignment2]

			# update new structure with new_chain1 and new_chain2
			new_model.add_chain(new_chain1)
			new_model.add_chain(new_chain2)

			# Add the new model to the new structure
			new_structure.add_model(new_model)

			# Write the new structure to a CIF file
			new_structure.make_mmcif_document().write_file(output_file)

		# else:
		# 	print('.cif file was already created for PDB chain-pair:', pdb, chains_to_keep)


	# updates self.converted_blast_alignments
	# accounts for structural gaps (gap in PDB); amino acid residue is in PDB SeqRes, but is not structurally resolved (not in .cif file)
	def convert_blast_alignments(self, uniprot, pdb_chain):
		
		# print('Converting BLAST alignment for:', uniprot, pdb_chain)
		
		# convert sstart and send from label_seq_id format to auth_seq_id format
		qstart, qend, label_sstart, label_send, qalignment, salignment = self.corrected_blast_best_alignments[(uniprot, pdb_chain)] 
		# get only label_seq_ids that exist in the mmCIF file
		label_residues_that_exist = [i for i in range(int(label_sstart), int(label_send) + 1) if i in self.label_auth_dict[(uniprot, pdb_chain)]] # label_seq_ids that exist in the mmCIF file

		# change residues in salignment that do not exist in the mmCIF file to '-'
		# change residues that are HETATM entries (in the mmCIF file, but not in the PDB SEQRES) to '.'
		all_hetatms_auth_seq_ids = []
		if (uniprot, pdb_chain) in self.hetatm_records:
			all_hetatms_auth_seq_ids = [t[0] for t in self.hetatm_records[(uniprot, pdb_chain)]] # t[0] is the auth_seq_id (int) and t[1] is the name of the HETATM residue (.e.g 'TPO')
		converted_salignment = ''
		i = int(label_sstart) # sanity check
		for residue in salignment:
			if residue != '-':
				if i in label_residues_that_exist:
					auth_seq_id = self.label_auth_dict[(uniprot, pdb_chain)][i]
					if auth_seq_id in all_hetatms_auth_seq_ids: # check whether this particular residue is a HETATM
						converted_salignment += '.' # MODELLER's HETATM representation
					else: 
						converted_salignment += residue
				else:
					converted_salignment += '-'
				i += 1
			else:
				converted_salignment += '-'

		if i != (int(label_send) + 1): # sanity check
			print('ERROR! Something went wrong with the label iterations!!', i, label_send)

		auth_sstart = self.label_auth_dict[(uniprot, pdb_chain)][label_residues_that_exist[0]]
		auth_send = self.label_auth_dict[(uniprot, pdb_chain)][label_residues_that_exist[len(label_residues_that_exist)-1]]
		self.converted_blast_alignments[(uniprot, pdb_chain)] = [qstart, qend, str(auth_sstart), str(auth_send), qalignment, converted_salignment]

	# Biopython version for downloading mmCIF files
	def get_pdb_structures_to_download(self):
		pdb_structures_to_download = set()
		for protein_pair in self.all_ppi_templates: # key = (uniprot1, uniprot2), value = (pdb_chain1, pdb_chain2)
			pdb_chain1 = self.all_ppi_templates[protein_pair][0]
			pdb, _ = pdb_chain1.split('_')
			pdb_structures_to_download.add(pdb)

		pdb_structures_to_manually_download = set()
		# download CIF file
		for pdb in pdb_structures_to_download:
			download_pdb_structure(pdb, self.pdb_download_dir)
			if not check_exist(osp.join(self.pdb_download_dir, pdb + '.cif')):
				pdb_structures_to_manually_download.add(pdb)
		for pdb in pdb_structures_to_manually_download:
			print(pdb)

	# requests version for downloading mmCIF files
	def download_pdb_structures_requests(self):
		pdb_structures_to_download = set()
		for protein_pair in self.all_ppi_templates: # key = (uniprot1, uniprot2), value = (pdb_chain1, pdb_chain2)
			pdb_chain1 = self.all_ppi_templates[protein_pair][0]
			pdb, _ = pdb_chain1.split('_')
			pdb_structures_to_download.add(pdb)
		print('Number of PDB structures to download:', len(pdb_structures_to_download))

		pdb_structures_to_manually_download = set()
		# download CIF file
		for pdb in pdb_structures_to_download:
			pdb_id = pdb.upper()
			url = 'https://files.rcsb.org/download/' + pdb_id + '.cif'
			response = requests.get(url)
			if response.status_code == 200:
				with open(osp.join(self.pdb_download_dir, pdb + '.cif'), 'wb') as f:
					f.write(response.content)
				print('mmCIF file', pdb, 'downloaded successfully.')
			else:
				print('Failed to download mmCIF file:', pdb)
				pdb_structures_to_manually_download.add(pdb)
		for pdb in pdb_structures_to_manually_download:
			print(pdb)

	def create_heterodimer_ali_file_for_modeller(self, uniprot1, uniprot2, pdb_chain1, pdb_chain2, fname):
		# create .ali file for a heterodimer
		# print('Creating .ali file for:', uniprot1, uniprot2, pdb_chain1, pdb_chain2)
		seq_name = '_'.join([uniprot1, uniprot2])
		qstart1, qend1, auth_sstart1, _, converted_qalignment1, converted_salignment1 = self.converted_blast_alignments[(uniprot1, pdb_chain1)]
		qstart2, qend2, _, auth_send2, converted_qalignment2, converted_salignment2 = self.converted_blast_alignments[(uniprot2, pdb_chain2)]
		pdb, chain1 = pdb_chain1.split('_')
		_, chain2 = pdb_chain2.split('_')
		pdb_chain1_chain2 = '_'.join([pdb, chain1, chain2])

		fname_ali = osp.join(self.alignments_dir, fname + '.ali')
		# cif_fname = osp.join(self.pdb_updated_dir, fname + '.cif') # unable to do this...MODELLER doesn't seem to recognize WSL file system
		cif_fname = osp.join('..', 'data', 'processed', 'pdb_cif', 'modified_cif', fname + '.cif')
		# if not check_exist(fname):
		with open(fname_ali, 'w') as f:
			f.write('>P1;' + pdb + '\n')
			f.write('structure:' + cif_fname + ':' + auth_sstart1 + ':' + chain1 + ':' + auth_send2 + ':' + chain2 + '::::' + '\n')
			f.write(converted_salignment1 + '/' + converted_salignment2 + '*' + '\n\n')
			f.write('>P1;' + seq_name + '\n')
			f.write('sequence:' + seq_name + ':' + qstart1 + '::' + qend1 + '::' + qstart2 + '::' + qend2 + ':' + '\n')
			f.write(converted_qalignment1 + '/' + converted_qalignment2 + '*' + '\n')


	def create_heterodimer_files(self):
		total_num = len(self.all_ppi_templates)
		curr = 1
		for (p1, p2) in self.all_ppi_templates: # key = (uniprot1, uniprot2), value = (pdb_chain1, pdb_chain2)
			
			print('Creating files for:', str(curr) + '/' + str(total_num))
			pdb_chain1, pdb_chain2 = self.all_ppi_templates[(p1, p2)]
			pdb, chain1 = pdb_chain1.split('_')
			_, chain2 = pdb_chain2.split('_')
			fname = '_'.join([p1, pdb_chain1, p2, pdb_chain2])
			chains_to_keep = [chain1, chain2]

			# create new .cif file with only the 2 structural template chains
			# update self.label_auth_dict with the mappings for the template chains
			self.create_cif_with_two_heterodimer_chains(p1, p2, pdb, fname, chains_to_keep)

			# print(self.converted_blast_alignments)
			# print(self.label_auth_dict)

			# get converted BLASTP alignment, need for MODELLER
			if (p1, pdb_chain1) not in self.converted_blast_alignments:
				self.convert_blast_alignments(p1, pdb_chain1)
			if (p2, pdb_chain2) not in self.converted_blast_alignments:
				self.convert_blast_alignments(p2, pdb_chain2)

			# create heterodimer ali file for input into MODELLER
			self.create_heterodimer_ali_file_for_modeller(p1, p2, pdb_chain1, pdb_chain2, fname)

			curr += 1

		# pickle dump
		pickle_dump(self.label_auth_dict, self.label_auth_dict_pickle_file)
		pickle_dump(self.converted_blast_alignments, self.converted_blast_alignments_pickle_file)
		pickle_dump(self.hetatm_records, self.hetatm_records_pickle_file)
		pickle_dump(self.CRO_hetatm, self.CRO_hetatm_pickle_file)

def main():
	script_dir = osp.dirname(__file__)
	interactome_data_dir = osp.join(script_dir, '..', 'data', 'processed', 'interactome')
	pdb_download_dir = osp.join(script_dir, '..', 'data', 'processed', 'pdb_cif')
	alignments_dir = osp.join(script_dir, '..', 'data', 'processed', 'alignments')

	# download_pdb_structure('4awl', '.')
	h = HeterodimerFiles(pdb_download_dir, interactome_data_dir, alignments_dir)
	# h.download_pdb_structures_requests() # FIRST MAKE SURE TO MANUALLY DOWNLOAD (FROM INTERNET) .CIF FILES THAT ARE OBSOLETE
	h.create_heterodimer_files()

if __name__ == '__main__':
	main()