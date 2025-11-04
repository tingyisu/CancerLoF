'''
Run MODELLER to build 3D PPI structural models (heterodimers)
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import argparse
import os.path as osp
import os
from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
from simple_tools import pickle_load, pickle_dump, check_exist

# Override MyModel methods
class MyModel (automodel):
	"""Overrides MyModel methods. Renumbers chain residues based on UniProt numbering

	Args:
		automodel (object): modeller environ object.

	"""
	def __init__(self, env, alnfile, knowns, sequence, assess_methods, root_name, uniprot1_res_start, uniprot2_res_start):
		# Initialize parent class (automodel) with required arguments
		super().__init__(env, alnfile, knowns, sequence)
		
		# Store additional parameters
		self.assess_methods = assess_methods
		self.root_name = root_name
		self.uniprot1_res_start = uniprot1_res_start
		self.uniprot2_res_start = uniprot2_res_start

	def special_patches (self, aln):
		# Renumber residues in each chain to UniProt numbering
		# segment_ids should be ['A', 'B'] (modeling dimers here)
		self.rename_segments (segment_ids = [c.name for c in self.chains],
							  renumber_residues = [self.uniprot1_res_start, self.uniprot2_res_start])


class RunModeller:
	def __init__(self, interactome_data_dir, alignments_dir, modeller_dir):
		self.interactome_data_dir = interactome_data_dir
		self.alignments_dir = alignments_dir
		self.modeller_dir = modeller_dir
		self.converted_blast_alignments = pickle_load(osp.join(interactome_data_dir, 'converted_blast_alignments.pickle'))
		self.all_ppi_templates = pickle_load(osp.join(interactome_data_dir, 'all_combined_structural_templates.pickle')) # key = (uniprot1, uniprot2), value = (pdb_chain1, pdb_chain2)
		self.hetatm_records = pickle_load(osp.join(interactome_data_dir, 'hetatm_records.pickle')) # only stores for residues within BLAST local alignments
		self.uncompleted_models = []
		self.uniprot_residues_over_9999 = [] # [(uniprot, pdb_chain), etc] list of uniprot & pdb chain-pairs with aligned UniProt residues over 9999 (.pdb files cannot store past 9999)
		self.uncompleted_models_pickle_file = osp.join(interactome_data_dir, 'uncompleted_models.pickle')
		self.uniprot_residues_over_9999_pickle_file = osp.join(interactome_data_dir, 'uniprot_residues_over_9999.pickle')

	def get_modulo_10000(self, uniprot, pdb_chain, original_res_num):
		new_res_num = str(int(original_res_num) % 10000)
		if (uniprot, pdb_chain) not in self.uniprot_residues_over_9999:
			self.uniprot_residues_over_9999.append((uniprot, pdb_chain))
		return new_res_num

	def build_dimer_model(self, pdb, uniprot1, uniprot2, chain1, chain2, uniprot1_res_start, uniprot2_res_start):
		# run if haven't found model already
		fname = osp.join('_'.join([uniprot1, pdb, chain1, uniprot2, pdb, chain2]))
		# root_name = osp.join(self.modeller_dir, fname)
		# root_fname = osp.join('..', 'data', 'processed', 'modeller', fname)
		# model_name = osp.join(self.modeller_dir, fname + '.B99990001.pdb')
		model_fname = osp.join(fname + '.B99990001.pdb')
		ali_fname = osp.join('..', 'data', 'processed', 'alignments', fname + '.ali')
		print('Finding model for:', fname)
		if not check_exist(model_fname):
			try:
				log.minimal()
				env = Environ()
				env.io.hetatm = True # (NOTE: MODELLER seems to take HETATM MSE regardless of whether this is set to True/False)! It converts MSE to MET! https://salilab.org/archives/modeller_usage/2009/msg00185.html
				a = MyModel(env, alnfile=ali_fname, knowns=pdb, sequence=uniprot1+'_'+uniprot2, assess_methods=(assess.DOPE, assess.GA341), root_name=fname, uniprot1_res_start=uniprot1_res_start, uniprot2_res_start=uniprot2_res_start)                             
				a.starting_model = 1
				a.ending_model = 1
				a.max_molpdf = 1e8 # default is 1e7, increase for some not so great alignments/PDB files
				a.make()
				# ['_'.join([pdb, chain1]), '_'.join([pdb, chain2])]
				
			except:
				self.uncompleted_models.append(fname)
				print('ERROR!', fname)
		# else:
		# 	print('Model has already been found for:', uniprot, pdb_chain + '...')

	def build_all_dimer_models(self):
		total_num = len(self.all_ppi_templates)
		# curr = 1
		total_remodelled = []
		# ali_files = ''
		for (p1, p2) in self.all_ppi_templates:
			pdb_chain1, pdb_chain2 = self.all_ppi_templates[(p1, p2)]
			pdb, chain1 = pdb_chain1.split('_')
			_, chain2 = pdb_chain2.split('_')
			
			uniprot1_res_start, _, _, _, _, _ = self.converted_blast_alignments[(p1, pdb_chain1)]
			uniprot2_res_start, _, _, _, _, _ = self.converted_blast_alignments[(p2, pdb_chain2)]

			if (p1, pdb_chain1) in self.hetatm_records or (p2, pdb_chain2) in self.hetatm_records:
				# print('Remodelling dimer with HETATM entries:', p1, p2, pdb_chain1, pdb_chain2)
				total_remodelled.append((p1, p2))
				# ali_files += '\'' + '_'.join([p1, pdb_chain1, p2, pdb_chain2 + '.B99990001.pdb']) + '\' '
				new_uniprot1_res_start, new_uniprot2_res_start = uniprot1_res_start, uniprot2_res_start

				# take % 10000 if uniprot residue num > 9999; .pdb files cannot store residues over 9999 (or over 4 chars)
				if int(uniprot1_res_start) > 9999:
					new_uniprot1_res_start = self.get_modulo_10000(p1, pdb_chain1, uniprot1_res_start)
				if int(uniprot2_res_start) > 9999:
					new_uniprot2_res_start = self.get_modulo_10000(p2, pdb_chain2, uniprot2_res_start)

				self.build_dimer_model(pdb, p1, p2, chain1, chain2, int(new_uniprot1_res_start), int(new_uniprot2_res_start))
				# break
		print('Num remodelled:', len(total_remodelled))
		# print(total_remodelled)
		# print(ali_files)

		# pickle_dump(self.uncompleted_models, self.uncompleted_models_pickle_file)
		# pickle_dump(self.uniprot_residues_over_9999, self.uniprot_residues_over_9999_pickle_file)

		# print(self.uniprot_residues_over_9999)

def main():
	script_dir = osp.dirname(__file__)
	interactome_data_dir = osp.join(script_dir, '..', 'data', 'processed', 'interactome')
	alignments_dir = osp.join(script_dir, '..', 'data', 'processed', 'alignments')
	modeller_dir = osp.join(script_dir, '..', 'data', 'processed', 'modeller')

	r = RunModeller(interactome_data_dir, alignments_dir, modeller_dir)
	r.build_all_dimer_models()

if __name__=='__main__':
	main()
