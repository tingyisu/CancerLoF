'''
Optional: Comparing interfacial residues in structural templates vs. those in modeller heterodimers (after undergoing FoldX RepairPDB)
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os.path as osp
from simple_tools import pickle_load, pickle_dump

def main():
	script_dir = osp.dirname(__file__)
	interactome_data_dir = osp.join(script_dir, '..', 'data', 'processed', 'interactome')
	structural_templates_interfacial_residues = pickle_load(osp.join(interactome_data_dir, 'all_combined_interfacial_residues_modeled.pickle'))
	modeller_heterodimer_interfacial_residues = pickle_load(osp.join(interactome_data_dir, 'modeller_heterodimer_interfacial_residues.pickle'))
	print(len(structural_templates_interfacial_residues), len(modeller_heterodimer_interfacial_residues))

	num_lost_interfacial_residues = 0
	num_larger = 0
	num_smaller = 0
	num_same = 0
	num_same_diff_irs = 0
	num_same_same_irs = 0

	ppis_lost_interfacial_residues = []
	ppis_lost_interfacial_residues_pickle_file = osp.join(interactome_data_dir, 'ppis_without_interfacial_residues_after_modeling.pickle')

	for pair in structural_templates_interfacial_residues:
		structural_template = structural_templates_interfacial_residues[pair]
		modeller_heterodimer = modeller_heterodimer_interfacial_residues[pair]

		# NOTE: ''.split(',') returns ['']
		st_residues1, st_residues2 = structural_template[0].split(','), structural_template[1].split(',')
		mh_residues1, mh_residues2 = modeller_heterodimer[0].split(','), modeller_heterodimer[1].split(',')

		if modeller_heterodimer[0] == '' or modeller_heterodimer[1] == '':
			print('NO IRs between modeller heterodimers:', pair, mh_residues1, mh_residues2, st_residues1, st_residues2)
			ppis_lost_interfacial_residues.append((pair))
			num_lost_interfacial_residues += 1

		elif len(mh_residues1) + len(mh_residues2) > len(st_residues1) + len(st_residues2):
			num_larger += 1

		elif len(mh_residues1) + len(mh_residues2) < len(st_residues1) + len(st_residues2):
			num_smaller += 1

		else: # same number of IRs
			if mh_residues1 != st_residues1 or mh_residues2 != st_residues2:
				num_same_diff_irs += 1
			else:
				num_same_same_irs += 1
			num_same += 1

	print('Number of PPIs with the same # of IRs after modeling:', num_same, num_same_diff_irs, num_same_same_irs)
	print('Number of PPIs that lost IRs after modeling:', num_lost_interfacial_residues)
	print('Number of PPIs that have more IRs after modeling:', num_larger)
	print('Number of PPIs that have less IRs after modeling:', num_smaller)

	print('Total number of PPIs:', num_same + num_lost_interfacial_residues + num_larger + num_smaller)

	pickle_dump(ppis_lost_interfacial_residues, ppis_lost_interfacial_residues_pickle_file)


if __name__ == '__main__':
	main()