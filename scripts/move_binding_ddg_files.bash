#!/bin/bash
# '''
# Move FoldX PSSM (binding DDG) output files
# ----------------------------------------------
# Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
# '''

binding_ddg_dir=/home/username/scratch/binding_ddg
mkdir $binding_ddg_dir
echo Created individual binding DDG folder in the scratch dir
# move PSSM files
cd /home/username/scratch/foldx_pssm_all
for i in {0..28}
do
	# echo ${dir}
	rm -r split_$i/molecules/
	echo Removed FoldX folder of molecules in split_$i
	for dir in split_$i/*/
	do
		f_name="$(basename ${dir%/})"
		dimer="${f_name%-*}"
		binding_ddg_fname="${dir}PSSM_${dimer}_Repair.txt"
		# echo "${dimer}" "${binding_ddg_fname}"
		cp $binding_ddg_fname $binding_ddg_dir/"PSSM_${f_name}.txt"
		# mv $binding_ddg_fname $binding_ddg_dir
	done
   	echo Moved binding DDG files in split_$i to own directories...
done