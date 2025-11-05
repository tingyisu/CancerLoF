#!/bin/bash
# '''
# Check whether BuildModel jobs were completely successfully
# ----------------------------------------------
# Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
# '''

cd /home/username/scratch/foldx_buildmodel_all
folding_ddg_dir=/home/username/scratch/folding_ddg
mkdir $folding_ddg_dir
echo Created individual folding DDG folder in the scratch dir
for i in {0..85}
do
	rm -r split_$i/molecules/
	echo Removed FoldX folder of molecules in split_$i
	for dir in split_$i/*/
	do
   		f_name="$(basename ${dir%/})"
   		dimer="${f_name%-*}"
   		folding_ddg_fname="${dir}Dif_${dimer}_Repair.fxout"
   		# echo "${rsa_fname}" "${folding_ddg_fname}"
   		# echo $folding_ddg_fname $folding_ddg_dir/"Dif_${f_name}.fxout"
   		cp $folding_ddg_fname $folding_ddg_dir/"Dif_${f_name}.fxout"
	done
   	echo Moved folding DDG files in split_$i to own directories...
done