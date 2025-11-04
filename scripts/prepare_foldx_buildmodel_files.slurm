#!/bin/bash
# '''
# Run each split FoldX BuildModel job for preparing configuration files
# ----------------------------------------------
# Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
# '''

#SBATCH -n 1
#SBATCH --mem-per-cpu 1G
#SBATCH -t 1-00:00
#SBATCH -o prepare_buildmodel.%N.%j.log
#SBATCH -e prepare_buildmodel.%N.%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email_address
#SBATCH --account=def-*
python3 -u prepare_foldx_buildmodel_files.py -s /home/username/projects/def-*/username/foldx -c /home/username/scratch -a def-* -f foldx_20251231
