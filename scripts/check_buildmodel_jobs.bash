#!/bin/bash
# '''
# Check whether BuildModel jobs were completed successfully
# ----------------------------------------------
# Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
# '''

cd /home/username/scratch/foldx_buildmodel_all # need to change username
for i in {0..n}; do
	dir="split_$i"
	if [[ -d "$dir" ]]; then
		cd "$dir" || continue
		for file in *.log; do
			[[ -e "$file" ]] || continue
			if [[ "$file" =~ ([0-9]{8,})\.log$ ]]; then
				num=${BASH_REMATCH[1]}
				if (( num >= 58875887 )); then # need to change 58875887 to the number of mutations you have
					count=$(grep -c 'Cleaning BuildModel...DONE' "$file") # same as cat | grep | wc -l
					echo "$dir $count" $file
				fi
			fi
		done
		cd ..
	fi
done