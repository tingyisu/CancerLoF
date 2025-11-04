# '''
# Splits and processes COSMIC genome screen mutations (in batches)
# ----------------------------------------------
# Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
# '''

#!/bin/bash
# split genome screen files into different files (b/c each has 50 million + mutations and will get seg fault)
split -l 800000 --numeric-suffixes Cosmic_GenomeScreensMutant_v100_GRCh38.tsv cosmic_genome_screen_mutants_split_files/cosmic_genome_screen_mutants_split_

# move the header from the first split file into its own header file
cd cosmic_genome_screen_mutants_split_files
head -n1 cosmic_genome_screen_mutants_split_00 >> header.tsv
sed -i -e '1,1 d' cosmic_genome_screen_mutants_split_00

# process all COSMIC genome screen split files
for i in {0..n} # n = # of splits-1
do
	python process_cosmic_genome_screen_mutants.py -s $i
	echo finished processing genome screen file number $i ...
done

# COSMIC GENOME SCREEN MUTANTS: concatenate all the processed split missense, synonymous, and nonsense mutation files together
echo concatenating the genome screen mutant files...
cat ../data/original/cosmic_genome_screen_mutants_split_files/processed_cosmic_genome_screen_mutants_missense_split_00.tsv > ../data/processed/mutations/cosmic_genome_screen_missense.tsv
cat ../data/original/cosmic_genome_screen_mutants_split_files/processed_cosmic_genome_screen_mutants_synonymous_split_00.tsv > ../data/processed/mutations/cosmic_genome_screen_synonymous.tsv
cat ../data/original/cosmic_genome_screen_mutants_split_files/processed_cosmic_genome_screen_mutants_nonsense_split_00.tsv > ../data/processed/mutations/cosmic_genome_screen_nonsense.tsv

num_to_divide=10
for i in {1..n} # n = # of splits-1
do
	result=$(expr $i / $num_to_divide)
	if (($result == 0))
	then
		# missense mutations (skip header line)
		awk FNR!=1 ../data/original/cosmic_genome_screen_mutants_split_files/processed_cosmic_genome_screen_mutants_missense_split_"$result$i".tsv >> ../data/processed/mutations/cosmic_genome_screen_missense.tsv
		# synonymous mutations (skip header line)
		awk FNR!=1 ../data/original/cosmic_genome_screen_mutants_split_files/processed_cosmic_genome_screen_mutants_synonymous_split_"$result$i".tsv >> ../data/processed/mutations/cosmic_genome_screen_synonymous.tsv
		# nonsense mutations (skip header line)
		awk FNR!=1 ../data/original/cosmic_genome_screen_mutants_split_files/processed_cosmic_genome_screen_mutants_nonsense_split_"$result$i".tsv >> ../data/processed/mutations/cosmic_genome_screen_nonsense.tsv
		echo finished concatenating file number "$result$i"
	else
		# missense mutations (skip header line)
		awk FNR!=1 ../data/original/cosmic_genome_screen_mutants_split_files/processed_cosmic_genome_screen_mutants_missense_split_"$i".tsv >> ../data/processed/mutations/cosmic_genome_screen_missense.tsv
		# synonymous mutations (skip header line)
		awk FNR!=1 ../data/original/cosmic_genome_screen_mutants_split_files/processed_cosmic_genome_screen_mutants_synonymous_split_"$i".tsv >> ../data/processed/mutations/cosmic_genome_screen_synonymous.tsv
		# nonsense mutations (skip header line)
		awk FNR!=1 ../data/original/cosmic_genome_screen_mutants_split_files/processed_cosmic_genome_screen_mutants_nonsense_split_"$i".tsv >> ../data/processed/mutations/cosmic_genome_screen_nonsense.tsv
		echo finished concatenating file number $i
	fi
done

# do the following for a sanity check...
# check number of mutations in the concatenated file (should be -n, = # of splits-1)
wc -l ../data/processed/mutations/cosmic_genome_screen_missense.tsv
wc -l ../data/processed/mutations/cosmic_genome_screen_synonymous.tsv
wc -l ../data/processed/mutations/cosmic_genome_screen_nonsense.tsv
# the above should be (should be -n, = # of splits-1) of the corresponding split files below (to account for repeated header lines)
wc -l ../data/original/cosmic_genome_screen_mutants_split_files/processed_cosmic_genome_screen_mutants_missense_split_*.tsv 
wc -l ../data/original/cosmic_genome_screen_mutants_split_files/processed_cosmic_genome_screen_mutants_synonymous_split_*.tsv 
wc -l ../data/original/cosmic_genome_screen_mutants_split_files/processed_cosmic_genome_screen_mutants_nonsense_split_*.tsv 
