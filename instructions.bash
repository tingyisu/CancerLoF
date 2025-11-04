#!/bin/bash

# '''
# Instructions for running scripts in CancerLoF
# This is the most up-to-date and complete pipeline with all mutation types (dbSNP, ClinVar, and COSMIC)
# ----------------------------------------------
# Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
# '''

# Please follow the steps below carefully
# Some steps are very memory/time consuming to run and should ideally be run in parallel on a server (scripts for submitting SLURM jobs are provided for these steps).

# first make sure you're in the MutDeleteriousness directory
# and create the following folders to store data
mkdir data
cd data
mkdir original
mkdir processed

# make subfolders in the data/processed dir
cd processed
mkdir interactome
mkdir mutations
mkdir mutations_final
mkdir edgotypes
mkdir edgotypes_final
mkdir nonsense_on_si


# ----------DOWNLOAD THE FOLLOWING DATA FILES:----------
cd ../original
# **1. UniProt mappings file**
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
gunzip HUMAN_9606_idmapping.dat.gz

# **2. RefSeqGene mappings (need for mapping mRNA accession NM_* to protein accession NP_* for mutations)** 
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/LRG_RefSeqGene

# **3. RefSeq protein transcript sequence files**
mkdir ref_seq
cd ref_seq
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.*.protein.faa.gz
gunzip human.*.protein.faa.gz
cat * > ../../processed/merged_ref_seq_protein.fa

# download Ensembl mRNA/protein transcripts to UniProt mappings
wget ftp://ftp.ensembl.org/pub/current_tsv/homo_sapiens/Homo_sapiens.GRCh38.112.uniprot.tsv.gz
gunzip Homo_sapiens.GRCh38.112.uniprot.tsv.gz

# download protein sequences form Ensembl protein acessions (ENSP*)
wget ftp://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
gunzip Homo_sapiens.GRCh38.pep.all.fa.gz

# **4. Uniprot reviewed (UniProtKB/SWISS-PROT) human reference proteome in list and FASTA format**
cd ../
wget 'https://rest.uniprot.org/uniprotkb/stream?format=list&query=%28%28taxonomy_id%3A9606%29%20AND%20%28reviewed%3Atrue%29%29%20AND%20%28model_organism%3A9606%29%20AND%20%28reviewed%3Atrue%29' -O uniprot_reviewed_human_proteome.list
wget 'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28taxonomy_id%3A9606%29%20AND%20%28reviewed%3Atrue%29%29%20AND%20%28model_organism%3A9606%29%20AND%20%28reviewed%3Atrue%29' -O uniprot_reviewed_human_proteome.fasta

# **5. HI-Union human binary protein-protein interaction dataset**
wget http://www.interactome-atlas.org/data/HI-union.tsv

# **6. IntAct human binary protein-protein interaction dataset**
wget ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip
unzip intact.zip
rm intact.zip intact_negative.txt

# **7. Protein Data Bank (PDB) SeqRes chain sequences**
wget https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt

# **8. File containing the resolutions of PDB structures**
wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/resolu.idx

# **9. ClinVar (Mendelian disease-causing) mutations
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
gunzip variant_summary.txt.gz

# **10. COSMIC (cancer-associated) mutations
# download COSMIC mutations
# (1) Mutations in Cancer Gene Census (Cosmic_MutantCensus_v100_GRCh38.tsv) --> mutations in cancer-driving genes
# (2) Genome screen mutants (coding mutations across entire cancer genome) (Cosmic_GenomeScreensMutant_v100_GRCh38.tsv) --> cancer-associated
# untar (tar -xvf) & unzip (gunzip)

# ----------SELECT PDB STRUCTURAL TEMPLATES FOR HI-UNION AND INTACT REFERENCE PPIS----------
cd ../../scripts
# **1. Process all mappings necessary for processing & mapping mutations to the structural interactomes**
python3 process_mappings.py # need to get from human_interactome dir

# **2. Process HI-union and IntAct data to get human reference PPIs**
python3 process_hiunion_data.py # also processes PDB SeqRes sequences for input into BLASTP
python3 process_intact_data.py

# **3. Perform BLASTP sequence alignment of SwissProt proteins (within HI-union and IntAct) against PDB chains**
cd ../data/processed/interactome
makeblastdb -in pdb_seqres_blast.fasta -dbtype prot -out pdb_seqres_db
# -----HI-union-----
blastp -db pdb_seqres_db -query hiunion_uniprot_sequences_blast.fasta -out hiunion_blast_results_with_seq.tsv -outfmt "6 std qseq sseq"
# keep all alignments with E-value < 10^-5 and pick the best alignment (listed first) for each pair of SwissProt protein and PDB chain
awk -F"\t" '$11<10^-5' "hiunion_blast_results_with_seq.tsv" > "hiunion_blast_results_with_seq_eval-5.tsv"
cut -d$'\t' -f 1-2 hiunion_blast_results_with_seq_eval-5.tsv | awk -F"\t" '!seen[$1, $2]++' > hiunion_blast_for_finding_interactions.tsv
# -----IntAct-----
blastp -db pdb_seqres_db -query intact_uniprot_sequences_blast.fasta -out intact_blast_results_with_seq.tsv -outfmt "6 std qseq sseq"
# keep all alignments with E-value < 10^-5 and pick the best alignment (listed first) for each pair of SwissProt protein and PDB chain
awk -F"\t" '$11<10^-5' "intact_blast_results_with_seq.tsv" > "intact_blast_results_with_seq_eval-5.tsv"
cut -d$'\t' -f 1-2 intact_blast_results_with_seq_eval-5.tsv | awk -F"\t" '!seen[$1, $2]++' > intact_blast_for_finding_interactions.tsv

# **4. Use a distance-based approach to find interfacial residues between interacting chain-pairs in PDB structures**
# NOTE: Please use a server for this step as it will take several months if you don't run multiple jobs in parallel on a server 
# The server we used was Compute Canada/Digital Research Alliance of Canada, which uses the SLURM job scheduler

# (1) CREATE A SCRIPT DIRECTORY called 'interfacial_residues' on your server and COPY the following scripts over (this is where you will run your scripts):
# If you're using Compute Canada, you should create the 'interfacial_residues' folder within your scratch directory (e.g. /home/username/scratch/interfacial_residues)
# The scratch directory can store up to 20TB of data, but they will be deleted after 60 days
# (1) process_blast_results.py
# (2) get_hiunion_pdb_chain_dict.py
# (3) get_intact_pdb_chain_dict.py
# (4) split_pdb_chain_dict.py
# (5) create_residues_slurm_files.py
# (6) create_bash_for_submitting_residues_slurm_jobs.py
# (7) hiunion_residues_split.py
# (8) intact_residues_split.py
# (9) memo_residues.py
# (10) simple_tools.py
# (11) hiunion_write_interactions.py
# (12) intact_write_interactions.py
# (13) get_auth_label_dict.py
# (14) get_auth_label_dict.slurm

# (2) CREATE AN OUTPUT DIRECTORY (-o --output_directory) called 'files' under your 'interfacial_residues' folder on your server and move the following files from data/processed/interactome there:
# If you're using Compute Canada, you should create the 'interfacial_residues/files' folder within your scratch directory (e.g. /home/username/scratch/interfacial_residues/files)
# (1) hiunion_blast_for_finding_interactions.tsv
# (2) HI-union_unique_uniprot_pairs.pickle
# (3) intact_blast_for_finding_interactions.tsv
# (4) intact_uniprot_pairs_physical.pickle
# (5) HI-union_interacting_pairs.pickle
# (6) intact_gene_dict_physical.pickle
# (7) intact_id_dict_physical.pickle

# (3) CREATE A PDB DIRECTORY (-d --pdb_download_path) called 'pdb_cif' to store downloaded PDB mmCIF files
# If you're using Compute Canada, you should create the 'pdb_cif' folder within your scratch directory (e.g. /home/username/scratch/pdb_cif)

# the following arguments are REQUIRED (please use ABSOLUTE PATHS):
# -o --output_directory (this is the output dir you created in step (2))
# -d --pdb_download_path (this is the pdb dir you created in step (3))
# -u --num_splits_for_hiunion (how many splits/jobs to create for HI-Union)
# -i --num_splits_for_intact (how many splits/jobs to create for IntAct)
# -n --num_days (max number of days to run each job)
# -e --email_address (input your email address if you would like to receive emails on when a job has completed or failed)
# -a --account (account name/compute canada group, e.g. def-*)

# run the following commands one-by-one in the script directory on your server:

# (4) Find possible homologs from BLASTP alignment results
python3 process_blast_results.py -o --output_directory 

# (5) Get HI-union and IntAct pdb_chain_dicts (dictionaries of pairs of PDB chains) & download PDB structures in .cif (mmCIF) format
python3 get_hiunion_pdb_chain_dict.py -d --pdb_download_path -o --output_directory
python3 get_intact_pdb_chain_dict.py -d --pdb_download_path -o --output_directory

# (6) Split the HI-union pdb_chain_dict into 20 jobs and the IntAct pdb_chain_dict into 50 jobs
python3 split_pdb_chain_dict.py -o --output_directory -u 20 -i 50 # 20 splits for HI-Union, 50 for IntAct

# (7) Create SLURM job files to find interfacial residues between pairs of PDB chains
python3 create_residues_slurm_files.py -d --pdb_download_path -o --output_directory -u 20 -i 50 -n 5 -e --email_address -a --account 
python3 create_bash_for_submitting_residues_slurm_jobs.py -u 20 -i 50 

# (8) Submit SLURM jobs
bash submit_residues_slurm_jobs.bash

# (9) If any of the SLURM jobs were not completed in the specified amount of time, 'JOB FAILED' messages will be sent to your email address
# Simply submit the uncompleted jobs again and they will continue to run from where they stopped at; resubmit as many times as needed
sbatch _name_of_uncompleted_job.slurm # replace _name_of_uncompleted_job with the name of the job (e.g. hiunion_1)

# (10) Once all SLURM jobs are completed, combine the memo files in your -o --output_directory
cd output_directory # replace with your output dir name
cat hiunion_memo_residues_split*.tsv > hiunion_memo_residues_combined.tsv
cat intact_memo_residues_split*.tsv > intact_memo_residues_combined.tsv # excluding overlapping HI-union chain pairs
cat hiunion_memo_residues_combined.tsv intact_memo_residues_combined.tsv > memo_residues_combined.tsv # all chain pairs, need to use for IntAct

# (11) Run the following python scripts to write interfacial residues found in *_memo_residues_combined.tsv to file
cd ..
python3 hiunion_write_interactions.py -o --output_directory
python3 intact_write_interactions.py -o --output_directory

# (12) Create a dictionary storing PDB auth_seq_id to label_seq_id mappings for converting auth residues to label residues
sbatch get_auth_label_dict.slurm # change the -d (--pdb_download_path), -o (--output_directory), and sbatch --mail_user and --account arguments

# (13) MOVE the following files back to the local ../data/processed/interactome directory, and exit out of your server
# 1. hiunion_uniprot_pdb_chains.pickle
# 2. intact_uniprot_pdb_chains.pickle
# 3. hiunion_auth_label_dict.pickle
# 4. intact_auth_label_dict.pickle
# 5. hiunion_interfacial_residues.tsv
# 6. intact_interfacial_residues_physical.tsv
# and the following additional (optional) files for backup...
# 7. hiunion_memo_residues_combined.tsv
# 8. intact_memo_residues_combined.tsv

# **5. Construct the HI-union and IntAct structural interactomes**

# (1) Create final interactions file (keep interactions with >=50% of interfacial residues mapping to their corresponding Uniprot proteins)
python3 build_structural_interactome.py

# (2) Select PDB structural templates for each binary interaction
# threshold for templates: 0.0 < PDB resolution <= 3.5 & BLASTP alignment bitscore >= 50
python3 select_pdb_structural_template.py

# CAN THIS BE DONE HERE? AFTER SELECTING PDB STRUCTURAL TEMPLATES
# WAS ORIGINALLY DONE BEFORE get_modeller_templates_foldx_mutations.py 
# **1. Combine BLASTP alignments from HI-Union and IntAct interactomes**
# the overlapping uniprot, pdb chain pairs should have the same BLASTP alignment
# will need this to run MODELLER (also reduces the number of computations needed)
python3 compare_hiunion_intact_blast_alignments.py


# ----------BUILD COMBINED STRUCTURAL INTERACTOME USING MODELLER----------

# **1. Combine HI-union and IntAct structural templates**
python3 get_combined_structural_templates.py

cd ../data/processed
mkdir alignments
mkdir pdb_cif
mkdir foldx_pdb
mkdir mutations
mkdir mutations_final
cd pdb_cif
mkdir modified_cif
# make sure that data/processed/pdb_cif/modified_cif within pdb_cif is CASE SENSITIVE 
# otherwise chains such as S and s will be confused as the same chain
# if using Windows OS, run Windows PowerShell as administrator
# and type the following (replace <path> with the absolute path of your modified_cif directory)
fsutil.exe file setCaseSensitiveInfo <path> enable

# NEED TO REWRITE THE FOLLOWING SCRIPTS TO RUN ON COMPUTE CANADA
# (1) edit_cif_files_create_heterodimer_ali_files.py
# (2) run_modeller_heterodimer.py
# and subsequent analyses/files to copy back

# **2. Create files needed for running MODELLER**
cd ../../../../scripts
python3 edit_cif_files_create_heterodimer_ali_files.py

# **3. Download and run MODELLER** 
# NOTE: Please use a server for this step

# (1) Install MODELLER (https://salilab.org/modeller/10.3/release.html#deb, choose either with conda or Linux distributions to run on compute canada)
# if you have Anaconda, you can install MODELLER as follows:
conda install -c salilab modeller

# (2) Make sure that your 'modeller' folder is CASE SENSITIVE 
# otherwise chains such as S and s will be confused as the same chain
cd ../data/processed
mkdir modeller
# if using Windows OS, run Windows PowerShell as administrator
# and type the following (replace <path> with the absolute path of your 'modeller' directory)
fsutil.exe file setCaseSensitiveInfo <path> enable

# (3) Run MODELLER
# if any uniprot alignments have residue indices over 9999
# take % 10000 of the indices (since .pdb files can't take over 4 chars)
python3 run_modeller_heterodimer.py

# (4) Remove all intermediate MODELLER files
# and move to ../data/processed/modeller
rm *.ini *.rsr *.sch *.D00000001 *.V99990001
mv *.pdb ../data/processed/modeller

# (5) Check for uncompleted models and remove them
# in ../data/processed/interactome/uncompleted_models.pickle
# these are uncompleted because of discrepancies between PDB SeqRes and mmCIF...
# can choose to fix .cif and/or .ali file, but since there are only 9, decided to just discard
# will pickle/write updated files with the suffix *_modeled*
python3 remove_unmodeled_ppis.py


#-----------PREPARE STRUCTURAL MODELS FOR FOLDX ENERGY CALCULATIONS---------------

# running on compute canada
# on compute canada create the following dir
# ~/projects/def-*/username/foldx
# underneath foldx, need the following directories & files (copy them over from local data dir)
# (1) interactome
# 	(i) all_combined_structural_templates_modeled.pickle
# 	(ii) all_blast_best_alignments.pickle
# (2) modeller 
# 	(i) all MODELLER constructed .B99990001.pdb files
# (3) foldx_pdb

# scripts need to be in ~/projects/def-*/username/foldx (copy them over from local scripts dir)
# (1) create_foldx_repair_pdb_job.py
# (2) simple_tools.py
# (3) find_interfacial_residues_between_modeller_heterodimers.py
# (4) create_repaired_pdb_interfacial_residues_split_slurm_jobs.py
# (5) combine_repaired_pdb_interfacial_residues_pickle_files.py 

# input arguments (please use absolute paths)
# -d (or --data_dir)
# -s (or --script_dir)
# -f (or --foldx_executable_name)
# -n (or --num_splits)
# -a (or --account)
# -e (or --email_address)


# **1. Make sure that your 'foldx_pdb' folder is CASE SENSITIVE** 
# if using Windows OS, run Windows PowerShell as administrator
# and type the following (replace <path> with the absolute path of your 'foldx_pdb' directory)
fsutil.exe file setCaseSensitiveInfo <path> enable

# **2. Create .bash script for running FoldX RepairPDB**
python3 create_foldx_repair_pdb_job.py -d <data_dir> -f <foldx_exceutable_name> -n <num_splits> -a <account> -e <email_address>

# **3. Download FoldX executable file**
# download FoldX in the foldx/foldx_pdb dir
# install FoldX (download from the website, need to sign up for an account to get a license), should get a .zip file (e.g. foldx5_1Linux64_0.zip, if using FoldX version 5.1)
# there should get an executable file called foldx_202x1231 (x is the year, e.g. foldx_20251231) inside the .zip file, this is needed as a parameter for the next step

# **4. Run FoldX RepairPDB**
# copy and unzip FoldX .zip file to each split_* directory
# run RepairPDB to repair MODELLER constructed homology models
# needs to be repaired as small issues may change mutagenesis energy calculations done by FoldX
cd ~/projects/def-*/username/foldx/foldx_pdb
for i in {0..n} # n = num_splits-1
do
	cp foldx5_1Linux64_0.zip split_$i
	echo copied zipped FoldX file to split_$i
	unzip split_$i/foldx5_1Linux64_0.zip -d split_$i
	echo extracted files zipped FoldX file to split_$i
	cd split_$i
	sbatch repair_pdb_$i.slurm
	echo submitted .slurm file for split_$i
	cd ../
done

# **4. Move all repaired PDB files to a single directory**
cd ~/projects/def-*/username/foldx
mkdir repaired_pdb
cd ~/projects/def-yxia/tingyisu/foldx/foldx_pdb
for i in {0..n} # n = num_splits-1
do
	for file in split_$i/*_Repair.pdb
	do
		mv $file ../repaired_pdb
	done
done

# **5. Activate virtural environment (on compute canada) if necessary to be able to use Biopython**
python_load
env_activate

# **6. Find interfacial residues between interacting chains of repaired homology models**
# uses find_interfacial_residues_between_modeller_heterodimers.py
python3 create_repaired_pdb_interfacial_residues_split_slurm_jobs.py -s <script_dir> -n <num_splits> -a <account> -e <email_address>

# **7. Submit SLURM jobs for finding interfacial residues between interacting chains of repaired homology models**
for i in {0..n} # n = num_splits-1
do
	sbatch interfacial_residues_$i.slurm
done

# **8. Combine file outputs of all SLURM jobs into a single .tsv file**
cd ~/projects/def-*/username/foldx/interactome
cat modeller_heterodimer_structural_interactome_0.tsv > modeller_heterodimer_structural_interactome.tsv
for i in {1..n} # n = num_splits-1
do
	awk FNR!=1 modeller_heterodimer_structural_interactome_$i.tsv >> modeller_heterodimer_structural_interactome.tsv 
done

# **9. Combine pickle files across all SLURM jobs into a single pickle file**
cd ../
python3 combine_repaired_pdb_interfacial_residues_pickle_files.py -s <script_dir> -n <num_splits>


# **10. Copy remote files back to local dir**
# copy the following files back to local directory
# under foldx/interactome to local ../data/processed/interactome
# (1) modeller_heterodimer_interfacial_residues.pickle
# (2) modeller_heterodimer_structural_interactome.tsv
# (3) all files in repaired_pdb (if have enough space?)

# **11. (OPTIONAL: for sanity check purposes)**
# check and compare interfacial residues from structural templates vs. from modeller heterodimers (after undergoing FoldX RepairPDB)
# save PPIs that lost interfacial residues after modelling
python3 compare_structural_templates_modeller_heterodimers_interfacial_residues.py

# ----------PROCESS AND MAP MISSENSE AND NONSENSE MUTATIONS ONTO STRUCTURAL TEMPLATES----------

# **1. Process ClinVar mutations**
python3 process_clinvar_mutations.py

# **2. Download and process dbSNP mutations on your server**
# NOTE: need to run on a server because dbSNP data is very large (each chromosome file is tens of GB)

# (1) Create a 'dbsnp' directory on your server; all downloaded & processed dbSNP data will be stored here
# If you're using Compute Canada, create the directory your scratch directory as the dbSNP data is hundreds of GB total
# replace, 'username' with your compute canada username
cd /home/username/scratch/
mkdir dbsnp
cd dbsnp

# (2) Copy the following script files from your local script directory to your 'dbsnp' directory on your server
# 1. download_dbsnp_json.slurm
# 2. process_dbsnp_json.py
# 3. simple_tools.py 

# (3) Copy all 'dbsnp_json_*.slurm' files under the local script/dbsnp_json dir to your 'dbsnp' dir on your server

# (4) Create a new 'data' dir under the 'dbsnp' dir on your server and copy the following data files from your local data/processed/ dir:
mkdir data
cd data
# 1. swissprot_ids_list.pickle
# 2. refseq_prot_gene_name_uniprot_dict.pickle
# 3. refseq_prot_seq_dict.pickle

# (5) Download dbSNP data in .json format (one .json file for each chromosome) 
# (need to change the directory names in the .slurm job)
cd ../
sbatch download_dbsnp_json.slurm # change the sbatch --mail_user and --account arguments

# (6) Process the dbSNP .json files
# Change the the -d (--data_dir) and sbatch --mail_user and --account arguments in each .slurm job
# These jobs will process the dbSNP .json files (4 chromsomes per job)
for i in {1..6}
do
	sbatch dbsnp_json_"$i".slurm
done

# (7) Combine the dbSNP output files
# missense mutations
cat data/dbsnp_missense_mutations_chr1.tsv > data/dbsnp_missense_mutations_all.tsv # want the header line at the beginning of the file
# nonsense mutations
cat data/dbsnp_nonsense_mutations_chr1.tsv > data/dbsnp_nonsense_mutations_all.tsv # want the header line at the beginning of the file
# nonstop mutations
cat data/dbsnp_nonstop_mutations_chr1.tsv > data/dbsnp_nonstop_mutations_all.tsv # want the header line at the beginning of the file
chroms=('2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X' 'Y')
for i in "${chroms[@]}"
do
	# missense mutations
	awk FNR!=1 data/dbsnp_missense_mutations_chr"$i".tsv >> data/dbsnp_missense_mutations_all.tsv # skip header line
	# nonsense mutations
	awk FNR!=1 data/dbsnp_nonsense_mutations_chr"$i".tsv >> data/dbsnp_nonsense_mutations_all.tsv # skip header line
	# nonstop mutations
	awk FNR!=1 data/dbsnp_nonstop_mutations_chr"$i".tsv >> data/dbsnp_nonstop_mutations_all.tsv # skip header line
done

# (8) Copy the output files back to your local data/processed/mutations directory
# create a 'dbsnp_unselected' folder under data/processed/mutations
cd data/processed/mutations
mkdir dbsnp_unselected
cd dbsnp_unselected
# move the *all.tsv files on compute canada to local dbsnp_unselected dir
# use scp to move; replace username with your compute canada username
scp username@graham.computecanada.ca:/home/username/scratch/dbsnp/data/*all.tsv .

# (9) Remove ambiguous dbSNP mutations
# i.e. mutations where the HGVSC (nucleotide change on mRNA/cDNA) doesn't correspond with the HGVSG (nucleotide change on chromosome)
# mutations are written to data/processed/mutations
python3 remove_ambiguous_dbsnp_mutations.py

# **3. Process COSMIC genome screen mutations**
# back to compute canada
# process COSMIC genome screen mutations
bash process_cosmic_genome_mutations.bash # uses process_genome_screen_mutations.py 

# **4. Process COSMIC mutations in CGC (Cancer Gene Census)**
python3 process_cosmic_cgc_mutations.py

# **5. Map dbSNP, ClinVar, and COSMIC (both genome screen and CGC) mutations**
# onto the UniProt proteins in the two structural interactomes (HI-union-SI and IntAct-SI)
python3 map_mutation_flanking_seq_to_uniprot.py

# **5. Remove redundant mutations**
# i.e. keep only one mutation with the same amino acid change (even if diff nucleotide change) at a given position on a UniProt protein
python3 remove_redundant_mutations.py

# **6. Find mutations that lie on interfacial residues (IR mutations)**
# do for both interfacial residues between structural templates and
# for interfacial residues between MODELLER constructed homology models (after repairing using RepairPDB)
python3 get_ir_mutations.py

# **6. Prepare mutations for running FoldX PSSM and BuildModel**
python3 get_foldx_mutations.py

# ----------CREATE CONFIGURATION FILES AND RUN FOLDX PSSM AND BUILDMODEL (ENERGY CALCULATIONS)----------
# go back to compute canada
# copy the following files over
# create a 'files' directory in your projects folder
cd ~/projects/def-*/username/foldx
mkdir files
# 1. Copy the following files in (/data/processed/edgotypes) onto compute canada under the 'files' (e.g. ~/projects/def-*/username/foldx/files) folder:
#	(1) all_foldx_pssm_mutations.pickle
#	(2) all_foldx_buildmodel_mutations.pickle
#	(3) *_ir_homology_models_with_modeller_foldx_mutations.tsv
#	(4) *_ir_structural_templates_with_modeller_foldx_mutations.tsv
# ***is essentially all of the files in /data/processed/edgotypes, so can just do the following
cd ../data/processed/edgotypes
scp * username@cedar.computecanada.ca:/home/username/projects/def-*/username/foldx/files

# input arguments (please use absolute paths)
# -s (or --script_dir)
# -c (or --scratch_dir)
# -a (or --account)
# -f (or --foldx_executable_name)

# **1. Prepare FoldX configuration files and SLURM jobs**
python3 prepare_foldx_pssm_files.py -s <script_dir> -c <scratch_dir> -a <account> -f <foldx_executable_name>

# **2. Run FoldX PSSM**
# copy foldx5_1Linux64_0.zip to the created PSSM dir in the scratch directory (/home/username/scratch/foldx_pssm_all)
# then submit each FoldX PSSM SLURM job
cd /home/username/scratch/foldx_pssm_all
for i in {0..28} # n = num_splits-1
do
	cp foldx5_1Linux64_0.zip split_$i
	echo copied compressed FoldX5.1 file to split_$i
	unzip split_$i/foldx5_1Linux64_0.zip -d split_$i
	echo extracted files from compressed FoldX5.1 file to split_$i
	cd split_$i
	sbatch foldx_pssm_$i.slurm
	echo submitted .slurm file for split_$i
	cd ../
done

# **3. Move FoldX PSSM outputs (binding DDG) out of scratch dir**
# if don't have much space in scratch dir, find edgetic mutations first
# need to combine binding DDG files together into a folder
bash move_binding_ddg_files.bash

# **4. Save (Pickle) FoldX PSSM outputs (binding DDG)**
# pickle all binding DDG calculations into a dictionary (key = FoldX PSSM mutation, value = binding DDG)
# will make the next step (finding edgetic mutations) much faster
python3 save_append_to_dict.py -f /home/username/scratch/binding_ddg -p /home/username/projects/def-*/username/foldx/files/foldx_pssm_binding_ddg_dict.pickle -d binding

# **5. Find edgetic (binding-destabilizing) mutations (using binding DDGs of IR mutations)**
# find edgetic mutations
python get_edgetic_mutations.py -s <script_dir> -c <scratch_dir>

# **6. Remove FoldX PSSM intermediate/output files to save space**
# need to rm foldx_pssm_all to save space
rm -r /home/username/scratch/foldx_pssm_all

# **7. Rrepare FoldX BuildModel configuration files and SLURM jobs**
sbatch prepare_foldx_buildmodel_files.slurm # uses prepare_foldx_buildmodel_files.py, need to change the dir arguments in the file

# **8. Run FoldX BuildModel**
# copy foldx5_1Linux64_0.zip to the created BuildModel dir in the scratch directory (/home/username/scratch/foldx_buildmodel_all)
# then submit each FoldX BuildModel SLURM job
cd /home/username/scratch/foldx_buildmodel_all
for i in {0..n} # n = num_splits-1
do
	cp foldx5_1Linux64_0.zip split_$i
	echo copied compressed FoldX5.1 file to split_$i
	unzip split_$i/foldx5_1Linux64_0.zip -d split_$i
	echo extracted files from compressed FoldX5.1 file to split_$i
	cd split_$i
	sbatch foldx_buildmodel_$i.slurm
	echo submitted .slurm file for split_$i
	cd ../
done

# **9. Check whether FoldX BuildModel jobs were completed successfully**
# if jobs were completed successfully
# should print 2000 for each split, other than the last split (which can range from 0-2000 mutations)
bash check_buildmodel_jobs.bash

# **10. Move FoldX BuildModel outputs (folding DDG) out of scratch dir**
# and into their own separate directory
bash move_folding_ddg_files.bash

# **11. Save (Pickle) FoldX BuildModel outputs (folding DDG)**
# pickle all folding DDG calculations into a dictionary (key = FoldX BuildModel mutation, value = folding DDG)
# will make the next step (finding quasi-null/quasi-wildtype mutations) much faster
python3 save_append_to_dict.py -f /home/username/scratch/folding_ddg -p /home/username/projects/def-*/username/foldx/files/foldx_buildmodel_folding_ddg_dict.pickle -d folding

# **12. Find quasi-null (folding-destabilizing) and quasi-wildtype (non-destabilizing) mutations based on folding DDG**
# categorize non-edgetic mutations into quasi-null or quasi-wildtype based on folding DDG
# copying files takes a long time, so submit job as a .slurm file
python3 get_quasi_null_wildtype.py -s <script_dir>

# **13. Remove FoldX BuildModel intermediate/output files to save space**
# need to rm foldx_buildmodel_all to save space
rm -r /home/username/scratch/foldx_buildmodel_all

# **14. Copy outputs back to local dir**
# copy the following files from ~/projects/def-*/username/foldx/files to the local data/processed/edgotypes dir
# (1) all_foldx_buildmodel_mutations_combined.pickle
# (2) additional_foldx_buildmodel_mutations.pickle
# (3) foldx_pssm_binding_ddg_dict.pickle
# (4) foldx_buildmodel_folding_ddg_dict.pickle
# (5) *_mutation_edgotypes_quasi_null_wildtype.tsv
# (6) *_mutation_edgotypes.tsv 
# back to local dir...

# **15. Print edgotype numbers (# of edgetic, quasi-null, and quasi-wildtyp mutations) and stats**
python3 get_edgotype_quasi_null_wildtype_nums.py

# **16. Map nonsense mutations onto structural interactomes**
python3 map_nonsense_mutations_onto_si.py