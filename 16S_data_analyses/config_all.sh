#!/bin/bash

# scripts for 16S data analysis
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# list of library IDs, separated by single spaces
l_list="1018_A 1018_B 1018_C 1018_D 1018_E 1416A_1A 1416B_1B 1416C_2A 1416D_2B 1416E_3A 1416F_3B"

# working_dir is the folder where (intermediate) results are stored
working_dir="/biodata/dep_psl/grp_psl/garridoo/lotus/454/results/all"

# data_dir is the folder where the raw data (and maping file) are deposited
data_dir="/biodata/dep_psl/grp_psl/garridoo/lotus/454/data"




# paths to relevant reference database files
refdata_dir="/projects/dep_psl/grp_psl/garridoo/16S/ref_data"
gg_core_db=$refdata_dir"/gg_13_8_otus/rep_set/97_otus.fasta"
gg_core_aligned_db=$refdata_dir"/gg_13_8_otus/rep_set_aligned/97_otus.fasta"
gg_taxonomy=$refdata_dir"/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt"
gold_db=$refdata_dir"/cs_gold.fa"

# path to folder containing usearch binary and python scripts
usearch_dir="/projects/dep_psl/grp_psl/garridoo/16S/usearch"

### debug

output=$working_dir"/output.txt"
logfile=$working_dir"/log.txt"

### parameters

# number of threads for parallel steps
n_cores=30
# min Phred score for quality control
phred=30
# min qual score for quality control (454)
qual=25
# trimming length (454)
t_length=323
# max. barcode errors
bc_length=14
# max. barcode errors
bc_err=0
# min. number of reads for a given OTU
min_size=2
# min. number of reads for a given sample
min_cov=500
# identity threshold for OTU clustering
id_threshold=0.97
# min length of alingment required to sequences in the database
min_id_aln=0.75
# subsampling depth used for alpha-diversity
subsampling_depth=1000
# OTU clustering method
cls_method="uclust"

