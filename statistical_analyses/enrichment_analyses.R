
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de
#

options(warn=-1)

rm(list=ls())

# load plotting functions

source("plotting_functions.R")

# directories

results.dir <- "/biodata/dep_psl/grp_psl/garridoo/lotus/454/results/all/"
figures.dir <- "/biodata/dep_psl/grp_psl/garridoo/lotus/454/figures/all/"

# files

design.file <- paste(results.dir, "design.txt", sep="")
taxonomy.file <- paste(results.dir, "taxonomy.txt", sep="")
otu_table.file <- paste(results.dir, "otu_table.txt", sep="")

# load data

design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, check.names=F)
taxonomy <- read.table(taxonomy.file, sep="\t", header=F, fill=T)

# re-order data matrices

idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

# remove non-bacterial (Chloroplast and Mitochondria)
# as well as Chloroflexi OTUs

taxonomy <- taxonomy[taxonomy[, 2]=="Bacteria", ]
taxonomy <- taxonomy[taxonomy[, 3]!="Chloroflexi", ]

idx <- rownames(otu_table) %in% taxonomy[, 1]
otu_table <- otu_table[idx, ]

# total sum normalization

otu_table_norm <- apply(otu_table, 2, function(x) x / sum(x))

# thresholding

threshold <- .05
idx <- rowSums(otu_table_norm * 100 > threshold) >= 1
otu_table <- otu_table[idx, ]
otu_table_norm <- otu_table_norm[idx, ]

# parameters

scale <- 1000               # scale for the ternary plots
alpha <- 0.05               # significance threshold
p.adj.method <- "fdr"       # FDR p-value adjustment method

# ternary plots - all soil batches

source("gifu_root_rhizo_soil.R")
source("mutant_root_rhizo_soil.R")

# ternary plots - CAS8 soil

source("gifu_root_rhizo_soil_cas8.R")
source("nfr5_2_root_rhizo_soil_cas8.R")

# ternary plots - CAS9 soil

source("gifu_root_rhizo_soil_cas9.R")
source("nfr5_2_root_rhizo_soil_cas9.R")

# ternary plots - CAS10 soil

source("gifu_root_rhizo_soil_cas10.R")
source("hit1_1_root_rhizo_soil_cas10.R")
source("nfr5_3_root_rhizo_soil_cas10.R")
source("nin2_root_rhizo_soil_cas10.R")

# WT v. mutant heatmaps

source("gifu_v_mutant_heatmaps.R")

# WT and mutant root and rhizpshere Manhattan plots

source("manhattan_gifu.R")
source("manhattan_gifu.R")

# WT v. mutant Manhattan plots

source("manhattan_gifu_v_mutant_root.R")
source("manhattan_gifu_v_mutant_rhizo.R")

