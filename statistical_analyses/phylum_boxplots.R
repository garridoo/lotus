
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de
#

options(warn=-1)

# cleanup

rm(list=ls())

# load plotting functions

source("plotting_functions.R")
library(reshape)

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

# remove non-bacterial and Chloroflexi OTUs

taxonomy <- taxonomy[taxonomy[, 2]=="Bacteria", ]
taxonomy <- taxonomy[taxonomy[, 3]!="Chloroflexi", ]

idx <- rownames(otu_table) %in% taxonomy[, 1]
otu_table <- otu_table[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

# normalize OTU table

otu_table_norm <- apply(otu_table, 2, function(x) x / sum(x))

# aggregate relative abundances to the phylum level

tax <- taxonomy[match(rownames(otu_table), taxonomy[, 1]), 3]
tax <- as.character(tax)
idx <- !tax %in% c("Proteobacteria", "Actinobacteria", "Bacteroidetes", "Firmicutes")
tax[idx] <- "Other"

l <- unique(tax)

tax_table <- matrix(nrow=length(l), ncol=dim(otu_table)[2])
colnames(tax_table) <- colnames(otu_table)
rownames(tax_table) <- l

for (i in 1:length(l)) {
    
    tax_table[i, ] <- colSums(otu_table_norm[tax==l[i], ])

}

# generate data frame for plotting

design$genotype <- as.character(design$genotype)
idx <- design$genotype %in% c("hit1_1", "nfr5_2", "nfr5_3", "nin2")
design$genotype[idx] <- "mutant"

design <- design[match(colnames(tax_table), design$SampleID), ]
df <- melt(cbind(design[, 1:3], t(tax_table)))

df$variable <- factor(df$variable, levels=rev(c("Proteobacteria", "Actinobacteria",
                                                "Bacteroidetes", "Firmicutes", "Other")))

df <- df[df$genotype %in% c("gifu", "mutant", "soil"), ]

# plot boxplots for root samples

df_root <- df[df$compartment %in% c("root"), ]

p1 <- ggplot(df_root, aes(x=variable, y=value, fill=genotype)) +
             geom_boxplot(alpha=1, outlier.size=0, size=0.5, width=.8,
                          position=position_dodge(width=0.7)) +
             scale_fill_manual(values=c("transparent", "darkgrey")) +
             scale_y_continuous(labels=percent, limits=c(0, 0.8)) +
             labs(x="", y="relative abundance") +
             coord_flip() +
             main_theme

ggsave(paste(figures.dir, "phylum_boxplots_root.pdf", sep=""), p1, height=4, width=5.5)

# plot boxplots for rhizosphere samples

df_rhizosphere <- df[df$compartment %in% c("rhizosphere"), ]

p2 <- ggplot(df_rhizosphere, aes(x=variable, y=value, fill=genotype)) +
             geom_boxplot(alpha=1, outlier.size=0, size=0.5, width=0.8,
                          position=position_dodge(width=0.7)) +
             scale_fill_manual(values=c("transparent", "darkgrey")) +
             scale_y_continuous(labels=percent, limits=c(0, 0.8)) +
             labs(x="", y="relative abundance") +
             coord_flip() +
             main_theme

ggsave(paste(figures.dir, "phylum_boxplots_rhizosphere.pdf", sep=""), p2, height=4, width=5.5)

# plot boxplots for soil samples

df_soil <- df[df$compartment %in% c("soil"), ]

p3 <- ggplot(df_soil, aes(x=variable, y=value, fill=genotype)) +
             geom_boxplot(alpha=1, outlier.size=0, size=0.5, width=0.8,
                          position=position_dodge(width=0.7)) +
             scale_fill_manual(values=c("transparent", "darkgrey")) +
             scale_y_continuous(labels=percent, limits=c(0, 0.8)) +
             labs(x="", y="relative abundance") +
             coord_flip() +
             main_theme

ggsave(paste(figures.dir, "phylum_boxplots_soil.pdf", sep=""), p3, height=4, width=5.5)

