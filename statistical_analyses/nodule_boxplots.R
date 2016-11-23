
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de
#

options(warn=-1)

# cleanup

rm(list=ls())

# load libraries

library(reshape)

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

# remove non-bacterial and Chloroflexi OTUs

taxonomy <- taxonomy[taxonomy[, 2]=="Bacteria", ]
taxonomy <- taxonomy[taxonomy[, 3]!="Chloroflexi", ]

idx <- rownames(otu_table) %in% taxonomy[, 1]
otu_table <- otu_table[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

# normalize OTU table

otu_table_norm <- apply(otu_table, 2, function(x) x / sum(x))

### order level

# aggregate relative abundances to the order level

tax <- taxonomy[match(rownames(otu_table), taxonomy[, 1]), 5]
tax <- as.character(tax)
tax <- gsub(".__", "", tax)

l <- unique(tax)
l <- l[l!=""]
l <- l[l!="1.00"]

tax_table <- matrix(nrow=length(l), ncol=dim(otu_table)[2])
colnames(tax_table) <- colnames(otu_table)
rownames(tax_table) <- l

empty <- rep(0, dim(otu_table)[2])

for (i in 1:length(l)) {
    
    tax_table[i, ] <- colSums(rbind(empty, otu_table_norm[tax==l[i], ]))

}

idx <- design$genotype %in% c("hit1_1", "nfr5_2", "nfr5_3", "nin2", "gifu") &
       design$compartment %in% c("pooled_nodules", "root")

samples <- design$SampleID[idx]
tax_table <- tax_table[, colnames(tax_table) %in% samples]

design <- design[match(colnames(tax_table), design$SampleID), ]
design$genotype <- as.character(design$genotype)
design$genotype[design$genotype %in% c("hit1_1", "nfr5_2", "nfr5_3", "nin2")] <- "mutant"

# rank abundance plots nodules only

idx <- design$genotype=="gifu" &
       design$compartment=="pooled_nodules"
samples <- design$SampleID[idx]
tax_table_nodules <- tax_table[, colnames(tax_table) %in% samples]

design_nodules <- design[match(colnames(tax_table_nodules), design$SampleID), ]

# generate data frame for plotting

df <- melt(cbind(design_nodules[, c(1, 2)], t(tax_table_nodules)))

abundant_orders <- names(sort(apply(tax_table_nodules, 1, mean), decreasing=T)[1:20])
df <- df[df$variable %in% abundant_orders, ]

df$variable <- factor(df$variable, levels=abundant_orders)

p1 <- ggplot(df, aes(x=variable, y=value, fill=genotype)) +
             geom_boxplot(alpha=1, outlier.size=0, size=0.5, width=.8,
                          position=position_dodge(width=0.7)) +
             scale_fill_manual(values=c("transparent", "darkgrey")) +
             scale_y_continuous(labels=percent, limits=c(0, 1)) +
             labs(x="", y="relative abundance") +
             main_theme +
             theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(paste(figures.dir, "order_boxplots_nodule.pdf", sep=""), p1, height=5.5, width=5)

# boxplots rhizobiales

tax <- taxonomy[taxonomy[, 5]=="Rhizobiales", ]
idx <- tax[, 7]=="Mesorhizobium"
# idx <- tax[, 6]=="Phyllobacteriaceae"
rhizobiales_otus <- tax[, 1]
meso_otus <- tax[idx, 1]
non_meso_otus <- tax[!idx, 1]

# generate data frame for plotting

idx <- colnames(otu_table_norm) %in% design$SampleID
otu_table_norm <- otu_table_norm[, idx]

meso_ra <- colSums(otu_table_norm[rownames(otu_table_norm) %in% meso_otus, ])
non_meso_ra <- colSums(otu_table_norm[rownames(otu_table_norm) %in% non_meso_otus, ])
rhizobiales_ra <- rbind(meso_ra, non_meso_ra)

df <- melt(cbind(design[, c(1, 2, 3)], t(rhizobiales_ra)))

df$group[df$genotype=="gifu" & df$compartment=="root"] <- "gifu_root"
df$group[df$genotype=="mutant" & df$compartment=="root"] <- "mutant_root"
df$group[df$genotype=="gifu" & df$compartment=="pooled_nodules"] <- "gifu_nodules"
df <- df[!(df$genotype=="mutant" & df$compartment=="pooled_nodules"),]

df$group <- factor(df$group, levels=c("gifu_root", "mutant_root", "gifu_nodules"))

p2 <- ggplot(df, aes(x=group, y=value, fill=variable)) +
             geom_boxplot(alpha=1, outlier.size=0, size=0.5, width=.7,
                          position=position_dodge(width=0.6)) +
             scale_fill_manual(values=c("darkgrey", "transparent")) +
             scale_y_continuous(labels=percent, limits=c(0, 1)) +
             labs(x="", y="relative abundance") +
             main_theme +
             theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(paste(figures.dir, "rhizobiales_boxplots.pdf", sep=""), p2, height=5.5, width=5)

