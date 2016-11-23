
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de
#

options(warn=-1)

# cleanup

rm(list=ls())

# load plotting functions

source("plotting_functions.R")

# directories

results.dir <- "/biodata/dep_psl/grp_psl/garridoo/lotus/454/results/all/"
figures.dir <- "/biodata/dep_psl/grp_psl/garridoo/lotus/454/figures/all/"

# files

design.file <- paste(results.dir, "design.txt", sep="")
taxonomy.file <- paste(results.dir, "taxonomy.txt", sep="")
otu_table.file <- paste(results.dir, "otu_table_norm.txt", sep="")

uw_unifrac.file <- paste(results.dir, "unweighted_unifrac.txt", sep="")
w_unifrac.file <- paste(results.dir, "weighted_unifrac.txt", sep="")
bray_curtis.file <- paste(results.dir, "bray_curtis.txt", sep="")

shannon.file <- paste(results.dir, "shannon.txt", sep="")
chao.file <- paste(results.dir, "chao.txt", sep="")
observed_otus.file <- paste(results.dir, "observed_otus.txt", sep="")

# load data

design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, check.names=F)
taxonomy <- read.table(taxonomy.file, sep="\t", header=F, fill=T)

uw_unifrac <- read.table(uw_unifrac.file, sep="\t", header=T, check.names=F)
w_unifrac <- read.table(w_unifrac.file, sep="\t", header=T, check.names=F)
bray_curtis <- read.table(bray_curtis.file, sep="\t", header=T, check.names=F)

shannon <- read.table(shannon.file, sep="\t", header=T, check.names=F)
chao <- read.table(chao.file, sep="\t", header=T, check.names=F)
observed_otus <- read.table(observed_otus.file, sep="\t", header=T, check.names=F)

# re-order data matrices

idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[idx, idx]
uw_unifrac <- uw_unifrac[idx, idx]
w_unifrac <- w_unifrac[idx, idx]
bray_curtis <- bray_curtis[idx, idx]

# remove non-bacterial and Chloroflexi OTUs

taxonomy <- taxonomy[taxonomy[, 2]=="k__Bacteria", ]
taxonomy <- taxonomy[taxonomy[, 3]!="p__Chloroflexi", ]

idx <- rownames(otu_table) %in% taxonomy[, 1]
otu_table <- otu_table[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

# remove negative control samples

neg_control_samples <- design$SampleID[design$compartment=="negative_control"]
idx <- !colnames(otu_table) %in% neg_control_samples
otu_table <- otu_table[, idx]
design <- design[idx, ]

# remove individual nodules samples

ind_nod_samples <- design$SampleID[design$compartment=="individual_nodule"]
idx <- !colnames(otu_table) %in% ind_nod_samples
otu_table <- otu_table[, idx]
design <- design[idx, ]
 
# subset samples of interest from distance matrices

idx <- rownames(uw_unifrac) %in% colnames(otu_table)

uw_unifrac <- uw_unifrac[idx, idx]
w_unifrac <- w_unifrac[idx, idx]
bray_curtis <- bray_curtis[idx, idx]

idx <- shannon[, 1] %in% design$SampleID

shannon <- shannon[idx, ]
chao <- chao[idx, ]
observed_otus <- observed_otus[idx, ]

### alpha diversity

colors <- data.frame(group=c("pooled_nodules", "root", "rhizosphere", "soil"),
                     color=c(c_red, c_very_dark_green, c_dark_red, c_dark_brown))

shapes <- data.frame(group=c("gifu", "hit1_1", "nfr5_2", "nfr5_3", "nin2", "soil"),
                     shape=c(19, 0, 24, 25, 18, 3))

# shannon index

index <- cbind(shannon[, 2], design[match(shannon[, 1], design$SampleID), ])
colnames(index)[1] <- "value"

index$genotype <- factor(index$genotype, levels=shapes$group)
index$compartment <- factor(index$compartment, levels=colors$group)

# reorder boxplots

l <- c("soil", "rhizosphere", "root", "pooled_nodules")
index$compartment <- factor(index$compartment, levels=l)
colors <- colors[match(l, colors$group), ]

p <- ggplot(index, aes(x=compartment, y=value, color=compartment)) +
            geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
            geom_jitter(aes(shape=genotype), position=position_jitter(0.17), size=1, alpha=0.7) +
            scale_colour_manual(values=as.character(colors$color)) +
            scale_shape_manual(values=shapes$shape) +
            labs(x="", y="shannon index") +
            main_theme

ggsave(paste(figures.dir, "shannon.pdf", sep=""), p)

# chao index

index <- cbind(chao[, 2], design[match(chao[, 1], design$SampleID), ])
colnames(index)[1] <- "value"

index$genotype <- factor(index$genotype, levels=shapes$group)
index$compartment <- factor(index$compartment, levels=colors$group)

# reorder boxplots

l <- c("soil", "rhizosphere", "root", "pooled_nodules")
index$compartment <- factor(index$compartment, levels=l)
colors <- colors[match(l, colors$group), ]

p <- ggplot(index, aes(x=compartment, y=value, color=compartment)) +
            geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
            geom_jitter(aes(shape=genotype), position=position_jitter(0.17), size=1, alpha=0.7) +
            scale_colour_manual(values=as.character(colors$color)) +
            scale_shape_manual(values=shapes$shape) +
            labs(x="", y="chao index") +
            main_theme

ggsave(paste(figures.dir, "chao.pdf", sep=""), p)

# observed species

index <- cbind(observed_otus[, 2], design[match(observed_otus[, 1], design$SampleID), ])
colnames(index)[1] <- "value"

index$genotype <- factor(index$genotype, levels=shapes$group)
index$compartment <- factor(index$compartment, levels=colors$group)

# reorder boxplots

l <- c("soil", "rhizosphere", "root", "pooled_nodules")
index$compartment <- factor(index$compartment, levels=l)
colors <- colors[match(l, colors$group), ]

p <- ggplot(index, aes(x=compartment, y=value, color=compartment)) +
            geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
            geom_jitter(aes(shape=genotype), position=position_jitter(0.17), size=1, alpha=0.7) +
            scale_colour_manual(values=as.character(colors$color)) +
            scale_shape_manual(values=shapes$shape) +
            labs(x="", y="observed OTUs") +
            main_theme

ggsave(paste(figures.dir, "observed_OTUs.pdf", sep=""), p)

