
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

# re-order data matrices

idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]
uw_unifrac <- uw_unifrac[idx, idx]
w_unifrac <- w_unifrac[idx, idx]
bray_curtis <- bray_curtis[idx, idx]

# remove non-bacterial and Chloroflexi OTUs

taxonomy <- taxonomy[taxonomy[, 2]=="Bacteria", ]
taxonomy <- taxonomy[taxonomy[, 3]!="Chloroflexi", ]

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
 
# subset soil batch

soil <- "CAS10"

soil_batch_samples <- design$SampleID[design$soil==soil]
idx <- colnames(otu_table) %in% soil_batch_samples
otu_table <- otu_table[, idx]
design <- design[idx, ]
soil <- gsub(" ", "_", Reduce(x=soil, f=paste))
           
# subset samples of interest from distance matrices

idx <- rownames(uw_unifrac) %in% colnames(otu_table)

uw_unifrac <- uw_unifrac[idx, idx]
w_unifrac <- w_unifrac[idx, idx]
bray_curtis <- bray_curtis[idx, idx]

### beta diversity

colors <- data.frame(group=c("pooled_nodules", "rhizosphere", "root", "soil"),
                     color=c("transparent", c_dark_red, c_very_dark_green, c_dark_brown))

shapes <- data.frame(group=c("gifu", "hit1_1", "nfr5_2", "nfr5_3", "nin2", "soil"),
                     shape=c(19, 0, 24, 25, 18, 3))

colors <- colors[colors$group %in% design$compartment, ]
shapes <- shapes[shapes$group %in% design$genotype, ]

# PCoA Bray-Curtis

k <- 2
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design[match(rownames(points), design$SampleID), ])

points$genotype <- factor(points$genotype, levels=shapes$group)
points$compartment <- factor(points$compartment, levels=colors$group)

# plot PCo 1 and 2

p <- ggplot(points, aes(x=x, y=y, color=compartment, shape=genotype)) +
     geom_point(alpha=.7, size=2) +
     scale_colour_manual(values=as.character(colors$color)) +
     scale_shape_manual(values=shapes$shape) +
     labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
     y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
     main_theme +
     theme(legend.position="top")

ggsave(paste(figures.dir, "PCoA_BC_", soil, ".pdf", sep=""), p)

