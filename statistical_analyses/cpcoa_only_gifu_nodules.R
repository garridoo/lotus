
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de
#

options(warn=-1)

# cleanup

rm(list=ls())

# load plotting functions

source("plotting_functions.R")
source("cpcoa.func.R")

# load plotting functions

library("ggplot2")
library("scales")
library("grid")
library("vegan")

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

# subset samples of interest

idx <- design$compartment %in% c("soil", "root", "rhizosphere", "pooled_nodules") &
       design$genotype %in% c("soil", "gifu")
samples <- design$SampleID[idx]
idx <- colnames(otu_table) %in% samples
otu_table <- otu_table[, idx]
design <- design[idx, ]

# subset samples of interest from distance matrices

idx <- rownames(uw_unifrac) %in% colnames(otu_table)

uw_unifrac <- uw_unifrac[idx, idx]
w_unifrac <- w_unifrac[idx, idx]
bray_curtis <- bray_curtis[idx, idx]

### CPCoA

colors <- data.frame(group=c("pooled_nodules", "rhizosphere", "root", "soil"),
                     color=c(c_red, c_dark_red, c_very_dark_green, c_dark_brown))

shapes <- data.frame(group=c("gifu", "soil"),
                     shape=c(19, 3))

sqrt_transform <- T

d <- design

capscale.gen <- capscale(t(otu_table) ~ compartment + Condition(soil), data=d, add=F, sqrt.dist=sqrt_transform, distance="bray")

# ANOVA-like permutation analysis

perm_anova.gen <- anova.cca(capscale.gen)
print(perm_anova.gen)
                                                    
# generate variability tables and calculate confidence intervals for the variance

var_tbl.gen <- variability_table(capscale.gen)

eig <- capscale.gen$CCA$eig

variance <- var_tbl.gen["constrained", "proportion"]
p.val <- perm_anova.gen[1, 4]

# extract the weighted average (sample) scores

points <- capscale.gen$CCA$wa[, 1:2]
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design[match(rownames(points), design$SampleID), ])

# plot CPCo 1 and 2

p <- ggplot(points, aes(x=x, y=y, color=compartment, shape=genotype)) +
     geom_point(alpha=.7, size=1.5) +
     scale_colour_manual(values=as.character(colors$color)) +
     scale_shape_manual(values=shapes$shape)+
     labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
     y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
     ggtitle(paste(format(100 * variance, digits=3), " % of variance; p=",
                   format(p.val, digits=2),
                   sep="")) +
     main_theme +
     theme(legend.position="top")
plot(p)
ggsave(paste(figures.dir, "CPCoA_BC_only_gifu_nodules.pdf", sep=""), p)

