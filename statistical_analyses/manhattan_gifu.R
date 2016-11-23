
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de
#

library(edgeR)

### gifu - root v. rhizosphere v. soil

# subset samples

idx <- design$compartment %in% c("root", "soil", "rhizosphere") & 
       design$genotype %in% c("gifu", "soil")
design_subset <- design[idx, ]
otu_table_subset <- otu_table[, idx]

### Generalized Linear Model (GLM)

# create DGE list

groups <- design_subset$compartment

d <- DGEList(counts=otu_table_subset, group=groups)
d <- calcNormFactors(d)

# fit the GLM

design.mat <- model.matrix(~ 0 + d$samples$group)
d2 <- estimateGLMCommonDisp(d, design.mat)
d2 <- estimateGLMTagwiseDisp(d2, design.mat)

fit <- glmFit(d2, design.mat)

lrt_rhizo_root <- glmLRT(fit, contrast=c(1, -1, 0))
lrt_rhizo_soil <- glmLRT(fit, contrast=c(1, 0, -1))
lrt_rhizo <- glmLRT(fit, contrast=c(1, -0.5, -0.5))

de_rhizo_root <- decideTestsDGE(lrt_rhizo_root, adjust.method=p.adj.method, p.value=alpha)
de_rhizo_soil <- decideTestsDGE(lrt_rhizo_soil, adjust.method=p.adj.method, p.value=alpha)
de_rhizo <- decideTestsDGE(lrt_rhizo, adjust.method=p.adj.method, p.value=alpha)

rhizo_otus <- rownames(otu_table_subset)[de_rhizo_root==1 & de_rhizo_soil==1]

# rhizo_pvals <- lrt_rhizo_root$table
rhizo_pvals <- lrt_rhizo$table

lrt_root_rhizo <- glmLRT(fit, contrast=c(-1, 1, 0))
lrt_root_soil <- glmLRT(fit, contrast=c(0, 1, -1))
lrt_root <- glmLRT(fit, contrast=c(-0.5, 1, -0.5))

de_root_rhizo <- decideTestsDGE(lrt_root_rhizo, adjust.method=p.adj.method, p.value=alpha)
de_root_soil <- decideTestsDGE(lrt_root_soil, adjust.method=p.adj.method, p.value=alpha)
de_root <- decideTestsDGE(lrt_root, adjust.method=p.adj.method, p.value=alpha)

root_otus <- rownames(otu_table_subset)[de_root_rhizo==1 & de_root_soil==1]

# root_pvals <- lrt_root_rhizo$table
root_pvals <- lrt_root$table

### get mean R.A. in gifu root samples

otu_table_norm_all <- apply(otu_table, 2, function(x) x / sum(x))
 
idx <- design$compartment %in% c("root") &
       design$genotype %in% c("gifu")

root_gifu_samples <- design$SampleID[idx]
otu_table_norm_root_gifu <- otu_table_norm_all[, colnames(otu_table_norm_all) %in% root_gifu_samples]

### get mean R.A. in gifu rhizosphere samples

idx <- design$compartment %in% c("rhizosphere") &
       design$genotype %in% c("gifu")

rhizo_gifu_samples <- design$SampleID[idx]
otu_table_norm_rhizo_gifu <- otu_table_norm_all[, colnames(otu_table_norm_all) %in% rhizo_gifu_samples]

### gifu root-enriched OTUs

### generate data frame for plotting

root_pvals$otu <- rownames(root_pvals)

# P values

root_pvals$neglogp <- -log(root_pvals$PValue)

# enrichment status

root_pvals$enrichment <- root_pvals$otu %in% root_otus

# log fold change

idx <- root_pvals$logFC < 0
root_pvals$neglogp[idx] <- 0

# order OTUs according to taxonomy

taxonomy <- taxonomy[order(taxonomy[, 3], taxonomy[, 4], taxonomy[, 5]), ]
idx <- taxonomy[, 1] %in% root_pvals$otu
taxonomy <- taxonomy[idx, ]

idx <- match(taxonomy[, 1], root_pvals$otu)
root_pvals <- root_pvals[idx, ]

root_pvals$tax <- taxonomy[, 5]

#~ root_pvals$tax <- as.character(root_pvals$tax)
#~ idx <- sort(root_pvals$tax, index.return=T)$ix
#~ root_pvals <- root_pvals[idx, ]

# root_pvals$tax <- as.character(root_pvals$tax)
# idx <- match(taxonomy[, 1], root_pvals$otu)
# root_pvals <- root_pvals[idx, ]

root_pvals$otu <- factor(root_pvals$otu, levels=root_pvals$otu)

# relative abundances

ra <- apply(otu_table_norm_root_gifu, 1, mean)
root_pvals$ra <- ra[match(root_pvals$otu, names(ra))]

# generate vector of colors

colors <- read.table("taxon_colors.txt", sep="\t", header=T)
taxon <- sort(unique(root_pvals$tax))
colors <- data.frame(taxon=taxon, colors=colors$color[match(taxon, colors$taxon)])
colors$colors <- as.character(colors$colors)
colors$colors[is.na(colors$colors)] <- "grey"
enriched_taxa <- unique(taxonomy[taxonomy[, 1] %in% c(root_otus), 5])
colors$colors[!colors$taxon %in% enriched_taxa] <- "grey"

# multiple testing correction thresholds

BF <- -log(0.05 / dim(root_pvals)[1])
FDR <- min(root_pvals$neglogp[root_pvals$enrichment==TRUE])

# plot

p1 <- ggplot(root_pvals, aes(x=otu, y=neglogp, color=tax, size=ra, shape=enrichment)) +
             geom_point(alpha=.8) +
             geom_hline(yintercept=FDR, linetype=2, color="lightgrey") +
             scale_color_manual(values=colors$colors) +
             scale_shape_manual(values=c(21, 19)) +
             scale_size(breaks=c(0, 0.04, 0.08, 0.12, 0.16)) +
             ylim(c(0, 250)) +
             labs(x="OTU", y="-log10(P)") +
             main_theme +
             theme(axis.ticks.x=element_blank(),
                   axis.text.x=element_blank(),
                   legend.position="none")

ggsave(paste(figures.dir, "manhattan_gifu_root.pdf", sep=""), p1, width=10, height=3, useDingbats=F)

### gifu rhizosphere-enriched OTUs

### generate data frame for plotting

rhizo_pvals$otu <- rownames(rhizo_pvals)

# P values

rhizo_pvals$neglogp <- -log(rhizo_pvals$PValue)

# enrichment status

rhizo_pvals$enrichment <- rhizo_pvals$otu %in% rhizo_otus

# log fold change

idx <- rhizo_pvals$logFC < 0
rhizo_pvals$neglogp[idx] <- 0

# order OTUs according to taxonomy

taxonomy <- taxonomy[order(taxonomy[, 3], taxonomy[, 4], taxonomy[, 5]), ]
idx <- taxonomy[, 1] %in% rhizo_pvals$otu
taxonomy <- taxonomy[idx, ]

idx <- match(taxonomy[, 1], rhizo_pvals$otu)
rhizo_pvals <- rhizo_pvals[idx, ]

rhizo_pvals$tax <- taxonomy[, 5]

# rhizo_pvals$tax <- as.character(rhizo_pvals$tax)
# idx <- sort(rhizo_pvals$tax, index.return=T)$ix
# rhizo_pvals <- rhizo_pvals[idx, ]

# root_pvals$tax <- as.character(root_pvals$tax)
# idx <- match(taxonomy[, 1], root_pvals$otu)
# root_pvals <- root_pvals[idx, ]

rhizo_pvals$otu <- factor(rhizo_pvals$otu, levels=rhizo_pvals$otu)

# relative abundances

ra <- apply(otu_table_norm_rhizo_gifu, 1, mean)
rhizo_pvals$ra <- ra[match(rhizo_pvals$otu, names(ra))]

# generate vector of colors

colors <- read.table("taxon_colors.txt", sep="\t", header=T)
taxon <- sort(unique(rhizo_pvals$tax))
colors <- data.frame(taxon=taxon, colors=colors$color[match(taxon, colors$taxon)])
colors$colors <- as.character(colors$colors)
colors$colors[is.na(colors$colors)] <- "grey"
enriched_taxa <- unique(taxonomy[taxonomy[, 1] %in% c(rhizo_otus), 5])
colors$colors[!colors$taxon %in% enriched_taxa] <- "grey"

# multiple testing correction thresholds

BF <- -log(0.05 / dim(rhizo_pvals)[1])
FDR <- min(rhizo_pvals$neglogp[rhizo_pvals$enrichment==TRUE])

# plot

p1 <- ggplot(rhizo_pvals, aes(x=otu, y=neglogp, color=tax, size=ra, shape=enrichment)) +
             geom_point(alpha=.8) +
             geom_hline(yintercept=FDR, linetype=2, color="lightgrey") +
             scale_color_manual(values=colors$colors) +
             scale_shape_manual(values=c(21, 19)) +
             scale_size(breaks=c(0, 0.04, 0.08, 0.12, 0.16)) +
             ylim(c(0, 95)) +
             labs(x="OTU", y="-log10(P)") +
             main_theme +
             theme(axis.ticks.x=element_blank(),
                   axis.text.x=element_blank(),
                   legend.position="none")
             
ggsave(paste(figures.dir, "manhattan_gifu_rhizo.pdf", sep=""), p1, width=10, height=3, useDingbats=F)

