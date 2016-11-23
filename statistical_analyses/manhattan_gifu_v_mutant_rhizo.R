
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de
#

library(edgeR)

### root compartment - gifu v. mutant

# subset samples

idx <- design$soil %in% c("CAS10") &
       design$compartment %in% c("rhizosphere") & 
       design$genotype %in% c("gifu", "hit1_1", "nfr5_2", "nfr5_3", "nin2")
design_subset <- design[idx, ]
otu_table_subset <- otu_table[, idx]

design_subset$genotype <- as.character(design_subset$genotype)
idx <- design_subset$genotype %in% c("hit1_1", "nfr5_2", "nfr5_3", "nin2")
design_subset$genotype[idx] <- "mutant"

### Generalized Linear Model (GLM)

# create DGE list

groups <- design_subset$genotype

d <- DGEList(counts=otu_table_subset, group=groups)
d <- calcNormFactors(d)

# fit the GLM

design.mat <- model.matrix(~ 0 + d$samples$group)
d2 <- estimateGLMCommonDisp(d, design.mat)
d2 <- estimateGLMTagwiseDisp(d2, design.mat)

fit <- glmFit(d2, design.mat)

lrt_gifu_mutant <- glmLRT(fit, contrast=c(1, -1))
de_gifu_mutant <- decideTestsDGE(lrt_gifu_mutant, adjust.method=p.adj.method, p.value=alpha)

gifu_otus_rhizo <- rownames(otu_table_subset)[de_gifu_mutant==1]

gifu_mutant_pvals <- lrt_gifu_mutant$table

lrt_mutant_gifu <- glmLRT(fit, contrast=c(-1, 1))
de_mutant_gifu <- decideTestsDGE(lrt_mutant_gifu, adjust.method=p.adj.method, p.value=alpha)

mutant_otus_rhizo <- rownames(otu_table_subset)[de_mutant_gifu==1]

mutant_gifu_pvals <- lrt_mutant_gifu$table

# write table of enriched OTUs

enriched_otus_rhizo <- data.frame(gifu_enriched=rownames(otu_table) %in% gifu_otus_rhizo,
                            mutant_enriched=rownames(otu_table) %in% mutant_otus_rhizo)
row.names(enriched_otus_rhizo) <- rownames(otu_table)

write.table(enriched_otus_rhizo, paste(results.dir, "gifu_v_mutant_rhizo_OTUs.txt", sep=""),
            quote=F, sep="\t", col.names=T, row.names=T)

### get mean R.A. in gifu rhizo samples

otu_table_norm_all <- apply(otu_table, 2, function(x) x / sum(x))
 
idx <- design$compartment %in% c("rhizosphere") &
       design$soil %in% c("CAS10") &
       design$genotype %in% c("gifu")
rhizo_gifu_samples <- design$SampleID[idx]
otu_table_norm_rhizo_gifu <- otu_table_norm_all[, colnames(otu_table_norm_all) %in% rhizo_gifu_samples]

### get mean R.A. in mutant rhizo samples

idx <- design$compartment %in% c("rhizosphere") &
       design$soil %in% c("CAS10") &
       design$genotype %in% c("hit1_1", "nfr5_2", "nfr5_3", "nin2")
rhizo_mutant_samples <- design$SampleID[idx]
otu_table_norm_rhizo_mutant <- otu_table_norm_all[, colnames(otu_table_norm_all) %in% rhizo_mutant_samples]

### gifu rhizosphere-enriched OTUs

### generate data frame for plotting

gifu_mutant_pvals$otu <- rownames(gifu_mutant_pvals)

# P values

gifu_mutant_pvals$neglogp <- -log(gifu_mutant_pvals$PValue)

# enrichment status

gifu_mutant_pvals$enrichment <- gifu_mutant_pvals$otu %in% gifu_otus_rhizo

# log fold change

idx <- gifu_mutant_pvals$logFC < 0
gifu_mutant_pvals$neglogp[idx] <- 0

# order OTUs according to taxonomy

taxonomy <- taxonomy[order(taxonomy[, 3], taxonomy[, 4], taxonomy[, 5]), ]
idx <- taxonomy[, 1] %in% gifu_mutant_pvals$otu
taxonomy <- taxonomy[idx, ]

idx <- match(taxonomy[, 1], gifu_mutant_pvals$otu)
gifu_mutant_pvals <- gifu_mutant_pvals[idx, ]

gifu_mutant_pvals$tax <- taxonomy[, 5]

# gifu_mutant_pvals$tax <- as.character(gifu_mutant_pvals$tax)
# idx <- sort(gifu_mutant_pvals$tax, index.return=T)$ix
# gifu_mutant_pvals <- gifu_mutant_pvals[idx, ]

gifu_mutant_pvals$otu <- factor(gifu_mutant_pvals$otu, levels=gifu_mutant_pvals$otu)

# relative abundances

ra <- apply(otu_table_norm_rhizo_gifu, 1, mean)
gifu_mutant_pvals$ra <- ra[match(gifu_mutant_pvals$otu, names(ra))]

# generate vector of colors

colors <- read.table("taxon_colors.txt", sep="\t", header=T)
taxon <- sort(unique(gifu_mutant_pvals$tax))
colors <- data.frame(taxon=taxon, colors=colors$color[match(taxon, colors$taxon)])
colors$colors <- as.character(colors$colors)
colors$colors[is.na(colors$colors)] <- "grey"
enriched_taxa <- unique(taxonomy[taxonomy[, 1] %in% c(gifu_otus_rhizo), 5])
colors$colors[!colors$taxon %in% enriched_taxa] <- "grey"

# multiple testing correction thresholds

BF <- -log(0.05 / dim(gifu_mutant_pvals)[1])
FDR <- min(gifu_mutant_pvals$neglogp[gifu_mutant_pvals$enrichment==TRUE])

# plot

p1 <- ggplot(gifu_mutant_pvals, aes(x=otu, y=neglogp, color=tax, size=ra, shape=enrichment)) +
             geom_point(alpha=.8) +
             geom_hline(yintercept=FDR, linetype=2, color="lightgrey") +
             scale_color_manual(values=colors$colors) +
             scale_shape_manual(values=c(21, 19)) +
             scale_size(breaks=c(0, 0.04, 0.08, 0.12, 0.16)) +
             ylim(c(0, 70)) +
             labs(x="OTU", y="-log10(P)") +
             main_theme +
             theme(axis.ticks.x=element_blank(),
                   axis.text.x=element_blank(),
                   legend.position="none")

ggsave(paste(figures.dir, "manhattan_gifu_v_mutant_rhizo.pdf", sep=""), p1, width=10, height=3, useDingbats=F)

### mutant rhizosphere-enriched OTUs

### generate data frame for plotting

mutant_gifu_pvals$otu <- rownames(mutant_gifu_pvals)

# P values

mutant_gifu_pvals$neglogp <- -log(mutant_gifu_pvals$PValue)

# enrichment status

mutant_gifu_pvals$enrichment <- mutant_gifu_pvals$otu %in% mutant_otus_rhizo

# log fold change

idx <- mutant_gifu_pvals$logFC < 0
mutant_gifu_pvals$neglogp[idx] <- 0

# order OTUs according to taxonomy

taxonomy <- taxonomy[order(taxonomy[, 3], taxonomy[, 4], taxonomy[, 5]), ]
idx <- taxonomy[, 1] %in% mutant_gifu_pvals$otu
taxonomy <- taxonomy[idx, ]

idx <- match(taxonomy[, 1], mutant_gifu_pvals$otu)
mutant_gifu_pvals <- mutant_gifu_pvals[idx, ]

mutant_gifu_pvals$tax <- taxonomy[, 5]

# mutant_gifu_pvals$tax <- as.character(mutant_gifu_pvals$tax)
# idx <- sort(mutant_gifu_pvals$tax, index.return=T)$ix
# mutant_gifu_pvals <- mutant_gifu_pvals[idx, ]

mutant_gifu_pvals$otu <- factor(mutant_gifu_pvals$otu, levels=mutant_gifu_pvals$otu)

# relative abundances

ra <- apply(otu_table_norm_rhizo_mutant, 1, mean)
mutant_gifu_pvals$ra <- ra[match(mutant_gifu_pvals$otu, names(ra))]

# generate vector of colors

colors <- read.table("taxon_colors.txt", sep="\t", header=T)
taxon <- sort(unique(mutant_gifu_pvals$tax))
colors <- data.frame(taxon=taxon, colors=colors$color[match(taxon, colors$taxon)])
colors$colors <- as.character(colors$colors)
colors$colors[is.na(colors$colors)] <- "grey"
enriched_taxa <- unique(taxonomy[taxonomy[, 1] %in% c(mutant_otus_rhizo), 5])
colors$colors[!colors$taxon %in% enriched_taxa] <- "grey"

# multiple testing correction thresholds

BF <- -log(0.05 / dim(mutant_gifu_pvals)[1])
FDR <- min(mutant_gifu_pvals$neglogp[mutant_gifu_pvals$enrichment==TRUE])

# plot

p1 <- ggplot(mutant_gifu_pvals, aes(x=otu, y=neglogp, color=tax, size=ra, shape=enrichment)) +
             geom_point(alpha=.8) +
             geom_hline(yintercept=FDR, linetype=2, color="lightgrey") +
             scale_color_manual(values=colors$colors) +
             scale_shape_manual(values=c(21, 19)) +
             scale_size(breaks=c(0, 0.04, 0.08, 0.12, 0.16)) +
             ylim(c(0, 70)) +
             labs(x="OTU", y="-log10(P)") +
             main_theme +
             theme(axis.ticks.x=element_blank(),
                   axis.text.x=element_blank(),
                   legend.position="none")
             
ggsave(paste(figures.dir, "manhattan_mutant_v_gifu_rhizo.pdf", sep=""), p1, width=10, height=3, useDingbats=F)

