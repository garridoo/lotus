
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de
#

library(edgeR)

### nin2 mutant: root/rhizosphere/soil enriched OTUs

# subset samples

idx <- design$compartment %in% c("root", "soil", "rhizosphere") & 
       design$genotype %in% c("nin2", "soil") &
       design$soil=="CAS10"
design_subset <- design[idx, ]
otu_table_subset <- otu_table[, idx]

### Generalized Linear Model (GLM)

# create DGE list

groups <- design_subset$compartment
groups <- droplevels(groups)

d <- DGEList(counts=otu_table_subset, group=groups)
d <- calcNormFactors(d)

# fit the GLM

design.mat <- model.matrix(~ 0 + d$samples$group)
d2 <- estimateGLMCommonDisp(d, design.mat)
d2 <- estimateGLMTagwiseDisp(d2, design.mat)

fit <- glmFit(d2, design.mat)

lrt_rhizo_root <- glmLRT(fit, contrast=c(1, -1, 0))
lrt_rhizo_soil <- glmLRT(fit, contrast=c(1, 0, -1))

de_rhizo_root <- decideTestsDGE(lrt_rhizo_root, adjust.method=p.adj.method, p.value=alpha)
de_rhizo_soil <- decideTestsDGE(lrt_rhizo_soil, adjust.method=p.adj.method, p.value=alpha)

rhizo_otus <- rownames(otu_table_subset)[de_rhizo_root==1 & de_rhizo_soil==1]

lrt_root_rhizo <- glmLRT(fit, contrast=c(-1, 1, 0))
lrt_root_soil <- glmLRT(fit, contrast=c(0, 1, -1))

de_root_rhizo <- decideTestsDGE(lrt_root_rhizo, adjust.method=p.adj.method, p.value=alpha)
de_root_soil <- decideTestsDGE(lrt_root_soil, adjust.method=p.adj.method, p.value=alpha)

root_otus <- rownames(otu_table_subset)[de_root_rhizo==1 & de_root_soil==1]

lrt_soil_rhizo <- glmLRT(fit, contrast=c(-1, 0, 1))
lrt_soil_root <- glmLRT(fit, contrast=c(0, -1, 1))

de_soil_rhizo <- decideTestsDGE(lrt_soil_rhizo, adjust.method=p.adj.method, p.value=alpha)
de_soil_root <- decideTestsDGE(lrt_soil_root, adjust.method=p.adj.method, p.value=alpha)

soil_otus <- rownames(otu_table_subset)[de_soil_rhizo==1 & de_soil_root==1]

enriched_otus <- data.frame(root_enriched=rownames(otu_table) %in% root_otus,
                            rhizo_enriched=rownames(otu_table) %in% rhizo_otus,
                            soil_enriched=rownames(otu_table) %in% soil_otus)
row.names(enriched_otus) <- rownames(otu_table)

### ternary plots

# normalize subsetted OTU table and apply log transform

otu_table_norm <- apply(otu_table_subset, 2, function(x) x / sum(x))
otu_table_norm_log <- log2(otu_table_norm + 1)

# create vectors of mean reltive abundances

idx <- design_subset$compartment=="rhizosphere"
rhizo_means <- apply(otu_table_norm[, idx], 1, mean)

idx <- design_subset$compartment=="root"
root_means <- apply(otu_table_norm[, idx], 1, mean)

idx <- design_subset$compartment=="soil"
soil_means <- apply(otu_table_norm[, idx], 1, mean)

# create matrix of average r.a. per group

df <- data.frame(root=root_means, rhizosphere=rhizo_means, soil=soil_means)
df <- df[rowSums(df)!=0, ]
df <- log2(df * scale + 1)

# sort the rows by decreasing abundance (looks better)

idx <- sort(rowSums(df), decreasing=F, index.return=T)$ix
df <- df[idx, ]

# create vector of colors according to enrichment

colors <- rep(c_grey, dim(df)[1])
colors[rownames(df) %in% root_otus] <- c_very_dark_green
colors[rownames(df) %in% soil_otus] <- c_dark_brown
colors[rownames(df) %in% rhizo_otus] <- c_dark_red

idx <- sort(colors==c_grey, decreasing=T, index.return=T)$ix
df <- df[idx, ]
colors <- colors[idx]

# plot colored by enrichment

pdf(file=paste(figures.dir, "nin2_root_rhizo_soil_enrichment_cas10.pdf", sep=""))

tern_e(df, prop_size=T, col=colors, grid_color="grey",
       labels_color="transparent", pch=19, main="nin2")

dev.off()

