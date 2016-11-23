
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de
#

library(edgeR)

# subset samples

idx <- design$compartment %in% c("root", "soil", "rhizosphere") & 
       design$genotype %in% c("hit1_1", "nfr5_2", "nfr5_3", "nin2", "soil")
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
rhizo_otus_mutant <- rhizo_otus

lrt_root_rhizo <- glmLRT(fit, contrast=c(-1, 1, 0))
lrt_root_soil <- glmLRT(fit, contrast=c(0, 1, -1))

de_root_rhizo <- decideTestsDGE(lrt_root_rhizo, adjust.method=p.adj.method, p.value=alpha)
de_root_soil <- decideTestsDGE(lrt_root_soil, adjust.method=p.adj.method, p.value=alpha)

root_otus <- rownames(otu_table_subset)[de_root_rhizo==1 & de_root_soil==1]
root_otus_mutant <- root_otus

lrt_soil_rhizo <- glmLRT(fit, contrast=c(-1, 0, 1))
lrt_soil_root <- glmLRT(fit, contrast=c(0, -1, 1))

de_soil_rhizo <- decideTestsDGE(lrt_soil_rhizo, adjust.method=p.adj.method, p.value=alpha)
de_soil_root <- decideTestsDGE(lrt_soil_root, adjust.method=p.adj.method, p.value=alpha)

soil_otus <- rownames(otu_table_subset)[de_soil_rhizo==1 & de_soil_root==1]

enriched_otus <- data.frame(root_enriched=rownames(otu_table) %in% root_otus,
                            rhizo_enriched=rownames(otu_table) %in% rhizo_otus,
                            soil_enriched=rownames(otu_table) %in% soil_otus)
row.names(enriched_otus) <- rownames(otu_table)

write.table(enriched_otus, paste(results.dir, "mutant_root_rhizo_soil_OTUs.txt", sep=""),
            quote=F, sep="\t", col.names=T, row.names=T)

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

pdf(file=paste(figures.dir, "mutant_root_rhizo_soil_enrichment.pdf", sep=""))

tern_e(df, prop_size=T, col=colors, grid_color="grey",
       labels_color="transparent", pch=19, main="mutant")

dev.off()

### boxplots of aggregated relative abundances

idx <- design_subset$compartment=="rhizosphere"

rhizo_rhizo <- colSums(otu_table_norm[rownames(otu_table_norm) %in% rhizo_otus, idx])
root_rhizo <- colSums(otu_table_norm[rownames(otu_table_norm) %in% root_otus, idx])
soil_rhizo <- colSums(otu_table_norm[rownames(otu_table_norm) %in% soil_otus, idx])

df_rhizo <- rbind(data.frame(otus="rhizo", compartment="rhizo", ra=rhizo_rhizo),
                  data.frame(otus="root", compartment="rhizo", ra=root_rhizo),
                  data.frame(otus="soil", compartment="rhizo", ra=soil_rhizo))

idx <- design_subset$compartment=="root"

rhizo_root <- colSums(otu_table_norm[rownames(otu_table_norm) %in% rhizo_otus, idx])
root_root <- colSums(otu_table_norm[rownames(otu_table_norm) %in% root_otus, idx])
soil_root <- colSums(otu_table_norm[rownames(otu_table_norm) %in% soil_otus, idx])

df_root <- rbind(data.frame(otus="rhizo", compartment="root", ra=rhizo_root),
                 data.frame(otus="root", compartment="root", ra=root_root),
                 data.frame(otus="soil", compartment="root", ra=soil_root))

idx <- design_subset$compartment=="soil"

rhizo_soil <- colSums(otu_table_norm[rownames(otu_table_norm) %in% rhizo_otus, idx])
root_soil <- colSums(otu_table_norm[rownames(otu_table_norm) %in% root_otus, idx])
soil_soil <- colSums(otu_table_norm[rownames(otu_table_norm) %in% soil_otus, idx])

df_soil <- rbind(data.frame(otus="rhizo", compartment="soil", ra=rhizo_soil),
                 data.frame(otus="root", compartment="soil", ra=root_soil),
                 data.frame(otus="soil", compartment="soil", ra=soil_soil))

df <- rbind(df_rhizo, df_root, df_soil)

df$compartment <- factor(df$compartment, levels=c("soil", "rhizo", "root"))
df$otus <- factor(df$otus, levels=c("soil", "rhizo", "root"))

p1 <- ggplot(df, aes(x=compartment, y=ra, color=otus)) +
             geom_boxplot(alpha=1, outlier.size=0, size=0.6, width=0.8,
                          position=position_dodge(width=0.7)) +
             scale_color_manual(values=c(c_dark_brown, c_dark_red, c_very_dark_green)) +
             scale_y_continuous(labels=percent, limits=c(0, 1)) +
             labs(x="", y="relative abundance") +
             main_theme +
             theme(legend.position="none")

ggsave(paste(figures.dir, "boxplots_mutant_root_rhizo_soil.pdf", sep=""), p1, width=5, height=3)

