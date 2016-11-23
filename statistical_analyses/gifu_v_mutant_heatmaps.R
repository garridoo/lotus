
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de
#

library(edgeR)

### root compartment - gifu v. mutant

# subset samples

idx <- design$soil %in% c("CAS10") &
       design$compartment %in% c("root") & 
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

gifu_otus_root <- rownames(otu_table_subset)[de_gifu_mutant==1]

lrt_mutant_gifu <- glmLRT(fit, contrast=c(-1, 1))
de_mutant_gifu <- decideTestsDGE(lrt_mutant_gifu, adjust.method=p.adj.method, p.value=alpha)

mutant_otus_root <- rownames(otu_table_subset)[de_mutant_gifu==1]

### rhizosphere compartment - gifu v. mutant

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

lrt_mutant_gifu <- glmLRT(fit, contrast=c(-1, 1))
de_mutant_gifu <- decideTestsDGE(lrt_mutant_gifu, adjust.method=p.adj.method, p.value=alpha)

mutant_otus_rhizo <- rownames(otu_table_subset)[de_mutant_gifu==1]

### root AND rhizosphere compartments - gifu v. mutant

# subset samples

idx <- design$compartment %in% c("root", "rhizosphere") & 
       design$genotype %in% c("gifu", "hit1_1", "nfr5_2", "nfr5_3", "nin2")
design_subset <- design[idx, ]
otu_table_subset <- otu_table[, idx]

design_subset$genotype <- as.character(design_subset$genotype)
idx <- design_subset$genotype %in% c("hit1_1", "nfr5_2", "nfr5_3", "nin2")
design_subset$genotype[idx] <- "mutant"

### Generalized Linear Model (GLM)

# create DGE list

groups <- design_subset$genotype
#groups <- droplevels(groups)

d <- DGEList(counts=otu_table_subset, group=groups)
d <- calcNormFactors(d)

# fit the GLM

design.mat <- model.matrix(~ 0 + d$samples$group)
d2 <- estimateGLMCommonDisp(d, design.mat)
d2 <- estimateGLMTagwiseDisp(d2, design.mat)

fit <- glmFit(d2, design.mat)

lrt_gifu_mutant <- glmLRT(fit, contrast=c(1, -1))
de_gifu_mutant <- decideTestsDGE(lrt_gifu_mutant, adjust.method=p.adj.method, p.value=alpha)

gifu_otus_rr <- rownames(otu_table_subset)[de_gifu_mutant==1]

lrt_mutant_gifu <- glmLRT(fit, contrast=c(-1, 1))
de_mutant_gifu <- decideTestsDGE(lrt_mutant_gifu, adjust.method=p.adj.method, p.value=alpha)

mutant_otus_rr <- rownames(otu_table_subset)[de_mutant_gifu==1]




dif_otus <- gifu_otus_root

idx <- design$soil %in% c("CAS10") &
       design$compartment %in% c("root") & 
       design$genotype %in% c("gifu", "hit1_1", "nfr5_2", "nfr5_3", "nin2")
design_subset <- design[idx, ]
otu_table_subset <- otu_table[, idx]

otu_table_norm <- apply(otu_table_subset, 2, function(x) x / sum(x))

df <- melt(otu_table_norm)
colnames(df) <- c("otu", "sample", "value")

df <- df[df$otu %in% dif_otus, ]

df$sample <- factor(df$sample, levels=design_subset$SampleID)
samples <- levels(df$sample)
l <- c(samples[grepl("gifu", samples)],
       samples[grepl("nf2", samples)],
       samples[grepl("nf2", samples)],
       samples[grepl("nf3", samples)],
       samples[grepl("nin", samples)],
       samples[grepl("hit", samples)])
df$sample <- factor(df$sample, levels=l)

means <- apply(otu_table_norm[, grepl("gifu", colnames(otu_table_norm))], 1, mean)
l <- rev(names(sort(means, decreasing=T)))
df$otu <- factor(df$otu, levels=l)

colors <- c("#030303", "#454409", "#7b7807", "#b0a609",
            "#ddd300", "#ffff04", "#ffcd00", "#ff9c0a",
            "#fd6e10", "#fe2f1d")
b <- c(0, .01, .05, .1, .5, 1, 2, 4, 6, 8, 100)
df$perc_bins <- cut(df$value * 100, breaks=b, right=F)

p1 <- ggplot(df, aes(x=sample, y=otu, fill=perc_bins)) +
      geom_tile() +
      scale_fill_manual(values=colors) +
      labs(x="", y="", title="") + 
      main_theme +
      theme(legend.title=element_blank(),
            legend.position="left",
            axis.text.x=element_text(angle=45, hjust=1, size=10),
            axis.ticks=element_blank(),
            axis.line=element_blank(),
            plot.margin=unit(c(0, 0, 0, 0), "mm"))

ggsave(file=paste(figures.dir, "heatmap_gifu_root.pdf", sep=""), p1, height=13, width=15)



dif_otus <- mutant_otus_root

idx <- design$soil %in% c("CAS10") &
       design$compartment %in% c("root") & 
       design$genotype %in% c("gifu", "hit1_1", "nfr5_2", "nfr5_3", "nin2")
design_subset <- design[idx, ]
otu_table_subset <- otu_table[, idx]

otu_table_norm <- apply(otu_table_subset, 2, function(x) x / sum(x))

df <- melt(otu_table_norm)
colnames(df) <- c("otu", "sample", "value")

df <- df[df$otu %in% dif_otus, ]

df$sample <- factor(df$sample, levels=design_subset$SampleID)
samples <- levels(df$sample)
l <- c(samples[grepl("gifu", samples)],
       samples[grepl("nf2", samples)],
       samples[grepl("nf2", samples)],
       samples[grepl("nf3", samples)],
       samples[grepl("nin", samples)],
       samples[grepl("hit", samples)])
df$sample <- factor(df$sample, levels=l)

means <- apply(otu_table_norm[, !grepl("gifu", colnames(otu_table_norm))], 1, mean)
l <- rev(names(sort(means, decreasing=F)))
df$otu <- factor(df$otu, levels=l)

colors <- c("#030303", "#454409", "#7b7807", "#b0a609",
            "#ddd300", "#ffff04", "#ffcd00", "#ff9c0a",
            "#fd6e10", "#fe2f1d")
b <- c(0, .01, .05, .1, .5, 1, 2, 4, 6, 8, 100)
df$perc_bins <- cut(df$value * 100, breaks=b, right=F)

p1 <- ggplot(df, aes(x=sample, y=otu, fill=perc_bins)) +
      geom_tile() +
      scale_fill_manual(values=colors) +
      labs(x="", y="", title="") + 
      main_theme +
      theme(legend.title=element_blank(),
            legend.position="left",
            axis.text.x=element_text(angle=45, hjust=1, size=10),
            axis.ticks=element_blank(),
            axis.line=element_blank(),
            plot.margin=unit(c(0, 0, 0, 0), "mm"))

ggsave(file=paste(figures.dir, "heatmap_mutant_root.pdf", sep=""), p1, height=13, width=15)




dif_otus <- gifu_otus_rhizo

idx <- design$soil %in% c("CAS10") &
       design$compartment %in% c("rhizosphere") & 
       design$genotype %in% c("gifu", "hit1_1", "nfr5_2", "nfr5_3", "nin2")
design_subset <- design[idx, ]
otu_table_subset <- otu_table[, idx]

otu_table_norm <- apply(otu_table_subset, 2, function(x) x / sum(x))

df <- melt(otu_table_norm)
colnames(df) <- c("otu", "sample", "value")

df <- df[df$otu %in% dif_otus, ]

df$sample <- factor(df$sample, levels=design_subset$SampleID)
samples <- levels(df$sample)
l <- c(samples[grepl("gifu", samples)],
       samples[grepl("nf2", samples)],
       samples[grepl("nf2", samples)],
       samples[grepl("nf3", samples)],
       samples[grepl("nin", samples)],
       samples[grepl("hit", samples)])
df$sample <- factor(df$sample, levels=l)

means <- apply(otu_table_norm[, grepl("gifu", colnames(otu_table_norm))], 1, mean)
l <- rev(names(sort(means, decreasing=T)))
df$otu <- factor(df$otu, levels=l)

colors <- c("#030303", "#454409", "#7b7807", "#b0a609",
            "#ddd300", "#ffff04", "#ffcd00", "#ff9c0a",
            "#fd6e10", "#fe2f1d")
b <- c(0, .1, 1, 1.5, 2, 4, 8, 5, 6, 7, 100)
b <- c(0, .01, .1, 1, 5, 10, 15, 20, 25, 30, 100)
b <- c(0, .01, .05, .1, .5, 1, 2, 4, 6, 8, 100)
df$perc_bins <- cut(df$value * 100, breaks=b, right=F)

p1 <- ggplot(df, aes(x=sample, y=otu, fill=perc_bins)) +
      geom_tile() +
      scale_fill_manual(values=colors) +
      labs(x="", y="", title="") + 
      main_theme +
      theme(legend.title=element_blank(),
            legend.position="left",
            axis.text.x=element_text(angle=45, hjust=1, size=10),
            axis.ticks=element_blank(),
            axis.line=element_blank(),
            plot.margin=unit(c(0, 0, 0, 0), "mm"))

ggsave(file=paste(figures.dir, "heatmap_gifu_rhizo.pdf", sep=""), p1, height=13, width=15)



dif_otus <- mutant_otus_rhizo

idx <- design$soil %in% c("CAS10") &
       design$compartment %in% c("rhizosphere") & 
       design$genotype %in% c("gifu", "hit1_1", "nfr5_2", "nfr5_3", "nin2")
design_subset <- design[idx, ]
otu_table_subset <- otu_table[, idx]

otu_table_norm <- apply(otu_table_subset, 2, function(x) x / sum(x))

df <- melt(otu_table_norm)
colnames(df) <- c("otu", "sample", "value")

df <- df[df$otu %in% dif_otus, ]

df$sample <- factor(df$sample, levels=design_subset$SampleID)
samples <- levels(df$sample)
l <- c(samples[grepl("gifu", samples)],
       samples[grepl("nf2", samples)],
       samples[grepl("nf2", samples)],
       samples[grepl("nf3", samples)],
       samples[grepl("nin", samples)],
       samples[grepl("hit", samples)])
df$sample <- factor(df$sample, levels=l)

means <- apply(otu_table_norm[, !grepl("gifu", colnames(otu_table_norm))], 1, mean)
l <- rev(names(sort(means, decreasing=F)))
df$otu <- factor(df$otu, levels=l)

colors <- c("#030303", "#454409", "#7b7807", "#b0a609",
            "#ddd300", "#ffff04", "#ffcd00", "#ff9c0a",
            "#fd6e10", "#fe2f1d")
b <- c(0, .1, 1, 1.5, 2, 4, 8, 5, 6, 7, 100)
b <- c(0, .01, .1, 1, 5, 10, 15, 20, 25, 30, 100)
b <- c(0, .01, .05, .1, .5, 1, 2, 4, 6, 8, 100)
df$perc_bins <- cut(df$value * 100, breaks=b, right=F)

p1 <- ggplot(df, aes(x=sample, y=otu, fill=perc_bins)) +
      geom_tile() +
      scale_fill_manual(values=colors) +
      labs(x="", y="", title="") + 
      main_theme +
      theme(legend.title=element_blank(),
            legend.position="left",
            axis.text.x=element_text(angle=45, hjust=1, size=10),
            axis.ticks=element_blank(),
            axis.line=element_blank(),
            plot.margin=unit(c(0, 0, 0, 0), "mm"))

ggsave(file=paste(figures.dir, "heatmap_mutant_rhizo.pdf", sep=""), p1, height=13, width=15)

