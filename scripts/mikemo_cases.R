library(ggplot2)
library(ggpubr)
library(phyloseq)
library(DESeq2)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(microViz)
library(dplyr)
#import data
ovs<-read.csv("ov-pos.csv",check.names = F)[,1:26]
colnames(ovs)
ovstax <- read.csv("ov-pos.csv",check.names = F)[,c(1,28:33)]
colnames(ovstax)
mgs <- read.csv("mg-pos.csv",check.names = F)[,1:24]
colnames(mgs)
mgstax <- read.csv("mg-pos.csv",check.names = F)[,c(1,26:31)]
colnames(mgstax)
#merge the two data sets by species column
mergeddata<-merge(x=ovs,y=mgs, by="genus",all=T)
mergeddata
mergeddata[is.na(mergeddata)]<-0
mergeddata
rownames(mergeddata)<-mergeddata$genus
mergeddata<-mergeddata[,-1]

#taxa table
mergetax <- rbind(ovstax, mgstax)
mergetax
mergetax <- distinct(mergetax,genus,.keep_all = T)
rownames(mergetax)<-mergetax$genus
mergetax<-mergetax[,-1]

#load metadata
metadata <- read.csv("casemetadata.csv",row.names = 1)
metadata

#check if metadata and row names are present and in same order
all(rownames(metadata)%in%colnames(mergeddata))
all(rownames(metadata)==colnames(mergeddata))

all(rownames(mergetax)%in%rownames(mergeddata))
all(rownames(mergetax)==rownames(mergeddata))

# create phyloseq object

otu1<-otu_table(mergeddata, taxa_are_rows = T)
sample<-phyloseq::sample_data(metadata)
tax <- tax_table(mergetax)
rownames(tax) <- rownames(mergetax)
colnames(tax) <- colnames(mergetax)
ps_object<-phyloseq(otu1, sample, tax)
ps_object

#estimate alpha diversity
alpha_meas = c("Shannon", "Simpson")
#?plot_richness
#MMB
(p <- plot_richness(ps_object, "tissue",measures=alpha_meas))
p$layers <- p$layers[-1]
#MMB plot
p + geom_boxplot(data=p$data, aes(x=tissue, y=value, fill = tissue), 
  alpha=0.2, show.legend = F, outlier.shape = NA) + geom_jitter(width = 0.2) + 
  theme_bw(base_size = 30) + 
  stat_compare_means(method="kruskal.test",label.y.npc = "bottom" ,label.x.npc = "center") +
  labs(x="Mosquito tissue")

# #species
# (p <- plot_richness(ps_object, "mspecies",measures=alpha_meas))
# p$layers <- p$layers[-1]
# 
# p + geom_boxplot(data=p$data, aes(x=mspecies, y=value, fill = mspecies), 
#                  alpha=0.2, show.legend = F, outlier.shape = NA) + geom_jitter(width = 0.2) + 
#   theme(axis.text.x = element_text(angle = 0, hjust=0.5)) + theme_pubclean() + 
#   stat_compare_means(method="kruskal.test", label.y.npc = "bottom", 
#                      label.x.npc = "left") + 
#   labs(title="Alpha diversity plots Shannon and Simpson",x="Mosquito Species")
# #site
# (p <- plot_richness(ps_object, "site",measures=alpha_meas))
# p$layers <- p$layers[-1]
# 
# p + geom_boxplot(data=p$data, aes(x=site, y=value, fill = site), 
#                  alpha=0.2, show.legend = F, outlier.shape = NA) + geom_jitter(width = 0.2) + 
#   theme(axis.text.x = element_text(angle = 0, hjust=0.5)) + theme_pubclean() + 
#   stat_compare_means(method="kruskal.test", label.y.npc = "bottom", 
#                      label.x.npc = "left") + 
#   labs(title="Alpha diversity plots Shannon and Simpson",x="Site")
# 
# #loc
# (p <- plot_richness(ps_object, "loc",measures=alpha_meas))
# p$layers <- p$layers[-1]
# 
# p + geom_boxplot(data=p$data, aes(x=loc, y=value, fill = loc), 
#                  alpha=0.2, show.legend = F, outlier.shape = NA) + geom_jitter(width = 0.2) + 
#   theme(axis.text.x = element_text(angle = 0, hjust=0.5)) + theme_pubclean() + 
#   stat_compare_means(method="kruskal.test", label.y.npc = "bottom", 
#                      label.x.npc = "left") + 
#   labs(title="Alpha diversity plots Shannon and Simpson",x="Location in-door or out-door")
# 
# #House
# (p <- plot_richness(ps_object, "house",measures=alpha_meas))
# p$layers <- p$layers[-1]
# 
# p + geom_boxplot(data=p$data, aes(x=house, y=value, fill = house), 
#                  alpha=0.2, show.legend = F, outlier.shape = NA) + geom_jitter(width = 0.2) + 
#   theme(axis.text.x = element_text(angle = 0, hjust=0.5)) + theme_pubclean() + 
#   stat_compare_means(method="kruskal.test", label.y.npc = "bottom", 
#                      label.x.npc = "left") + 
#   labs(title="Alpha diversity plots Shannon and Simpson",x="Household ID")

#t.test, wilcox.test, anova, kruskal.test
# Calculate distances
DistBC = phyloseq::distance(ps_object, method = "bray") #jaccard, unifrac, wunifrac
ordBC = ordinate(ps_object, method = "PCoA", distance = DistBC)
plot_scree(ordBC, "Scree Plot: Bray-Curtis MDS")
plot_ordination(ps_object, ordBC, color = "tissue") +
  geom_point(aes(color=tissue), size = 4) + 
  theme_bw(base_size = 20) + labs(color="Tissue")

# plot_ordination(ps_object, ordBC, color = "tissue", shape = "site") +
#   geom_point(aes(color=tissue), size = 4) +
#   ggtitle("PCoA: Bray-Curtis") + labs(color="MMB Positivity",shape="Site")
# 
# plot_ordination(ps_object, ordBC, color = "tissue", shape = "loc") +
#   geom_point(aes(color=tissue), size = 4) +
#   ggtitle("PCoA: Bray-Curtis") + labs(color="MMB Positivity",shape="Location")
# 
# 
# 
# DistJD = phyloseq::distance(ps_object, method = "manhattan") #jaccard, unifrac, wunifrac
# ordJD = ordinate(ps_object, method = "PCoA", distance = DistJD)
# plot_scree(ordJD, "Scree Plot: Jaccard MDS")
# plot_ordination(ps_object, ordJD, color = "tissue", shape = "tissue") +
#   geom_point(aes(color=tissue), size = 4) +
#   ggtitle("PCoA: Jaccard Distance") + labs(color="MMB Positivity",shape="Mosquito Spp")
# 
# plot_ordination(ps_object, ordJD, color = "tissue", shape = "site") +
#   geom_point(aes(color=tissue), size = 4) +
#   ggtitle("PCoA: Jaccard Distance") + labs(color="MMB Positivity",shape="Site")
# 
# plot_ordination(ps_object, ordJD, color = "tissue", shape = "loc") +
#   geom_point(aes(color=tissue), size = 4) +
#   ggtitle("PCoA: Jaccard Distance") + labs(color="MMB Positivity",shape="Location")
# #ordBC$vectors

# df<-as.data.frame(ordBC$vectors[,1:2])
# df$condition<-metadata$condition
# df
# 
# ggplot(df, aes(Axis.1, Axis.2, color=condition))+geom_point(size=4)#PERMANOVA

#############################################################################
diagdds = phyloseq_to_deseq2(ps_object, ~ tissue)
#diagdds = phyloseq_to_deseq2(ps_object, ~ site)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds$sizeFactor



######################################################################################
diagdds = DESeq(diagdds,
                fitType="local",
                test = "LRT",
                reduced = ~1)
resultsNames(diagdds)

ovs_vs_mgs = results(diagdds,#filterFun = ihw,
                        independentFiltering = FALSE,
                        test = "Wald",
                        name = "tissue_ovary_vs_midgut")
summary(ovs_vs_mgs)

ovs_vs_mgs<-ovs_vs_mgs[order(ovs_vs_mgs$padj),]
ovs_vs_mgs
# ovs_vs_mgs_sig<-subset(ovs_vs_mgs, ovs_vs_mgs$padj<0.05)
# ovs_vs_mgs_sig

# 
# rpm<-fpm(diagdds)
# rpm

# #heatmap(rpm[rownames(ovs_vs_mgs_sig),])
# 
# Heatmap(t(scale(t(log2(rpm[rownames(ovs_vs_mgs_sig),]+2)))), column_split = metadata$tissue)

#BiocManager::install("ComplexHeatmap")
#BiocManager::install("EnhancedVolcano")
# plot volcano, MA plots

ovs_vs_mgs.df <- as.data.frame(ovs_vs_mgs)
# EnhancedVolcano(ovs_vs_mgs.df,x="log2FoldChange",y="padj",lab = row.names(ovs_vs_mgs.df),
#                 FCcutoff = 2,pCutoffCol = "padj", pCutoff = 0.05)
ovs_vs_mgs.df$sig <- ifelse(ovs_vs_mgs.df$padj<0.05,"yes","no")
ovs_vs_mgs.df$labs <- ifelse(ovs_vs_mgs.df$sig=="yes",rownames(ovs_vs_mgs.df),NA)
ggmaplot(ovs_vs_mgs.df, fdr = 0.05, fc = 2, size = 2,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(ovs_vs_mgs.df$labs),label.rectangle = TRUE,
         font.label = c(19, "plain", "black"),
         xlab = "Log2 mean abundance",ggtheme = theme_bw(base_size = 19))

# ggplot(ovs_vs_mgs.df, aes(log2FoldChange, -log10(padj), color=sig)) + 
#   geom_point()

#barplot top 20

# top20 <- names(sort(taxa_sums(ps_object), decreasing = T))[1:300]
# pstop20 <- transform_sample_counts(ps_object, function(OTU) OTU/sum(OTU))
# pstop20 <- prune_taxa(top20, pstop20)
# pstop20 <- prune_samples(sample_sums(pstop20) > 10, pstop20)
# plot_bar(pstop20, facet_grid = ~condition, fill="class")
# plot_bar(ps_objectm, facet_grid = ~condition, fill="class")
# plot_heatmap(pstop20, sample.label = "condition",sample.order="condition",
#              taxa.label = "class", taxa.order = "class")
ps.df <- tax_transform(ps_object,"identity", rank = "class") #identity compositional
# ps.df
# ps.df <- tax_filter(ps.df, min_total_abundance = 10)
# ps.df
# ps.dff <- otu_table(ps.df)
# ps.dff <- as.matrix(ps.dff)
# Heatmap(ps.dff,name="Composition",row_km=2,row_km_repeats=100,
#         column_split = metadata$tissue, left_annotation = row_ha) 
# row_ha <- rowAnnotation(count=taxa_sums(ps.df),bars=anno_barplot(taxa_sums(ps.df)),
#                         )
# Heatmap(ps.dff, name="Counts", column_split = metadata$tissue,
#         show_column_dend = FALSE, show_row_dend = FALSE)

(pp <- comp_barplot(ps.df,
                    n_taxa = 41, tax_order = sum, tax_level = "class",
                    bar_outline_colour = "black", merge_other = FALSE,
                    sample_order = "bray", facet_by = "tissue"
) + theme_bw(base_size = 20) + theme(axis.text.x = element_text(angle = 90),
                                     legend.position = "bottom"))

