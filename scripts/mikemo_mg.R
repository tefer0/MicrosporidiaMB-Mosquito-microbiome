library(DESeq2)
library(ggplot2)
library(ggpubr)
library(phyloseq)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(dplyr)
library(microViz)
#import data
posm<-read.csv("mg-pos.csv",check.names = F)[,1:24]
colnames(posm)
postaxm<-read.csv("mg-pos.csv",check.names = F)[,c(1,26:31)]
colnames(postaxm)
negm <- read.csv("mg-neg.csv",check.names = F)[,1:24]
colnames(negm)
negtaxm <- read.csv("mg-neg.csv",check.names = F)[,c(1,26:31)]
colnames(negtaxm)

#merge the two data sets by species column
mergeddatam<-merge(x=posm,y=negm, by="genus",all=T)
mergeddatam
mergeddatam[is.na(mergeddatam)]<-0
mergeddatam
rownames(mergeddatam)<-mergeddatam$genus
mergeddatam<-mergeddatam[,-1]

#merge taxa
# mergetaxm<-merge(x=negtaxm,y=postaxm, by="species",all=T)
# mergetaxm
# mergetaxm[is.na(mergetaxm)]<-0
# mergetaxm
mergetaxm <- rbind(negtaxm, postaxm)
mergetaxm
mergetaxm <- distinct(mergetaxm,genus,.keep_all = T)
rownames(mergetaxm)<-mergetaxm$genus
mergetaxm<-mergetaxm[,-1]

#load metadata
metadatam <- read.csv("mgmetadata.csv",row.names = 1)
metadatam

#check if metadata and row names are present and in same order
all(rownames(metadatam)%in%colnames(mergeddatam))
all(rownames(metadatam)==colnames(mergeddatam))
#metadatam<-metadatam[colnames(mergeddatam),drop=F,]
#all(rownames(metadatam)==colnames(mergeddatam))

all(rownames(mergetaxm)%in%rownames(mergeddatam))
all(rownames(mergetaxm)==rownames(mergeddatam))
#mergetaxm<-mergetaxm[rownames(mergeddatam),drop=F,]

# create phyloseq object
#library(phyloseq)
otu1m<-otu_table(mergeddatam, taxa_are_rows = T)
samplem<-phyloseq::sample_data(metadatam)
taxm <- tax_table(mergetaxm)
rownames(taxm) <- rownames(mergetaxm)
colnames(taxm) <- colnames(mergetaxm)
ps_objectm<-phyloseq(otu1m, samplem, taxm)
ps_objectm
#any(taxa_sums(ps_objectm)==0)

#estimate alpha diversity
alpha_meas = c("Shannon", "Simpson")
#?plot_richness
#MMB
(pm <- plot_richness(ps_objectm, "condition",measures=alpha_meas))
pm$layers <- pm$layers[-1]
#MMB plot
pm + geom_boxplot(data=pm$data, aes(x=condition, y=value, fill = condition), 
  alpha=0.2, show.legend = F, outlier.shape = NA) + geom_jitter(width = 0.2) + 
  theme_bw(base_size = 30) + 
  stat_compare_means(method="kruskal.test",label.y.npc = "bottom",label.x.npc = "center") +
  labs(x="Microsporidia MB Positivity")
#MMB species
# pm + geom_boxplot(data=pm$data, aes(x=condition, y=value, fill = mspecies), 
#                  alpha=0.2, show.legend = T, outlier.shape = NA) + #geom_jitter(width = 0.2) + 
#   theme(axis.text.x = element_text(angle = 0, hjust=0.5)) + theme_pubclean() + 
#   stat_compare_means(method="kruskal.test",label.y.npc = "bottom", 
#                      label.x.npc = "center") + 
#   labs(title="Alpha diversity plots Shannon and Simpson", 
#        fill="Mosquito spp", x="Microsporidia MB Positivity")
# 
# #MMB Site
# pm + geom_boxplot(data=pm$data, aes(x=condition, y=value, fill = site), 
#                  alpha=0.2, show.legend = T, outlier.shape = NA) + #geom_jitter(width = 0.2) + 
#   theme(axis.text.x = element_text(angle = 0, hjust=0.5)) + theme_pubclean() + 
#   stat_compare_means(method="kruskal.test", label.y.npc = "bottom", 
#                      label.x.npc = "center") + 
#   labs(title="Alpha diversity plots Shannon and Simpson", 
#        fill="Site", x="Microsporidia MB Positivity")
# 
# #MMB loc
# pm + geom_boxplot(data=pm$data, aes(x=condition, y=value, fill = loc), 
#                  alpha=0.2, show.legend = T, outlier.shape = NA) + #geom_jitter(width = 0.2) + 
#   theme(axis.text.x = element_text(angle = 0, hjust=0.5)) + theme_pubclean() + 
#   stat_compare_means(method="kruskal.test", label.y.npc = "bottom", 
#                      label.x.npc = "center") + 
#   labs(title="Alpha diversity plots Shannon and Simpson", 
#        fill="Location", x="Microsporidia MB Positivity")
# 
# #MMB house
# pm + geom_boxplot(data=pm$data, aes(x=condition, y=value, fill = house), 
#                  alpha=0.2, show.legend = T, outlier.shape = NA) + #geom_jitter(width = 0.2) + 
#   theme(axis.text.x = element_text(angle = 0, hjust=0.5)) + theme_pubclean() + 
#   stat_compare_means(method="kruskal.test", label.y.npc = "bottom", 
#                      label.x.npc = "center") + 
#   labs(title="Alpha diversity plots Shannon and Simpson", 
#        fill="House Hold ID", x="Microsporidia MB Positivity")

#species
# (pm <- plot_richness(ps_objectm, "mspecies",measures=alpha_meas))
# pm$layers <- pm$layers[-1]
# 
# pm + geom_boxplot(data=pm$data, aes(x=mspecies, y=value, fill = mspecies), 
#                  alpha=0.2, show.legend = T, outlier.shape = NA) + 
#   geom_jitter(width = 0.2) + theme_pubclean() + 
#   stat_compare_means(method="kruskal.test", label.y.npc = "bottom", 
#                      label.x.npc = "left") + 
#   labs(title="Alpha diversity plots Shannon and Simpson",x="Mosquito Species",
#        fill="Mosquito species")
# #site
# (pm <- plot_richness(ps_objectm, "site",measures=alpha_meas))
# pm$layers <- pm$layers[-1]
# 
# pm + geom_boxplot(data=pm$data, aes(x=site, y=value, fill = site), 
#                  alpha=0.2, show.legend = F, outlier.shape = NA) + geom_jitter(width = 0.2) + 
#   theme(axis.text.x = element_text(angle = 0, hjust=0.5)) + theme_pubclean() + 
#   stat_compare_means(method="kruskal.test", label.y.npc = "bottom", 
#                      label.x.npc = "left") + 
#   labs(title="Alpha diversity plots Shannon and Simpson",x="Sample Collection Sites")
# 
# #loc
# (pm <- plot_richness(ps_objectm, "loc",measures=alpha_meas))
# pm$layers <- pm$layers[-1]
# 
# pm + geom_boxplot(data=pm$data, aes(x=loc, y=value, fill = loc), 
#                  alpha=0.2, show.legend = F, outlier.shape = NA) + geom_jitter(width = 0.2) + 
#   theme(axis.text.x = element_text(angle = 0, hjust=0.5)) + theme_pubclean() + 
#   stat_compare_means(method="kruskal.test", label.y.npc = "bottom", 
#                      label.x.npc = "center") + 
#   labs(title="Alpha diversity plots Shannon and Simpson",x="Location indoor or outdoor")
# 
# #House
# (pm <- plot_richness(ps_objectm, "house",measures=alpha_meas))
# pm$layers <- pm$layers[-1]
# 
# pm + geom_boxplot(data=pm$data, aes(x=house, y=value, fill = house), 
#                  alpha=0.2, show.legend = F, outlier.shape = NA) + geom_jitter(width = 0.2) + 
#   theme(axis.text.x = element_text(angle = 0, hjust=0.5)) + theme_pubclean() + 
#   stat_compare_means(method="kruskal.test", label.y.npc = "bottom", 
#                      label.x.npc = "center") + 
#   labs(title="Alpha diversity plots Shannon and Simpson",x="House hold ID")

#t.test, wilcox.test, anova, kruskal.test
# Calculate distances
#Bray cutis
DistBC = phyloseq::distance(ps_objectm, method = "bray") #jaccard, unifrac, wunifrac
ordBC = ordinate(ps_objectm, method = "PCoA", distance = DistBC)
plot_scree(ordBC, "Scree Plot: Bray-Curtis MDS")
plot_ordination(ps_objectm, ordBC, color = "condition") +
  geom_point(aes(color=condition), size = 4) +
  theme_bw(base_size = 20) + labs(color="MMB Positivity")

# plot_ordination(ps_objectm, ordBC, color = "condition", shape = "site") +
#   geom_point(aes(color=condition), size = 4) +
#   ggtitle("PCoA: Bray-Curtis") + labs(color="MMB Positivity",shape="Site")
# 
# plot_ordination(ps_objectm, ordBC, color = "condition", shape = "loc") +
#   geom_point(aes(color=condition), size = 4) +
#   ggtitle("PCoA: Bray-Curtis") + labs(color="MMB Positivity",shape="Location")
# 
# # plot_ordination(ps_objectm, ordBC, color = "house", shape = "condition") +
# #   geom_point(aes(color=house), size = 4) +
# #   ggtitle("PCoA: Bray-Curtis") + labs(shape="MMB Positivity",color="House hold")
# 
# 
# #Jaccard Distance
# DistJD = phyloseq::distance(ps_objectm, method = "jaccard") #jaccard, unifrac, wunifrac
# ordJD = ordinate(ps_objectm, method = "PCoA", distance = DistJD)
# plot_scree(ordJD, "Scree Plot: Jaccard MDS")
# plot_ordination(ps_objectm, ordJD, color = "condition", shape = "mspecies") +
#   geom_point(aes(color=condition), size = 4) +
#   ggtitle("PCoA: Jaccard Distance") + labs(color="MMB Positivity",shape="Mosquito Spp")
# 
# plot_ordination(ps_objectm, ordJD, color = "condition", shape = "site") +
#   geom_point(aes(color=condition), size = 4) +
#   ggtitle("PCoA: Jaccard Distance") + labs(color="MMB Positivity",shape="Site")
# 
# plot_ordination(ps_objectm, ordJD, color = "condition", shape = "loc") +
#   geom_point(aes(color=condition), size = 4) +
#   ggtitle("PCoA: Jaccard Distance") + labs(color="MMB Positivity",shape="Location")
# #ordBC$vectors

# df<-as.data.frame(ordBC$vectors[,1:2])
# df$condition<-metadata$condition
# df
# 
# ggplot(df, aes(Axis.1, Axis.2, color=condition))+geom_point(size=4)#PERMANOVA

#############################################################################
diagddsm = phyloseq_to_deseq2(ps_objectm, ~ condition)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
library(DESeq2)
geoMeansm = apply(counts(diagddsm), 1, gm_mean)
diagddsm = estimateSizeFactors(diagddsm, geoMeans = geoMeansm)
diagddsm$sizeFactor

# diagddsm = phyloseq_to_deseq2(pstop20, ~ condition)
# geoMeansm = apply(counts(diagddsm), 1, gm_mean)
# diagddsm = estimateSizeFactors(diagddsm, geoMeans = geoMeansm)
# diagddsm$sizeFactor

######################################################################################
diagddsm = DESeq(diagddsm,
                fitType="local",
                test = "LRT",
                reduced = ~1)
resultsNames(diagddsm)

pos_vs_negm = results(diagddsm,#filterFun = ihw,
                        independentFiltering = FALSE,
                        test = "Wald",
                        name = "condition_positive_vs_negative")
summary(pos_vs_negm)

pos_vs_negm<-pos_vs_negm[order(pos_vs_negm$padj),]
pos_vs_negm
# pos_vs_neg_sigm<-subset(pos_vs_negm, pos_vs_negm$padj<0.05)
# pos_vs_neg_sigm
# 
# 
# rpm<-fpm(diagddsm)
# rpm
# 
# #heatmap(rpm[rownames(pos_vs_neg_sigm),])
# library(ComplexHeatmap)
# Heatmap(t(scale(t(log2(rpm[rownames(pos_vs_neg_sigm),]+2)))), 
#         column_split = metadatam$condition, show_column_dend = FALSE, 
#         show_row_dend = FALSE) #condition, site
# # Heatmap(t(scale(t(log2(rpm[rownames(pos_vs_negm),]+2)))), 
#         column_split = metadatam$condition)

#BiocManager::install("ComplexHeatmap")
#BiocManager::install("EnhancedVolcano")
# plot volcano, MA plots
# library(EnhancedVolcano)
pos_vs_negm.df <- as.data.frame(pos_vs_negm)
# EnhancedVolcano(pos_vs_negm.df,x="log2FoldChange",y="padj",lab = row.names(pos_vs_negm.df),
#                 FCcutoff = 2,pCutoffCol = "padj", pCutoff = 0.1)
pos_vs_negm.df$sig <- ifelse(pos_vs_negm.df$padj<0.05,"yes","no")
pos_vs_negm.df$labs <- ifelse(pos_vs_negm.df$sig=="yes",rownames(pos_vs_negm.df),NA)
ggmaplot(pos_vs_negm.df, fdr = 0.05, fc = 2, size = 2,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(pos_vs_negm.df$labs),label.rectangle = TRUE,
         font.label = c(19, "plain", "black"),
         xlab = "Log2 mean abundance",ggtheme = theme_bw(base_size = 19))

# ggplot(pos_vs_negm.df, aes(log2FoldChange, -log10(padj), color=sig)) + 
#   geom_point()

#barplot top 20
# 
# top20 <- names(sort(taxa_sums(ps_objectm), decreasing = T))[1:300]
# pstop20 <- transform_sample_counts(ps_objectm, function(OTU) OTU/sum(OTU))
# pstop20 <- prune_taxa(top20, pstop20)
# #pstop20 <- prune_samples(sample_sums(pstop20) > 10, pstop20)
# plot_bar(pstop20, facet_grid = ~condition, fill="class")
# plot_bar(ps_objectm, facet_grid = ~condition, fill="genus")
# (p <- plot_heatmap(pstop20, sample.label = "condition",sample.order="condition",
#              taxa.label = "class", taxa.order = "class", title = ))
# p$scales$scales[[1]]$name <- "Microsporidia MB positivity"
# p$scales$scales[[2]]$name <- "Class of bacteria"
# p
# Heatmap(pstab, column_title=c("Negative samples","Positive samples"),
#         row_title_side=c("left", "right"),row_labels = mergetaxm$class,
#         column_split = metadatam$condition)
# pstab <- otu_table(pstop20)
# pstab <- as.data.frame(pstab)
# heatmap(pstab)
# install.packages("remotes")
# remotes::install_github("david-barnett/microViz")
# 
# ps_objectm %>% tax_transform("compositional", rank = "class") %>% 
#   comp_heatmap(sample_names_show = T,sample_anno =
#   sampleAnnotation(MMB_Positivity = anno_sample_cat("condition", legend_title = 
#                                                "MMB positivity")))

ps.df <- tax_transform(ps_objectm,"identity", rank = "class") #compositional
# ps.df1 <- tax_transform(ps_objectm, "compositional", rank = "class")
# ps.df1 <- tax_filter(ps.df1, min_prevalence = 5)
# subset_samples(ps.df,condition == "negative")
# #ps.df <- tax_filter(ps.df, min_prevalence = 3)
# ps.dff <- otu_table(ps.df)
# ps.dff <- as.matrix(ps.dff)
# # Heatmap(ps.dff,name="Composition",row_km=2,row_km_repeats=100,
# #                 column_split = metadatam$condition)
# Heatmap(ps.dff, name="Counts per sample", column_split = metadatam$condition,
# show_column_dend = FALSE, show_row_dend = FALSE,row_dend_reorder = FALSE,
# left_annotation = row_ha)
# 
# row_ha <- rowAnnotation("Total Counts"=row_anno_barplot(taxa_sums(ps.dff)), 
#                         width = unit(4,"cm"))
# row_ha
# plot_bar(ps.df1,fill = "class", facet_grid = ~condition)
# (p <- plot_heatmap(ps.df, sample.label = "condition",sample.order="condition",
#                                  taxa.label = "class", taxa.order = NULL))
 (pp <- comp_barplot(ps.df,
                  n_taxa = 41, tax_order = sum, tax_level = "class",
                  bar_outline_colour = "black", merge_other = FALSE,
                  sample_order = "bray", facet_by = "condition"
) + theme_bw(base_size = 20) + theme(axis.text.x = element_text(angle = 90),
                                     legend.position = "bottom")) 
# library(patchwork)
# patch <- patchwork::wrap_plots(pp, ncol = 2, guides = "collect")
# patch
