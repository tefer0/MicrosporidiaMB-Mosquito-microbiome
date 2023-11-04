library(DESeq2)
library(ggplot2)
library(ggpubr)
library(phyloseq)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(microViz)
library(dplyr)
#import data
pos<-read.csv("ov-pos.csv",check.names = F)[,1:26]
colnames(pos)
postax <- read.csv("ov-pos.csv",check.names = F)[,c(1,28:33)]
colnames(postax)
neg <- read.csv("ov-neg.csv",check.names = F)[,1:24]
colnames(neg)
negtax <- read.csv("ov-neg.csv",check.names = F)[,c(1,26:31)]
colnames(negtax)
#merge the two data sets by species column
mergeddata<-merge(x=pos,y=neg, by="genus",all=T)
mergeddata
mergeddata[is.na(mergeddata)]<-0
mergeddata
rownames(mergeddata)<-mergeddata$genus
mergeddata<-mergeddata[,-1]

#taxa table
mergetax <- rbind(postax, negtax)
mergetax
mergetax <- distinct(mergetax,genus,.keep_all = T)
rownames(mergetax)<-mergetax$genus
mergetax<-mergetax[,-1]

#load metadata
metadata <- read.csv("ovmetadata.csv",row.names = 1)
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
(p <- plot_richness(ps_object, "condition",measures=alpha_meas))
p$layers <- p$layers[-1]
#MMB plot
p + geom_boxplot(data=p$data, aes(x=condition, y=value, fill = condition), 
  alpha=0.2, show.legend = F, outlier.shape = NA) + geom_jitter(width = 0.2) + 
  theme_bw(base_size = 30) + 
  stat_compare_means(method="kruskal.test",label.y.npc = "bottom" ,label.x.npc = "left") +
  labs(x="Microsporidia MB Positivity")
# #MMB species
# p + geom_boxplot(data=p$data, aes(x=condition, y=value, fill = mspecies), 
#                  alpha=0.2, show.legend = T, outlier.shape = NA) + #geom_jitter(width = 0.2) + 
#   theme(axis.text.x = element_text(angle = 0, hjust=0.5)) + theme_pubclean() + 
#   stat_compare_means(method="kruskal.test", label.y = 1.1, 
#                      label.x.npc = "center") + 
#   labs(title="Alpha diversity plots Shannon and Simpson", 
#        fill="Mosquito spp", x="Microsporidia MB Positivity")
# 
# #MMB Site
# p + geom_boxplot(data=p$data, aes(x=condition, y=value, fill = site), 
#                  alpha=0.2, show.legend = T, outlier.shape = NA) + #geom_jitter(width = 0.2) + 
#   theme(axis.text.x = element_text(angle = 0, hjust=0.5)) + theme_pubclean() + 
#   stat_compare_means(method="kruskal.test", label.y = 1.1, 
#                      label.x.npc = "center") + 
#   labs(title="Alpha diversity plots Shannon and Simpson", 
#        fill="Site", x="Microsporidia MB Positivity")
# 
# #MMB loc
# p + geom_boxplot(data=p$data, aes(x=condition, y=value, fill = loc), 
#                  alpha=0.2, show.legend = T, outlier.shape = NA) + #geom_jitter(width = 0.2) + 
#   theme(axis.text.x = element_text(angle = 0, hjust=0.5)) + theme_pubclean() + 
#   stat_compare_means(method="kruskal.test", label.y = 1.1, 
#                      label.x.npc = "center") + 
#   labs(title="Alpha diversity plots Shannon and Simpson", 
#        fill="Location", x="Microsporidia MB Positivity")
# 
# #MMB house
# p + geom_boxplot(data=p$data, aes(x=condition, y=value, fill = house), 
#                  alpha=0.2, show.legend = T, outlier.shape = NA) + #geom_jitter(width = 0.2) + 
#   theme(axis.text.x = element_text(angle = 0, hjust=0.5)) + theme_pubclean() + 
#   stat_compare_means(method="kruskal.test", label.y = 1.1, 
#                      label.x.npc = "center") + 
#   labs(title="Alpha diversity plots Shannon and Simpson", 
#        fill="House Hold ID", x="Microsporidia MB Positivity")

#species
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
plot_ordination(ps_object, ordBC, color = "condition") +
  geom_point(aes(color=condition), size = 4) + 
  theme_bw(base_size = 20) + labs(color="MMB Positivity")

# plot_ordination(ps_object, ordBC, color = "condition", shape = "site") +
#   geom_point(aes(color=condition), size = 4) +
#   ggtitle("PCoA: Bray-Curtis") + labs(color="MMB Positivity",shape="Site")
# 
# plot_ordination(ps_object, ordBC, color = "condition", shape = "loc") +
#   geom_point(aes(color=condition), size = 4) +
#   ggtitle("PCoA: Bray-Curtis") + labs(color="MMB Positivity",shape="Location")
# 
# 
# 
# DistJD = phyloseq::distance(ps_object, method = "jaccard") #jaccard, unifrac, wunifrac
# ordJD = ordinate(ps_object, method = "PCoA", distance = DistJD)
# plot_scree(ordJD, "Scree Plot: Jaccard MDS")
# plot_ordination(ps_object, ordJD, color = "condition", shape = "mspecies") +
#   geom_point(aes(color=condition), size = 4) +
#   ggtitle("PCoA: Jaccard Distance") + labs(color="MMB Positivity",shape="Mosquito Spp")
# 
# plot_ordination(ps_object, ordJD, color = "condition", shape = "site") +
#   geom_point(aes(color=condition), size = 4) +
#   ggtitle("PCoA: Jaccard Distance") + labs(color="MMB Positivity",shape="Site")
# 
# plot_ordination(ps_object, ordJD, color = "condition", shape = "loc") +
#   geom_point(aes(color=condition), size = 4) +
#   ggtitle("PCoA: Jaccard Distance") + labs(color="MMB Positivity",shape="Location")
# #ordBC$vectors

# df<-as.data.frame(ordBC$vectors[,1:2])
# df$condition<-metadata$condition
# df
# 
# ggplot(df, aes(Axis.1, Axis.2, color=condition))+geom_point(size=4)#PERMANOVA

#############################################################################
diagdds = phyloseq_to_deseq2(ps_object, ~ condition)
#diagdds = phyloseq_to_deseq2(ps_object, ~ loc)
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
#plotDispEsts(diagdds)

pos_vs_neg = results(diagdds,#filterFun = ihw,
                        independentFiltering = FALSE,
                        test = "Wald",
                        name = "condition_positive_vs_negative")
summary(pos_vs_neg)

pos_vs_neg<-pos_vs_neg[order(pos_vs_neg$padj),]
pos_vs_neg
# pos_vs_neg_sig<-subset(pos_vs_neg, pos_vs_neg$padj<0.05)
# pos_vs_neg_sig
# 
# 
# rpm<-fpm(diagdds)
# rpm
# 
# #heatmap(rpm[rownames(pos_vs_neg_sig),])

# Heatmap(t(scale(t(log2(rpm[rownames(pos_vs_neg_sig),]+2)))), 
#         column_split = metadata$condition, show_column_dend = FALSE, 
#         show_row_dend = FALSE ) #condition, loc

#BiocManager::install("ComplexHeatmap")
#BiocManager::install("EnhancedVolcano")
# plot volcano, MA plots

pos_vs_neg.df <- as.data.frame(pos_vs_neg)
# EnhancedVolcano(pos_vs_neg.df,x="log2FoldChange",y="padj",lab = row.names(pos_vs_neg.df),
#                 FCcutoff = 2,pCutoffCol = "padj", pCutoff = 0.05)
pos_vs_neg.df$sig <- ifelse(pos_vs_neg.df$padj<0.05,"yes","no")
pos_vs_neg.df$labs <- ifelse(pos_vs_neg.df$sig=="yes",rownames(pos_vs_neg.df),NA)
# ggplot(pos_vs_neg.df, aes(log2FoldChange, -log10(padj), color=sig, label=labs)) + 
#   geom_point() + geom_text(col="red")


ggmaplot(pos_vs_neg.df, fdr = 0.05, fc = 2, size = 2,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(pos_vs_neg.df$labs),label.rectangle = TRUE,
         font.label = c(19, "plain", "black"),
         xlab = "Log2 mean abundance",ggtheme = theme_bw(base_size = 19))



ps.df <- tax_transform(ps_object,"identity", rank = "class") #compositional


(pp <- comp_barplot(ps.df,
                    n_taxa = 41, tax_order = sum, tax_level = "class",
                    bar_outline_colour = "black", merge_other = FALSE,
                    sample_order = "bray", facet_by = "condition"
) + theme_bw(base_size = 20) + theme(axis.text.x = element_text(angle = 90),
                                      legend.position = "bottom")) 
