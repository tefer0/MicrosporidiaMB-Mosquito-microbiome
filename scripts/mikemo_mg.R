library(DESeq2)
library(ggplot2)
library(ggpubr)
library(phyloseq)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(dplyr)
library(microViz)
#import data
posm<-read.csv("mg-pos.csv",check.names = F)[,1:25]
colnames(posm)
postaxm<-read.csv("mg-pos.csv",check.names = F)[,c(1,27:32)]
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

#t.test, wilcox.test, anova, kruskal.test
# Calculate distances
#Bray cutis beta diversity
DistBC = phyloseq::distance(ps_objectm, method = "bray") #jaccard, unifrac, wunifrac
ordBC = ordinate(ps_objectm, method = "PCoA", distance = DistBC)
plot_scree(ordBC, "Scree Plot: Bray-Curtis MDS")
plot_ordination(ps_objectm, ordBC, color = "condition") +
  geom_point(aes(color=condition), size = 4) +
  theme_bw(base_size = 20) + labs(color="MMB Positivity")

#Jaccard Distance beta diversity
DistJD = phyloseq::distance(ps_objectm, method = "jaccard") #jaccard, unifrac, wunifrac
ordJD = ordinate(ps_objectm, method = "PCoA", distance = DistJD)
plot_scree(ordJD, "Scree Plot: Jaccard MDS")
plot_ordination(ps_objectm, ordJD, color = "condition") +
  geom_point(aes(color=condition), size = 4) +
  theme_bw(base_size = 20) + labs(color="MMB Positivity")


#############################################################################
#Differential abundance 

diagddsm = phyloseq_to_deseq2(ps_objectm, ~ condition)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
library(DESeq2)
geoMeansm = apply(counts(diagddsm), 1, gm_mean)
diagddsm = estimateSizeFactors(diagddsm, geoMeans = geoMeansm)
diagddsm$sizeFactor

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

pos_vs_negm.df <- as.data.frame(pos_vs_negm)

pos_vs_negm.df$sig <- ifelse(pos_vs_negm.df$padj<0.05,"yes","no")
pos_vs_negm.df$labs <- ifelse(pos_vs_negm.df$sig=="yes",rownames(pos_vs_negm.df),NA)
ggmaplot(pos_vs_negm.df, fdr = 0.05, fc = 2, size = 2,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(pos_vs_negm.df$labs),label.rectangle = TRUE,
         font.label = c(19, "plain", "black"),
         xlab = "Log2 mean abundance",ggtheme = theme_bw(base_size = 19))

#Relative abundance barplot
ps.df <- tax_transform(ps_objectm,"identity", rank = "class") #compositional

 (pp <- comp_barplot(ps.df,
                  n_taxa = 8, tax_order = sum, tax_level = "class",
                  bar_outline_colour = "black", merge_other = FALSE,
                  sample_order = "bray", facet_by = "condition"
) + theme_bw(base_size = 20) + theme(axis.text.x = element_text(angle = 90),
                                     legend.position = "bottom")) 

