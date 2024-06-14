library(ggplot2)
library(ggpubr)
library(phyloseq)
library(DESeq2)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(microViz)
library(dplyr)
#import data
ovs<-read.csv("ov-neg.csv",check.names = F)[,1:24]
colnames(ovs)
ovstax <- read.csv("ov-neg.csv",check.names = F)[,c(1,26:31)]
colnames(ovstax)
mgs <- read.csv("mg-neg.csv",check.names = F)[,1:24]
colnames(mgs)
mgstax <- read.csv("mg-neg.csv",check.names = F)[,c(1,26:31)]
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
metadata <- read.csv("ctrlmetadata.csv",row.names = 1)
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


#t.test, wilcox.test, anova, kruskal.test
# Calculate distances
#Bray cutis beta diversity
DistBC = phyloseq::distance(ps_object, method = "bray") #jaccard, unifrac, wunifrac, euclidean
ordBC = ordinate(ps_object, method = "PCoA", distance = DistBC)
plot_scree(ordBC, "Scree Plot: Bray-Curtis MDS")
plot_ordination(ps_object, ordBC, color = "tissue") +
  geom_point(aes(color=tissue), size = 4) + 
  theme_bw(base_size = 20) + labs(color="Tissue")
# Jaccard distance beta diversity
DistJD = phyloseq::distance(ps_object, method = "jaccard") #jaccard, unifrac, wunifrac
ordJD = ordinate(ps_object, method = "PCoA", distance = DistJD)
plot_scree(ordJD, "Scree Plot: Jaccard MDS")
plot_ordination(ps_object, ordJD, color = "tissue") +
  geom_point(aes(color=tissue), size = 4) +
  theme_bw(base_size = 20) + labs(color="Tissue")

#############################################################################
#Differential abundance

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

ovs_vs_mgs.df <- as.data.frame(ovs_vs_mgs)

ovs_vs_mgs.df$sig <- ifelse(ovs_vs_mgs.df$padj<0.2,"yes","no") #0.2
ovs_vs_mgs.df$labs <- ifelse(ovs_vs_mgs.df$sig=="yes",rownames(ovs_vs_mgs.df),NA)
ggmaplot(ovs_vs_mgs.df, fdr = 0.05, fc = 2, size = 2,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(ovs_vs_mgs.df$labs),label.rectangle = TRUE,
         font.label = c(19, "plain", "black"),
         top = 4, label.select = c("Asaia","Niallia","Zymobacter","Massilia"),
         xlab = "Log2 mean abundance",ggtheme = theme_bw(base_size = 19))

#        genenames = as.vector(ovs_vs_mgs.df$labs),top = 0,label.rectangle = TRUE,
#         label.select = c("Asaia","Niallia","Zymobacter","Massilia"))


ps.df <- tax_transform(ps_object,"identity", rank = "class") #compositional

#Draw relative abundance plot
(pp <- comp_barplot(ps.df,
                    n_taxa = 8, tax_order = sum, tax_level = "class",
                    bar_outline_colour = "black", merge_other = FALSE,
                    sample_order = "bray", facet_by = "tissue"
) + theme_bw(base_size = 20) + theme(axis.text.x = element_text(angle = 90),
                                     legend.position = "bottom")) 
