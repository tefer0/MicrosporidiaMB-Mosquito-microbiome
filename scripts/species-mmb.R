#genera
library(tidyverse)
library(normalize)
library(reshape2)

gen <- read.csv("specs.csv")
colnames(gen)
tissue <- read.csv("specs-tissue.csv")
tissue

#ggplot(aes())
gen[,c(2:5)] <- round(normalize(gen[,c(2:5)], center = FALSE, scale =TRUE), 
                      digits = 2)
gen
mgen <- melt(gen,id = "taxa")
mgen
all(mgen$taxa%in%tissue$taxa)
all(mgen$taxa==tissue$taxa)

mgen$tissue <- tissue$tissue
mgen

theme_set(theme_bw(base_size = 20) + theme(legend.position = "bottom"))
# ggplot(mgen, aes(taxa, value, fill=variable, label=value)) +
#   geom_bar(stat = "identity", position = "dodge") + 
#   theme_bw(base_size = 20) + 
#   labs(x="Bacteria Genera",y="Normalised Counts",fill="Mosquito Tissue(MMB +/-)") + 
#   geom_text(position = position_dodge(width = 1), vjust = -0.5, hjust=0.5, 
#              size = 5) 

ggplot(mgen, aes(x=taxa, y=value, fill=variable, label=value)) +
  geom_bar(stat = "identity", position = "dodge") + facet_wrap(~tissue) +
  labs(x="Bacteria Genera",y="Normalised Counts",fill="Mosquito Tissue(MMB positivity)") + 
  geom_text(position = position_dodge(width = 1), vjust = -0.5, hjust=0.5, 
            size = 4)
###################################################################################

# gen <- read.csv("specs-mg.csv")
# colnames(gen)
# gen[,c(2:3)] <- round(normalize(gen[,c(2:3)], center = FALSE, scale =TRUE), 
#                       digits = 2)
# gen
# mgen <- melt(gen,id = "taxa")
# mgen
# 
# ogen <- read.csv("specs-ov.csv")
# colnames(ogen)
# ogen[,c(2:3)] <- round(normalize(ogen[,c(2:3)], center = FALSE, scale =TRUE), 
#                        digits = 2)
# ogen
# omgen <- melt(ogen,id = "taxa")
# omgen
# 
# par(mfrow = c(2,2))
# theme_set(theme_bw(base_size = 20) + theme(axis.text.x = element_text(angle = 90),legend.position = "bottom"))
# 
# ggplot(mgen, aes(x=taxa, y=value, fill=variable, label=value)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(x="Bacteria Genera",y="Normalised Counts",fill="Mosquito Tissue(MMB positivity)") + 
#   geom_text(position = position_dodge(width = 1), vjust = -0.5, hjust=0.5, 
#             size = 5)
# 
# 
# 
# 
# 
# #theme_set(theme_bw(base_size = 20) + theme(legend.position = "bottom"))
# 
# ggplot(mgen, aes(x=taxa, y=value, fill=variable, label=value)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(x="Bacteria Genera",y="Normalised Counts",fill="Mosquito Tissue(MMB positivity)") + 
#   geom_text(position = position_dodge(width = 1), vjust = -0.5, hjust=0.5, 
#             size = 5)
