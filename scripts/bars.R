library(tidyverse)

seqs <- read.csv("reads.csv", header = TRUE)
reads <- seqs[seqs$group == "ovpos",]

ggplot(reads, aes(sample_id, count, label=count)) + 
  geom_bar(stat = "identity") + theme_classic(base_size = 20) +
  labs(y="Number of reads", x="Samples") + 
  geom_text(position = position_dodge(width = 1),vjust=-0.5,hjust=0.5) + 
  theme(axis.text.x = element_text(angle = 90))
