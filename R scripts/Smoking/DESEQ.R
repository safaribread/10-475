setwd("~/Desktop/Smoking/all_phyloseq")

library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(microbiome)
library(ggVennDiagram)
library(indicspecies)
library(DESeq2)


#### Load in RData ####
load("smk_rare.RData")
load("smk_final.RData")


#DESEQ

smk_plus1 <- transform_sample_counts(smk_final, function(x) x+1)
smk_deseq <- phyloseq_to_deseq2(smk_plus1, ~`Smoking_status`)
DESEQ_smk <- DESeq(smk_deseq)

#Smoking vs non smoking
res_C_NS <- results(DESEQ_smk, tidy=TRUE, 
                    contrast = c("Smoking_status","C","NS"))

res_EC_NS <- results(DESEQ_smk, tidy=TRUE, 
                     contrast = c("Smoking_status","EC","NS"))

vol_plot_C_NS <- res_C_NS %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

vol_plot_C_NS
ggsave(filename="vol_plot_C_NS.png",vol_plot_C_NS)


vol_plot_EC_NS <- res_EC_NS %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

vol_plot_EC_NS
ggsave(filename="vol_plot_EC_NS.png",vol_plot_EC_NS)


#bar plots_C_NS
sigASVs_C_NS <- res_C_NS %>% 
  filter(padj<0.001 & abs(log2FoldChange)>2.5) %>%
  dplyr::rename(ASV=row)
View(sigASVs_C_NS)
# Get only asv names
sigASVs_vec_C_NS <- sigASVs_C_NS %>%
  pull(ASV)

smk_DESeq_C_NS <- prune_taxa(sigASVs_vec_C_NS, smk_final)
sigASVs_C_NS <- tax_table(smk_DESeq_C_NS) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_C_NS) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))


Deseqbar_C_NS<-ggplot(sigASVs_C_NS) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

Deseqbar_C_NS
ggsave(filename="Deseqbar_C_NS.png",Deseqbar_C_NS)
#___

sigASVs_EC_NS <- res_EC_NS %>% 
  filter(padj<0.00000000000001 & abs(log2FoldChange)>3.5) %>%
  dplyr::rename(ASV=row)
View(sigASVs_EC_NS)
# Get only asv names
sigASVs_vec_EC_NS <- sigASVs_EC_NS %>%
  pull(ASV)

smk_DESeq_EC_NS <- prune_taxa(sigASVs_vec_EC_NS, smk_final)
sigASVs_EC_NS <- tax_table(smk_DESeq_EC_NS) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_EC_NS) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))


Deseqbar_EC_NS<-ggplot(sigASVs_EC_NS) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

Deseqbar_EC_NS

ggsave(filename="Deseqbar_EC_NS.png",Deseqbar_EC_NS)

