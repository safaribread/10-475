setwd("~/Desktop/Smoking/all_phyloseq")
library(phyloseq)
library(ape)
library(tidyverse)
library(picante)

#### Load in RData ####
load("smk_rare_saliva.RData")
smk_rare_saliva_uk <- subset_samples(smk_rare_saliva, QuickCode %in% c("Houston"))


#### Taxonomy bar plots ####

# Plot bar plot of taxonomy
plot_bar(smk_rare_saliva_uk, fill="Phylum") 

# Convert to relative abundance
smk_RA <- transform_sample_counts(smk_rare_saliva_uk, function(x) x/sum(x))

# To remove black bars, "glom" by phylum first
smk_phylum <- tax_glom(smk_RA, taxrank = "Phylum", NArm=FALSE)

plot_bar(smk_phylum, fill="Phylum") + 
  facet_wrap(.~Smoking_status, scales = "free_x")

gg_taxa <- plot_bar(smk_phylum, fill="Phylum") + 
  facet_wrap(.~Smoking_status, scales = "free_x")
gg_taxa

ggsave("Metadata_Full_Results/plot_taxonomy.png"
       , gg_taxa
       , height=8, width =49)
