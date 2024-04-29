setwd("~/Desktop/Smoking/all_phyloseq")

library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(microbiome)
library(ggVennDiagram)
library(indicspecies)

#### Load in RData ####
load("smk_rare.RData")
load("smk_final.RData")

#here, the UK dataset is referred to as "america" due to the data being from Houston, Texas
smk_final_saliva_us <- subset_samples(smk_final, Sample_type %in% c("Saliva") & QuickCode %in% c("Houston")) 



#core microbiome analysis
# Convert absolute abundance to relative abundance
smk_RA_us <- transform_sample_counts(smk_final_saliva_us, fun=function(x) x/sum(x))#each counts/total ASV

c_us <- subset_samples(smk_RA_us, `Smoking_status` == "C")
ec_us <- subset_samples(smk_RA_us, `Smoking_status` == "EC")
ns_us <- subset_samples(smk_RA_us, `Smoking_status` == "NS")

c_ASVs_us <- core_members(c_us, detection=0.001, prevalence = 0.05) 
ec_ASVs_us <- core_members(ec_us, detection=0.001, prevalence = 0.05)
ns_ASVs_us <- core_members(ns_us, detection=0.001, prevalence = 0.05)
first_venn_us <- ggVennDiagram(x=list(Cigarette = c_ASVs_us, ECigarette = ec_ASVs_us, Nonsmoking = ns_ASVs_us)) +
  labs(col = "Smoking Status")


# Adjust the plot with theme modifications for font size

first_venn_us <- first_venn_us + 
  theme(text = element_text(size = 20), # Adjusts global text size, change 20 as needed
        legend.text = element_text(size = 20)) # Specifically adjusts legend text size, change 16 as needed

first_venn_us
# Save the modified plot with adjusted font sizes
ggsave(filename = "us_venn_NS_vs_S.png",
       plot = first_venn_us,
       height = 25, width = 25, dpi = 300)


