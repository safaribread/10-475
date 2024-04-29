setwd("~/Desktop/Smoking/all_phyloseq")
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(microbiome)
library(ggVennDiagram)
library(indicspecies)

#### Load in RData ####
load("smk_final_saliva.RData")


smk_unrare_saliva_uk <- subset_samples(smk_final_saliva, QuickCode %in% c("Houston"))


#Indicator Species Analysis

smk_genus <- tax_glom(smk_unrare_saliva_uk, "Species", NArm = FALSE)
smk_RA_uk <- transform_sample_counts(smk_unrare_saliva_uk, fun=function(x) x/sum(x))

ISA_smk_uk <- multipatt(t(otu_table(smk_RA_uk)), cluster = sample_data(smk_RA_uk)$'Smoking_status')
summary(ISA_smk_uk)

taxtable <- tax_table(smk_final_saliva) %>% as.data.frame() %>% rownames_to_column(var="ASV")

ISA_uk_res <- ISA_smk_uk$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) 

View(ISA_uk_res)
