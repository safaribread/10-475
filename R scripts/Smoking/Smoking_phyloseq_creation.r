setwd("~/Desktop/Smoking/all_phyloseq")

#Importing required R libraries
library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)


#Importing files from qiime into current environment
metafp <- "metadataFull.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "all_merged_unrooted_tree.nwk"
phylotree <- read.tree(phylotreefp)



#creating phyloseq objects out of imports
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

samp_df <- as.data.frame(meta[,-1])
rownames(samp_df)<- meta$'Sample_ID'
SAMP <- sample_data(samp_df)
class(SAMP)

tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix()
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature ID`
TAX <- tax_table(tax_mat)
class(TAX)

#merge into a single phyloseq object
smk <- phyloseq(OTU, SAMP, TAX, tree = phylotree)

sample_variables(smk)
get_variable(smk, c("Smoking_status"))

# Remove non-bacterial sequences, if any; only keep bacteria samples 
smk_filt <- subset_taxa(smk,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# (optional) Remove ASVs that have less than 5 counts total (low abundance samples, ASV<5). If sum >5, keep
smk_filt_nolow <- filter_taxa(smk_filt, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads, quality control as samples <100 reads have poor sampling run 
smk_final <- prune_samples(sample_sums(smk_filt_nolow)>100, smk_filt_nolow)

#rarefaction 
rarecurve(t(as.data.frame(otu_table(smk_final))), cex=0.1) #use a transposed otu column, cex = font size
smk_rare <- rarefy_even_depth(smk_final, rngseed = 1, sample.size = 9500) 

save(smk_final, file="smk_final.RData")
save(smk_rare, file="smk_rare.RData")


#Creating phyloseq objects that only contain saliva and buccal swab samples (manly for beta diversity)
smk_final_saliva_buccal <- subset_samples(smk_final, Sample_type %in% c("Saliva", "Buccal_Swab"))
smk_rare_saliva_buccal <- subset_samples(smk_rare,Sample_type %in% c("Saliva", "Buccal_Swab"))

unique(get_variable(smk_final_saliva_buccal, c("Sample_type")))

save(smk_final_saliva_buccal, file="smk_final_saliva_buccal.RData")
save(smk_rare_saliva_buccal, file="smk_rare_saliva_buccal.RData")


#Creating phyloseq objects that only contain saliva 
smk_final_saliva <- subset_samples(smk_final, Sample_type == "Saliva")
smk_rare_saliva <- subset_samples(smk_rare, Sample_type == "Saliva")

unique(get_variable(smk_rare_saliva, c("Sample_type")))

save(smk_final_saliva, file="smk_final_saliva.RData")
save(smk_rare_saliva, file="smk_rare_saliva.RData")

