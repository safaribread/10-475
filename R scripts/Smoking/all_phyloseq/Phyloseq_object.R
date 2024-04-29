library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)

metafp <- "metadataFull.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")


phylotreefp <- "all_merged_unrooted_tree.nwk"
phylotree <- read.tree(phylotreefp)

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
TAX <- tax_table(tax_mat[,-1])
rownames(TAX)<- tax$`Feature ID`
class(TAX)
smk <- phyloseq(OTU, SAMP, TAX, tree = phylotree)

View(smk)
sample_variables(smk)
otu_table(smk)
sample_data(smk)
tax_table(smk)
phy_tree(smk)