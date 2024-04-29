setwd("~/Desktop/Smoking/all_phyloseq")
library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggsignif)
library(vegan)


#### Load in RData ####
load("smk_rare.RData")
load("smk_final.RData")

#here, the UK dataset is referred to as "america" due to the data being from Houston, Texas
smk_rare_saliva_america <- subset_samples(smk_rare, Sample_type %in% c("Saliva") & QuickCode %in% c("Houston"))


###### ALPHA  ######################################################
america_gg_richness_o <- plot_richness(smk_rare_saliva_america, x = "Smoking_status", measures = c("Observed")) +
  xlab("Smoking Status") +
  geom_boxplot(aes(fill = Smoking_status)) +
  theme_classic() +
  ylim(0, 350) +
  scale_fill_manual(values = c("darkred","darkblue","darkgreen")) +
  theme(axis.title = element_text(face = "bold", size = rel(1.2)),
        axis.text = element_text(face = "bold", size = rel(1.2)),
        legend.position = "none")
america_gg_richness_o

america_gg_richness_s <- plot_richness(smk_rare_saliva_america, x = "Smoking_status", measures = c("Shannon")) +
  xlab("Smoking Status") +
  geom_boxplot(aes(fill = Smoking_status)) +
  theme_classic() +
  ylim(0, 6) +
  scale_fill_manual(values = c("darkred","darkblue","darkgreen")) +
  theme(axis.title = element_text(face = "bold", size = rel(1.2)),
        axis.text = element_text(face = "bold", size = rel(1.2)),
        legend.position = "none")
america_gg_richness_s

america_gg_richness_c <- plot_richness(smk_rare_saliva_america, x = "Smoking_status", measures = c("Chao1")) +
  xlab("Smoking Status") +
  geom_boxplot(aes(fill = Smoking_status)) +
  theme_classic() +
  ylim(0, 350) +
  scale_fill_manual(values = c("darkred","darkblue","darkgreen")) +
  theme(axis.title = element_text(face = "bold", size = rel(1.2)),
        axis.text = element_text(face = "bold", size = rel(1.2)),
        legend.position = "none")
america_gg_richness_c

america_alphadiv <- estimate_richness(smk_rare_saliva_america)
samp_dat <- sample_data(smk_rare_saliva_america)
samp_dat_wdiv <- data.frame(samp_dat, america_alphadiv)
# wilcox.test(Observed ~ Smoking_status, data=samp_dat_wdiv, exact = FALSE) --> no significant difference in all 3 metrics 
krt_obs_america <- kruskal.test(Observed ~ Smoking_status, data=samp_dat_wdiv) 
lm_obs_america <- lm(log(Observed) ~ `Smoking_status`, data=samp_dat_wdiv)
anova_obs_america <- aov(lm_obs_america)
summary(anova_obs_america)
tukey_sum_obs_america <- TukeyHSD(anova_obs_america)

krt_shannon_america <- kruskal.test(Shannon ~ Smoking_status, data=samp_dat_wdiv) 
lm_shannon_america <- lm(log(Shannon) ~ `Smoking_status`, data=samp_dat_wdiv)
anova_shannon_america <- aov(lm_shannon_america)
summary(anova_shannon_america)
tukey_sum_shannon_america <- TukeyHSD(anova_shannon_america)

krt_chao_america <- kruskal.test(Chao1 ~ Smoking_status, data=samp_dat_wdiv) 
lm_chao_america <- lm(log(Chao1) ~ `Smoking_status`, data=samp_dat_wdiv)
anova_chao_america <- aov(lm_chao_america)
summary(anova_chao_america)
tukey_sum_chao_america <- TukeyHSD(anova_chao_america)

# phylogenetic diversity
# calculate Faith's phylogenetic diversity as PD
america_phylo_dist <- pd(t(otu_table(smk_rare_saliva_america)), phy_tree(smk_rare_saliva_america),
                    include.root=F) 
# add PD to metadata table
sample_data(smk_rare_saliva_america)$PD <- america_phylo_dist$PD
# plot any metadata category against the PD
america_plot.pd <- ggplot(sample_data(smk_rare_saliva_america), aes(Smoking_status, PD)) + 
  geom_boxplot(aes(fill = Smoking_status)) +
  theme_classic() +
  ylim(0, 20) +
  scale_fill_manual(values = c("darkred","darkblue","darkgreen")) +
  theme(axis.title = element_text(face = "bold", size = rel(1.2)),
        axis.text = element_text(face = "bold", size = rel(1.2)),
        legend.position = "none") +
  xlab("Smoking status") +
  ylab("Phylogenetic Diversity")

#significance test for PD --> no significance 
samp_dat_pd <- data.frame(samp_dat, america_phylo_dist)
krt_pd_america <- kruskal.test(PD ~ Smoking_status, data=samp_dat_pd) 
lm_pd_america <- lm(log(PD) ~ `Smoking_status`, data=samp_dat_pd)
anova_pd_america <- aov(lm_pd_america)
summary(anova_pd_america)
tukey_sum_pd_america <- TukeyHSD(anova_pd_america)

# Convert to relative abundance
smk_rare_saliva_america_RA <- transform_sample_counts(smk_rare_saliva_america, function(x) x/sum(x))
# To remove black bars, "glom" by phylum first
smk_rare_saliva_america_RA_phylum <- tax_glom(smk_rare_saliva_america_RA, taxrank = "Phylum", NArm=FALSE)
smk_rare_saliva_america_RA_phylum_taxa <- plot_bar(smk_rare_saliva_america_RA_phylum, fill="Phylum") + 
  facet_wrap(.~Smoking_status, scales = "free_x")

smk_rare_saliva_america_RA_family <- tax_glom(smk_rare_saliva_america_RA, taxrank = "Family", NArm=FALSE)
smk_rare_saliva_america_RA_family_taxa <- plot_bar(smk_rare_saliva_america_RA_family, fill="Family") + 
  facet_wrap(.~Smoking_status, scales = "free_x")

smk_rare_saliva_america_RA_order <- tax_glom(smk_rare_saliva_america_RA, taxrank = "Order", NArm=FALSE)
smk_rare_saliva_america_RA_order_taxa <- plot_bar(smk_rare_saliva_america_RA_order, fill="Order") + 
  facet_wrap(.~Smoking_status, scales = "free_x")



### BETA colored based on smoking status  ######################################################
bc_dm <- distance(smk_rare_saliva_america, method = "bray")
pcoa_bc <- ordinate(smk_rare_saliva_america, method = "PCoA", distance = bc_dm)
gg_pcoa_saliva_america_bcd <- plot_ordination(smk_rare_saliva_america, pcoa_bc, color = "Smoking_status") +
  labs(col = "Smoking Status") + stat_ellipse(type = "norm")
gg_pcoa_saliva_america_bcd

ord_unweighted_unifrac <- ordinate(smk_rare_saliva_america, method="PCoA", distance="unifrac", weighted=FALSE)
gg_pcoa_saliva_america_ununi <- plot_ordination(smk_rare_saliva_america, ord_unweighted_unifrac, color = "Smoking_status") +
  labs(col = "Smoking Status") + stat_ellipse(type = "norm")
gg_pcoa_saliva_america_ununi

ord_weighted_unifrac <- ordinate(smk_rare_saliva_america, method="PCoA", distance="unifrac", weighted=TRUE)
gg_pcoa_saliva_america_weuni <- plot_ordination(smk_rare_saliva_america, ord_weighted_unifrac, color = "Smoking_status") +
  labs(col = "Smoking Status") + stat_ellipse(type = "norm")
gg_pcoa_saliva_america_weuni
