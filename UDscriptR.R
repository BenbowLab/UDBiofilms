library(vegan)
library(MASS)
library(ggplot2)
library(plyr)
library(dplyr)
library(magrittr)
#library(mctoolsr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(randomForest)
library(knitr)
library(ape)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000","#CC79A7")
theme_set(theme_bw(base_size = 20)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))



biome=import_biom(file.choose(),parseFunction= parse_taxonomy_greengenes)#,file.choose()

physeq=biome
physeq
sample_variables(physeq)
rank_names(physeq)
plot_richness(physeq, x="Treatment",color="Treatment",shape="Treatment", measures=c("Simpson"),)+geom_boxplot(aes(x=Treatment, y=value, color=Treatment), alpha=0.05)


GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
Biofilm = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE) #filter out any taxa lower tha 0.1%
Biofilm  = transform_sample_counts(Biofilm, function(x) x / sum(x) )#transform samples so they are based on relative abundance
Biofilm <- subset_taxa(Biofilm, Family != "mitochondria" & Class != "Chloroplast")

ord=ordinate(physeq,"PCoA", "bray")
ordplot=plot_ordination(physeq, ord,"samples", color="Treatment", shape="Treatment")+geom_point(size=4)+scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = Treatment))+ theme(legend.justification=c(1,0), legend.position=c(1,0))
  facet_wrap(~Date)

GPdist=phyloseq::distance(physeq, "bray")
adonis(GPdist ~ Treatment, as(sample_data(physeq), "data.frame"))
adonis(GPdist ~Sampling_Date*Treatment, as(sample_data(physeq), "data.frame"))
plot_bar(physeq, "Phylum","Abundance")#, "Sampling_Date",facet_grid="Treatment~.")

  
  