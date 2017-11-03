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
theme_set(theme_bw(base_size = 15)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
file.choose()

##ITS import
metadata=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\UDBiofilms\\UDBF_ITS_mapping.txt",header=TRUE)
otufull=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\UDBiofilms\\UDBF_ITS_mapping.txt",header=TRUE)#L6otu.txt
taxmatrixfull=as.matrix(read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\UDBiofilms\\UDTAXAITS.txt"))#TAXAtableL6.txt
OTU=otu_table(otufull, taxa_are_rows=TRUE)
#head(OTU) #should be 6 taxa and 40 samples
TAX=tax_table(taxmatrixfull)
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$NAME
#head(TAX) #will say six taxa by six ranks, theres more in the file it is just counting the header
#rownames(OTU)
#(TAX)
colnames(TAX) <- c("Domain", "Phylum", "Class", "Order", "Family","Genus")#assigns names to the taxa levels
#taxa_names(TAX)
#row.names(OTU)
taxa_names(TAX)=row.names(OTU)
physeq=phyloseq(OTU,TAX,sampdat)#joins together OTU,TAX, and metadata into a 4D  object



######################################
#########16S
######################################
head(ITSmeta)
file.choose()
biome=import_biom("C:\\Users\\Joe Receveur\\Documents\\MSU data\\UDBiofilms\\UDBF_meta_ns_no_mito.biom",parseFunction= parse_taxonomy_greengenes)
biome
physeq=biome
sample_names(biome)
ITSmeta$SampleID
sample_data(physeq)

physeq
samdat=sample_data(ITSmeta)
head(samdat)
sample_names(samdat)=sample_data(samdat)$SampleID
sample_names(samdat)
ITSnames=sample_names(ITS)
MetaNames

sample_data(physeq)$Treatment2=c("Modified","Reference","Modified","Modified","Modified","Reference","Modified","Modified","Reference","Modified","Modified","Modified","Modified","Modified","Modified","Modified","Modified","Modified","Modified","Modified","Modified","Modified","Modified","Modified","Modified","Modified","Reference","Reference","Reference","Modified","Reference","Modified","Modified","Modified","Modified","Modified","Modified","Reference","Modified","Reference","Modified","Modified","Reference","Modified","Reference","Reference","Modified","Reference","Modified","Reference","Reference","Modified","Modified","Reference","Modified","Modified","Reference","Reference","Reference","Modified","Reference","Modified","Modified","Modified","Modified","Modified","Reference","Modified","Reference","Reference","Reference","Modified","Modified","Reference","Modified","Reference","Modified","Reference","Modified","Reference","Reference","Modified","Modified","Modified","Reference","Modified","Modified","Modified","Modified","Modified","Modified","Reference","Modified","Reference","Modified","Modified","Modified","Modified","Modified","Modified","Modified","Modified","Modified","Reference","Reference","Modified","Modified","Reference","Modified","Modified","Modified","Modified","Modified","Reference","Modified","Reference","Reference","Reference","Reference","Modified","Modified","Modified","Modified","Reference","Reference")
sample_data(physeq)$DateTreat= paste0(sample_data(physeq)$Treatment2," ",sample_data(physeq)$Date) #creates DateTreat variable

physeq
rank_names(physeq)
sample_variables(physeq)
rank_names(physeq)
plot_richness(physeq, x="Date",color="Date", measures=c("Simpson"),)+geom_boxplot(aes(x=Date, y=value, color=Date), alpha=0.05)

GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
GPr=tax_glom(GPr, "Family")
Biofilm = filter_taxa(GPr, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 0.1%
#Biofilm  = transform_sample_counts(Biofilm, function(x) x / sum(x) )#transform samples so they are based on relative abundance
Biofilm

ord=ordinate(physeq,"PCoA", "bray")
ordplot=plot_ordination(physeq, ord,"samples", color="Treatment")+geom_point(size=4)+scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = Treatment))+ theme(legend.justification=c(1,0), legend.position=c(1,0))#+ facet_wrap(~Date)

GPdist=phyloseq::distance(physeq, "jaccard")
adonis(GPdist ~ Treatment, as(sample_data(physeq), "data.frame"))
adonis(GPdist ~Date*Treatment, as(sample_data(physeq), "data.frame"))

Biofilmplot=plot_bar(Biofilm, "Name","Abundance", fill='Phylum') +ylab("Relative Bacterial Abundance (> 0.1%)")+ facet_grid(Treatment ~ .,scale="free")#+scale_fill_manual(values=cbPalette)
Biofilmplot+ theme(axis.text.x = element_text(angle = 0, hjust = 0.5))#+ theme(legend.justification=c(0.05,0.95), legend.position=c(0.05,0.95))
sample_data(physeq)$Name=sample_names(physeq)

#PSMELT
df <- psmelt(Biofilm)
df$Abundance=df$Abundance*100
cdata <- ddply(df, c("Family", "Treatment"), summarise,
               N    = length(Abundance),
               mean = mean(Abundance),
               sd   = sd(Abundance),
               se   = sd / sqrt(N)
)
head(df)
head(cdata)
#cdata$mean
# plot bar graph with standard deviation as error bars

head(cdata)

# barplots! (replace SampleType with treatment)
Plot=ggplot(cdata, aes(x=Treatment, y=mean,fill=Treatment))+facet_wrap(~Family)
Plot+geom_bar(stat="identity")+ geom_errorbar(aes(ymin=mean-se, ymax=mean+se),col="black")+
 theme(axis.text.x = element_text(angle = 45, hjust = 1))+ylab("Relative Abundance (> 1%)")+
library(ggpubr)
ggboxplot(df, x = "Treatment", y = "Abundance",color = "Treatment", palette = "jco",xlab= "Treatment")+ facet_wrap(~Phylum)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))







###################33

GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance

GPrFamily=tax_glom(GPr,"Family")
FamilyLevel = filter_taxa(GPrFamily, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 1%
FamilyLevel


df <- psmelt(FamilyLevel)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Family", "Treatment2","Date"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)



TrtDate=merge_samples(FamilyLevel,"Trt_Date")
sample_data(TrtDate)
sample_data(TrtDate)$Treatment=c("Heavy","Heavy","Heavy","Heavy","Heavy","Heavy","Heavy","Heavy","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Reference","Reference","Reference","Reference","Reference","Reference","Reference","Reference")
sample_data(TrtDate)$Date=c("11/16","11/30","5/17","5/31","7/14","7/28","9/16","9/30","11/16","11/30","5/17","5/31","7/14","7/28","9/16","9/30","11/16","11/30","5/17","5/31","7/14","7/28","9/16","9/30")
sample_data(TrtDate)$Date = factor(sample_data(TrtDate)$Date, levels = c("5/17", "5/31","7/14","7/28","9/16","9/30","11/16","11/30")) #fixes x-axis labels

sample_data(TrtDate)$Instar=sample_names(TrtDate)
#sample_data(Instar)
TrtDate=transform_sample_counts(TrtDate, function(x) 100*x/sum(x)) #merging samples #(averaging)
TrtDate
TrtDateFamilyplot=plot_bar(TrtDate, "Date","Abundance", fill='Family')  +ylab("Relative Bacterial Abundance (> 1%)")+facet_grid(Treatment ~ .)#+scale_fill_manual(values=cbPalette)
TrtDateFamilyplot+ theme(axis.text.x = element_text(angle = 45, hjust = 1))#+ theme(legend.justification=c(0.05,0.95), legend.position=c(0.05,0.95))


#Richness
factor(sample_data(physeq)$Date)
sample_data(physeq)$Date = factor(sample_data(physeq)$Date, levels = c("17-May-16","31-May-16","14-Jul-16","28-Jul-16","16-Sep-16","30-Sep-16","16-Nov-16","30-Nov-16")) #fixes x-axis labels

plot_richness(physeq, x="Date",color="Date", measures=c("Simpson"),)+geom_boxplot(aes(x=Date, y=value, color=Date), alpha=0.05)+facet_wrap(~Treatment)+ylab("Simpson Diversity")
plot_richness(physeq, x="Date",color="Date", measures=c("Shannon"),)+geom_boxplot(aes(x=Date, y=value, color=Date), alpha=0.05)+facet_wrap(~Treatment)+ylab("Shannon Diversity")

file.choose()
Evenness=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\UDBiofilms\\UDBFevennessPielou.txt", header=TRUE)
Evenness$Date = factor(Evenness$Date, levels = c("17-May-16","31-May-16","14-Jul-16","28-Jul-16","16-Sep-16","30-Sep-16","16-Nov-16","30-Nov-16")) #fixes x-axis labels


EvenSum <- ddply(Evenness, c( "Treatment","Date"), summarise,
                   N    = length(Evenness),
                   mean = mean(Evenness),
                   sd   = sd(Evenness),
                   se   = sd / sqrt(N)
)

##Evenness boxplot
ggplot(Evenness, aes(x=Date, y=Evenness,color=Date))+facet_wrap(~Treatment)+geom_point()+
geom_boxplot(aes(x=Date, y=Evenness, color=Date))+ylab("Species Evenness")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

Plot=ggplot(EvenSum, aes(x=Date, y=mean))+facet_grid(Treatment ~ .)###Evenness se bars
Plot+geom_point()+ geom_errorbar(aes(ymin=mean-se, ymax=mean+se),col="black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
HeavyVRef=subset_samples(physeq, Treatment!= "f")
GPdist=phyloseq::distance(physeq, "jaccard")
#sample_variables(physeq)
beta <- betadisper(GPdist, sample_data(physeq)$Treatment)
permutest(beta, pairwise= FALSE, permutations= 999)
beta.HSD=TukeyHSD(beta)
beta.HSD$group
plot(beta, label=TRUE,ellipse=TRUE, hull=FALSE, conf=0.95, segments = FALSE)
boxplot(beta)

sample_variables(physeq)



plot_richness(physeq, x="Date",color="Date", measures=c("Simpson"),)+geom_boxplot(aes(x=Date, y=value, color=Date), alpha=0.05)+facet_wrap(~Treatment2)+ylab("Simpson Diversity")
