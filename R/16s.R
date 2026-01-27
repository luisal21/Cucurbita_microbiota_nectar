### Análisis de las secuencias de 16s
# Luis Alberto Villanueva Espino

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("phyloseq")


# Cargando paquete phyloseq y paquetes adicionales
library(phyloseq)
library(ggplot2)
#library(VennDiagram) # no instalado
library(picante)
library(ape)
library("microbiome")
library("dplyr")
#library("hrbrthemes") # no instalado
#library("gcookbook") # no instalado
#library("tidyverse")
library("RColorBrewer")
library("ggpubr")
#library("breakaway") # no instalado
#library(MicrobiotaProcess)
#library(reltools)
library(microbiomeutilities)
library(vegan)
library(gridExtra)
library(tidyr)
library(reshape2)
#library(UpSetR) # no instalado
#library(minpack.lm) # no instalado
library(Hmisc)
library(stats4)
#library(patchwork) # no instalado


#### Set up directory in which the databases are located ####
setwd("/home/luis/Documentos/papers_publicados/paper2_tesisDoctorado/Cucurbita_microbiota_nectar/data")


### Primero transponer la tabla original "seqtab_nochim_16s_2024.csv, despues
# deshabilitar este codigo
### Transponer tabla de ASV
#asvtable <- read.csv("seqtab_nochim_16s_2024.csv", header = T, row.names = 1)
#asv2 = t(asvtable)
#head(asv2, 1)
#write.csv(asv2, file = "asv_table.csv")

# cargar primero el archivo de metadatos
meta_16s <- read.csv("meta_16s_2024.csv", header = T)
head(meta_16s)
# despues cargar la tabla de ASV
seqtab.16s <- read.csv("asv_table.csv", header = T, row.names = 1)
#head(seqtab.16s)
#tail(seqtab.16s)
str(seqtab.16s)
# Debido a que al cargar el archivo agrega una X al principio del nombre de
# cada columna, hay que decirle que los nombres de las columnas son iguales
# a las filas de sample_id de metadatos para corregir este error
colnames(seqtab.16s) = meta_16s$sample_id
str(seqtab.16s)
# Volver a cargar archivo de metadatos
meta_16s2 <- read.csv("meta_16s_2024.csv", header = T, row.names = 1)
head(meta_16s2)
# cargando taxonomia
taxo_16s <-read.csv("taxa_16s_species_2024.csv", header = T,row.names = 1)
head(taxo_16s)
# corrobando match entre archivos
setdiff(rownames(taxo_16s), rownames(seqtab.16s))
all(rownames(seqtab.16s) == rownames(taxo_16s))
all(colnames(seqtab.16s) == unlist(meta_16s2[,'sample.id']))

# Construct phyloseq object (straightforward from dada2 outputs)
s16_1 <- phyloseq(otu_table(as.matrix(seqtab.16s), taxa_are_rows=T),
                  tax_table(as.matrix(taxo_16s)), 
                  sample_data(meta_16s2))
s16_1
summarize_phyloseq(s16_1) #13,732,428 reads, 4160 ASVs

# deleting bad quality samples with low reads
s16_1 = subset_samples(s16_1, !sampleid2 %in% c("281022-EMo22-SA014",
                      "281022-EMo22-SA035", "281022-EMo22-SD064",
                      "281022-EMo22-SD085", "EMo1223S-SC273",
                      "ENESMo24-bLVEv416s82"))
s16_1
# Building phylogenetic tree with the function rtree of package "ape"
arbol_16s = rtree(ntaxa(s16_1), rooted=TRUE, tip.label=taxa_names(s16_1))
#plot(arbol_16s)
# Uniendo el arbol con el objeto phyloseq anterior
s16_2 = merge_phyloseq(s16_1, arbol_16s)
s16_2
summarize_phyloseq(s16_2) #13,731,996 reads, 4160 ASVs in 356 samples
                          

# Remove contaminants using the frecuency method of the package decontam
library(decontam)

head(sample_data(s16_2))

df <- as.data.frame(sample_data(s16_2)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(s16_2)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=sample_or_control)) + geom_point()

# Identify Contaminants -Frequency
contamdf.freq <- isContaminant(s16_2, method="frequency", conc="conce_mg")
head(contamdf.freq)
View(contamdf.freq)
table(contamdf.freq$contaminant)
head(which(contamdf.freq$contaminant))

# Identify Contaminants -Prevalence
sample_data(s16_2)$is.neg <- sample_data(s16_2)$sample_or_control == "Control Sample"
contamdf.prev <- isContaminant(s16_2, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant), 20)
View(contamdf.prev)

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(s16_2, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$sample_or_control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$sample_or_control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

### Remove this contaminant ASVs
s16_3 <- prune_taxa(!contamdf.prev$contaminant, s16_2)
summarize_phyloseq(s16_3) #13,343,174 reads, 4112 ASVs in 356 samples
sum(sample_sums(s16_2)) - sum(sample_sums(s16_3))
# 388822 read identify as contaminants of 48 ASVs

################### Extracting only contaminants ##############################
contaminantes <- prune_taxa(contamdf.prev$contaminant, s16_2) 
contaminantes2 = prune_samples(sample_sums(contaminantes)>=0, contaminantes)
View(tax_table(contaminantes2))
# Plotting total abundance to phylum level of contaminants
plot_bar(contaminantes2, fill = "Phylum")+ 
  geom_bar(aes(color = Phylum), stat = "identity", position = "stack")
###############################################################################

### A useful next step is to explore feature prevalence in the dataset,
# which we will define here as the number of samples in which a taxa appears
# at least once.

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(s16_3),
               MARGIN = ifelse(taxa_are_rows(s16_3), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(s16_3),
                    tax_table(s16_3))
#View(prevdf)
#write.csv(prevdf, file = "prevalencia.csv")

plyr::ddply(prevdf, "Phylum", function(df1)
{cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Filtering phyla whose prevalence is 0
filterPhyla = c("Campylobacterota", "Crenarchaeota", "Euryarchaeota",
                "Fibrobacterota", "Patescibacteria", "Synergistota")

# Filtering phyla whose prevalence is 1
filterPhyla = c("Euryarchaeota", "FCPU426", "RCP2-54",
                "Synergistota", "Thermotogota", "WPS-2")
# Filter entries with low prevalence Phylum.
s16_4 = subset_taxa(s16_3, !Phylum %in% filterPhyla)
s16_4 # Final ASV data set: 4106 ASVs in 227 samples


# Filtering samples of Morelia site
s16_5 = subset_samples(s16_4, sitio=="morelia")
summarize_phyloseq(s16_5)  #9,453,658 reads, 4106 ASVs
# Removing Chloroplast and Mitochondria sequences and na samples in condition
s16_5 = subset_samples(s16_5, condicion !="BLANCO")
s16_5 = subset_taxa(s16_5, Order != "Chloroplast")
s16_5 = subset_taxa(s16_5, Family != "Mitochondria")
# Hay algunas muestras que quedaron sin secuencias, ya que todas sus secuencias
# pertenecian al Reino Viridiplantae. Por lo que voy a excluir esas 
# muestras que quedaron con 0 secuencias:
s16_5 = prune_samples(sample_sums(s16_5)>1, s16_5)
s16_5
summarize_phyloseq(s16_5) #146,201 reads, 3089 ASVs

#Remove taxa with 0 observation
s16_6 <- prune_taxa(taxa_sums(s16_5) > 0, s16_5)
s16_6
summarize_phyloseq(s16_6) #146,201 reads, 1635 ASVs in 227 samples

#View(otu_table(s16_6))
#write.csv(otu_table(s16_6), file = "s16_6.csv")

# Primero resumimos el contenido del objeto phyloseq
summarize_phyloseq(s16_6) # reads totales 146,201
sort(sample_sums(s16_6))

# información del número de reads
summary(sample_sums(s16_6))
range(sample_sums(s16_6))
#Library size spans 2 - 24079 reads 

# Obteniendo el número de read por muestra
reads_sample = readcount(s16_6)
reads_sample[1:5]
#write.csv(reads_sample, file = "reads_sample_16s_todas.csv")

#s16_7 = subset_samples(s16_6, sampleid2 != "281022-EMo22-SA084")
#View(otu_table(s16_7))
#write.csv(otu_table(s16_7), file = "outlier.csv")

# Plotting total abundance to phylum level
plot_bar(s16_6, fill = "Phylum")+ 
  geom_bar(aes(color = Phylum), stat = "identity", position = "stack")

# Plotting relative abundance to phylum level 
s16_6_re = transform_sample_counts(s16_6, function(x) x / sum(x))
plot_bar(s16_6_re, fill = "Phylum")+ 
  geom_bar(aes(color = Phylum), stat = "identity", position = "stack")

# checking correlation between reads and ASVs observados
ggplot(data = data.frame("total_reads" =  phyloseq::sample_sums(s16_6),
                         "observed" = phyloseq::estimate_richness(s16_6,
                                      measures = "Observed")[, 1]),
       aes(x = total_reads, y = observed)) +
  geom_point() +
  geom_smooth(method="lm", se = FALSE) +
  labs(x = "\nTotal Reads", y = "Observed Richness\n")


###############################################################################
#### ANALYSIS WITH PROPORTION DATA ####
### Calculating relative abundances
mergedGP_pro <- merge_samples(s16_6, "especies")
# Use psmelt to obtain a long-format data.frame
phy_pro <- mergedGP_pro %>% tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
#View(phy_pro)
#write.csv(phy_pro, file = "abundancias_family_propor.csv")

#mergedtop = merge_samples(topfam, "especies")
#phy2 <- mergedtop %>% tax_glom(taxrank = "Family") %>%
#  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
#View(phy2)

#### HEATMAPs of most abundant bacteria ####
### Collapse dataset taxonomy to 10 most abundant families
topfam_pro <- aggregate_top_taxa2(s16_6, "Family", top = 10)
#View(otu_table(topfam))
topfam_pro

# transform families dataset to relative abundance
transform <- microbiome::transform
famrel_pro <- transform(topfam_pro, "compositional")
#View(otu_table(famrel))

### Heatmap for 10 most abundant bacterial families between species
heatmap_k_pro <- simple_heatmap(famrel_pro,
                                group.facet = "especies",
                                group.order = NULL,
                                abund.thres = 0.01,
                                prev.thres = 0.01,
                                level = "Family",
                                scale.color = "log10",
                                na.fill = "white",
                                color.fill = NULL,
                                taxa.arrange=TRUE,
                                remove.other=TRUE,
                                panel.arrange="grid",
                                ncol=NULL,
                                nrow=NULL)
heatmap_k_pro


# Separando por sexos florales
hembras_pro <- subset_samples(s16_6, sexo == "hembra")
machos_pro <- subset_samples(s16_6, sexo == "macho")

### Calculating relative abundances for Female flowers
merged_hembras_pro <- merge_samples(hembras_pro, "especies")
# Use psmelt to obtain a long-format data.frame
phy_hembras_pro <- merged_hembras_pro %>% tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
#View(phy_hembras_pro)
#write.csv(phy_hembras_pro, file = "abundancias_hembras_pro.csv")

### Calculating relative abundances for Male flowers
merged_machos_pro <- merge_samples(machos_pro, "especies")
# Use psmelt to obtain a long-format data.frame
phy_machos_pro <- merged_machos_pro %>% tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
#View(phy_machos_pro)
#write.csv(phy_machos_pro, file = "abundancias_machos_pro.csv")


# Flores hembra
topfam_hempro <- aggregate_top_taxa2(hembras_pro, "Family", top = 10)
# transform families dataset to relative abundance
transform <- microbiome::transform
famrel_hempro <- transform(topfam_hempro, "compositional")
# Heatmap de flores hembra
heatmap_hempro <- simple_heatmap(famrel_hempro,
                                 group.facet = "especies",
                                 group.order = NULL,
                                 abund.thres = 0.01,
                                 prev.thres = 0.01,
                                 level = "Family",
                                 scale.color = "log10",
                                 na.fill = "white",
                                 color.fill = NULL,
                                 taxa.arrange=TRUE,
                                 remove.other=TRUE,
                                 panel.arrange="grid",
                                 ncol=NULL,
                                 nrow=NULL)
heatmap_hempro

# Flores machos
topfam_macpro <- aggregate_top_taxa2(machos_pro, "Family", top = 10)
# transform families dataset to relative abundance
transform <- microbiome::transform
famrel_macpro <- transform(topfam_macpro, "compositional")
# Heatmap de flores hembra
heatmap_macpro <- simple_heatmap(famrel_macpro,
                                 group.facet = "especies",
                                 group.order = NULL,
                                 abund.thres = 0.01,
                                 prev.thres = 0.01,
                                 level = "Family",
                                 scale.color = "log10",
                                 na.fill = "white",
                                 color.fill = NULL,
                                 taxa.arrange=TRUE,
                                 remove.other=TRUE,
                                 panel.arrange="grid",
                                 ncol=NULL,
                                 nrow=NULL)
heatmap_macpro

### Arrange heatmaps in one graphic
library(ggpubr)
library(gridExtra)
heatmaps <- grid.arrange(heatmap_hempro, heatmap_macpro, nrow = 2)

# Most dominant families en todas las muestras
domipro <- dominant_taxa(s16_6,level = "Family") #group = "especies"
#View(domipro$dominant_overview)
write.csv(domipro$dominant_overview, file = "dominant_overview_pro.csv")


#### CORE MICROBIOTA ####

# PREVALENCE OF TAXONOMIC GROUPS
ps.relpro = microbiome::transform(s16_6, "compositional")
# Relative population frequencies; at 1% compositional abundance threshold:
head(prevalence(ps.relpro, detection = 10/100, sort = TRUE))
# Absolute population frequencies (sample count)
head(prevalence(ps.relpro, detection = 1/100, sort = TRUE, count = TRUE))

# Core tresholds: minimum ASV relative abundance 0.001 % minimum prevalence
# among samples 10%
core.taxapro = core_members(ps.relpro, detection = .001, prevalence = .1,
                            include.lowest = T)
print(core.taxapro)
# A full phyloseq object of the core microbiota is obtained as follows:
full_corepro <- core(ps.relpro, detection = 0.001, prevalence = .1)
View(tax_table(full_corepro))
# Retrieving the core taxa names from the phyloseq object:
core.taxa2 = taxa(full_corepro)

# CORE ABUNDANCE AND DIVERSITY
# total core abundance in each sample (sum of abundance of the core members):
core.abundance <- sample_sums(core(ps.relpro, detection = .001, prevalence = .1))


prevdf = apply(X = otu_table(s16_6),
               MARGIN = ifelse(taxa_are_rows(s16_6), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(s16_6),
                    tax_table(s16_6))
# plotting prevalence
table(prevdf$Prevalence)
ggplot(prevdf, aes(Prevalence)) + geom_histogram(bins = 50) + theme_bw() + 
  ggtitle("Histogram of Taxa Prevalence")

View(head(otu_table(s16_6)))
head(prevdf)

# CORE VISUALIZATION
# With compositional (relative) abundances
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)
#ggplot(d) + geom_point(aes(x, y)) + scale_x_continuous(trans="log10", limits=c(NA,1))

plot_core(ps.rel, 
          prevalences = prevalences, 
          detections = det, 
          plot.type = "lineplot") + 
  xlab("Relative Abundance (%)")

### CORE HEATMAPS
library(RColorBrewer)
library(reshape)

prevalences <- seq(.05, 1, .05)

detections <- round(10^seq(log10(0.01), log10(.2), length = 9), 3)

# Also define gray color palette
gray <- gray(seq(0,1,length=5))

#Added pseq.rel, I thin... must be checked if it was in the the rednred version,; where it is initialized
#pseq.rel<- microbiome::transform(pseq, 'compositional')
#min-prevalence gets the 100th highest prevalence
p <- plot_core(ps.rel,
               plot.type = "heatmap", 
               colours = gray,
               prevalences = prevalences, 
               detections = detections, 
               min.prevalence = prevalence(ps.rel, sort = TRUE)[100]) +
  labs(x = "Detection Threshold\n(Relative Abundance (%))") +
  
  #Adjusts axis text size and legend bar height
  theme(axis.text.y= element_text(size=8, face="italic"),
        axis.text.x.bottom=element_text(size=8),
        axis.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))

print(p)

# Plotting with upsetR
library(MicrobiotaProcess)
library(UpSetR)


uppenpro <- get_upset(s16_6, factorNames = "especies")
View(uppenpro)
# función para gráficar los ASVs compartidos entre especies
?upset

upset(uppenpro, nsets = 6, mainbar.y.label = "Number of ASVs shared",
      order.by = "freq",
      queries = list(list(query = intersects,
                          params = list("6CM", "5CAA", "4CPP"),
                          color = "orange",
                          active = T),
                     list(query = intersects,
                          params = list("6CM", "5CAA", "4CPP", 
                                        "3CAS", "2CPF", "1CF"),
                          color = "green",
                          active = T)))

#list(query = intersects,
#    params = list("3CAS", "2CPF", "1CF"),
#   color = "pink",
#  active = T)))
#text.scale = 1.1, 
#point.size = 2.8, 
#line.size = 1




###############################################################################
#### ALPHA DIVERSITY WITH PROPORTIONS DATA  ####
#### ANALYSIS WITH PROPORTIONS DATA ####
# Primero nos aseguramos que estamos usando las funciones del paquete correcto
transform <- microbiome::transform
# Proportion Normalization
pseq_16s <- transform(s16_6, "compositional")


library(knitr)

# Obteniendo los valores de las metricas de diversidad alfa
div2 = microbiome::alpha(pseq_16s, index = "all")
kable(head(div2))
View(div2)

# get the metadata out as separate object
ps2.meta <- meta(pseq_16s)
kable(head(ps2.meta))
# Add the rownames as a new colum for easy integration later.
ps2.meta$sam_name2 <- rownames(ps2.meta)

# Add the rownames to diversity table
div2$sam_name2 <- rownames(div2)

# merge these two data frames into one
div.df2 <- merge(div2,ps2.meta, by = "sam_name2")

# check the tables
colnames(div.df2)

### Phylogenetic diversity (Faith)
library(picante)

ps0.rar.asvtab2 <- as.data.frame(pseq_16s@otu_table)

ps0.rar.tree2 <- pseq_16s@phy_tree

# hmp.meta from previous code chunks

# We first need to check if the tree is rooted or not 
pseq_16s@phy_tree

# it is a rooted tree
df.pd2 <- pd(t(ps0.rar.asvtab2), ps0.rar.tree2,include.root=T) 
# t(ou_table) transposes the table for use in picante and the tre file comes
# from the first code chunck we used to read tree file (see making a phyloseq 
# object section).
#datatable(df.pd2)

# We add the results of PD to this file and then plot
div.df2$Phylogenetic_Diversity <- df.pd2$PD
head(div.df2)


#### Diversidad alpha en flores hembra ####
head(div.df2)
hembras_alfa = subset(div.df2, sexo == "hembra")
str(hembras_alfa)

# data
sort(hembras_alfa$observed)
hembras_alfa$diversity_shannon
hembras_alfa$Phylogenetic_Diversity
# Histograma
hist(hembras_alfa$observed)
hist(hembras_alfa$diversity_shannon) 
hist(hembras_alfa$Phylogenetic_Diversity)
# Shapiro test
shapiro.test(hembras_alfa$observed)
shapiro.test(hembras_alfa$diversity_shannon)
shapiro.test(hembras_alfa$Phylogenetic_Diversity)

#### STATISTICAL ANALYSIS OF ALPHA DIVERSITY ####
# Analysis with kruskal-Wallis and Wilcox text in female flowers #
# Observed by species
kruskal.test(hembras_alfa$observed~hembras_alfa$especies)
pairwise.wilcox.test(hembras_alfa$observed, hembras_alfa$especies, exact=F)

# shannon
kruskal.test(hembras_alfa$diversity_shannon~hembras_alfa$especies)
pairwise.wilcox.test(hembras_alfa$diversity_shannon, hembras_alfa$especies, exact=F)
# PD
kruskal.test(hembras_alfa$Phylogenetic_Diversity~hembras_alfa$especies)
pairwise.wilcox.test(hembras_alfa$Phylogenetic_Diversity, hembras_alfa$especies, exact=F)


#### Diversidad alpha en flores macho ####
head(div.df2)
machos_alfa = subset(div.df2, sexo == "macho")
str(machos_alfa)

# data
sort(machos_alfa$observed)
machos_alfa$diversity_shannon
machos_alfa$Phylogenetic_Diversity
# Histograma
hist(machos_alfa$observed)
hist(machos_alfa$diversity_shannon) 
hist(machos_alfa$Phylogenetic_Diversity)
# Shapiro test
shapiro.test(machos_alfa$observed)
shapiro.test(machos_alfa$diversity_shannon)
shapiro.test(machos_alfa$Phylogenetic_Diversity)

# Analysis with kruskal-Wallis and Wilcox text in male flowers #
# Observed
kruskal.test(machos_alfa$observed~machos_alfa$especies)
pairwise.wilcox.test(machos_alfa$observed, machos_alfa$especies, exact=F)
# shannon
kruskal.test(machos_alfa$diversity_shannon~machos_alfa$especies)
pairwise.wilcox.test(machos_alfa$diversity_shannon, machos_alfa$especies, exact=F)
# PD
kruskal.test(machos_alfa$Phylogenetic_Diversity~machos_alfa$especies)
pairwise.wilcox.test(machos_alfa$Phylogenetic_Diversity, machos_alfa$especies,
                     exact=F)


### STATISTICAL ANALYSIS OF ALPHA DIVERSITY
library(emmeans)
mean(div.df2$observed)
var(div.df2$observed)
# Richness
mod4 = glm(observed ~ especies*sexo, data = div.df2,
           family = quasipoisson)
summary(mod4)
anova(mod4)
# post-hoc comparisons
compa4 = emmeans(mod4, 
                 specs = pairwise ~ especies | sexo, 
                 type = "response")
compa4
compa4.1 = emmeans(mod4,
                   specs = pairwise ~ especies,
                   type = "response")
compa4.1
#write.csv(compa1$contrasts, file = "observed_compa_rarefied.csv")

# Shannon
mod5 = glm(diversity_shannon ~ especies*sexo, data = div.df2)
summary(mod5)
anova(mod5)
# post-hoc comparisons
compa5 = emmeans(mod5, 
                 specs = pairwise ~ especies | sexo, 
                 type = "response")
compa5
#modx = aov(diversity_shannon ~ especies*sexo*tratamiento, data = div.df)
#summary(modx)

# Phylogenetic Diversity index
mod6 = glm(Phylogenetic_Diversity ~ especies*sexo, data = div.df2)
summary(mod6)
anova(mod6)
# post-hoc comparisons
compa6 = emmeans(mod6, 
                 specs = pairwise ~ especies | sexo, 
                 type = "response")
compa6

##### VISUALIZATION  OF ALPHA DIVERSITY WITH PROPORTION DATA####
# Stripchart of observed ASVs of female flowers
obs_hem <- ggplot(hembras_alfa, aes(x=especies, y=observed, color=especies)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 2)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Richness")+
  scale_color_discrete(name = "Species", 
                       labels = c("2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  scale_x_discrete(labels=c("2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  #theme_bw()
  theme_classic()
#Exportamos grafico de alfa div (observed) a pdf
ggsave("observed_hembras.pdf", 
       plot = obs_hem,
       width = 25, 
       height = 20, 
       units = "cm")

# Stripchart of observed ASVs of male flowers
obs_machos <- ggplot(machos_alfa, aes(x=especies, y=observed, color=especies)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 2)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Richness")+
  scale_color_discrete(name = "Species", 
                       labels = c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  scale_x_discrete(labels=c("1CF" = "CF", "2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  #theme_bw()
  theme_classic()
#Exportamos grafico de alfa div (observed) a pdf
ggsave("observed_machos.pdf", 
       plot = obs_machos,
       width = 25, 
       height = 20, 
       units = "cm")

# combinando plots
# cargando paqueteria gridExtra
library(gridExtra)
library(cowplot)
obs_pro <- grid.arrange(obs_hem, obs_machos, 
                        ncol = 2, nrow = 1)

obs_pro <- ggarrange(obs_hem, obs_machos, ncol = 2, nrow = 1,
                     labels = "AUTO") 
# exporting plots
ggsave("alpha_obser_nov06.pdf", 
       plot = obs_pro,
       width = 25, 
       height = 20, 
       units = "cm")

### stripchart de Shannon diversity index in female flowers
sha_hem = ggplot(hembras_alfa, aes(x=especies, y=diversity_shannon,
                                   color=especies)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 1)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Shannon diversity index")+
  scale_color_discrete(name = "Species", 
                       labels = c("2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  scale_x_discrete(labels=c("2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  theme_classic()
sha_hem    


### stripchart de Shannon diversity index in male flowers
sha_mac = ggplot(machos_alfa, aes(x=especies, y=diversity_shannon,
                                  color=especies)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 1)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Shannon diversity index")+
  scale_color_discrete(name = "Species", 
                       labels = c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  scale_x_discrete(labels=c("1CF" = "CF", "2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  theme_classic()
sha_mac

#facet_wrap(~sexo, labeller = labeller(sexo = sexos))
#stat_summary(fun=mean, geom="point", shape=18,
#            size=6, color="red")
sha_pro <- ggarrange(sha_hem, sha_mac, ncol = 2, nrow = 1,
                     labels = "AUTO") 
# exporting plots
ggsave("alpha_shan_nov06.pdf", 
       plot = sha_pro,
       width = 25, 
       height = 20, 
       units = "cm")


### stripchart de phylogenetic diversity index in female flowers
pd_hem = ggplot(hembras_alfa, aes(x=especies, y=Phylogenetic_Diversity,
                                  color=especies)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 1)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Shannon diversity index")+
  scale_color_discrete(name = "Species", 
                       labels = c("2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  scale_x_discrete(labels=c("2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  theme_classic()
pd_hem    


### stripchart de Shannon diversity index in male flowers
pd_mac = ggplot(machos_alfa, aes(x=especies, y=Phylogenetic_Diversity,
                                 color=especies)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 1)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Shannon diversity index")+
  scale_color_discrete(name = "Species", 
                       labels = c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  scale_x_discrete(labels=c("1CF" = "CF", "2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  theme_classic()
pd_mac

#facet_wrap(~sexo, labeller = labeller(sexo = sexos))
#stat_summary(fun=mean, geom="point", shape=18,
#            size=6, color="red")
pd_pro <- ggarrange(pd_hem, pd_mac, ncol = 2, nrow = 1,
                    labels = "AUTO") 
# exporting plots
ggsave("PD_shan_nov06.pdf", 
       plot = pd_pro,
       width = 25, 
       height = 20, 
       units = "cm")

### stripchart Shannon diversity index facet with 4 panels
pro16_2 = ggplot(div.df2, aes(x=especies, y=diversity_shannon, color=especies))+ 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 2)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Shannon diversity index",
       title = "")+
  scale_color_discrete(name = "Species", 
                       labels = c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  scale_x_discrete(labels=c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  theme_bw()+
  facet_wrap(~sexo, labeller = labeller(sexo = sexos))
#stat_summary(fun=mean, geom="point", shape=18,
#            size=6, color="red")
# shannon_16s_1 = shannon_16s_1 + 
#stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = 3)

### stripchart phylogenetic diversity index facet with 2 panels
pro16_3 = ggplot(div.df2, aes(x=especies, y=Phylogenetic_Diversity, color=especies))+ 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 2)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Phylogenetic Diversity index",
       title = "")+
  scale_color_discrete(name = "Species", 
                       labels = c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  scale_x_discrete(labels=c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  theme_bw()+
  facet_wrap(~sexo, labeller = labeller(sexo = sexos))
#stat_summary(fun=mean, geom="point", shape=18,
#            size=6, color="red")
# shannon_16s_1 = shannon_16s_1 + 
#stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = 3)

# cargando paqueteria gridExtra
library(gridExtra)
library(cowplot)
gt_pro <- grid.arrange(pro16_1, pro16_2, pro16_3, 
                       ncol = 1, nrow = 3)

### gráfica del indice de Shannon entre especies y sexo
pmult <- ggboxplot(hembras_alfa, x = "especies", y = "observed", alpha = .5,
                   fill = "especies", add = "jitter",
                   short.panel.labs = FALSE)
pmult + stat_compare_means(aes(label = paste0("p = ", ..p.format..)))

labs(x = "Especies", y = "Índice de diversidad de Shannon")+
  scale_fill_discrete(name = "Especies") +
  scale_fill_discrete(name = "Especies", 
                      labels = c("CF","CPF","CPP","CAS", "CAA", "CM"))+
  scale_x_discrete(labels=c("1CF" = "CF", "2CPF" = "CPF",
                            "3CPP" = "CPP", "4CAS" = "CAS",
                            "5CAA" = "CAA", "6CM" = "CM"))

pmult + stat_compare_means()
pmult + stat_compare_means(aes(label = paste0("p = ", ..p.format..)))


# create a list of pairwise comparisons
spp <- base::levels(as.factor(hembras_alfa$especie)) # get the variables
# make a pairwise list that we want to compare.
spp.pairs <- combn(seq_along(spp), 2, simplify = FALSE, FUN = function(i)spp[i])
print(spp.pairs)
pmult <- pmult + stat_compare_means(comparisons = spp.pairs) 

p2 <- p2 + stat_compare_means(comparisons = spp.pairs) 
print(p2) # Wilcox test


########## BETA DIVERSITY  WITH PROPORTION DATA ##########
library(vegan)
# phyloseq object filtered but no rarefied
s16_6 
# divided by sex
hembras_pro
machos_pro

#### Unweighted UniFrac HEMBRAS ####
# Primero nos aseguramos que estamos usando las funciones del paquete correcto
transform <- microbiome::transform
# Proportion Normalization
pseq_16s_hem <- transform(hembras_pro, "compositional")


# PCoA plot using the unweighted UniFrac as distance
ordination_16s_hem = ordinate(pseq_16s_hem, method="PCoA", "unifrac",
                              weighted = F)

baseplot_16s_hem <- plot_ordination(pseq_16s_hem, ordination_16s_hem,
                                    color="especies",
                                    type = "samples",
) + #shape = "sexo"
  ggtitle("")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 2)
#stat_ellipse()

baseplot_16s_hem

#baseplot_16s_hem +
#  geom_point(size = 1)
#  geom_label(aes(label = clave), size = 3.5, color = "black")

#### Unweighted UniFrac MACHOS ####
# Primero nos aseguramos que estamos usando las funciones del paquete correcto
transform <- microbiome::transform
# Proportion Normalization
pseq_16s_mac <- transform(machos_pro, "compositional")

# PCoA plot using the unweighted UniFrac as distance
ordination_16s_mac = ordinate(pseq_16s_mac, method="PCoA", "unifrac",
                              weighted = F)

baseplot_16s_mac <- plot_ordination(pseq_16s_mac, ordination_16s_mac,
                                    color="especies",
                                    type = "samples",
) + #shape = "sexo"
  ggtitle("")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 2)
#stat_ellipse()
baseplot_16s_mac


#### Weighted UniFrac HEMBRAS####
# PCoA plot using the weighted UniFrac as distance
ordination_16s_hem_wei = ordinate(pseq_16s_hem,
                                  method="PCoA", "wunifrac")
# Plot
baseplot_16s_hem_wei <- plot_ordination(pseq_16s_hem, 
                                        ordination_16s_hem_wei,
                                        color="especies",
                                        type = "samples") + #shape = "sexo" shape = "tratamiento"
  ggtitle("")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 2)
#labs(x = "PCoa1 (31.3%)", y = "PCoA2 (10.4%)")+
#scale_color_discrete(name = "Condition", 
#                    labels = c("Domesticated","Wild"))+
#scale_color_discrete(name = "Treatment",
#                     labels = c("No-visited","Visited"))

baseplot_16s_hem_wei

#baseplot_16s_hem_wei + stat_ellipse()
#
#baseplot2 +
#  geom_point(size = 1) +
#  geom_label(aes(label = clave), size = 3.5, color = "black")
#

#### Weighted UniFrac MACHOS####
# PCoA plot using the weighted UniFrac as distance
ordination_16s_mac_wei = ordinate(pseq_16s_mac,
                                  method="PCoA", "wunifrac")
# Plot
baseplot_16s_mac_wei <- plot_ordination(pseq_16s_mac, 
                                        ordination_16s_mac_wei,
                                        color="especies",
                                        type = "samples") + #shape = "sexo" shape = "tratamiento"
  ggtitle("")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 2)
#labs(x = "PCoa1 (31.3%)", y = "PCoA2 (10.4%)")+
#scale_color_discrete(name = "Condition", 
#                    labels = c("Domesticated","Wild"))+
#scale_color_discrete(name = "Treatment",
#                     labels = c("No-visited","Visited"))

baseplot_16s_mac_wei

#baseplot_16s_mac_wei + stat_ellipse()
#
#baseplot2 +
#  geom_point(size = 1) +
#  geom_label(aes(label = clave), size = 3.5, color = "black")


#### PERMANOVA by floral gender####
## female flowers
metadf_hem = data.frame(sample_data(hembras_pro))

# calculating Weighted-UniFrac
wunifrac.dist_hem <- UniFrac(pseq_16s_hem, 
                             weighted = TRUE, 
                             normalized = TRUE,  
                             parallel = FALSE, 
                             fast = TRUE)

#Test dispersion among body sites
disp_wun <- betadisper(wunifrac.dist_hem, sample_data(hembras_pro)$especies)
permutest(disp_wun)

# Run PERMANOVA
set.seed(100)   # Set the seed to make the analysis reproducible
perm_hem = adonis2(wunifrac.dist_hem ~ especies, data = metadf_hem,
                   by = "terms", permutations = 999)
perm_hem

## Male flowers
metadf_mac = data.frame(sample_data(machos_pro))

# calculating Weighted-UniFrac
wunifrac.dist_mac <- UniFrac(pseq_16s_mac, 
                             weighted = TRUE, 
                             normalized = TRUE,  
                             parallel = FALSE, 
                             fast = TRUE)

#Test dispersion among body sites
disp_wun_mac <- betadisper(wunifrac.dist_mac, sample_data(machos_pro)$especies)
permutest(disp_wun_mac)

# Run PERMANOVA
set.seed(100)   # Set the seed to make the analysis reproducible
perm_mac = adonis2(wunifrac.dist_mac ~ especies, data = metadf_mac,
                   by = "terms", permutations = 999)
perm_mac


###############################################################################
# calculating Unweighted-UniFrac
unifrac.dist <- UniFrac(pseq_16s, 
                        weighted = FALSE, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)
# Run PERMANOVA
set.seed(100)   # Set the seed to make the analysis reproducible
perm4 = adonis2(unifrac.dist ~ especies*sexo, data = metadf,
                by = "terms", permutations = 999)
perm4  

?betadisper
ps.disper <- betadisper(unifrac.dist, metadf$especies)
permutest(ps.disper, pairwise = TRUE)
anova(ps.disper)
TukeyHSD(ps.disper)










####################### METODOS DE NORMALIZACION ##############################
                         #### RAREFIED DATA ####
otu_tab = t(abundances(s16_6))
p = vegan::rarecurve(otu_tab,
                     step = 50, label = FALSE,
                     sample = min(rowSums(otu_tab),
                                  col = "blue, cex = 0.6"))

#vegan::rarecurve(t(otu_table(s16_6)), step=50, cex=0.5)

### Rarefaction plots
#Set seed
set.seed(1)
subsamples <- seq(0, 20000, by=1000)[-1]

#Repeat for each data set, here exemplified with king penguin clean data set
sample_data(s16_6)
p <- plot_alpha_rcurve(s16_6, index="observed", 
                       #subsamples = subsamples,
                       lower.conf = 0.025, 
                       upper.conf = 0.975,
                       group="especies",
                       label.color = "brown3",
                       label.size = 3,
                       label.min = TRUE) 
p

### Normalize library size to 1000 reads
set.seed(111)
ps.rarefied1000 <- rarefy_even_depth(s16_6, 
                                    1000, rngseed = 1,
                                    replace = FALSE)             
ps.rarefied1000
# phyloseq object with normalized library size has
# 754 ASVs in 39 samples
# 688 ASVs in 36 samples 20-03-2025

### Normalize library size to 100 reads
set.seed(111)
ps.rarefied100 <- rarefy_even_depth(s16_6, 
                                   100, rngseed = 1,
                                   replace = FALSE)             
ps.rarefied100
# phyloseq object with normalized library size has
# 906 ASVs in 134 samples
# 506 ASVs in 132 samples
# 898 ASVs in 132 samples 20/03/2025
# 901 ASVs in 132 samples 21/03/2025
reads_sample2 = readcount(ps.rarefied100)
reads_sample2[1:20]
summarize_phyloseq(ps.rarefied100) # 13200 reads
summary(sample_sums(ps.rarefied100))

### Normalize library size to 50 reads
set.seed(111)
ps.rarefied50 <- rarefy_even_depth(s16_6, 
                             sample.size =50, rngseed = 1,
                             replace = FALSE)             
ps.rarefied50 # use this rarefied depth
# phyloseq object with normalized library size has
# 785 ASVs in 168 samples
# 757 asvs in 168 samples 26/02/2025
# 784 ASVs in 168 samples 03/03/2025
# 771 ASVs in 167 samples 20/03/2025
reads_sample2 = readcount(ps.rarefied50)
reads_sample2[1:20]
summarize_phyloseq(ps.rarefied50)
# total number of reads= 8400 in rarefied data
summary(sample_sums(ps.rarefied50))
abundances(ps.rarefied50)

### Normalize library size to 10 reads
set.seed(111)
ps.rarefied10 <- rarefy_even_depth(s16_6, 
                                   10, rngseed = 1,
                                   replace = FALSE)             
ps.rarefied10
# phyloseq object with normalized library size has
# 392 ASVs in 220 samples

# sampling depth
barplot(sample_sums(ps.rarefied100), las = 2)

# taxa prevalence
plot_taxa_prevalence(s16_6, "Phylum")
plot_taxa_prevalence(ps.rarefied100, "Phylum")

                #### ALPHA DIVERSITY WITH RAREFIED DATA  ####
library(knitr)

div1 = microbiome::alpha(ps.rarefied100, index = "all")
kable(head(div1))
View(div1)

# get the metadata out as separate object
ps1.meta <- meta(ps.rarefied100)
kable(head(ps1.meta))
# Add the rownames as a new colum for easy integration later.
ps1.meta$sam_name <- rownames(ps1.meta)

# Add the rownames to diversity table
div1$sam_name <- rownames(div1)

# merge these two data frames into one
div.df <- merge(div1,ps1.meta, by = "sam_name")

# check the tables
colnames(div.df)

### Phylogenetic diversity (Faith)
library(picante)

ps0.rar.asvtab <- as.data.frame(ps.rarefied100@otu_table)

ps0.rar.tree <- ps.rarefied100@phy_tree

# hmp.meta from previous code chunks

# We first need to check if the tree is rooted or not 
ps.rarefied100@phy_tree

# it is a rooted tree
df.pd <- pd(t(ps0.rar.asvtab), ps0.rar.tree,include.root=T) 
# t(ou_table) transposes the table for use in picante and the tre file comes
# from the first code chunck we used to read tree file (see making a phyloseq 
# object section).
#datatable(df.pd)

# We add the results of PD to this file and then plot
div.df$Phylogenetic_Diversity <- df.pd$PD
head(div.df)


# data
sort(div.df$observed)
div.df$diversity_shannon
div.df$Phylogenetic_Diversity
# Histograma
hist(div.df$observed)
hist(div.df$diversity_shannon) 
hist(div.df$Phylogenetic_Diversity)
# Shapiro test
shapiro.test(div.df$observed)
shapiro.test(div.df$diversity_shannon)
shapiro.test(div.df$Phylogenetic_Diversity)


### STATISTICAL ANALYSIS OF ALPHA DIVERSITY
library(emmeans)
# Richness
mean(div.df$observed)
var(div.df$observed)

mod1 = glm(observed ~ especies*sexo, data = div.df, family = quasipoisson)
summary(mod1)
anova(mod1)
shapiro.test(residuals(mod1))
plot(residuals(mod1))

# post-hoc comparisons
compa1 = emmeans(mod1, 
        specs = pairwise ~ especies, 
        type = "response")
compa1
compa1.1 = emmeans(mod1,
                   specs = pairwise ~ especies | sexo,
                   type = "response")
compa1.1
#write.csv(compa1$contrasts, file = "observed_compa_rarefied.csv")

# Shannon
mod2 = glm(diversity_shannon ~ especies*sexo, data = div.df)
summary(mod2)
anova(mod2)
# post-hoc comparisons
compa2 = emmeans(mod2, 
                   specs = pairwise ~ especies, 
                   type = "response")
compa2
compa2.1 = emmeans(mod2, 
                 specs = pairwise ~ especies | sexo, 
                 type = "response")
compa2.1
#modx = aov(diversity_shannon ~ especies*sexo*tratamiento, data = div.df)
#summary(modx)

# Phylogenetic Diversity index
mod3 = glm(Phylogenetic_Diversity ~ especies*sexo, data = div.df)
summary(mod3)
anova(mod3)
# post-hoc comparisons
compa3 = emmeans(mod3, 
                 specs = pairwise ~ especies, 
                 type = "response")
compa3
compa3.1 = emmeans(mod3, 
                 specs = pairwise ~ especies | sexo, 
                 type = "response")
compa3.1
######## Analysis with kruskal-Wallis and Wilcox text
# Observed by species
kruskal.test(div.df$observed~div.df$especies)
pairwise.wilcox.test(div.df$observed, div.df$especies, exact=F)
# shannon
kruskal.test(div.df$diversity_shannon~div.df$especies)
pairwise.wilcox.test(div.df$diversity_shannon, div.df$especies, exact=F)
# PD
kruskal.test(div.df$Phylogenetic_Diversity~div.df$especies)
pairwise.wilcox.test(div.df$Phylogenetic_Diversity, div.df$especies, exact=F)

# by sex
# Observed
kruskal.test(div.df$observed~div.df$sexo)
pairwise.wilcox.test(div.df$observed, div.df$sexo, exact=F)
# shannon
kruskal.test(div.df$diversity_shannon~div.df$sexo)
pairwise.wilcox.test(div.df$diversity_shannon, div.df$sexo, exact=F)
# PD
kruskal.test(div.df$Phylogenetic_Diversity~div.df$sexo)
pairwise.wilcox.test(div.df$Phylogenetic_Diversity, div.df$sexo, exact=F)



##### VISUALIZATION  OF ALPHA DIVERSITY with rarefied data###
# Stripchart of observed ASVs with 1 panel
ggplot(div.df, aes(x=especies, y=observed, color=especies)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 2)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Especies", y = "Riqueza de ASVs")+
  scale_color_discrete(name = "Especies", 
                       labels = c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  scale_x_discrete(labels=c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  theme_bw()

### stripchart de Observed ASVs facet with 2 panels (divided by sex)
sexos <- c("Hembras", "Machos")
names(sexos) <- c("hembra", "macho")

riqueza_16 <- ggplot(div.df, aes(x=especies, y=observed, color=especies)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 1)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Especies", y = "Riqueza de ASVs")+
  scale_color_discrete(name = "Especies", 
                       labels = c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  scale_x_discrete(labels=c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  theme_bw()+
  facet_wrap(~sexo, labeller = labeller(sexo = sexos))
#stat_summary(fun=mean, geom="point", shape=18,
#            size=6, color="red")
riqueza_16
### stripchart Shannon diversity index facet with 2 panels (divided by sex)
shannon_16 <- ggplot(div.df, aes(x=especies, y=diversity_shannon, color=especies))+ 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 2)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Especies", y = "Indice de diversidad de Shannon",
       title = "")+
  scale_color_discrete(name = "Especies", 
                       labels = c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  scale_x_discrete(labels=c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  theme_bw()+
  facet_wrap(~sexo,
             labeller = labeller(sexo = sexos))
  #stat_summary(fun=mean, geom="point", shape=18,
  #            size=6, color="red")
  # shannon_16s_1 = shannon_16s_1 + 
  #stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = 3)
shannon_16

cartel1 = ggarrange(riqueza_16, shannon_16,
                             nrow = 2, labels = "AUTO")
cartel1
### stripchart phylogenetic diversity index facet with 2 panels (divided by sex)
pd_16 <- ggplot(div.df, aes(x=especies, y=Phylogenetic_Diversity,
                            color=especies))+ 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 2)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Phylogenetic Diversity index",
       title = "")+
  scale_color_discrete(name = "Species", 
                       labels = c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  scale_x_discrete(labels=c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  theme_bw()+
  facet_wrap(~sexo,
             labeller = labeller(sexo = sexos))
#stat_summary(fun=mean, geom="point", shape=18,
#            size=6, color="red")
# shannon_16s_1 = shannon_16s_1 + 
#stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = 3)

# cargando paqueteria gridExtra
library(gridExtra)
library(cowplot)

gt_s16 = ggarrange(riqueza_16, shannon_16, pd_16,
                   nrow = 3, labels = "AUTO")
gt_s16


gt.its = ggarrange(riqueza_its, shannon_its, pd_its,
                     nrow = 3, labels = "AUTO")
gt.its



###### OTRO TIPO DE GRÁFICA ####
### gráfica de riqueza por especie y sexo 1 panel
riqueza_16_group <- ggplot(div.df, aes(x=especies, y=observed, fill=sexo)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 1)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Richness")+
  scale_fill_discrete(name = "Species", 
                       labels = c("hembra" = "Female","macho" = "Male"))+
  scale_x_discrete(labels=c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  theme_bw()

#  facet_wrap(~sexo, labeller = labeller(sexo = sexos))
#stat_summary(fun=mean, geom="point", shape=18,
#            size=6, color="red")

### grafica de shannon por especie y sexo 1 panel
shannon_16_group <- ggplot(div.df, aes(x=especies, y=diversity_shannon,
                                       fill=sexo)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 1)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Richness")+
  scale_fill_discrete(name = "Species", 
                      labels = c("hembra" = "Female","macho" = "Male"))+
  scale_x_discrete(labels=c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  theme_bw()

#  facet_wrap(~sexo, labeller = labeller(sexo = sexos))
#stat_summary(fun=mean, geom="point", shape=18,
#            size=6, color="red")

### grafica de PD por especie y sexo 1 panel
pd_16_group <- ggplot(div.df, aes(x=especies, y=Phylogenetic_Diversity,
                                       fill=sexo)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 1)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Richness")+
  scale_fill_discrete(name = "Species", 
                      labels = c("hembra" = "Female","macho" = "Male"))+
  scale_x_discrete(labels=c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  theme_bw()

#  facet_wrap(~sexo, labeller = labeller(sexo = sexos))
#stat_summary(fun=mean, geom="point", shape=18,
#            size=6, color="red")

gt2 <- grid.arrange(riqueza_16_group, shannon_16_group, pd_16_group, 
                   ncol = 1, nrow = 3)

            ########## BETA DIVERSITY  WITH RAREFIED DATA ##########
library(vegan)
#### Unweighted UniFrac ####
# PCoA plot using the unweighted UniFrac as distance
ordination_16s = ordinate(ps.rarefied100, method="PCoA", "unifrac",
                          weighted = F)

baseplot_16s <- plot_ordination(ps.rarefied100, ordination_16s,
                                color="especies",
                                type = "samples",
) + #shape = "sexo"
  ggtitle("")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 2)
baseplot_16s + stat_ellipse()
#baseplot1 +
#  geom_point(size = 1) +
#  geom_label(aes(label = clave), size = 3.5, color = "black")
#
baseplot_16s +
  geom_point(size = 2)+
  ggtitle("")+
  facet_wrap(~sexo, labeller = labeller(sexo = sexos))+
  stat_ellipse()

#### Weighted UniFrac ####
# Transform data to relative abundance
ps.rarefied100
s16_re = microbiome::transform(ps.rarefied100, "compositional")
# PCoA plot using the weighted UniFrac as distance
ordination_16s_2 = ordinate(s16_re, method="PCoA", "wunifrac")
# Plot
baseplot_16s_2 <- plot_ordination(s16_re, ordination_16s_2,
                                  color="especies",
                                  type = "samples") + #shape = "sexo" shape = "tratamiento"
  ggtitle("Bacterias")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 3)+
  scale_color_discrete(name = "Especies", 
                     labels = c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                                "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))


  #scale_x_discrete(labels=c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
   #                         "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))
  #labs(x = "PCoa1 (31.3%)", y = "PCoA2 (10.4%)")+
  #scale_color_discrete(name = "Condition", 
  #                    labels = c("Domesticated","Wild"))+
  #scale_color_discrete(name = "Treatment",
  #                     labels = c("No-visited","Visited"))

baseplot_16s_2
#
wei_16s = baseplot_16s_2 + stat_ellipse(type = "norm", linetype = 2)
#+ coord_fixed()

#+ stat_ellipse(type = "t")
baseplot_16s_2 + stat_ellipse(type = "t")

#
#baseplot2 +
#  geom_point(size = 1) +
#  geom_label(aes(label = clave), size = 3.5, color = "black")
#
baseplot_16s_2 +
  geom_point(size = 2)+
  ggtitle("")+
  facet_wrap(~sexo, labeller = labeller(sexo = sexos))+
  stat_ellipse(type = "norm", linetype = 2)


ggarrange(wei_16s,
          wei_its, nrow = 2)

                         #### PERMANOVA ####
s16_6_rel <- microbiome::transform(ps.rarefied100, "compositional")
metadf = data.frame(sample_data(ps.rarefied100))
# calculating Weighted-UniFrac
wunifrac.dist <- UniFrac(s16_6_rel, 
                        weighted = TRUE, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)
# Run PERMANOVA
set.seed(100)   # Set the seed to make the analysis reproducible
perm1 = adonis2(wunifrac.dist ~ especies*sexo, data = metadf,
                by = "terms", permutations = 999)
perm1        

# calculating Unweighted-UniFrac
unifrac.dist <- UniFrac(s16_6_rel, 
                        weighted = FALSE, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)
# Run PERMANOVA
set.seed(100)   # Set the seed to make the analysis reproducible
perm2 = adonis2(unifrac.dist ~ especies*sexo, data = metadf,
                by = "terms", permutations = 999)
perm2  

?betadisper
ps.disper <- betadisper(unifrac.dist, metadf$especies)
permutest(ps.disper, pairwise = TRUE)
anova(ps.disper)
TukeyHSD(ps.disper)

#### ANALYSIS WITH RAREFIED DATA ####
### Calculating relative abundances
mergedGP <- merge_samples(ps.rarefied100, "especies")
# Use psmelt to obtain a long-format data.frame
phy <- mergedGP %>% tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
View(phy)
#write.csv(phy, file = "abundancias.csv")

#mergedtop = merge_samples(topfam, "especies")
#phy2 <- mergedtop %>% tax_glom(taxrank = "Family") %>%
#  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
#View(phy2)

#### HEATMAPs of most abundant bacteria ####
### Collapse dataset taxonomy to 10 most abundant families
topfam <- aggregate_top_taxa2(ps.rarefied100, "Family", top = 10)
#View(otu_table(topfam))
topfam

# transform families dataset to relative abundance
transform <- microbiome::transform
famrel <- transform(topfam, "compositional")
#View(otu_table(famrel))

### Define colors for heatmap
heat.cols <- c("#F4D166","#146C36")
### Create a gradient color palette for abundance
grad_ab <- colorRampPalette(c("#F4D166","#146C36"))
heat.cols <- grad_ab(10)

paletteer_c("ggthemes::Green-Gold", 30)

### Heatmap for 10 most abundant bacterial families between species
heatmap_k <- simple_heatmap(famrel,
                            group.facet = "especies",
                            group.order = NULL,
                            abund.thres = 0.01,
                            prev.thres = 0.01,
                            level = "Family",
                            scale.color = "log10",
                            na.fill = "white",
                            color.fill = heat.cols,
                            taxa.arrange=TRUE,
                            remove.other=TRUE,
                            panel.arrange="grid",
                            ncol=NULL,
                            nrow=NULL)
heatmap_k


# Separando por sexos florales
hembras <- subset_samples(ps.rarefied100, sexo == "hembra")
machos <- subset_samples(ps.rarefied100, sexo == "macho")

### Calculating relative abundances for Female flowers
merged_hembras <- merge_samples(hembras, "especies")
# Use psmelt to obtain a long-format data.frame
phy_hembras <- merged_hembras %>% tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
#View(phy_hembras)
write.csv(phy_hembras, file = "abundancias_hembras.csv")

### Calculating relative abundances for Male flowers
merged_machos <- merge_samples(machos, "especies")
# Use psmelt to obtain a long-format data.frame
phy_machos <- merged_machos %>% tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
#View(phy_machos)
write.csv(phy_machos, file = "abundancias_machos.csv")


# Flores hembra
topfam_hem <- aggregate_top_taxa2(hembras, "Family", top = 10)
# transform families dataset to relative abundance
transform <- microbiome::transform
famrel_hem <- transform(topfam_hem, "compositional")
# Heatmap de flores hembra
heatmap_hem <- simple_heatmap(famrel_hem,
                            group.facet = "especies",
                            group.order = NULL,
                            abund.thres = 0.01,
                            prev.thres = 0.01,
                            level = "Family",
                            scale.color = "log10",
                            na.fill = "white",
                            color.fill = heat.cols,
                            taxa.arrange=TRUE,
                            remove.other=TRUE,
                            panel.arrange="grid",
                            ncol=NULL,
                            nrow=NULL)
heatmap_hem

# Flores machos
topfam_mac <- aggregate_top_taxa2(machos, "Family", top = 10)
# transform families dataset to relative abundance
transform <- microbiome::transform
famrel_mac <- transform(topfam_mac, "compositional")
# Heatmap de flores hembra
heatmap_mac <- simple_heatmap(famrel_mac,
                              group.facet = "especies",
                              group.order = NULL,
                              abund.thres = 0.01,
                              prev.thres = 0.01,
                              level = "Family",
                              scale.color = "log10",
                              na.fill = "white",
                              color.fill = heat.cols,
                              taxa.arrange=TRUE,
                              remove.other=TRUE,
                              panel.arrange="grid",
                              ncol=NULL,
                              nrow=NULL)
heatmap_mac

### Arrange heatmaps in one graphic
library(ggpubr)
library(gridExtra)
heatmaps_16s <- grid.arrange(heatmap_hem, heatmap_mac, nrow = 2,
                             labels="AUTO")

heatmaps_16s = ggarrange(heatmap_hem, heatmap_mac, nrow = 2,
                         labels=c("Hembras", "Machos"))
heatmaps_16s

# para cartel
heatmaps_cartel <- grid.arrange(heatmap_hem, heatmap_mac,
                                heatmap_hemits, heatmap_macits,
                                nrow = 2, ncol = 2)


# Most dominant families en todas las muestras
domi <- dominant_taxa(ps.rarefied100,level = "Family") #group = "especies"
View(domi$dominant_overview)
write.csv(domi$dominant_overview, file = "dominant_overview.csv")


#### CORE MICROBIOTA ####

# PREVALENCE OF TAXONOMIC GROUPS
ps.rel = microbiome::transform(ps.rarefied100, "compositional")
# Relative population frequencies; at 1% compositional abundance threshold:
head(prevalence(ps.rel, detection = 10/100, sort = TRUE))
# Absolute population frequencies (sample count)
head(prevalence(ps.rel, detection = 1/100, sort = TRUE, count = TRUE))

# Core tresholds: minimum ASV relative abundance 0.001 % minimum prevalence
# among samples 10%
core.taxa = core_members(ps.rel, detection = .001, prevalence = .1,
                         include.lowest = T)
print(core.taxa)
# A full phyloseq object of the core microbiota is obtained as follows:
full_core <- core(ps.rel, detection = 0.001, prevalence = .1)
View(tax_table(full_core))
# Retrieving the core taxa names from the phyloseq object:
core.taxa2 = taxa(full_core)

# CORE ABUNDANCE AND DIVERSITY
# total core abundance in each sample (sum of abundance of the core members):
core.abundance <- sample_sums(core(ps.rel, detection = .001, prevalence = .1))


prevdf = apply(X = otu_table(s16_6),
               MARGIN = ifelse(taxa_are_rows(s16_6), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(s16_6),
                    tax_table(s16_6))
# plotting prevalence
table(prevdf$Prevalence)
ggplot(prevdf, aes(Prevalence)) + geom_histogram(bins = 50) + theme_bw() + 
  ggtitle("Histogram of Taxa Prevalence")

View(head(otu_table(s16_6)))
head(prevdf)

# CORE VISUALIZATION
# With compositional (relative) abundances
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)
#ggplot(d) + geom_point(aes(x, y)) + scale_x_continuous(trans="log10", limits=c(NA,1))

plot_core(ps.rel, 
          prevalences = prevalences, 
          detections = det, 
          plot.type = "lineplot") + 
  xlab("Relative Abundance (%)")

### CORE HEATMAPS
library(RColorBrewer)
library(reshape)

prevalences <- seq(.05, 1, .05)

detections <- round(10^seq(log10(0.01), log10(.2), length = 9), 3)

# Also define gray color palette
gray <- gray(seq(0,1,length=5))

#Added pseq.rel, I thin... must be checked if it was in the the rednred version,; where it is initialized
#pseq.rel<- microbiome::transform(pseq, 'compositional')
#min-prevalence gets the 100th highest prevalence
p <- plot_core(ps.rel,
               plot.type = "heatmap", 
               colours = gray,
               prevalences = prevalences, 
               detections = detections, 
               min.prevalence = prevalence(ps.rel, sort = TRUE)[100]) +
  labs(x = "Detection Threshold\n(Relative Abundance (%))") +
  
  #Adjusts axis text size and legend bar height
  theme(axis.text.y= element_text(size=8, face="italic"),
        axis.text.x.bottom=element_text(size=8),
        axis.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))

print(p)

# Plotting with upsetR
#install.packages("UpSetR")
library(MicrobiotaProcess)
library(UpSetR)


uppen <- get_upset(ps.rarefied100, factorNames = "especies")
View(uppen)
write.csv(uppen, file = "shared.csv")
# función para gráficar los ASVs compartidos entre especies
?upset

upset(uppen, nsets = 6, mainbar.y.label = "Number of ASVs shared",
      order.by = "freq",
      queries = list(list(query = intersects,
                         params = list("6CM", "5CAA", "4CPP"),
                          color = "orange3",
                         active = T),
                     list(query = intersects,
                          params = list("6CM", "5CAA", "4CPP", 
                                        "3CAS", "2CPF", "1CF"),
                          color = "green",
                          active = T),
                     list(query = intersects,
                          params = list("3CAS", "2CPF", "1CF"),
                          color = "blue",
                          active = T)))
                        #text.scale = 1.1, 
                        #point.size = 2.8, 
                        #line.size = 1

###############################################################################
###############################################################################
library(ggalluvial)
library(pals)

# objeto phyloseq
s16_5

# To plot 'alluvium = Genus' later
# Merges ASVs that have the same taxonomy rank (Genus)
genus2 = tax_glom(s16_5, taxrank = "Genus")

# Merge by SubjectID
merged2 = merge_samples(genus2, "condicion")
merged2 #taxa are columns


# Calculate relative abundance
genus.rb2 = transform_sample_counts(genus2, function(x) {x/sum(x)} )
merged.rb2 = transform_sample_counts(merged2, function(x) {x/sum(x)} )

# Retrieve top5 taxa names per species
mat2 = otu_table(merged.rb2)
top5 = list()
for(i in 1:nrow(mat2)) {
  # Sort by relative abundance and select top 5
  top5[[rownames(mat2)[i]]] = as.data.frame(t(mat2[i,])) %>%
    arrange(desc(.[1])) %>% rownames %>% head(10)
}

length(top5)
head(top5)

# Build genus.rb data.frame
df2 = psmelt(merged.rb2)
head(df2)
df2$Abundance

# Example 1: Combine top5 from all especies
df2$Taxa1 = df2$Genus #modificar nombre 
df2$Taxa1[!df2$OTU %in% unique(unlist(top5))] = "Others"

alluvia_data3 = df2 %>% group_by(Sample, especies, Taxa1) %>%
  summarise(Abundance = sum(Abundance)) %>% as.data.frame
#> `summarise()` has grouped output by 'SampleID', 'SubjectID', 'Sampledetails'.
#> You can override using the `.groups` argument.
head(alluvia_data3)
alluvia_data3$Abundance
#alluvia_data2 = subset(alluvia_data1, Abundance != "NaN")
#alluvia_data2$Abundance

ggplot(data = alluvia_data3,
       aes(x = Sample, y = Abundance, alluvium = Taxa1)) +
  geom_alluvium(aes(fill = Taxa1, color = Taxa1), alpha = 0.5)+
  #facet_wrap(~ especie, scales = "free_x") +
  scale_fill_manual("Genus", values = pals::cols25()) +
  scale_color_manual("Genus", values = pals::cols25())+ 
  labs(x = "Condition")+
  scale_x_discrete(labels = c("Domesticated","Wild"))
#ggtitle("Top5 from all subjects")

#############################################################################
#install.packages("eulerr") # If not installed
library(eulerr)
library(microbiome)
#devtools::install_github('microsud/microbiomeutilities')
library(microbiomeutilities)


# Primero nos aseguramos que estamos usando las funciones del paquete correcto
transform <- microbiome::transform
# Proportion Normalization
pseq_16s <- transform(s16_5, "compositional")

#
trata_states <- unique(as.character(meta(pseq_16s)$tratamiento))
print(trata_states)

# Write a for loop to go through each of the disease_states one by one and 
# combine identified core taxa into a list.
list_core <- c() # an empty object to store information

for (n in trata_states){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq_16s, tratamiento == n) 
  # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.0001, # 0.001 in at least 90% samples 
                         prevalence = .1)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) 
  # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

print(list_core)

# Specify colors and plot venn
# supplying colors in the order they appear in list_core
mycols <- c(nonCRC="#d6e2e9", CRC="#cbf3f0", H="#fcf5c7") 

plot(venn(list_core),
     fills = mycols)
###############################################################################



