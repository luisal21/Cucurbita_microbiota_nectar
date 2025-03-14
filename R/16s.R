### Análisis de las secuencias de 16s
# Luis Alberto Villanueva Espino


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
setwd("/home/luis/Documents/papers_publicados/paper2_tesisDoctorado/Cucurbita_microbiota_nectar/data")

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
head(seqtab.16s)
tail(seqtab.16s)
str(seqtab.16s)
# Debido a que al cargar el archivo agrega una X al principio del nombre de
# cada columna, hay que decirle que los nombres de las columnas son iguales
# a las filas de sample_id de metadatos para corregir este error
colnames(seqtab.16s) = meta_16s$sample_id
str(seqtab.16s)
# Volver a cargar archivo de metadatos
meta_16s <- read.csv("meta_16s_2024.csv", header = T, row.names = 1)
head(meta_16s)
# cargando taxonomia
taxo_16s <-read.csv("taxa_16s_species_2024.csv", header = T,row.names = 1)
head(taxo_16s)
# corrobando match entre archivos
setdiff(rownames(taxo_16s), rownames(seqtab.16s))
all(rownames(seqtab.16s) == rownames(taxo_16s))
all(colnames(seqtab.16s) == unlist(meta_16s[,'sample.id']))

# Construct phyloseq object (straightforward from dada2 outputs)
s16_1 <- phyloseq(otu_table(as.matrix(seqtab.16s), taxa_are_rows=T),
                  tax_table(as.matrix(taxo_16s)), 
                  sample_data(meta_16s))
s16_1

# Building phylogenetic tree with the function rtree of package "ape"
arbol_16s = rtree(ntaxa(s16_1), rooted=TRUE, tip.label=taxa_names(s16_1))
#plot(arbol_16s)
# Uniendo el arbol con el objeto phyloseq anterior
s16_2 = merge_phyloseq(s16_1, arbol_16s)
s16_2

# Filtering samples of Morelia site
s16_3 = subset_samples(s16_2, sitio=="morelia")
# Removing Chloroplast and Mitochondria sequences and na samples in condition
s16_4 = subset_samples(s16_3, condicion !="na")
s16_4 = subset_taxa(s16_4, Order != "Chloroplast")
s16_4 = subset_taxa(s16_4, Family != "Mitochondria")
# Hay algunas muestras que quedaron sin secuencias, ya que todas sus secuencias
# pertenecian al Reino Viridiplantae. Por lo que voy a excluir esas 
# muestras que quedaron con 0 secuencias:
s16_5 = prune_samples(sample_sums(s16_4)>1, s16_4)
s16_5
# A useful next step is to explore feature prevalence in the dataset,
# which we will define here as the number of samples in which a taxa appears
# at least once.

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(s16_5),
               MARGIN = ifelse(taxa_are_rows(s16_5), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(s16_5),
                    tax_table(s16_5))
#View(prevdf)
#write.csv(prevdf, file = "prevalencia.csv")

plyr::ddply(prevdf, "Phylum", function(df1)
{cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Filtering phyla whose prevalence is 0
filterPhyla = c("Campylobacterota", "Crenarchaeota", "Euryarchaeota",
                "Fibrobacterota", "Patescibacteria", "Synergistota")
# Filter entries with low prevalence Phylum.
s16_6 = subset_taxa(s16_5, !Phylum %in% filterPhyla)
s16_6 # Final ASV data set: 3116 ASVs

# Primero resumimos el contenido del objeto phyloseq
summarize_phyloseq(s16_6) # reads totales 155,244
sort(sample_sums(s16_6))

# información del número de reads
summary(sample_sums(s16_6))
range(sample_sums(s16_6))
#Library size spans 2 - 24079 reads 

# Obteniendo el número de read por muestra
reads_sample = readcount(s16_6)
reads_sample[1:5]
#write.csv(reads_sample, file = "reads_sample_16s_todas.csv")

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


                    #### METODOS DE NORMALIZACION ####
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

### Normalize library size to 100 reads
set.seed(111)
ps.rarefied100 <- rarefy_even_depth(s16_6, 
                                   100, rngseed = 1,
                                   replace = FALSE)             
ps.rarefied100
# phyloseq object with normalized library size has
# 906 ASVs in 134 samples
reads_sample2 = readcount(ps.rarefied100)
reads_sample2[1:20]
summarize_phyloseq(ps.rarefied100)
summary(sample_sums(ps.rarefied100))
abundances(ps.rarefied100)

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
barplot(sample_sums(ps.rarefied50), las = 2)

# taxa prevalence
plot_taxa_prevalence(s16_6, "Phylum")
plot_taxa_prevalence(ps.rarefied50, "Phylum")

                #### ALPHA DIVERSITY WITH RAREFIED DATA  ####
library(knitr)

div1 = microbiome::alpha(ps.rarefied50, index = "all")
kable(head(div1))
View(div1)

# get the metadata out as separate object
ps1.meta <- meta(ps.rarefied50)
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

ps0.rar.asvtab <- as.data.frame(ps.rarefied50@otu_table)

ps0.rar.tree <- ps.rarefied50@phy_tree

# hmp.meta from previous code chunks

# We first need to check if the tree is rooted or not 
ps.rarefied50@phy_tree

# it is a rooted tree
df.pd <- pd(t(ps0.rar.asvtab), ps0.rar.tree,include.root=T) 
# t(ou_table) transposes the table for use in picante and the tre file comes
# from the first code chunck we used to read tree file (see making a phyloseq 
# object section).
datatable(df.pd)

# We add the results of PD to this file and then plot
div.df$Phylogenetic_Diversity <- df.pd$PD
head(div.df)


# data
div.df$observed
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
mod1 = glm(observed ~ especies*sexo, data = div.df,
           family = poisson)
summary(mod1)
anova(mod1)
# post-hoc comparisons (all variables)
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
                 specs = pairwise ~ especies | sexo, 
                 type = "response")
compa2
compa2.1 = emmeans(mod2, 
                 specs = pairwise ~ especies, 
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
                 specs = pairwise ~ especies | sexo, 
                 type = "response")
compa3
compa3.1 = emmeans(mod3, 
                 specs = pairwise ~ especies, 
                 type = "response")
compa3.1
######## Analysis with kruskal-Wallis and Wilcox text
# by species
kruskal.test(div.df$observed~div.df$especies)
pairwise.wilcox.test(div.df$observed, div.df$especies, exact=F)
# by sex
kruskal.test(div.df$observed~div.df$sexo)
pairwise.wilcox.test(div.df$observed, div.df$sexo, exact=F)

##### VISUALIZATION  OF ALPHA DIVERSITY with rarefied data####
# Now use this data frame to plot 
p.shannon <- boxplot_alpha(ps.rarefied50, 
                           index = "shannon",
                           x_var = "especies")
p.shannon
#fill.colors = c(female="cyan4", male="deeppink4"))

p.shannon <- p.shannon + theme_minimal() + 
  labs(x="Species", y="Shannon diversity") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16))
p.shannon

# Stripchart of observed ASVs with 1 panel
ggplot(div.df, aes(x=especies, y=observed, color=especies)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 2)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Richness")+
  scale_color_discrete(name = "Species", 
                       labels = c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  scale_x_discrete(labels=c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  theme_bw()

### stripchart de Observed ASVs facet with 2 panels (divided by sex)
sexos <- c("Female", "Male")
names(sexos) <- c("hembra", "macho")

riqueza_16 <- ggplot(div.df, aes(x=especies, y=observed, color=especies)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 1)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Richness")+
  scale_color_discrete(name = "Species", 
                       labels = c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  scale_x_discrete(labels=c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  theme_bw()+
  facet_wrap(~sexo, labeller = labeller(sexo = sexos))
#stat_summary(fun=mean, geom="point", shape=18,
#            size=6, color="red")

### stripchart Shannon diversity index facet with 2 panels (divided by sex)
shannon_16 <- ggplot(div.df, aes(x=especies, y=diversity_shannon, color=especies))+ 
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
  theme_bw()+
  facet_wrap(~sexo,
             labeller = labeller(sexo = sexos))
  #stat_summary(fun=mean, geom="point", shape=18,
  #            size=6, color="red")
  # shannon_16s_1 = shannon_16s_1 + 
  #stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = 3)

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
gt <- grid.arrange(riqueza_16, shannon_16, pd_16, 
                   ncol = 1, nrow = 3)

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

citation()

            ########## BETA DIVERSITY  WITH RAREFIED DATA ##########
library(vegan)
#### Unweighted UniFrac ####
# PCoA plot using the unweighted UniFrac as distance
ordination_16s = ordinate(ps.rarefied50, method="PCoA", "unifrac",
                          weighted = F)

baseplot_16s <- plot_ordination(ps.rarefied50, ordination_16s,
                                color="especies",
                                type = "samples",
) + #shape = "sexo"
  ggtitle("")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 2)
baseplot_16s
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
ps.rarefied50
s16_re = microbiome::transform(ps.rarefied50, "compositional")
# PCoA plot using the weighted UniFrac as distance
ordination_16s_2 = ordinate(s16_re, method="PCoA", "wunifrac")
# Plot
baseplot_16s_2 <- plot_ordination(s16_re, ordination_16s_2,
                                  color="especies",
                                  type = "samples") + #shape = "sexo" shape = "tratamiento"
  ggtitle("Bacteria")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 3)
  #labs(x = "PCoa1 (31.3%)", y = "PCoA2 (10.4%)")+
  #scale_color_discrete(name = "Condition", 
  #                    labels = c("Domesticated","Wild"))+
  #scale_color_discrete(name = "Treatment",
  #                     labels = c("No-visited","Visited"))

baseplot_16s_2
#
baseplot_16s_2 + stat_ellipse()
#
#baseplot2 +
#  geom_point(size = 1) +
#  geom_label(aes(label = clave), size = 3.5, color = "black")
#
baseplot_16s_2 +
  geom_point(size = 2)+
  ggtitle("")+
  facet_wrap(~sexo, labeller = labeller(sexo = sexos))+
  stat_ellipse()


# por condicion
# Weighted UniFrac distance matrix
wunifrac_dist_16s = phyloseq::distance(pseq_16s, method="wunifrac")
# PCoA plot using the weighted UniFrac as distance
ordination_16s_2 = ordinate(pseq_16s, method="PCoA", distance=wunifrac_dist_16s)

baseplot_16s_2 <- plot_ordination(pseq_16s, ordination_16s_2,
                                  color="condicion",
                                  type = "samples") + #shape = "sexo" shape = "tratamiento"
  ggtitle("Bacteria")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 3)+
  labs(x = "PCoa1 (31.3%)", y = "PCoA2 (10.4%)")+
  scale_color_discrete(name = "Condition", 
                       labels = c("Domesticated","Wild"))
#scale_color_discrete(name = "Treatment",
#                    labels = c("No-visited","Visited"))

baseplot_16s_2
#
baseplot_16s_2 + stat_ellipse()


# por sexo
# Weighted UniFrac distance matrix
wunifrac_dist_16s = phyloseq::distance(pseq_16s, method="wunifrac")
# PCoA plot using the weighted UniFrac as distance
ordination_16s_2 = ordinate(pseq_16s, method="PCoA", distance=wunifrac_dist_16s)

baseplot_16s_2 <- plot_ordination(pseq_16s, ordination_16s_2,
                                  color="sexo",
                                  type = "samples") + #shape = "sexo" shape = "tratamiento"
  ggtitle("Bacteria")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 3)+
  labs(x = "PCoa1 (31.3%)", y = "PCoA2 (10.4%)")+
  scale_color_discrete(name = "Floral sex", 
                       labels = c("Female","Male"))
#scale_color_discrete(name = "Treatment",
#                    labels = c("No-visited","Visited"))

baseplot_16s_2
#
baseplot_16s_2 + stat_ellipse()


baseplot_16s_2
baseplot_its_2
# Juntando gráficas
gt <- grid.arrange(baseplot_its_2_fungi, baseplot_16s_2, 
                   ncol = 1, nrow = 2)


baseplot_its_2_fungi
baseplot_16s_2

ordu.unwt.uni <- ordinate(s16_6_re, "PCoA", "unifrac", weighted=F)

                         #### PERMANOVA ####
s16_6_rel <- microbiome::transform(ps.rarefied50, "compositional")
metadf = data.frame(sample_data(ps.rarefied50))
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
mergedGP <- merge_samples(ps.rarefied50, "especies")
# Use psmelt to obtain a long-format data.frame
phy <- mergedGP %>% tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
View(phy)
write.csv(phy, file = "abundancias.csv")

#mergedtop = merge_samples(topfam, "especies")
#phy2 <- mergedtop %>% tax_glom(taxrank = "Family") %>%
#  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
#View(phy2)

#### HEATMAPs of most abundant bacteria ####
### Collapse dataset taxonomy to 10 most abundant families
topfam <- aggregate_top_taxa2(ps.rarefied50, "Family", top = 10)
#View(otu_table(topfam))
topfam

# transform families dataset to relative abundance
transform <- microbiome::transform
famrel <- transform(topfam, "compositional")
#View(otu_table(famrel))

### Heatmap for 10 most abundant bacterial families between species
heatmap_k <- simple_heatmap(famrel,
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
heatmap_k


# Separando por sexos florales
hembras <- subset_samples(ps.rarefied50, sexo == "hembra")
machos <- subset_samples(ps.rarefied50, sexo == "macho")

### Calculating relative abundances for Female flowers
merged_hembras <- merge_samples(hembras, "especies")
# Use psmelt to obtain a long-format data.frame
phy_hembras <- merged_hembras %>% tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
View(phy_hembras)
write.csv(phy_hembras, file = "abundancias_hembras.csv")

### Calculating relative abundances for Male flowers
merged_machos <- merge_samples(machos, "especies")
# Use psmelt to obtain a long-format data.frame
phy_machos <- merged_machos %>% tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
View(phy_machos)
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
                            color.fill = NULL,
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
                              color.fill = NULL,
                              taxa.arrange=TRUE,
                              remove.other=TRUE,
                              panel.arrange="grid",
                              ncol=NULL,
                              nrow=NULL)
heatmap_mac

### Arrange heatmaps in one graphic
library(ggpubr)
library(gridExtra)
heatmaps <- grid.arrange(heatmap_hem, heatmap_mac, nrow = 2)

### Save heatmaps

ggsave("heatmaps_june.pdf", 
       plot = heatmaps,
       width = 20, 
       height = 20, 
       units = "cm")

# Most dominant families en todas las muestras
domi <- dominant_taxa(ps.rarefied50,level = "Family") #group = "especies"
View(domi$dominant_overview)
write.csv(domi$dominant_overview, file = "dominant_overview.csv")


#### CORE MICROBIOTA ####

# PREVALENCE OF TAXONOMIC GROUPS
ps.rel = microbiome::transform(ps.rarefied50, "compositional")
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
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MicrobiotaProcess")

#install.packages("UpSetR")
library(MicrobiotaProcess)
library(UpSetR)


uppen <- get_upset(ps.rarefied50, factorNames = "especies")
View(uppen)
# función para gráficar los ASVs compartidos entre especies
?upset

upset(uppen, nsets = 6, mainbar.y.label = "Number of ASVs shared",
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
###############################################################################
                  #### ANALYSIS WITH PROPORTIONS DATA ####
# Primero nos aseguramos que estamos usando las funciones del paquete correcto
transform <- microbiome::transform
# Proportion Normalization
pseq_16s <- transform(s16_6, "compositional")

#### ALPHA DIVERSITY WITH PROPORTIONS DATA  ####
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


# data
div.df2$observed
div.df2$diversity_shannon
div.df2$Phylogenetic_Diversity
# Histograma
hist(div.df2$observed)
hist(div.df2$diversity_shannon) 
hist(div.df2$Phylogenetic_Diversity)
# Shapiro test
shapiro.test(div.df2$observed)
shapiro.test(div.df2$diversity_shannon)
shapiro.test(div.df2$Phylogenetic_Diversity)


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

##### VISUALIZATION  OF ALPHA DIVERSITY with rarefied data####
# Stripchart of observed ASVs with 1 panel
ggplot(div.df2, aes(x=especies, y=observed, color=especies)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 2)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Richness")+
  scale_color_discrete(name = "Species", 
                       labels = c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  scale_x_discrete(labels=c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  theme_bw()

### stripchart de Observed ASVs facet with 2 panels
pro16_1 = ggplot(div.df2, aes(x=especies, y=observed, color=especies)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 1)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Richness")+
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



head(div.df)
### gráfica del indice de Shannon entre especies y sexo
pmult <- ggboxplot(div.df, x = "especies", y = "observed", alpha = .5,
                   fill = "especies", add = "jitter", facet.by = "sexo",
                   panel.labs = list(sexo = c("female",
                                              "male")),
                   short.panel.labs = FALSE)
pmult + stat_compare_means(aes(label = paste0("p = ", ..p.format..)))

labs(x = "Especies", y = "Índice de diversidad de Shannon")+
  scale_fill_discrete(name = "Especies") +
  scale_fill_discrete(name = "Especies", 
                      labels = c("CF","CPF","CPP","CAS", "CAA", "CM"))+
  scale_x_discrete(labels=c("1CF" = "CF", "2CPF" = "CPF",
                            "3CPP" = "CPP", "4CAS" = "CAS",
                            "5CAA" = "CAA", "6CM" = "CM"))

pmult + stat_compare_means(aes(label = paste0("p = ", ..p.format..)))


########## BETA DIVERSITY  WITH PROPORTION DATA ##########
pseq_16s # proportion data
s16_6 # phyloseq object filtered but no rarefied
library(vegan)
#### Unweighted UniFrac ####
# PCoA plot using the unweighted UniFrac as distance
ordination_16s = ordinate(s16_6, method="PCoA", "unifrac",
                          weighted = F)

baseplot_16s <- plot_ordination(s16_6, ordination_16s,
                                color="especies",
                                type = "samples",
) + #shape = "sexo"
  ggtitle("")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 2)
baseplot_16s
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
# PCoA plot using the weighted UniFrac as distance
ordination_16s_2 = ordinate(pseq_16s, method="PCoA", "wunifrac")
# Plot
baseplot_16s_2 <- plot_ordination(pseq_16s, ordination_16s_2,
                                  color="especies",
                                  type = "samples") + #shape = "sexo" shape = "tratamiento"
  ggtitle("Bacteria")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 3)
#labs(x = "PCoa1 (31.3%)", y = "PCoA2 (10.4%)")+
#scale_color_discrete(name = "Condition", 
#                    labels = c("Domesticated","Wild"))+
#scale_color_discrete(name = "Treatment",
#                     labels = c("No-visited","Visited"))

baseplot_16s_2
#
baseplot_16s_2 + stat_ellipse()
#
#baseplot2 +
#  geom_point(size = 1) +
#  geom_label(aes(label = clave), size = 3.5, color = "black")
#
baseplot_16s_2 +
  geom_point(size = 2)+
  ggtitle("")+
  facet_wrap(~sexo, labeller = labeller(sexo = sexos))+
  stat_ellipse()

#### PERMANOVA ####
metadf = data.frame(sample_data(pseq_16s))

# calculating Weighted-UniFrac
wunifrac.dist <- UniFrac(pseq_16s, 
                         weighted = TRUE, 
                         normalized = TRUE,  
                         parallel = FALSE, 
                         fast = TRUE)
# Run PERMANOVA
set.seed(100)   # Set the seed to make the analysis reproducible
perm3 = adonis2(wunifrac.dist ~ especies*sexo, data = metadf,
                by = "terms", permutations = 999)
perm3        

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



#### ANALYSIS WITH PROPORTION DATA ####
### Calculating relative abundances
mergedGP_pro <- merge_samples(s16_6, "especies")
# Use psmelt to obtain a long-format data.frame
phy_pro <- mergedGP_pro %>% tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
View(phy_pro)
write.csv(phy, file = "abundancias_family_propor.csv")

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
View(phy_hembras_pro)
write.csv(phy_hembras, file = "abundancias_hembras.csv")

### Calculating relative abundances for Male flowers
merged_machos_pro <- merge_samples(machos_pro, "especies")
# Use psmelt to obtain a long-format data.frame
phy_machos_pro <- merged_machos_pro %>% tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
View(phy_machos_pro)
write.csv(phy_machos, file = "abundancias_machos.csv")


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

### Save heatmaps

ggsave("heatmaps_june.pdf", 
       plot = heatmaps,
       width = 20, 
       height = 20, 
       units = "cm")

# Most dominant families en todas las muestras
domipro <- dominant_taxa(s16_6,level = "Family") #group = "especies"
View(domipro$dominant_overview)
write.csv(domi$dominant_overview, file = "dominant_overview.csv")


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
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MicrobiotaProcess")

#install.packages("UpSetR")
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
###############################################################################
###############################################################################
### Código que sobra
##### DESeq2 #### No correr por ahora
# install package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)

kostic <- prune_samples(sample_sums(its5) > 10, its5)
kostic

sample_data(kostic)$condicion <- as.factor(sample_data(kostic)$condicion)

# Convert to DESeq2 object
DESeq2_Object <- phyloseq_to_deseq2(kostic, ~ condicion)
#> converting counts to integer mode
#> Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
#> design formula are characters, converting to factors

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(DESeq2_Object), 1, gm_mean)
DESeq2_dds = estimateSizeFactors(DESeq2_Object, geoMeans = geoMeans)
DESeq2_dds = DESeq(DESeq2_dds, fitType="local")
#> using pre-existing size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing

# Explore results
res = results(DESeq2_dds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"),
               as(tax_table(kostic)[rownames(sigtab), ], "matrix"))
head(sigtab)

posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj",
                        "Phylum", "Class", "Family", "Genus")]



#### NORMALIZATION METHODS ####
# Proportion Normalization
# Primero nos aseguramos que estamos usando las funciones del paquete correcto
# Paquete microbiome
transform <- microbiome::transform
s16_prop <- transform(s16_6, "compositional")
otu_table(s16_prop)[1:10, 1:10]
head(sample_sums(s16_prop))

### Abundance value transformation
plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("Bacteroidota"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "sexo",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}

# Transform to relative abundance. Save as new object.
ps3ra = transform_sample_counts(s16_6, function(x){x / sum(x)})

#Now plot the abundance values before and after transformation.
plotBefore = plot_abundance(s16_6,"")
plotAfter = plot_abundance(ps3ra,"")
# Combine each plot into one graphic.
grid.arrange(nrow = 2, plotBefore, plotAfter)


write.csv(otu_table(s16_prop), file = "prop_16s.csv") 

