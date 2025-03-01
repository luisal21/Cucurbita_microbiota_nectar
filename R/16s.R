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


                    #### METODOS DE NORMALIZACION ####
                         #### RAREFIED DATA ####
otu_tab = t(abundances(s16_6))
p = vegan::rarecurve(otu_tab,
                     step = 50, label = FALSE,
                     sample = min(rowSums(otu_tab),
                                  col = "blue, cex = 0.6"))

#vegan::rarecurve(t(otu_table(s16_6)), step=50, cex=0.5)

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
                             50, rngseed = 1,
                             replace = FALSE)             
ps.rarefied50 # use this rarefied depth
# phyloseq object with normalized library size has
# 785 ASVs in 168 samples
# 757 asvs in 168 samples 26/02/2025
reads_sample2 = readcount(ps.rarefied50)
reads_sample2[1:20]
summarize_phyloseq(ps.rarefied50)
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
div.df$diversity_shannon
div.df$observed
div.df$Phylogenetic_Diversity
# Histograma
hist(div.df$diversity_shannon) 
hist(div.df$observed)
hist(div.df$Phylogenetic_Diversity)
# Shapiro test
shapiro.test(div.df$diversity_shannon)
shapiro.test(div.df$observed)
shapiro.test(div.df$Phylogenetic_Diversity)


### STATISTICAL ANALYSIS OF ALPHA DIVERSITY
library(emmeans)
# Richness
mod1 = glm(observed ~ especies*sexo*tratamiento, data = div.df,
           family = poisson)
summary(mod1)
anova(mod1)
# post-hoc comparisons
compa1 = emmeans(mod1, 
        specs = pairwise ~ especies | sexo | tratamiento, 
        type = "response")
#write.csv(compa1$contrasts, file = "observed_compa_rarefied.csv")

# Shannon
mod2 = glm(diversity_shannon ~ especies*sexo*tratamiento, data = div.df)
summary(mod2)
anova(mod2)
# post-hoc comparisons
compa2 = emmeans(mod2, 
                 specs = pairwise ~ especies | sexo | tratamiento, 
                 type = "response")
#modx = aov(diversity_shannon ~ especies*sexo*tratamiento, data = div.df)
#summary(modx)

# Phylogenetic Diversity index
mod3 = glm(Phylogenetic_Diversity ~ especies*sexo*tratamiento, data = div.df)
summary(mod3)
anova(mod3)
# post-hoc comparisons
compa3 = emmeans(mod3, 
                 specs = pairwise ~ especies | sexo | tratamiento, 
                 type = "response")

######## Analysis with kruskal-Wallis and Wilcox text
# by species
kruskal.test(div.df$observed~div.df$especies)
pairwise.wilcox.test(div.df$observed, div.df$especies, exact=F)
# by sex
kruskal.test(div.df$observed~div.df$sexo)
pairwise.wilcox.test(div.df$observed, div.df$sexo, exact=F)
# by treatment
kruskal.test(div.df$observed~div.df$tratamiento)
pairwise.wilcox.test(div.df$observed, div.df$tratamiento, exact=F)


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
  geom_jitter(position=position_jitter(0.2), size = 3)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Especies", y = "ASV observados")+
  scale_color_discrete(name = "Especies", 
                       labels = c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  scale_x_discrete(labels=c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  theme_bw()

### stripchart de Observed ASVs facet with 4 panels
ggplot(div.df, aes(x=especies, y=observed, color=especies)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 3)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Especies", y = "ASV observados")+
  scale_color_discrete(name = "Especies", 
                       labels = c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  scale_x_discrete(labels=c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  theme_bw()+
  facet_wrap(~sexo+tratamiento)
#stat_summary(fun=mean, geom="point", shape=18,
#            size=6, color="red")

### stripchart Shannon diversity index by species
ggplot(div.df, aes(x=especies, y=diversity_shannon, color=especies))+ 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 3)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Especies", y = "Shannon",
       title = "")+
  scale_color_discrete(name = "Especies", 
                       labels = c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  scale_x_discrete(labels=c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  theme_bw()+
  facet_wrap(~sexo+tratamiento)
  #stat_summary(fun=mean, geom="point", shape=18,
  #            size=6, color="red")
  # shannon_16s_1 = shannon_16s_1 + 
  #stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = 3)

### stripchart phylogenetic diversity index
ggplot(div.df, aes(x=especies, y=Phylogenetic_Diversity, color=especies))+ 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 3)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Especies", y = "PD",
       title = "")+
  scale_color_discrete(name = "Especies", 
                       labels = c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  scale_x_discrete(labels=c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  theme_bw()+
  facet_wrap(~sexo+tratamiento)
#stat_summary(fun=mean, geom="point", shape=18,
#            size=6, color="red")
# shannon_16s_1 = shannon_16s_1 + 
#stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = 3)

# cargando paqueteria gridExtra
library(gridExtra)
library(cowplot)
gt <- grid.arrange(shannon_hongos2,shannon_16s, 
                   ncol = 2, nrow = 1)


gt <- grid.arrange(shannon_hongos2, shannon_16s, 
                   ncol = 2, nrow = 1)


gt3 <- grid.arrange(shannon_hongos2, shannon_hongos4, shannon_hongos3,
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


                         ##### BETA DIVERSITY #####
library(vegan)
# PCoA plot using the unweighted UniFrac as distance
ordination_16s = ordinate(ps.rarefied50, method="PCoA", "unifrac",
                          weighted = F)

baseplot_16s <- plot_ordination(ps.rarefied50, ordination_16s,
                                color="especies",
                                type = "samples",
) + #shape = "sexo"
  ggtitle("")+
  theme(aspect.ratio=1)+
  theme_bw()
baseplot_16s
#
baseplot_16s + geom_point(size = 2)
#
#baseplot1 +
#  geom_point(size = 1) +
#  geom_label(aes(label = clave), size = 3.5, color = "black")
#
baseplot_16s +
  geom_point(size = 2)+
  ggtitle("")+
  facet_wrap(~sexo+tratamiento)

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


  labs(x = "PCoa1 (31.3%)", y = "PCoA2 (10.4%)")+
  #scale_color_discrete(name = "Condition", 
  #                    labels = c("Domesticated","Wild"))+
  scale_color_discrete(name = "Treatment",
                       labels = c("No-visited","Visited"))

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
  facet_wrap(~tratamiento,2)


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
s16_6_rel <- microbiome::transform(s16_6, "compositional")
metadf = data.frame(sample_data(s16_6))
# calculating Weighted-UniFrac
wunifrac.dist <- UniFrac(s16_6_rel, 
                        weighted = TRUE, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)
# Run PERMANOVA
set.seed(100)   # Set the seed to make the analysis reproducible
perm1 = adonis2(wunifrac.dist ~ especies*sexo*tratamiento, data = metadf,
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
perm2 = adonis2(unifrac.dist ~ especies*sexo*tratamiento, data = metadf,
                by = "terms", permutations = 999)
perm2  

?betadisper
ps.disper <- betadisper(unifrac.dist, metadf$especies)
permutest(ps.disper, pairwise = TRUE)
anova(ps.disper)
TukeyHSD(ps.disper)


#### ANALYSIS WITH RAREFIED DATA ####
#### HEATMAPs of most abundant bacteria ####

### Collapse dataset taxonomy to 10 most abundant families

topfam <- aggregate_top_taxa2(ps.rarefied100, "Family", top = 10)
View(otu_table(topfam))
# transform families dataset to relative abundance
transform <- microbiome::transform
famrel <- transform(topfam, "compositional")
View(otu_table(famrel))

### Define colors for heatmap
heat.cols <- c("#a8dadc","#457b9d", "#1d3557")
### Create a gradient color palette for abundance
grad_ab <- colorRampPalette(c("#96d4ca","#d3f3f1", "#7c65a9"))
heat.cols <- grad_ab(10)

### Heatmap for 10 most abundant bacterial families in king penguin body sites

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

#### TOP FIVE FAMILIES ####
library(microbiomeutilities)

kcols <- c("steelblue", "grey50", "red2", "green4", "yellow4", "pink4")

pn <- plot_taxa_boxplot(ps.rarefied50,
                        taxonomic.level = "Family",
                        top.otu = 5, 
                        group = "especies",
                        add.violin= FALSE,
                        title = "Top five families", 
                        keep.other = FALSE,
                        #group.order = c("back","chest","foot"),
                        group.colors = kcols,
                        dot.size = 1)

print(pn + theme_biome_utils())

pn$data
View(pn$data)










#### ALPHA DIVERSITY WITH PROPORTIONS DATA  ####
library(knitr)

# Primero nos aseguramos que estamos usando las funciones del paquete correcto
transform <- microbiome::transform
# Proportion Normalization
pseq_16s <- transform(s16_6, "compositional")

# Obteniendo los valores de las metricas de diversidad alfa
div1 = microbiome::alpha(pseq_16s, index = "all")
kable(head(div1))
View(div1)

tab_16 <- microbiome::alpha(pseq_16s, index = c("observed","shannon", "simpson",
                                                "absolute", "low_abundance",
                                                "rare_abundance"))
kable(head(tab_16))
#write.csv(tab, file = "stats_alfa.csv")

#metadatos
ps1.meta <- meta(pseq_16s)
kable(head(ps1.meta))


#Add the diversity table to metadata
ps1.meta$Shannon <- tab_16$diversity_shannon 
ps1.meta$observed <- tab_16$observed
ps1.meta$even_simp <- tab_16$evenness_simpson
ps1.meta$dom_simp <- tab_16$dominance_simpson
ps1.meta$dom_abs <- tab_16$dominance_absolute
ps1.meta$low <- tab_16$rarity_low_abundance
ps1.meta$rare <- tab_16$rarity_rare_abundance
head(ps1.meta)




# Riqueza
tab <- richness(pseq_16s)
kable(head(tab))
levels(as.factor(ps1.meta$especie))

# Dominancia
tab <- dominance(pseq_16s, index = "all")
kable(head(tab))
tab$absolute

# Most dominant
dominant(pseq_16s)

# rarity
tab <- rarity(pseq_16s, index = "all")
kable(head(tab))
hist(tab$low_abundance)
hist(tab$rare_abundance)
tab$low_abundance











                        
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

