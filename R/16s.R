### Análisis de las secuencias de 16s

devtools::install_github("microsud/microbiomeutilities")


# Cargando paquete phyloseq y paquetes adicionales
library(phyloseq); packageVersion("phyloseq")
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
library(UpSetR)
library(minpack.lm)
library(Hmisc)
library(stats4)
library(patchwork)


#### Set up directory in which the databases are located ####
setwd("/home/luis/Documents/papers_publicados/paper2_tesisDoctorado/Cucurbita_microbiota_nectar/data")

### Primero transponer la tabla original "seqtab_nochim_16s_2024.csv, despues
# deshabilitar este codigo
### Transponer tabla de ASV
#asvtable <- read.csv("seqtab_nochim_16s_2024.csv", header = T, row.names = 1)
#head(asvtable)
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

# obteniendo el arbol filogenetico
arbol_16s = rtree(ntaxa(s16_1), rooted=TRUE, tip.label=taxa_names(s16_1))
#plot(arbol_16s)
# Uniendo el arbol con el objeto phyloseq anterior
s16_2 = merge_phyloseq(s16_1, arbol_16s)
s16_2
# Si solo filtro morelia
s16_3 = subset_samples(s16_2, sitio=="morelia")
s16_3

# Primero resumimos el contenido del objeto phyloseq
summarize_phyloseq(s16_3) # reads totales 5,757,299

# información del número de reads
summary(sample_sums(s16_3))

##### Extraemos todas aquellas secuencias de bacterias #####
s16_4 = subset_samples(s16_3, condicion!="na")
s16_4
s16_4 = subset_taxa(s16_4, Order != "Chloroplast")
s16_4
s16_4 = subset_taxa(s16_4, Family != "Mitochondria")
s16_4

# Hay algunas muestras que quedaron sin secuencias, ya que todas sus secuencias
# pertenecian al Reino Viridiplantae. Por lo que voy a excluir esas 
# muestras que quedaron con 0 secuencias:
s16_5 = prune_samples(sample_sums(s16_4)>1, s16_4)
s16_5
#Final ASV data set: 4655 ASVs
summary(sample_sums(s16_5))
summarize_phyloseq(s16_5)

# Obteniendo el número de read por muestra
reads_sample = readcount(s16_5)
reads_sample[1:5]
write.csv(reads_sample, file = "reads_sample_16s_todas.csv")

range(sample_sums(s16_5))
#Library size spans 2 - 5691 reads 

# abundancia total a nivel de Phylum
plot_bar(s16_5, fill = "Phylum")+ 
  geom_bar(aes(color = Phylum), stat = "identity", position = "stack")

# Abundancia relativa de secuencias 
s16_5_re = transform_sample_counts(s16_5, function(x) x / sum(x))
s16_5_re
# Abundancia relativa a nivel de Phylum
plot_bar(s16_5_re, fill = "Phylum")+ 
  geom_bar(aes(color = Phylum), stat = "identity", position = "stack")



blancos = subset_samples(s16_1, condicion=="na")
blancos
blancos = subset_taxa(blancos, Order != "Chloroplast")
blancos
blancos = subset_taxa(blancos, Family != "Mitochondria")
blancos
# abundancia total a nivel de Phylum
plot_bar(blancos, fill = "Class")+ 
  geom_bar(aes(color = Class), stat = "identity", position = "stack")
reads_sample = readcount(blancos)
reads_sample[1:5]
write.csv(reads_sample, file = "reads_blancos.csv")

#### NORMALIZATION METHODS ####
# Proportion Normalization
# Primero nos aseguramos que estamos usando las funciones del paquete correcto
# Paquete microbiome
transform <- microbiome::transform
s16_prop <- transform(s16_5, "compositional")
otu_table(s16_prop)[1:10, 1:10]
head(sample_sums(s16_prop))

write.csv(otu_table(s16_prop), file = "prop_16s.csv") 

### OTROS METODOS DE NORMALIZACION
# Curvas de Rarefacción
otu_tab = t(abundances(s16_5))
p = vegan::rarecurve(otu_tab,
                     step = 50, label = FALSE,
                     sample = min(rowSums(otu_tab),
                                  col = "blue, cex = 0.6"))

rarecurve(t(otu_table(s16_5)), step=50, cex=0.5)

### Normalize library size to 1000 reads
ps.rarefied1000 <- rarefy_even_depth(s16_5, 
                                    1000, rngseed = 1,
                                    replace = FALSE)             
ps.rarefied1000
# phyloseq object with normalized library size has
# 784 ASVs in 18 samples

### Normalize library size to 100 reads
ps.rarefied100 <- rarefy_even_depth(s16_5, 
                                   100, rngseed = 1,
                                   replace = FALSE)             
ps.rarefied100
# phyloseq object with normalized library size has
# 1600 ASVs in 123 samples

### Normalize library size to 50 reads
ps.rarefied50 <- rarefy_even_depth(s16_5, 
                             50, rngseed = 1,
                             replace = FALSE)             
ps.rarefied50
# phyloseq object with normalized library size has
# 1421 ASVs in 149 samples
reads_sample2 = readcount(ps.rarefied50)
reads_sample2[1:20]
summarize_phyloseq(ps.rarefied50)
summary(sample_sums(ps.rarefied50))
abundances(ps.rarefied50)

### Normalize library size to 10 reads
ps.rarefied10 <- rarefy_even_depth(s16_5, 
                                   10, rngseed = 1,
                                   replace = FALSE)             
ps.rarefied10
# phyloseq object with normalized library size has
# 824 ASVs in 209 samples

#### HEATMAPs of most abundant bacteria ####

### Collapse dataset taxonomy to 10 most abundant families

topfam <- aggregate_top_taxa2(ps.rarefied50, "Family", top = 20)
View(otu_table(topfam))
# transform families dataset to relative abundance

famrel <- transform(topfam, "compositional")
View(otu_table(famrel))

### Define colors for heatmap
heat.cols <- c("#a8dadc","#457b9d", "#1d3557")
### Create a gradient color palette for abundance
grad_ab <- colorRampPalette(c("#96d4ca","#d3f3f1", "#7c65a9"))
heat.cols <- grad_ab(20)

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















#### ALPHA DIVERSITY ####
library(knitr)

# Primero nos aseguramos que estamos usando las funciones del paquete correcto
transform <- microbiome::transform
# Proportion Normalization
pseq_16s <- transform(s16_5, "compositional")

# Obteniendo los valores de las metricas de diversidad alfa
#div1 = microbiome::alpha(pseq_16s, index = "all")
#kable(head(div1))

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


hist(tab_16$diversity_shannon) 
hist(tab_16$observed)
hist(tab_16$evenness_simpson)
hist(tab_16$dominance_simpson)
hist(tab_16$dominance_absolute)

shapiro.test(ps1.meta$Shannon)
shapiro.test(ps1.meta$observed)
shapiro.test(ps1.meta$even_simp)
shapiro.test(ps1.meta$dom_simp)
shapiro.test(ps1.meta$dom_abs)
length(ps1.meta$Shannon)

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


######## analisis estadistica de diversidad alfa
pairwise.wilcox.test(tab$observed, sample_data(pseq_16s)$especies)
pairwise.wilcox.test(tab$diversity_shannon, sample_data(pseq)$condicion)
pairwise.wilcox.test(tab$observed, sample_data(pseq)$especie)

pairwise.wilcox.test(tab$diversity_shannon, sample_data(pseq)$especie)

# análisis estadístico del tamaño del polen con una prueba de kruskal-wallis
kruskal.test(ps1.meta$Shannon ~ sample_data(ps1.meta)$especies)
# comparaciones multiples
pairwise.wilcox.test(ps1.meta$Shannon, sample_data(ps1.meta)$especies)


pairwise.wilcox.test(tab$observed, sample_data(pseq)$sexo)
pairwise.wilcox.test(tab$diversity_shannon, sample_data(pseq)$sexo)
pairwise.wilcox.test(tab$observed, sample_data(pseq)$tratamiento)
pairwise.wilcox.test(tab$diversity_shannon, sample_data(pseq)$tratamiento)



# VISUALIZATION de diversidad alfa
### stripchart de Observed
ggplot(ps1.meta, aes(x=tratamiento, y=observed, color=tratamiento)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 3)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Especies", y = "ASV observados")+
  scale_color_discrete(name = "Especies", 
                       labels = c("3CPP"="CPP", "4CAS" = "CAS", "5CAA"="CAA",
                                  "6CM"="CM"))+
  scale_x_discrete(labels=c("3CPP"="CPP", "4CAS" = "CAS", "5CAA"="CAA",
                            "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  theme_bw()

#scale_color_discrete(name = "Especies", 
#                     labels = c("1CF" = "CF","2CPF" = "CPF", "3CPP"="CPP",
#                                "4CAS" = "CAS", "5CAA"="CAA", "6CM"="CM"))+
#scale_x_discrete(labels=c("1CF" = "CF","2CPF" = "CPF", "3CPP"="CPP",
#                          "4CAS" = "CAS", "5CAA"="CAA", "6CM"="CM"))+
#scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
#                                                       "Visitadas")) +
# theme_bw()
#  stat_summary(fun=mean, geom="point", shape=18,
#               size=6, color="red")

### stripchart de Shannon por condición
shannon_16s_1 = ggplot(ps1.meta, aes(x=condicion, y=Shannon,
                                     color=condicion)) + 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 3)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Condition", y = "Alpha diversity (Shannon)",
       title = "")+
  scale_color_discrete(name = "Especies", 
                       labels = c("3CPP"="CPP", "4CAS" = "CAS", "5CAA"="CAA",
                                  "6CM"="CM"))+
  scale_x_discrete(labels=c("Domesticated","Wild"))+
  theme_bw()+
  scale_y_continuous(limits=c(0, 3))+
  theme(legend.position = "none")
shannon_16s_1 = shannon_16s_1 + stat_compare_means(method = "wilcox.test", 
                                                   label.x = 1.5, label.y = 3)


shannon_16s_1

# cargando paqueteria gridExtra
library(gridExtra)
library(cowplot)
gt <- grid.arrange(shannon_hongos2,shannon_16s, 
                   ncol = 2, nrow = 1)


gt <- grid.arrange(shannon_hongos2, shannon_16s, 
                   ncol = 2, nrow = 1)


gt3 <- grid.arrange(shannon_hongos2, shannon_hongos4, shannon_hongos3,
                    ncol = 1, nrow = 3)

#### Shannon por tratamientos ####
shannon_16s_2 = ggplot(ps1.meta, aes(x=tratamiento, y=Shannon,
                                     color=tratamiento)) + 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 3)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Treatment", y = "Alpha diversity (Shannon)",
       title = "")+
  scale_color_discrete(name = "Especies", 
                       labels = c("3CPP"="CPP", "4CAS" = "CAS", "5CAA"="CAA",
                                  "6CM"="CM"))+
  scale_x_discrete(labels=c("3CPP"="CPP", "4CAS" = "CAS", "5CAA"="CAA",
                            "6CM"="CM"))+
  theme_bw()+
  scale_y_continuous(limits=c(0, 3))+
  scale_x_discrete(labels = c("No-visited","Visited"))+
  theme(legend.position = "none")
shannon_16s_2 = shannon_16s_2 + stat_compare_means(method = "wilcox.test", 
                                                   label.x = 1.5, label.y = 3)
shannon_16s_2
#########
#### Shannon por sexo ####
shannon_16s_3 = ggplot(ps1.meta, aes(x=sexo, y=Shannon,
                                     color=sexo)) + 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 3)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Sex", y = "Alpha diversity (Shannon)",
       title = "")+
  scale_color_discrete(name = "Especies", 
                       labels = c("3CPP"="CPP", "4CAS" = "CAS", "5CAA"="CAA",
                                  "6CM"="CM"))+
  scale_x_discrete(labels=c("3CPP"="CPP", "4CAS" = "CAS", "5CAA"="CAA",
                            "6CM"="CM"))+
  theme_bw()+
  scale_y_continuous(limits=c(0, 3))+
  scale_x_discrete(labels = c("Female","Male"))+
  theme(legend.position = "none")

shannon_16s_3 = shannon_16s_3 + stat_compare_means(method = "wilcox.test", 
                                                   label.x = 1.5, label.y = 3)

shannon_16s_3

# cargando paqueteria gridExtra
gt4 <- grid.arrange(shannon_16s_1, shannon_16s_3, shannon_16s_2, 
                    ncol = 1, nrow = 3)

#### shannon por especie ####
# primero haces el análisis
compare_means(Shannon ~ especies,  data = ps1.meta)
# Visualize: Specify the comparisons you want
#my_comparisons <- list( c("5CAA", "3CPP"), c("5CAA", "4CAS"), c("5CAA", "6CM"),
#                        c("3CPP", "4CAS"), c("3CPP", "6CM"), c("4CAS", "6CM"))
my_comparisons <- list( c("5CAA", "4CAS"), c("5CAA", "6CM"),
                        c("3CPP", "6CM"))

shannon_16s_3 = ggplot(ps1.meta, aes(x=especies, y=Shannon,
                                     color=especies)) + 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 3)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Especies", y = "Índice de diversidad de Shannon",
       title = "Bacterias")+
  scale_color_discrete(name = "Especies", 
                       labels = c("3CPP"="CPP", "4CAS" = "CAS", "5CAA"="CAA",
                                  "6CM"="CM"))+
  scale_x_discrete(labels=c("3CPP"="CPP", "4CAS" = "CAS", "5CAA"="CAA",
                            "6CM"="CM"))+
  theme_bw()+
  #scale_y_continuous(limits=c(0, 3))+
  theme(legend.position = "none")+
  ## Add pairwise comparisons p-value
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.x = .7 , label.y = 4.5)    # Add global p-value
shannon_16s_3

shannon_16s_3 = shannon_16s_3 + stat_compare_means(method = "wilcox.test", 
                                                   label.x = 1.5, label.y = 3)

shannon_hongos3
shannon_16s_2

# cargando paqueteria gridExtra
gt2 <- grid.arrange(shannon_hongos3, shannon_16s_2, 
                    ncol = 2, nrow = 1)


# graficas por sexo
shannon_hongos4
shannon_16s_4
# cargando paqueteria gridExtra
gt3 <- grid.arrange(shannon_hongos4, shannon_16s_4, 
                    ncol = 2, nrow = 1)


# otra gráfica
p.shannon <- boxplot_alpha(pseq, 
                           index = "shannon",
                           x_var = "especie",
                           fill.colors = c("1CF" = "yellow",
                                           "2CPF" = "blue",
                                           "3CPP"="green",
                                           "4CAS"="cyan4",
                                           "5CAA"="deeppink4",
                                           "6CM"="red"
                           ))
p.shannon

p.shannon <- p.shannon + theme_minimal() + 
  labs(x="Especies", y="Indice de diversidad de Shannon") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16))
p.shannon


# create a list of pairwise comparisons
spp <- base::levels(as.factor(ps1.meta$especie)) # get the variables
# make a pairwise list that we want to compare.
spp.pairs <- combn(seq_along(spp), 2, simplify = FALSE, FUN = function(i)spp[i])
print(spp.pairs)

#
library(ggpubr)
# gráfica del indice de diversida de Shannon por especie
#p1 <- ggboxplot(ps1.meta, x = "especie", y = "Shannon",
#                add = c("boxplot", "jitter"), fill = "especie")+
# labs(x = "Especies", y = "Índice de diversidad de Shannon")
#  scale_fill_discrete(name = "Especies", 
#                    labels = c("CF","CPF","CPP","CAS", "CAA", "CM"))+
#  scale_x_discrete(labels=c("1CF" = "CF", "2CPF" = "CPF",
#                            "3CPP" = "CPP", "4CAS" = "CAS",
#                            "5CAA" = "CAA", "6CM" = "CM"))
#print(p1)
#p1 <- p1 + stat_compare_means(comparisons = spp.pairs) 
#print(p1) # Wilcox test

# Gráfica del número de taxones observados por especie
# p2 <- ggboxplot(ps1.meta, x = "especie", y = "observed",alpha =.5,
#                add = c("boxplot","jitter"), fill = "especie")+
#  labs(x = "Especies", y = "Riqueza observada")+
#  scale_fill_discrete(name = "Especies", 
#                      labels = c("CAS","CAA", "CM", "CPP"))+
#  scale_x_discrete(labels=c("1CAS" = "CAS", "2CAA" = "CAA",
#                            "3CM" = "CM", "4CPP" = "CPP"))
# print(p2)
# p2 <- p2 + stat_compare_means(comparisons = spp.pairs) 
# print(p2) # Wilcox test


### gráfica del indice de Shannon entre especies
pmult <- ggboxplot(ps1.meta, x = "especie", y = "Shannon", alpha = .5,
                   fill = "especie", add = "jitter")+ 
  labs(x = "Especies", y = "Índice de diversidad de Shannon")+
  scale_fill_discrete(name = "Especies", 
                      labels = c("CPP","CAS", "CAA", "CM"))+
  scale_x_discrete(labels=c("2CAA" = "CAA", "4CPP" = "CPP",
                            "1CAS" = "CAS", "3CM" = "CM" ))
pmult
ggpar(pmult, ylim = c(0,3))

pmult + stat_compare_means(aes(label = paste0("p = ", ..p.format..)))

### gráfica del indice de Shannon entre especies y tratamiento
pmult <- ggboxplot(ps1.meta, x = "especie", y = "Shannon", alpha = .5,
                   fill = "especie", add = "jitter", facet.by = "tratamiento",
                   panel.labs = list(tratamiento = c("No visitadas",
                                                     "Visitadas")),
                   short.panel.labs = FALSE)+ 
  labs(x = "Especies", y = "Índice de diversidad de Shannon")+
  scale_fill_discrete(name = "Especies") +
  scale_fill_discrete(name = "Especies", 
                      labels = c("CF","CPF","CPP","CAS", "CAA", "CM"))+
  scale_x_discrete(labels=c("1CF" = "CF", "2CPF" = "CPF",
                            "3CPP" = "CPP", "4CAS" = "CAS",
                            "5CAA" = "CAA", "6CM" = "CM"))

pmult + stat_compare_means(aes(label = paste0("p = ", ..p.format..)))

### gráfica del número de taxones observados entre especies y tratamientos
pmult <- ggboxplot(ps1.meta, x = "especie", y = "observed", alpha = .5,
                   fill = "especie", add = "jitter", facet.by = "tratamiento",
                   panel.labs = list(tratamiento = c("t1", "t2")),
                   short.panel.labs = FALSE)+ 
  labs(x = "Especies", y = "Riqueza observada")+
  scale_fill_discrete(name = "Especies")

pmult + stat_compare_means(aes(label = paste0("p = ", ..p.format..)))


##### gráfica del indice de Shannon por especie y sexo 
pmult <- ggboxplot(ps1.meta2, x = "especie", y = "Shannon", alpha = .5,
                   fill = "especie", add = "jitter", 
                   facet.by = "sexo",
                   panel.labs = list(sexo = c("Hembra",
                                              "Macho")),
                   short.panel.labs = FALSE)+ 
  labs(x = "Especies", y = "Indice de diversidad de Shannon")+
  scale_fill_discrete(name = "Especies",
                      labels = c("CPP","CAS", "CAA", "CM"))+
  scale_x_discrete(labels=c("1CAS" = "CAS", "2CAA" = "CAA",
                            "3CM" = "CM", "4CPP" = "CPP"))+
  theme(legend.position = "top")

pmult + stat_compare_means(aes(label = paste0("p = ", ..p.format..)))

##### número de otus observados por especie, sexo y tratamiento
pmult <- ggboxplot(ps1.meta, x = "especie", y = "observed", alpha = .5,
                   color = "sexo", fill = "especie", add = "jitter", 
                   facet.by = "tratamiento",
                   panel.labs = list(tratamiento = c("Sin visita",
                                                     "Con visita")),
                   short.panel.labs = FALSE)+ 
  labs(x = "Especies", y = "Riqueza observada")+
  scale_fill_discrete(name = "Especies",
                      labels = c("CAS","CAA", "CM", "CPP"))+
  scale_x_discrete(labels=c("1CAS" = "CAS", "2CAA" = "CAA",
                            "3CM" = "CM", "4CPP" = "CPP"))+
  theme(legend.position = "bottom")

pmult + stat_compare_means(aes(label = paste0("p = ", ..p.format..)))


#### BETA DIVERSITY ####
# unweighted UniFrac distance matrix
unwunifrac_dist_16s = phyloseq::distance(pseq_16s, method="uunifrac")
as.matrix(unwunifrac_dist_16s)[1:6, 1:6]

# PCoA plot using the unweighted UniFrac as distance
ordination_16s = ordinate(pseq_16s, method="PCoA", distance=unwunifrac_dist_16s)

baseplot_16s <- plot_ordination(pseq_16s, ordination_16s,
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
  facet_wrap(~tratamiento,2)

# Weighted UniFrac distance matrix
wunifrac_dist_16s = phyloseq::distance(pseq_16s, method="wunifrac")
# PCoA plot using the weighted UniFrac as distance
ordination_16s_2 = ordinate(pseq_16s, method="PCoA", distance=wunifrac_dist_16s)

baseplot_16s_2 <- plot_ordination(pseq_16s, ordination_16s_2,
                                  color="tratamiento",
                                  type = "samples") + #shape = "sexo" shape = "tratamiento"
  ggtitle("Bacteria")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 3)+
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

# PERMANOVA
adonis2(unwunifrac_dist_16s ~ sample_data(pseq_16s)$especies*sample_data(pseq_16s)$sexo*
          sample_data(pseq_16s)$tratamiento, permutations = 100)

adonis2(wunifrac_dist_16s ~ sample_data(pseq_16s)$especies*sample_data(pseq_16s)$sexo*
          sample_data(pseq_16s)$tratamiento, permutations = 100)

# Metadata and estimate Unifrac
metadata <- as(sample_data(pseq_16s), "data.frame")
unifrac_dist_16s <- UniFrac(pseq_16s, weighted = TRUE)

set.seed(100)  # Set the seed to make the analysis reproducible

permanova_16s <- adonis2(unifrac_dist_16s ~ especies*sexo*tratamiento,
                         data = metadata,
                         permutations = 999)
permanova_16s
# por condicion
permanova_16s_2 <- adonis2(unifrac_dist_16s ~ condicion*sexo*tratamiento,
                           data = metadata,
                           permutations = 999)

permanova_16s_2

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





