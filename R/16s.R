################### Análisis de las secuencias de 16s #########################
##### Luis Alberto Villanueva Espino
#### Marzo 2026

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("phyloseq")

#### Cargando paquete phyloseq y paquetes adicionales ####
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
all(colnames(seqtab.16s) == rownames(meta_16s2))

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
#### ANALISIS DE COMPOSICIÓN TAXONÓMICA DE BACTERIAS SEPARADAS POR SEXO FLORAL 

# Separando por sexos florales
hembras_pro <- subset_samples(s16_6, sexo == "hembra")
machos_pro <- subset_samples(s16_6, sexo == "macho")

##### HEMBRAS ####
#### Obteniendo las abundancias relativas por especie
# 1. Preparación del objeto Phyloseq
# A. Aglomerar a nivel de Familia. Esto suma los reads de todos los ASVs/OTUs
# que pertenecen a la misma familia. Reduce la complejidad drásticamente.
ps_fam <- tax_glom(hembras_pro, taxrank = "Family", NArm = FALSE)

# B. Transformar a abundancia relativa (escala 0-1)
# Esto estandariza las muestras independientemente de su profundidad de secuenciación.
ps_rel <- transform_sample_counts(ps_fam, function(x) x / sum(x))

# C. Convertir el objeto phyloseq a un dataframe largo ('melted' dataframe)
df_melt <- psmelt(ps_rel)

# D. Agrupar por Especie y Familia, y calcular la abundancia relativa PROMEDIO
# para esa especie de planta.
df_grouped <- df_melt %>%
  dplyr::group_by(especies, Family) %>% 
  dplyr::summarize(MeanRelAbund = mean(Abundance), .groups = "drop")

# E. Identificar el Top 10 de familias más abundantes POR ESPECIE.
# Usamos slice_max para quedarnos con las 10 filas con mayor abundancia promedio por grupo.
top10_per_species <- df_grouped %>%
  group_by(especies) %>%
  slice_max(order_by = MeanRelAbund, n = 5, with_ties = FALSE) %>%
  mutate(Classification = Family) %>% # Guardamos el nombre para el cruce
  ungroup()

# F. Crear la categoría "Otras".
# Cruzamos la tabla agrupada completa con la tabla del Top 10.
df_final <- df_grouped %>%
  left_join(select(top10_per_species, especies, Family, Classification), 
            by = c("especies", "Family")) %>%
  # Si la familia NO estaba en el top10 (Classification es NA), la llamamos "Otras"
  dplyr::mutate(Taxonomy_Display = ifelse(is.na(Classification), 
                                   "Other", 
                                   as.character(Family))) %>%
  # Volvemos a agrupar y sumar para colapsar todas las "Otras" en una sola fila por especie
  dplyr::group_by(especies, Taxonomy_Display) %>%
  dplyr::summarize(TotalMeanAbund = sum(MeanRelAbund), .groups = 'drop') %>%
  # Convertir abundancia relativa (0-1) a porcentaje (0-100) para el gráfico
  dplyr::mutate(Percentage = TotalMeanAbund * 100)
View(df_final)
#write.csv(df_final, file = "abundancias_family_hembras_marzo2026_10.csv")

# 3. Preparación para Graficar
# ==============================================================================

# G. Ordenar los factores para que "Otras" siempre aparezca abajo en la barra y leyenda.
# Primero obtenemos todos los nombres únicos excepto "Otras" y los ordenamos alfabéticamente
taxa_names <- sort(setdiff(unique(df_final$Taxonomy_Display), "Other"))
# Creamos el factor con "Otras" al final
df_final$Taxonomy_Display <- factor(df_final$Taxonomy_Display, 
                                    levels = c(taxa_names, "Other"))

# H. Definir una paleta de colores. El Top 10 por especie puede generar
# fácilmente entre 15 y 30 familias únicas en total para graficar.
# Necesitamos una paleta grande.
n_colors <- length(unique(df_final$Taxonomy_Display))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
custom_colors <- getPalette(n_colors)

# Asignar explícitamente el color gris a "Otras"
names(custom_colors) <- levels(df_final$Taxonomy_Display)
custom_colors["Other"] <- "#D3D3D3" # Gris claro

# 4. Visualización con GGPLOT2
# ==============================================================================

grafico_taxonomia_ps_H <- ggplot(df_final, aes(x = especies, y = Percentage,
                                               fill = Taxonomy_Display)) +
  # Barras apiladas al 100% (position = "stack")
  geom_bar(stat = "identity", position = "stack", color = "black",
           linewidth = 0.2) +
  
  # Usar la paleta de colores personalizada
  scale_fill_manual(values = custom_colors) +
  
  # Etiquetas y formato
  labs(x = "Species",
    y = "Relative Abundance (%)",
    fill = "Family"
  ) +
  scale_x_discrete(labels=c("2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  
  # Tema limpio para tesis
  theme_classic() +
  theme(
    # Nombres de especies de Cucurbita en cursiva y rotados
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 9)
  )

# Mostrar gráfico
grafico_taxonomia_ps_H

#### MACHOS #####
# 1. Preparación del objeto Phyloseq
# A. Aglomerar a nivel de Familia. Esto suma los reads de todos los ASVs/OTUs
# que pertenecen a la misma familia. Reduce la complejidad drásticamente.
ps_fam_m <- tax_glom(machos_pro, taxrank = "Family", NArm = FALSE)

# B. Transformar a abundancia relativa (escala 0-1)
# Esto estandariza las muestras independientemente de su profundidad de secuenciación.
ps_rel_m <- transform_sample_counts(ps_fam_m, function(x) x / sum(x))

# C. Convertir el objeto phyloseq a un dataframe largo ('melted' dataframe)
df_melt_m <- psmelt(ps_rel_m)

# D. Agrupar por Especie y Familia, y calcular la abundancia relativa PROMEDIO
# para esa especie de planta.
df_grouped_m <- df_melt_m %>%
  dplyr::group_by(especies, Family) %>% 
  dplyr::summarize(MeanRelAbund = mean(Abundance), .groups = "drop")

# E. Identificar el Top 10 de familias más abundantes POR ESPECIE.
# Usamos slice_max para quedarnos con las 10 filas con mayor abundancia promedio por grupo.
top10_per_species_m <- df_grouped_m %>%
  group_by(especies) %>%
  slice_max(order_by = MeanRelAbund, n = 5, with_ties = FALSE) %>%
  mutate(Classification = Family) %>% # Guardamos el nombre para el cruce
  ungroup()

# F. Crear la categoría "Otras".
# Cruzamos la tabla agrupada completa con la tabla del Top 10.
df_final_m <- df_grouped_m %>%
  left_join(select(top10_per_species_m, especies, Family, Classification), 
            by = c("especies", "Family")) %>%
  # Si la familia NO estaba en el top10 (Classification es NA), la llamamos "Otras"
  dplyr::mutate(Taxonomy_Display = ifelse(is.na(Classification), 
                                          "Other", 
                                          as.character(Family))) %>%
  # Volvemos a agrupar y sumar para colapsar todas las "Otras" en una sola fila por especie
  dplyr::group_by(especies, Taxonomy_Display) %>%
  dplyr::summarize(TotalMeanAbund = sum(MeanRelAbund), .groups = 'drop') %>%
  # Convertir abundancia relativa (0-1) a porcentaje (0-100) para el gráfico
  dplyr::mutate(Percentage = TotalMeanAbund * 100)
View(df_final_m)
#write.csv(df_final_m, file = "abundancias_family_machos_marzo2026_10.csv")
# 3. Preparación para Graficar
# ==============================================================================

# G. Ordenar los factores para que "Otras" siempre aparezca abajo en la barra y leyenda.
# Primero obtenemos todos los nombres únicos excepto "Otras" y los ordenamos alfabéticamente
taxa_names_m <- sort(setdiff(unique(df_final_m$Taxonomy_Display), "Other"))
# Creamos el factor con "Otras" al final
df_final_m$Taxonomy_Display <- factor(df_final_m$Taxonomy_Display, 
                                    levels = c(taxa_names_m, "Other"))

# H. Definir una paleta de colores. El Top 10 por especie puede generar
# fácilmente entre 15 y 30 familias únicas en total para graficar.
# Necesitamos una paleta grande.
n_colors_m <- length(unique(df_final_m$Taxonomy_Display))
getPalette_m = colorRampPalette(brewer.pal(12, "Paired"))
custom_colors_m <- getPalette(n_colors_m)

# Asignar explícitamente el color gris a "Otras"
names(custom_colors_m) <- levels(df_final_m$Taxonomy_Display)
custom_colors_m["Other"] <- "#D3D3D3" # Gris claro

# 4. Visualización con GGPLOT2
# ==============================================================================

grafico_taxonomia_ps_M <- ggplot(df_final_m, aes(x = especies, y = Percentage,
                                               fill = Taxonomy_Display)) +
  # Barras apiladas al 100% (position = "stack")
  geom_bar(stat = "identity", position = "stack", color = "black",
           linewidth = 0.2) +
  
  # Usar la paleta de colores personalizada
  scale_fill_manual(values = custom_colors_m) +
  
  # Etiquetas y formato
  labs(x = "Species",
       y = "Relative Abundance (%)",
       fill = "Family"
  ) +
  scale_x_discrete(labels=c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  # Tema limpio para tesis
  theme_classic() +
  theme(
    # Nombres de especies de Cucurbita en cursiva y rotados
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 9)
  )

# Mostrar gráfico
grafico_taxonomia_ps_M

# Uniendo graficas
ggarrange(grafico_taxonomia_ps_H, grafico_taxonomia_ps_M, ncol = 2, nrow = 1,
          labels = "AUTO")



###############################################################################
# ANALIZANDO UNICAMENTE LAS MUESTRAS DE FLORES VIRGENES (SIN VISITA) PARA
# EVALUAR EL EFECTO DE LA DOMESTICACIÓN SOBRE LA COMPOSICIÓN DE LA 
# COMUNIDAD DE MICROORGANISMOS.
#s16_7 = subset_samples(s16_6, tratamiento=="t1")
#s16_7
#summarize_phyloseq(s16_7)  #79771 reads, 1635 ASVs
s16_6

#### ANALYSIS WITH PROPORTION DATA ####
### Calculating relative abundances
mergedGP_pro <- merge_samples(s16_7, "especies")
# Use psmelt to obtain a long-format data.frame
phy_pro <- mergedGP_pro %>% tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
View(phy_pro)
#write.csv(phy_pro, file = "abundancias_family_propor_NV.csv")

#mergedtop = merge_samples(topfam, "especies")
#phy2 <- mergedtop %>% tax_glom(taxrank = "Family") %>%
#  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
#View(phy2)

#### HEATMAPs of most abundant bacteria ####
### Collapse dataset taxonomy to 10 most abundant families
topfam_pro <- aggregate_top_taxa2(s16_6, "Family", top = 10)
View(otu_table(topfam_pro))
topfam_pro

# transform families dataset to relative abundance
transform <- microbiome::transform
famrel_pro <- transform(topfam_pro, "compositional")
View(otu_table(famrel_pro))

### Heatmap for 10 most abundant bacterial families between species
heatmap_k_pro <- simple_heatmap(famrel_pro,
                                group.facet = "tratamiento",
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
write.csv(phy_hembras_pro, file = "abundancias_hembras_NV_pro.csv")

### Calculating relative abundances for Male flowers
merged_machos_pro <- merge_samples(machos_pro, "especies")
# Use psmelt to obtain a long-format data.frame
phy_machos_pro <- merged_machos_pro %>% tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
#View(phy_machos_pro)
write.csv(phy_machos_pro, file = "abundancias_machos_NV_pro.csv")


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

heatmaps_16s = ggarrange(heatmap_hempro, heatmap_macpro, nrow = 2,
                         labels=c("A", "B"))
heatmaps_16s



#### CORE MICROBIOTA ####
domipro <- dominant_taxa(s16_6,level = "Family")
# Most dominant families en todas las muestras
domipro_H <- dominant_taxa(hembras_pro,level = "Family") #group = "especies"
domipro_M <- dominant_taxa(machos_pro,level = "Family") #group = "especies"

#View(domipro$dominant_overview)
write.csv(domipro$dominant_overview, file = "dominant_overview.csv")
write.csv(domipro_H$dominant_overview, file = "dominant_overview_H.csv")
write.csv(domipro_M$dominant_overview, file = "dominant_overview_M.csv")

# PREVALENCE OF TAXONOMIC GROUPS global
ps.relpro = microbiome::transform(s16_6, "compositional")
# Relative population frequencies; at 1% compositional abundance threshold:
head(prevalence(ps.relpro, detection = 10/100, sort = TRUE))
# Absolute population frequencies (sample count)
head(prevalence(ps.relpro, detection = 1/100, sort = TRUE, count = TRUE))

# CORE MICROBIOME BY SEX
# Verificar cuántas muestras quedaron en cada grupo
cat("Muestras hembras:", nsamples(hembras_pro), "\n")
cat("Muestras machos: ", nsamples(machos_pro), "\n")

# Verificar distribución por especie dentro de cada sexo
cat("\nHembras por especie:\n")
table(sample_data(hembras_pro)$especie)

cat("\nMachos por especie:\n")
table(sample_data(machos_pro)$especie)

# Luego recalcular prevalencia por sexo
for(objeto in list(hembras_pro, machos_pro)){
  
  nombre <- ifelse(identical(objeto, hembras_pro), "HEMBRAS", "MACHOS")
  
  ps_rel <- microbiome::transform(objeto, "compositional")
  prev   <- prevalence(ps_rel, detection = 0.00001)
  
  cat("\n─────", nombre, "─────\n")
  cat("N muestras:", nsamples(objeto), "\n")
  cat("Taxa >50%:", sum(prev > 0.50), "\n")
  cat("Taxa >25%:", sum(prev > 0.25), "\n")
  cat("Taxa >10%:", sum(prev > 0.10), "\n")
  cat("Taxa >5%: ", sum(prev > 0.05), "\n")
}



# Verificar distribución de abundancias en tus datos
summary(otu_table(ps.relpro))

# ¿Cuántos taxa tienen abundancia > 0.00001 en al menos el 90% de muestras?
cat("Taxa en core con esos umbrales:", 
    ntaxa(core(ps.relpro, detection = 0.00001, prevalence = 0.90)), "\n")

# Ver cómo cambia el tamaño del core según la prevalencia
# con la abundancia fija en 0.001%
prevalencias <- c(0.50, 0.75, 0.90, 0.95, 0.99)

for(prev in prevalencias){
  n <- ntaxa(core(ps.relpro, 
                  detection  = 0.00001, 
                  prevalence = prev))
  cat("Prevalencia", prev*100, "%  →", n, "taxa en el core\n")
}

# Ver cómo se distribuyen las abundancias en tus datos
abundancias <- as.vector(otu_table(ps.relpro))
abundancias_nonzero <- abundancias[abundancias > 0]

summary(abundancias_nonzero)

# Histograma de abundancias
hist(log10(abundancias_nonzero),
     main = "Distribución de abundancias (log10)",
     xlab = "log10(Abundancia relativa)",
     col = "steelblue",
     breaks = 50)

# Ver la prevalencia real de los taxa más comunes
prev_real <- prevalence(ps.relpro, detection = 0.00001)
summary(prev_real)

# ¿Cuántos taxa tienen prevalencia > 50%?
cat("Taxa con prevalencia >50%:", sum(prev_real > 0.50), "\n")
cat("Taxa con prevalencia >25%:", sum(prev_real > 0.25), "\n")
cat("Taxa con prevalencia >10%:", sum(prev_real > 0.10), "\n")
cat("Taxa con prevalencia >5%: ", sum(prev_real > 0.05), "\n")

# Ver los 20 taxa más prevalentes
sort(prev_real, decreasing = TRUE)[1:20]



# Ver el core completo con ambos umbrales simultáneamente
plot_core(ps.relpro,
          prevalences = c(0.5, 0.75, 0.9, 0.99),
          detections = c(0.01, 0.05, 0.10, 0.20),
          plot.type = "lineplot")

nsamples(s16_6)

# Gráfica con umbrales biológicamente relevantes para 227 muestras
plot_core(ps.relpro,
          prevalences = c(0.10, 0.25, 0.50, 0.75, 0.90),
          detections  = c(0.01, 0.05, 0.10, 0.20),
          plot.type   = "lineplot") +
  xlab("Abundancia relativa") +
  ylab("Tamaño del core (N taxa)") +
  ggtitle("Core microbiome - Bacterias 16S")


# Con 227 muestras, convertir proporciones a conteos absolutos
# 0.10 * 227 = 23 muestras
# 0.25 * 227 = 57 muestras
# 0.50 * 227 = 114 muestras
# 0.75 * 227 = 170 muestras
# 0.90 * 227 = 204 muestras

plot_core(ps.relpro,
          prevalences = c(23, 57, 114, 170, 204),
          detections  = c(0.01, 0.05, 0.10, 0.20),
          plot.type   = "lineplot") +
  xlab("Abundancia relativa") +
  ylab("Tamaño del core (N taxa)") +
  ggtitle("Core microbiome - Bacterias 16S")


# Ver qué rango de prevalencias tienen tus datos actualmente
summary(prevalence(ps.relpro, detection = 0.01))

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

plot_core(ps.relpro, 
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
p <- plot_core(ps.relpro,
               plot.type = "heatmap", 
               colours = gray,
               prevalences = prevalences, 
               detections = detections, 
               min.prevalence = prevalence(ps.relpro, sort = TRUE)[100]) +
  labs(x = "Detection Threshold\n(Relative Abundance (%))") +
  
  #Adjusts axis text size and legend bar height
  theme(axis.text.y= element_text(size=8, face="italic"),
        axis.text.x.bottom=element_text(size=8),
        axis.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))

print(p)

# Plotting with upsetR
#remotes::install_github("YuLab-SMU/MicrobiotaProcess")
library(MicrobiotaProcess)
library(UpSetR)

# General
uppenpro <- get_upset(s16_6, factorNames = "especies")
View(uppenpro)
# función para gráficar los ASVs compartidos entre especies
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

# HEMBRAS
uppenpro_H <- get_upset(hembras_pro, factorNames = "especies")
View(uppenpro_H)
# función para gráficar los ASVs compartidos entre especies
uphembra = upset(uppenpro_H, nsets = 5, mainbar.y.label = "Number of ASVs shared",
      order.by = "freq",
      queries = list(list(query = intersects,
                          params = list("6CM", "5CAA", "4CPP"),
                          color = "orange",
                          active = T),
                     list(query = intersects,
                          params = list("6CM", "5CAA", "4CPP", 
                                        "3CAS", "2CPF"),
                          color = "green",
                          active = T),
                     list(query = intersects,
                         params = list("3CAS", "2CPF"),
                        color = "pink",
                       active = T)))


# MACHOS
uppenpro_M <- get_upset(machos_pro, factorNames = "especies")
View(uppenpro_M)
# función para gráficar los ASVs compartidos entre especies
upmacho = upset(uppenpro_M, nsets = 6, mainbar.y.label = "Number of ASVs shared",
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

# juntando gráficas
upsetgraf_16s <- grid.arrange(uphembra, upmacho, nrow = 2,
                             labels="AUTO")

upsetgraf_16s = ggarrange(heatmap_hem, heatmap_mac, nrow = 2,
                         labels=c("Hembras", "Machos"))
upsetgraf_16s

sha_pro <- ggarrange(uphembra, upmacho, ncol = 2, nrow = 1,
                     labels = "AUTO")
sha_pro

###############################################################################
###############################################################################
################################ ALPHA DIVERSITY  #############################
s16_6
# Primero nos aseguramos que estamos usando las funciones del paquete correcto
transform <- microbiome::transform
# Proportion Normalization
pseq16s_to <- transform(s16_6, "compositional")

# Cargamos librería
library(knitr)

# Obteniendo los valores de las metricas de diversidad alfa
div5 = microbiome::alpha(pseq16s_to, index = "all")
kable(head(div5))
View(div5)

# get the metadata out as separate object
ps2.meta2 <- meta(pseq16s_to)
kable(head(ps2.meta2))
# Add the rownames as a new colum for easy integration later.
ps2.meta2$sam_name2 <- rownames(ps2.meta2)

# Add the rownames to diversity table
div5$sam_name2 <- rownames(div5)

# merge these two data frames into one
div.df5 <- merge(div5,ps2.meta2, by = "sam_name2")

# check the tables
colnames(div.df5)

### Phylogenetic diversity (Faith)
library(picante)

ps0.rar.asvtab22 <- as.data.frame(pseq16s_to@otu_table)

ps0.rar.tree22 <- pseq16s_to@phy_tree

# hmp.meta from previous code chunks

# We first need to check if the tree is rooted or not 
pseq16s_to@phy_tree

# it is a rooted tree
df.pd5 <- pd(t(ps0.rar.asvtab22), ps0.rar.tree22,include.root=T) 
# t(ou_table) transposes the table for use in picante and the tre file comes
# from the first code chunck we used to read tree file (see making a phyloseq 
# object section).
#datatable(df.pd2)

# We add the results of PD to this file and then plot
div.df5$Phylogenetic_Diversity <- df.pd5$PD
head(div.df5)

### DIVERSIDAD ALPHA TOTAL
head(div.df5)
sort(div.df5$observed)
sort(div.df5$diversity_shannon)
sort(div.df5$Phylogenetic_Diversity)
#Histograma
hist(div.df5$observed)
hist(div.df5$diversity_shannon)
hist(div.df5$Phylogenetic_Diversity)
# #shapiro
shapiro.test(div.df5$observed)
shapiro.test(div.df5$diversity_shannon)
shapiro.test(div.df5$Phylogenetic_Diversity)

## STATISTICAL ANALYSIS OF ALPHA DIVERSITY ####
# Ver cuántas observaciones tienes por celda
div.df5 %>%
  group_by(condicion, tratamiento, sexo) %>%
  summarise(n = n()) %>%
  print(n = Inf)

# Resumen del desbalance por sexo
div.df5 %>%
  group_by(especies, sexo) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = sexo, 
                     values_from = n,
                     values_fill = 0)

# Excluyendo CF de machos
div.df6 <- div.df5 %>%
  filter(especies != "1CF")

# Este es tu modelo final
# modelo sin efecto aleatorio
modelo_total <- glm(diversity_shannon ~ condicion * sexo +
                       condicion * tratamiento +
                       sexo * tratamiento, 
                     data = div.df5, family = gaussian())

# Resultados
library(car)
Anova(modelo_total)

summary(modelo_total)     # resumen completo
# Verificar supuestos
plot(modelo_total)
qqnorm(resid(modelo_total))
qqline(resid(modelo_total))
shapiro.test(residuals(modelo_total))


# Modelo con efecto aleatorio = especies
library(lmerTest)
modelo_mixto <- lmer(diversity_shannon ~ condicion * sexo +
                       condicion * tratamiento +
                       sexo * tratamiento +
                       (1 | especies),
                     data = div.df5) # %>% filter(especies != "1CF")

# Resultados
anova(modelo_mixto) 

summary(modelo_mixto)     # resumen completo
# Verificar supuestos
plot(modelo_mixto)
qqnorm(resid(modelo_mixto))
qqline(resid(modelo_mixto))
shapiro.test(residuals(modelo_mixto))

# Post-hoc
library(emmeans)
emmeans(modelo_total, pairwise ~ condicion * sexo)

# Descomponer la interacción condicion × sexo
emm_condicion_sexo <- emmeans(modelo_total, 
                              pairwise ~ sexo | condicion)
summary(emm_condicion_sexo)

# También verlo al revés
emmeans(modelo_total, pairwise ~ condicion | sexo)

# Ver las medias estimadas
emmeans(modelo_total, ~ condicion * sexo)



#### Diversidad alpha en flores hembra ####
head(div.df5)
hembras_alfa_to = subset(div.df5, sexo == "hembra")
str(hembras_alfa_to)
head(hembras_alfa_to)

# data
sort(hembras_alfa_to$observed)
hembras_alfa_to$diversity_shannon
hembras_alfa_to$Phylogenetic_Diversity
# Histograma
hist(hembras_alfa_to$observed)
hist(hembras_alfa_to$diversity_shannon) 
hist(hembras_alfa_to$Phylogenetic_Diversity)
# Shapiro test
shapiro.test(hembras_alfa_to$observed)
shapiro.test(hembras_alfa_to$diversity_shannon)
shapiro.test(hembras_alfa_to$Phylogenetic_Diversity)

#### Diversidad alpha en flores macho ####
head(div.df5)
machos_alfa_to = subset(div.df5, sexo == "macho")
str(machos_alfa_to)

# data
sort(machos_alfa_to$observed)
machos_alfa_to$diversity_shannon
machos_alfa_to$Phylogenetic_Diversity
# Histograma
hist(machos_alfa_to$observed)
hist(machos_alfa_to$diversity_shannon) 
hist(machos_alfa_to$Phylogenetic_Diversity)
# Shapiro test
shapiro.test(machos_alfa_to$observed)
shapiro.test(machos_alfa_to$diversity_shannon)
shapiro.test(machos_alfa_to$Phylogenetic_Diversity)



# Flores Hembra
modglm3 <- glm(diversity_shannon ~ condicion / especies * tratamiento, 
               family = gaussian(link = "identity"), data = hembras_alfa_to)
summary(modglm3)
anova(modglm3)

modglm4 = glm(diversity_shannon ~ especies * tratamiento, 
            family = gaussian(link = "identity"), data = hembras_alfa_to)
summary(modglm4)
anova(modglm4)

library(emmeans)

# Asumiendo que tu modelo se llama modelo_shannon
# La barra "|" le dice a R: "Compara las Especies, pero hazlo POR SEPARADO 
# dentro de cada Estado"
# comparaciones con el modelo anidado
comp1 = emmeans(modglm1,
                pairwise ~ especies | condicion)
comp1
# comparando todas las especies sin anidar
comp2 <- emmeans(modglm1, pairwise ~ especies)
comp2

# hembras
comp3 <- emmeans(modglm3,
                 pairwise ~ especies | condicion)
comp3
# comparando todas las especies sin anidar
comp4 <- emmeans(modglm3, pairwise ~ especies)
comp4



##### VISUALIZATION  OF ALPHA DIVERSITY BY CONDITION####
sha_con = ggplot(div.df5, aes(x=condicion, y=diversity_shannon,
                                   fill=condicion)) +
  geom_boxplot(alpha = 0.8)+
  #geom_jitter(position=position_jitter(0.2), size = 2)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Condition", y = "Shannon diversity index")+
  scale_fill_discrete(name = "Condition", 
                      labels = c("domesticada" = "Domesticated",
                                 "silvestre" = "Wild"),
                      palette = c("#E6AB02","#4D9221"))+
  scale_x_discrete(labels=c("domesticada" = "Domesticated","silvestre" = "Wild"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  theme_classic()

sha_con

##### VISUALIZATION  OF ALPHA DIVERSITY BY SEX ####
sha_sex = ggplot(div.df5, aes(x=sexo, y=diversity_shannon,
                              fill=sexo)) +
  geom_boxplot(alpha = 0.8)+
  #geom_jitter(position=position_jitter(0.2), size = 2)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Floral sex", y = "Shannon diversity index")+
  scale_fill_discrete(name = "Floral sex", 
                      labels = c("hembra" = "Pistillate",
                                 "macho" = "Staminate"),
                      palette = c("#E6AB02","#4D9221"))+
  scale_x_discrete(labels=c("hembra" = "Pistillate","macho" = "Staminate"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  theme_classic()

sha_sex

##### VISUALIZATION  OF ALPHA DIVERSITY BY TREATMENT #####
library(ggsignif)

sha_tra <- ggplot(div.df5, aes(x = tratamiento, y = diversity_shannon,
                               fill = tratamiento)) +
  geom_boxplot(alpha = 0.8) +
  geom_signif(
    comparisons = list(c("t1", "t2")),  # grupos a comparar
    annotation = "**",                   # texto: *, **, ***, ns, o p-value
    tip_length = 0.02,                   # longitud de las patitas
    y_position = max(div.df5$diversity_shannon) * .9  # altura de la barra
  ) +
  scale_fill_manual(name = "Treatment",
                    values = c("t1" = "#E6AB02", "t2" = "#4D9221"),
                    labels = c("t1" = "Unvisited", "t2" = "Visited")) +
  scale_x_discrete(labels = c("t1" = "Unvisited", "t2" = "Visited")) +
  labs(x = "Treatment", y = "Shannon diversity index") +
  theme_classic()

sha_tra

# combinando plots
# cargando paqueteria gridExtra
library(gridExtra)
library(cowplot)

shato2 <- ggarrange(sha_con, sha_sex, sha_tra, 
                    ncol = 1, nrow = 3,
                     labels = "AUTO")
shato2

######## BOXPLOT CONDITION X SEX  ################
sha_conxsex <- ggplot(div.df5, aes(x = condicion, y = diversity_shannon,
                                   fill = sexo)) +
  geom_boxplot(alpha = 0.8) +
  
  # Asterisco SOLO en domesticadas (hembra vs macho)
  geom_signif(
    xmin = 0.8,          # borde izquierdo de la barra (caja hembra en domesticada)
    xmax = 1.2,          # borde derecho de la barra (caja macho en domesticada)
    y_position = max(div.df5$diversity_shannon, na.rm = TRUE) * 1.05,
    annotation = "**",   # p = 0.009
    tip_length = 0.02
  ) +
  
  scale_fill_manual(name = "Floral sex",
                    values = c("hembra" = "#E6AB02", "macho" = "#4D9221"),
                    labels = c("hembra" = "Pistillate",
                               "macho"  = "Staminate")) +
  scale_x_discrete(labels = c("domesticada" = "Domesticated",
                              "silvestre"   = "Wild")) +
  labs(x = "Condition", y = "Shannon diversity index") +
  theme_classic()

sha_conxsex


####### CONDITION X SPECIES ########
mod_spp = glm(diversity_shannon ~ especies,
              data = div.df5, family = gaussian)
Anova(mod_spp)
## Post-hoc
library(emmeans)
emmeans(mod_spp, pairwise ~ especies)

# Gráfica
# Dataframe con las letras y posiciones
# y_pos debe estar encima de cada boxplot
letras_df <- data.frame(
  especies = c("1CF", "2CPF", "3CAS", "4CPP", "5CAA", "6CM"),
  letra    = c("ab",  "ab",   "b",    "b",    "a",    "b"),
  # calcular posición y encima del máximo de cada especie
  y_pos = tapply(div.df5$diversity_shannon,
                 div.df5$especies,
                 function(x) max(x, na.rm = TRUE)) + 0.15
)

sha_spp <- ggplot(div.df5, aes(x = especies, y = diversity_shannon,
                               fill = condicion)) +
  geom_boxplot(alpha = 0.8) +
  
  # Letras de significancia
  geom_text(data = letras_df,
            aes(x = especies, y = y_pos, label = letra),
            inherit.aes = FALSE,
            size = 4) + #fontface = "bold"
  
  scale_fill_manual(name = "Condition",
                    values = c("domesticada" = "#E6AB02",
                               "silvestre"   = "#4D9221"),
                    labels = c("domesticada" = "Domesticated",
                               "silvestre"   = "Wild")) +
  scale_x_discrete(labels = c("1CF"  = "CF",
                              "2CPF" = "CPF",
                              "3CAS" = "CAS",
                              "4CPP" = "CPP",
                              "5CAA" = "CAA",
                              "6CM"  = "CM")) +
  labs(x = "Species",
       y = "Shannon diversity index") +
  theme_classic()
  #theme(axis.text.x = element_text(face = "italic"))

sha_spp


####### SEX X SPECIES ########
mod_sexXspp = glm(diversity_shannon ~ especies * sexo,
              data = div.df5, family = gaussian)
Anova(mod_sexXspp)
## Post-hoc
library(emmeans)
emmeans(mod_sexXspp, pairwise ~ sexo | especies)

# Gráfica
sha_sppxsex = ggplot(div.df5, aes(x=especies, y=diversity_shannon,
                              fill=sexo)) +
  geom_boxplot(alpha = 0.8)+
  #geom_jitter(position=position_jitter(0.2), size = 2)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Shannon diversity index")+
  scale_fill_discrete(name = "Floral sex", 
                      labels = c("hembra" = "Pistillate",
                                 "macho" = "Staminate"),
                      palette = c("#E6AB02","#4D9221"))+
  scale_x_discrete(labels=c("1CF"="CF","2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  theme_classic()

sha_sppxsex

####### TREATMENT X SPECIES ########
mod_traXspp = glm(diversity_shannon ~ especies * tratamiento,
                  data = div.df5, family = gaussian)
Anova(mod_traXspp)
## Post-hoc
library(emmeans)
emmeans(mod_traXspp, pairwise ~ tratamiento | especies)
# Gráfica
sha_sppxtreat = ggplot(div.df5, aes(x=especies, y=diversity_shannon,
                                  fill=tratamiento)) +
  geom_boxplot(alpha = 0.8)+
  #geom_jitter(position=position_jitter(0.2), size = 2)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Shannon diversity index")+
  scale_fill_discrete(name = "Treatment", 
                      labels = c("t1" = "Unvisited",
                                 "t2" = "Visited"),
                      palette = c("#E6AB02","#4D9221"))+
  scale_x_discrete(labels=c("1CF"="CF", "2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  theme_classic()

sha_sppxtreat

# Calcular posiciones y encima de cada especie
y_max_spp <- div.df5 %>%
  group_by(especies) %>%
  summarise(y_pos = max(diversity_shannon, na.rm = TRUE) + 0.15)

sha_sppxtreat <- ggplot(div.df5, aes(x = especies, y = diversity_shannon,
                                     fill = tratamiento)) +
  geom_boxplot(alpha = 0.8) +
  
  # Asterisco en 4CPP (significativa p = 0.024)
  # 4CPP es la 4ta especie → posición x = 4
  # t1 (unvisited) = 4 - 0.2 = 3.8
  # t2 (visited)   = 4 + 0.2 = 4.2
  geom_signif(
    xmin       = 3.8,
    xmax       = 4.2,
    y_position = y_max_spp$y_pos[y_max_spp$especies == "4CPP"],
    annotation = "*",
    tip_length = 0.02,
    size       = 0.5,
    textsize   = 5
  ) +
  
  scale_fill_manual(name   = "Treatment",
                    values = c("t1" = "#E6AB02", "t2" = "#4D9221"),
                    labels = c("t1" = "Unvisited", "t2" = "Visited")) +
  scale_x_discrete(labels = c("1CF"  = "CF",
                              "2CPF" = "CPF",
                              "3CAS" = "CAS",
                              "4CPP" = "CPP",
                              "5CAA" = "CAA",
                              "6CM"  = "CM")) +
  labs(x       = "Species",
       y       = "Shannon diversity index") +
  theme_classic() +
  theme(plot.caption = element_text(hjust = 0, size = 9))

sha_sppxtreat



### Shannon flores macho
# Crear el gráfico base
grafico_shannon1 <- ggplot(machos_alfa_to, aes(x = especies,
                          y = diversity_shannon, fill = condicion)) +
  #, fill = tratamiento
  # Añadir las cajas (boxplots)
  geom_boxplot(color = "black", alpha = 0.8, outlier.shape = 16) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1),
             alpha = 0.4, color = "black") +
  # Dividir el gráfico en dos paneles (Silvestre y Domesticada)
  # scales = "free_x" y space = "free_x" son CLAVES para diseños anidados
  facet_grid(. ~ condicion, scales = "free_x", space = "free_x") +
  
  # Personalizar los colores de Visitadas y No visitadas (puedes cambiar los colores)
  scale_fill_manual(values = c("t2" = "#4DAF4A", "t1" = "#377EB8")) +
  
  # Añadir etiquetas claras
  labs(x = "Species",
    y = "Shannon Diversity Index",
    fill = "Treatment:"
  ) +
  
  # Usar un tema limpio y minimalista, ideal para publicaciones o tesis
  theme_bw() +
  theme(
    # Rotar un poco los nombres de las especies si son muy largos
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "italic"),
    axis.title = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"), # Texto de los paneles (Estado)
    strip.background = element_rect(fill = "lightgray"),
    legend.position = "bottom"
  )

# Mostrar el gráfico
grafico_shannon1

### Shannon flores macho
# Crear el gráfico base
grafico_shannon2 <- ggplot(hembras_alfa_to, aes(x = especies,
                           y = diversity_shannon, fill = tratamiento)) +
  
  # Añadir las cajas (boxplots)
  geom_boxplot(color = "black", alpha = 0.8, outlier.shape = 16) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1),
             alpha = 0.4, color = "black") +
  # Dividir el gráfico en dos paneles (Silvestre y Domesticada)
  # scales = "free_x" y space = "free_x" son CLAVES para diseños anidados
  facet_grid(. ~ condicion, scales = "free_x", space = "free_x") +
  
  # Personalizar los colores de Visitadas y No visitadas (puedes cambiar los colores)
  scale_fill_manual(values = c("t2" = "#4DAF4A", "t1" = "#377EB8")) +
  
  # Añadir etiquetas claras
  labs(x = "Species",
       y = "Shannon Diversity Index",
       fill = "Treatment:"
  ) +
  
  # Usar un tema limpio y minimalista, ideal para publicaciones o tesis
  theme_bw() +
  theme(
    # Rotar un poco los nombres de las especies si son muy largos
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "italic"),
    axis.title = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"), # Texto de los paneles (Estado)
    strip.background = element_rect(fill = "lightgray"),
    legend.position = "bottom"
  )

# Mostrar el gráfico
print(grafico_shannon2)










# Stripchart of observed ASVs of female flowers
obs_hem <- ggplot(hembras_alfa, aes(x=especies, y=observed, color=condicion)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 2)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Richness")+
  scale_color_discrete(name = "Condition", 
                       labels = c("domesticada" = "Domesticated",
                                  "silvestre" = "Wild"))+
  scale_x_discrete(labels=c("2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  #theme_bw()
  theme_classic()
obs_hem

head(hembras_alfa)
#Exportamos grafico de alfa div (observed) a pdf
ggsave("observed_hembras.pdf", 
       plot = obs_hem,
       width = 25, 
       height = 20, 
       units = "cm")

# Stripchart of observed ASVs of male flowers
obs_machos <- ggplot(machos_alfa, aes(x=especies, y=observed, color=condicion)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 2)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Richness")+
  scale_color_discrete(name = "Condition", 
                       labels = c("domesticada" = "Domesticated",
                                  "silvestre" = "Wild"))+
  scale_x_discrete(labels=c("1CF" = "CF", "2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  #theme_bw()
  theme_classic()
obs_machos
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

#### stripchart de Shannon diversity index in female flowers
sha_hem = ggplot(hembras_alfa, aes(x=especies, y=diversity_shannon,
                                   fill=condicion)) +
  geom_boxplot(alpha = 0.8)+
  geom_jitter(position=position_jitter(0.2), size = 2)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Shannon diversity index")+
  scale_fill_discrete(name = "Condition", 
                      labels = c("domesticada" = "Domesticated",
                                 "silvestre" = "Wild"),
                      palette = c("#E6AB02","#4D9221"))+
  scale_x_discrete(labels=c("2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  theme_classic()

sha_hem    

#display.brewer.pal(n = 8, name = 'Greens')
#brewer.pal(n = 8, name = "PiYG")

### stripchart de Shannon diversity index in male flowers
sha_mac = ggplot(machos_alfa, aes(x=especies, y=diversity_shannon,
                                  fill=condicion)) +
  geom_boxplot(alpha = 0.8)+
  geom_jitter(position=position_jitter(0.2), size = 2)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Species", y = "Shannon diversity index")+
  scale_fill_discrete(name = "Condition", 
                      labels = c("domesticada" = "Domesticated",
                                 "silvestre" = "Wild"),
                      palette = c("#E6AB02","#4D9221"))+
  scale_x_discrete(labels=c("1CF" = "CF", "2CPF" = "CPF","3CAS" = "CAS",
                            "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  #scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
  #                                                       "Visitadas")) +
  theme_classic() +
  theme(legend.position = "none")
sha_mac

#facet_wrap(~sexo, labeller = labeller(sexo = sexos))
#stat_summary(fun=mean, geom="point", shape=18,
#            size=6, color="red")
sha_pro <- ggarrange(sha_hem, sha_mac, ncol = 2, nrow = 1,
                     labels = "AUTO") 
# exporting plots
ggsave("alpha_shan_mar26.pdf", 
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



###############################################################################
#### DIVERSIDAD BETA ####
########## BETA DIVERSITY  WITH PROPORTION DATA ##########
library(vegan)
# phyloseq object filtered but no rarefied
s16_6 

# divided by sex
#hembras_pro
#machos_pro

#### Unweighted UniFrac ####
# Primero nos aseguramos que estamos usando las funciones del paquete correcto
transform <- microbiome::transform
# Proportion Normalization
pseq_16s <- transform(s16_6, "compositional")


# PCoA plot using the unweighted UniFrac as distance
ordination_16s = ordinate(pseq_16s, method="PCoA", "unifrac",
                              weighted = F)

baseplot_16s <- plot_ordination(pseq_16s, ordination_16s,
                                    color="condicion",
                                    #shape = "condicion",
                                    type = "samples",
) + #shape = "sexo"
  ggtitle("")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 5)


  scale_color_discrete(name = "Species", 
                       labels = c("2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  scale_shape_discrete(name = "Condition",
                       labels = c("domesticada" = "Domesticated",
                                  "silvestre" = "Wild"))
#stat_ellipse()
baseplot_16s

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
                                    shape = "condicion",
                                    type = "samples",
) + #shape = "sexo"
  ggtitle("")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 5)+
  scale_color_discrete(name = "Species", 
                       labels = c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))+
  scale_shape_discrete(name = "Condition",
                       labels = c("domesticada" = "Domesticated",
                                  "silvestre" = "Wild"))
#stat_ellipse()
baseplot_16s_mac


#### Weighted UniFrac GLOBAL####
# PCoA plot using the weighted UniFrac as distance
ordination_16s_wei = ordinate(pseq_16s,
                                  method="PCoA", "wunifrac")
# Plot
baseplot_16s_wei <- plot_ordination(pseq_16s, 
                                        ordination_16s_wei,
                                        color="condicion",
                                        #shape = "condicion",
                                        type = "samples") + #shape = "sexo" shape = "tratamiento"
  ggtitle("")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 4)+
  scale_color_discrete(name = "Condition",
                       labels = c("domesticada" = "Domesticated",
                                  "silvestre" = "Wild"),
                       palette = c("#E6AB02","#4D9221"))
  #scale_color_discrete(name = "Species", 
   #                    labels = c("2CPF" = "CPF","3CAS" = "CAS",
    #                              "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"),
               #        palette = "Set1")+
  #scale_shape_discrete(name = "Condition",
   #                    labels = c("domesticada" = "Domesticated",
    #                              "silvestre" = "Wild"))
#scale_shape_manual(values=c(16,3))
#labs(x = "PCoa1 (31.3%)", y = "PCoA2 (10.4%)")+
#scale_color_discrete(name = "Condition", 
#                    labels = c("Domesticated","Wild"))+
#scale_color_discrete(name = "Treatment",
#                     labels = c("No-visited","Visited"))

#baseplot_16s_wei
#baseplot_16s_hem_wei + geom_polygon(aes(fill = especies))
baseplot_16s_wei + stat_ellipse(linetype = 2)
#
#baseplot2 +
#  geom_point(size = 1) +
#  geom_label(aes(label = clave), size = 3.5, color = "black")
#
#wei_16s = baseplot_16s_hem_wei + stat_ellipse(type = "norm", linetype = 2)
#wei_16s
#+ coord_fixed()

#+ stat_ellipse(type = "t")
#baseplot_16s_hem_wei + stat_ellipse(type = "t")

# Plot
baseplot_16s_wei2 <- plot_ordination(pseq_16s, 
                                    ordination_16s_wei,
                                    color="tratamiento",
                                    #shape = "condicion",
                                    type = "samples") + 
  ggtitle("")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 4)+
  scale_color_discrete(name = "Treatment",
                       labels = c("t1" = "Unvisited", "t2" = "Visited"),
                       palette = c("#E6AB02","#4D9221"))

baseplot_16s_wei2 + stat_ellipse(linetype = 2)




#### Weighted UniFrac MACHOS####
# PCoA plot using the weighted UniFrac as distance
ordination_16s_mac_wei = ordinate(pseq_16s_mac,
                                  method="PCoA", "wunifrac")
# Plot
baseplot_16s_mac_wei <- plot_ordination(pseq_16s_mac, 
                                        ordination_16s_mac_wei,
                                        color="especies",
                                        shape = "condicion",
                                        type = "samples") + #shape = "sexo" shape = "tratamiento"
  ggtitle("")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 4)+
  scale_color_discrete(name = "Species", 
                       labels = c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                                  "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"),
                       palette = "Set1")+
  scale_shape_discrete(name = "Condition",
                       labels = c("domesticada" = "Domesticated",
                                  "silvestre" = "Wild"))
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


################ Run PERMANOVA global ###################
transform <- microbiome::transform
# Proportion Normalization
pseq_16s <- transform(s16_6, "compositional")
metadf = data.frame(sample_data(s16_6))

# calculating Weighted-UniFrac
wunifrac.dist <- UniFrac(pseq_16s, 
                             weighted = TRUE, 
                             normalized = TRUE,  
                             parallel = FALSE, 
                             fast = TRUE)

set.seed(100)   # Set the seed to make the analysis reproducible
perm = adonis2(wunifrac.dist ~ condicion * sexo * tratamiento,
                   data = metadf, by = "terms", permutations = 999)
perm

# Condicion (el más importante dado que es significativo en beta)
bd_condicion <- betadisper(wunifrac.dist, metadf$condicion)
permutest(bd_condicion, permutations = 999)

# Tratamiento
bd_trat <- betadisper(wunifrac.dist, metadf$tratamiento)
permutest(bd_trat, permutations = 999)


###############################################################################
#### PERMANOVA by floral gender####
## female flowers
metadf_hem = data.frame(sample_data(hembras_pro))

# calculating Weighted-UniFrac
wunifrac.dist_hem <- UniFrac(pseq_16s_hem, 
                             weighted = TRUE, 
                             normalized = TRUE,  
                             parallel = FALSE, 
                             fast = TRUE)

# Run PERMANOVA
set.seed(100)   # Set the seed to make the analysis reproducible
perm_hem = adonis2(wunifrac.dist_hem ~ condicion / especies * tratamiento,
                   data = metadf_hem, by = "terms", permutations = 999)
perm_hem


#Test dispersion among treatments and species
disp_wun1 <- betadisper(wunifrac.dist_hem, sample_data(hembras_pro)$condicion)
permutest(disp_wun1)
disp_wun2 <- betadisper(wunifrac.dist_hem, sample_data(hembras_pro)$tratamiento)
permutest(disp_wun2)
disp_wun3 <- betadisper(wunifrac.dist_hem, sample_data(hembras_pro)$especies)
permutest(disp_wun3) # especies heterogeneas


#### Male flowers ####
metadf_mac = data.frame(sample_data(machos_pro))

# calculating Weighted-UniFrac
wunifrac.dist_mac <- UniFrac(pseq_16s_mac, 
                             weighted = TRUE, 
                             normalized = TRUE,  
                             parallel = FALSE, 
                             fast = TRUE)

# Run PERMANOVA of MALE FLOWERS
set.seed(100)   # Set the seed to make the analysis reproducible
perm_mac = adonis2(wunifrac.dist_mac ~ condicion / especies * tratamiento,
                   data = metadf_mac, by = "terms", permutations = 999)
perm_mac



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
# A. Aglomerar a nivel de Familia
ps_fam <- tax_glom(s16_6, taxrank = "Family", NArm = FALSE)

# B. Transformar a abundancia relativa
ps_rel <- transform_sample_counts(ps_fam, function(x) x / sum(x))

# C. Psmelt - cada fila es una muestra × familia
df_melt <- psmelt(ps_rel) %>%
  mutate(Family = ifelse(is.na(Family), "Unclassified", 
                         as.character(Family)))

# D. Top 10 familias más abundantes globalmente
top10_familias <- df_melt %>%
  group_by(Family) %>%
  summarise(AbundMedia = mean(Abundance)) %>%
  slice_max(order_by = AbundMedia, n = 10) %>%
  pull(Family)

# E. Colapsar el resto en "Other" por muestra
df_final <- df_melt %>%
  mutate(Taxonomy_Display = ifelse(Family %in% top10_familias,
                                   Family, "Other")) %>%
  group_by(Sample, especies, Taxonomy_Display) %>%
  summarise(Percentage = sum(Abundance) * 100, .groups = "drop")

# F. Ordenar factores
taxa_names <- sort(setdiff(unique(df_final$Taxonomy_Display), "Other"))
df_final$Taxonomy_Display <- factor(df_final$Taxonomy_Display,
                                    levels = c(taxa_names, "Other"))

# G. Ordenar muestras por especie para que queden agrupadas
df_final <- df_final %>%
  arrange(especies, Sample)
df_final$Sample <- factor(df_final$Sample, 
                          levels = unique(df_final$Sample))

# H. Paleta de colores
n_colors <- length(unique(df_final$Taxonomy_Display))
getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
custom_colors <- getPalette(n_colors)
names(custom_colors) <- levels(df_final$Taxonomy_Display)
custom_colors["Other"] <- "#D3D3D3"

# I. Gráfica por muestra
grafico_muestras <- ggplot(df_final,
                           aes(x = Sample,
                               y = Percentage,
                               fill = Taxonomy_Display)) +
  geom_bar(stat = "identity", position = "stack",
           color = NA) +  # sin borde para 227 muestras
  scale_fill_manual(values = custom_colors) +
  labs(x    = "Sample",
       y    = "Relative Abundance (%)",
       fill = "Family") +
  theme_classic() +
  theme(axis.text.x   = element_blank(),  # ocultar nombres de muestras
        axis.ticks.x  = element_blank(),  # son 227, no caben
        axis.title    = element_text(size = 12, face = "bold"),
        legend.text   = element_text(size = 9))

print(grafico_muestras)

# Calcular posiciones de corte entre especies
# Agregar columna con posición numérica de cada muestra
df_final <- df_final %>%
  mutate(sample_num = as.numeric(factor(Sample, 
                                        levels = levels(Sample))))

# Ahora calcular cortes usando sample_num
cortes <- df_final %>%
  group_by(especies) %>%
  summarise(
    inicio = min(sample_num),
    fin    = max(sample_num),
    centro = (min(sample_num) + max(sample_num)) / 2
  )

print(cortes)

###
grafico_muestras <- ggplot(df_final,
                           aes(x = Sample,
                               y = Percentage,
                               fill = Taxonomy_Display)) +
  geom_bar(stat = "identity", position = "stack",
           color = NA) +
  scale_fill_manual(values = custom_colors) +
  
  # Líneas divisorias entre especies
  geom_vline(xintercept = c(4.5, 15.5, 73.5, 109.5, 171.5),
             linetype = "dashed",
             color    = "gray40",
             linewidth = 0.5) +
  
  # Etiquetas de especie debajo del eje x
  annotate("text",
           x        = cortes$centro,
           y        = -5,
           label    = c("CF", "CPF", "CAS", "CPP", "CAA", "CM"),
           size     = 3.5,
           fontface = "italic") +
  
  labs(x    = "Samples",
       y    = "Relative Abundance (%)",
       fill = "Family") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.02))) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title   = element_text(size = 12, face = "bold"),
        legend.text  = element_text(size = 9),
        plot.margin  = margin(t = 5, r = 5, b = 30, l = 5))

print(grafico_muestras)



#######################33#####
#### Obteniendo las abundancias relativas por especie
# 1. Preparación del objeto Phyloseq
# A. Aglomerar a nivel de Familia. Esto suma los reads de todos los ASVs/OTUs
# que pertenecen a la misma familia. Reduce la complejidad drásticamente.
ps_fam <- tax_glom(s16_6, taxrank = "Family", NArm = FALSE)

# B. Transformar a abundancia relativa (escala 0-1)
# Esto estandariza las muestras independientemente de su profundidad de secuenciación.
ps_rel <- transform_sample_counts(ps_fam, function(x) x / sum(x))

# C. Convertir el objeto phyloseq a un dataframe largo ('melted' dataframe)
df_melt <- psmelt(ps_rel)

# D. Agrupar por Especie y Familia, y calcular la abundancia relativa PROMEDIO
# para esa especie de planta.
df_grouped <- df_melt %>%
  dplyr::group_by(condicion, Family) %>% 
  dplyr::summarize(MeanRelAbund = mean(Abundance), .groups = "drop")

# E. Identificar el Top 10 de familias más abundantes POR ESPECIE.
# Usamos slice_max para quedarnos con las 10 filas con mayor abundancia promedio por grupo.
top10_per_species <- df_grouped %>%
  group_by(condicion) %>%
  slice_max(order_by = MeanRelAbund, n = 10, with_ties = FALSE) %>%
  mutate(Classification = Family) %>% # Guardamos el nombre para el cruce
  ungroup()

# F. Crear la categoría "Otras".
# Cruzamos la tabla agrupada completa con la tabla del Top 10.
df_final <- df_grouped %>%
  left_join(select(top10_per_species, condicion, Family, Classification), 
            by = c("condicion", "Family")) %>%
  # Si la familia NO estaba en el top10 (Classification es NA), la llamamos "Otras"
  dplyr::mutate(Taxonomy_Display = ifelse(is.na(Classification), 
                                          "Other", 
                                          as.character(Family))) %>%
  # Volvemos a agrupar y sumar para colapsar todas las "Otras" en una sola fila por especie
  dplyr::group_by(condicion, Taxonomy_Display) %>%
  dplyr::summarize(TotalMeanAbund = sum(MeanRelAbund), .groups = 'drop') %>%
  # Convertir abundancia relativa (0-1) a porcentaje (0-100) para el gráfico
  dplyr::mutate(Percentage = TotalMeanAbund * 100)
View(df_final)
#write.csv(df_final, file = "abundancias_family_hembras_marzo2026_10.csv")

# 3. Preparación para Graficar
# ==============================================================================

# G. Ordenar los factores para que "Otras" siempre aparezca abajo en la barra y leyenda.
# Primero obtenemos todos los nombres únicos excepto "Otras" y los ordenamos alfabéticamente
taxa_names <- sort(setdiff(unique(df_final$Taxonomy_Display), "Other"))
# Creamos el factor con "Otras" al final
df_final$Taxonomy_Display <- factor(df_final$Taxonomy_Display, 
                                    levels = c(taxa_names, "Other"))

# H. Definir una paleta de colores. El Top 10 por especie puede generar
# fácilmente entre 15 y 30 familias únicas en total para graficar.
# Necesitamos una paleta grande.
n_colors <- length(unique(df_final$Taxonomy_Display))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
custom_colors <- getPalette(n_colors)

# Asignar explícitamente el color gris a "Otras"
names(custom_colors) <- levels(df_final$Taxonomy_Display)
custom_colors["Other"] <- "#D3D3D3" # Gris claro

# 4. Visualización con GGPLOT2
# ==============================================================================

grafico_taxonomia <- ggplot(df_final, aes(x = condicion, y = Percentage,
                                               fill = Taxonomy_Display)) +
  # Barras apiladas al 100% (position = "stack")
  geom_bar(stat = "identity", position = "stack", color = "black",
           linewidth = 0.2) +
  
  # Usar la paleta de colores personalizada
  scale_fill_manual(values = custom_colors) +
  
  # Etiquetas y formato
  labs(x = "Condition",
       y = "Relative Abundance (%)",
       fill = "Family"
  ) +
  scale_x_discrete(labels=c("domesticada" = "Domesticated",
                            "silvestre" = "Wild"))+
  # Tema limpio para tesis
  theme_classic() +
  theme(
    # Nombres de especies de Cucurbita en cursiva y rotados
   # axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 11),
    #axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 9)
  )

# Mostrar gráfico
grafico_taxonomia

###############################################
install.packages("tinytex")
tinytex::install_tinytex()
tinytex::install_tinytex()



