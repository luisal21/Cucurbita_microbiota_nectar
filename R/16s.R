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
#head(asvtable, 1)
#asv2 = t(asvtable)
#head(asv2, 1)
#write.csv(asv2, file = "asv_table.csv")

# cargar primero el archivo de metadatos
meta_16s <- read.csv("meta_16s_2024.csv", header = T)
head(meta_16s)
# despues cargar la tabla de ASV
seqtab.16s <- read.csv("asv_table.csv", header = T, row.names = 1)
head(seqtab.16s)
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
#arbol_16s = rtree(ntaxa(s16_1), rooted=TRUE, tip.label=taxa_names(s16_1))
#plot(arbol_16s)
# Uniendo el arbol con el objeto phyloseq anterior
#s16_2 = merge_phyloseq(s16_1, arbol_16s)
s16_2 = s16_1
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
#View(contamdf.freq)
table(contamdf.freq$contaminant)
head(which(contamdf.freq$contaminant))

# Identify Contaminants -Prevalence
sample_data(s16_2)$is.neg <- sample_data(s16_2)$sample_or_control == "Control Sample"
contamdf.prev <- isContaminant(s16_2, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant), 20)
#View(contamdf.prev)

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


#### CONTAR TAXONES POR NIVEL TAXONÓMICO #####
# Extraer la tabla de taxonomía
taxo <- as.data.frame(tax_table(s16_6))

# Contar taxones únicos por nivel
n_species <- length(unique(na.omit(taxo$Species)))
n_genera  <- length(unique(na.omit(taxo$Genus)))
n_families <- length(unique(na.omit(taxo$Family)))
n_orders  <- length(unique(na.omit(taxo$Order)))
n_classes <- length(unique(na.omit(taxo$Class)))
n_phyla   <- length(unique(na.omit(taxo$Phylum)))

# Ver todos juntos
data.frame(
  Nivel   = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
  N_taxones = c(n_phyla, n_classes, n_orders, n_families, n_genera, n_species)
)

###############################################################################
###############################################################################
### ANÁLISIS DE NÚMEROS DE HILL

#install.packages("iNEXT")       # rarefacción y extrapolación
#install.packages("hillR")       # Hill a partir de phyloseq

library(hillR)
library(phyloseq)

# Extraer tabla de OTU/ASV (muestras en filas)
asv_table <- as.data.frame(t(otu_table(s16_6)))

# Calcular Hill para cada q
hill_q0 <- hill_taxa(asv_table, q = 0)  # riqueza
hill_q1 <- hill_taxa(asv_table, q = 1)  # Shannon
hill_q2 <- hill_taxa(asv_table, q = 2)  # Simpson

# Unir con metadatos
hill_df <- data.frame(
  sample_id  = names(hill_q0),
  q0         = hill_q0,
  q1         = hill_q1,
  q2         = hill_q2
)

# Agregar variables del experimento
meta <- as.data.frame(sample_data(s16_6))
hill_df <- merge(hill_df, meta, by = "row.names")
View(hill_df)





###   CURVAS DE RAREFACCIÓN
### POR CONDICIÓN
library(iNEXT)

# Extraer tabla de ASVs y metadatos del objeto phyloseq
asv_table <- as.data.frame(otu_table(s16_6))# ASVs en filas, muestras en columnas
metadata <- as.data.frame(as.matrix(sample_data(s16_6)))

# Ver los niveles de condición
unique(metadata$condicion)

# Separar muestras por condición y sumar abundancias
muestras_dom <- metadata %>% 
  filter(condicion == "domesticada") %>% 
  rownames()

muestras_sil <- metadata %>% 
  filter(condicion == "silvestre") %>% 
  rownames()

# Sumar todas las muestras de cada condición (un vector por condición)
asv_domesticada <- rowSums(asv_table[, muestras_dom])
asv_silvestre   <- rowSums(asv_table[, muestras_sil])

# Crear lista para iNEXT (eliminar ASVs con abundancia 0)
asv_por_condicion <- list(
  domesticada = asv_domesticada[asv_domesticada > 0],
  silvestre   = asv_silvestre[asv_silvestre > 0]
)

# Verificar
str(asv_por_condicion)

# Números de Hill por condición
out_condicion <- iNEXT(asv_por_condicion,
                       q        = c(0, 1, 2),
                       datatype = "abundance",
                       nboot    = 100,
                       knots    = 40)
out_condicion$DataInfo
out_condicion$AsyEst

write.csv(out_condicion$DataInfo, file = "Hill_condicion_datainfo.csv")
write.csv(out_condicion$AsyEst, file = "Hill_condicion_asyest.csv")

coverage_incidence_con <- estimateD(asv_por_condicion,
                                q = c(0, 1, 2),
                                datatype = "abundance",
                                base = "coverage",
                                conf = 0.95, nboot = 100)
coverage_incidence_con
unique(coverage_incidence_con$SC)

### Estimando los numeros de Hill con la función estimateD por cobertura
coverage_incidence <- estimateD(asv_table,
                                q = c(0, 1, 2),
                                datatype = "abundance",
                                base = "coverage",
                                conf = 0.95, nboot = 100)
coverage_incidence
unique(coverage_incidence$SC)
View(coverage_incidence)

# Uniendo tabla de estimateD con metadatos
library(dplyr)
library(tidyr)
# Uniendo tablas
hill_wide_ic <- coverage_incidence %>%
  filter(Method == "Observed") %>%
  dplyr::select(Assemblage, Order.q, qD, qD.LCL, qD.UCL) %>%
  pivot_wider(names_from  = Order.q,
              values_from = c(qD, qD.LCL, qD.UCL),
              names_sep   = "_q") %>%    # genera qD_q0, qD.LCL_q0, etc.
  rename(sampleid2 = Assemblage) %>%
  left_join(metadata, by = "sampleid2")
# Verificando
colnames(hill_wide_ic)
head(hill_wide_ic)

# Paso 1: histogramas de Hill numbers
hist(hill_wide_ic$qD_q0)
hist(hill_wide_ic$qD_q1)
hist(hill_wide_ic$qD_q2)

shapiro.test(hill_wide_ic$qD_q0)
shapiro.test(hill_wide_ic$qD_q1)
shapiro.test(hill_wide_ic$qD_q2)

# Paso 2: Análisis estadístico con tus variables
# Para cada orden de Hill
# Modelo con q0
mod_q0 <- glm(qD_q0 ~ condicion * sexo * tratamiento, 
              data = hill_wide_ic,
              family = poisson())
summary(mod_q0)
# Calcular ratio de dispersión (debe ser ~1 si no hay sobredispersión)
E2 <- resid(mod_q0, type = "pearson")
sum(E2^2) / df.residual(mod_q0)
# Modelo con q0 binomial negativo
library(MASS)
# q0 modelo bueno
modnb_q0 <- glm.nb(qD_q0 ~ condicion * sexo * tratamiento, 
                 data = hill_wide_ic,
                 link = "log")
summary(modnb_q0)
# Modelo con q1
mod_q1 <- glm(qD_q1 ~ condicion * sexo * tratamiento, 
              data = hill_wide_ic,
              family = Gamma(link = "log"))
summary(mod_q1)
# Modelo con q2
mod_q2 <- glm(qD_q2 ~ condicion * sexo * tratamiento, 
              data = hill_wide_ic,
              family = Gamma(link = "log"))
summary(mod_q2)

# ANOVA tipo II
library(car)
Anova(modnb_q0)
Anova(mod_q1)
Anova(mod_q2)

# Summary
#q0
tapply(hill_wide_ic$qD_q0, list(hill_wide_ic$condicion,hill_wide_ic$sexo,
                                hill_wide_ic$tratamiento), mean)
tapply(hill_wide_ic$qD_q0, list(hill_wide_ic$condicion,hill_wide_ic$sexo,
                        hill_wide_ic$tratamiento), sd)
#q1
tapply(hill_wide_ic$qD_q1, list(hill_wide_ic$condicion,hill_wide_ic$sexo,
                        hill_wide_ic$tratamiento), mean)
tapply(hill_wide_ic$qD_q1, list(hill_wide_ic$condicion,hill_wide_ic$sexo,
                        hill_wide_ic$tratamiento), sd)
#q2
tapply(hill_wide_ic$qD_q2, list(hill_wide_ic$condicion,hill_wide_ic$sexo,
                        hill_wide_ic$tratamiento), mean)
tapply(hill_wide_ic$qD_q2, list(hill_wide_ic$condicion,hill_wide_ic$sexo,
                        hill_wide_ic$tratamiento), sd)


# Gráfica de curvas de cobertura de muestra
# Ver qué tan completo está cada muestra/grupo
# ¿Mis muestras están bien secuenciadas? (cobertura)
ggiNEXT(out_condicion, type = 2)+
  scale_color_manual(name = "Condition",
                     values = c("domesticada" = "#E6AB02", "silvestre" = "#4D9221"),
                     labels = c("domesticada" = "Domesticated", "silvestre" = "Wild"))+
  scale_fill_manual(name = "Condition",
                    values = c("domesticada" = "#E6AB02", "silvestre" = "#4D9221"),
                    labels = c("domesticada" = "Domesticated", "silvestre" = "Wild"))+
  guides(
    shape    = "none"   # ← quita la leyenda negra de formas
    #linetype = "none"    # ← quita la leyenda de líneas si aparece duplicada
  ) + 
  theme_classic()+
  #theme_bw()
  labs(#title  = "Curvas de rarefacción por condición",
    x      = "Number of sequences")
#y      = "Diversidad (Hill)",
#color  = "Condición",
#fill   = "Condición")

# Curvas de diversidad por cobertura
# Comparar diversidad a igual completitud  
# ¿Difieren los grupos a igual completitud? (comparación justa)  
ggiNEXT(out_condicion, type = 3, facet.var = "Order.q")+
  scale_color_manual(name = "Condition",
                     values = c("domesticada" = "#E6AB02", "silvestre" = "#4D9221"),
                     labels = c("domesticada" = "Domesticated", "silvestre" = "Wild"))+
  scale_fill_manual(name = "Condition",
                    values = c("domesticada" = "#E6AB02", "silvestre" = "#4D9221"),
                    labels = c("domesticada" = "Domesticated", "silvestre" = "Wild"))+
  guides(
    shape    = "none"   # ← quita la leyenda negra de formas
    #linetype = "none"    # ← quita la leyenda de líneas si aparece duplicada
  )+
  theme_classic()
#facet_wrap(~Order.q, scales = "free_y",
#           labeller = labeller(Order.q = c("0" = "q = 0  (Riqueza)",
#                                          "1" = "q = 1  (Shannon)",
#                                         "2" = "q = 2  (Simpson)")))

### CURVAS DE RAREFACCIÓN POR SEXO
# Separar muestras por sexo y sumar abundancias
muestras_hem <- metadata %>% 
  filter(sexo == "hembra") %>% 
  rownames()

muestras_mac <- metadata %>% 
  filter(sexo == "macho") %>% 
  rownames()

# Sumar todas las muestras de cada sexo (un vector por sexo)
asv_hem <- rowSums(asv_table[, muestras_hem])
asv_mac   <- rowSums(asv_table[, muestras_mac])

# Crear lista para iNEXT (eliminar ASVs con abundancia 0)
asv_por_sexo <- list(
  hembra = asv_hem[asv_hem > 0],
  macho   = asv_mac[asv_mac > 0]
)

# Verificar
str(asv_por_sexo)

# Números de Hill por sexo
out_sexo <- iNEXT(asv_por_sexo,
                  q        = c(0, 1, 2),
                  datatype = "abundance",
                  nboot    = 100,
                  knots    = 40)
out_sexo$DataInfo
out_sexo$AsyEst

write.csv(out_sexo$DataInfo, file = "Hill_sexo_datainfo.csv")
write.csv(out_sexo$AsyEst, file = "Hill_sexo_asyest.csv")


coverage_incidence_sex <- estimateD(asv_por_sexo,
                                q = c(0, 1, 2),
                                datatype = "abundance",
                                base = "coverage",
                                conf = 0.95, nboot = 100)
coverage_incidence_sex

# Gráfica de curvas de cobertura de muestra
# Ver qué tan completo está cada muestra/grupo
# ¿Mis muestras están bien secuenciadas? (cobertura)
ggiNEXT(out_sexo, type = 2) +
  scale_color_manual(name = "Floral sex",
                     values = c("hembra" = "#E6AB02", "macho" = "#4D9221"),
                     labels = c("hembra" = "Pistillate", "macho" = "Staminate"))+
  scale_fill_manual(name = "Floral sex",
                    values = c("hembra" = "#E6AB02", "macho" = "#4D9221"),
                    labels = c("hembra" = "Pistillate", "macho" = "Staminate"))+
  guides(
    shape    = "none"   # ← quita la leyenda negra de formas
    #linetype = "none"    # ← quita la leyenda de líneas si aparece duplicada
  ) + 
  theme_classic()+
  #theme_bw()
  labs(#title  = "Curvas de rarefacción por condición",
    x      = "Number of sequences")
#y      = "Diversidad (Hill)",
#color  = "Condición",
#fill   = "Condición")

# Curvas de diversidad por cobertura
# Comparar diversidad a igual completitud  
# ¿Difieren los grupos a igual completitud? (comparación justa)  
ggiNEXT(out_sexo, type = 3, facet.var = "Order.q")+
  scale_color_manual(name = "Floral sex",
                     values = c("hembra" = "#E6AB02", "macho" = "#4D9221"),
                     labels = c("hembra" = "Pistillate", "macho" = "Staminate"))+
  scale_fill_manual(name = "Floral sex",
                    values = c("hembra" = "#E6AB02", "macho" = "#4D9221"),
                    labels = c("hembra" = "Pistillate", "macho" = "Staminate"))+
  guides(
    shape    = "none"   # ← quita la leyenda negra de formas
    #linetype = "none"    # ← quita la leyenda de líneas si aparece duplicada
  )+
  theme_classic()
#facet_wrap(~Order.q, scales = "free_y",
#           labeller = labeller(Order.q = c("0" = "q = 0  (Riqueza)",
#                                          "1" = "q = 1  (Shannon)",
#                                         "2" = "q = 2  (Simpson)")))

### CURVAS DE RAREFACCIÓN POR TRATAMIENTO
# Separar muestras por tratamiento y sumar abundancias
muestras_t1 <- metadata %>% 
  filter(tratamiento == "t1") %>% 
  rownames()

muestras_t2 <- metadata %>% 
  filter(tratamiento == "t2") %>% 
  rownames()

# Sumar todas las muestras de cada tratamiento (un vector por tratamiento)
asv_t1 <- rowSums(asv_table[, muestras_t1])
asv_t2   <- rowSums(asv_table[, muestras_t2])

# Crear lista para iNEXT (eliminar ASVs con abundancia 0)
asv_por_tra <- list(
  t1 = asv_t1[asv_t1 > 0],
  t2   = asv_t2[asv_t2 > 0]
)

# Verificar
str(asv_por_tra)

# Números de Hill por tratamiento
out_trata <- iNEXT(asv_por_tra,
                   q        = c(0, 1, 2),
                   datatype = "abundance",
                   nboot    = 100,
                   knots    = 40)

out_trata$DataInfo
out_trata$AsyEst

write.csv(out_trata$DataInfo, file = "Hill_tratamiento_datainfo.csv")
write.csv(out_trata$AsyEst, file = "Hill_tratamiento_asyest.csv")


coverage_incidence_tra <- estimateD(asv_por_tra,
                                q = c(0, 1, 2),
                                datatype = "abundance",
                                base = "coverage",
                                conf = 0.95, nboot = 100)
coverage_incidence_tra

# Gráfica de curvas de cobertura de muestra
# Ver qué tan completo está cada muestra/grupo
# ¿Mis muestras están bien secuenciadas? (cobertura)
ggiNEXT(out_trata, type = 2) +
  scale_color_manual(name = "Treatment",
                     values = c("t1" = "#E6AB02", "t2" = "#4D9221"),
                     labels = c("t1" = "Unvisited", "t2" = "Visited"))+
  scale_fill_manual(name = "Floral sex",
                    values = c("t1" = "#E6AB02", "t2" = "#4D9221"),
                    labels = c("t1" = "Unvisited", "t2" = "Visited"))+
  guides(
    shape    = "none"   # ← quita la leyenda negra de formas
    #linetype = "none"    # ← quita la leyenda de líneas si aparece duplicada
  ) + 
  theme_classic()+
  #theme_bw()
  labs(#title  = "Curvas de rarefacción por condición",
    x      = "Number of sequences")
#y      = "Diversidad (Hill)",
#color  = "Condición",
#fill   = "Condición")

# Curvas de diversidad por cobertura
# Comparar diversidad a igual completitud  
# ¿Difieren los grupos a igual completitud? (comparación justa)  
ggiNEXT(out_trata, type = 3, facet.var = "Order.q") +
  scale_color_manual(name = "Treatment",
                     values = c("t1" = "#E6AB02", "t2" = "#4D9221"),
                     labels = c("t1" = "Unvisited", "t2" = "Visited"))+
  scale_fill_manual(name = "Floral sex",
                    values = c("t1" = "#E6AB02", "t2" = "#4D9221"),
                    labels = c("t1" = "Unvisited", "t2" = "Visited"))+
  guides(
    shape    = "none"   # ← quita la leyenda negra de formas
    #linetype = "none"    # ← quita la leyenda de líneas si aparece duplicada
  )+
  theme_classic()
#facet_wrap(~Order.q, scales = "free_y",
#           labeller = labeller(Order.q = c("0" = "q = 0  (Riqueza)",
#                                          "1" = "q = 1  (Shannon)",
#                                         "2" = "q = 2  (Simpson)")))

#### CURVAS DE RAREFACCIÓN POR GRUPOS
# Separar muestras por condición × sexo × tratamiento
muestras_dom_hem_t1 <- metadata %>% filter(condicion == "domesticada", sexo == "hembra", tratamiento == "t1") %>% rownames()
muestras_dom_hem_t2 <- metadata %>% filter(condicion == "domesticada", sexo == "hembra", tratamiento == "t2") %>% rownames()
muestras_dom_mac_t1 <- metadata %>% filter(condicion == "domesticada", sexo == "macho",  tratamiento == "t1") %>% rownames()
muestras_dom_mac_t2 <- metadata %>% filter(condicion == "domesticada", sexo == "macho",  tratamiento == "t2") %>% rownames()
muestras_sil_hem_t1 <- metadata %>% filter(condicion == "silvestre",   sexo == "hembra", tratamiento == "t1") %>% rownames()
muestras_sil_hem_t2 <- metadata %>% filter(condicion == "silvestre",   sexo == "hembra", tratamiento == "t2") %>% rownames()
muestras_sil_mac_t1 <- metadata %>% filter(condicion == "silvestre",   sexo == "macho",  tratamiento == "t1") %>% rownames()
muestras_sil_mac_t2 <- metadata %>% filter(condicion == "silvestre",   sexo == "macho",  tratamiento == "t2") %>% rownames()


muestras_t1_dom_hem <- metadata %>% filter(tratamiento == "t1", condicion == "domesticada", sexo == "hembra") %>% rownames()
muestras_t2_dom_hem <- metadata %>% filter(tratamiento == "t2", condicion == "domesticada", sexo == "hembra") %>% rownames()
muestras_t1_dom_mac <- metadata %>% filter(tratamiento == "t1", condicion == "domesticada", sexo == "macho") %>% rownames()
muestras_t2_dom_mac <- metadata %>% filter(tratamiento == "t2", condicion == "domesticada", sexo == "macho") %>% rownames()
muestras_t1_sil_hem <- metadata %>% filter(tratamiento == "t1", condicion == "silvestre",   sexo == "hembra") %>% rownames()
muestras_t2_sil_hem <- metadata %>% filter(tratamiento == "t2", condicion == "silvestre",   sexo == "hembra") %>% rownames()
muestras_t1_sil_mac <- metadata %>% filter(tratamiento == "t1", condicion == "silvestre",   sexo == "macho") %>% rownames()
muestras_t2_sil_mac <- metadata %>% filter(tratamiento == "t2", condicion == "silvestre",   sexo == "macho") %>% rownames()

# Verificar cuántas muestras hay en cada grupo
cat("t1-dom-hembra:", length(muestras_t1_dom_hem), "\n")
cat("t2-dom-hembra:", length(muestras_t2_dom_hem), "\n")
cat("t1-dom-macho: ", length(muestras_t1_dom_mac), "\n")
cat("t2-dom-macho: ", length(muestras_t2_dom_mac), "\n")
cat("t1-sil-hembra:", length(muestras_t1_sil_hem), "\n")
cat("t2-sil-hembra:", length(muestras_t2_sil_hem), "\n")
cat("t1-sil-macho: ", length(muestras_t1_sil_mac), "\n")
cat("t2-sil-macho: ", length(muestras_t2_sil_mac), "\n")

# Sumar abundancias por grupo
asv_t1_dom_hem <- rowSums(asv_table[, muestras_t1_dom_hem])
asv_t2_dom_hem <- rowSums(asv_table[, muestras_t2_dom_hem])
asv_t1_dom_mac <- rowSums(asv_table[, muestras_t1_dom_mac])
asv_t2_dom_mac <- rowSums(asv_table[, muestras_t2_dom_mac])
asv_t1_sil_hem <- rowSums(asv_table[, muestras_t1_sil_hem])
asv_t2_sil_hem <- rowSums(asv_table[, muestras_t2_sil_hem])
asv_t1_sil_mac <- rowSums(asv_table[, muestras_t1_sil_mac])
asv_t2_sil_mac <- rowSums(asv_table[, muestras_t2_sil_mac])

# Crear lista para iNEXT (eliminar ASVs con abundancia 0)
asv_por_grupo <- list(
  "UDP" = asv_t1_dom_hem[asv_t1_dom_hem > 0],
  "VDP" = asv_t2_dom_hem[asv_t2_dom_hem > 0],
  "UDS" = asv_t1_dom_mac[asv_t1_dom_mac > 0],
  "VDS" = asv_t2_dom_mac[asv_t2_dom_mac > 0],
  "UWP" = asv_t1_sil_hem[asv_t1_sil_hem > 0],
  "VWP" = asv_t2_sil_hem[asv_t2_sil_hem > 0],
  "UWS" = asv_t1_sil_mac[asv_t1_sil_mac > 0],
  "VWS" = asv_t2_sil_mac[asv_t2_sil_mac > 0]
)

# Correr iNEXT
out_grupos <- iNEXT(asv_por_grupo,
                    q        = c(0, 1, 2),
                    datatype = "abundance",
                    nboot    = 100,
                    knots    = 40)
out_grupos$DataInfo
out_grupos$AsyEst

#write.csv(out_grupos$DataInfo, file = "Hill_grupos_datainfo.csv")
#write.csv(out_grupos$AsyEst, file = "Hill_grupos_asyest.csv")


coverage_incidence_grupo <- estimateD(asv_por_grupo,
                                      q = c(0, 1, 2),
                                      datatype = "abundance",
                                      base = "coverage",
                                      conf = 0.95, nboot = 100)
coverage_incidence_grupo
write.csv(coverage_incidence_grupo, file = "Hill_numbers.csv")


colores_8 <- c(
  "UDP" = "#E74C3C",
  "VDP" = "#E74C3C",
  "UDS"  = "#E67E22",
  "VDS"  = "#E67E22",
  "UWP" = "#2980B9",
  "VWP" = "#2980B9",
  "UWS"  = "#27AE60",
  "VWS"  = "#27AE60"
)

# Gráfica tipo 1
p <- ggiNEXT(out_grupos, type = 1, facet.var = "Order.q")

# Modificar el grosor de las líneas en los layers
p$layers[[1]]$aes_params$linewidth <- 0.5  # líneas
p$layers[[2]]$aes_params$linewidth <- 0.5  # líneas extrapoladas

p +
  scale_color_manual(values = colores_8) +
  scale_fill_manual(values  = colores_8) +
  scale_shape_manual(values = c(
    "UDP" = 16,   # círculo relleno
    "VDP" = 21,   # círculo con borde
    "UDS"  = 17,   # triángulo relleno
    "VDS"  = 24,   # triángulo con borde
    "UWP" = 15,   # cuadrado relleno
    "VWP" = 22,   # cuadrado con borde
    "UWS"  = 18,   # diamante relleno
    "VWS"  = 23))+ # diamante con borde
  #guides(shape = "none", linetype = "none") +
  theme_classic() +
  labs(#title = "Números de Hill por condición × sexo × tratamiento",
    x     = "Number of sequences")
#y     = "Diversidad (Hill)",
#color = "Grupo")

# Gráfica tipo 2
p2 = ggiNEXT(out_grupos, type = 2) 

# Modificar el grosor de las líneas en los layers
p2$layers[[1]]$aes_params$linewidth <- 0.5  # líneas
p2$layers[[2]]$aes_params$linewidth <- 0.5  # líneas extrapoladas

p2 +
  scale_color_manual(values = colores_8) +
  scale_fill_manual(values  = colores_8) +
  scale_shape_manual(values = c(
    "UDP" = 16,   # círculo relleno
    "VDP" = 21,   # círculo con borde
    "UDS"  = 17,   # triángulo relleno
    "VDS"  = 24,   # triángulo con borde
    "UWP" = 15,   # cuadrado relleno
    "VWP" = 22,   # cuadrado con borde
    "UWS"  = 18,   # diamante relleno
    "VWS"  = 23)) + # diamante con borde
  theme_classic() 


# Gráfica tipo 3
p3 = ggiNEXT(out_grupos, type = 3, facet.var = "Order.q")

# Modificar el grosor de las líneas en los layers
p3$layers[[1]]$aes_params$linewidth <- 0.5  # líneas
p3$layers[[2]]$aes_params$linewidth <- 0.5  # líneas extrapoladas

p3 +
  scale_color_manual(values = colores_8) +
  scale_fill_manual(values  = colores_8) +
  scale_shape_manual(values = c(
    "UDP" = 16,   # círculo relleno
    "VDP" = 21,   # círculo con borde
    "UDS"  = 17,   # triángulo relleno
    "VDS"  = 24,   # triángulo con borde
    "UWP" = 15,   # cuadrado relleno
    "VWP" = 22,   # cuadrado con borde
    "UWS"  = 18,   # diamante relleno
    "VWS"  = 23))+ # diamante con borde
  theme_classic() 


### Extraer valores y graficar con intervalos de confianza

# Extraer valores en el punto observado (tamaño real de muestra)
hill_valores <- coverage_incidence_grupo %>% # solo el punto observado
  select(Assemblage, Order.q, qD, qD.LCL, qD.UCL) %>%
  mutate(Order.q = case_when(
    Order.q == 0 ~ "q = 0 (Riqueza)",
    Order.q == 1 ~ "q = 1 (Shannon)",
    Order.q == 2 ~ "q = 2 (Simpson)"
  ))


# Verificar
head(hill_valores)

# Paso 2: Separar las variables del nombre del grupo
hill_valores <- hill_valores %>%
  mutate(
    tratamiento = substr(Assemblage, 1, 1),  # U = unvisited, V = visited
    condicion   = substr(Assemblage, 2, 2),  # D = domesticated, W = wild
    sexo        = substr(Assemblage, 3, 3)   # P = pistillate, S = staminate
  )

head(hill_valores)

# Reordenar los niveles para que aparezcan en orden lógico en el eje X
#hill_valores$Assemblage <- factor(hill_valores$Assemblage, 
 #                                 levels = c("t1-Dom-Hembra", "t2-Dom-Hembra:",
  #                                           "t1-Dom-Macho", "t2-Dom-Macho",
    #                                         "t1-Sil-Hembra", "t2-Sil-Hembra",
      #                                         "t1-Sil-Macho", "t2-Sil-Macho"))

ggplot(hill_valores, aes(x = Assemblage, y = qD)) +
  # Intervalo de confianza
  geom_errorbar(aes(ymin = qD.LCL, ymax = qD.UCL),
                width     = 0.2,
                linewidth = 0.3) +
  # Puntos
  geom_point(size = 2.5, shape = 21,
             fill = "grey70", color = "black") +
  facet_wrap(~Order.q, scales = "free_y",
             labeller = labeller(Order.q = c(
               "q = 0 (Riqueza)"  = "Richness",
               "q = 1 (Shannon)"  = "Hill-Shannon",
               "q = 2 (Simpson)"  = "Hill-Simpson"))) +
  theme_bw(base_size = 13) +
  theme(axis.text.x    = element_text(angle = 45, hjust = 1, face = "italic"),
        strip.background = element_rect(fill = "grey90"),
        strip.text       = element_text(face = "bold"),
        legend.position  = "none") +
  labs(x = NULL,
       y = "Diversity (Hill)")


  scale_x_discrete(labels= c("Dom-Hembra-t1" = "DPU",
                             "Dom-Hembra-t2" = "DPV",
                             "Dom-Macho-t1" = "DSU",
                             "Dom-Macho-t2" = "DSV",
                             "Sil-Hembra-t1" = "WPU",
                             "Sil-Hembra-t2" = "WPV",
                             "Sil-Macho-t1" = "WSU",
                             "Sil-Macho-t2" = "WSV"))

###############################################################################
################ BETA DIVERSITY OF BACTERIA
BiocManager::install("DECIPHER")

# Cargando paquetes
library(DECIPHER)
library(phangorn)
library(Biostrings)


ntaxa(s16_6)

# 1. extraer directo de los rownames
seqs <- rownames(otu_table(s16_6))

# 2. Verificar que son secuencias y no IDs
head(seqs)  # deben verse como "ATCGGCTA..." 
nchar(seqs[1])  # debe tener longitud de ~250 pb para ITS2

# 3. Convertir a DNAStringSet y nombrar
seqs_dna <- DNAStringSet(seqs)
names(seqs_dna) <- seqs

# 4. Verificar
class(seqs_dna)
length(seqs_dna)  # debe coincidir con número de ASVs

# 5. Ahora sí correr el alineamiento
alineamiento <- AlignSeqs(seqs_dna)

# 6. Guardar el alineamiento en formato FASTA
writeXStringSet(alineamiento, 
                filepath = "alineamiento_its2.fasta")

# 7. Verificar que se guardó
file.exists("alineamiento_its2.fasta")

library(ape)

##### Importar árbol ####
arbol_ft <- read.tree("arbol_its2.nwk")

# Verificar
arbol_ft
plot(arbol_ft, cex = 0.2)  # grafica rápida para verificar

# Incorporar a phyloseq
ps_con_arbol <- merge_phyloseq(s16_6, phy_tree(arbol_ft))

# Verificar que quedó bien
ps_con_arbol
phy_tree(ps_con_arbol)

# UniFrac no ponderado
unifrac_nw <- UniFrac(ps_con_arbol, weighted = FALSE, normalized = TRUE)

# UniFrac ponderado
unifrac_w  <- UniFrac(ps_con_arbol, weighted = TRUE, normalized = TRUE)

# PERMANOVA y gráficas para UniFrac
library(vegan)

metadata <- as.data.frame(as.matrix(sample_data(ps_con_arbol)))

# UniFrac no ponderado
set.seed(1)
adonis2(unifrac_nw ~ condicion * sexo * tratamiento,
        data         = metadata,
        permutations = 999, by = "terms")

# UniFrac ponderado
set.seed(1)
adonis2(unifrac_w ~ condicion * sexo * tratamiento,
        data         = metadata,
        permutations = 999, by = "terms")

### Verificar homogeneidad de dispersión
# UniFrac no ponderado
disp_nw1 <- betadisper(unifrac_nw, metadata$sexo)
permutest(disp_nw1, permutations = 999)

disp_nw2 <- betadisper(unifrac_nw, metadata$tratamiento)
permutest(disp_nw2, permutations = 999)

# UniFrac ponderado
disp_w1 <- betadisper(unifrac_w, metadata$sexo)
permutest(disp_w1, permutations = 999)

disp_w2 <- betadisper(unifrac_w, metadata$tratamiento)
permutest(disp_w2, permutations = 999)


### Gráficas PCoA
#### Weighted UniFrac GLOBAL####
# PCoA plot using the weighted UniFrac as distance
ordination_16s_wei = ordinate(ps_con_arbol,
                              method="PCoA", "wunifrac")
# Plot
baseplot_16s_wei <- plot_ordination(ps_con_arbol, 
                                    ordination_16s_wei,
                                    color="tratamiento",
                                    #shape = "sexo",
                                    type = "samples") + #shape = "sexo" shape = "tratamiento"
  ggtitle("")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 4)+
  scale_color_discrete(name = "Treatment",
                       labels = c("t1" = "Unvisited",
                                  "t2" = "Visited"),
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




graficar_pcoa_unifrac <- function(dist_matrix, titulo) {
  
  pcoa <- cmdscale(dist_matrix, eig = TRUE, k = 2)
  
  df_pcoa <- data.frame(
    PC1         = pcoa$points[, 1],
    PC2         = pcoa$points[, 2],
    condicion   = metadata$condicion,
    sexo        = metadata$sexo,
    tratamiento = metadata$tratamiento
  )
  
  var_exp <- round(pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100, 1)
  
  ggplot(df_pcoa, aes(x = PC1, y = PC2,
                      color = condicion,
                      shape = tratamiento)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(aes(group = condicion,
                     color = condicion),
                 linewidth = 0.7) +
    scale_color_manual(values = c("domesticada" = "#C0392B",
                                  "silvestre"   = "#2980B9")) +
    scale_shape_manual(values = c("t1" = 16,
                                  "t2" = 17)) +
    theme_bw(base_size = 13) +
    theme(legend.position  = "right",
          strip.background = element_rect(fill = "grey90")) +
    labs(title = titulo,
         x     = paste0("PC1 (", var_exp[1], "%)"),
         y     = paste0("PC2 (", var_exp[2], "%)"),
         color = "Condición",
         shape = "Tratamiento")
}

# Generar las gráficas
p_nw <- graficar_pcoa_unifrac(unifrac_nw, "PCoA - UniFrac no ponderado")
p_w  <- graficar_pcoa_unifrac(unifrac_w,  "PCoA - UniFrac ponderado")

# Ver por separado
p_nw
p_w

# O juntas con patchwork
library(patchwork)
p_nw + p_w

######################### Bray-Curtis y Sørensen en R #########################
library(vegan)

# Extraer tabla de ASVs (muestras en filas, ASVs en columnas)
asv_matrix <- t(as.data.frame(otu_table(ps_con_arbol)))

# Bray-Curtis (basada en abundancias)
dist_bray <- vegdist(asv_matrix, method = "bray")

# Sørensen (basada en presencia/ausencia)
dist_sorensen <- vegdist(asv_matrix, method = "bray", binary = TRUE)
# binary = TRUE convierte a presencia/ausencia antes de calcular

# PERMANOVA para probar diferencias entre grupos
metadata <- as.data.frame(as.matrix(sample_data(ps_con_arbol)))

# Con Bray-Curtis
set.seed(1)
adonis2(dist_bray ~ condicion * sexo * tratamiento,
        data         = metadata,
        permutations = 999, by = "terms")

# Con Sørensen
set.seed(1)
adonis2(dist_sorensen ~ condicion * sexo * tratamiento,
        data         = metadata,
        permutations = 999, by = "terms")


# Bray-Curtis BETADISPER
disp_bray1 <- betadisper(dist_bray, metadata$condicion)
permutest(disp_bray1, permutations = 999)

disp_bray3 <- betadisper(dist_bray, metadata$tratamiento)
permutest(disp_bray3, permutations = 999)

# Sørensen
disp_sor1 <- betadisper(dist_sorensen, metadata$condicion)
permutest(disp_sor1, permutations = 999)

disp_sor3 <- betadisper(dist_sorensen, metadata$tratamiento)
permutest(disp_sor3, permutations = 999)


# Función para graficar PCoA
graficar_pcoa <- function(dist_matrix, titulo){
  
  pcoa <- cmdscale(dist_matrix, eig = TRUE, k = 2)
  
  df_pcoa <- data.frame(
    PC1       = pcoa$points[, 1],
    PC2       = pcoa$points[, 2],
    condicion  = metadata$condicion,
    sexo       = metadata$sexo,
    tratamiento = metadata$tratamiento
  )
  
  # Varianza explicada
  var_exp <- round(pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100, 1)
  
  ggplot(df_pcoa, aes(x = PC1, y = PC2,
                      color = condicion,
                      shape = tratamiento)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(aes(group = tratamiento,
                     color = tratamiento),
                 linewidth = 0.7) +
    theme_bw(base_size = 13) +
    labs(title = titulo,
         x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)"),
         color = "Condición",
         shape = "Tratamiento")
}

# Graficar ambas
graficar_pcoa(dist_bray,     "PCoA - Bray-Curtis")
graficar_pcoa(dist_sorensen, "PCoA - Sørensen")

#
p_bc  <- graficar_pcoa(dist_bray,      "Bray-Curtis")
p_sor <- graficar_pcoa(dist_sorensen,  "Sørensen")
p_nw  <- graficar_pcoa_unifrac(unifrac_nw, "UniFrac no ponderado")
p_w   <- graficar_pcoa_unifrac(unifrac_w,  "UniFrac ponderado")

(p_bc + p_sor) / (p_nw + p_w) +
  plot_annotation(title = "Diversidad beta - Microbiota de néctar",
                  theme = theme(plot.title = element_text(face = "bold",
                                                          size = 14)))


# Beta total = Recambio (turnover) + Anidamiento (nestedness)

#Recambio   → comunidades con especies completamente diferentes
#Anidamiento → una comunidad es subconjunto de otra

# En R con el paquete betapart
install.packages('rcdd')
install.packages("betapart")

library(betapart)


# Extraer tabla de ASVs (muestras en filas, ASVs en columnas)
comunidades <- t(as.data.frame(otu_table(ps)))

# Verificar dimensiones
dim(asv_matrix)  # debe ser: número de muestras × número de ASVs
head(rownames(asv_matrix))  # nombres de muestras
head(colnames(asv_matrix))  # nombres de ASVs

# Convertir a presencia/ausencia (0 y 1)
comunidades <- ifelse(asv_matrix > 0, 1, 0)

beta_total <- beta.pair(comunidades, index.family = "jaccard")
beta_total$beta.jtu   # recambio (turnover)
beta_total$beta.jne   # anidamiento (nestedness)
beta_total$beta.jac   # beta total


# Con abundancias usar bray en lugar de jaccard
beta_bray <- beta.pair.abund(asv_matrix, index.family = "bray")

beta_bray$beta.bray.bal   # recambio (balancedness)
beta_bray$beta.bray.gra   # anidamiento (gradient)
beta_bray$beta.bray       # beta total

boxplot(beta_bray$beta.bray.bal)
boxplot(beta_bray$beta.bray.gra)
boxplot(beta_bray$beta.bray)

### Gráficar
# Correr betapart (presencia/ausencia)
beta_total <- beta.pair(comunidades, index.family = "jaccard")

# Convertir cada matriz de distancia a data.frame
turnover   <- as.matrix(beta_total$beta.jtu)
nestedness <- as.matrix(beta_total$beta.jne)
beta_jac   <- as.matrix(beta_total$beta.jac)

# Convertir a formato largo
mat_to_long <- function(mat, nombre) {
  as.data.frame(mat) %>%
    mutate(sample1 = rownames(.)) %>%
    pivot_longer(-sample1, names_to = "sample2", values_to = nombre) %>%
    filter(sample1 != sample2)  # quitar diagonal
}

df_turn <- mat_to_long(turnover,   "turnover")
df_nest <- mat_to_long(nestedness, "nestedness")
df_jac  <- mat_to_long(beta_jac,   "beta_total")

# Unir los tres componentes
df_beta <- df_turn %>%
  left_join(df_nest, by = c("sample1", "sample2")) %>%
  left_join(df_jac,  by = c("sample1", "sample2"))

head(df_beta)

# Agregar metadatos
metadata <- as.data.frame(as.matrix(sample_data(ps_con_arbol)))

# Agregar metadatos de sample1
df_beta <- df_beta %>%
  left_join(metadata %>% 
              mutate(sample1 = rownames(.)) %>%
              select(sample1, 
                     condicion1   = condicion,
                     sexo1        = sexo,
                     tratamiento1 = tratamiento),
            by = "sample1") %>%
  # Agregar metadatos de sample2
  left_join(metadata %>%
              mutate(sample2 = rownames(.)) %>%
              select(sample2,
                     condicion2   = condicion,
                     sexo2        = sexo,
                     tratamiento2 = tratamiento),
            by = "sample2") %>%
  # Crear columna de comparación entre grupos
  mutate(
    comp_condicion   = paste(pmin(condicion1,   condicion2),
                             pmax(condicion1,   condicion2),   sep = " vs "),
    comp_sexo        = paste(pmin(sexo1,        sexo2),
                             pmax(sexo1,        sexo2),        sep = " vs "),
    comp_tratamiento = paste(pmin(tratamiento1, tratamiento2),
                             pmax(tratamiento1, tratamiento2), sep = " vs ")
  )

# Convertir a formato largo para graficar turnover y nestedness juntos
df_long <- df_beta %>%
  pivot_longer(cols      = c(turnover, nestedness),
               names_to  = "componente",
               values_to = "valor") %>%
  mutate(componente = factor(componente,
                             levels = c("turnover", "nestedness"),
                             labels = c("Recambio (turnover)",
                                        "Anidamiento (nestedness)")))

View(df_long)
###
# Crear grupo combinado de comparación
df_long <- df_long %>%
  mutate(grupo = paste0(comp_condicion, "\n",
                        comp_sexo, "\n",
                        comp_tratamiento))

ggplot(df_long, aes(x = comp_condicion, y = valor,
                    fill = comp_tratamiento)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5,
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("t1 vs t1" = "#E74C3C",
                               "t1 vs t2" = "#E67E22",
                               "t2 vs t2" = "#27AE60")) +
  facet_grid(componente ~ comp_sexo, scales = "free_y") +
  theme_bw(base_size = 12) +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey90"),
        strip.text       = element_text(face = "bold"),
        legend.position  = "bottom") +
  labs(title = "Recambio y anidamiento por condición × sexo × tratamiento",
       x     = "Condición",
       y     = "Distancia de Jaccard",
       fill  = "Comparación de tratamiento")


###############################################################################
###############################################################################
# install microViz
#remotes::install_github("david-barnett/microViz")
library("microViz")
library(phyloseq)
library(ggplot2)
library(patchwork) # for arranging groups of plots


# Paquete microbiome
# 16S
# Paquete microbiome
transform <- microbiome::transform
s16_prop <- transform(s16_6, "compositional")
otu_table(s16_prop)[1:10, 1:10]
head(sample_sums(s16_prop))

s16_prop %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  comp_barplot(tax_level = "Family", n_taxa = 15, facet_by = "tratamiento")+
  coord_flip() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# modificada
s16_prop %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  comp_barplot(tax_level = "Family", n_taxa = 15)+
  facet_wrap(~ condicion + sexo + tratamiento, scales = "free")+ #, scales = "free"
  coord_flip() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())


# ITS (fungi)
transform <- microbiome::transform
ps_prop <- transform(s16_6, "compositional")
otu_table(ps_prop)[1:10, 1:10]
head(sample_sums(ps_prop))

ps_prop %>%
  tax_fix() %>%
  comp_barplot(tax_level = "Genus", n_taxa = 15, facet_by = "tratamiento",
               group_by = "")+
  coord_flip() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

#facet_grid(labeller = labeller(tratamiento = trata))

###############################################################################
###############################################################################
### ANÁLISIS CON DESEQ2
# Paso 1. Instalar y cargar paquetes
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)
library(phyloseq)

# 2. Convertir phyloseq a DESeq2
# primero agrupar a nivel de familia
s16_familia <- tax_glom(s16_6, taxrank = "Family")
# Para la comparación por condición (domesticated vs wild)
ds_condicion <- phyloseq_to_deseq2(s16_familia, ~ condicion)
# 3. Estimar factores de tamaño con conteos cero
ds_condicion <- estimateSizeFactors(ds_condicion, 
                                    type = "poscounts")
# 4. Correr DESeq2
ds_condicion <- DESeq(ds_condicion, 
                      test = "Wald", 
                      fitType = "parametric")
# 5. Extraer resultados
res_condicion <- results(ds_condicion,
                         contrast = c("condicion", "domesticada", "silvestre"),
                         alpha = 0.05)
# Ver resumen
summary(res_condicion)

# Para la comparación por SEXO (hembra vs macho)
ds_sexo <- phyloseq_to_deseq2(s16_familia, ~ sexo)
# 3. Estimar factores de tamaño con conteos cero
ds_sexo <- estimateSizeFactors(ds_sexo, 
                                    type = "poscounts")
# 4. Correr DESeq2
ds_sexo <- DESeq(ds_sexo, 
                      test = "Wald", 
                      fitType = "parametric")
# 5. Extraer resultados
res_sexo <- results(ds_sexo,
                         contrast = c("sexo", "hembra", "macho"),
                         alpha = 0.05)
# Ver resumen
summary(res_sexo)
# Ver qué familia es significativa en sexo
res_sexo_df <- as.data.frame(res_sexo) %>%
  rownames_to_column("ASV") %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(desc(log2FoldChange))

tax_tabla_fam <- as.data.frame(tax_table(s16_familia)) %>%
  rownames_to_column("ASV")

res_sexo_df <- left_join(res_sexo_df, tax_tabla_fam, by = "ASV")

res_sexo_df %>% dplyr::select(Family, log2FoldChange, padj)

# Comparación por tratamiento (efecto más fuerte en PERMANOVA)
s16_familia <- tax_glom(s16_6, taxrank = "Family")

ds_tratamiento <- phyloseq_to_deseq2(s16_familia, ~ tratamiento)
ds_tratamiento <- estimateSizeFactors(ds_tratamiento, type = "poscounts")
ds_tratamiento <- DESeq(ds_tratamiento, test = "Wald", fitType = "parametric")

res_tratamiento <- results(ds_tratamiento,
                           contrast = c("tratamiento", "t2", "t1"),
                           alpha = 0.05)
summary(res_tratamiento)

# Extraer resultados como dataframe
res_trat_df <- as.data.frame(res_tratamiento) %>%
  rownames_to_column("ASV") %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(desc(log2FoldChange))

# Agregar taxonomía
tax_tabla <- as.data.frame(tax_table(s16_familia)) %>%
  rownames_to_column("ASV")

res_trat_df <- left_join(res_trat_df, tax_tabla, by = "ASV")

# Ver qué familias son
res_trat_df %>% dplyr::select(Family, log2FoldChange, padj)

# Graficando
ggplot(res_trat_df, aes(x = reorder(Family, log2FoldChange), 
                        y = log2FoldChange)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(x = "Family", 
       y = "log2 Fold Change (visited vs. unvisited)",
       title = "Differentially abundant bacterial families by treatment") +
  theme_bw()

### VOLCANO PLOT
library(ggplot2)
library(ggrepel)  # para etiquetas que no se sobrepongan

# Preparar dataframe completo (incluyendo no significativos)
res_trat_volcano <- as.data.frame(res_tratamiento) %>%
  rownames_to_column("ASV") %>%
  filter(!is.na(padj)) %>%
  left_join(tax_tabla_fam, by = "ASV") %>%
  mutate(
    significativo = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "More abundant in visited",
      padj < 0.05 & log2FoldChange < 0 ~ "More abundant in unvisited",
      TRUE ~ "Not significant"
    )
  )

# Volcano plot
ggplot(res_trat_volcano, aes(x = log2FoldChange, 
                             y = -log10(padj),
                             color = significativo,
                             label = Family)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), 
             linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 0, 
             linetype = "dashed", color = "gray40") +
  geom_text_repel(data = filter(res_trat_volcano, padj < 0.05),
                  size = 3.5,
                  fontface = "italic",
                  max.overlaps = 20) +
  scale_color_manual(values = c(
    "More abundant in visited"   = "steelblue",
    "More abundant in unvisited" = "tomato",
    "Not significant"            = "gray70"
  )) +
  labs(x = "log2 Fold Change (visited vs. unvisited)",
       y = "-log10 (adjusted p-value)",
       color = NULL) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid.minor = element_blank())

##############################################################################
###  differential abundance (DA) analysis of microbial taxa with ANCOMBC2 
#devtools::install_github("FrederickHuangLin/ANCOMBC")
citation("ANCOMBC")
library(ANCOMBC)
library(ggplot2)

# CONDICION
ancom_con_16s <- ancombc2(
  data = s16_6,
  tax_level = "Family",
  fix_formula = "condicion",
  p_adj_method = "BH",
  alpha = 0.05,
  prv_cut = 0.10
)

# Ver resultados significativos
res_con_16s <- ancom_con_16s$res
res_con_16s_sig <- res_con_16s[res_con_16s$diff_condicionsilvestre == TRUE, 
                                 c("taxon", "lfc_condicionsilvestre", 
                                   "q_condicionsilvestre")]
print(res_con_16s_sig)

# Gráfica de Barras por condicion
# Preparar datos
res_con_16s_plot <- res_con_16s_sig %>%
  mutate(
    direccion = ifelse(lfc_condicionsilvestre > 0, 
                       "More abundant in domesticated", 
                       "More abundant in wild"),
    taxon = gsub("f__", "", taxon)  # limpiar prefijo f__
  )

# Grafica
ggplot(res_con_16s_plot, 
       aes(x = reorder(taxon, lfc_condicionsilvestre),
           y = lfc_condicionsilvestre,
           fill = direccion)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  labs(x = "Family",
       y = "log2 Fold Change (Domesticated vs. Wild)",
       fill = NULL) +
  theme_bw() +
  theme(legend.position = "top")

# SEXO FLORAL
ancom_sex_16s <- ancombc2(
  data = s16_6,
  tax_level = "Family",
  fix_formula = "sexo",
  p_adj_method = "BH",
  alpha = 0.05,
  prv_cut = 0.10
)

# Ver resultados significativos
res_sex_16s <- ancom_sex_16s$res
res_sex_16s_sig <- res_sex_16s[res_sex_16s$diff_sexomacho == TRUE, 
                               c("taxon", "lfc_sexomacho", 
                                 "q_sexomacho")]
print(res_sex_16s_sig)

# Gráfica de Barras por sexo
# Preparar datos
res_sex_16s_plot <- res_sex_16s_sig %>%
  mutate(
    direccion = ifelse(lfc_sexomacho > 0, 
                       "More abundant in staminate", 
                       "More abundant in pistillate"),
    taxon = gsub("f__", "", taxon)  # limpiar prefijo f__
  )

# Grafica
ggplot(res_sex_16s_plot, 
       aes(x = reorder(taxon, lfc_sexomacho),
           y = lfc_sexomacho,
           fill = direccion)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  labs(x = "Family",
       y = "log2 Fold Change (Staminate vs. Pistillate)",
       fill = NULL) +
  theme_bw() +
  theme(legend.position = "top")


# Tratamiento
ancom_trat_16s <- ancombc2(
  data = s16_6,
  tax_level = "Family",
  fix_formula = "tratamiento",
  p_adj_method = "BH",
  alpha = 0.05,
  prv_cut = 0.10
)

# Ver resultados significativos
res_trat_16s <- ancom_trat_16s$res
res_trat_16s_sig <- res_trat_16s[res_trat_16s$diff_tratamientot2 == TRUE, 
                                 c("taxon", "lfc_tratamientot2", 
                                   "q_tratamientot2")]
print(res_trat_16s_sig)

# Gráfica de Barras de tratamiento
# Preparar datos
res_trat_16s_plot <- res_trat_16s_sig %>%
  mutate(
    direccion = ifelse(lfc_tratamientot2 > 0, 
                       "More abundant in visited", 
                       "More abundant in unvisited"),
    taxon = gsub("f__", "", taxon)  # limpiar prefijo f__
  )

# Grafica
ggplot(res_trat_16s_plot, 
       aes(x = reorder(taxon, lfc_tratamientot2),
           y = lfc_tratamientot2,
           fill = direccion)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  labs(x = "Family",
       y = "log2 Fold Change (visited vs. unvisited)",
       fill = NULL) +
  theme_bw() +
  theme(legend.position = "top")

###############################################################################
### Código para Top 10 Familias Global (Sin Agrupar por Variables)
# A. Aglomerar a nivel de Familia
ps_fam1 <- tax_glom(s16_6, taxrank = "Family", NArm = FALSE)

# B. Transformar a abundancia relativa
ps_rel1 <- transform_sample_counts(ps_fam1, function(x) x / sum(x))

# C. Psmelt
df_melt1 <- psmelt(ps_rel1)

# D. Agrupar SOLO por Familia (sin condicion ni sexo)
df_grouped_global <- df_melt1 %>%
  dplyr::group_by(Family) %>%                          # ← solo Family
  dplyr::summarize(MeanRelAbund = mean(Abundance),
                   .groups = "drop") %>%
  arrange(desc(MeanRelAbund))                          # ordenar de mayor a menor

# Ver el resultado
print(df_grouped_global, n = 20)

# E. Top 10 familias más abundantes globalmente
top10_global <- df_grouped_global %>%
  slice_max(order_by = MeanRelAbund, n = 10,
            with_ties = FALSE)

cat("Top 10 familias más abundantes:\n")
print(top10_global)
write.csv(top10_global, file = "top10_familias_bacterias.csv")

# F. Crear categoría "Other" para el resto
df_final_global <- df_grouped_global %>%
  mutate(Taxonomy_Display = ifelse(
    Family %in% top10_global$Family,
    as.character(Family),
    "Other"
  )) %>%
  group_by(Taxonomy_Display) %>%
  dplyr::summarize(TotalMeanAbund = sum(MeanRelAbund),
            .groups = "drop") %>%
  mutate(Percentage = TotalMeanAbund * 100) %>%
  arrange(desc(TotalMeanAbund))

# Verificar que suma 100%
cat("Suma total:", sum(df_final_global$Percentage), "%\n")

View(df_final_global)
write.csv(df_final_global, "abundancias_top10_global.csv")

###############################################################################
#### EXTRAYENDO LAS MUESTRAS DE FLORES VISITADAS PARA CONSTRUIR REDES
### DE INTERACCIÓN
visitadas <- subset_samples(ps_con_arbol, tratamiento == "t2")
visitadas

#### Obteniendo las abundancias relativas por especie
# 1. Preparación del objeto Phyloseq
# A. Aglomerar a nivel de Familia. Esto suma los reads de todos los ASVs/OTUs
# que pertenecen a la misma familia. Reduce la complejidad drásticamente.
ps_fam <- tax_glom(visitadas, taxrank = "Family", NArm = FALSE)

# B. Transformar a abundancia relativa (escala 0-1)
# Esto estandariza las muestras independientemente de su profundidad de secuenciación.
ps_rel <- transform_sample_counts(ps_fam, function(x) x / sum(x))

# C. Convertir el objeto phyloseq a un dataframe largo ('melted' dataframe)
df_melt <- psmelt(ps_rel)

# D. Agrupar por Especie y Familia, y calcular la abundancia relativa PROMEDIO
# para esa especie de planta.
df_grouped <- df_melt %>%
  dplyr::group_by(condicion, sexo, Family) %>% 
  dplyr::summarize(MeanRelAbund = mean(Abundance), .groups = "drop")

# E. Identificar el Top 10 de familias más abundantes POR ESPECIE.
# Usamos slice_max para quedarnos con las 10 filas con mayor abundancia promedio por grupo.
top10_per_species <- df_grouped %>%
  group_by(condicion, sexo) %>%
  slice_max(order_by = MeanRelAbund, n = 10, with_ties = FALSE) %>%
  mutate(Classification = Family) %>% # Guardamos el nombre para el cruce
  ungroup()

# F. Crear la categoría "Otras".
# Cruzamos la tabla agrupada completa con la tabla del Top 10.
df_final <- df_grouped %>%
  left_join(select(top10_per_species, condicion, sexo, Family, Classification), 
            by = c("condicion","sexo", "Family")) %>%
  # Si la familia NO estaba en el top10 (Classification es NA), la llamamos "Otras"
  dplyr::mutate(Taxonomy_Display = ifelse(is.na(Classification), 
                                          "Other", 
                                          as.character(Family))) %>%
  # Volvemos a agrupar y sumar para colapsar todas las "Otras" en una sola fila por especie
  dplyr::group_by(condicion, sexo, Taxonomy_Display) %>%
  dplyr::summarize(TotalMeanAbund = sum(MeanRelAbund), .groups = 'drop') %>%
  # Convertir abundancia relativa (0-1) a porcentaje (0-100) para el gráfico
  dplyr::mutate(Percentage = TotalMeanAbund * 100)
View(df_final)
write.csv(df_final, file = "abundancias_visitadas_abril.csv")












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



###############################################################################
###############################################################################
#### Redes de interacción DE POLINIZADORES EN PARCELA DE MORELIA ####
## Estableciendo directorio en el cual se encuentran las bases de datos ##
setwd("/home/luis/Documentos/3_Doctorado_UNAM/a-proyecto/analisis_florales_tesis_2022")
more2021 <- read.csv("visitas_morelia_2021.csv", header = T)
head(more2021)
str(more2021)
more2021$ID.planta
more2021$codigo.video

# Cargando librerias
library(fossil)
library (bipartite)


### Horas totales de grabación entre 2021 y 2022: 128:52:28
# Horas de grabación por especie
# CPF: 7:31:30
# CPP: 30:41:42
# CAS: 31:25:02
# CAA: 22:50:56
# CM: 13:18:00

30+23+13
7.5+31.5
66+39


#To calculate the visitation rate for each floral visitor group we divided the
# number of visits by the total duration of the video recording 
# (decimal time = time x 1440) for each flower.


# numero de visitas por especie
tapply(more2021$especie, more2021$condicion, length)
tapply(more2021$especie, more2021$Especie.planta, length)

#
visitas <- more2021 %>%
  dplyr::group_by(condicion, sexo, especie) %>%
  dplyr::summarize(num_visitas = length(especie))

## frecuencia de visitas
frec <- more2021 %>%
  dplyr::group_by(año, fecha, condicion, ID.planta, sexo, especie,
                  duración.decimal) %>%
  dplyr::summarize(num_visitas = length(especie)) %>% 
  dplyr::mutate(frecu = num_visitas/duración.decimal)

frec2 <- frec %>%
  dplyr::group_by(condicion, sexo, especie) %>%
  dplyr::summarize(mean_frec = mean(frecu))


View(frec2)
#
ggplot(frec, aes(condicion, fill = especie))+
  geom_bar()
#
ggplot(more2021, aes(condicion, fill = especie))+
  geom_bar(position = "fill")
#
ggplot(more2021, aes(condicion, fill = especie))+
  geom_bar()+
  facet_grid(.~ sexo)
#
ggplot(more2021, aes(condicion, fill = especie))+
  geom_bar(position = "fill")+
  facet_grid(.~ sexo)

### Realizando redes de interacción ###
visitas2 = visitas %>%
  unite("plant_sex", condicion:sexo)
df_vis <- create.matrix(as.data.frame(visitas2), tax.name = "especie",
                        locality = "plant_sex",
                        abund.col = "num_visitas",
                        abund = TRUE)
df1_vis <- t(df_vis)

# frecuencias de visita
frec3 = frec2 %>%
  unite("plant_sex", condicion:sexo)
#write.csv(frec3, file = "frecuencias_visitas.csv")

df <- create.matrix(as.data.frame(frec3), tax.name = "especie",
                    locality = "plant_sex",
                    abund.col = "mean_frec",
                    abund = TRUE)
df1 <- t(df)

#cargar la libreria bipartite


#Dibujar la red de interacciones
# numero de visita
plotweb(df1_vis)
plotweb(df1_vis,method="normal", y.lim=c(-1,2.5), labsize=1.2, text.rot=90,
        col.high = "brown4", col.low = "#D9A404")

# frecuencia de visita
plotweb(df1,method="normal", y.lim=c(-1,2.5), labsize=1.2, text.rot=90,
        col.high = "brown4", col.low = "#D9A404")
plotweb(df1, sorting = "normal", higher_color = "#E6AB02",
        lower_color = "#4D9221")


specieslevel(df1)
beeR<-specieslevel(df1)
summary(beeR)
str(beeR)
class(beeR)

colors()

#Calcular la m?tricas de red y guardarlos en un archivo de texto
bee_strength<-beeR$'higher level'
write.table(bee_strength, "bee_strength")

plant_strength<-beeR$'lower level'
write.table(plant_strength, "plant_strength")

visweb(df1) 
visweb(df1, type = "diagonal")
networklevel(df1) 


#### Cargando las 2 bases de datos, de frecuencia de visita y de 
### bacterias 
frec2
top10_per_species
View(frec2)
View(top10_per_species)

# Corregir la columna sexo en visitantes
frec2 <- frec2 %>%
  mutate(sexo = case_when(
    sexo == "H" ~ "hembra",
    sexo == "M" ~ "macho",
    TRUE        ~ sexo
  ))

# Verificar que quedó bien
unique(frec2$sexo)
# Debe decir: "hembra" "macho"

# Volver a unir las tablas
tabla_integrada <- left_join(
  top10_per_species,
  frec2,
  by = c("condicion", "sexo")
)

# Verificar que ya no hay NAs
head(tabla_integrada)
View(tabla_integrada)
sum(is.na(tabla_integrada$mean_frec))

## Tabla que integra visitantes y bacterias
integrada = tabla_integrada %>%
  unite("plant_sex", condicion:sexo)
View(integrada)
# Matríz de bacterias
df2 <- create.matrix(as.data.frame(integrada), tax.name = "Family",
                    locality = "plant_sex",
                    abund.col = "MeanRelAbund",
                    abund = TRUE)
df3 <- t(df2)
# 
specieslevel(df3)
beeR<-specieslevel(df3)
summary(beeR)
str(beeR)
class(beeR)

#Calcular la m?tricas de red y guardarlos en un archivo de texto
bee_strength <- beeR$'higher level'
write.table(bee_strength, "bee_strength")

plant_strength<-beeR$'lower level'
write.table(plant_strength, "plant_strength")

visweb(df3) 
visweb(df3, type = "diagonal")
networklevel(df3) 

# Matriz de visitantes
df4 <- create.matrix(as.data.frame(integrada), tax.name = "especie",
                     locality = "plant_sex",
                     abund.col = "mean_frec",
                     abund = TRUE)
df5 <- t(df4)

# Visualizar redes lado a lado
par(mfrow = c(1, 2),
    mar   = c(8, 2, 2, 2))
# Red 1: plantas → bacterias
plotweb(df3, sorting = "normal", higher_color = "#E6AB02",
        lower_color = "#4D9221", horizontal = TRUE)

# Red 2: visitantes → plantas
plotweb(df5, sorting = "normal", higher_color = "#E6AB02",
        lower_color = "#4D9221",horizontal = TRUE)

# Verificando orientación de las redes
# web1: filas = nivel inferior, columnas = nivel medio
# web2: filas = nivel medio,    columnas = nivel superior

# Verificar orientación
dim(df5)  # visitantes × plantas
dim(df3)  # plantas × bacterias

# El nivel medio (plantas) debe coincidir:
# colnames(web1) == rownames(web2)
colnames(df5)  # deben ser las plantas
rownames(df3)  # deben ser las plantas

# Verificar que coinciden exactamente
all(sort(colnames(df5)) == sort(rownames(df3)))

# Transponer si es necesario
# web1 necesita: filas=visitantes, columnas=plantas
red_visitantes_ok <- t(df5)  # o t(red_visitantes)

# web2 necesita: filas=plantas, columnas=bacterias
red_microbiota_ok <- df3  # o t(red_microbiota)

# colnames(web1) == rownames(web2)
colnames(red_visitantes_ok)  # deben ser las plantas
rownames(red_microbiota_ok)  # deben ser las plantas
# Verificar que coinciden exactamente
all(sort(colnames(red_visitantes_ok)) == sort(rownames(red_microbiota_ok)))


plotweb2(red_microbiota_ok, red_visitantes_ok, method = "normal")



# Normalizar cada matriz entre 0 y 1
normalizar <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

red_vis_norm  <- normalizar(red_visitantes_ok)
red_micro_norm <- normalizar(red_microbiota_ok)

plotweb2(red_vis_norm, red_micro_norm, method = "normal")

######################## Datos de volumen de néctar ############################
volumen  = read.csv("flores_morelia.csv", header = T)
head(volumen)
str(volumen)
volumen2 = subset(volumen, NV != 0)
str(volumen2)



#### Análisis de nectar entre tratamientos ####
setwd("/home/luis/Documentos/3_Doctorado_UNAM/a-proyecto/analisis_florales_tesis_2022/datos_completos/")

azu = read.csv("azucares_2022.csv", header = T)
head(azu)
str(azu)
azu2 =subset(azu, sitio == "morelia")
head(azu2)

tapply(azu2$Suma.total, list(azu2$condicion, azu2$sexo, azu2$trat), length)

tapply(azu2$Suma.total, list(azu2$spp, azu2$sexo, azu2$trat), length)
tapply(azu2$Suma.total, list(azu2$spp, azu2$trat), length)
tapply(azu2$Suma.total, azu2$trat, length)

# Pruebas de normalidad de cada azucar
shapiro.test(azu2$fructosa)
shapiro.test(azu2$glucosa)
shapiro.test(azu2$sacarosa)
shapiro.test(azu2$Suma.total)

### Histogramas de cada azucar
hist(azu2$fructosa)
hist(azu2$glucosa)
hist(azu2$sacarosa)
hist(azu2$Suma.total)

azu2$fructosa
azu2$glucosa

# Modelo 1: análisis de varianza de Concentración de Sacarosa
mod1 = glm(fructosa ~ condicion * sexo * trat, 
           data = azu2, family = Gamma(link = "log"))
summary(mod1)
library(car)
Anova(mod1)

# Modelo 2: análisis de varianza de Concentración de Sacarosa
mod2 = glm(glucosa ~ condicion * sexo * trat, 
           data = azu2, family = Gamma(link = "log"))
summary(mod2)
Anova(mod1)

# Modelo 3: análisis de varianza de Concentración de Sacarosa
mod3 = aov(azu2$sacarosa ~ azu2$condicion * azu2$sexo * azu2$trat)
summary(mod3)

# Modelo 4: análisis de varianza de Concentración total
mod4 = aov(azu2$Suma.total ~ azu2$condicion * azu2$sexo * azu2$trat)
summary(mod4)

TukeyHSD(mod1)

# media
media2 = azu2 %>%
  dplyr::group_by(condicion, sexo, trat) %>%
  dplyr::summarize(fruc_m = mean(fructosa), glu_m = mean(glucosa),
                   sac_m = mean(sacarosa),
                   tot_m = mean(Suma.total))
media2
# sd
sd2 = azu2 %>%
  dplyr::group_by(condicion, sexo, trat) %>%
  dplyr::summarize(fruc_sd = sd(fructosa), glu_sd = sd(glucosa),
                   sac_sd = sd(sacarosa),
                   tot_sd = sd(Suma.total))
sd2


### 
boxplot(azu2$Suma.total ~ azu2$condicion)
boxplot(azu2$Suma.total ~ azu2$sexo)
boxplot(azu2$Suma.total ~ azu2$trat)

boxplot(azu2$Suma.total ~ azu2$condicion + azu2$trat)



library(dplyr)
resultado <- left_join(media2, sd2, by = c("condicion", "sexo", "trat"))
resultado

library(ggplot2)

## Gráfica de concentración de fructosa
fru1 = ggplot(resultado, aes(x = interaction(trat, sexo, condicion), 
                             y = fruc_m,
                             color = interaction(trat, sexo, condicion))) + #fill = condicion
  geom_col(position = "dodge", width = 0.7,
           fill = NA,
           linewidth = 0.8) +
  geom_errorbar(aes(ymin = fruc_m - fruc_sd, 
                    ymax = fruc_m + fruc_sd),
                width = 0.25) +
  labs(
    #title = "Azúcares totales por grupo experimental",
    x = "",
    y = "Fructose",
    #fill = "Condición"
  ) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
scale_color_manual(
  name   = "Group",   # título de la leyenda
  values = c(         # colores para cada combinación
    "t1.h.domesticada" = "#E6AB02",
    "t2.h.domesticada"  = "#4D9221",
    "t1.m.domesticada"   = "#1F78B4",
    "t2.m.domesticada"    = "#E31A1C",
    "t1.h.silvestre" = "#FF7F00",
    "t2.h.silvestre"  = "#6A3D9A",
    "t1.m.silvestre"   = "#B15928",
    "t2.m.silvestre"    = "#33A02C"
  ),
  labels = c(         # nombres que aparecen en la leyenda
    "t1.h.domesticada" = "UDP",
    "t2.h.domesticada"  = "VDP",
    "t1.m.domesticada"   = "UDS",
    "t2.m.domesticada"    = "VDS",
    "t1.h.silvestre" = "UWP",
    "t2.h.silvestre"  = "VWP",
    "t1.m.silvestre"   = "UWS",
    "t2.m.silvestre"    = "VWS"
  )
) + 
  scale_x_discrete(labels=c("t1.h.domesticada" = "UDP",
                            "t2.h.domesticada"  = "VDP",
                            "t1.m.domesticada"   = "UDS",
                            "t2.m.domesticada"    = "VDS",
                            "t1.h.silvestre" = "UWP",
                            "t2.h.silvestre"  = "VWP",
                            "t1.m.silvestre"   = "UWS",
                            "t2.m.silvestre"    = "VWS"))
fru1

## Gráfica de concentración de glucosa
glu1 = ggplot(resultado, aes(x = interaction(trat, sexo, condicion), 
                             y = glu_m,
                             color = interaction(trat, sexo, condicion))) + #fill = condicion
  geom_col(position = "dodge", width = 0.7,
           fill = NA, linewidth = .8) +
  geom_errorbar(aes(ymin = glu_m - glu_sd, 
                    ymax = glu_m + glu_sd),
                width = 0.25) +
  labs(
    #title = "Azúcares totales por grupo experimental",
    x = "",
    y = "Glucose",
    #fill = "Condición"
  ) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(
    name   = "Group",   # título de la leyenda
    values = c(         # colores para cada combinación
      "t1.h.domesticada" = "#E6AB02",
      "t2.h.domesticada"  = "#4D9221",
      "t1.m.domesticada"   = "#1F78B4",
      "t2.m.domesticada"    = "#E31A1C",
      "t1.h.silvestre" = "#FF7F00",
      "t2.h.silvestre"  = "#6A3D9A",
      "t1.m.silvestre"   = "#B15928",
      "t2.m.silvestre"    = "#33A02C"
    ),
    labels = c(         # nombres que aparecen en la leyenda
      "t1.h.domesticada" = "UDP",
      "t2.h.domesticada"  = "VDP",
      "t1.m.domesticada"   = "UDS",
      "t2.m.domesticada"    = "VDS",
      "t1.h.silvestre" = "UWP",
      "t2.h.silvestre"  = "VWP",
      "t1.m.silvestre"   = "UWS",
      "t2.m.silvestre"    = "VWS"
    )
  ) + 
  scale_x_discrete(labels=c("t1.h.domesticada" = "UDP",
                            "t2.h.domesticada"  = "VDP",
                            "t1.m.domesticada"   = "UDS",
                            "t2.m.domesticada"    = "VDS",
                            "t1.h.silvestre" = "UWP",
                            "t2.h.silvestre"  = "VWP",
                            "t1.m.silvestre"   = "UWS",
                            "t2.m.silvestre"    = "VWS"))
glu1
## Gráfica de concentración de sacarosa
sac1 = ggplot(resultado, aes(x = interaction(trat, sexo, condicion), 
                             y = sac_m,
                             color = interaction(trat, sexo, condicion))) + #fill = condicion
  geom_col(position = "dodge", width = 0.7,
           fill = NA, linewidth = .8) +
  geom_errorbar(aes(ymin = sac_m - sac_sd, 
                    ymax = sac_m + sac_sd),
                width = 0.25) +
  labs(
    #title = "Azúcares totales por grupo experimental",
    x = "",
    y = "Sacarose",
    #fill = "Condición"
  ) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(
    name   = "Group",   # título de la leyenda
    values = c(         # colores para cada combinación
      "t1.h.domesticada" = "#E6AB02",
      "t2.h.domesticada"  = "#4D9221",
      "t1.m.domesticada"   = "#1F78B4",
      "t2.m.domesticada"    = "#E31A1C",
      "t1.h.silvestre" = "#FF7F00",
      "t2.h.silvestre"  = "#6A3D9A",
      "t1.m.silvestre"   = "#B15928",
      "t2.m.silvestre"    = "#33A02C"
    ),
    labels = c(         # nombres que aparecen en la leyenda
      "t1.h.domesticada" = "UDP",
      "t2.h.domesticada"  = "VDP",
      "t1.m.domesticada"   = "UDS",
      "t2.m.domesticada"    = "VDS",
      "t1.h.silvestre" = "UWP",
      "t2.h.silvestre"  = "VWP",
      "t1.m.silvestre"   = "UWS",
      "t2.m.silvestre"    = "VWS"
    )
  ) + 
  scale_x_discrete(labels=c("t1.h.domesticada" = "UDP",
                            "t2.h.domesticada"  = "VDP",
                            "t1.m.domesticada"   = "UDS",
                            "t2.m.domesticada"    = "VDS",
                            "t1.h.silvestre" = "UWP",
                            "t2.h.silvestre"  = "VWP",
                            "t1.m.silvestre"   = "UWS",
                            "t2.m.silvestre"    = "VWS"))

sac1
## Gráfica de concentración total de azucar
total = ggplot(azu2, aes(x = trat,
                         y = Suma.total,
                         color = trat))+
  geom_boxplot()

total1 = ggplot(resultado, aes(x = interaction(trat, sexo, condicion), 
                               y = tot_m,
                               color = interaction(trat, sexo, condicion))) + #fill = condicion
  geom_col(position = "dodge", width = 0.7,
           fill = NA, linewidth = .8) +
  geom_errorbar(aes(ymin = tot_m - tot_sd, 
                    ymax = tot_m + tot_sd),
                width = 0.25) +
  labs(
    #title = "Azúcares totales por grupo experimental",
    x = "",
    y = "Total Sugar Concentration (w/v)",
    #fill = "Condición"
  ) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(
    name   = "Groups",   # título de la leyenda
    values = c(         # colores para cada combinación
      "t1.h.domesticada" = "#E6AB02",
      "t2.h.domesticada"  = "#4D9221",
      "t1.m.domesticada"   = "#1F78B4",
      "t2.m.domesticada"    = "#E31A1C",
      "t1.h.silvestre" = "#FF7F00",
      "t2.h.silvestre"  = "#6A3D9A",
      "t1.m.silvestre"   = "#B15928",
      "t2.m.silvestre"    = "#33A02C"
    ),
    labels = c(         # nombres que aparecen en la leyenda
      "t1.h.domesticada" = "UDP",
      "t2.h.domesticada"  = "VDP",
      "t1.m.domesticada"   = "UDS",
      "t2.m.domesticada"    = "VDS",
      "t1.h.silvestre" = "UWP",
      "t2.h.silvestre"  = "VWP",
      "t1.m.silvestre"   = "UWS",
      "t2.m.silvestre"    = "VWS"
    )
  ) + 
  scale_x_discrete(labels=c("t1.h.domesticada" = "UDP",
                            "t2.h.domesticada"  = "VDP",
                            "t1.m.domesticada"   = "UDS",
                            "t2.m.domesticada"    = "VDS",
                            "t1.h.silvestre" = "UWP",
                            "t2.h.silvestre"  = "VWP",
                            "t1.m.silvestre"   = "UWS",
                            "t2.m.silvestre"    = "VWS"))

total1

library(gridExtra)
library(cowplot)
library(ggpubr)

# joining plots (2 plots by spp)
gt <- ggarrange(fru1, glu1, sac1, total1,
                ncol = 2, nrow = 2, common.legend = TRUE, legend = "right",
                labels="AUTO")
gt

gt <- ggarrange(fru1, glu1, sac1, total1,
                ncol = 2, nrow = 2, legend = "none",
                labels="AUTO")
gt


