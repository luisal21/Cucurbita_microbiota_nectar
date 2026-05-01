############ ASIGNACIÓN TAXONÓMICA CON BASE UNITE ALL COMPLETA ################
#### Analysis post-taxonomic assignation
### change directory
setwd("/home/luis/Documentos/papers_publicados/paper2_tesisDoctorado/Cucurbita_microbiota_nectar/data_its")

# Cargando paquete phyloseq y paquetes adicionales
library(phyloseq); packageVersion("phyloseq")
library(ggplot2)
library(picante)
library(ape)
library("microbiome")
library("dplyr")
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
#library(minpack.lm) # no instalado
library(Hmisc)
library(stats4)

# cargar primero el archivo de metadatos
meta_its <- read.csv("its_metadata_2024.csv", header = T)
head(meta_its)
# despues cargar la tabla de ASV
seqtab.its <- read.csv("asv_table_2026.csv", header = T, row.names = 1)
head(seqtab.its)
tail(seqtab.its)
str(seqtab.its)
# Debido a que al cargar el archivo agrega una X al principio del nombre de
# cada columna, hay que decirle que los nombres de las columnas son iguales
# a las filas de sample_id de metadatos para corregir este error
colnames(seqtab.its) = meta_its$sample.id
str(seqtab.its)
# Volver a cargar archivo de metadatos
meta_its2 <- read.csv("its_metadata_2024.csv", header = T, row.names = 1)
head(meta_its2)
# cargando taxonomia
taxo_its <-read.csv("taxa_fungi_all_abril2026_completa_all.csv", header = T,row.names = 1)
head(taxo_its)
# corrobando match entre archivos
setdiff(rownames(taxo_its), rownames(seqtab.its))
all(rownames(seqtab.its) == rownames(taxo_its))
all(colnames(seqtab.its) == rownames(meta_its2))

# Construct phyloseq object (straightforward from dada2 outputs)
its1 <- phyloseq(otu_table(as.matrix(seqtab.its), taxa_are_rows=T),
                 tax_table(as.matrix(taxo_its)), 
                 sample_data(meta_its2))
its1

# obteniendo el arbol filogenetico
#arbol_its = rtree(ntaxa(its1), rooted=TRUE, tip.label=taxa_names(its1))
#plot(arbol_its)
# Uniendo el arbol con el objeto phyloseq anterior
#its2 = merge_phyloseq(its1, arbol_its)
its2 = its1
# Primero resumimos el contenido del objeto phyloseq
summarize_phyloseq(its2) # reads totales 9,027,900, 1257 ASVs in 364 samples
# información del número de reads
summary(sample_sums(its2)) #range of 0 to 53465

# A useful next step is to explore feature prevalence in the dataset,
# which we will define here as the number of samples in which a taxa appears
# at least once.
# Compute prevalence of each feature, store as data.frame
#prevdf_its = apply(X = otu_table(its2),
#                   MARGIN = ifelse(taxa_are_rows(its2), yes = 1, no = 2),
#                   FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
#prevdf_its = data.frame(Prevalence = prevdf_its,
#                        TotalAbundance = taxa_sums(its2),
#                        tax_table(its2))
# Revisando prevalencia por phylum
#plyr::ddply(prevdf_its, "Phylum", function(df1)
#{cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Filtering phyla whose prevalence is 1 (only present in one sample)
#filterPhyla_its = c(,"p__Chytridiomycota", "p__Mortierellomycota",
#                    "p__Rozellomycota")

# Filter entries with low prevalence Phylum.
#its3 = subset_taxa(its2, !Phylum %in% filterPhyla_its)
#its3 # Final ASV data set: 1294 ASVs in 364 samples


# Filtering samples of Morelia site
its4 = subset_samples(its2, sitio == "morelia")
summarize_phyloseq(its4)  #5,538,466 reads, 1294 ASVs de 236 samples
##### Extraemos todas aquellas secuencias del reino Fungi #####
its4 = subset_taxa(its4, Kingdom == "k__Fungi")
its4
# Hay algunas muestras que quedaron sin secuencias, ya que todas sus secuencias
# pertenecian al Reino Viridiplantae. Por lo que voy a excluir esas 
# muestras que quedaron con 0 secuencias:
its5 = prune_samples(sample_sums(its4)>1, its4)
its5 # 1294 taxa and 235 samples

#Remove taxa with 0 observation
its6 <- prune_taxa(taxa_sums(its5) > 0, its5)
its6 #235 ASVs in 186 samples
summarize_phyloseq(its6) #24188 reads 
sort(sample_sums(its6))

# información del número de reads
summary(sample_sums(its6))
range(sample_sums(its6))
#Library size spans 2 - 3460 reads 

# Obteniendo el número de read por muestra
reads_sample = readcount(its6)
reads_sample[1:5]
#write.csv(reads_sample, file = "reads_sample_its_todas_2026.csv")

#### grafica de abundancia con phyloseq
# creamos una primera visualización de la abundancia total de secuencias por
# muestra a nivel de Phylum
plot_bar(its6, fill = "Phylum") +
  geom_bar(aes(color = Phylum), stat = "identity",
           position = "stack")


#### CONTAR TAXONES POR NIVEL TAXONÓMICO #####
# Extraer la tabla de taxonomía
taxo <- as.data.frame(tax_table(its6))

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

#############################################################################
###############################################################################
###############################################################################
### ANÁLISIS DE NÚMEROS DE HILL DE HONGOS

#install.packages("iNEXT")       # rarefacción y extrapolación
#install.packages("hillR")       # Hill a partir de phyloseq

library(hillR)
library(phyloseq)

# Extraer tabla de OTU/ASV (muestras en filas)
asv_table <- as.data.frame(t(otu_table(its6)))

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

View(hill_df)

# Agregar variables del experimento
meta <- as.data.frame(sample_data(its6))
hill_df <- merge(hill_df, meta, by = "row.names")
View(hill_df)

hist(hill_df$q0)
hist(hill_df$q1)
hist(hill_df$q2)

# Paso 2: Análisis estadístico con tus variables
# Para cada orden de Hill
mod_q0 <- glm(q0 ~ condicion * sexo * tratamiento, 
              data = hill_df,
              family = poisson())
summary(mod_q0)
# Calcular ratio de dispersión (debe ser ~1 si no hay sobredispersión)
E2 <- resid(mod_q0, type = "pearson")
sum(E2^2) / df.residual(mod_q0)

library(MASS)
# q0 modelo bueno
mod_q0 <- glm.nb(q0 ~ condicion * sexo * tratamiento, 
              data = hill_df,
              link = "log")
summary(mod_q0)

#q1
mod_q1 <- glm(q1 ~ condicion * sexo * tratamiento, 
              data = hill_df,
              family = Gamma(link = "log"))

summary(mod_q1)

# q2
mod_q2 <- glm(q2 ~ condicion * sexo * tratamiento, 
              data = hill_df,
              family = Gamma(link = "log"))
summary(mod_q2)


# ANOVA tipo II
library(car)
Anova(mod_q0)
Anova(mod_q1)
Anova(mod_q2)

# Summary
#q1
tapply(hill_df$q1, list(hill_df$condicion,hill_df$sexo,
                        hill_df$tratamiento), mean)
tapply(hill_df$q1, list(hill_df$condicion,hill_df$sexo,
                        hill_df$tratamiento), sd)
#q2
tapply(hill_df$q2, list(hill_df$condicion,hill_df$sexo,
                        hill_df$tratamiento), mean)
tapply(hill_df$q2, list(hill_df$condicion,hill_df$sexo,
                        hill_df$tratamiento), sd)

###   CURVAS DE RAREFACCIÓN
### POR CONDICIÓN
library(iNEXT)

# Extraer tabla de ASVs y metadatos del objeto phyloseq
asv_table <- as.data.frame(otu_table(its6))# ASVs en filas, muestras en columnas
metadata <- as.data.frame(as.matrix(sample_data(its6)))

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

coverage_incidence <- estimateD(asv_table,
                                    q = c(0, 1, 2),
                                    datatype = "abundance",
                                    base = "coverage",
                                    conf = 0.95, nboot = 100)
coverage_incidence
warnings()
# Unir con metadatos
hill_df2 <- data.frame(
  sample_id  = names(hill_q0),
  q0         = hill_q0,
  q1         = hill_q1,
  q2         = hill_q2
)

View(hill_df2)





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
###############################################################################
###  differential abundance (DA) analysis of microbial taxa 
#devtools::install_github("FrederickHuangLin/ANCOMBC")
### ANALISIS CON ANCOM-BC2
library(ANCOMBC)

# Condición
ancom_cond_its <- ancombc2(
  data = its6,
  tax_level = "Family",
  fix_formula = "condicion",
  p_adj_method = "BH",
  alpha = 0.05,
  prv_cut = 0.10
)

res_cond_its <- ancom_cond_its$res
res_cond_its_sig <- res_cond_its[res_cond_its$diff_condicionsilvestre == TRUE |
                                   res_cond_its$diff_condicionsilvestre == TRUE,
                                 c("taxon", "lfc_condicionsilvestre",
                                   "q_condicionsilvestre")]
print(res_cond_its_sig)

# Sexo
ancom_sexo_its <- ancombc2(
  data = its6,
  tax_level = "Family",
  fix_formula = "sexo",
  p_adj_method = "BH",
  alpha = 0.05,
  prv_cut = 0.10
)

res_sexo_its <- ancom_sexo_its$res
res_sexo_its_sig <- res_sexo_its[res_sexo_its$diff_sexomacho == TRUE |
                                   res_sexo_its$diff_sexomacho == TRUE,
                                 c("taxon", "lfc_sexomacho",
                                   "q_sexomacho")]
print(res_sexo_its_sig)
# Tratamiento
ancom_trat_its <- ancombc2(
  data = its6,
  tax_level = "Genus",
  fix_formula = "tratamiento",
  p_adj_method = "BH",
  alpha = 0.05,
  prv_cut = 0.10
)

# Ver resultados significativos
res_trat_its <- ancom_trat_its$res
res_trat_its_sig <- res_trat_its[res_trat_its$diff_tratamientot2 == TRUE, 
                                 c("taxon", "lfc_tratamientot2", 
                                   "q_tratamientot2")]
print(res_trat_its_sig)




############################################################################
### ANÁLISIS CON DESEQ2
# Paso 1. Instalar y cargar paquetes
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)
library(phyloseq)

# 2. Convertir phyloseq a DESeq2
# primero revisar cuantos NA hay
sum(is.na(tax_table(its6)[, "Family"]))
# Eliminar ASVs sin asignación a nivel de Familia y agrupar 
its_familia <- its6%>%
  subset_taxa(!is.na(Family)) %>%
  tax_glom(taxrank = "Family")
# Para la comparación por condición (domesticated vs wild)
ds_condicion_its <- phyloseq_to_deseq2(its_familia, ~ condicion)


# Método alternativo
geoMeans <- apply(counts(ds_condicion_its), 1, function(row) {
  if (all(row == 0)) 0 
  else exp(mean(log(row[row != 0])))
})

ds_condicion_its <- estimateSizeFactors(ds_condicion_its, 
                                      geoMeans = geoMeans)
        



# 3. Estimar factores de tamaño con conteos cero
ds_condicion_its <- estimateSizeFactors(ds_condicion_its, 
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




###############################################################################
### Código para Top 10 Familias Global (Sin Agrupar por Variables)
# A. Aglomerar a nivel de Familia
ps_fam1 <- tax_glom(its6, taxrank = "Family", NArm = FALSE)

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
write.csv(top10_global, file = "top10_familias_hongos.csv")

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
###############################################################################
###############################################################################
################ BETA DIVERSITY OF BACTERIA
BiocManager::install("DECIPHER")

# Cargando paquetes
library(DECIPHER)
library(phangorn)
library(Biostrings)


ntaxa(its6)

# 1. extraer directo de los rownames
seqs <- rownames(otu_table(its6))

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


##### Importar árbol ####
library(ape)
arbol_ft <- read.tree("arbol_fasttree.nwk")

# Verificar
arbol_ft
plot(arbol_ft, cex = 0.2)  # grafica rápida para verificar

# Incorporar a phyloseq
ps_con_arbol <- merge_phyloseq(its6, phy_tree(arbol_ft))

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
                                    shape = "sexo",
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

disp_bray2 <- betadisper(dist_bray, metadata$sexo)
permutest(disp_bray2, permutations = 999)

disp_bray3 <- betadisper(dist_bray, metadata$tratamiento)
permutest(disp_bray3, permutations = 999)

# Sørensen
disp_sor1 <- betadisper(dist_sorensen, metadata$condicion)
permutest(disp_sor1, permutations = 999)

disp_sor2 <- betadisper(dist_sorensen, metadata$sexo)
permutest(disp_sor2, permutations = 999)

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

