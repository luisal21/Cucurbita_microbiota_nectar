############## Taxonomic and diversity analysis with R of 16s amplicon
############## sequencing from nectar samples.
##### Luis Alberto Villanueva Espino
#### November 2026

#### Loading packages ####
library(phyloseq)
library(ggplot2)
library(picante)
library(ape)
library("microbiome")
library("dplyr")
library("RColorBrewer")
library("ggpubr")
library(vegan)
library(gridExtra)
library(tidyr)
library(reshape2)
library(Hmisc)
library(stats4)

#### Set up directory in which the databases are located ####
setwd("/home/luis/Documentos/papers_publicados/paper2_tesisDoctorado/Cucurbita_microbiota_nectar/data_16s")

### First, import the original table "seqtab_nochim_16s_2024.csv, then
# transpose the table. NOTE: Disable this code at the end
#asvtable <- read.csv("seqtab_nochim_16s_2024.csv", header = T, row.names = 1)
#head(asvtable, 1)
#asv2 = t(asvtable)
#head(asv2, 1)
#write.csv(asv2, file = "asv_table.csv")

# Load metadata file:
meta_16s <- read.csv("meta_16s_2024.csv", header = T)
head(meta_16s)

# Then, Load the new ASV table with the rows and columns tranposed:
seqtab.16s <- read.csv("asv_table.csv", header = T, row.names = 1)
head(seqtab.16s)
str(seqtab.16s)

# Because uploading the file adds an “X” to the beginning of each row's name,
# you need to tell it that the column names are the same as the “sample_id” 
# rows in the metadata to correct this error: 
colnames(seqtab.16s) = meta_16s$sample_id
str(seqtab.16s)

# Load again metadata file:
meta_16s2 <- read.csv("meta_16s_2024.csv", header = T, row.names = 1)
head(meta_16s2)

# Now, load taxonomy file:
taxo_16s <-read.csv("taxa_16s_species_2024.csv", header = T,row.names = 1)
head(taxo_16s)

# Verifying file matches:
setdiff(rownames(taxo_16s), rownames(seqtab.16s))
all(rownames(seqtab.16s) == rownames(taxo_16s))
all(colnames(seqtab.16s) == rownames(meta_16s2))

# Construct phyloseq object (straightforward from dada2 outputs):
s16_1 <- phyloseq(otu_table(as.matrix(seqtab.16s), taxa_are_rows=T),
                  tax_table(as.matrix(taxo_16s)), 
                  sample_data(meta_16s2))

# Check phyloseq object:
s16_1

# Check summary with microbiome package
microbiome::summarize_phyloseq(s16_1) #13,732,428 reads, 4160 ASVs

# deleting bad quality samples with low reads
s16_1 = subset_samples(s16_1, !sampleid2 %in% c("281022-EMo22-SA014",
                      "281022-EMo22-SA035", "281022-EMo22-SD064",
                      "281022-EMo22-SD085", "EMo1223S-SC273",
                      "ENESMo24-bLVEv416s82"))

# Check summary after removing bad quality samples
summarize_phyloseq(s16_1) #13,731,996 reads, 4160 ASVs in 356 samples
                          
### Remove contaminants using the prevalence method of the package decontam
#BiocManager::install("decontam")
library(decontam)

head(sample_data(s16_1))
df <- as.data.frame(sample_data(s16_1)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(s16_1)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=sample_or_control)) + geom_point()

# Identify Contaminants -Frequency
contamdf.freq <- isContaminant(s16_1, method="frequency", conc="conce_mg")
head(contamdf.freq)
#View(contamdf.freq)
table(contamdf.freq$contaminant)
head(which(contamdf.freq$contaminant))

# Identify Contaminants -Prevalence
sample_data(s16_1)$is.neg <- sample_data(s16_1)$sample_or_control == "Control Sample"
contamdf.prev <- isContaminant(s16_1, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant), 20)
#View(contamdf.prev)

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(s16_1, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$sample_or_control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$sample_or_control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

### Remove this contaminant ASVs
s16_2 <- prune_taxa(!contamdf.prev$contaminant, s16_1)
summarize_phyloseq(s16_2) #13,343,174 reads, 4112 ASVs in 356 samples
sum(sample_sums(s16_1)) - sum(sample_sums(s16_2))
# 388822 read identify as contaminants of 48 ASVs

################### Extracting only contaminants ##############################
contaminantes <- prune_taxa(contamdf.prev$contaminant, s16_1) 
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
prevdf = apply(X = otu_table(s16_2),
               MARGIN = ifelse(taxa_are_rows(s16_2), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(s16_2),
                    tax_table(s16_2))
#View(prevdf)

plyr::ddply(prevdf, "Phylum", function(df1)
{cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Filtering phyla whose prevalence is 1
filterPhyla = c("Euryarchaeota", "FCPU426", "RCP2-54",
                "Synergistota", "Thermotogota", "WPS-2")
# Filter entries with low prevalence Phylum.
s16_3 = subset_taxa(s16_2, !Phylum %in% filterPhyla)
s16_3 # Final ASV data set: 4106 ASVs in 356 samples

# Filtering samples of Morelia site
s16_4 = subset_samples(s16_3, sitio=="morelia")
summarize_phyloseq(s16_4)  #9,453,658 reads, 4106 ASVs in 231 samples

# Removing Chloroplast and Mitochondria sequences and na samples in condition
s16_5 = subset_samples(s16_4, condicion !="BLANCO")
s16_5 = subset_taxa(s16_5, Order != "Chloroplast")
s16_5 = subset_taxa(s16_5, Family != "Mitochondria")

# There are some samples that ended up with no sequences, since all of their
# sequences belonged to the Viridiplantae kingdom. Therefore, I am going to 
# exclude those samples that ended up with 0 sequences:
s16_5 = prune_samples(sample_sums(s16_5)>1, s16_5)
s16_5

#Remove taxa with 0 observation
s16_6 <- prune_taxa(taxa_sums(s16_5) > 0, s16_5)
s16_6
View(otu_table(s16_6))

# Check summary of phyloseq object
summarize_phyloseq(s16_6) #146,201 reads, 1635 ASVs in 227 samples

# Summary of number of read
summary(sample_sums(s16_6))
range(sample_sums(s16_6)) #Library size spans 2 - 24079 reads 

# Calculating the number of reads per sample
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


#### COUNT TAXONS BY TAXONOMIC LEVEL #####
# Extract the taxonomy table
taxo <- as.data.frame(tax_table(s16_6))

# Count taxons by taxonomic level
n_species <- length(unique(na.omit(taxo$Species)))
n_genera  <- length(unique(na.omit(taxo$Genus)))
n_families <- length(unique(na.omit(taxo$Family)))
n_orders  <- length(unique(na.omit(taxo$Order)))
n_classes <- length(unique(na.omit(taxo$Class)))
n_phyla   <- length(unique(na.omit(taxo$Phylum)))

# Make a table
data.frame(
  Nivel   = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
  N_taxones = c(n_phyla, n_classes, n_orders, n_families, n_genera, n_species)
)

###############################################################################
###############################################################################
### Analysis of alpha diversity calculating Hill numbers

# Extract data from phyloseq object:
otu_data <- otu_table(s16_6, taxa_are_rows = TRUE)
metadata <- as(sample_data(s16_6), "data.frame")
head(metadata)


# 01. Calculating hill numbers
library(hilldiv3)

q0 <- hilldiv(otu_data, q = 0)
q1 <- hilldiv(otu_data, q = 1)
q2 <- hilldiv(otu_data, q = 2)

# Merge metadata with Hill numbers
library(tidyr)

q012_ancho <- bind_rows(q0, q1, q2) %>%
  pivot_wider(names_from = q,
              values_from = value,
              names_prefix = "q")

head(q012_ancho)

# join
metadata_with_hill <- q012_ancho %>%
  inner_join(metadata, by = c("sample"="sampleid2"))

# show
head(metadata_with_hill)

# Save table
#write.csv(metadata_with_hill, file = "metadata_with_hill.csv")


### 02. Get Hill means
# Reorder q values
meta_qs <- metadata_with_hill %>%
  pivot_longer(cols = q0:q2, names_to = "q", values_to = "value") %>%
  filter(q %in% c("q0", "q1", "q2")) %>%
  mutate(
    qs = case_when(
      q == "q0" ~ "q0=Observed",
      q == "q1" ~ "q1=Exp Shannon",
      q == "q2" ~ "q2=Inv Simpson",
    ))

#Get means of Hill numbers
means <- meta_qs %>% group_by(condicion, sexo, tratamiento, qs) %>%
  summarise(means = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE),
            .groups = 'drop')

print(means)
#write.csv(means, file = "Hill_means.csv")

### 03. Get sequencing effort
### Sample effort
library(vegan)
mat <- as(t(otu_table(s16_6)), "matrix")
raremax <- min(rowSums(mat))

# plot
system.time(rarecurve(mat, step = 100, sample = raremax, col = "green4", 
                      label = FALSE))


### 04. Check alpha diversity correlation with sequencing depth
#Get effective reads
Effective_reads <- colSums(otu_data)
metadata_with_hill$Effective_reads <- Effective_reads
View(metadata_with_hill)

# Alpha diversity depth correlation to order q=0
library(ggpubr)

q0_vs_depth <- ggscatter(metadata_with_hill,
               x = "Effective_reads", y = "q0", 
               xlab= "Sequencing depth (number of reads)",
               ylab = "q=0", add = "reg.line",
               add.params = list(color = "#114e9d", fill = "lightgray"),
               conf.int = TRUE, # Add confidence interval
               cor.coef = TRUE, # Add correlation coefficient. 
               cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)+
  theme(legend.title = element_blank(), legend.position = "none")

# Alpha diversity depth correlation to order q=1
q1_vs_depth <- ggscatter(metadata_with_hill, 
               x = "Effective_reads", y = "q1", 
               xlab= "Sequencing depth (number of reads)",
               ylab = "q=1", add = "reg.line",
               add.params = list(color = "#114e9d", fill = "lightgray"),
               conf.int = TRUE, # Add confidence interval
               cor.coef = TRUE, # Add correlation coefficient
               cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
  theme(legend.title = element_blank(), 
        legend.position = "none")

# Alpha diversity depth correlation to order q=2
q2_vs_depth <- ggscatter(metadata_with_hill, 
               x = "Effective_reads", y = "q2", 
               xlab= "Sequencing depth (number of reads)",
               ylab = "q=2", add = "reg.line",
               add.params = list(color = "#114e9d", fill = "lightgray"),
               conf.int = TRUE, # Add confidence interval
               cor.coef = TRUE, # Add correlation coefficient
               cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
  theme(legend.title = element_blank(), 
        legend.position = "none")

#Combine plot
library(cowplot)
title_corr_plot <- ggdraw() + 
  draw_label("Alpha diversity depth correlation")

correlation_plot_q012 <- plot_grid(title_corr_plot,
                                   q0_vs_depth, 
                                   q1_vs_depth, 
                                   q2_vs_depth,
                                   ncol = 1, 
                                   rel_heights = c(0.15, 1, 1, 1))
# Print plot
correlation_plot_q012


### 05. Statistical analysis of Hill numbers
# histograms of Hill numbers
hist(metadata_with_hill$q0)
hist(metadata_with_hill$q1)
hist(metadata_with_hill$q2)


# Generalized linear model of q0 with negative binomial distribution
library(MASS)
q0glm <- glm.nb(q0 ~ log(Effective_reads) + condicion * sexo * tratamiento,
                data = metadata_with_hill)

# summary
summary(q0glm)

#effects
library(car)
Anova(q0glm)

# check model reliability
library(DHARMa)

sim <- simulateResiduals(q0glm)
plot(sim)              # QQ-plot of uniformity and residuals versus predicted 
# values
testDispersion(sim)    # Is there overdispersion or underdispersion?
testZeroInflation(sim) # Are there more zeros than expected?
testOutliers(sim)      # Are there any atypical ones?

#get differences with emmeans
library(emmeans)
emm <- emmeans(q0glm, ~ condicion | sexo, type = "response")
emm                               # estimated means by group

pairs(emm)                        # condition within each sex

# plot
emmip(q0glm, condicion ~ sexo, type = "response", CIs = TRUE)


### q1 Diversity Exp Shannon
mod_q1 <- glm(q1 ~ log(Effective_reads) + condicion * sexo * tratamiento, 
              data = metadata_with_hill,
              family = Gamma(link = "log"))

# summary
summary(mod_q1)

#effects
Anova(mod_q1)

### q2 Dominance Inv Simpson
mod_q2 <- glm(q2 ~ condicion * sexo * tratamiento, 
              data = metadata_with_hill,
              family = Gamma(link = "log"))

# summary
summary(mod_q2)

#effects
Anova(mod_q2)

###############################################################################
###############################################################################
################ BETA DIVERSITY OF BACTERIA
#BiocManager::install("DECIPHER")
  
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



####################### Importar árbol ##########################
library(ape)
library(phangorn)
setwd("/home/luis/Documentos/papers_publicados/paper2_tesisDoctorado/Cucurbita_microbiota_nectar/data")

# Cargando arbol
arbol_ft <- read.tree("arbol_its2.nwk")

# Verificar
arbol_ft
#plot(arbol_ft, cex = 0.2)  # grafica rápida para verificar
# Enraizando arbol
arbol_ft_root = midpoint(arbol_ft)
is.rooted(arbol_ft_root)

# Incorporar a phyloseq
ps_con_arbol <- merge_phyloseq(s16_6, phy_tree(arbol_ft_root))

# Verificar que quedó bien
ps_con_arbol
phy_tree(ps_con_arbol)

######### UniFrac no ponderado #########
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

### Plot separated by treatment
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
baseplot_16s_wei
pcoa_16s = baseplot_16s_wei + stat_ellipse(linetype = 2)
pcoa_16s +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key.size = unit(1, "cm")
  )
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

### Plot separated by condition
baseplot_16s_wei <- plot_ordination(ps_con_arbol, 
                                    ordination_16s_wei,
                                    color="condicion",
                                    #shape = "sexo",
                                    type = "samples") + #shape = "sexo" shape = "tratamiento"
  ggtitle("")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 4)+
  scale_color_discrete(name = "Condition",
                       labels = c("domesticada" = "Domesticated",
                                  "silvestre" = "Wild"),
                       palette = c("#E6AB02","#4D9221"))

baseplot_16s_wei
pcoa_16s = baseplot_16s_wei + stat_ellipse(linetype = 2)
pcoa_16s +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key.size = unit(1, "cm")
  )

#scale_shape_manual(values=c(16,3))
#labs(x = "PCoa1 (31.3%)", y = "PCoA2 (10.4%)")+
#scale_color_discrete(name = "Condition", 
#                    labels = c("Domesticated","Wild"))+
#scale_color_discrete(name = "Treatment",
#                     labels = c("No-visited","Visited"))

  
#scale_color_discrete(name = "Species", 
#                    labels = c("2CPF" = "CPF","3CAS" = "CAS",
#                              "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"),
#        palette = "Set1")+


#baseplot_16s_wei
#baseplot_16s_hem_wei + geom_polygon(aes(fill = especies))
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
                       "More abundant in wild", 
                       "More abundant in domesticated"),
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


### Tratamiento ###
ancom_trat_16s <- ancombc2(
  data = s16_6,
  tax_level = "Genus",
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
ancom_16s_trat = ggplot(res_trat_16s_plot, 
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


ancom_16s_trat

# Grafica
ancom_16s_trat = ggplot(res_trat_16s_plot, 
                        aes(x = reorder(taxon, lfc_tratamientot2),
                            y = lfc_tratamientot2,
                            fill = direccion)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_fill_manual(values = "tomato") + #"steelblue"
  labs(x = "Genus",
       y = "log2 Fold Change (visited vs. unvisited)",
       fill = NULL) +
  theme_bw() +
  theme(legend.position = "top")


ancom_16s_trat

citation()
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

etiquetas = c("t1" = "Unvisited", "t2" = "Visited")

# gráfica modificada para bacterias
s16_prop %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  comp_barplot(tax_level = "Family", n_taxa = 15)+
  facet_wrap(~ tratamiento, scales = "free",
             labeller = as_labeller(etiquetas))+ #, scales = "free"
  coord_flip() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"))

# Gráfica para hongos
its_prop %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  comp_barplot(tax_level = "Family", n_taxa = 15, facet_by = "tratamiento")+
  coord_flip() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())




# ITS (fungi)
transform <- microbiome::transform
its_prop <- transform(its6, "compositional")
otu_table(its_prop)[1:10, 1:10]
head(sample_sums(its_prop))

its_prop %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  comp_barplot(tax_level = "Family", n_taxa = 15)+
  facet_wrap(~ condicion + sexo + tratamiento, scales = "free")+ #, scales = "free"
  coord_flip() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())



ps_prop %>%
  tax_fix() %>%
  comp_barplot(tax_level = "Family", n_taxa = 15, facet_by = "tratamiento",
               group_by = "")+
  coord_flip() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

#facet_grid(labeller = labeller(tratamiento = trata))

###############################################################################
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
#write.csv(top10_global, file = "top10_familias_bacterias.csv")

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
### Usando codigo de metagenomics carprentries
percentages <- transform_sample_counts(s16_6, function(x) x*100 / sum(x) )
head(percentages@otu_table@.Data)


percentages_glom <- tax_glom(percentages, taxrank = 'Genus')

View(percentages_glom@tax_table@.Data)
write.csv(percentages_glom@tax_table@.Data, "genus_taxtable.csv")

percentages_df <- psmelt(percentages_glom)
str(percentages_df)


percentages_df$Family <- as.factor(percentages_df$Family)

family_colors_rel<- colorRampPalette(brewer.pal(8,"Dark2")) (length(levels(percentages_df$Family)))

relative_plot <- ggplot(data=percentages_df, aes(x=Sample, y=Abundance,
                                                 fill=Family))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  scale_fill_manual(values = family_colors_rel)

relative_plot

# Agrupando las familias menos abundantes
percentages_df$Family <- as.character(percentages_df$Family) # Return the Phylum column to be of type character
percentages_df$Family[percentages_df$Abundance < 30] <- "Families < 30% abund."
unique(percentages_df$Family)

# Gráfica a nivel de familia
percentages_df$Family <- as.factor(percentages_df$Family)

family_colors_rel<- colorRampPalette(brewer.pal(8,
                    "Dark2")) (length(levels(percentages_df$Family)))

relative_plot <- ggplot(data=percentages_df, aes(x=Sample, y=Abundance,
                                                 fill=Family))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  scale_fill_manual(values = family_colors_rel)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
relative_plot

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
View(more2021)
more2021$ID.planta
more2021$codigo.video
sort(more2021$Hora.visita.real)
sort(more2021$fecha)
more2021$duración.decimal

tapply(more2021$ID.planta, list(more2021$condicion, more2021$sexo), length)


# Cargando librerias
library(fossil)
library (bipartite)
library(dplyr)
library(ggplot2)


### Horas totales de grabación entre 2021 y 2022: 106 horas
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
tapply(more2021$especie, list(more2021$condicion,more2021$sexo), length)

tapply(more2021$especie, more2021$Especie.planta, length)

#
visitas <- more2021 %>%
  dplyr::group_by(condicion, Especie.planta, sexo, especie) %>%
  dplyr::summarize(num_visitas = length(especie))

head(visitas)
## frecuencia de visitas
frec <- more2021 %>%
  dplyr::group_by(año, fecha, condicion, Especie.planta, ID.planta, sexo, especie,
                  duración.decimal) %>%
  dplyr::summarize(num_visitas = n()) %>% 
  dplyr::mutate(
    duracion_horas = duración.decimal / 60,  # convertir a horas
    frecu = num_visitas / duracion_horas      # visitas por hora
  )
frec
View(frec)

sum(frec$num_visitas)
tapply(frec$num_visitas, frec$especie, sum)
tapply(frec$num_visitas, list(frec$especie,frec$condicion, frec$sexo), sum)
(587/1145)*100 # apis
(387/1145)*100 # pepo macho
(88/1145)*100 # xenoglossa
(73/1145)*100 # pepo hembra
(7/1145)*100 # halictidae
(2/1145)*100 # xilocopa
(1/1145)*100 # syrphidae
(3/1145)*100 # syrphidae

frec$especie
tapply(frec$frecu, frec$especie, length)
hist(frec$frecu)
sum(frec$frecu == 0)

# contabilizando cuantos datos tengo de cada especie
frec %>%
  group_by(especie) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

# agrupando las 2 especies menos frecuentes en other
frec <- frec %>%
  mutate(especie_agrupada = case_when(
    especie %in% c("Xilocopa_tabaniformis", 
                   "Syrphidae") ~ "Other",
    TRUE ~ as.character(especie)
  ))

# Verificar
frec %>%
  group_by(especie_agrupada) %>%
  summarise(n = n()) %>%
  arrange(desc(n))



### Modelos estadisticos
modvis = glm(frecu ~ condicion * sexo,
    data = frec, family = Gamma(link = "log"))
summary(modvis)

library(car)
Anova(modvis)

modvis2 = glm(frecu ~ condicion * sexo * especie_agrupada,
              data = frec, family = Gamma(link = "log"))
summary(modvis2)
Anova(modvis2)

library(emmeans)
library(multcomp)
#install.packages("multcompView")
library(multcompView)

# Letras de significancia por especie (efecto más fuerte, p < 0.001)
emm <- emmeans(modvis2, ~ especie | sexo,
               type = "response")
pares <- pairs(emm, adjust = "BH")
summary(pares)


emm_esp <- emmeans(modvis2, ~ especie_agrupada, type = "response")
emm_esp
# Comparaciones entre especies dentro de cada grupo
pares <- pairs(emm_esp, adjust = "BH")
summary(pares)

letras <- cld(emm_esp, adjust = "BH", Letters = letters)
print(letras)


######### plot
library(ggplot2)
library(dplyr)

# Calcular media y error estándar por especie agrupada y condición
frec_summary <- frec %>%
  group_by(especie_agrupada, condicion) %>%
  summarise(
    media = mean(frecu),
    se    = sd(frecu) / sqrt(n()),
    .groups = "drop"
  )

# Ordenar por media descendente para mejor visualización
orden_especies <- frec_summary %>%
  group_by(especie_agrupada) %>%
  summarise(media_total = mean(media)) %>%
  arrange(desc(media_total)) %>%
  pull(especie_agrupada)

frec_summary <- frec_summary %>%
  mutate(especie_agrupada = factor(especie_agrupada, 
                                   levels = orden_especies))

# Gráfica
ggplot(frec_summary, aes(x = especie_agrupada,
                         y = media,
                         fill = condicion)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7,
           color = "black",
           linewidth = 0.3) +
  geom_errorbar(aes(ymin = media - se,
                    ymax = media + se),
                position = position_dodge(width = 0.8),
                width = 0.25,
                linewidth = 0.3) +
  scale_fill_manual(values = c("domesticada" = "white",
                               "silvestre"   = "black"),
                    labels = c("domesticada" = "Domesticated",
                               "silvestre"   = "Wild")) +
  scale_x_discrete(labels = c(
    "Peponapis_macho"  = "Peponapis sp. ♂",
    "Peponapis_hembra" = "Peponapis sp. ♀",
    "Xenoglossa"       = "Xenoglossa sp.",
    "Apis_mellifera"   = "Apis mellifera",
    "Halictidae"       = "Halictidae",
    "Other"            = "Other"
  )) +
  labs(x = NULL,
       y = expression("Mean visitation rate (visits hour"^-1*")"),
       fill = NULL) +
  theme_classic() +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1,
                                 face = "italic", size = 9),
    axis.text.y   = element_text(size = 10),
    axis.title.y  = element_text(size = 11),
    legend.position = "top",
    legend.text   = element_text(size = 10)
  )

##########
# Calcular media y error estándar por especie agrupada, condición y sexo
frec_summary <- frec %>%
  group_by(especie_agrupada, condicion, sexo) %>%
  summarise(
    media = mean(frecu),
    se    = sd(frecu) / sqrt(n()),
    .groups = "drop"
  )

# Ordenar por media descendente
orden_especies <- frec_summary %>%
  group_by(especie_agrupada) %>%
  summarise(media_total = mean(media)) %>%
  arrange(desc(media_total)) %>%
  pull(especie_agrupada)

frec_summary <- frec_summary %>%
  mutate(especie_agrupada = factor(especie_agrupada,
                                   levels = orden_especies))

# Gráfica con facet por sexo
ggplot(frec_summary, aes(x = especie_agrupada,
                         y = media,
                         fill = condicion)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7,
           color = "black",
           linewidth = 0.3) +
  geom_errorbar(aes(ymin = media - se,
                    ymax = media + se),
                position = position_dodge(width = 0.8),
                width = 0.25,
                linewidth = 0.3) +
  scale_fill_manual(values = c("domesticada" = "white",
                               "silvestre"   = "black"),
                    labels = c("domesticada" = "Domesticated",
                               "silvestre"   = "Wild")) +
  scale_x_discrete(labels = c(
    "Peponapis_macho"  = "Peponapis sp. ♂",
    "Peponapis_hembra" = "Peponapis sp. ♀",
    "Xenoglossa"       = "Xenoglossa sp.",
    "Apis_mellifera"   = "Apis mellifera",
    "Halictidae"       = "Halictidae",
    "Other"            = "Other"
  )) +
  facet_wrap(~ sexo,
             ncol = 1,
             labeller = labeller(sexo = c(
               "H" = "Pistillate flowers",
               "M" = "Staminate flowers"
             ))) +
  labs(x = NULL,
       y = expression("Mean visitation rate (visits hour"^-1*")"),
       fill = NULL) +
  theme_classic() +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1,
                                 face = "italic", size = 9),
    axis.text.y   = element_text(size = 10),
    axis.title.y  = element_text(size = 11),
    legend.position = "top",
    legend.text   = element_text(size = 10),
    strip.text    = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = "black")
  )

############################
library(ggplot2)
library(dplyr)

# Calcular media y error estándar por especie
frec_summary_esp <- frec %>%
  group_by(especie_agrupada) %>%
  summarise(
    media = mean(frecu),
    se    = sd(frecu) / sqrt(n()),
    .groups = "drop"
  )

# Agregar letras del cld
letras_df <- as.data.frame(letras) %>%
  dplyr::select(especie_agrupada, .group) %>%
  mutate(.group = trimws(.group))

frec_summary_esp <- left_join(frec_summary_esp, letras_df, by = "especie_agrupada")

# Ordenar por media descendente
frec_summary_esp <- frec_summary_esp %>%
  mutate(especie_agrupada = factor(especie_agrupada, 
                          levels = especie_agrupada[order(media, decreasing = FALSE)]))

# Gráfica
ggplot(frec_summary_esp, aes(x = especie_agrupada, y = media)) +
  geom_col(width = 0.7,
           color = "black",
           fill  = "gray80",
           linewidth = 0.3) +
  geom_errorbar(aes(ymin = media - se,
                    ymax = media + se),
                width = 0.25,
                linewidth = 0.3) +
  geom_text(aes(y = media + se + 0.3,
                label = .group),
            size = 4) +
  scale_x_discrete(labels = c(
    "Peponapis_macho"      = "Peponapis sp. ♂",
    "Peponapis_hembra"     = "Peponapis sp. ♀",
    "Xenoglossa"           = "Xenoglossa sp.",
    "Apis_mellifera"       = "Apis mellifera",
    "Halictidae"           = "Halictidae",
    "Other"                = "Other"
  )) +
  coord_flip() +
  labs(x = NULL,
       y = expression("Mean visitation rate (visits hour"^-1*")")) +
  theme_classic() +
  theme(
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(face = "italic", size = 9),
    axis.title.x = element_text(size = 11)
  )


###########3
frec_summary_int <- frec %>%
  group_by(especie_agrupada, condicion, sexo) %>%
  summarise(
    media = mean(frecu),
    se    = sd(frecu) / sqrt(n()),
    .groups = "drop"
  )

ggplot(frec_summary_int, 
       aes(x = reorder(especie_agrupada, media),
           y = media,
           fill = condicion)) +
  geom_col(position = position_dodge(0.8),
           width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = media - se, ymax = media + se),
                position = position_dodge(0.8),
                width = 0.25, linewidth = 0.3) +
  scale_fill_manual(values = c("domesticada" = "white",
                               "silvestre"   = "black"),
                    labels = c("domesticada" = "Domesticated",
                               "silvestre"   = "Wild")) +
  facet_wrap(~ sexo, ncol = 1,
             labeller = labeller(sexo = c("H" = "Pistillate",
                                          "M" = "Staminate"))) +
  coord_flip() +
  labs(x = NULL,
       y = expression("Mean visitation rate (visits hour"^-1*")"),
       fill = NULL) +
  theme_classic() +
  theme(
    axis.text.y  = element_text(face = "italic"),
    legend.position = "top",
    strip.text = element_text(face = "bold")
  )
##############################################################################
### Promediando la frecuencia de visita para la red
frec2 <- frec %>%
  dplyr::group_by(condicion, Especie.planta, sexo, especie_agrupada) %>%
  dplyr::summarize(mean_frec = mean(frecu))
frec2

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

library(bipartite)
library(fossil)
library(tidyr)
### Realizando redes de interacción ###
visitas2 = visitas %>%
  unite("plant_sex", Especie.planta:sexo)
visitas2

df_vis <- create.matrix(as.data.frame(visitas2), tax.name = "especie",
                        locality = "plant_sex",
                        abund.col = "num_visitas",
                        abund = TRUE)
df1_vis <- t(df_vis)

# frecuencias de visita
frec3 = frec2 %>%
  unite("plant_sex", condicion:sexo)
#write.csv(frec3, file = "frecuencias_visitas.csv")
df <- create.matrix(as.data.frame(frec3), tax.name = "especie_agrupada",
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
plotweb(df1, sorting = "dec", higher_color = "#E6AB02",
        lower_color = "#4D9221")

# Orden de nombres
colnames(df1)
rownames(df1)

# renombrando etiquetas
mapeo_especies <- c(
  "Apis_mellifera"   = "Apis mellifera",
  "Halictidae"       = "Halictidae",
  "Other"            = "Otro",
  "Peponapis_hembra" = "Xenoglossa (Peponapis) ♀",
  "Peponapis_macho"  = "Xenoglossa (Peponapis) ♂",
  "Xenoglossa"       = "Xenoglossa sp."
)

colnames(df1) <- mapeo_especies[colnames(df1)]

mapeo_grupos <- c(
  "domesticada_CAA_H" = "Domesticada CAA ♀",
  "domesticada_CAA_M" = "Domesticada CAA ♂",
  "domesticada_CM_H" = "Domesticada CM ♀",
  "domesticada_CM_M" = "Domesticada CM ♂",
  "domesticada_CPP_H" = "Domesticada CPP ♀",
  "domesticada_CPP_M" = "Domesticada CPP ♂",
  "silvestre_CAS_H" = "Silvestrre CAS ♀",
  "silvestre_CAS_M"   = "Silvestre CAS ♂",
  "silvestre_CPF_M"   = "Silvestre CPF ♂"
)

rownames(df1) <- mapeo_grupos[rownames(df1)]

# Verificando cambio de nombres
colnames(df1)
rownames(df1)


# frecuencia de visita
plotweb(df1, sorting = "dec",
        higher_color = "#964B00",
        lower_color = c("#4D9221","#4D9221","#E6AB02","#E6AB02","#E6AB02",
                        "#E6AB02","#E6AB02","#4D9221","#E6AB02"))



### Metricas básicas de la red
# Calcular métricas de red
networklevel(df1, index = c(
  "connectance",           # proporción de interacciones realizadas vs posibles
  "linkage density",       # número promedio de interacciones por especie
  "weighted connectance",  # conectancia ponderada por frecuencia
  "H2"                     # especialización a nivel de red (0=generalista, 1=especialista)
))

# Métricas a nivel de especie
specieslevel(df1, index = c(
  "degree",           # número de interacciones por especie
  "normalised degree", # grado normalizado
  "d"                 # especialización individual
))

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
str(azu2)

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

# Modelo 1: GLM de concentración de fructosa
mod1 = glm(fructosa ~ condicion * sexo * tratamiento, 
           data = azu2, family = Gamma(link = "log"))
summary(mod1)
library(car)
Anova(mod1)

# Modelo 2: análisis de varianza de Concentración de Sacarosa
mod2 = glm(glucosa ~ condicion * sexo * tratamiento, 
           data = azu2, family = Gamma(link = "log"))
summary(mod2)
Anova(mod2)

# Modelo 3: análisis de varianza de Concentración de Sacarosa
mod3 = glm(sacarosa ~ condicion * sexo * tratamiento,
           data = azu2, family = gaussian(link = "identity"))
summary(mod3)
Anova(mod3)

# Modelo 4: análisis de varianza de Concentración total
mod4 = glm(Suma.total ~ condicion * sexo * tratamiento,
           data = azu2, family = gaussian(link = "identity"))
summary(mod4)
Anova(mod4)


# media
media2 = azu2 %>%
  dplyr::group_by(tratamiento) %>%
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

### Gráfica de Hill numbers: q0
head(azu2)
nombres_condicion <- c("domesticada" = "Domesticated",
                       "silvestre"   = "Wild")

nombres_sexo <- c("hembra" = "Pistillate",
                  "macho"  = "Staminate")
# plot
nectar = ggplot(azu2, aes(y = sacarosa, x = tratamiento, fill = tratamiento))+
  geom_boxplot(alpha = .8)+
  #geom_jitter(position = position_jitter(0.1), size = 1, alpha = .8)+
  facet_grid(sexo ~ condicion,
             labeller = labeller(
               condicion = nombres_condicion,
               sexo = nombres_sexo
             ))+
  labs(x = "Treatment",
       y = "Sugar concentration (w/v)")+
  scale_x_discrete(label = c("t1" = "Unvisited", "t2" = "Visited"))+
  scale_fill_discrete(name = "tratamiento", 
                      labels = c("t1" = "Unvisited",
                                 "t2" = "Visited"),
                      palette = c("#E6AB02","#4D9221"))+
  theme_bw()+
  theme(legend.position = "none")

nectar

# plot 1 : fructosa
plot1 = ggplot(azu2, aes(y = fructosa, x = tratamiento, fill = tratamiento))+
  geom_boxplot(alpha = .8)+
  #geom_jitter(position = position_jitter(0.1), size = 1, alpha = .8)+
  labs(x = "Treatment",
       y = "Sugar concentration (w/v)")+
  scale_x_discrete(label = c("t1" = "Unvisited", "t2" = "Visited"))+
  scale_fill_discrete(name = "tratamiento", 
                      labels = c("t1" = "Unvisited",
                                 "t2" = "Visited"),
                      palette = c("#E6AB02","#4D9221"))+
  theme_bw()+
  theme(legend.position = "none")
plot1
# plot 2 : Glucose
plot2 = ggplot(azu2, aes(y = glucosa, x = tratamiento, fill = tratamiento))+
  geom_boxplot(alpha = .8)+
  #geom_jitter(position = position_jitter(0.1), size = 1, alpha = .8)+
  labs(x = "Treatment",
       y = "Sugar concentration (w/v)")+
  scale_x_discrete(label = c("t1" = "Unvisited", "t2" = "Visited"))+
  scale_fill_discrete(name = "tratamiento", 
                      labels = c("t1" = "Unvisited",
                                 "t2" = "Visited"),
                      palette = c("#E6AB02","#4D9221"))+
  theme_bw()+
  theme(legend.position = "none")
plot2
# plot 3 : Glucose
plot3 = ggplot(azu2, aes(y = sacarosa, x = tratamiento, fill = tratamiento))+
  geom_boxplot(alpha = .8)+
  #geom_jitter(position = position_jitter(0.1), size = 1, alpha = .8)+
  labs(x = "Treatment",
       y = "Sugar concentration (w/v)")+
  scale_x_discrete(label = c("t1" = "Unvisited", "t2" = "Visited"))+
  scale_fill_discrete(name = "tratamiento", 
                      labels = c("t1" = "Unvisited",
                                 "t2" = "Visited"),
                      palette = c("#E6AB02","#4D9221"))+
  theme_bw()+
  theme(legend.position = "none")
plot3
# plot 4 : total
plot4 = ggplot(azu2, aes(y = Suma.total, x = tratamiento, fill = tratamiento))+
  geom_boxplot(alpha = .8)+
  geom_signif(comparisons  = list(c("t1", "t2")),  # grupos a comparar
              annotations  = "*",                   # ← asterisco
              y_position   = max(azu2$Suma.total, na.rm = TRUE) * .9,
              tip_length   = 0.02,
              textsize     = 5) +
  #geom_jitter(position = position_jitter(0.1), size = 1, alpha = .8)+
  labs(x = "Treatment",
       y = "Sugar concentration (w/v)")+
  scale_x_discrete(label = c("t1" = "Unvisited", "t2" = "Visited"))+
  scale_fill_discrete(name = "tratamiento", 
                      labels = c("t1" = "Unvisited",
                                 "t2" = "Visited"),
                      palette = c("#E6AB02","#4D9221"))+
  theme_bw()+
  theme(legend.position = "none")
plot4

plot_azucar = ggarrange(plot1, plot2, plot3, plot4,
                        nrow = 2, ncol = 2, labels = "AUTO")
plot_azucar



str(azu2)
head(azu2)
sort(azu2$id_planta)
sort(azu2$sexo)

str(hill_wide_ic)
head(hill_wide_ic)
sort(hill_wide_ic$sexo)
sort(hill_wide_ic$ano)
tapply(hill_wide_ic$qD_q0, hill_wide_ic$ano, length)

### Paso 1: Verificar los IDs de cada tabla
# Ver IDs únicos en cada tabla
sort(unique(azu2$id_planta))
sort(unique(hill_wide_ic$id_planta))

# Ver cuántos coinciden
ids_azu  <- unique(tolower(azu2$id_planta))         # convertir a minúsculas
ids_hill <- unique(tolower(hill_wide_ic$id_planta))  # convertir a minúsculas

# IDs que están en azu2 pero NO en hill
setdiff(ids_azu, ids_hill)

# IDs que están en hill pero NO en azu2
setdiff(ids_hill, ids_azu)

# IDs que coinciden en ambas
intersect(ids_azu, ids_hill)

### Paso 2: Estandarizar la capitalización antes de unir
library(dplyr)

# Convertir id_planta a minúsculas en ambas tablas
azu2 <- azu2 %>%
  mutate(ano = as.character(ano),
         id_planta = tolower(id_planta))

hill_wide_ic <- hill_wide_ic %>%
  mutate(id_planta = tolower(id_planta))



### Paso 3: Elegir el tipo de join según tu objetivo
# Inner join: solo conserva filas que coinciden en ambas tablas
datos_unidos <- hill_wide_ic %>%
  inner_join(azu2 %>%
               dplyr::select(ano, id_planta, sexo, tratamiento, fructosa, glucosa,
                             sacarosa, Suma.total),
             by = c("ano","id_planta","sexo", "tratamiento"))

# Verificar cuántas filas quedaron
dim(datos_unidos)
datos_unidos
View(datos_unidos)

# datos de hongos
datos_unidos_its <- hill_wide_ic_its %>%
  inner_join(azu2 %>%
               dplyr::select(ano, sexo, tratamiento, fructosa, glucosa,
                             sacarosa, Suma.total),
             by = c("ano","sexo", "tratamiento"))

# Verificar cuántas filas quedaron
dim(datos_unidos_its)
datos_unidos_its
View(datos_unidos_its)






# left_join
datos_unidos_left <- hill_wide_ic %>%
  left_join(azu2 %>% 
              dplyr::select(año, id_planta, sexo, tratamiento, fructosa, glucosa, 
                            sacarosa, Suma.total),
            by = c("año","id_planta", "sexo", "tratamiento"))
View(datos_unidos_left)

# Ver cuántos NAs quedan por columna
colSums(is.na(datos_unidos_left))

# Ver cuántas filas completas tienes
sum(complete.cases(datos_unidos_left))

# Ver cuántas filas tienen NA en azúcares
sum(is.na(datos_unidos_left$fructosa))

# Correlaciones (na.rm automático)
fruc1 = cor.test(datos_unidos_left$qD_q1, datos_unidos_left$fructosa,
                 method = "spearman", use    = "complete.obs")

fruc1

# Ver cuántas muestras hay por grupo
table(datos_unidos$condicion)
table(datos_unidos$tratamiento)
table(datos_unidos$sexo)

# Ver cuántas muestras por año
table(datos_unidos$año)

# Ver combinaciones
table(datos_unidos$condicion,datos_unidos$sexo, datos_unidos$tratamiento)
table(datos_unidos$condicion, datos_unidos$sexo)

boxplot(datos_unidos$qD_q2 ~ datos_unidos$condicion)
boxplot(datos_unidos$qD_q2 ~ datos_unidos$sexo)
boxplot(datos_unidos$qD_q2 ~ datos_unidos$tratamiento)


boxplot(datos_unidos$fructosa ~ datos_unidos$condicion)
boxplot(datos_unidos$fructosa ~ datos_unidos$sexo)
boxplot(datos_unidos$fructosa ~ datos_unidos$tratamiento)


boxplot(datos_unidos$glucosa ~ datos_unidos$condicion)
boxplot(datos_unidos$glucosa ~ datos_unidos$sexo)
boxplot(datos_unidos$glucosa ~ datos_unidos$tratamiento)

boxplot(datos_unidos$sacarosa ~ datos_unidos$condicion)
boxplot(datos_unidos$sacarosa ~ datos_unidos$sexo)
boxplot(datos_unidos$sacarosa ~ datos_unidos$tratamiento)

boxplot(datos_unidos$Suma.total ~ datos_unidos$condicion)
boxplot(datos_unidos$Suma.total ~ datos_unidos$sexo)
boxplot(datos_unidos$Suma.total ~ datos_unidos$tratamiento)


cor_fru_q0 <- cor.test(datos_unidos$qD_q0, datos_unidos$fructosa, method = "spearman")
cor_fru_q1 <- cor.test(datos_unidos$qD_q1, datos_unidos$fructosa, method = "spearman")
cor_fru_q2 <- cor.test(datos_unidos$qD_q2, datos_unidos$fructosa, method = "spearman")

cor_glu_q0 <- cor.test(datos_unidos$qD_q0, datos_unidos$glucosa, method = "spearman")
cor_glu_q1 <- cor.test(datos_unidos$qD_q1, datos_unidos$glucosa, method = "spearman")
cor_glu_q2 <- cor.test(datos_unidos$qD_q2, datos_unidos$glucosa, method = "spearman")

cor_sac_q0 <- cor.test(datos_unidos$qD_q0, datos_unidos$sacarosa, method = "spearman")
cor_sac_q1 <- cor.test(datos_unidos$qD_q1, datos_unidos$sacarosa, method = "spearman")
cor_sac_q2 <- cor.test(datos_unidos$qD_q2, datos_unidos$sacarosa, method = "spearman")

cor_tot_q0 <- cor.test(datos_unidos$qD_q0, datos_unidos$Suma.total, method = "spearman")
cor_tot_q1 <- cor.test(datos_unidos$qD_q1, datos_unidos$Suma.total, method = "spearman")
cor_tot_q2 <- cor.test(datos_unidos$qD_q2, datos_unidos$Suma.total, method = "spearman")

# Resumir todos en una tabla
tabla_cor <- data.frame(
  Azucar     = rep(c("Fructosa", "Glucosa", "Sacarosa", "Suma.total"), each = 3),
  Hill_q     = rep(c("q0", "q1", "q2"), times = 4),
  rho        = c(cor_fru_q0$estimate, cor_fru_q1$estimate, cor_fru_q2$estimate,
                 cor_glu_q0$estimate, cor_glu_q1$estimate, cor_glu_q2$estimate,
                 cor_sac_q0$estimate, cor_sac_q1$estimate, cor_sac_q2$estimate,
                 cor_tot_q0$estimate, cor_tot_q1$estimate, cor_tot_q2$estimate),
  p_valor    = c(cor_fru_q0$p.value, cor_fru_q1$p.value, cor_fru_q2$p.value,
                 cor_glu_q0$p.value, cor_glu_q1$p.value, cor_glu_q2$p.value,
                 cor_sac_q0$p.value, cor_sac_q1$p.value, cor_sac_q2$p.value,
                 cor_tot_q0$p.value, cor_tot_q1$p.value, cor_tot_q2$p.value)
)


tabla_cor$significancia <- ifelse(tabla_cor$p_valor <= 0.001, "***",
                                  ifelse(tabla_cor$p_valor <= 0.01,  "**",
                                         ifelse(tabla_cor$p_valor <= 0.05,  "*",
                                                ifelse(tabla_cor$p_valor <= 0.10,  ".", "ns"))))

print(tabla_cor)


mod_azu1 = glm(qD_q0 ~ fructosa + glucosa + sacarosa +
                 condicion * sexo * tratamiento,
               data = datos_unidos,
               family = poisson(link = "log"))
summary(mod_azu1)
Anova(mod_azu1)

# Calcular ratio de dispersión (debe ser ~1 si no hay sobredispersión)
E2 <- resid(mod_azu1, type = "pearson")
sum(E2^2) / df.residual(mod_azu1)
# Modelo con q0 binomial negativo
library(MASS)
# q0 modelo fructosa
modnb_azu1 <- glm.nb(qD_q0 ~ fructosa +
                       tratamiento,
                     data = datos_unidos,
                   link = "log")
summary(modnb_azu1)
Anova(modnb_azu1)
# q0 con glucosa
modnb_azu1 <- glm.nb(qD_q0 ~ glucosa +
                       tratamiento,
                     data = datos_unidos,
                     link = "log")
summary(modnb_azu1)
Anova(modnb_azu1)
# q0 con sacarosa
modnb_azu1 <- glm.nb(qD_q0 ~ sacarosa +
                       tratamiento,
                     data = datos_unidos,
                     link = "log")
summary(modnb_azu1)
Anova(modnb_azu1)
# q0 con concentración total
modnb_azu1 <- glm.nb(qD_q0 ~ Suma.total +
                       tratamiento,
                     data = datos_unidos,
                     link = "log")
summary(modnb_azu1)
Anova(modnb_azu1)







### Test de Mantel
library(vegan)
library(phyloseq)

# Filtrar el objeto phyloseq a las 53 muestras
muestras_53 <- datos_unidos$sampleid2
ps_53 <- prune_samples(muestras_53, ps_con_arbol)
nsamples(ps_53)
dim(asv_53)                        # ¿cuántas filas tiene la matriz?
length(muestras_53)                # ¿cuántos IDs hay en el vector?

dim(asv_53)                        # ¿cuántas filas tiene la matriz?
length(muestras_53)                # ¿cuántos IDs hay en el vector?

# Comparar IDs
ids_datos   <- datos_unidos$sampleid2
ids_phyloseq <- sample_names(ps_con_arbol)

# ¿Cuántos coinciden?
sum(ids_datos %in% ids_phyloseq)   # debe ser 53

# Verificar dimensiones detalladas
dim(otu_table(ps_53))        # ¿cuántas filas y columnas?
taxa_are_rows(ps_53)         # TRUE = ASVs en filas, FALSE = ASVs en columnas

# Ver las primeras filas y columnas
otu_table(ps_53)[1:3, 1:3]


# Recalcular beta diversidad
asv_53    <- t(as.data.frame(otu_table(ps_53)))
dist_bray_53 <- vegdist(asv_53, method = "bray")

# Matriz de distancias de azúcares
dist_azucares <- dist(datos_unidos[, c("fructosa", "glucosa", "sacarosa",
                                       "Suma.total")],
                      method = "euclidean")

# Mantel test
set.seed(123)
mantel(dist_bray_53, dist_azucares,
       method       = "spearman",
       permutations = 999)

####
library(vegan)

# Matriz de ASVs de las 53 muestras
asv_53 <- t(as.data.frame(otu_table(ps_53)))

# Variables ambientales (azúcares)
azucares <- datos_unidos %>%
  dplyr::select(fructosa, glucosa, sacarosa)

# dbRDA con Bray-Curtis
dbrda_bac <- dbrda(asv_53 ~ fructosa + glucosa + sacarosa,
                   data   = azucares,
                   method = "bray")

# Prueba de significancia global
set.seed(123)
anova(dbrda_bac, permutations = 999)

# Prueba por cada variable
anova(dbrda_bac, by = "terms", permutations = 999)

# Varianza explicada
RsquareAdj(dbrda_bac)


### NMDS
library(vegan)

# Matriz de ASVs de las 53 muestras completas
asv_53 <- t(as.data.frame(otu_table(ps_53)))
#asvs <- t(as.data.frame(otu_table(s16_6)))
unifrac_w  <- UniFrac(ps_con_arbol, weighted = TRUE, normalized = TRUE)

dim(asv_53)

# NMDS
set.seed(123)
nmds <- metaMDS(asv_53,
                distance = "bray",
                k        = 2,       # 2 dimensiones
                trymax   = 100)     # intentos para converger

# Verificar el stress (debe ser < 0.2)
nmds$stress
names(datos_unidos)
# Variables ambientales (azúcares)
azucares <- datos_unidos %>%
  dplyr::select(fructosa, glucosa, sacarosa, Suma.total)

# envfit: correlación de azúcares con el NMDS
set.seed(123)
ef <- envfit(nmds, azucares,
             permutations = 999,
             na.rm        = TRUE)

# Ver resultados
print(ef)

# Graficar todas las flechas aunque no sean significativas
# pero indicando cuáles son significativas
# Extraer coordenadas del NMDS
# 1. Construir nmds_coords
nmds_coords <- as.data.frame(scores(nmds, display = "sites"))
nmds_coords <- cbind(nmds_coords, 
                     condicion   = datos_unidos$condicion,
                     tratamiento = datos_unidos$tratamiento,
                     sexo        = datos_unidos$sexo)

# 2. Construir ef_coords desde el envfit
ef_coords <- as.data.frame(scores(ef, display = "vectors"))
ef_coords$variable <- rownames(ef_coords)
ef_coords$r2       <- ef$vectors$r
ef_coords$pval     <- ef$vectors$pvals
ef_coords$sig      <- ifelse(ef_coords$pval < 0.05, "p < 0.05", "p > 0.05")

# 3. Verificar que los objetos tienen el formato correcto
head(nmds_coords)
head(ef_coords)

# 4. Gráfica
ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = condicion, shape = tratamiento),
             size = 3, alpha = 0.8) +
  stat_ellipse(aes(color = condicion, group = condicion),
               linewidth = 0.7) +
  geom_segment(data = ef_coords,
               aes(x = 0, y = 0,
                   xend    = NMDS1 * 0.3,
                   yend    = NMDS2 * 0.3,
                   linetype = sig),
               arrow     = arrow(length = unit(0.3, "cm")),
               color     = "black",
               linewidth = 0.7) +
  geom_text(data = ef_coords,
            aes(x     = NMDS1 * 0.35,
                y     = NMDS2 * 0.35,
                label = paste0(variable, "\n(r²=",
                               round(r2, 3), ", p=",
                               round(pval, 3), ")")),
            size = 3, fontface = "italic") +
  scale_color_manual(values  = c("domesticada" = "#C0392B",
                                 "silvestre"   = "#2980B9")) +
  scale_shape_manual(values  = c("t1" = 16, "t2" = 17),
                     labels  = c("t1" = "Unvisited",
                                 "t2" = "Visited")) +
  scale_linetype_manual(values = c("p < 0.05" = "solid",
                                   "p > 0.05" = "dashed")) +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text       = element_text(face = "bold")) +
  labs(title    = paste0("NMDS Bray-Curtis (stress = ",
                         round(nmds$stress, 3), ")"),
       color    = "Condición",
       shape    = "Tratamiento",
       linetype = "Significancia")




muestras_problematicas <- c("ENESMo24-cLVEv416s130", "ENESMo24-cLVEv416s138")

# Suma de secuencias por muestra
seq_por_muestra <- rowSums(asv_53)
seq_por_muestra[muestras_problematicas]

# Ver el mínimo y distribución
summary(seq_por_muestra)
hist(seq_por_muestra, main = "Secuencias por muestra")

# Ver las distancias de esas muestras al resto
dist_mat <- as.matrix(dist_bray_53)
dist_mat[muestras_problematicas, 1:5]

# Si tienen distancia = 1 con todas las demás muestras
# significa que no comparten ningún ASV con el resto

# Eliminar las muestras problemáticas
muestras_ok <- setdiff(rownames(asv_53), muestras_problematicas)

asv_53_limpio       <- asv_53[muestras_ok, ]
datos_completos_ok  <- datos_unidos %>%
  filter(!sampleid2 %in% muestras_problematicas)

# Verificar
dim(asv_53_limpio)        # debe ser 51 × ASVs
nrow(datos_completos_ok)  # debe ser 51

# Recalcular NMDS sin las muestras problemáticas
set.seed(123)
nmds_limpio <- metaMDS(asv_53_limpio,
                       distance = "bray",
                       k        = 2,
                       trymax   = 100)

nmds_limpio$stress  # ahora debe ser > 0 y < 0.2

# Nuevas coordenadas sin muestras problemáticas
nmds_coords_ok <- as.data.frame(scores(nmds_limpio, display = "sites"))
nmds_coords_ok <- cbind(nmds_coords_ok,
                        condicion   = datos_completos_ok$condicion,
                        tratamiento = datos_completos_ok$tratamiento,
                        sexo        = datos_completos_ok$sexo)

# Recalcular envfit
azucares_ok <- datos_completos_ok %>%
  dplyr::select(fructosa, glucosa, sacarosa, Suma.total)

set.seed(123)
ef_ok <- envfit(nmds_limpio, azucares_ok,
                permutations = 999,
                na.rm        = TRUE)
print(ef_ok)

# 1. Construir nmds_coords con muestras limpias
nmds_coords_ok <- as.data.frame(scores(nmds_limpio, display = "sites"))
nmds_coords_ok <- cbind(nmds_coords_ok,
                        condicion   = datos_completos_ok$condicion,
                        tratamiento = datos_completos_ok$tratamiento,
                        sexo        = datos_completos_ok$sexo)

# 2. Construir ef_coords desde el envfit limpio
ef_coords_ok <- as.data.frame(scores(ef_ok, display = "vectors"))
ef_coords_ok$variable <- rownames(ef_coords_ok)
ef_coords_ok$r2       <- ef_ok$vectors$r
ef_coords_ok$pval     <- ef_ok$vectors$pvals
ef_coords_ok$sig      <- ifelse(ef_coords_ok$pval < 0.05, "p < 0.05", "p > 0.05")

# 3. Gráfica
ggplot(nmds_coords_ok, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = condicion, shape = tratamiento),
             size = 3, alpha = 0.8) +
  stat_ellipse(aes(color = condicion, group = condicion),
               linewidth = 0.7) +
  geom_segment(data = ef_coords_ok,
               aes(x = 0, y = 0,
                   xend     = NMDS1 * 0.3,
                   yend     = NMDS2 * 0.3,
                   linetype = sig),
               arrow     = arrow(length = unit(0.3, "cm")),
               color     = "black",
               linewidth = 0.7) +
  geom_text(data = ef_coords_ok,
            aes(x     = NMDS1 * 0.35,
                y     = NMDS2 * 0.35,
                label = paste0(variable, "\n(r²=",
                               round(r2, 3), ", p=",
                               round(pval, 3), ")")),
            size     = 3,
            fontface = "italic") +
  scale_color_manual(values  = c("domesticada" = "#C0392B",
                                 "silvestre"   = "#2980B9"),
                     labels  = c("domesticada" = "Domesticated",
                                 "silvestre"   = "Wild")) +
  scale_shape_manual(values  = c("t1" = 16, "t2" = 17),
                     labels  = c("t1" = "Unvisited",
                                 "t2" = "Visited")) +
  scale_linetype_manual(values = c("p < 0.05" = "solid",
                                   "p > 0.05" = "dashed")) +
  theme_bw(base_size = 13) +
  theme(legend.position  = "right",
        strip.background = element_rect(fill = "grey90"),
        strip.text       = element_text(face = "bold")) +
  labs(title    = paste0("NMDS Bray-Curtis (stress = ",
                         round(nmds_limpio$stress, 3), ")"),
       x        = "NMDS1",
       y        = "NMDS2",
       color    = "Condición",
       shape    = "Tratamiento",
       linetype = "Significancia")






nmds_coords <- as.data.frame(scores(nmds, display = "sites"))
nmds_coords <- cbind(nmds_coords, azucares)

# Extraer vectores de envfit (solo los significativos)
ef_coords <- as.data.frame(scores(ef, display = "vectors"))
ef_coords$variable <- rownames(ef_coords)
ef_coords$r2       <- ef$vectors$r
ef_coords$pval     <- ef$vectors$pvals



ef_coords <- ef_coords %>%
  mutate(sig = ifelse(pval < 0.05, "p < 0.05", "p > 0.05"))

ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = condicion, shape = tratamiento),
             size = 3, alpha = 0.8) +
  stat_ellipse(aes(color = condicion, group = condicion),
               linewidth = 0.7) +
  geom_segment(data = ef_coords,
               aes(x = 0, y = 0,
                   xend = NMDS1 * 0.3,
                   yend = NMDS2 * 0.3,
                   linetype = sig),
               arrow     = arrow(length = unit(0.3, "cm")),
               color     = "black",
               linewidth = 0.7) +
  geom_text(data = ef_coords,
            aes(x     = NMDS1 * 0.35,
                y     = NMDS2 * 0.35,
                label = paste0(variable, "\n(r²=",
                               round(r2, 3), ", p=",
                               round(pval, 3), ")")),
            size = 3, fontface = "italic") +
  scale_color_manual(values   = c("domesticada" = "#C0392B",
                                  "silvestre"   = "#2980B9")) +
  scale_shape_manual(values   = c("t1" = 16, "t2" = 17)) +
  scale_linetype_manual(values = c("p < 0.05" = "solid",
                                   "p > 0.05" = "dashed")) +
  theme_bw(base_size = 13) +
  labs(title    = paste0("NMDS Bray-Curtis (stress = ",
                         round(nmds$stress, 3), ")"),
       color    = "Condición",
       shape    = "Tratamiento",
       linetype = "Significancia")

###############################################################################
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


