# Pipeline for bioinformatic analysis of ITS sequences in Cucurbita species
# Luis Alberto Villanueva Espino

# Modified from: https://benjjneb.github.io/dada2/ITS_workflow.html

# Address of LANASE cluster
# IP: 132.247.186.43
# Port: 60307

# Cluster connection
#ssh villanueva@132.247.186.43 -p 60307

# First activate conda environment "Bioconductor"
# Then, activate R in Bioconductor environment

# Loading R packages
library(dada2)
library(ShortRead)
library(Biostrings)

# Change to directory where there are files:
path <- "/home/.lvillanueva/work/analisis_2024/its"

list.files(path)

# Forward and reverse fastq filenames have format: 
# SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
head(fnFs)
head(fnRs)

#some_files <- sample(1:375, 16, replace = FALSE)
#plotQualityProfile(ITS2_files[c(some_files)])


# Primers usados
fwd1 <- "CATCGATGAAGAACGCAG"
fwd2 <- "CAACGATGAAGAACGCAG"
fwd3 <- "CACCGATGAAGAACGCAG"
fwd5 <- "CATCGATGAAGAACGTGG"
rev <- "TCCTCCGCTTATTGATATGC"

#FWD <- c("CATCGATGAAGAACGCAG", "CAACGATGAAGAACGCAG", "CACCGATGAAGAACGCAG",
#         "CATCGATGAAGAACGTGG")

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather
  #than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), 
               Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
fwd1.orients <- allOrients(fwd1);fwd1.orients
fwd2.orients <- allOrients(fwd2);fwd2.orients
fwd3.orients <- allOrients(fwd3);fwd3.orients
fwd5.orients <- allOrients(fwd5);fwd5.orients
rev.orients <- allOrients(rev); rev.orients



# The presence of ambiguous bases (Ns) in the sequencing reads makes accurate 
# mapping of short primer sequences difficult. Next we are going to “pre-filter”
# the sequences just to remove those with Ns, but perform no other filtering.
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) 
# Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
# 
filtarN = filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0,
                        multithread = TRUE)


# We are now ready to count the number of times the primers appear in the
# forward and reverse read, while considering all possible primer orientations.
# Identifying and counting the primers on one set of paired end FASTQ files is 
# sufficient, assuming all the files were created using the same library 
# preparation, so we’ll just process the first sample.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

primer1_count <- rbind(FWD.ForwardReads = sapply(fwd1.orients, primerHits,
                fn = fnFs.filtN[[1]]), FWD.ReverseReads = sapply(fwd1.orients,
                primerHits, fn = fnRs.filtN[[1]]), 
                       REV.ForwardReads = sapply(rev.orients, primerHits,
                fn = fnFs.filtN[[1]]), REV.ReverseReads = sapply(rev.orients,
                primerHits, fn = fnRs.filtN[[1]]))


primer2_count <- rbind(FWD.ForwardReads = sapply(fwd2.orients, primerHits,
                fn = fnFs.filtN[[1]]), FWD.ReverseReads = sapply(fwd2.orients,
                primerHits, fn = fnRs.filtN[[1]]), 
                       REV.ForwardReads = sapply(rev.orients, primerHits,
                fn = fnFs.filtN[[1]]), REV.ReverseReads = sapply(rev.orients,
                primerHits, fn = fnRs.filtN[[1]]))

primer3_count <- rbind(FWD.ForwardReads = sapply(fwd3.orients, primerHits,
                fn = fnFs.filtN[[1]]), FWD.ReverseReads = sapply(fwd3.orients,
                primerHits, fn = fnRs.filtN[[1]]), 
                       REV.ForwardReads = sapply(rev.orients, primerHits,
                fn = fnFs.filtN[[1]]), REV.ReverseReads = sapply(rev.orients,
                primerHits, fn = fnRs.filtN[[1]]))

primer5_count <- rbind(FWD.ForwardReads = sapply(fwd5.orients, primerHits,
                fn = fnFs.filtN[[1]]), FWD.ReverseReads = sapply(fwd5.orients,
                primerHits, fn = fnRs.filtN[[1]]), 
                       REV.ForwardReads = sapply(rev.orients, primerHits,
                fn = fnFs.filtN[[1]]), REV.ReverseReads = sapply(rev.orients,
                primerHits, fn = fnRs.filtN[[1]]))

primer1_count
primer2_count
primer3_count
primer5_count

### Removiendo primers con cutadapt
#cutadapt <- "/bin/cutadapt" # CHANGE ME to the cutadapt path on your
#cutadapt <- "/opt/Anaconda/envs/bioconductor/bin/cutadapt"
# Ruta de cutadapt en el nuevo servidor
cutadapt <- "/usr/local/sbin/cutadapt"


# machine
system2(cutadapt, args = "--version") # Run shell commands from R

# We now create output filenames for the cutadapt-ed files, and define the 
# parameters we are going to give the cutadapt command. The critical parameters
# are the primers, and they need to be in the right orientation, i.e. the FWD 
# primer should have been matching the forward-reads in its forward orientation,
# and the REV primer should have been matching the reverse-reads in its forward 
# orientation. Warning: A lot of output will be written to the screen by 
# cutadapt!
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs.filtN))
fnRs.cut <- file.path(path.cut, basename(fnRs.filtN))


# Complemento reverso del forward
FWD1.RC <- dada2:::rc(fwd1)
FWD2.RC <- dada2:::rc(fwd2)
FWD3.RC <- dada2:::rc(fwd3)
FWD5.RC <- dada2:::rc(fwd5)
# Complemento reverso del Reverse
REV.RC <- dada2:::rc(rev)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", c(fwd1,fwd2,fwd3,fwd5), "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", rev, "-A", c(FWD1.RC,FWD2.RC,FWD3.RC,FWD5.RC)) 


# Run Cutadapt
#for(i in seq_along(fnFs.filtN)) {
#  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to 
#                             # remove FWD and REV from reads
#                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
#                             fnFs.filtN[i], fnRs.filtN[i])) # input files
#}

for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}


for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("--trim-n", 
                             fnFs[i], fnRs[i]))
}
# As a sanity check, we will count the presence of primers in the first 
# cutadapt-ed sample:
#rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
#      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
#      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
#      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

rbind(FWD.ForwardReads = sapply(fwd1.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(fwd1.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(rev.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(rev.orients, primerHits, fn = fnRs.cut[[1]]))

rbind(FWD.ForwardReads = sapply(fwd2.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(fwd2.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(rev.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(rev.orients, primerHits, fn = fnRs.cut[[1]]))

rbind(FWD.ForwardReads = sapply(fwd3.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(fwd3.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(rev.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(rev.orients, primerHits, fn = fnRs.cut[[1]]))

rbind(FWD.ForwardReads = sapply(fwd5.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(fwd5.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(rev.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(rev.orients, primerHits, fn = fnRs.cut[[1]]))


# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# Visualizando la calidad de los 2 primeros archivos Forward:
plotQualityProfile(cutFs[1:2])

gra1_F <- plotQualityProfile(cutFs[1:2])

# Visualizando la calidad de los 2 primeros archivos Reverse:
plotQualityProfile(cutRs[1:2])

list.files(path.cut)

### Filtrado y corte (trim):
# Colocar archivos filtrados en un subdirectorio
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


# For this dataset, we will use standard filtering paraments: maxN=0 
# (DADA2 requires sequences contain no Ns), truncQ = 2, rm.phix = TRUE and 
# maxEE=2. The maxEE parameter sets the maximum number of “expected errors” 
# allowed in a read, which is a better filter than simply averaging quality 
# scores. Note: We enforce a minLen here, to get rid of spurious very low-length
# sequences. This was not needed in the 16S Tutorial Workflow because truncLen
# already served that purpose.

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, minLen = 20, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, minLen = 50, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

# filterAndTrim relajando los parametros de maxEE y truncQ
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = 8,
                     truncQ = 8, minLen = 50, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

# muestra cuantas secuencias pasaron el filtrado
head(out)
str(out)


# Aprendiendo las tasas de error
errF <- learnErrors(filtFs, multithread=TRUE)
# 104317402 total bases in 446788 reads from 20 samples will be used for 
# learning the error rates.
#104321907 total bases in 446861 reads from 20 samples will be used for
## analisis 07/03/2024
#103640251 total bases in 444070 reads from 28 samples will be used for
#learning the error rates.
## maxEE=8
#102289011 total bases in 437812 reads from 11 samples will be used for
#learning the error rates.

errR <- learnErrors(filtRs, multithread=TRUE)
# 102674345 total bases in 446788 reads from 20 samples will be used for 
# learning the error rates.
#102678066 total bases in 446861 reads from 20 samples will be used for
#learning the error rates.

## analisis 07/03/2024
# 102057436 total bases in 444070 reads from 28 samples will be used for
# learning the error rates.

#102065542 total bases in 437071 reads from 16 samples will be used for
#learning the error rates.
#100273541 total bases in 437071 reads from 16 samples will be used for 
#learning the error rates.

## maxEE=8
#101175643 total bases in 437812 reads from 11 samples will be used for
#learning the error rates.

# Graficando tasas de error
plotErrors(errF, nominalQ=TRUE)

# Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)


# Ensamble de reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construcción de tabla de secuencias
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# muestras: 21, secuencias: 513
# muestras> 153, secuencuas> 1129
#364 3241


# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Removiendo quimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# 364 1297
#View(seqtab.nochim)
write.csv(seqtab.nochim, file = "seqtab_nochim_minlen20.csv")
head(seqtab.nochim)


sum(seqtab.nochim)/sum(seqtab)
# 0.798244
# 0.7890555
# 0.7928839

# Inspect distribution of sequence lengths without quimeras
table(nchar(getSequences(seqtab.nochim)))
# muestras: 21, secuencias: 110
# muestras: 21, secuencias: 111
# muestras: 153, secuencias: 393

# Estadisticas del filtrado
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace 
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged",
                     "nonchim")
rownames(track) <- sample.names
head(track)

write.csv(track, file = "resumen_filtrado.csv")


### Asignando taxonomía
unite.ref <- "~/all_dynamics_unite/sh_general_release_dynamic_all_29.11.2022.fasta"
#unite_fungi <- "/home/luis/fungi_unite/sh_general_release_dynamic_s_29.11.2022.fasta"

unite.ref <- "/home/.lvillanueva/its2_2021/sh_general_release_dynamic_s_all_29.11.2022.fasta"
unite.ref <- "/home/.lvillanueva/work/analisis_2024/sh_general_release_dynamic_s_all_25.07.2023.fasta"

taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE,
                       tryRC = TRUE, verbose = TRUE)

head(taxa)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
#write.csv(taxa.print, file = "taxa_fungi_all.csv")
write.csv(taxa.print, file = "taxa_fungi_all_8.csv")




#library(patchwork) # no instalado


### Transponer tabla de ASV
#seqtab.nochim2 <- read.csv("asv_tabla.csv", header = T, row.names = 1)
#asvtable <- read.csv("seqtab_nochim_minlen50.csv", header = T, row.names = 1)
#head(asvtable)
#asv2 = t(asvtable)
#head(asv2, 1)
#write.csv(asv2, file = "asv_table.csv")

# cargar primero el archivo de metadatos
meta_its <- read.csv("its_metadata_2024.csv", header = T)
head(meta_its)
# despues cargar la tabla de ASV
seqtab.its <- read.csv("asv_table.csv", header = T, row.names = 1)
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
taxo_its <-read.csv("taxa_fungi_all_2024.csv", header = T,row.names = 1)
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
arbol_its = rtree(ntaxa(its1), rooted=TRUE, tip.label=taxa_names(its1))
plot(arbol_its)
# Uniendo el arbol con el objeto phyloseq anterior
its2 = merge_phyloseq(its1, arbol_its)
its2
# Primero resumimos el contenido del objeto phyloseq
summarize_phyloseq(its2) # reads totales 6,060,868, 1257 ASVs in 364 samples
# información del número de reads
summary(sample_sums(its2)) #range of 0 to 31014

# A useful next step is to explore feature prevalence in the dataset,
# which we will define here as the number of samples in which a taxa appears
# at least once.
# Compute prevalence of each feature, store as data.frame
prevdf_its = apply(X = otu_table(its2),
                   MARGIN = ifelse(taxa_are_rows(its2), yes = 1, no = 2),
                   FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_its = data.frame(Prevalence = prevdf_its,
                        TotalAbundance = taxa_sums(its2),
                        tax_table(its2))
# Revisando prevalencia por phylum
plyr::ddply(prevdf_its, "Phylum", function(df1)
{cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Filtering phyla whose prevalence is 1 (only present in one sample)
filterPhyla_its = c("p__Bryophyta","p__Chytridiomycota", "p__Monilophyta"
                    ,"p__Mortierellomycota")
# Filter entries with low prevalence Phylum.
its3 = subset_taxa(its2, !Phylum %in% filterPhyla_its)
its3 # Final ASV data set: 1253 ASVs in 364 samples


# Filtering samples of Morelia site
its4 = subset_samples(its3, sitio == "morelia")
summarize_phyloseq(its4)  #3,796,606 reads, 1253 ASVs de 236 samples
##### Extraemos todas aquellas secuencias del reino Fungi #####
its4 = subset_taxa(its4, Kingdom == "k__Fungi")
its4
# Hay algunas muestras que quedaron sin secuencias, ya que todas sus secuencias
# pertenecian al Reino Viridiplantae. Por lo que voy a excluir esas 
# muestras que quedaron con 0 secuencias:
its5 = prune_samples(sample_sums(its4)>1, its4)
its5

#Remove taxa with 0 observation
its6 <- prune_taxa(taxa_sums(its5) > 0, its5)
its6
summarize_phyloseq(its6) #16692 reads, 210 ASVs in 178 samples
sort(sample_sums(its6))

# información del número de reads
summary(sample_sums(its6))
range(sample_sums(its6))
#Library size spans 2 - 2873 reads 

# Obteniendo el número de read por muestra
reads_sample = readcount(its6)
reads_sample[1:5]
#write.csv(reads_sample, file = "reads_sample_its_todas.csv")


# Get count of phyla
table(phyloseq::tax_table(its6)[,"Phylum"])
# Convert to relative abundance
ps_rel_abund = phyloseq::transform_sample_counts(its6, function(x)
  {x / sum(x)})
phyloseq::otu_table(its6)[1:5, 1:5]
phyloseq::otu_table(ps_rel_abund)[1:5, 1:5]

#### grafica de abundancia con phyloseq
# creamos una primera visualización de la abundancia total de secuencias por
# muestra a nivel de Phylum
plot_bar(its6, fill = "Phylum") +
  geom_bar(aes(color = Phylum), stat = "identity",
           position = "stack")
# abundancia total de secuencias a nivel de Order
plot_bar(its6, fill = "Order") +
  geom_bar(aes(color = Order), stat = "identity",
           position = "stack")
# abundancia total de secuencias a nivel de Order
plot_bar(its6, fill = "Family") +
  geom_bar(aes(color = Family), stat = "identity",
           position = "stack")

# Abundancia relativa a nivel de Phylum
plot_bar(ps_rel_abund, fill = "Phylum") +
  geom_bar(aes(color = Phylum), stat = "identity", position = "stack")
# Abundancia relativa a nivel de Class
plot_bar(ps_rel_abund, fill = "Class") +
  geom_bar(aes(color = Class), stat = "identity", position = "stack")

# Elimando NAs a nivel de phylum
#its2_4 = subset_taxa(its2_3, Phylum != "NA")
#its2_4
#
#its2_4_re  = transform_sample_counts(its2_4, function(x) x / sum(x) )
#its2_4_re

#plot_bar(its2_4_re, fill="Phylum") +
#  geom_bar(aes(color = Phylum), stat = "identity", position = "stack")

# eliminando muestras sin secuencias
#its2_5 = prune_samples(sample_sums(its2_4)>1, its2_4)
#its2_5
#otu_table(its2_5)
#View(otu_table(its2_5))

#its2_5_re  = transform_sample_counts(its2_5, function(x) x / sum(x) )
#its2_5_re
#plot_bar(its2_5_re, fill="Phylum") +
#  geom_bar(aes(color = Phylum), stat = "identity", position = "stack")

# checking correlation between reads and ASVs observados
ggplot(data = data.frame("total_reads" =  phyloseq::sample_sums(its6),
                         "observed" = phyloseq::estimate_richness(its6,
                           measures = "Observed")[, 1]),
       aes(x = total_reads, y = observed)) +
  geom_point() +
  geom_smooth(method="lm", se = FALSE) +
  labs(x = "\nTotal Reads", y = "Observed Richness\n")

##############################################################################
ps_fam_its <- tax_glom(its6, taxrank = "Family", NArm = FALSE)

# B. Transformar a abundancia relativa
ps_rel_its <- transform_sample_counts(ps_fam_its, function(x) x / sum(x))

# C. Psmelt - cada fila es una muestra × familia
df_melt_its <- psmelt(ps_rel_its) %>%
  mutate(Family = ifelse(is.na(Family), "Unclassified", 
                         as.character(Family)))

# D. Top 10 familias más abundantes globalmente
top10_familias_its <- df_melt_its %>%
  group_by(Family) %>%
  summarise(AbundMedia = mean(Abundance)) %>%
  slice_max(order_by = AbundMedia, n = 10) %>%
  pull(Family)

# E. Colapsar el resto en "Other" por muestra
df_final_its <- df_melt_its %>%
  mutate(Taxonomy_Display = ifelse(Family %in% top10_familias_its,
                                   Family, "Other")) %>%
  group_by(Sample, especie, Taxonomy_Display) %>%
  summarise(Percentage = sum(Abundance) * 100, .groups = "drop")

# F. Ordenar factores
taxa_names_its <- sort(setdiff(unique(df_final_its$Taxonomy_Display), "Other"))
df_final_its$Taxonomy_Display <- factor(df_final_its$Taxonomy_Display,
                                    levels = c(taxa_names_its, "Other"))

# G. Ordenar muestras por especie para que queden agrupadas
df_final_its <- df_final_its %>%
  arrange(especie, Sample)
df_final_its$Sample <- factor(df_final_its$Sample, 
                          levels = unique(df_final_its$Sample))

# H. Paleta de colores
n_colors_its <- length(unique(df_final_its$Taxonomy_Display))
getPalette_its <- colorRampPalette(brewer.pal(12, "Paired"))
custom_colors_its <- getPalette_its(n_colors_its)
names(custom_colors_its) <- levels(df_final_its$Taxonomy_Display)
custom_colors_its["Other"] <- "#D3D3D3"

# I. Gráfica por muestra
grafico_muestras_its <- ggplot(df_final_its,
                           aes(x = Sample,
                               y = Percentage,
                               fill = Taxonomy_Display)) +
  geom_bar(stat = "identity", position = "stack",
           color = NA) +  # sin borde para 227 muestras
  scale_fill_manual(values = custom_colors_its) +
  labs(x    = "Sample",
       y    = "Relative Abundance (%)",
       fill = "Family") +
  theme_classic() +
  theme(axis.text.x   = element_blank(),  # ocultar nombres de muestras
        axis.ticks.x  = element_blank(),  # son 227, no caben
        axis.title    = element_text(size = 12, face = "bold"),
        legend.text   = element_text(size = 9))

print(grafico_muestras_its)
##############################################################################
# Separando por sexos florales
hembras_its <- subset_samples(its6, sexo == "hembra")
hembras_its
machos_its <- subset_samples(its6, sexo == "macho")
machos_its

###### HEMBRAS ####
# Ver si hay familias duplicadas con diferente contexto taxonómico
tax_table(hembras_its) %>%
  as.data.frame() %>%
  group_by(Family) %>%
  summarise(
    n_phyla = n_distinct(Phylum),
    n_class = n_distinct(Class),
    n_order = n_distinct(Order)
  ) %>%
  filter(n_order > 1 | n_class > 1) %>%
  print(n = Inf)

# Reemplazar NAs en TODA la tax_table, no solo en Family
tax <- tax_table(hembras_its) %>% as.data.frame()

# Reemplazar NA en cada nivel taxonómico
tax <- tax %>%
  mutate(across(everything(), ~ifelse(is.na(.), "Unclassified", .)))

# Para los que tienen Family = "Unclassified", estandarizar
# toda su jerarquía superior para que tax_glom los fusione
tax <- tax %>%
  mutate(
    Phylum = ifelse(Family == "Unclassified", "Unclassified", Phylum),
    Class  = ifelse(Family == "Unclassified", "Unclassified", Class),
    Order  = ifelse(Family == "Unclassified", "Unclassified", Order)
  )

# Devolver la tax_table al objeto phyloseq
tax_table(hembras_its) <- tax_table(as.matrix(tax))

# Verificar que ahora todos los Unclassified tienen jerarquía homogénea
tax_table(hembras_its) %>%
  as.data.frame() %>%
  group_by(Family) %>%
  summarise(
    n_phyla = n_distinct(Phylum),
    n_class = n_distinct(Class),
    n_order = n_distinct(Order)
  ) %>%
  filter(n_order > 1 | n_class > 1)
# Debería devolver tibble vacío (0 filas)

# Visualizar tabla modificada
View(tax_table(hembras_its))

# 1. Preparación del objeto Phyloseq
# A. Aglomerar a nivel de Familia. Esto suma los reads de todos los ASVs/OTUs
# que pertenecen a la misma familia. Reduce la complejidad drásticamente.
ps_fam_itsh <- tax_glom(hembras_its, taxrank = "Family", NArm = FALSE)

# B. Transformar a abundancia relativa (escala 0-1)
# Esto estandariza las muestras independientemente de su profundidad de secuenciación.
ps_rel_itsh <- transform_sample_counts(ps_fam_itsh, function(x) x / sum(x))

# C. Convertir el objeto phyloseq a un dataframe largo ('melted' dataframe)
df_melt_itsh <- psmelt(ps_rel_itsh)
View(df_melt_itsh)

# D. Agrupar por Especie y Familia, y calcular la abundancia relativa PROMEDIO
# para esa especie de planta.
#library(dplyr)
df_grouped_itsh <- df_melt_itsh %>%
  dplyr::group_by(especie, Family) %>% 
  dplyr::summarize(MeanRelAbund = mean(Abundance), .groups = "drop")
# Visualizar tabla 
View(df_grouped_itsh)
write.csv(df_grouped_itsh, file = "abundancias_its_ABRIL-2026.csv")

# Verificar que ahora suma 100%
df_grouped_itsh %>%
  group_by(especie) %>%
  summarise(total = sum(MeanRelAbund) * 100)

# E. Identificar el Top 5 de familias más abundantes POR ESPECIE.
# Usamos slice_max para quedarnos con las 5 filas con mayor abundancia promedio por grupo.
top10_per_species_itsh <- df_grouped_itsh %>%
  group_by(especie) %>%
  slice_max(order_by = MeanRelAbund, n = 5, with_ties = FALSE) %>%
  mutate(Classification = Family) %>% # Guardamos el nombre para el cruce
  ungroup()

# F. Crear la categoría "Otras".
# Cruzamos la tabla agrupada completa con la tabla del Top 10.
df_final_itsh <- df_grouped_itsh %>%
  left_join(select(top10_per_species_itsh, especie, Family, Classification), 
            by = c("especie", "Family")) %>%
  # Si la familia NO estaba en el top10 (Classification es NA), la llamamos "Otras"
  dplyr::mutate(Taxonomy_Display = ifelse(is.na(Classification), 
                                          "Otras", 
                                          as.character(Family))) %>%
  # Volvemos a agrupar y sumar para colapsar todas las "Otras" en una sola fila por especie
  dplyr::group_by(especie, Taxonomy_Display) %>%
  dplyr::summarize(TotalMeanAbund = sum(MeanRelAbund), .groups = 'drop') %>%
  # Convertir abundancia relativa (0-1) a porcentaje (0-100) para el gráfico
  dplyr::mutate(Percentage = TotalMeanAbund * 100)
View(df_final_itsh)
write.csv(df_final_itsh, file = "abundancias_family_hembrasITS_abril2026.csv")

# 3. Preparación para Graficar
# ==============================================================================

# G. Ordenar los factores para que "Otras" siempre aparezca abajo en la barra y leyenda.
# Primero obtenemos todos los nombres únicos excepto "Otras" y los ordenamos alfabéticamente
taxa_names_itsh <- sort(setdiff(unique(df_final_itsh$Taxonomy_Display), "Otras"))
# Creamos el factor con "Otras" al final
df_final_itsh$Taxonomy_Display <- factor(df_final_itsh$Taxonomy_Display, 
                                    levels = c(taxa_names_itsh, "Otras"))

# H. Definir una paleta de colores. El Top 10 por especie puede generar
# fácilmente entre 15 y 30 familias únicas en total para graficar.
# Necesitamos una paleta grande.
n_colors_itsh <- length(unique(df_final_itsh$Taxonomy_Display))
getPalette_itsh = colorRampPalette(brewer.pal(12, "Paired"))
custom_colors_itsh <- getPalette_itsh(n_colors_itsh)


# Asignar explícitamente el color gris a "Otras"
names(custom_colors_itsh) <- levels(df_final_itsh$Taxonomy_Display)
custom_colors_itsh["Otras"] <- "#D3D3D3" # Gris claro
custom_colors_itsh["Unclassified"] <- "#404040" # Gris oscuro/Casi negro para las desconocidas


# 4. Visualización con GGPLOT2
# ==============================================================================

grafico_taxonomia_its_H <- ggplot(df_final_itsh, aes(x = especie, y = Percentage,
                                               fill = Taxonomy_Display)) +
  # Barras apiladas al 100% (position = "stack")
  geom_bar(stat = "identity", position = "stack", color = "black",
           linewidth = 0.2) +
  
  # Usar la paleta de colores personalizada
  scale_fill_manual(values = custom_colors_itsh) +
  
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
print(grafico_taxonomia_its_H)


#### MACHOS #####
# Ver si hay familias duplicadas con diferente contexto taxonómico
tax_table(machos_its) %>%
  as.data.frame() %>%
  group_by(Family) %>%
  summarise(
    n_phyla = n_distinct(Phylum),
    n_class = n_distinct(Class),
    n_order = n_distinct(Order)
  ) %>%
  filter(n_order > 1 | n_class > 1) %>%
  print(n = Inf)

# Reemplazar NAs en TODA la tax_table, no solo en Family
tax <- tax_table(machos_its) %>% as.data.frame()

# Reemplazar NA en cada nivel taxonómico
tax <- tax %>%
  mutate(across(everything(), ~ifelse(is.na(.), "Unclassified", .)))

# Para los que tienen Family = "Unclassified", estandarizar
# toda su jerarquía superior para que tax_glom los fusione
tax <- tax %>%
  mutate(
    Phylum = ifelse(Family == "Unclassified", "Unclassified", Phylum),
    Class  = ifelse(Family == "Unclassified", "Unclassified", Class),
    Order  = ifelse(Family == "Unclassified", "Unclassified", Order)
  )

# Devolver la tax_table al objeto phyloseq
tax_table(machos_its) <- tax_table(as.matrix(tax))

# Verificar que ahora todos los Unclassified tienen jerarquía homogénea
tax_table(machos_its) %>%
  as.data.frame() %>%
  group_by(Family) %>%
  summarise(
    n_phyla = n_distinct(Phylum),
    n_class = n_distinct(Class),
    n_order = n_distinct(Order)
  ) %>%
  filter(n_order > 1 | n_class > 1) # Debería devolver tibble vacío (0 filas)

# Visualizar tabla modificada
View(tax_table(machos_its))

# 1. Preparación del objeto Phyloseq
# A. Aglomerar a nivel de Familia. Esto suma los reads de todos los ASVs/OTUs
# que pertenecen a la misma familia. Reduce la complejidad drásticamente.
ps_fam_itsm <- tax_glom(machos_its, taxrank = "Family", NArm = FALSE)

# B. Transformar a abundancia relativa (escala 0-1)
# Esto estandariza las muestras independientemente de su profundidad de secuenciación.
ps_rel_itsm <- transform_sample_counts(ps_fam_itsm, function(x) x / sum(x))

# C. Convertir el objeto phyloseq a un dataframe largo ('melted' dataframe)
df_melt_itsm <- psmelt(ps_rel_itsm)

# D. Agrupar por Especie y Familia, y calcular la abundancia relativa PROMEDIO
# para esa especie de planta.
df_grouped_itsm <- df_melt_itsm %>%
  dplyr::group_by(especie, Family) %>% 
  dplyr::summarize(MeanRelAbund = mean(Abundance), .groups = "drop")
View(df_grouped_itsm)

# Verificar que ahora suma 100%
df_grouped_itsm %>%
  group_by(especie) %>%
  summarise(total = sum(MeanRelAbund) * 100)

# E. Identificar el Top 5 de familias más abundantes POR ESPECIE.
# Usamos slice_max para quedarnos con las 5 filas con mayor abundancia promedio por grupo.
top10_per_species_itsm <- df_grouped_itsm %>%
  group_by(especie) %>%
  slice_max(order_by = MeanRelAbund, n = 5, with_ties = FALSE) %>%
  mutate(Classification = Family) %>% # Guardamos el nombre para el cruce
  ungroup()

# F. Crear la categoría "Otras".
# Cruzamos la tabla agrupada completa con la tabla del Top 10.
df_final_itsm <- df_grouped_itsm %>%
  left_join(select(top10_per_species_itsm, especie, Family, Classification), 
            by = c("especie", "Family")) %>%
  # Si la familia NO estaba en el top10 (Classification es NA), la llamamos "Otras"
  dplyr::mutate(Taxonomy_Display = ifelse(is.na(Classification), 
                                          "Otras", 
                                          as.character(Family))) %>%
  # Volvemos a agrupar y sumar para colapsar todas las "Otras" en una sola fila por especie
  dplyr::group_by(especie, Taxonomy_Display) %>%
  dplyr::summarize(TotalMeanAbund = sum(MeanRelAbund), .groups = 'drop') %>%
  # Convertir abundancia relativa (0-1) a porcentaje (0-100) para el gráfico
  dplyr::mutate(Percentage = TotalMeanAbund * 100)
View(df_final_itsm)
write.csv(df_final_itsm, file = "abundancias_family_machos_abril2026.csv")
# 3. Preparación para Graficar
# ==============================================================================

# G. Ordenar los factores para que "Otras" siempre aparezca abajo en la barra y leyenda.
# Primero obtenemos todos los nombres únicos excepto "Otras" y los ordenamos alfabéticamente
taxa_names_itsm <- sort(setdiff(unique(df_final_itsm$Taxonomy_Display), "Otras"))
# Creamos el factor con "Otras" al final
df_final_itsm$Taxonomy_Display <- factor(df_final_itsm$Taxonomy_Display, 
                                      levels = c(taxa_names_itsm, "Otras"))

# H. Definir una paleta de colores. El Top 10 por especie puede generar
# fácilmente entre 15 y 30 familias únicas en total para graficar.
# Necesitamos una paleta grande.
n_colors_itsm <- length(unique(df_final_itsm$Taxonomy_Display))
getPalette_itsm = colorRampPalette(brewer.pal(12, "Paired"))
custom_colors_itsm <- getPalette_itsm(n_colors_itsm)

# Asignar explícitamente el color gris a "Otras"
names(custom_colors_itsm) <- levels(df_final_itsm$Taxonomy_Display)
custom_colors_itsm["Otras"] <- "#D3D3D3" # Gris claro

# 4. Visualización con GGPLOT2
# ==============================================================================

grafico_taxonomia_its_M <- ggplot(df_final_itsm, aes(x = especie, y = Percentage,
                                                 fill = Taxonomy_Display)) +
  # Barras apiladas al 100% (position = "stack")
  geom_bar(stat = "identity", position = "stack", color = "black",
           linewidth = 0.2) +
  
  # Usar la paleta de colores personalizada
  scale_fill_manual(values = custom_colors_itsm) +
  
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
print(grafico_taxonomia_its_M)

####################### ANALYSIS WITH PROPORTION DATA #########################
# ANALIZANDO UNICAMENTE LAS MUESTRAS DE FLORES VIRGENES (SIN VISITA) PARA
# EVALUAR EL EFECTO DE LA DOMESTICACIÓN SOBRE LA COMPOSICIÓN DE LA 
# COMUNIDAD DE MICROORGANISMOS.
#its7 = subset_samples(its6, tratamiento=="t1")
#its7
#summarize_phyloseq(its7)  #79771 reads, 1635 ASVs


### Calculating relative abundances
#mergedGP_its <- merge_samples(its7, "especie")
# Use psmelt to obtain a long-format data.frame
#phy_its <- mergedGP_its %>% tax_glom(taxrank = "Family") %>%
#  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
#View(phy_its)
#write.csv(phy_its, file = "abundancias_family_its_NV.csv")

#mergedtop = merge_samples(topfam, "especies")
#phy2 <- mergedtop %>% tax_glom(taxrank = "Family") %>%
#  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
#View(phy2)

#### HEATMAPs of most abundant bacteria ####
### Collapse dataset taxonomy to 10 most abundant families
library(microbiomeutilities)
topfam_its <- aggregate_top_taxa2(its6, "Family", top = 10)
View(otu_table(topfam_its))
topfam_its

# transform families dataset to relative abundance
transform <- microbiome::transform
famrel_its <- transform(topfam_its, "compositional")
View(otu_table(famrel_its))
#View(otu_table(famrel))

### Heatmap for 10 most abundant bacterial families between species
heatmap_its <- simple_heatmap(famrel_its,
                              group.facet = "especie",
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
heatmap_its


### Calculating relative abundances for Female flowers
merged_hembras_its <- merge_samples(hembras_its, "especie")
# Use psmelt to obtain a long-format data.frame
phy_hembras_its <- merged_hembras_its %>% tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
#View(phy_hembras_its)
write.csv(phy_hembras_its, file = "abundancias_hembras_its_NV.csv")

### Calculating relative abundances for Male flowers
merged_machos_its <- merge_samples(machos_its, "especie")
# Use psmelt to obtain a long-format data.frame
phy_machos_its <- merged_machos_its %>% tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()



#View(phy_machos_its)
write.csv(phy_machos_its, file = "abundancias_machos_its_NV.csv")


##### Flores HEMBRA
topfam_hemits <- aggregate_top_taxa2(hembras_its, "Family", top = 10)
# transform families dataset to relative abundance
transform <- microbiome::transform
famrel_hemits <- transform(topfam_hemits, "compositional")
# Heatmap de flores hembra
heatmap_hemits <- simple_heatmap(famrel_hemits,
                                 group.facet = "especie",
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
heatmap_hemits

# Flores machos
topfam_macits <- aggregate_top_taxa2(machos_its, "Family", top = 10)
# transform families dataset to relative abundance
transform <- microbiome::transform
famrel_macits <- transform(topfam_macits, "compositional")
# Heatmap de flores hembra
heatmap_macits <- simple_heatmap(famrel_macits,
                                 group.facet = "especie",
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
heatmap_macits

### Arrange heatmaps in one graphic
library(ggpubr)
library(gridExtra)
heatmaps.its <- grid.arrange(heatmap_macits, heatmap_hemits, nrow = 2)

heatmaps_its = ggarrange(heatmap_hemits, heatmap_macits, nrow = 2,
                         labels=c("A", "B"))
heatmaps_its

### Save heatmaps

ggsave("heatmaps_abril2026.pdf", 
       plot = heatmaps,
       width = 20, 
       height = 20, 
       units = "cm")

# Most dominant families en todas las muestras
domiits <- dominant_taxa(its6, level = "Family") #group = "especies"
View(domiits$dominant_overview)
write.csv(domiits$dominant_overview, file = "dominant_overview_its.csv")


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


uppenits <- get_upset(its6, factorNames = "especie")
View(uppenits)
write.csv(uppenits, file = "matrix_its.csv")
# función para gráficar los ASVs compartidos entre especies
?upset

upset(uppenits, nsets = 6, mainbar.y.label = "Number of ASVs shared",
      order.by = "freq",
      queries = list(list(query = intersects,
                          params = list("6CM", "5CAA", "4CPP"),
                          color = "orange",
                          active = T)))

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
pseq_its <- transform(its6, "compositional")

library(knitr)

# Obteniendo los valores de las metricas de diversidad alfa
alpha_its = microbiome::alpha(pseq_its, index = "all")
tab <- microbiome::alpha(pseq_its, index = c("observed","shannon", "simpson"));tab
kable(head(tab))
View(tab)

# get the metadata out as separate object
its.meta <- meta(its6)
kable(head(its.meta))
# Add the rownames as a new colum for easy integration later.
its.meta$sam_name2 <- rownames(its.meta)

# Add the rownames to diversity table
tab$sam_name2 <- rownames(tab)

# merge these two data frames into one
div.df3 <- merge(tab,its.meta, by = "sam_name2")

# check the tables
colnames(div.df3)

### Phylogenetic diversity (Faith)
library(picante)

ps0.rar.asvtab3 <- as.data.frame(its6@otu_table)

ps0.rar.tree3 <- its6@phy_tree

# hmp.meta from previous code chunks

# We first need to check if the tree is rooted or not 
its6@phy_tree

# it is a rooted tree
df.pd3 <- pd(t(ps0.rar.asvtab3), ps0.rar.tree3,include.root=T) 
# t(ou_table) transposes the table for use in picante and the tre file comes
# from the first code chunck we used to read tree file (see making a phyloseq 
# object section).
#datatable(df.pd2)

# We add the results of PD to this file and then plot
div.df3$Phylogenetic_Diversity <- df.pd3$PD
head(div.df3)

#### Diversidad alpha ####
head(div.df3)
str(div.df3)

# data
sort(div.df3$observed)
div.df3$diversity_shannon
div.df3$Phylogenetic_Diversity
# Histograma
hist(div.df3$observed)
hist(div.df3$diversity_shannon) 
hist(div.df3$Phylogenetic_Diversity)
# Shapiro test
shapiro.test(div.df3$observed)
shapiro.test(div.df3$diversity_shannon)
shapiro.test(div.df3$Phylogenetic_Diversity)

### STATISTICAL ANALYSIS ###
modits <- glm(diversity_shannon ~ condicion * sexo * tratamiento, 
               family = gaussian(link = "identity"), data = div.df3)
summary(modits)

library(car)
Anova(modits)

# POST HOC
library(emmeans)
# Comparar tratamiento dentro de cada combinación condicion:sexo
emmeans(modits, pairwise ~ tratamiento | condicion * sexo)

# Esto te permitirá ver en qué combinaciones específicas el tratamiento tiene
# efecto sobre la diversidad de hongos.
ggplot(div.df3, aes(x = tratamiento, y = diversity_shannon, 
                      color = sexo, group = sexo)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun = mean, geom = "point", size = 3) +
  facet_wrap(~ condicion) +
  theme_bw() +
  labs(title = "Interacción triple: condición × sexo × tratamiento",
       y = "Índice de Shannon")

################################################################################
################# GRÁFICAS DE INDICE DE SHANNON EN HONGOS ######################
##### VISUALIZATION  OF ALPHA DIVERSITY BY CONDITION####
sha_con_its = ggplot(div.df3, aes(x=condicion, y=diversity_shannon,
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

sha_con_its

##### VISUALIZATION  OF ALPHA DIVERSITY BY SEX ####
sha_sex_its = ggplot(div.df3, aes(x=sexo, y=diversity_shannon,
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

sha_sex_its

##### VISUALIZATION  OF ALPHA DIVERSITY BY TREATMENT #####
library(ggsignif)

sha_tra_its <- ggplot(div.df3, aes(x = tratamiento, y = diversity_shannon,
                               fill = tratamiento)) +
  geom_boxplot(alpha = 0.8) +
  geom_signif(
    comparisons = list(c("t1", "t2")),  # grupos a comparar
    annotation = "*",                   # texto: *, **, ***, ns, o p-value
    tip_length = 0.001,                   # longitud de las patitas
    y_position = max(div.df3$diversity_shannon) * .96  # altura de la barra
  ) +
  scale_fill_manual(name = "Treatment",
                    values = c("t1" = "#E6AB02", "t2" = "#4D9221"),
                    labels = c("t1" = "Unvisited", "t2" = "Visited")) +
  scale_x_discrete(labels = c("t1" = "Unvisited", "t2" = "Visited")) +
  labs(x = "Treatment", y = "Shannon diversity index") +
  theme_classic()

sha_tra_its

# combinando plots
# cargando paqueteria gridExtra
library(gridExtra)
library(cowplot)

shato3 <- ggarrange(sha_con_its, sha_sex_its, sha_tra_its, 
                    ncol = 1, nrow = 3,
                    labels = "AUTO")
shato3

######## BOXPLOT CONDITION X SEX  ################

sha_conxsex_its <- ggplot(div.df3, aes(x = condicion, y = diversity_shannon,
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
mod_spp_its = glm(diversity_shannon ~ especie,
              data = div.df3, family = gaussian)
Anova(mod_spp_its)
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
###############################################################################
###############################################################################

#### Diversidad alpha en flores macho ####
head(div.df3)
machos_alfa_its = subset(div.df3, sexo == "macho")
str(machos_alfa_its)
head(machos_alfa_its)
# data
sort(machos_alfa_its$observed)
machos_alfa_its$diversity_shannon
machos_alfa_its$Phylogenetic_Diversity
# Histograma
hist(machos_alfa_its$observed)
hist(machos_alfa_its$diversity_shannon) 
hist(machos_alfa_its$Phylogenetic_Diversity)
# Shapiro test
shapiro.test(machos_alfa_its$observed)
shapiro.test(machos_alfa_its$diversity_shannon)
shapiro.test(machos_alfa_its$Phylogenetic_Diversity)

#### STATISTICAL ANALYSIS OF ALPHA DIVERSITY ####
# Flores Hembra
modglm5 <- glm(diversity_shannon ~ condicion / especie * tratamiento, 
               family = gaussian(link = "identity"), data = hembras_alfa_its)
summary(modglm5)
anova(modglm5)

# Flores Macho
modglm6 <- glm(diversity_shannon ~ condicion / especie * tratamiento, 
               family = gaussian(link = "identity"), data = machos_alfa_its)
summary(modglm6)
anova(modglm6)

library(emmeans)

# Asumiendo que tu modelo se llama modelo_shannon
# La barra "|" le dice a R: "Compara las Especies, pero hazlo POR SEPARADO 
# dentro de cada Estado"
# comparaciones con el modelo anidado
comp5 = emmeans(modglm5,
                pairwise ~ especie | condicion)
comp5
# comparando todas las especies sin anidar
comp6 <- emmeans(modglm5, pairwise ~ especie)
comp6

# hembras
comp7 <- emmeans(modglm6,
                 pairwise ~ especie | condicion)
comp7
# comparando todas las especies sin anidar
comp8 <- emmeans(modglm6, pairwise ~ especie)
comp8


##### VISUALIZATION  OF ALPHA DIVERSITY WITH PROPORTION DATA####
# Stripchart of observed ASVs of female flowers
obs_hem_its <- ggplot(hembras_alfa_its, aes(x=especie, y=observed,
                                            color=condicion)) +
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
obs_hem_its
#Exportamos grafico de alfa div (observed) a pdf
ggsave("observed_hembras_its.pdf", 
       plot = obs_hem_its,
       width = 25, 
       height = 20, 
       units = "cm")

# Stripchart of observed ASVs of Male flowers
obs_machos_its <- ggplot(machos_alfa_its, aes(x=especie, y=observed,
                                              color=especie)) +
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
ggsave("observed_machos_its.pdf", 
       plot = obs_machos_its,
       width = 25, 
       height = 20, 
       units = "cm")

# combinando plots
# cargando paqueteria gridExtra
library(gridExtra)
library(cowplot)
obs_pro_its <- grid.arrange(obs_hem_its, obs_machos_its, 
                        ncol = 2, nrow = 1)

obs_pro_its <- ggarrange(obs_hem_its, obs_machos_its, ncol = 2, nrow = 1,
                     labels = "AUTO") 
# exporting plots
ggsave("alpha_obser_its_ene26.pdf", 
       plot = obs_pro_its,
       width = 25, 
       height = 20, 
       units = "cm")


### stripchart de Shannon diversity index in female flowers
sha_hem_its = ggplot(hembras_alfa_its, aes(x=especie, y=diversity_shannon,
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
sha_hem_its    


### stripchart de Shannon diversity index in male flowers
sha_mac_its = ggplot(machos_alfa_its, aes(x=especie, y=diversity_shannon,
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
  theme_classic()+
  theme(legend.position = "none")
sha_mac_its

#facet_wrap(~sexo, labeller = labeller(sexo = sexos))
#stat_summary(fun=mean, geom="point", shape=18,
#            size=6, color="red")
sha_pro_its <- ggarrange(sha_hem_its, sha_mac_its, ncol = 2, nrow = 1,
                     labels = "AUTO") 
# exporting plots
ggsave("alpha_shan_its_nov06.pdf", 
       plot = sha_pro_its,
       width = 25, 
       height = 20, 
       units = "cm")


shannon <- ggarrange(sha_mac, sha_hem, sha_mac_its, sha_hem_its,
                     ncol = 2, nrow = 2, labels = "AUTO")
shannon
### stripchart de Phylogenetic Diversity index in female flowers
pd_hem_its = ggplot(hembras_alfa_its, aes(x=especie, y=Phylogenetic_Diversity,
                                  color=especie)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 2)+
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
pd_hem_its    


### stripchart de Shannon diversity index in male flowers
pd_mac_its = ggplot(machos_alfa_its, aes(x=especie, y=Phylogenetic_Diversity,
                                 color=especie)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size = 2)+
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
pd_mac_its

#facet_wrap(~sexo, labeller = labeller(sexo = sexos))
#stat_summary(fun=mean, geom="point", shape=18,
#            size=6, color="red")
pd_pro_its <- ggarrange(pd_hem_its, pd_mac_its, ncol = 2, nrow = 1,
                    labels = "AUTO") 
# exporting plots
ggsave("PD_shan_its_nov06.pdf", 
       plot = pd_pro_its,
       width = 25, 
       height = 20, 
       units = "cm")

###############################################################################
############################## BETA DIVERSITY #################################
library(vegan)
# phyloseq object filtered but no rarefied
its6
# divided by sex
hembras_its
machos_its

#### Unweighted UniFrac HEMBRAS ####
# Primero nos aseguramos que estamos usando las funciones del paquete correcto
transform <- microbiome::transform
# Proportion Normalization
its_beta <- transform(its6, "compositional")


# PCoA plot using the unweighted UniFrac as distance
ordination_its = ordinate(its_beta, method="PCoA", "unifrac",
                              weighted = F)

baseplot_its <- plot_ordination(its_beta, ordination_its,
                                    color="condicion",
                                    type = "samples",
) + #shape = "sexo"
  ggtitle("")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 2)
#stat_ellipse()

baseplot_its

#baseplot_16s_hem +
#  geom_point(size = 1)
#  geom_label(aes(label = clave), size = 3.5, color = "black")


#### Unweighted UniFrac MACHOS ####
# Primero nos aseguramos que estamos usando las funciones del paquete correcto
transform <- microbiome::transform
# Proportion Normalization
its_mac <- transform(machos_its, "compositional")

# PCoA plot using the unweighted UniFrac as distance
ordination_its_mac = ordinate(its_mac, method="PCoA", "unifrac",
                              weighted = F)

baseplot_its_mac <- plot_ordination(its_mac, ordination_its_mac,
                                    color="especie",
                                    type = "samples",
) + #shape = "sexo"
  ggtitle("")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 2)
#stat_ellipse()
baseplot_its_mac


#### Weighted UniFrac####
# PCoA plot using the weighted UniFrac as distance
ordination_its_wei = ordinate(its_beta,
                                  method="PCoA", "wunifrac")
# Plot
baseplot_its_wei <- plot_ordination(its_beta, 
                                        ordination_its_wei,
                                        color="tratamiento",
                                        #shape = "condicion",
                                        type = "samples") + #shape = "sexo" shape = "tratamiento"
  ggtitle("")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 4)+
  #scale_color_discrete(name = "Species", 
   #                    labels = c("2CPF" = "CPF","3CAS" = "CAS",
    #                              "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"),
     #                  palette = "Set1")+
  scale_color_discrete(name = "Treatment",
                       labels = c("t1" = "Unvisited",
                                  "t2" = "Visited"),
                       palette = c("#E6AB02","#4D9221"))
#labs(x = "PCoa1 (31.3%)", y = "PCoA2 (10.4%)")+
#scale_color_discrete(name = "Condition", 
#                    labels = c("Domesticated","Wild"))+
#scale_color_discrete(name = "Treatment",
#                     labels = c("No-visited","Visited"))

baseplot_its_wei + stat_ellipse(linetype = 2)

#baseplot_its_hem_wei + stat_ellipse()
#
#baseplot2 +
#  geom_point(size = 1) +
#  geom_label(aes(label = clave), size = 3.5, color = "black")
#

#### Weighted UniFrac MACHOS####
# PCoA plot using the weighted UniFrac as distance
ordination_its_mac_wei = ordinate(its_mac,
                                  method="PCoA", "wunifrac")
# Plot
baseplot_its_mac_wei <- plot_ordination(its_mac, 
                                        ordination_its_mac_wei,
                                        color="especie",
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

baseplot_its_mac_wei

#baseplot_its_mac_wei + stat_ellipse()
#
#baseplot2 +
#  geom_point(size = 1) +
#  geom_label(aes(label = clave), size = 3.5, color = "black")

baseplot_16s_hem_wei
baseplot_16s_mac_wei
baseplot_its_hem_wei
baseplot_its_mac_wei

beta_analisis <- ggarrange(baseplot_16s_mac_wei, baseplot_16s_hem_wei,
                           baseplot_its_mac_wei, baseplot_its_hem_wei,
                           ncol = 2, nrow = 2, labels = "AUTO") 
beta_analisis

#### PERMANOVA by floral gender####
## female flowers
metadf_hem_its = data.frame(sample_data(hembras_its))

# calculating Weighted-UniFrac
wunifrac.dist_its_hem <- UniFrac(its_hem, 
                             weighted = TRUE, 
                             normalized = TRUE,  
                             parallel = FALSE, 
                             fast = TRUE)

# Run PERMANOVA
set.seed(100)   # Set the seed to make the analysis reproducible
perm_hem_its = adonis2(wunifrac.dist_its_hem ~ condicion / especie * tratamiento,
                       data = metadf_hem_its,
                   by = "terms", permutations = 999)
perm_hem_its


#Test dispersion among body sites
disp_wun_its <- betadisper(wunifrac.dist_its_hem,
                           sample_data(hembras_its)$tratamiento)
permutest(disp_wun_its)



## Male flowers
metadf_its = data.frame(sample_data(its_beta))

# calculating Weighted-UniFrac
wunifrac.dist_its <- UniFrac(its_beta, 
                             weighted = TRUE, 
                             normalized = TRUE,  
                             parallel = FALSE, 
                             fast = TRUE)

# Run PERMANOVA
set.seed(100)   # Set the seed to make the analysis reproducible
perm_its = adonis2(wunifrac.dist_its ~ condicion * sexo * tratamiento,
                       data = metadf_its,
                   by = "terms", permutations = 999)
perm_its


#Test dispersion among body sites
disp_wun_its <- betadisper(wunifrac.dist_its,
                               sample_data(metadf_its)$tratamiento)
permutest(disp_wun_its)

disp_wun_mac_its2 <- betadisper(wunifrac.dist_its_mac,
                                sample_data(machos_its)$especie)
permutest(disp_wun_mac_its2)








# PCoA plot using the unweighted UniFrac as distance
wunifrac_dist = phyloseq::distance(ps.rarefied, method="unifrac", weighted=F)
ordination = ordinate(ps.rarefied, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied, ordination, color="especie") + theme(aspect.ratio=1)


# PERMANOVA
adonis(wunifrac_dist ~ sample_data(ps.rarefied)$especie)
########

# VISUALIZATION
p.shannon <- boxplot_alpha(ps.rarefied, 
                           index = "shannon",
                           x_var = "especie",
                           fill.colors = c("4CAS"="cyan4",
                                           "5CAA"="deeppink4",
                                           "6CM"="red",
                                           "3CPP"="green",
                                           "2CPF" = "blue"))


p.shannon <- p.shannon + theme_minimal() + 
  labs(x="Especies", y="Shannon diversity") +
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
p1 <- ggboxplot(ps1.meta, x = "especie", y = "Shannon",
                add = c("boxplot", "jitter"), fill = "especie")+
  labs(x = "Especies", y = "Índice de diversidad de Shannon")+
  scale_fill_discrete(name = "Especies", 
                      labels = c("CPP","CAS", "CAA", "CM"))+ 
  scale_x_discrete(labels=c("1CAS" = "CAS", "2CAA" = "CAA",
                            "3CM" = "CM", "4CPP" = "CPP"))
print(p1)

p1 <- p1 + stat_compare_means(comparisons = spp.pairs) 
print(p1) # Wilcox test

# Gráfica del número de taxones observados por especie
p2 <- ggboxplot(ps1.meta, x = "especie", y = "observed",alpha =.5,
                add = c("boxplot","jitter"), fill = "especie")+
  labs(x = "Especies", y = "Riqueza observada")+
  scale_fill_discrete(name = "Especies", 
                      labels = c("CAS","CAA", "CM", "CPP"))+
  scale_x_discrete(labels=c("1CAS" = "CAS", "2CAA" = "CAA",
                            "3CM" = "CM", "4CPP" = "CPP"))
print(p2)

p2 <- p2 + stat_compare_means(comparisons = spp.pairs) 
print(p2) # Wilcox test


### gráfica del indice de Shannon por sexo
ggviolin(ps1.meta2, x = "sexo", y = "Shannon",
         add = "boxplot", fill = "sexo")+
  stat_compare_means()

### gráfica del nuḿero de taxones por tratamiento
ggviolin(ps1.meta2, x = "tratamiento", y = "Shannon",
         add = "boxplot", fill = "tratamiento")+
  stat_compare_means()

### gráfica del indice de Shannon entre especies
pmult <- ggboxplot(ps1.meta2, x = "especie", y = "Shannon", alpha = .5,
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
pmult <- ggboxplot(ps1.meta2, x = "especie", y = "Shannon", alpha = .5,
                   fill = "especie", add = "jitter", facet.by = "tratamiento",
                   panel.labs = list(tratamiento = c("No visitadas",
                                                     "Visitadas")),
                   short.panel.labs = FALSE)+ 
  labs(x = "Especies", y = "Índice de diversidad de Shannon")+
  scale_fill_discrete(name = "Especies") +
  scale_fill_discrete(name = "Especies", 
                      labels = c("CPP","CAS", "CAA", "CM"))+
  scale_x_discrete(labels=c("2CAA" = "CAA", "4CPP" = "CPP",
                            "1CAS" = "CAS", "3CM" = "CM" ))

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

### stripchart
ggplot(ps1.meta, aes(x=especie, y=Shannon, color=especie,
                     shape=tratamiento)) + 
  geom_jitter(position=position_jitter(0.2), size = 3)+
  scale_shape_manual(values=c(1,19))+
  labs(x="Especies", y = "Índice de diversidad de Shannon")+
  scale_color_discrete(name = "Especies", 
                       labels = c("CAS","CAA", "CM", "CPP"))+
  scale_x_discrete(labels=c("1CAS" = "CAS", "2CAA" = "CAA", "3CM" = "CM",
                            "4CPP" = "CPP"))+
  scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
                                                         "Visitadas")) +
  theme_bw()+
  stat_summary(fun=mean, geom="point", shape=18,
               size=6, color="red")



ggplot(iris, aes(Sepal.Length, Sepal.Width)) +
  geom_point(aes(shape = Species, color = Species), size = 3) +
  scale_shape_manual(values = c(5, 16, 17)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  theme_minimal() +
  theme(legend.position = "top")

###############################################################################
###############################################################################
########### DIVERSIDAD BETA #####
ps.rarefied
its2_5_re <- microbiome::transform(ps.rarefied, "compositional")


plot_bar(its2_5_re, fill = "Family") +
  geom_bar(aes(color = Family), stat = "identity", position = "stack")

ordu3.wt.uni1 <- ordinate(its2_5_re , "NMDS", "bray")
ordu3.wt.uni2 <- ordinate(its2_5_re , "NMDS", "jaccard")

library(vegan)
BC.nmds <- metaMDS(otu_table(ps.rarefied), distance="bray", k=2, trymax=20)

wUF.ordu = ordinate(ps.rarefied, method="NMDS", distance="unifrac",
                    weighted=TRUE)
uwUF.ordu = ordinate(ps.rarefied, method="NMDS", distance="unifrac",
                     weighted=FALSE)

#ordu3.wt.uni3 <- ordinate(its2_5_re , "PCoA", "unifrac", weighted=T)  
#ordu3.wt.uni4 <- ordinate(its2_5_re , "PCoA", "unifrac", weighted=F)


# PCoA por especie
plot_ordination(ps.rarefied, wUF.ordu, color="especie") + 
  ggtitle("Weighted Unifrac") +
  geom_point(size = 3) +
  theme_classic() + 
  scale_color_brewer("especie", palette = "Set2")+
  scale_color_discrete(name = "Especies",
                       labels = c("2CPF","3CPP", "4CAS", "5CAA", "6CM"))
# print(wt.unifrac6 + stat_ellipse(linetype = 1, cex = 1.5)) 

# PCoA por tratamiento
plot_ordination(ps.rarefied, wUF.ordu, color="tratamiento") +
  ggtitle("Weighted UniFrac")+
  geom_point(size = 3) +
  theme_classic() + 
  scale_color_brewer("tratamiento", palette = "Set2")+
  scale_color_discrete(name = "Tratamiento",
                       labels = c("No visitadas", "Visitadas"))
# print(wt.unifrac6 + stat_ellipse(linetype = 1, cex = 1.5))

# PCoA por especie y tratamiento
plot_ordination(its2_5, wUF.ordu, color="especie",
                shape = "tratamiento") + 
  ggtitle("Weighted UniFrac")+
  geom_point(size = 3) +
  theme_classic() + 
  scale_color_brewer("especie", palette = "Set2")+
  scale_color_discrete(name = "Especies",
                       labels = c("CAS","CAA", "CM", "CPP"))+
  scale_shape_discrete(name = "Tratamientos",
                       labels = c("No visitadas", "Visitadas"))
# print(wt.unifrac6 + stat_ellipse(linetype = 1, cex = 1.5))

###############################################################################
###############################################################################
############################### PERMANOVA #####################################
library(microbiome)
library(ggplot2)
library(dplyr)

# Make sure we use functions from correct package
transform <- microbiome::transform
# Pick relative abundances (compositional) and sample metadata 
pseq.rel <- microbiome::transform(ps.rarefied, "compositional")
otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)

plot_landscape(pseq.rel, method = "NMDS", distance = "bray",
               col = "especie", size = 3)

plot_landscape(pseq.rel, method = "PCoA", distance = "bray",
               col = "sexo", size = 3)
plot_landscape(pseq.rel, method = "PCoA", distance = "bray",
               col = "tratamiento", size = 3)

# haciendo la prueba de permanova
library(vegan)
set.seed(1)
permanova1 <- adonis2(t(otu) ~ especie*sexo*tratamiento,
                      data = meta, permutations=99, method = "bray")
permanova1

#### PERMANOVAS con las 4 metricas ####
### PERMANOVA con bray-curtis ###
BC.dist=vegdist(t(otu), distance="bray")
adonis2(BC.dist ~ condicion*especie*sexo*tratamiento, data = meta(ps.rarefied),
        permutations = 1000)
# PERMANOVA con jaccard
J.dist=vegdist(t(otu), distance="jaccard")
adonis2(J.dist ~ condicion*especie*sexo*tratamiento, data = meta(ps.rarefied),
        permutations = 1000)
# PERMANOVA con weighted unifrac
wUF.dist = UniFrac(ps.rarefied, weighted=TRUE, normalized=TRUE)
adonis2(wUF.dist ~ condicion*especie*sexo*tratamiento, data=meta(ps.rarefied),
        permutations = 1000)
# PERMANOVA con unweighted unifrac
uwUF.dist = UniFrac(ps.rarefied, weighted=FALSE, normalized=TRUE)
adonis2(uwUF.dist ~ condicion*especie*sexo*tratamiento, data=meta(ps.rarefied),
        permutations = 1000)


# Beta-dispersion por especie
disp.spp = betadisper(wUF.dist, meta$especie)
permutest(disp.spp, pairwise=TRUE, permutations=1000)
# Beta-dispersion por sexo
disp.sexo = betadisper(wUF.dist, meta$sexo)
permutest(disp.sexo, pairwise=TRUE, permutations=1000)
# Beta-dispersion por tratamiento
disp.trat = betadisper(wUF.dist, meta$tratamiento)
permutest(disp.trat, pairwise=TRUE, permutations=1000)


# checando homogeneidad
dist1 <- vegdist(t(otu))
anova(betadisper(dist1, meta$especie))
dist2 <- vegdist(t(otu))
anova(betadisper(dist2, meta$sexo))
dist3 <- vegdist(t(otu))
anova(betadisper(dist3, meta$tratamiento))

?betadisper
#####################################################################33
# Diversidad beta de todas las muestras sin separar por tratamiento
ordu.unwt.uni <- ordinate(phys4, "PCoA", "unifrac", weighted=F)

# check for Eigen values 
# barplot(ordu.unwt.uni$values$Eigenvalues[1:10])
# por especie  
unwt.unifrac <- plot_ordination(phys4, 
                                ordu.unwt.uni, color="especie", 
                                shape = "tratamiento") 
unwt.unifrac <- unwt.unifrac + ggtitle("Unweighted UniFrac") +
  geom_point(size = 2)
unwt.unifrac <- unwt.unifrac + theme_classic() +
  scale_color_brewer("Location", palette = "Set2")
print(unwt.unifrac)  

# por sexo
unwt.unifrac2 <- plot_ordination(phys4, 
                                 ordu.unwt.uni, color="sexo") 
unwt.unifrac2 <- unwt.unifrac2 + ggtitle("Unweighted UniFrac") +
  geom_point(size = 2)
unwt.unifrac2 <- unwt.unifrac2 + theme_classic() +
  scale_color_brewer("Location", palette = "Set2")
print(unwt.unifrac2)    

# por tratamiento  
unwt.unifrac3 <- plot_ordination(phys4, 
                                 ordu.unwt.uni, color="tratamiento") 
unwt.unifrac3 <- unwt.unifrac3 + ggtitle("Unweighted UniFrac") +
  geom_point(size = 2)
unwt.unifrac3 <- unwt.unifrac3 + theme_classic() +
  scale_color_brewer("Location", palette = "Set2")
print(unwt.unifrac3)  

############## Weighted Unifrac
ps1.rel <- microbiome::transform(phys4, "compositional")

ordu.wt.uni <- ordinate(ps1.rel , "PCoA", "unifrac", weighted=T)

# check for Eigen values 
# barplot(ordu.unwt.uni$values$Eigenvalues[1:10])

# PCoA de la especie
wt.unifrac <- plot_ordination(ps1.rel, 
                              ordu.wt.uni, color="especie") 
wt.unifrac <- wt.unifrac + ggtitle("Weighted UniFrac") +
  geom_point(size = 3)
wt.unifrac <- wt.unifrac + theme_classic() + 
  scale_color_brewer("Especies", palette = "Set2")+
  theme(legend.key.size = unit(1.5, 'cm'))
print(wt.unifrac)
print(wt.unifrac + stat_ellipse(linetype = 1, cex = 1.5))


# por sexo  
wt.unifrac2 <- plot_ordination(ps1.rel, 
                               ordu.wt.uni, color = "sexo") 
wt.unifrac2 <- wt.unifrac2 + ggtitle("Weighted UniFrac") +
  geom_point(size = 2)
wt.unifrac2 <- wt.unifrac2 + theme_classic() + 
  scale_color_brewer("Sexo", palette = "Set2",
                     labels = c("Hembras", "Machos"))+
  theme(legend.key.size = unit(1.5, 'cm'))
print(wt.unifrac2)  
print(wt.unifrac2 + stat_ellipse(linetype = 1, cex = 1.5))
# por tratamiento
wt.unifrac3 <- plot_ordination(ps1.rel, 
                               ordu.wt.uni, color="tratamiento") 
wt.unifrac3 <- wt.unifrac3 + ggtitle("Weighted UniFrac") +
  geom_point(size = 3)
wt.unifrac3 <- wt.unifrac3 + theme_classic() + 
  scale_color_brewer("Tratamientos", palette = "Set2",
                     labels = c("No visitadas", "Visitadas"))+
  theme(legend.key.size = unit(1.5, 'cm'))
print(wt.unifrac3) 
print(wt.unifrac3 + stat_ellipse(linetype = 1, cex = 1.5))


### PCoA combinando variables: la especie por el tratamiento
wt.unifrac4 <- plot_ordination(ps1.rel, 
                               ordu.wt.uni, color="especie", 
                               shape ="tratamiento") 
wt.unifrac4 <- wt.unifrac4 + ggtitle("Weighted UniFrac") +
  geom_point(size = 3)
wt.unifrac4 <- wt.unifrac4 + theme_classic() + 
  scale_color_brewer("Especies", palette = "Set2")+
  scale_shape_discrete(name = "Tratamientos", labels = c("No visitadas",
                                                         "Visitadas"))+
  theme(legend.key.size = unit(1.5, 'cm'))
print(wt.unifrac4)
print(wt.unifrac4 + stat_ellipse(linetype = 1, cex = 1.1))

### PCoA combinando variables: la especie por el tratamiento
wt.unifrac5 <- plot_ordination(ps1.rel, 
                               ordu.wt.uni, color="especie", 
                               shape ="sexo") 
wt.unifrac5 <- wt.unifrac5 + ggtitle("Weighted UniFrac") +
  geom_point(size = 3)
wt.unifrac5 <- wt.unifrac5 + theme_classic() + 
  scale_color_brewer("Especies", palette = "Set2")+
  scale_shape_discrete(name = "Sexos", labels = c("Hembras",
                                                  "Machos"))
print(wt.unifrac5)
print(wt.unifrac5 + stat_ellipse(linetype = 1, cex = 1.1))

### PCoA combinando variables: el sexo por el tratamiento
wt.unifrac7 <- plot_ordination(ps1.rel, 
                               ordu.wt.uni, color="sexo", 
                               shape ="tratamiento") 
wt.unifrac7 <- wt.unifrac7 + ggtitle("Weighted UniFrac") +
  geom_point(size = 3)
wt.unifrac7 <- wt.unifrac7 + theme_classic()

+ 
  scale_color_brewer("Especies", palette = "Set2")+
  scale_shape_discrete(name = "Sexos", labels = c("Hembras",
                                                  "Machos"))
print(wt.unifrac7)
print(wt.unifrac7 + stat_ellipse(linetype = 1, cex = 1.1))


library(vegan)
metadf <- data.frame(sample_data(ps1.rel))
unifrac.dist <- UniFrac(ps1.rel, 
                        weighted = TRUE, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)

permanova <- adonis2(unifrac.dist ~ especie, data = metadf)
permanova

metadata <- sample_data(food) %>% as("data.frame")

model <- adonis2(unifrac.dist ~ especie * tratamiento, data = metadf,
                 permutations = 999)
model

ps.disper <- betadisper(unifrac.dist, metadf$tratamiento)
permutest(ps.disper, pairwise = TRUE)  

ps.disper <- betadisper(unifrac.dist, metadf$especie)
permutest(ps.disper, pairwise = TRUE)





###############################################################################
##################### METODOS DE NORMALIZACION ################################
#### RAREFIED DATA ####
# Curvas de rarefacción
otu_tab_its = t(abundances(its6))
p = vegan::rarecurve(otu_tab_its,
                     step = 50, label = FALSE,
                     sample = min(rowSums(otu_tab_its),
                                  col = "blue, cex = 0.6"))

### Normalize library size to 1000 reads
set.seed(111)
rare.its1000 <- rarefy_even_depth(its6, 
                                     1000, rngseed = 1,
                                     replace = FALSE)             
rare.its1000
# phyloseq object with normalized library size has
# 17 ASVs in 4 samples

### Normalize library size to 100 reads
set.seed(111)
rare.its100 <- rarefy_even_depth(its6, 
                                  100, rngseed = 1,
                                  replace = FALSE)             
rare.its100
# phyloseq object with normalized library size has
# 72 ASVs in 24 samples

### Normalize library size to 50 reads
set.seed(111)
rare.its50 <- rarefy_even_depth(its6, 
                                 50, rngseed = 1,
                                 replace = FALSE)             
rare.its50
# phyloseq object with normalized library size has
# 85 ASVs in 36 samples

### Normalize library size to 10 reads
set.seed(111)
rare.its10 <- rarefy_even_depth(its6, 
                                10, rngseed = 1,
                                replace = FALSE)             
rare.its10
# phyloseq object with normalized library size has
# 109 ASVs in 100 samples
reads_sampleits = readcount(rare.its10)
reads_sampleits[1:20]
summarize_phyloseq(rare.its10)
# total number of reads= 1000 in rarefied data
summary(sample_sums(rare.its10))
abundances(rare.its10)

### TAXA PREVALENCE
plot_taxa_prevalence(its6, "Class")
plot_taxa_prevalence(rare.its10, "Class")



           #### ALPHA DIVERSITY WITHOUT RAREFIED DATA  ####
library(knitr)

div.its = microbiome::alpha(its6, index = c("observed", "shannon"))
kable(head(div.its))
View(div.its)
str(div.its)

# get the metadata out as separate object
its.meta <- meta(its6)
kable(head(its.meta))
# Add the rownames as a new colum for easy integration later.
its.meta$sam_name <- rownames(its.meta)

# Add the rownames to diversity table
div.its$sam_name <- rownames(div.its)

# merge these two data frames into one
div.its.df <- merge(div.its,its.meta, by = "sam_name")

# check the tables
colnames(div.its.df)

### Phylogenetic diversity (Faith)
library(picante)

ps0.its.asvtab <- as.data.frame(its6@otu_table)

ps0.its.tree <- its6@phy_tree

# hmp.meta from previous code chunks

# We first need to check if the tree is rooted or not 
its6@phy_tree

# it is a rooted tree
its.pd <- pd(t(ps0.its.asvtab), ps0.its.tree,include.root=T) 
# t(ou_table) transposes the table for use in picante and the tre file comes
# from the first code chunck we used to read tree file (see making a phyloseq 
# object section).
# datatable(df.pd)

# We add the results of PD to this file and then plot
div.its.df$Phylogenetic_Diversity <- its.pd$PD
head(div.its.df)


# data
div.its.df$observed
div.its.df$diversity_shannon
div.its.df$Phylogenetic_Diversity
# Histograma
hist(div.its.df$observed)
hist(div.its.df$diversity_shannon) 
hist(div.its.df$Phylogenetic_Diversity)
# Shapiro test
shapiro.test(div.its.df$observed)
shapiro.test(div.its.df$diversity_shannon)
shapiro.test(div.its.df$Phylogenetic_Diversity)

### STATISTICAL ANALYSIS OF ALPHA DIVERSITY
library(emmeans)
# Richness
mod7 = glm(observed ~ especie*sexo, data = div.its.df,
           family = poisson)
summary(mod7)
anova(mod7)
#mean(div.its.df$observed)
#var(div.its.df$observed)
# post-hoc comparisons (all variables)
compa7 = emmeans(mod7, 
                 specs = pairwise ~ especie, 
                 type = "response")
compa7
compa7.1 = emmeans(mod7,
                   specs = pairwise ~ especie | sexo,
                   type = "response")
compa7.1
#write.csv(compa1$contrasts, file = "observed_compa_rarefied.csv")

# Shannon
mod8 = glm(diversity_shannon ~ especie*sexo, data = div.its.df)
summary(mod8)
anova(mod8)
# post-hoc comparisons
compa8 = emmeans(mod8, 
                 specs = pairwise ~ especie, 
                 type = "response")
compa8
compa2.1 = emmeans(mod2, 
                   specs = pairwise ~ especies, 
                   type = "response")
compa2.1
#modx = aov(diversity_shannon ~ especies*sexo*tratamiento, data = div.df)
#summary(modx)

# Phylogenetic Diversity index
mod9 = glm(Phylogenetic_Diversity ~ especie*sexo, data = div.its.df)
summary(mod9)
anova(mod9)
# post-hoc comparisons
compa9 = emmeans(mod9, 
                 specs = pairwise ~ especie, 
                 type = "response")
compa9
compa3.1 = emmeans(mod3, 
                   specs = pairwise ~ especies, 
                   type = "response")
compa3.1

##### VISUALIZATION  OF ALPHA DIVERSITY without rarefied data####
# Now use this data frame to plot 
p.shannon.its <- boxplot_alpha(its6, 
                           index = "shannon",
                           x_var = "especie")
p.shannon.its
#fill.colors = c(female="cyan4", male="deeppink4"))

p.shannon <- p.shannon + theme_minimal() + 
  labs(x="Species", y="Shannon diversity") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16))
p.shannon

# Stripchart of observed ASVs with 1 panel
ggplot(div.its.df, aes(x=especie, y=observed, color=especie)) +
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

riqueza_its <- ggplot(div.its.df, aes(x=especie, y=observed, color=especie)) +
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
riqueza_its
### stripchart Shannon diversity index facet with 2 panels (divided by sex)
shannon_its <- ggplot(div.its.df, aes(x=especie, y=diversity_shannon, 
                                      color=especie))+ 
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
shannon_its

cartel2 =  ggarrange(riqueza_its, shannon_its,
                     nrow = 2, labels = "AUTO")
cartel2
### stripchart phylogenetic diversity index facet with 2 panels (divided by sex)
pd_its <- ggplot(div.its.df, aes(x=especie, y=Phylogenetic_Diversity,
                            color=especie))+ 
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
pd_its
# cargando paqueteria gridExtra
library(gridExtra)
library(cowplot)
gt.its <- grid.arrange(riqueza_its, shannon_its, pd_its, 
                   ncol = 1, nrow = 3)


########## BETA DIVERSITY  WITH RAREFIED DATA ##########
library(vegan)
#### Unweighted UniFrac ####
# PCoA plot using the unweighted UniFrac as distance
ordination_its = ordinate(its6, method="PCoA", "unifrac",
                          weighted = F)

baseplot_its <- plot_ordination(its6, ordination_its,
                                color="especie",
                                type = "samples",
) + #shape = "sexo"
  ggtitle("")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 2)
baseplot_its
#baseplot1 +
#  geom_point(size = 1) +
#  geom_label(aes(label = clave), size = 3.5, color = "black")
#
baseplot_its +
  geom_point(size = 2)+
  ggtitle("")+
  facet_wrap(~sexo, labeller = labeller(sexo = sexos))+
  stat_ellipse()

#### Weighted UniFrac ####
# Transform data to relative abundance
its6
its_re = microbiome::transform(its6, "compositional")
# PCoA plot using the weighted UniFrac as distance
ordination_its_2 = ordinate(its_re, method="PCoA", "wunifrac")
# Plot
baseplot_its_2 <- plot_ordination(its_re, ordination_its_2,
                                  color="especie",
                                  type = "samples") + #shape = "sexo" shape = "tratamiento"
  ggtitle("Hongos")+
  theme(aspect.ratio=1)+
  theme_bw()+
  geom_point(size = 3)+
scale_color_discrete(name = "Especies", 
                     labels = c("1CF" = "CF","2CPF" = "CPF","3CAS" = "CAS",
                                "4CPP"="CPP", "5CAA"="CAA", "6CM"="CM"))
#labs(x = "PCoa1 (31.3%)", y = "PCoA2 (10.4%)")+
#scale_color_discrete(name = "Condition", 
#                    labels = c("Domesticated","Wild"))+
#scale_color_discrete(name = "Treatment",
#                     labels = c("No-visited","Visited"))

baseplot_its_2
#
wei_its = baseplot_its_2 + stat_ellipse(type = "norm", linetype = 2)
#
#baseplot2 +
#  geom_point(size = 1) +
#  geom_label(aes(label = clave), size = 3.5, color = "black")
#
baseplot_its_2 +
  geom_point(size = 2)+
  ggtitle("")+
  facet_wrap(~sexo, labeller = labeller(sexo = sexos))+
  stat_ellipse(type = "norm", linetype = 2)


#### PERMANOVA ####
its_re
metadf = data.frame(sample_data(its_re))
# calculating Weighted-UniFrac
wunifrac.dist.its <- UniFrac(its_re, 
                         weighted = TRUE, 
                         normalized = TRUE,  
                         parallel = FALSE, 
                         fast = TRUE)
# Run PERMANOVA
set.seed(100)   # Set the seed to make the analysis reproducible
perm5 = adonis2(wunifrac.dist.its ~ especie*sexo, data = metadf,
                by = "terms", permutations = 999)
perm5        

# calculating Unweighted-UniFrac
unifrac.dist.its <- UniFrac(its_re, 
                        weighted = FALSE, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)
# Run PERMANOVA
set.seed(100)   # Set the seed to make the analysis reproducible
perm6 = adonis2(unifrac.dist.its ~ especie*sexo, data = metadf,
                by = "terms", permutations = 999)
perm6  

?betadisper
ps.disper <- betadisper(unifrac.dist, metadf$especies)
permutest(ps.disper, pairwise = TRUE)
anova(ps.disper)
TukeyHSD(ps.disper)






###############################################################################
#### Ahora vamos a utilizar el paquete "microbiome" ####
# Primero nos aseguramos que estamos usando las funciones del paquete correcto
transform <- microbiome::transform

# Primero resumimos el contenido del objeto phyloseq
summarize_phyloseq(ps.rarefied)
# información del número de reads
summary(sample_sums(ps.rarefied))

#### Fusionamos taxones poco abundantes para facilitar la visualización ####
ps.rarefied
pseq <- transform(ps.rarefied, "compositional")
pseq2 <- aggregate_rare(pseq, level = "Family", 
                        detection = 1/31, prevalence = 1/31)
pseq2

# Visualizamos las 14 familias mas abundantes más "Otros" y "Unknown"
nb.cols <- 14
mycolors <- colorRampPalette(brewer.pal(14, "Paired"))(nb.cols)

plot_composition(pseq2,
                 taxonomic.level = "Family",
                 sample.sort = "especie",
                 x.label = "Especies") +
  scale_fill_brewer("Familia", palette = "Paired") +
  guides(fill = guide_legend(ncol = 1)) +
  scale_y_percent() +
  labs(x = "Muestras", y = "Abundancia Relativa (%)",
       title = "Datos de abundancia relativa") + 
  theme_ipsum(grid="Y") +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.text = element_text(face = "italic"))#+
#  scale_fill_manual(values = mycolors)

plot_composition(pseq2,
                 taxonomic.level = "Family",
                 sample.sort = "especies",
                 x.label = "Especies") +
  scale_fill_brewer("Familia", palette = "Paired") +
  guides(fill = guide_legend(ncol = 1)) +
  scale_y_percent() +
  labs(x = "Muestras", y = "Abundancia Relativa (%)",
       title = "Datos de abundancia relativa") + 
  theme_ipsum(grid="Y") +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.text = element_text(face = "italic"))#+
#  scale_fill_manual(values = mycolors)


# Agrupadas por especie
plot_composition(pseq2,
                 average_by = "especie", 
                 transform = "compositional") +
  scale_fill_brewer("Familia", palette = "Paired") +
  theme_ipsum(grid="Y") +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.text = element_text(face = "italic"))#+
  #scale_fill_manual(values = mycolors)

# Agrupadas por sexo
plot_composition(pseq2,
                 average_by = "sexo", 
                 transform = "compositional") +
  scale_fill_brewer("Familia", palette = "Paired") +
  theme_ipsum(grid="Y") +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.text = element_text(face = "italic"))#+
#  scale_fill_manual(values = mycolors)

# Agrupadas por tratamiento
plot_composition(pseq2,
                 average_by = "tratamiento", 
                 transform = "compositional") +
  scale_fill_brewer("Familia", palette = "Paired") +
  theme_ipsum(grid="Y") +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.text = element_text(face = "italic"))#+
#  scale_fill_manual(values = mycolors)

#### Grafica de abundancia relativa por especie ####
pseq2 %>%  plot_composition(average_by = "especie")+
  scale_fill_brewer("Familias", palette = "Paired")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.text = element_text(face = "italic"))+
  labs(x = "Especies", y = "Abundancia Relativa")+
  theme_bw()+
  scale_x_discrete(labels=c("1CAS" = "CAS", "2CAA" = "CAA",
                            "3CM" = "CM", "4CPP" = "CPP"))+
  guides(fill = guide_legend(ncol = 1))+
  theme_minimal()

#### Grafica de abundancia relativa por sexo ####
pseq2 %>%  plot_composition(average_by = "sexo")+
  scale_fill_brewer("Familia", palette = "Paired")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.text = element_text(face = "italic"))+
  labs(x = "Sexo floral", y = "Abundancia Relativa")+
  theme_bw()+
  scale_x_discrete(labels=c("Hembras","Machos"))+
  guides(fill = guide_legend(ncol = 1))+
  theme_minimal()

#### Grafica de abundancia relativa por tratamiento ####
pseq2 %>%  plot_composition(average_by = "tratamiento")+
  scale_fill_brewer("Familia", palette = "Paired")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.text = element_text(face = "italic"))+
  labs(x = "Tratamientos", y = "Abundancia Relativa")+
  theme_bw()+
  scale_x_discrete(labels=c("No visitadas","Visitadas"))+
  guides(fill = guide_legend(ncol = 1))+
  theme_minimal()

#+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#        panel.background = element_blank(), axis.line = element_line(colour = "black"))

###############################################################################
#### Creando una gráfica que combina especies, sexo y tratamientos
# filtrando secuencias entre sexos y tratamientos
t1_m <- subset_samples(pseq2, sexo == "macho" & tratamiento == "t1");t1_m  
t2_m <- subset_samples(pseq2, sexo == "macho" & tratamiento == "t2");t2_m

t1_h <- subset_samples(pseq2, sexo == "hembra" & tratamiento == "t1");t1_h  
t2_h <- subset_samples(pseq2, sexo == "hembra" & tratamiento == "t2");t2_h

nb.cols <- 10
mycolors <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols)

# Flores macho del tratamiento 1
plot_t1_m  <- plot_composition(t1_m,
                               average_by = "especie") + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(x = "Flores sin visita", y = "")+ 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 0.5),
        legend.position = "none",
        axis.text.y = element_blank())+
  scale_fill_manual(values = mycolors)+
  theme(axis.text.x = element_text(face = "bold"))+ 
  #        color = c("#9E6334", "#005200", "#9E6334", "#9E6334"))) 
  scale_x_discrete(labels=c("1CAS" = "CAS", "2CAA" = "CAA",
                            "3CM" = "CM", "4CPP" = "CPP"))
plot_t1_m

# Flores macho del tratamiento 2
plot_t2_m  <- plot_composition(t2_m,
                               average_by = "especie") + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(x = "Flores con visita", y = "")+ 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 0.5),
        legend.position = "none",
        axis.text.y = element_blank())+
  scale_fill_manual(values = mycolors)+
  theme(axis.text.x = element_text(face = "bold"))+
  scale_x_discrete(labels=c("1CAS" = "CAS", "2CAA" = "CAA",
                            "3CM" = "CM", "4CPP" = "CPP"))
plot_t2_m

#### Flores hembra del tratamiento 1 #####
plot_t1_h  <- plot_composition(t1_h,
                               average_by = "especie") + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(x = "Flores sin visita", y = "")+ 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 0.5),
        legend.position = "none",
        axis.text.y = element_blank())+
  scale_fill_manual(values = mycolors)+
  theme(axis.text.x = element_text(face = "bold"))+
  scale_x_discrete(labels=c("1CAS" = "CAS", "2CAA" = "CAA",
                            "3CM" = "CM", "4CPP" = "CPP"))
plot_t1_h

# Flores hembra del tratamiento 2
plot_t2_h  <- plot_composition(t2_h,
                               average_by = "especie") + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(x = "Flores con visita", y = "")+ 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 0.5),
        legend.position = "none",
        axis.text.y = element_blank())+
  scale_fill_manual(values = mycolors)+
  theme(axis.text.x = element_text(face = "bold"))+
  scale_x_discrete(labels=c("1CAS" = "CAS", "2CAA" = "CAA",
                            "3CM" = "CM", "4CPP" = "CPP"))
plot_t2_h

# Combinando las gráficas con el paquete "ggpubr"
m <- ggarrange(plot_t1_m, plot_t2_m, nrow = 1, widths = c(1.2,1,1))
h <- ggarrange(plot_t1_h, plot_t2_h, nrow = 1, widths = c(1.2,1,1))

ggarrange(m, h, nrow = 2, labels = c("Flores masculinas","Flores femeninas"),
          vjust = 1,font.label = list(size = 12, face = "bold"))

# Flores macho tratamiento 1
plot_t1_m  <- plot_composition(t1_m,
                               average_by = "especie") + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(x = "Tratamiento 1",  y = "Abundancia relativa")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 0.5),
        axis.text.y = element_blank())+
  scale_fill_manual(values = mycolors, name = "Género")

# Extract the legend. Returns a gtable
leyenda <- get_legend(plot_t1_m)

# Convert to a ggplot and print
as_ggplot(leyenda)

### Heatmap
plot_heatmap(ps1)
plot_heatmap(its2_5)

p <- 
  microbiome::transform(its2_5, "compositional") %>% 
  plot_composition(plot.type = "heatmap",
                   sample.sort = "neatmap", 
                   otu.sort = "neatmap") +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.x = element_text(size = 9, hjust=1),
        legend.text = element_text(size = 8)) +
  ylab("Samples")
print(p)

### Heatmap 2
sample_variables(its2_3)
gpt = prune_taxa(names(sort(taxa_sums(its2_3), TRUE) [1:20]), its2_3)
gpt = prune_samples(sample_sums(gpt)>1, gpt)
plot_heatmap(gpt, sample.label = "especie")


### Plot taxa prevalence
its2_5
# Define detection and prevalence thresholds to filter out rare taxa
p0 <- core(its2_5, detection = 0.1/100, prevalence = 1/100)

# For the available taxonomic levels
plot_taxa_prevalence(p0, "Phylum", detection = 0.1/100)
plot_taxa_prevalence(p0, "Class", detection = 0.1/100)
plot_taxa_prevalence(p0, "Order", detection = 0.1/100)
plot_taxa_prevalence(p0, "Family", detection = 0.1/100)

###############################################################################
library(rentrez)

# Ver cuántas secuencias existen en total
query_all_fungi <- 'ITS2[All Fields] AND Fungi[Organism]'

busqueda_total <- entrez_search(
  db          = "nucleotide",
  term        = query_all_fungi,
  retmax      = 0,
  use_history = TRUE
)

cat("Total secuencias ITS2 de hongos en GenBank:", busqueda_total$count, "\n")



library(rentrez)

# Ver cuántas hay por phylum principal
phylums <- list(
  Ascomycota    = 'ITS2[All Fields] AND Ascomycota[Organism]',
  Basidiomycota = 'ITS2[All Fields] AND Basidiomycota[Organism]',
  Mucoromycota  = 'ITS2[All Fields] AND Mucoromycota[Organism]',
  Chytridiomycota = 'ITS2[All Fields] AND Chytridiomycota[Organism]',
  Glomeromycota = 'ITS2[All Fields] AND Glomeromycota[Organism]'
)

for(nombre in names(phylums)){
  r <- entrez_search(db = "nucleotide", term = phylums[[nombre]],
                     retmax = 0, use_history = TRUE)
  cat(nombre, ":", r$count, "\n")
  Sys.sleep(0.3)
}



###########
# txid4751 = Fungi (ID taxonómico directo)
# equivalente a Fungi[Organism] pero más preciso

query_its_hongos <- 'txid4751[ORGN] AND 
                     (internal transcribed spacer ITS2[Title] OR 
                      ITS2[Title] OR 
                      internal transcribed spacer 2[Title])'

busqueda_hongos <- entrez_search(
  db          = "nucleotide",
  term        = query_its_hongos,
  retmax      = 0,
  use_history = TRUE
)
cat("Total ITS2 hongos (query con txid):", busqueda_hongos$count, "\n")

# Comparar con tu query anterior
query_anterior <- 'ITS2[All Fields] AND Fungi[Organism]'
busqueda_anterior <- entrez_search(
  db     = "nucleotide",
  term   = query_anterior,
  retmax = 0,
  use_history = TRUE
)
cat("Total ITS2 hongos (query anterior):", busqueda_anterior$count, "\n")



##############################################################################
library(rentrez)

# Query mejorado con txid
query_hongos_completo <- 'txid4751[ORGN] AND 
                          (ITS2[Title] OR 
                           "internal transcribed spacer ITS2"[Title] OR
                           "internal transcribed spacer 2"[Title])'

# Ver total
busqueda <- entrez_search(
  db          = "nucleotide",
  term        = query_hongos_completo,
  retmax      = 0,
  use_history = TRUE
)
cat("Total:", busqueda$count, "\n")


## 2
# cambiar directorio
setwd("/home/luis/")

library(Biostrings)
library(dplyr)

# Paso 1: Leer
cat("Leyendo archivo...\n")
seqs_raw <- readDNAStringSet("hongos_ITS2_completo.fasta")
cat("Total descargadas:", length(seqs_raw), "\n")

# Paso 2: Eliminar duplicados por ID
seqs_unicas <- seqs_raw[!duplicated(names(seqs_raw))]
cat("Secuencias únicas:", length(seqs_unicas), "\n")

# Paso 3: Filtrar por longitud
# guardar amplio para ITSx (200-900 bp)
seqs_long <- seqs_unicas[
  width(seqs_unicas) >= 200 &
    width(seqs_unicas) <= 900
]
cat("Después de filtrar por longitud:", length(seqs_long), "\n")

# Paso 4: Filtrar Ns excesivos
prop_n <- letterFrequency(seqs_long, "N") / width(seqs_long)
seqs_limpias <- seqs_long[prop_n < 0.10]
cat("Después de filtrar Ns:", length(seqs_limpias), "\n")

# Paso 5: Ver distribución de longitudes
summary(width(seqs_limpias))
hist(width(seqs_limpias),
     breaks = 50,
     main   = "Distribución de longitudes ITS hongos",
     xlab   = "Longitud (bp)",
     col    = "steelblue",
     border = "white")
abline(v = c(200, 500), col = "red", lty = 2)

# Paso 6: Guardar versión limpia para ITSx
writeXStringSet(seqs_limpias, "hongos_ITS2_genbank_clean.fasta")
cat("✓ Archivo listo para ITSx\n")

##############################################################################
# Prueba para detectar secuencias iguales con nombres diferentes
library(Biostrings)
library(dplyr)

### change directory
setwd("/home/luis/Documentos/papers_publicados/paper2_tesisDoctorado/Cucurbita_microbiota_nectar")

# Leer la base combinada
base <- readDNAStringSet("base_ITS2_UNITE_GenBank_DADA2.fasta")
cat("Total secuencias:", length(base), "\n")

# Extraer secuencias como strings
seqs    <- as.character(base)
headers <- names(base)

# Crear dataframe con secuencia y taxonomía
df_base <- data.frame(
  header   = headers,
  taxonomy = headers,   # en DADA2 el header ES la taxonomía
  sequence = seqs,
  stringsAsFactors = FALSE
)

head(df_base)

### Paso 2: Detectar secuencias duplicadas con diferente taxonomía
# Agrupar por secuencia y contar cuántas taxonomías diferentes tiene
duplicados_tax <- df_base %>%
  group_by(sequence) %>%
  summarise(
    n_taxonomias    = n_distinct(taxonomy),
    taxonomias      = paste(unique(taxonomy), collapse = " | "),
    n_ocurrencias   = n(),
    .groups = "drop"
  ) %>%
  filter(n_taxonomias > 1) %>%  # solo las que tienen >1 taxonomía
  arrange(desc(n_taxonomias))

cat("Secuencias con taxonomía conflictiva:", nrow(duplicados_tax), "\n")
print(duplicados_tax)
View(duplicados_tax)
write.csv(duplicados_tax, file = "duplicados_prueba.csv")


### Paso 3:
# Ver cuántos conflictos hay por nivel taxonómico
if(nrow(duplicados_tax) > 0) {
  
  # Separar taxonomías para comparar nivel por nivel
  conflictos_detalle <- duplicados_tax %>%
    rowwise() %>%
    mutate(
      tax_list = list(strsplit(taxonomias, " \\| ")[[1]])
    ) %>%
    ungroup()
  
  # Ver los primeros conflictos
  cat("\nPrimeros conflictos encontrados:\n")
  for(i in 1:min(5, nrow(duplicados_tax))){
    cat("\n--- Conflicto", i, "---\n")
    cat("Secuencia:", substr(duplicados_tax$sequence[i], 1, 50), "...\n")
    cat("Taxonomías:\n")
    tax_split <- strsplit(duplicados_tax$taxonomias[i], " \\| ")[[1]]
    for(t in tax_split) cat(" →", t, "\n")
  }
}


CAACCCTCGAACCCCTCCGGGGGGTCGGCGTTGGGGATCGGCCCTCCCTTAGCGGGTGGCCGTCTCCGAAATACAGTGGC

GGTCTCGCCGCAGCCTCTCCTGCGCAGTAGTTTGCACACTCGCATCGGGAGCGCGGCGCGTCCACAGCCGTTAAACACCC
AACTTCTGAAA
