# Pipeline for bioinformatic analysis of 16s sequences of Cucurbita species

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

# Change to directory where there are files
path <- "/home/.lvillanueva/work/analisis_2024/16s"
list.files(path)

# Forward and reverse fastq filenames have format: 
# SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
head(fnFs)
head(fnRs)

# Primers of 16s microbiome earth project
fwd <- "GTGYCAGCMGCCGCGGTAA"
rev <- "GGACTACNVGGGTWTCTAAT"

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# For this dataset, we will use standard filtering paraments: maxN=0 
# (DADA2 requires sequences contain no Ns), truncQ = 2, rm.phix = TRUE and 
# maxEE=2. The maxEE parameter sets the maximum number of “expected errors” 
# allowed in a read, which is a better filter than simply averaging quality 
# scores. Note: We enforce a minLen here, to get rid of spurious very low-length
# sequences. This was not needed in the 16S Tutorial Workflow because truncLen
# already served that purpose.

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,200),
                     trimLeft=c(20,20), maxN=0, maxEE=c(2,2), truncQ=2, 
                     rm.phix=TRUE, compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

# If I a relax parameters maxEE and truncQ
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen = c(240,200),
                     maxN = 0, maxEE = 8,
                     truncQ = 8, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)

# Show how many sequences passed filtering
head(out)

# Learning the error rates:
errF <- learnErrors(filtFs, multithread=TRUE)
# 103181760 total bases in 469008 reads from 15 samples will be used for
# learning the error rates.
errR <- learnErrors(filtRs, multithread=TRUE)
# 105716520 total bases in 587314 reads from 18 samples will be used for
# learning the error rates.

# Plotting the error rates:
#plotErrors(errF, nominalQ=TRUE)

# Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)

# Distribution of sequence lengths after removing chimeras
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))

# Exporting sequence table
write.csv(seqtab.nochim, file = "seqtab_nochim_16s_2024.csv")
head(seqtab.nochim)

# Abundances of quimeras of the merged sequences reads
sum(seqtab.nochim)/sum(seqtab)
# 0.9757 only 2.5% of sequences were quimeras

# Inspect distribution of sequence lengths without quimeras
table(nchar(getSequences(seqtab.nochim)))

# Statistics of filtering
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace 
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged",
                     "nonchim")
rownames(track) <- sample.names
head(track)

# Exporting statistics of filtering
write.csv(track, file = "resumen_filtrado_16s_2024.csv")

## cuando lo hago en la compu
#setwd("/home/luis/datos_16s_setC_may24")
#seqtab.nochim = read.csv("seqtab_nochim_16s_2024.csv", row.names = 1)
#seqtab.nochim2 =  as.matrix(seqtab.nochim)
#head(seqtab.nochim2)

### Assign taxonomy
silva <- "/home/.lvillanueva/work/analisis_2024/silva_nr99_v138.1_train_set.fa.gz"
taxa <- assignTaxonomy(seqtab.nochim, silva, multithread=TRUE)

# extension
taxa2 <- addSpecies(taxa, "/home/.lvillanueva/work/analisis_2024/silva_species_assignment_v138.1.fa.gz")

# Taxonomic assignments
head(taxa)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
write.csv(taxa.print, file = "taxa_16s_2024.csv")

# Taxonomic assignments by species
head(taxa2)
taxa2.print <- taxa2 # Removing sequence rownames for display only
rownames(taxa2.print)
head(taxa2.print)
write.csv(taxa2.print, file = "taxa_16s_species_2024.csv")

