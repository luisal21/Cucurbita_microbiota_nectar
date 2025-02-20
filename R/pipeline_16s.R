# Pipeline análisis bioinformatico de ITS2 de néctar floral de cucurbita

# Luis Alberto Villanueva Espino
# Modified from: https://benjjneb.github.io/dada2/ITS_workflow.html

# Entrar al cluster del LANASE
# dirección
# IP: 132.247.186.43
# Puerto: 60307

# Conexión al cluster
#ssh villanueva@132.247.186.43 -p 60307

# path en el cluster
# /home/.lvillanueva/work/analisis_2024
# /home/.lvillanueva/work/analisis_2024/16s

# Activar ambiente de conda - Bioconductor
# Activar R dentro del ambiente Bioconductor


require(parallel)
#
cl  <- makeCluster(4)

# primero activar R en bioconductor

# Cargamos los paquetes requeridos
library(dada2)
library(ShortRead)
library(Biostrings)

# Moverme al directorio donde se encuentran las muestras:
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

# Primers usados-primers 16s microbiome earth project
fwd <- "GTGYCAGCMGCCGCGGTAA"
rev <- "GGACTACNVGGGTWTCTAAT"

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
fwd.orients <- allOrients(fwd);fwd.orients
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

# Conteo de primers
primers_count <- rbind(FWD.ForwardReads = sapply(fwd.orients, primerHits,
                 fn = fnFs.filtN[[1]]), FWD.ReverseReads = sapply(fwd.orients,
                 primerHits, fn = fnRs.filtN[[1]]), 
                       REV.ForwardReads = sapply(rev.orients, primerHits,
                 fn = fnFs.filtN[[1]]), REV.ReverseReads = sapply(rev.orients,
                 primerHits, fn = fnRs.filtN[[1]]))
# Revisar conteo de primers
primers_count

### Removiendo primers con cutadapt
cutadapt <- "/usr/local/bin/cutadapt" # CHANGE ME to the cutadapt path on your
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
FWD.RC <- dada2:::rc(fwd)
# Complemento reverso del Reverse
REV.RC <- dada2:::rc(rev)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", fwd, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", rev, "-A", FWD.RC)

# Run Cutadapt
for(i in seq_along(fnFs.filtN)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to 
                             # remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# As a sanity check, we will count the presence of primers in the first 
# cutadapt-ed sample:
rbind(FWD.ForwardReads = sapply(fwd.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(fwd.orients, primerHits, fn = fnRs.cut[[1]]), 
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
#plotQualityProfile(cutFs[1:2])

### Filtrado y corte (trim):
# Colocar archivos filtrados en un subdirectorio
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

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
# filterAndTrim relajando los parametros de maxEE y truncQ
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen = c(240,200),
                     maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE

# Relajando parametros maxEE and truncQ
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen = c(240,200),
                     maxN = 0, maxEE = 8,
                     truncQ = 8, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE


# muestra cuantas secuencias pasaron el filtrado
head(out)

# Aprendiendo las tasas de error
errF <- learnErrors(filtFs, multithread=TRUE)
# 101668613 total bases in 438036 reads from 15 samples will be used for 
# learning the error rates.
errR <- learnErrors(filtRs, multithread=TRUE)
# 01224846 total bases in 438036 reads from 15 samples will be used for
# learning the error rates.

# Graficando tasas de error
#plotErrors(errF, nominalQ=TRUE)

# Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

# Ensamble de reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construcción de tabla de secuencias
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Removiendo quimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))

write.csv(seqtab.nochim, file = "seqtab_nochim_16s_2024.csv")
head(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)
# 0.798244
# 0.584593 19/02/2025
# Inspect distribution of sequence lengths without quimeras
table(nchar(getSequences(seqtab.nochim)))

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

write.csv(track, file = "resumen_filtrado_16s_2024.csv")

## cuando lo hago en la compu
setwd("/home/luis/datos_16s_setC_may24")
seqtab.nochim = read.csv("seqtab_nochim_16s_2024.csv", row.names = 1)
seqtab.nochim2 =  as.matrix(seqtab.nochim)
head(seqtab.nochim2)

### Asignando taxonomía
silva <- "/home/.lvillanueva/work/analisis_2024/silva_nr99_v138.1_train_set.fa.gz"
#silva <- "silva_nr99_v138.1_train_set.fa.gz"

taxa <- assignTaxonomy(seqtab.nochim, silva, multithread=TRUE)

# extension
taxa2 <- addSpecies(taxa, "/home/.lvillanueva/work/analisis_2024/silva_species_assignment_v138.1.fa.gz")
#taxa2 <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz")

head(taxa)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
write.csv(taxa.print, file = "taxa_16s_2024.csv")

head(taxa2)
taxa2.print <- taxa2 # Removing sequence rownames for display only
rownames(taxa2.print)
head(taxa2.print)
write.csv(taxa2.print, file = "taxa_16s_species_2024.csv")

stopCluster(cl)
###############################################################################
