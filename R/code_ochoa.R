#ASV non chimeric sequences, taxonomic annotation and phylogenetic tree were 
# obtained followin: 
#Callahan BJ, Sankaran K, Fukuyama JA et al. 
#Bioconductor Workflow for Microbiome 
#Data Analysis: from raw reads to community analyses 
#[version 2; peer review: 3 approved]. F1000Research 2016,
# 5:1492 (https://doi.org/10.12688/f1000research.8986.2)

#######         Phyloseq raw object

#Magellan penguin meta data
site <- c(rep("back", 11), rep("foot", 10), rep("nest", 9), rep("chest", 10))
sp <- rep("S. magellanicus", 40)
regs <- data.frame(Regions = site,
                   Species = sp)
#King peguin meta data
kreg <- c(rep("back",8), rep("foot", 8), rep("chest", 8))
Species <- rep("A. patagonicus", 8)
king_meta <- data.frame(Regions = kreg,
                        Species)

### Merge metadata objects by rows

Meta <- rbind(king_meta, regs)

#Assing sample names to rownames of metadata data frame
row.names(Meta) <- row.names(pen)


###### Phyloseq object
#ASV table - non chimeric sequences from DADA2
#tax table - non chimeric sequences annotation file 
#Sample data - described above
#phy tree - phylogenetic tree constructed with non chimeric sequences

RAW <- phyloseq(otu_table(pen, taxa_are_rows = FALSE),
                tax_table(tax),
                sample_data(Meta),
                phy_tree(tree$tree))

# ADD DNA concentration data in DNA micrograms

DNA_quant_ug <- c(0.002,
                  0.004,
                  0.004,
                  0.004,
                  0.003,
                  0.002,
                  0.002,
                  0.009,
                  0.002,
                  0.002,
                  0.008,
                  0.002,
                  0.002,
                  0.002,
                  0.004,
                  0.006,
                  0.124,
                  0.004,
                  0.002,
                  0.027,
                  0.016,
                  0.005,
                  0.005,
                  0.059,
                  0.007,
                  0.002,
                  0.002,
                  0.005,
                  0.002,
                  0.002,
                  0.004,
                  0.004,
                  0.006,
                  0.004,
                  0.002,
                  0.009,
                  0.004,
                  0.006,
                  0.002,
                  0.003,
                  0.002,
                  0.001,
                  0.002,
                  0.003,
                  0.004,
                  0.029,
                  0.009,
                  0.035,
                  0.018,
                  0.025,
                  0.004,
                  0.004,
                  0.036,
                  0.017,
                  0.004,
                  0.006,
                  0.002,
                  0.013,
                  0.007,
                  0.002,
                  0.006,
                  0.005,
                  0.004,
                  0.02)

length(DNA_quant_ug)

# Add metadata field to RAW phyloseq object
sample_data(RAW)$DNA_quant_ug <- DNA_quant_ug

RAW

View(tax_table(RAW))

##### Filters

# Remove contaminants using the frecuency method

library(decontam)

contam_df_freq <- isContaminant(RAW, method = "frequency", conc = "DNA_quant_ug")

# Check number of ASVs identified as contaminants
table(contam_df_freq$contaminant)
# 20654 total ASVs
# 117 AsvS identified as contaminants 

# Remove this contaminant ASVs

RAW_nocontam <- prune_taxa(!contam_df_freq$contaminant, RAW)
sum(sample_sums(RAW)) - sum(sample_sums(RAW_nocontam))
#Remove taxa that only appears in one sample

#Compute prevalence
prevdf = apply(X = otu_table(RAW_nocontam),
               MARGIN = ifelse(taxa_are_rows(RAW_nocontam), yes = 1, no = 2),
               FUN = function(x) {sum(x>0)})

#Add taxonomy to prevalence data

prevdf = data.frame(prevalence = prevdf,
                    TotalAbundance = taxa_sums(RAW_nocontam),
                    tax_table(RAW_nocontam))


#Compute total and average prevalence for each phylum
plyr::ddply(prevdf, "Phylum", function(df1)
  {cbind(mean(df1$prevalence), sum(df1$prevalence))})


#Phyla whose prevalence is 1, that is, only occurs in one sample
filterPhyla = c("Zixibacteria", "WS2", "Cloacimonadota", "FCPU426")
#FIlter low prevalence phyla
filtpen <- subset_taxa(RAW_nocontam, !Phylum %in% filterPhyla)
sum(sample_sums(filtpen))
sum(sample_sums(RAW_nocontam))

#Remove taxa with 0 observation
penps <- prune_taxa(taxa_sums(filtpen) > 0, filtpen)

#Remove taxa annotated as Chloroplast or Mithocondria
penclean <- penps %>% subset_taxa(Order!="Chloroplast")
penclean2 <- penclean %>% subset_taxa(Family!= "Mitochondria")

sum(sample_sums(RAW_nocontam)) - sum(sample_sums(penclean))

sum(sample_sums(penclean)) - sum(sample_sums(penclean2))

#Final ASV data set: 17 227 ASVs
range(sample_sums(penclean2))
#Library size spans 30 104 - 94 631 reads

### Normalize library size

penrare <- rarefy_even_depth(penclean2, 30100, rngseed = 123, replace = FALSE)             

# phyloseq object with normalized library size has
# 17 005 ASVs 

range(sample_sums(penrare))


######          COMMUNITY ANALYSIS        #####

######### LOAD REQUIRED LIBRARIES

library(MicrobiotaProcess)
library(reltools)
library(microbiomeutilities)
library(microbiome)
library(vegan)
library(picante)
library(gridExtra)
library(dplyr)
library(tidyr)
library(reshape2)
library(UpSetR)
library(minpack.lm)
library(Hmisc)
library(stats4)
library(patchwork)
library(phyloseq)

#Phyloseq object
penrare

### Add metadata file (useful for Upset plot)

sampleData(penrare)$Penguin_niche <- c(rep("Ap back", 8), 
                                       rep("Ap foot", 8), 
                                       rep("Ap chest", 8),
                                       rep("Sm back", 11), 
                                       rep("Sm foot", 10 ),
                                       rep("Sm nest soil", 9),
                                       rep("Sm chest", 10))


#### HEATMAPs of most abundant bacteria

### Collapse dataset taxonomy to 10 most abundant families

topfam <- aggregate_top_taxa2(penrare, "Family", top = 10)
View(otu_table(topfam))
# transform families dataset to relative abundance

famrel <- transform(topfam, "compositional")
View(otu_table(famrel))
#### Subset by species

cf <- subset_samples(famrel, especies == "1CF")

mag <- subset_samples(famrel, Species == "S. magellanicus")

#### Subset by species from original phyloseq object

king_rare <- subset_samples(penrare, Species == "A. patagonicus")
mag_rare <- subset_samples(penrare, Species == "S. magellanicus")


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


#### TOP FIVE FAMILIES
library(microbiomeutilities)
#KING PENGUIN

kcols <- c("steelblue", "grey50", "red2")

pn <- plot_taxa_boxplot(king_rare,
                        taxonomic.level = "Family",
                        top.otu = 5, 
                        group = "Regions",
                        add.violin= FALSE,
                        title = "Top five families", 
                        keep.other = FALSE,
                        group.order = c("back","chest","foot"),
                        group.colors = kcols,
                        dot.size = 1)

print(pn + theme_biome_utils())

pn$data
View(pn$data)

#### MAGELLANIC PENGUIN

mcols <- c("steelblue", "grey50", "red2", "black")

magpn <- plot_taxa_boxplot(mag_rare,
                           taxonomic.level = "Family",
                           top.otu = 5, 
                           group = "Regions",
                           add.violin= FALSE,
                           title = "Top five families", 
                           keep.other = FALSE,
                           group.order = c("back","chest","foot", "nest"),
                           group.colors = mcols,
                           dot.size = 1)

print(magpn + theme_biome_utils())

View(magpn$data)


### Heatmap for 10 most abundant bacterial families in Magellan penguin body sites and nest soil

heatmap_ma <- simple_heatmap(mag,
                             group.facet = "Regions",
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
heatmap_ma

### Arrange heatmaps in one graphic

heatmaps <- grid.arrange(heatmap_k, heatmap_ma, nrow = 2)

### Save heatmaps

ggsave("heatmaps_june.pdf", 
       plot = heatmaps,
       width = 20, 
       height = 20, 
       units = "cm")


#######       CORE BACTERIA     #####

#### Get data with binary format to perfom upset graphic between all penguin sample types

uppen <- get_upset(penrare, factorNames = "Penguin_niche")

### Upset plot with all penguin sample types directly exported as a PDF


pdf(file="pensppsharedASV_july.pdf", 
    onefile=FALSE,
    width = 10) 

upset(uppen,
      sets = c("Ap back", "Ap chest", "Ap foot", 
               "Sm back", "Sm chest", "Sm foot", "Sm nest soil"),
      order.by = "freq",
      mainbar.y.label = "Number of ASVs shared",
      queries = list(list(query = intersects,
                          params = list("Ap back", "Sm back"),
                          color = "lightgoldenrod",
                          active = T),
                     list(query = intersects,
                          params = list("Ap foot", "Sm foot"),
                          color = "lightblue",
                          active = T),
                     list(query = intersects, 
                          params = list("Ap chest", "Sm chest"),
                          color = "indianred2",
                          active = T),
                     list(query = intersects, 
                          params = list("Ap back", 
                                        "Ap chest",
                                        "Ap foot"),
                          color = "orange",
                          active = T),
                     list(query = intersects, 
                          params = list("Sm back",
                                        "Sm chest",
                                        "Sm foot"),
                          color = "pink",
                          active = T),
                     list(query = intersects, 
                          params = list("Sm back",
                                        "Sm chest",
                                        "Sm foot",
                                        "Sm nest soil"),
                          color = "green",
                          active = T),
                     list(query = intersects, 
                          params = list("Ap back", 
                                        "Ap chest",
                                        "Ap foot",
                                        "Sm back",
                                        "Sm chest",
                                        "Sm foot",
                                        "Sm nest soil"),
                          color = "purple",
                          active = T),
                     list(query = intersects, 
                          params = list("Ap back", 
                                        "Ap chest",
                                        "Ap foot",
                                        "Sm back",
                                        "Sm chest",
                                        "Sm foot"),
                          color = "blue",
                          active = T)))
dev.off()

### CORE BACTERIA

###Transform original phyloseq object to relative abundance
penrel <- transform(penrare, "compositional")
### ADD "ASV" tag to each ASV
allpenASV <- add_refseq(penrel, tag = "ASV")
### Add best taxonomy level possible to each ASV
allpenASV <- format_to_besthit(allpenASV)

# Core tresholds: minimum ASV relative abundance 0.001 % minimum prevalence among samples 90%

#Among all samples
full_core <- core(allpenASV, detection = 0.001, prevalence = .9)
View(tax_table(full_core))

#Among king penguin samples

#Add ASV tag and best taxonomic hit to king penguin samples
allking
kingASV <- add_refseq(allking, tag = "ASV")
kingASV2 <- format_to_besthit(kingASV)
king_ASV2 <- transform(kingASV2, "compositional")
#Core bacteria in king penguin

king_core <- core(king_ASV2, detection = 0.001, prevalence = .9)
View(tax_table(king_core))

#Among Magellan penguin samples

allmag
magASV <- add_refseq(allmag, tag = "ASV")
magASV2 <- format_to_besthit(magASV)
mag_ASV2 <- transform(magASV2, "compositional")

mag_core <- core(mag_ASV2, detection = 0.001, prevalence = .9)

View(tax_table(mag_core))



#######               ALPHA DIVERSITY ANALYSIS       #####

#### Subset by species from original phyloseq object
penrare

allking <- subset_samples(penrare, Species == "A. patagonicus")
allmag <- subset_samples(penrare, Species == "S. magellanicus")

###                          KING PENGUIN

allking

#DATA FRAME with alpha diversity (PD and Shannon) across body sites
div <- data.frame(
  "Shannon" = estimate_richness(allking, measures = "Shannon"),
  "PD" = pd(samp = otu_table(allking), tree = phy_tree(allking), include.root = FALSE),
  "Region" = sample_data(allking)$Region,
  "Species" = sample_data(allking)$Species)

colnames(div)

# Change column names

colnames(div) <- c("Shannon", "PD", "PD.SR", "Region", "Species")

### King penguin alpha diversty (Shannon & PD) trends across body sites

king_alpha <- div %>%
  gather(key = metric, value = value, c("Shannon", "PD")) %>%
  mutate(metric = factor(metric, levels = c("Shannon", "PD"))) %>%
  ggplot(aes(x = Region, y = value)) +
  geom_jitter(aes(colour = Region), size = 3, alpha = 1) +
  scale_color_manual(values = c("lightgoldenrod", "indianred2", "lightblue", "tan4")) +
  theme_classic() +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme(legend.position="none") +
  labs(title = "King penguin alpha diversity measures")

king_alpha

#Exportamos grafico de alfa div a pdf

ggsave("king_alpha.pdf", 
       plot = king_alpha,
       width = 25, 
       height = 20, 
       units = "cm")

#### King penguin alpha diversity statistics

### Kruskal test between PD across king penguin body sites
kruskal.test(PD~Region, data=div)

#### Pairwise wilcoxon test to detect specific differences
pairwise.wilcox.test(div$PD, div$Region, exact=F)

### Kruskal test between Shannon diversity across king penguin body sites
kruskal.test(Shannon~Region, data=div) # NO DIFFERENCES

###             MAGELLAN PENGUIN

allmag

#DATA FRAME with alpha diversity (PD and Shannon) across body sites

magdiv <- data.frame(
  "Shannon" = estimate_richness(allmag, measures = "Shannon"),
  "PD" = pd(samp = otu_table(allmag), tree = phy_tree(allmag), include.root = FALSE),
  "Region" = sample_data(allmag)$Region,
  "Species" = sample_data(allmag)$Species)

### Change column names
colnames(magdiv) <- c("Shannon", "PD", "PD.SR", "Region", "Species")

### Magellan penguin alpha diversty (Shannon & PD) trends across body sites and nest

mag_alpha <- magdiv %>%
  gather(key = metric, value = value, c("Shannon", "PD")) %>%
  mutate(metric = factor(metric, levels = c("Shannon", "PD"))) %>%
  ggplot(aes(x = Region, y = value)) +
  geom_jitter(aes(colour = Region), size = 3, alpha = 1) +
  scale_color_manual(values = c("lightgoldenrod", "indianred2", "lightblue", "tan4")) +
  theme_classic()+
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme(legend.position="none") +
  labs(title = "Magellanic penguin alpha diversity measures")

mag_alpha

#### Magellan penguin alpha diversity statistics

### Kruskal test between PD across Magellan penguin body sites and nest soil
kruskal.test(PD~Region, data=magdiv)

#### Pairwise wilcoxon test to detect specific differences
pairwise.wilcox.test(magdiv$PD, magdiv$Region, exact=F)

### Kruskal test between Shannon across Magellan penguin body sites and nest soil
kruskal.test(Shannon~Region, data=magdiv) # NO DIFFERENCES


### Combine alpha diversity plots from each species and export to pdf

alpha <- king_alpha/mag_alpha

ggsave("penguin_alpha_july26.pdf", 
       plot = alpha,
       width = 25, 
       height = 20, 
       units = "cm")

##############      BETA DIVERSITY ANALYSIS     ##########

###   King penguin

# Transform king penguin phyloseq to relative abundances and compute weighted unifrac distance

allkrel <- transform(allking, "compositional")
wun_dist <- distance(allkrel, method = "wunifrac")

#Test dispersion among body sites
disp_wun <- betadisper(wun_dist, sample_data(allking)$Regions)
permutest(disp_wun)

#PCoA Ordination

kwuni <- ordinate(allkrel, method = "PCoA", distance = "wunifrac")

kingord <- plot_ordination(
  allkrel,
  kwuni, type="samples", color= "Regions") + 
  scale_colour_manual(values = c("lightgoldenrod", "indianred2", "lightblue")) +
  geom_point(size = 4) +
  stat_ellipse(aes(group = Regions), linetype = 2) +
  theme_classic() +                                                      #changes theme, removes grey background
  ggtitle("King penguin body sites")

### PERMANOVA global test among King penguin body sites

permall <- adonis(wun_dist ~ sample_data(allkrel)$Regions)

#Anova-like table for PERMANOVA results
permall$aov.tab

### Paired PERMANOVA for King penguin body site comparisons 


#Paired permanova function

pairwise.adonis.physeq <- function(x,factors,dist.m="wunifrac", perm=999,p.adjust.m ='BH'){
  co <- combn(unique(as.character(sample_data(x)[,factors][[1]])),2)
  pairs <- c()
  F.Model <- c()
  R2 <- c()
  p.value <- c()
  for (elem in 1:ncol(co)){ 
    cosa<-c(co[1,elem],co[2,elem])
    filtered<-prune_samples(sample_names(x)[which( sample_data(x)[,factors][[1]] %in% cosa)],x)
    filtered<- filter_taxa(filtered, function(x) sum(x) > 0, TRUE)
    x1=distance(filtered,dist.m)
    form<-paste("x1","~",factors)
    ad = adonis(as.formula(form),as(sample_data(filtered),"data.frame"),perm=perm);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
  return(pairw.res)
}

#Perform PERMANOVA paired comparisons

allkrel %>% filter_taxa(function(x) sum(x) > 0, T) %>% 
  pairwise.adonis.physeq(.,factors = "Regions", dist.m = "wunifrac")

###       SIMPER analysis

#We use phyloseq object with best taxonomy hit and ASV tags
kingASV2
kv2 <- transform(kingASV2, "compositional")

#Save as data frame sample data and ASV table
mdk <- data.frame(sample_data(kingASV2))
otuk2 <- data.frame(otu_table(kv2))

#Run simper analysis and save it
simking <- simper(otuk2, mdk$Regions, permutations = 999)
sumsk <- summary(simking)

#Access particular combinations and store them in csv format for manipulation in Excel

###back vs foot
back_foot <- sumsk$back_foot
head(back_foot)
write.csv(back, "back_foot_june.csv")
####back vs chest

back_chest <- sumsk$back_chest
head(back_chest)
write.csv(chest_back, "back_chest_june.csv")

####foot vs chest

chest_foot <- sumsk$foot_chest
head(chest_foot)
write.csv(back, "foot_chest_june.csv")

###       Magellanic penguin

# Transform Magellan penguin phyloseq to relative abundances and compute weighted unifrac distance

allmrel <- transform(allmag, "compositional")
mwundist <- distance(allmag, method = "wunifrac")

#Test dispersion among body sites
magdisp_wun <- betadisper(mwundist, sample_data(allmag)$Regions)
permutest(magdisp_wun)

#PCoA ordination

magwuni <- ordinate(allmrel, method = "PCoA", distance = "wunifrac")

magord <- plot_ordination(
  allmrel,
  magwuni, type="samples", color= "Regions") + 
  scale_colour_manual(values = c("lightgoldenrod", "indianred2", "lightblue","tan4")) +
  geom_point(size = 4) +
  stat_ellipse(aes(group = Regions), linetype = 2) +
  theme_classic() +                                                      #changes theme, removes grey background
  ggtitle("Magellanic penguin body sites and nest")

### PERMANOVA global test among Magellanic penguin body sites

per_mag <- adonis(mwundist ~ sample_data(allmrel)$Regions)

#Anova-like table for PERMANOVA results
per_mag$aov.tab

### Paired PERMANOVA for Magellanic penguin body site comparisons 

allmrel %>% filter_taxa(function(x) sum(x) > 0, T) %>% 
  pairwise.adonis.physeq(.,factors = "Regions", dist.m = "wunifrac")

###   SIMPER Analysis

#We use phyloseq object with best taxonomy hit and ASV tags

magASV2
mv2 <- transform(magASV2, "compositional")

#Save as data frame sample data and ASV table
magsam <- data.frame(sample_data(magASV2))
otum <- data.frame(otu_table(mv2))

#Run simper analysis and save it
simmag <- simper(otum, magsam$Regions, permutations = 999)
sumsm <- summary(simmag)

#Access particular combinations that render significant results in paired permanova comparisons and store them in csv format for manipulation in Excel

###back vs nest

back <- sumsm$back_nest
head(back)
write.csv(back, "back_nest_june.csv")

#chest vs nest
chest <- sumsm$nest_chest
head(chest)
write.csv(chest, "chest_nest_june.csv")

#foot vs nest
pata <- sumsm$foot_nest
head(pata)
write.csv(chest, "foot_nest_june.csv")

#### Intraspecific ordination array and export as PDF

pens <- kingord + magord
ggsave("penord2_july26.pdf",
       plot = pens,
       width = 30, 
       height = 20, 
       units = "cm")

####          Global ordination and interspecific beta diversity comparisons

#Normalized phyloseq object
penrare
#Remove nest, to display only body sites comparisons
bodies <- subset_samples(penrare, Penguin_niche != "Sm nest soil")

# Transform penguin spp phyloseq to relative abundances and compute weighted unifrac distance

bodrel <- transform(bodies, "compositional")
bodwund <- distance(bodies, method = "wunifrac")

#Test dispersions among sample types
boddisp_wun <- betadisper(bodwund, sample_data(bodies)$Species)
permutest(boddisp_wun) #Dispersion differs between species

#PCoA Species ordination

bodwuni <- ordinate(bodrel, method = "PCoA", distance = "wunifrac")

bodord <- plot_ordination(
  bodrel,
  bodwuni, type="samples", color= "Species") + 
  scale_colour_manual(values = c("lightgoldenrod", "indianred2")) +
  geom_point(size = 3) +
  stat_ellipse(linetype = 2) +
  theme_classic() +                                                     
  ggtitle("Penguin species comparison")

### PERMANOVA global test among Magellan and King penguin body sites

perm_species <- adonis(bodwund ~ sample_data(bodrel)$Species)

#Anova-like table for PERMANOVA results

perm_species$aov.tab 

############# Interspecific body site ordinations

##################    BACK

#Subset backs from both species
back <- subset_samples(bodrel, Regions == "back")

#Interspecific back ordination
backwuni <- ordinate(back, method = "PCoA", distance = "wunifrac")


backord <- plot_ordination(back, backwuni, type="samples", color= "Species") + 
  scale_colour_manual(values = c("lightgoldenrod", "indianred2")) +
  geom_point(size = 3) +
  stat_ellipse(linetype = 2) +
  theme_classic() +                                                      
  ggtitle("Back comparison")

#Test dispersions among penguin backs

backw <- distance(back, method = "wunifrac")
backdisp <- betadisper(backw, sample_data(back)$Species)
permutest(backdisp) #Dispersions differ among penguin backs

#PERMANOVA test between penguins backs

permback <- adonis(backw ~ sample_data(back)$Species)
permback$aov.tab 


##################    CHEST


#Subset chests from both species
chest <- subset_samples(bodrel, Regions == "chest")

#Interspecific chest ordination

chestwuni <- ordinate(chest, method = "PCoA", distance = "wunifrac")

chestord <- plot_ordination(
  chest,
  chestwuni, type="samples", color= "Species") + 
  scale_colour_manual(values = c("lightgoldenrod", "indianred2")) +
  geom_point(size = 3) +
  stat_ellipse(linetype = 2) +
  theme_classic() +                     
  ggtitle("Back comparison")

#Test dispersions among penguin chests
chestw <- distance(chest, method = "wunifrac")
chestdisp <- betadisper(chestw, sample_data(chest)$Species)
permutest(chestdisp) #Dispersions differ among penguin chests

#PERMANOVA test between penguins chests
permchest <- adonis(chestw ~ sample_data(chest)$Species)
permchest$aov.tab 

##################    FOOT

#Subset foot from both species
foot <- subset_samples(bodrel, Regions == "foot")

#Interspecific foot ordination

footwuni <- ordinate(foot, method = "PCoA", distance = "wunifrac")

footord <- plot_ordination(
  foot,
  footwuni, type="samples", color= "Species") + 
  scale_colour_manual(values = c("lightgoldenrod", "indianred2")) +
  geom_point(size = 4) +
  stat_ellipse(linetype = 2) +
  theme_classic() +                                                      #changes theme, removes grey background
  ggtitle("Foot comparison")

#Test dispersions among penguin feet
footw <- distance(foot, method = "wunifrac")
dispfo <- betadisper(footw, sample_data(foot)$Species)
permutest(dispfo) #Dispersions differ among penguin feet

#PERMANOVA test between penguins feet
permfoot <- adonis(footw ~ sample_data(foot)$Species)
permfoot$aov.tab 

### Arrange interspecific ordinations and export it as PDF

bodiesords <- bodord + backord/chestord/footord

ggsave("interpen2_june.pdf",
       plot = bodiesords,
       width = 39, 
       height = 20, 
       units = "cm")


#######         SLOAN NEUTRAL COMMUNITY MODEL          #######

#Code from:
#Burns, A., Stephens, W., Stagaman, K. et al. 
#Contribution of neutral processes to the assembly of 
#gut microbial communities in the zebrafish over host 
#development. ISME J 10, 655-664 (2016). 
#https://doi.org/10.1038/ismej.2015.142

# Normalized ASV table (here exemplified with king penguin)
# These steps must be repeated with Magellanic penguin (without nest samples)

allking

onlymag <- subset_samples(bodies, Species == "S. magellanicus")

kspp <- otu_table(allking)

# FIT SLOAN NEUTRAL COMMUNITY MODEL

fit_sncm <- function(spp, pool=NULL, taxon=NULL){
  
  options(warn=-1)
  
  #Calculate the number of individuals per community
  N <- mean(apply(spp, 1, sum))
  
  #Calculate the average relative abundance of each taxa across communities
  if(is.null(pool)){
    p.m <- apply(spp, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  } else {
    p.m <- apply(pool, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  }
  
  #Calculate the occurrence frequency of each taxa across communities
  spp.bi <- 1*(spp>0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]
  
  #Combine
  C <- merge(p, freq, by=0)
  C <- C[order(C[,2]),]
  C <- as.data.frame(C)
  #Removes rows with any zero (absent in either source pool or local communities)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
  p <- C.0[,2]
  freq <- C.0[,3]
  names(p) <- C.0[,1]
  names(freq) <- C.0[,1]
  
  #Calculate the limit of detection
  d = 1/N
  
  ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
  m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.001))
  m.ci <- confint(m.fit, 'm', level=0.95)
  
  ##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
  freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
  Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
  
  pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Calculate AIC for Poisson model
  pois.LL <- function(mu, sigma){
    R = freq - ppois(d, N*p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  
  aic.pois <- AIC(pois.mle, k=2)
  bic.pois <- BIC(pois.mle)
  
  ##Goodness of fit for Poisson model
  pois.pred <- ppois(d, N*p, lower.tail=FALSE)
  Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))
  
  pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Results
  fitstats <- data.frame(
    m=as.numeric(coef(m.fit)),
    m.ci=as.numeric(coef(m.fit)-m.ci[1]),
    poisLL=as.numeric(pois.mle@details$value),
    Rsqr=as.numeric(Rsqr), # measuring fit, #comparing fit differing datasets to the same model
    Rsqr.pois=as.numeric(Rsqr.pois),
    RMSE=as.numeric(RMSE), # measuring fit #comparing fit differing datasets to the same model
    RMSE.pois=as.numeric(RMSE.pois),
    AIC.pois=as.numeric(aic.pois),  #comparing differing models to the dataset
    BIC.pois=as.numeric(bic.pois), #comparing differing models to the dataset
    N=as.numeric(N),
    Samples=as.numeric(nrow(spp)),
    Richness=as.numeric(length(p)),
    Detect=as.numeric(d))
  
  A <- cbind(p, freq, freq.pred, pred.ci[,2:3])
  A <- as.data.frame(A)
  colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr')
  if(is.null(taxon)){
    B <- A[order(A[,1]),]
  } else {
    B <- merge(A, taxon, by=0, all=TRUE)
    row.names(B) <- B[,1]
    B <- B[,-1]
    B <- B[order(B[,1]),]
  }
  B <- B[!is.na(B$freq),]
  # fit_class for graphing
  B$fit_class <-"As predicted"
  B[which(B$freq < B$pred.lwr),"fit_class"]<- "Below prediction"
  B[which(B$freq > B$pred.upr),"fit_class"]<- "Above prediction"
  B[which(is.na(B$freq)),"fit_class"]<- "NA"
  
  # combine fit stats and predicitons into list
  i <- list(fitstats, B)
  names(i) <- c("fitstats", "predictions")
  return(i)
}

#Run neutral community model

unloadNamespace("MicrobiotaProcess") #Unload microbiota process, since has a redundant function for extracting taxonomy

kspp.out <- fit_sncm(kspp, pool=NULL, taxon=data.frame(tax_table(allking)))
kspp.out$fitstats

### Plot results

#Code extracted from https://rdrr.io/github/DanielSprockett/reltools/man/plot_sncm_fit.html

plot_sncm_fit <- function(spp.out, fill = NULL, title = NULL){
  
  tax_levels <- colnames(spp.out$predictions)[7:length(colnames(spp.out$predictions))-1]
  
  if(is.null(fill)){
    fill <- "fit_class"
  }
  
  r2_val <- paste("r^2 ==", round(spp.out$fitstats$Rsqr,4))
  m_val <- paste("m ==", round(spp.out$fitstats$m,4))
  df <- data.frame(t(table(spp.out$predictions$fit_class)))
  df <- df[,c(2,3)]
  colnames(df) <- c("Prediction", "ASV Abundance")
  
  p <- ggplot(data=spp.out$predictions)
  
  if(fill == "fit_class"){
    p <- p + geom_point(aes(x = log(p), y = freq, fill=eval(parse(text=fill))), shape =21, color="black", size =2, alpha=0.75)
    p <- p + scale_fill_manual(
      name = "Prediction",
      values = c("Above prediction" = "cyan4", "As predicted" = "blue", "Below prediction" = "goldenrod", "NA" = "white"),
      breaks = c("Above prediction", "As predicted", "Below prediction", "NA"),
      labels = c(paste0("Above prediction (",round((df[1,2]/spp.out$fitstats$Richness)*100, 1),"%)"),
                 paste0("As predicted (",round((df[2,2]/spp.out$fitstats$Richness)*100, 1),"%)"),
                 paste0("Below Prediction (",round((df[3,2]/spp.out$fitstats$Richness)*100, 1),"%)"),
                 paste0("NA (",df[4,2],")")))
    
  }else if (fill %in% tax_levels){
    p <- p + geom_point(aes(x = log(p), y = freq, fill=eval(parse(text=fill))), shape =21, color="black", size =2, alpha=0.75)
    p <- p + scale_fill_discrete(name = "Taxon")
    
  } else{
    print(paste0("fill variable: ", fill, " is not a valid taxonomic level or fit_class"))
  }
  
  p <- p + geom_line(aes(x = log(p), y = freq.pred), color = "blue")
  p <- p + geom_line(aes(x = log(p), y = pred.lwr), color = "blue", linetype="dashed")
  p <- p + geom_line(aes(x = log(p), y = pred.upr), color = "blue", linetype="dashed")
  p <- p + xlab("log(Mean Relative Abundance)")
  p <- p + ylab("Frequency")
  p <- p + ggtitle(title)
  p <- p + theme_classic()
  p <- p + annotate("text", x=mean(log(spp.out$predictions$p), na.rm = TRUE), y=0.95, size=5, label = r2_val, parse=TRUE)
  p <- p + annotate("text", x=mean(log(spp.out$predictions$p), na.rm = TRUE), y=0.9, size=5, label = m_val, parse=TRUE)
  
  return(p)
}

?geom_jitter
# Plot 

neutralking <- plot_sncm_fit(kspp.out, title = "Neutral model in King penguin body sites")

neutralking

### Extract data to obtain the taxa with most ASV in above and below prediction categories

res_king <- neutralking$data
write.csv(res_king, file = "neutral_king_june.csv")

# King penguin
# 
#Load taxa with most ASVs in each category

AB <- data_frame(Taxa = c("Psychrobacter",
                          "Chryseobacterium",
                          "Flavobacterium",
                          "Nocardioides",
                          "Deinococcus"),
                 ASV.number = c(10,
                                7,
                                4,
                                4,
                                3))

#Barplot of 5 taxa with most ASVs above prediction

above_k <-  ggplot(AB, aes(x = reorder(Taxa, ASV.number), y = ASV.number)) +
  geom_bar(stat = 'identity', fill = "cyan4") +
  xlab("Genus") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "italic")) +
  ggtitle("Genera with more ASVs above predicted frequency")


#Barplot of 5 taxa with most ASVs below prediction


BE <- data.frame(Taxa = c("Psychrobacter",
                          "Ornithobacterium",
                          "Fusobacterium",
                          "Suttonella",
                          "Corynebacterium"),
                 ASV.number = c(26,
                                11,
                                8,
                                8,
                                7))

#Barplot of 5 taxa with most ASVs below prediction

below_k <-  ggplot(BE, aes(x = reorder(Taxa, ASV.number), y = ASV.number)) +
  geom_bar(stat = 'identity', fill = "goldenrod") +
  xlab("Genus") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "italic")) +
  ggtitle("Genera with more ASVs below predicted frequency")




### Merge neutral model and barplots and export it to pdf

k_neutral_grobs <- neutralking + above_k/below_k

###EXPORTAMOS PDF

ggsave("neutral_king_july26.pdf",
       plot = k_neutral_grobs,
       width = 39, 
       height = 20, 
       units = "cm")

# Magellanic penguin
# 
#Load taxa with most ASVs in each category

ABmag <- data_frame(Taxa = c("Nocardioides",
                             "Sporosarcina",
                             "Pedobacter",
                             "Rhodococcus",
                             "Sphingomonas"),
                    ASV.number = c(11,
                                   6,
                                   4,
                                   4,
                                   4))

#Barplot of 5 taxa with most ASVs above prediction

above_mag <-  ggplot(ABmag, aes(x = reorder(Taxa, ASV.number), y = ASV.number)) +
  geom_bar(stat = 'identity', fill = "cyan4") +
  xlab("Genus") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "italic")) +
  ggtitle("Genera with more ASVs above predicted frequency")

above_mag

#Barplot of 5 taxa with most ASVs below prediction


BE_mag <- data.frame(Taxa = c("Psychrobacter",
                              "Fusobacterium",
                              "Bacteroides",
                              "Ornithobacterium",
                              "Staphylococcus"),
                     ASV.number = c(28,
                                    10,
                                    7,
                                    7,
                                    7))

#Barplot of 5 taxa with most ASVs below prediction

below_m <-  ggplot(BE_mag, aes(x = reorder(Taxa, ASV.number), y = ASV.number)) +
  geom_bar(stat = 'identity', fill = "goldenrod") +
  xlab("Genus") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "italic")) +
  ggtitle("Genera with more ASVs below predicted frequency")




### Merge neutral model and barplots and export it to pdf
### neutralking could have either data from king penguin or magellanic penguin!

k_neutral_grobs <- neutralking + above_k/below_k

mag_neutral_grobs <- neutralking + above_mag/below_m

###EXPORTAMOS PDF

ggsave("neutral_mag_july26.pdf",
       plot = mag_neutral_grobs,
       width = 39, 
       height = 20, 
       units = "cm")


##############      Supplemental figures/data     #########

### Rarefaction plots

cleanking <- subset_samples(penclean2, Species == "A. patagonicus")
cleanmag <- subset_samples(penclean2, Species == "S. magellanicus")


#Set seed
set.seed(1)
subsamples <- seq(0, 50000, by=1000)[-1]

#Repeat for each data set, here exemplified with king penguin clean data set
sample_data(cleanking)
p <- plot_alpha_rcurve(cleanmag, index="observed", 
                       subsamples = subsamples,
                       lower.conf = 0.025, 
                       upper.conf = 0.975,
                       group="Regions",
                       label.color = "brown3",
                       label.size = 3,
                       label.min = TRUE) 


### Most abundant genera

# In King penguin body site samples
dom_k <- dominant_taxa(allking,level = "Genus", group="Regions")
write.csv(dom_k$dominant_overview, "dom_genera_king_june.csv")

### In Magellan penguin body site and nest samples
dom_m <- dominant_taxa(allmag,level = "Genus", group="Regions")
write.csv(dom_m$dominant_overview, "dom_genera_mag_june.csv")


##### Psychrobacter phylogenetic tree sequences

allking <- subset_samples(penrare, Species == "A. patagonicus")
allmag <- subset_samples(penrare, Species == "S. magellanicus")

### Extract Psychrobacter sequences that only occur in King penguin
k_psy <- subset_taxa(allking, Genus == "Psychrobacter")

### Prune sequences with 0 abundance
k_psy <- prune_taxa(taxa_sums(k_psy) > 0, k_psy)

### TREE
k_tree <- plot_tree(k_psy,
                    color="Regions") +
  scale_color_manual(values = c("lightgoldenrod", "indianred2", "lightblue")) +
  ggtitle("Psychrobacter sequences associated to King penguin body sites")

### Magellan penguin phylogenetic tree

### Extract Psychrobacter sequences that only occur in King penguin
m_psy <- subset_taxa(allmag, Genus == "Psychrobacter")

### Prune sequences with 0 abundance
m_psy <- prune_taxa(taxa_sums(m_psy) > 0, m_psy)

### TREE
magtree <- plot_tree(m_psy,
                     color="Regions") +
  scale_color_manual(values = c("lightgoldenrod", "indianred2", "lightblue", "tan4")) +
  ggtitle("Psychrobacter sequences associated to Magellanic penguin samples")

### Arrange phylogenetic trees in one graphic and export it to PDF
trees <- k_tree + magtree

ggsave("psy_penspp_tree_july26.pdf", 
       plot = trees,
       width = 30, 
       height = 20, 
       units = "cm")

#####     Interspecific diversity comparisons

######## Shannon between penguin body sites

##Get specific columns from previous alpha diversity data frames
div[,c(1,4,5)]
magdiv[,c(1,4,5)]

#Shannon data frame 
Shdf <- rbind(div[,c(1,4,5)], magdiv[,c(1,4,5)])

### Remove nest soil samples to have only body site alpha diversity
Sh_bird <- subset(Shdf, Region != "nest")

#Prepare data to compare shannon diversity between body sites
Shdf <- data.frame(melt(Sh_bird))

cols <- c("lightgoldenrod", "indianred2")

SHintbox <- ggplot(data=Shdf) +
  geom_boxplot(aes(x=Region, y=value, fill=Species)) +
  theme_classic() + scale_fill_manual(values=cols) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13)) +
  theme(axis.text.y = element_text(hjust = 1, size=13)) +
  ggtitle("Shannon diversity") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=18))


#Specific tests comparing Shannon values between body sites
#foot
SH_foot <- subset(Sh_bird, Region == "foot")
kruskal.test(Shannon~Species, data=SH_foot)

#chest
SH_chest <- subset(Sh_bird, Region == "chest")
kruskal.test(Shannon~Species, data=SH_chest)
#back
SH_back <- subset(Sh_bird, Region == "back")
kruskal.test(Shannon~Species, data=SH_back)

########   PD between penguin body sites

##Get specific columns from previous alpha diversity data frames

div[,c(2,4,5)]
magdiv[,c(2,4,5)]

#PD data frame
PDdf <- rbind(div[,c(2,4,5)], magdiv[,c(2,4,5)])

### Remove nest soil samples to have only body site alpha diversity
PD_bird <- subset(PDdf, Region != "nest")

#Prepare data to compare PD diversity between body sites
df <- data.frame(melt(PD_bird))

cols <- c("lightgoldenrod", "indianred2")

PDintbox <- ggplot(data=df) +
  geom_boxplot(aes(x=Region, y=value, fill=Species)) +
  theme_classic() + scale_fill_manual(values=cols) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13)) +
  theme(axis.text.y = element_text(hjust = 1, size=13)) +
  ggtitle("PD diversity") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=18))

### PD diversity statistics between penguin samples

#Kruskal test comparing PD values between penguin species
kruskal.test(PD~Species, data=PD_bird)

#Specific tests comparing PD values between body sites
#foot
PD_foot <- subset(PD_bird, Region == "foot")
kruskal.test(PD~Species, data=PD_foot)
#Chest
PD_chest <- subset(PD_bird, Region == "chest")
kruskal.test(PD~Species, data=PD_chest)
#Back
PD_back <- subset(PD_bird, Region == "back")
kruskal.test(PD~Species, data=PD_back)


#Arrange plots and export to pdf
inter <- grid.arrange(PDintbox, SHintbox, ncol = 2)

ggsave("inter_alpha.pdf",
       plot = inter,
       width = 30, 
       height = 20, 
       units = "cm")