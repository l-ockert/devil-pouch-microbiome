# 1.0 Set-Up
## 1.1 Settings
```{r}
rm(list=ls()) # clears environment
theme_bw()
```

## 1.2 Load Packages
```{r}
library(qiime2R)
library(phyloseq)
library(tidyverse)
library(microbiome)
library(ggplot2)
library(decontam)
library(viridis)
library(cowplot)
library(ggpubr)
library(picante)
library(vegan)
```

# 2.0 Import Data
```{r}
#create base phyloseq object
##read in unfiltered table and taxonomy files from qiime2 output
pseq <- qza_to_phyloseq(features = "import-files/unfiltered-table.qza", taxonomy = "import-files/taxonomy.qza")

#read in metadata
dev.metadata <- read.table("import-files/metadata.tsv", header = T, sep = "\t", row.names = 1) #double check data was imported correctly
pseq <- merge_phyloseq(pseq, sample_data(dev.metadata))
rm(dev.metadata)

#read in tree
tree.import <- read_tree("import-files/rootedtree.nwk") #rooted tree needed for UniFrac
pseq <-merge_phyloseq(pseq, phy_tree(tree.import))
rm(tree.import)
```

## 2.1 Clean Data
```{r}
#rename taxa names to ASV1, ASV2, ASV3, etc.
taxa_names(pseq)
n_seqs <- seq(ntaxa(pseq))
len_n_seqs <- nchar(max(n_seqs))
taxa_names(pseq) <- paste("ASV", formatC(n_seqs, 
                                            width = len_n_seqs, 
                                            flag = "0"), sep = "_")
taxa_names(pseq) #check changes were made

#clean taxonomy
taxonomy_table <- data.frame(tax_table(pseq))
taxonomy_table <- taxonomy_table %>%
  rownames_to_column(var = "OTU") %>%
  mutate(Kingdom = gsub("d__", "", Kingdom),
         Phylum = gsub("p__", "", Phylum),
         Phylum = gsub("-", "_", Phylum),
         Phylum = gsub("_clade(Marine_group_B)", "", Phylum, fixed = T),
         Class = gsub("c__", "", Class),
         Class = gsub("_clade(Marine_group_B)", "", Class, fixed = T),
         Class = gsub("-", "_", Class),
         Order = gsub("o__", "", Order),
         Order = gsub("_(SR1)", "", Order, fixed = T),
         Order = gsub("_clade(Marine_group_B)", "", Order, fixed = T),
         Order = gsub("-", "_", Order),
         Family = gsub("f__", "", Family),
         Family = gsub("_(SR1)", "", Family, fixed = T),
         Family = gsub("_clade(Marine_group_B)", "", Family, fixed = T),
         Family = gsub("_(Subgroup_1)", "", Family, fixed = T),
         Family = gsub("-", "_", Family),
         Genus = gsub("g__", "", Genus),
         Genus = gsub("_(SR1)", "", Genus, fixed = T),
         Genus = gsub("_clade(Marine_group_B)", "", Genus, fixed = T),
         Genus = gsub("-", "_", Genus),
         Species = gsub("s__", "", Species),
         Species = gsub("-", "_", Species))

empty_as_na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
}

taxonomy_table <- taxonomy_table %>% mutate_each(list(empty_as_na)) 
taxonomy_table <- data.frame(taxonomy_table, row.names = 1)
taxonomy_table <- as.matrix(taxonomy_table)

tax_table(pseq) <- NULL
pseq <- merge_phyloseq(pseq, tax_table(taxonomy_table))

summarize_phyloseq(pseq)
```

# 3.0 Decontaminate and Filter Data
## 3.1 View Data
```{r}
#visualise library size
libsize <- as.data.frame(sample_data(pseq))
libsize$LibrarySize <- sample_sums(pseq)
libsize <- libsize [order(libsize$LibrarySize),]
libsize$Index <- seq(nrow(libsize))

#dictate order of samples to maintain consistency with colour theme
libsize$sample.or.control <- factor(libsize$sample.or.control, levels = c("Negative Control", "Biological Sample"))

#plot
ggplot(data = libsize, aes(x=Index, y=LibrarySize, color=sample.or.control)) + geom_point() + ylab("Library Size (Reads per Sample)") + xlab("Sample Index") + labs(colour = 'Sample Type') + theme_bw()

#save
ggsave("./Figures/decontam-libsize.png", height = 6, width = 10)
```

## 3.2 Ordination of Negative Control Sequences
```{r}
#ordinate and plot
pseq.ord <- ordinate(pseq, "PCoA", "bray")
ord <- plot_ordination(pseq, pseq.ord, type = "sample.or.control", color = "sample.or.control") + labs(colour='Sample Type') + theme_bw()

# Set colours
ord + scale_color_manual(values = c("#01BFC4", "#F8766D"))

# Save
ggsave("./Figures/decontam-PCoA.png", height = 5, width = 10)
```

## 3.3 Plot Decontam Scores
```{r}
sample_data(pseq)$is.neg <- sample_data(pseq)$sample.or.control == "Negative Control"
contamfull.prev.05 <- isContaminant(pseq, method="prevalence", neg="is.neg", threshold = 0.5) #threshold of 0.5 removes small samples without interfering with valuable data but this threshold is different for every dataset
table(contamfull.prev.05$contaminant) 

#plot
qplot(contamfull.prev.05$p, geom = "histogram", bins = 100, ylab = "Number of Species", xlab = "Decontam Score")

#save
ggsave("./Figures/decontam-score.png", height = 5, width = 10)
```

## 3.4 Visualise Prevalence
```{r}
#set variables
devil.contam <- transform_sample_counts(pseq, function(abund) 1*(abund>0)) #transfporm to relative abundance
devil.contam.neg <- prune_samples(sample_data(devil.contam)$sample.or.control == "Negative Control", devil.contam) #select negative controls
devil.contam.pos <- prune_samples(sample_data(devil.contam)$sample.or.control == "Biological Sample", devil.contam) #select actual samples

#make data.frame of prevalence in positive and negative samples
df.devil.contam <- data.frame(pa.pos=taxa_sums(devil.contam.pos),
                              pa.neg=taxa_sums(devil.contam.neg),
                      contaminant=contamfull.prev.05$contaminant)
ggplot(data=df.devil.contam, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (Biological Samples)") + labs(colour = "Contaminating ASV") + scale_color_manual(labels = c("False", "True"), values = c("#01BFC4", "#F8766D")) + theme_bw()

ggsave("./Figures/decontam-bar.png", height = 5, width = 10)
```

## 3.5 Remove Contaminants
```{r}
#you can remove contaminants as a group or individually per sample. Here we will remove them as a group as the same regents were used and so differences between samples are negligable. 
pseq.nocontam <- prune_taxa(!contamfull.prev.05$contaminant, pseq)
pseq.no.contam.rm <- prune_samples(sample_data(pseq.nocontam)$sample.or.control == "Biological Sample", pseq.nocontam)
pseq.clean <- pseq.no.contam.rm

#compare
pseq
pseq.clean #when comparing the two phyloseq objects, you can see the reduction in ASVs (reflected in otu_table, tax_table, and phy_tree), as well as the reduction in number of samples (7 samples removed here).

summarize_phyloseq(pseq.clean)

contaminating.taxa <- pseq.clean@tax_table[(which(contamfull.prev.05$contaminant))]
write.csv(contaminating.taxa, "contaminating-taxa.txt", row.names = T) #table of taxa identified as contaminating
```

## 3.6 Remove Mitochondrial, Chloroplast, Archael, and Kingdom Unassigned ASVs
```{r}
pseq <- pseq.clean %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family != "Mitochondria" &
      Class != "Chloroplast"
  )
# pseq.df <- as(tax_table(pseq), "matrix") # To check that correct taxa were removed
pseq
```

## 3.7 Filter to Remove Null Taxa
```{r}
# Sort samples on total read count
sort(sample_sums(pseq)) 

# Remove Null Taxa
pseq <- prune_taxa(taxa_sums(pseq) > 0, pseq) 
pseq

summarize_phyloseq(pseq) 

save(pseq, file ="devil.phyloseq")
load("devil.phyloseq")
pseq
```

# 4.0 Composition Plots
## 4.1 Transform Counts to Relative Abundance
```{r}
pseq.rel = transform_sample_counts(pseq, function(x) x/sum(x)*100)
```

## 4.2 Relative Abudance per Sample - Phylum
```{r}
glom.phy <- tax_glom(pseq.rel, taxrank = "Phylum", NArm = F)
pseq.melt.phy <- psmelt(glom.phy)
pseq.melt.phy$Phylum <- as.character(pseq.melt.phy$Phylum)

pseq.melt.phy <- pseq.melt.phy %>%
  group_by(reproductive.status, Phylum) %>%
  mutate(median = median(Abundance))

keep <- unique(pseq.melt.phy$Phylum[pseq.melt.phy$median > 0.25])
pseq.melt.phy$Phylum[!(pseq.melt.phy$Phylum %in% keep)] <- "Other (< 0.25%)"

pseq.melt.sum.phy <- pseq.melt.phy %>%
  group_by(Sample, reproductive.status, location, Phylum) %>%
  summarise(Abundance = sum(Abundance))

#move "Other" to end of list
pseq.melt.sum.phy$Phylum <- factor(pseq.melt.sum.phy$Phylum, levels = c("Actinobacteriota", "Bacteroidota", "Firmicutes", "Fusobacteriota", "Proteobacteria", "Other (< 0.25%)"))

#rename to order based on location
pseq.melt.sum.phy$new_Sample_ID <- sapply(strsplit(pseq.melt.sum.phy$Sample, "-"), function(x) paste(rev(x), collapse = "-")) 

#colour palette
colour_palette_phy <- c("#5254a3", "#c7c7c7", "#f7b6d2", "#dbdb8d", "#9edae5",
  "#6b6ecf")

#plot
ggplot(pseq.melt.sum.phy, aes(x = new_Sample_ID, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", aes(fill=Phylum)) + 
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="Reproductive Status", y="Relative Frequency (%)") +
  facet_wrap(~reproductive.status, scales= "free_x", nrow=1, strip.position = "bottom") +
  theme_classic() + 
    guides(fill = "none")+
  scale_fill_manual(values=colour_palette_phy) + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_rect(size = 1),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

ggsave("./Figures/sample-tax-phylum.png", height = 8, width = 12)
```

## 4.3 Relative Abudance per Sample - Family
```{r}
glom.fam <- tax_glom(pseq.rel, taxrank = "Family", NArm = F)
pseq.melt.fam <- psmelt(glom.fam)
pseq.melt.fam$Family <- as.character(pseq.melt.fam$Family)

pseq.melt.fam <- pseq.melt.fam %>%
  group_by(reproductive.status, Family) %>%
  mutate(median = median(Abundance))

keep <- unique(pseq.melt.fam$Family[pseq.melt.fam$median > 0.25])
pseq.melt.fam$Family[!(pseq.melt.fam$Family %in% keep)] <- "Other (< 0.25%)"

pseq.melt.sum.fam <- pseq.melt.fam %>%
  group_by(Sample, reproductive.status, Family) %>%
  summarise(Abundance = sum(Abundance))

#move "Other" to end of list
pseq.melt.sum.fam$Family <- factor(pseq.melt.sum.fam$Family, levels = c("Bacillaceae", "Carnobacteriaceae", "Clostridiaceae", "Corynebacteriaceae","Enterobacteriaceae", "Erysipelotrichaceae", "Family_XI", "Flavobacteriaceae", "Fusobacteriaceae", "Gemellaceae", "Lachnospiraceae", "Leptotrichiaceae", "Micrococcaceae", "Moraxellaceae", "Morganellaceae", "Mycoplasmataceae", "Neisseriaceae", "Pasteurellaceae", "Peptostreptococcaceae", "Planococcaceae", "Porphyromonadaceae", "Pseudomonadaceae", "Staphylococcaceae", "Streptococcaceae", "Vagococcaceae", "Weeksellaceae", "Wohlfahrtiimonadaceae", "Other (< 0.25%)"))

# order based on location
pseq.melt.sum.fam$new_Sample_ID <- sapply(strsplit(pseq.melt.sum.fam$Sample, "-"), function(x) paste(rev(x), collapse = "-")) 

# colour palette
colour_palette_fam <- c(
"#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#aec7e8", "#ffbb78", "#393b79", "#98df8a", "#ff9896", "#c5b0d5",
  "#c49c94", "#f7b6d2", "#5254a3", "#c7c7c7", "#dbdb8d", "#9edae5",
  "#6b6ecf"
)
#plot
ggplot(pseq.melt.sum.fam, aes(x = new_Sample_ID, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity", aes(fill=Family)) + 
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="Reproductive Status", y="Relative Frequency (%)") +
  facet_wrap(~reproductive.status, scales= "free_x", nrow=1, strip.position = "bottom") +
  theme_classic() + 
  scale_fill_manual(values=colour_palette_fam) + 
  guides(fill = "none")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_rect(size = 1),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

ggsave("./Figures/sample-tax-phylum.png", height = 8, width = 12)

ggsave("./Figures/sample-tax-family.png", height = 8, width = 12)
```

## 4.4 Relative Abudance per Sample - Genus
```{r}
glom.gen <- tax_glom(pseq.rel, taxrank = "Genus", NArm = F)
pseq.melt.gen <- psmelt(glom.gen)
pseq.melt.gen$Genus <- as.character(pseq.melt.gen$Genus)

pseq.melt.gen <- pseq.melt.gen %>%
  group_by(reproductive.status, Genus) %>%
  mutate(median = median(Abundance))

keep <- unique(pseq.melt.gen$Genus[pseq.melt.gen$median > 0.25])
pseq.melt.gen$Genus[!(pseq.melt.gen$Genus %in% keep)] <- "Other (< 0.25%)"

pseq.melt.sum.gen <- pseq.melt.gen %>%
  group_by(Sample, reproductive.status, Genus) %>%
  summarise(Abundance = sum(Abundance))

#move "Other" to end of list
pseq.melt.sum.gen$Genus <- factor(pseq.melt.sum.gen$Genus, levels = c("Acinetobacter", "Alysiella", "Bacillus", "Carnobacterium", "Cetobacterium", "Clostridium_sensu_stricto_1", "Clostridium_sensu_stricto_13", "Corynebacterium", "Fusobacterium", "Gemella", "Macrococcus", "Paeniclostridium", "Peptostreptococcus", "Plesiomonas", "Porphyromonas", "Pseudomonas", "Psychrobacter", "Romboutsia", "Staphylococcus", "Streptobacillus", "Streptococcus", "Terrisporobacter", "Vagococcus", "Other (< 0.25%)"))

# order based on location
pseq.melt.sum.gen$new_Sample_ID <- sapply(strsplit(pseq.melt.sum.gen$Sample, "-"), function(x) paste(rev(x), collapse = "-")) 

# colour palette
colour_palette_gen <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#5254a3", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#aec7e8", "#ffbb78", "#393b79", "#98df8a", "#ff9896", "#c5b0d5",
  "#c49c94", "#f7b6d2", "#dbdb8d", "#c7c7c7","#9edae5",
  "#6b6ecf"
)

#plot
ggplot(pseq.melt.sum.gen, aes(x = new_Sample_ID, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", aes(fill=Genus)) + 
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="Reproductive Status", y="Relative Frequency (%)") +
  facet_wrap(~reproductive.status, scales= "free_x", nrow=1, strip.position = "bottom") +
  theme_classic() + 
  scale_fill_manual(values=colour_palette_gen)+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_rect(size = 1),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

ggsave("./Figures/sample-tax-phylum.png", height = 8, width = 12)

ggsave("./Figures/sample-tax-genus.png", height = 8, width = 12)
```

# 5.0 Exploratory Analyses
## 5.1 Age
```{r}
# Generate metadata dataframe
clean.metadata <- data.frame(sample_data(pseq))

# Statistical test
t.test(devil.age ~ reproductive.status, data = clean.metadata)
## t = 6.5841, df = 31.819, p-value = 2.079e-07
## mean in group Lactating = 2.24, mean in group Non-Lactating = 1.13

## Standard deviations
lact.subset <- subset(clean.metadata, reproductive.status == "Lactating")
sd(lact.subset$devil.age)
nonlact.subset <- subset(clean.metadata, reproductive.status == "Non-Lactating")
sd(nonlact.subset$devil.age)

# Plot
age.plot <- ggplot(clean.metadata, aes(x=reproductive.status, y=devil.age)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() + geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = reproductive.status), height = 0) + xlab("Reproductive Status") + ylab("Age (years)") + theme_bw() + theme(legend.position = "none")
age.plot + scale_color_manual(values=c("darkorchid2", "darkorange1"))

ggsave("./Figures/age-barplot.png", height = 5, width = 7)
```

# 6.0 Alpha Diversity Analyses
## 6.1 Rarefy Reads
```{r}
pseq_rar <- rarefy_even_depth(pseq, rngseed = 123, replace = FALSE)

# Check
head(sample_sums(pseq_rar))
```

## 6.2 Generate Values
```{r}
shannon.index <- estimate_richness(pseq_rar, measures = "Shannon")
faiths.phylogenetic.diversity <- pd(samp = data.frame(t(data.frame(otu_table(pseq_rar)))), tree = phy_tree(pseq_rar), include.root = T) #pd already provides species richness (i.e. Observed ASVs)

clean.metadata02 <- cbind(clean.metadata, shannon.index, faiths.phylogenetic.diversity)
colnames(clean.metadata02)[11] <- "shannon.index"
colnames(clean.metadata02)[12] <- "faiths.phylogenetic.diversity"
colnames(clean.metadata02)[13] <- "observed.asvs"

write.csv(clean.metadata02, "./import-files/metadata-v2")

#subset data
lact.subset02 <- subset(clean.metadata02, reproductive.status == "Lactating")
nonlact.subset02 <- subset(clean.metadata02, reproductive.status == "Non-Lactating")

#mean and sd using data subsets
mean(lact.subset02$observed.asvs) #388.96
sd(lact.subset02$observed.asvs) #163.11
mean(nonlact.subset02$observed.asvs) #618.47
sd(nonlact.subset02$observed.asvs) #392.71

mean(lact.subset02$shannon.index) #4.51
sd(lact.subset02$shannon.index) #0.83
mean(nonlact.subset02$shannon.index) #5.05
sd(nonlact.subset02$shannon.index) #0.99

mean(lact.subset02$faiths.phylogenetic.diversity) #18.06
sd(lact.subset02$faiths.phylogenetic.diversity) #7.17
mean(nonlact.subset02$faiths.phylogenetic.diversity) #28.76
sd(nonlact.subset02$faiths.phylogenetic.diversity) #20.15
```

## 6.3 Plot
```{r}
observed <- ggplot(clean.metadata02, aes(x = reproductive.status, y = observed.asvs)) + stat_boxplot(geom ='errorbar') + geom_boxplot(outlier.colour = NA) + geom_jitter(aes(color = reproductive.status), height = 0) + labs(x = "\n", y = "Number of Observed ASVs") + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + ggtitle("Observed ASV") + scale_color_manual(values=c("darkorchid2", "darkorange1"))
observed

observed.nl <- nonlact.subset02 %>% ggplot(aes(x = location, y = faiths.phylogenetic.diversity)) + geom_boxplot() + geom_point(aes(x = location, y = observed.asvs))  # + facet_wrap(vars(nonlact.subset02$location))
observed.nl

shannon <- ggplot(clean.metadata02, aes(x = reproductive.status, y = shannon.index)) + stat_boxplot(geom ='errorbar') + geom_boxplot(outlier.colour = NA) + geom_jitter(aes(color = reproductive.status), height = 0) + labs(x = "\nReproductive Status", y = "Shannon Diversity Index") + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + ggtitle("Shannon Diversity") + scale_color_manual(values=c("darkorchid2", "darkorange1"))
shannon

faiths <- ggplot(clean.metadata02, aes(x = reproductive.status, y = faiths.phylogenetic.diversity)) + stat_boxplot(geom ='errorbar') + geom_boxplot(outlier.colour = NA) + geom_jitter(aes(color = reproductive.status), height = 0) + labs(x = "\n", y = "Faith's PD Index") + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + ggtitle("Faith's Phylogenetic Diversity") + scale_color_manual(values=c("darkorchid2", "darkorange1"))
faiths

alpha.diversity <- plot_grid(observed, shannon, faiths, labels = c("A", "B", "C"), ncol = 3, nrow = 1)
alpha.diversity

ggsave("./Figures/alpha_grid.png", width = 8, height = 6)
```

## 6.4 Check Significance
```{r}
# check repro status vs alpha metric
glm.observed.age <- glm(observed.asvs ~ reproductive.status + location, data = clean.metadata02)
summary(glm.observed)
glm.shannon <- glm(shannon.index ~ reproductive.status + location, data = clean.metadata02)
summary(glm.shannon)
glm.faiths <- glm(faiths.phylogenetic.diversity ~ reproductive.status + location, data = clean.metadata02)
summary(glm.faiths)

# confirm that significant effects are due to rep status and not age
glm.observed.age.l <- glm(observed.asvs ~ devil.age + location, data = lact.subset02)
summary(glm.observed.age.l)
glm.observed.age.nl <- glm(observed.asvs ~ devil.age + location, data = nonlact.subset02)
summary(glm.observed.age.nl)

glm.shannon.age.l <- glm(shannon.index ~ devil.age + location, data = lact.subset02)
summary(glm.shannon.age.l)
glm.shannon.age.nl <- glm(shannon.index ~ devil.age + location, data = nonlact.subset02)
summary(glm.shannon.age.nl)

glm.faiths.age.l <- glm(faiths.phylogenetic.diversity ~ devil.age + location, data = lact.subset02)
summary(glm.faiths.age.l)
glm.faiths.age.nl <- glm(faiths.phylogenetic.diversity ~ devil.age + location, data = nonlact.subset02)
summary(glm.faiths.age.nl)
```

# 7.0 Beta Diversity Analyses
## 7.1 Unweighted UniFrac
```{r}
set.seed(123)

uw.unifrac.distance <- phyloseq::distance(pseq_rar, method = "uunifrac")
```

### 7.1.1 Ordinate and Plot
```{r}
# NMDS
ordination.uw.nmds <- ordinate(pseq_rar, method = "NMDS", distance = uw.unifrac.distance)
uw.nmds <- plot_ordination(pseq_rar, ordination.uw.nmds, color = "reproductive.status", shape = "location") + theme_bw() + labs(shape="Location", colour="Reproductive Status") + annotate("text", x = -0.4, y = 0.3, label = paste0("Stress value: ", format(ordination$stress, digits = 4)), hjust = 0) + scale_color_manual(values=c("darkorchid2", "darkorange1")) + geom_point(size = 2.5)
uw.nmds
ggsave("./Figures/unweighted-NMDS.png", height = 6, width = 8)

# PCoA
ordination.uw.pcoa <- ordinate(pseq_rar, method = "PCoA", distance = uw.unifrac.distance)
uw.pcoa.1.2 <- plot_ordination(pseq_rar, ordination.uw.pcoa, color = "reproductive.status", shape = "location", axes = c(1,2)) + theme_bw() + labs(shape="Location", colour="Reproductive Status") + scale_color_manual(values=c("darkorchid2", "darkorange1")) + geom_point(size = 2.5)
uw.pcoa.1.2
uw.pcoa.1.3
uw.pcoa.2.3

# Scree Plot
p_scree.uw <- plot_ordination(pseq_rar, ordination.uw.pcoa, type = "scree")
p_scree.uw
ggsave("./Figures/scree-plot-uw-unifrac.png", height = 6, width = 8)
```

### 7.1.2 Significance Testing
```{r}
fb.subset <- subset_samples(pseq_rar, location == "Fentonbury")
uw.unifrac.distance.fb <- phyloseq::distance(fb.subset, method = "uunifrac")
w.unifrac.distance.fb <- phyloseq::distance(fb.subset, method = "wunifrac")

kp.subset <- subset_samples(pseq_rar, location == "Kempton")
uw.unifrac.distance.kp <- phyloseq::distance(kp.subset, method = "uunifrac")
w.unifrac.distance.kp <- phyloseq::distance(kp.subset, method = "wunifrac")

np.subset <- subset_samples(pseq_rar, location == "Narawntapu")
uw.unifrac.distance.np <- phyloseq::distance(np.subset, method = "uunifrac")
w.unifrac.distance.np <- phyloseq::distance(np.subset, method = "wunifrac")

# get new metadata values: 
fb.subset02 <- subset(clean.metadata02, location == "Fentonbury")
kp.subset02 <- subset(clean.metadata02, location == "Kempton")
np.subset02 <- subset(clean.metadata02, location == "Narawntapu")

#location permanova
##fentonbury
uw.uni.perm.fb <- adonis2(uw.unifrac.distance.fb ~ reproductive.status, data = fb.subset02)
uw.uni.perm.fb
w.uni.perm.fb <- adonis2(w.unifrac.distance.fb ~ reproductive.status, data = fb.subset02)
w.uni.perm.fb

ordination.uw.nmds.np <- ordinate(pseq_rar, method = "NMDS", distance = uw.unifrac.distance.np)
uw.nmds.np <- plot_ordination(pseq_rar, ordination.uw.nmds.np, color = "reproductive.status", title = "Unweighted UniFrac Distance", )
uw.nmds.np

ordination.w.nmds.np <- ordinate(pseq_rar, method = "NMDS", distance = w.unifrac.distance.np)
w.nmds.np <- plot_ordination(pseq_rar, ordination.w.nmds.np, color = "reproductive.status", title = "Weighted UniFrac Distance")
w.nmds.np

ggarrange(uw.nmds.np, w.nmds.np, ncol = 2, nrow = 1, common.legend = T)

##kempton
uw.uni.perm.kp <- adonis2(uw.unifrac.distance.kp ~ reproductive.status, data = kp.subset02)
uw.uni.perm.kp
w.uni.perm.kp <- adonis2(w.unifrac.distance.kp ~ reproductive.status, data = kp.subset02)
w.uni.perm.kp

##narawntapu
uw.uni.perm.np <- adonis2(uw.unifrac.distance.np ~ reproductive.status, data = np.subset02)
uw.uni.perm.np
w.uni.perm.np <- adonis2(w.unifrac.distance.np ~ reproductive.status, data = np.subset02)
w.uni.perm.np

# PERMANOVA
uw.unifrac.permanova <- adonis2(uw.unifrac.distance ~ reproductive.status + location, data = clean.metadata02)
uw.unifrac.permanova # repro.status and location significant
write.table(data.frame(uw.unifrac.permanova), "permanova-unweighted.txt", sep = "\t", dec = ".", col.names = NA)

## Pairwise test to determine significant groups
pairwise.uw <- pairwise.adonis2(uw.unifrac.distance ~ location, nperm = 9999, data = clean.metadata)
pairwise.uw
write.table(data.frame(pairwise.uw), "pairwise-unweighted.txt", sep = "\t", dec = ".", col.names = NA)

## Post-Hoc Testing
### Reproductive Status
uw.dispersion_repro <- betadisper(uw.unifrac.distance, clean.metadata$repro.status)
permutest(uw.dispersion_repro) # significant so null hypothesis that groups have same disperions cannot be rejected - not confident adonis result is a real result rather than due to differences in group dispersions
centroid.repro.uw <- plot(uw.dispersion_repro, main = NULL, sub = NULL)
distance.repro.uw <- boxplot(uw.dispersion_repro, xlab = "Reproductive Status", ylab = "Distance to Centroid")

# repro HSD not needed as only two groups

### Location
uw.dispersion_location <- betadisper(uw.unifrac.distance, clean.metadata$location)
permutest(uw.dispersion_location) # significant - as above
centroid.location.uw <- plot(uw.dispersion_location, main = NULL, sub = NULL)
distance.location.uw <- boxplot(uw.dispersion_location, xlab = "Reproductive Status", ylab = "Distance to Centroid")

uw.location.HSD <- TukeyHSD(uw.dispersion_location)
par(mar = c(5, 15, 5, 5))
plot(uw.location.HSD, las = 2)
```

## 7.2 Weighted UniFrac
```{r}
set.seed(123)

w.unifrac.distance <- phyloseq::distance(pseq_rar, method = "wunifrac")
```

```{r}
# NMDS
ordination.w = ordinate(pseq_rar, method = "NMDS", distance = w.unifrac.distance)
w.nmds <- plot_ordination(pseq_rar, ordination.w, color = "reproductive.status", shape = "location") + theme_bw() + labs(shape="Location", colour="Reproductive Status") + annotate("text", x = -0.01, y = 0.05, label = paste0("Stress value: ", format(ordination.w$stress, digits = 4))) + scale_color_manual(values=c("darkorchid2", "darkorange1")) + geom_point(size = 2.5)
w.nmds
ggsave("./Figures/weighted-nmds.png", height = 6, width = 8)

# PCoA
ordination.w.pcoa <- ordinate(pseq_rar, method = "PCoA", distance = w.unifrac.distance)
w.pcoa.2.3 <- plot_ordination(pseq_rar, ordination.w.pcoa, color = "reproductive.status", shape = "location", axes = c(2,3)) + theme_bw() + labs(shape="Location", colour="Reproductive Status") + scale_color_manual(values=c("darkorchid2", "darkorange1")) + geom_point(size = 2.5)
w.pcoa.1.2
w.pcoa.1.3
w.pcoa.2.3

## Arrange PCoA
ggarrange(w.pcoa.1.2, w.pcoa.1.3, w.pcoa.2.3,  uw.pcoa.1.2, uw.pcoa.1.3, uw.pcoa.2.3, nrow = 2, ncol = 3, labels = c("A", "B", "C", "D", "E", "F"), label.y = c("Weighted UniFrac", "Unweighted UniFrac"), widths=c(1,1), common.legend = T, legend = "right")

ggsave("./Figures/weighted-PCoA-2-3.png", height = 6, width = 8)

# Scree Plot
p_scree.w <- plot_ordination(pseq_rar, ordination.w.pcoa, type = "scree")
p_scree.w

ggsave("./Figures/scree-plot-uw-unifrac.png", height = 6, width = 8)


# PERMANOVA
w.unifrac.permanova <- adonis2(w.unifrac.distance ~ reproductive.status + location, data = clean.metadata)
w.unifrac.permanova # location significant
write.table(data.frame(w.unifrac.permanova), "permanova-weighted.txt", sep = "\t", dec = ".", col.names = NA)

## Pairwise test to determine significant groups
pairwise.adonis2(w.unifrac.distance ~ location, nperm = 999, data = clean.metadata)

## Post-Hoc Testing
### Reproductive Status
w.dispersion_repro <- betadisper(w.unifrac.distance, clean.metadata$repro.status)
permutest(w.dispersion_repro) # significant so null hypothesis that groups have same disperions cannot be rejected - not confident adonis result is a real result rather than due to differences in group dispersions

centroid.repro.w <- plot(w.dispersion_repro, main = NULL, sub = NULL, xlim = c(-0.020, 0.020))
distance.repro.w <- boxplot(w.dispersion_repro, xlab = "Reproductive Status", ylab = "Distance to Centroid")

# repro HSD not needed as only two groups

### Location
w.dispersion_location <- betadisper(w.unifrac.distance, clean.metadata$location)
permutest(w.dispersion_location) # not significant - dispersions are equal
centroid.location.w <- plot(w.dispersion_location, main = NULL, sub = NULL, xlim = c(-0.020, 0.020))
distance.location.w <- boxplot(w.dispersion_location, xlab = "Reproductive Status", ylab = "Distance to Centroid")

# location HSD not needed as dispersion equal
```

### 7.2.1 Ordinate and Plot
```{r}
# NMDS
ordination.uw.nmds <- ordinate(pseq_rar, method = "NMDS", distance = uw.unifrac.distance)
uw.nmds <- plot_ordination(pseq_rar, ordination.uw.nmds, color = "reproductive.status", shape = "location") + theme_bw() + labs(shape="Location", colour="Reproductive Status") + annotate("text", x = -0.4, y = 0.3, label = paste0("Stress value: ", format(ordination$stress, digits = 4)), hjust = 0) + scale_color_manual(values=c("darkorchid2", "darkorange1")) + geom_point(size = 2.5)
uw.nmds
ggsave("./Figures/unweighted-NMDS.png", height = 6, width = 8)

# PCoA
ordination.uw.pcoa <- ordinate(pseq_rar, method = "PCoA", distance = uw.unifrac.distance)
uw.pcoa <- plot_ordination(pseq_rar, ordination.uw.pcoa, color = "reproductive.status", shape = "location") + theme_bw() + labs(shape="Location", colour="Reproductive Status") + scale_color_manual(values=c("darkorchid2", "darkorange1")) + geom_point(size = 2.5)
uw.pcoa
ggsave("./Figures/unweighted-PCoA.png", height = 6, width = 8)

# Scree Plot
p_scree.uw <- plot_ordination(pseq_rar, ordination.uw.pcoa, type = "scree")
p_scree.uw

ggsave("./Figures/scree-plot-uw-unifrac.png", height = 6, width = 8)
```


## 7.3 UniFrac Combined Plot
```{r}
# Group Plot
ggarrange(uw.nmds, w.nmds, ncol = 1, nrow = 2, labels = c("A", "B"), widths=c(1,1), common.legend = T, legend = "right")
ggsave("./Figures/group-unifrac.png", height = 9, width = 8)

plot_grid(centroid.location.uw, centroid.location.w, distance.location.uw, distance.location.w, labels = c("A", "B", "C", "D"))
```
