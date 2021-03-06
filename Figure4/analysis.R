install.packages("ggplot2")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("phyloseq")
install.packages("ape")
install.packages("microbiome")
install.packages("vegan")
install.packages("ggpubr")
library(ggplot2)
library(dplyr)
library(tidyverse)
library(phyloseq)
library(ape)
library(microbiome)
library(vegan)
library(ggpubr)

# https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/multivariate-comparisons-of-microbial-community-composition.html

directory <- getwd()
setwd(directory)

### Physicochemical data
chem_data <- read.csv(file = "chem_data.csv",header = TRUE)

# Import tables and format for PhylloSeq
bracken <- read_csv("count_table.csv") 

# Import taxonomic ID table for PhylloSeq
tax_table <- read_csv("taxonomy_table.csv") %>% rename("kingdom"="superkingdom")

### Removing contaminate species
remove_taxa <- read_csv("remove_taxa.csv") 
contaminate_id_list = as.character(remove_taxa$taxID)
bracken = filter(bracken, !taxonomy_id %in% contaminate_id_list)
# Remove contaminates from Tax Table to fix index issue for figure
tax_table = filter(tax_table, !taxID %in% contaminate_id_list) 

## create TaxTable
taxMat <- tax_table %>% tibble::column_to_rownames("taxID")
TAX = tax_table(as.matrix(taxMat))

## create countTable
countMat <- bracken %>% tibble::column_to_rownames("taxonomy_id")
OTU = otu_table(as.matrix(countMat), taxa_are_rows = TRUE)

chem_data <- chem_data %>%tibble::column_to_rownames("File.Name")
SAMPLE = sample_data(as.data.frame(chem_data))

physeq = phyloseq(OTU, TAX, SAMPLE)

# summarize the table
microbiome::summarize_phyloseq(physeq)

# remove almost low abundance samples
reads_sample <- readcount(physeq)
sample_data(physeq)$reads_sample <- reads_sample
pseq.subset <- subset_samples(physeq, reads_sample > 100000)

# get the normalized counts
pseq.compo <- microbiome::transform(pseq.subset, "compositional")

## PERMANOVA analysis
counts <- abundances(pseq.compo)
SpeciesNames <- taxMat[rownames(taxMat), "species"]
row.names(counts) <- SpeciesNames
meta <- meta(pseq.compo)
set.seed(1)

# Calculate bray curtis distance matrix
B1_bray <- phyloseq::distance(pseq.compo, method = "bray")

# Permanova for Depth
# Adonis test 
permanova <- adonis(t(counts) ~ Depth,
                    data = meta, permutations=999, method = "bray")

# P-value
print(as.data.frame(permanova$aov.tab)["Depth", "Pr(>F)"])

# Homogeneity of dispersion test
anova(betadisper(B1_bray, meta$Depth))

# Select top 20 coefficients
coef <- coefficients(permanova)["Depth",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
cf <- as.data.frame(top.coef)
cf <- tibble::rownames_to_column(cf, "species")
cf <- left_join(cf, tax_table, by = "species")

# Use this one for figure
plt <- ggdotchart(cf, x = "species", y = "top.coef",
           color = "phylum",                                # Color by groups
           palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           rotate = TRUE,                                # Rotate vertically
           dot.size = 3.5,                                 # Large dot size
           y.text.col = TRUE,                            # Color y text by groups
           ggtheme = theme_pubr()                        # ggplot2 theme
)+
  theme_cleveland() 

plot(plt)

path_string = paste(directory,"permanova_figure.pdf",sep="/")

#save plots
ggsave(filename = path_string, 
       plot =  plt, width = 29.7, height = 21, units = "cm")

