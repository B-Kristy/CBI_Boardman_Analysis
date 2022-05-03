# Load Required Packages
library(compositions)
library(cowplot)
library(cluster)
library(data.table)
library(devtools)
library(plyr)
library(dplyr)
library(doParallel)
library(factoextra)
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(hillR)
library(microbiome)
library(multcomp)
library(nlme)
library(phyloseq)
library(phylosmith)
library(qiime2R)
library(tidyverse)
library(vegan)
library(zCompositions)

source("~/OneDrive/Documents/R/win-library/4.1/ANCOM-master/programs/ancom.R")

#Set-up ------------------------------------------------------------------------
# Set Working Directory

setwd("~/Dropbox/Cregger_CBI_drought/data/")

# Read QZA files into dataframe, re-format taxonomic tables, and re-upload them as .csv files

SVs16S <- read_qza("../qiime_output/16S-merged-table.qza")
SVs16Stable <- SVs16S$data
write.csv(SVs16Stable, file = "16S-merged-table.csv")
taxonomy16S <- read_qza("../qiime_output/16S-taxonomy.qza")
tax16S<-taxonomy16S$data %>% as.tibble() %>%
  mutate(Taxon=gsub("D_[0-9]__", "", Taxon)) %>% 
  separate(Taxon, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))%>%
  mutate(Phylum=replace_na(Phylum,"empty"))
write.csv(tax16S, file = "16S-taxonomy.csv", row.names =F)

SVsITS <- read_qza("../qiime_output/ITS-merged-table.qza")
SVsITStable <- SVsITS$data
write.csv(SVsITStable, file = "ITS-merged-table.csv")
taxonomyITS <- read_qza("../qiime_output/ITS-taxonomy.qza")
taxITS<-taxonomyITS$data %>% as.tibble() %>%
  mutate(Taxon=gsub("[a-z]__", "", Taxon)) %>% 
  separate(Taxon, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) %>%
  mutate(Phylum=replace_na(Phylum,"unidentified"))
write.csv(taxITS, file= "ITS-taxonomy.csv", row.names=F)



# Filtering --------------------------------------------------------------------

# 16S Library Filtering

# Write unfiltered data into a phyloseq object

data_16S_unfiltered <- read_csv2phyloseq(otu.file = "16S-merged-table.csv",
                                         taxonomy.file = "16S-taxonomy.csv",
                                         metadata.file = "16S_sample_metadata.csv")
summarize_phyloseq(data_16S_unfiltered)
# Filter out Plant, Mi1ochondrial, Unassigned DNA from 16S Dataset

data_16S_uf1 <- subset_taxa(data_16S_unfiltered, Kingdom == "Bacteria" | Kingdom == "Archaea")
data_16S_uf2 <- subset_taxa(data_16S_uf1, Kingdom != "Eukaryota")
data_16S_uf3 <- subset_taxa(data_16S_uf2, Order != "Chloroplast")
data_16S_uf4 <- subset_taxa(data_16S_uf3, Family != "Mitochondria")
summarize_phyloseq(data_16S_uf4)

# Remove samples with extremely low read depth 

data_16S_uf5 <- prune_samples(sample_sums(data_16S_uf4)>=1000, data_16S_uf4)

# Export filtered data

write.csv(data_16S_uf5@otu_table, "16S-merged-table-filtered.csv")
write.csv(data_16S_uf5@tax_table, "16S-taxonomy.csv")

# ITS Library Filtering

# Write unfiltered data into a phyloseq object

data_ITS_unfiltered <- read_csv2phyloseq(otu.file = "ITS-merged-table.csv",
                                         taxonomy.file = "ITS-taxonomy.csv",
                                         metadata.file = "ITS_sample_metadata.csv")
summarize_phyloseq(data_ITS_unfiltered)

# Filter out Plant, Micochondrial, or Unassigned DNA from ITS Dataset 

data_ITS_uf1 <- subset_taxa(data_ITS_unfiltered, Kingdom == "Fungi")

# Remove samples with extremely low read depth

data_ITS_uf2 <- prune_samples(sample_sums(data_ITS_uf1)>=1000, data_ITS_uf1)

# Export filtered data

write.csv(data_ITS_uf2@otu_table, "ITS-merged-table-filtered.csv")
write.csv(data_ITS_uf2@tax_table, "ITS-taxonomy.csv")

# Create phyloseq objects of filtered data (16S and ITS)
data_16S_filtered <- read_csv2phyloseq(otu.file = "16S-merged-table-filtered.csv",
                                       taxonomy.file = "16S-taxonomy.csv",
                                       metadata.file = "16S_sample_metadata.csv")
summarize_phyloseq(data_16S_filtered)


sample_data(data_16S_filtered)$habitat <- ordered(sample_data(data_16S_filtered)$habitat, c("soil", "root_endosphere", "rhizosphere"))
sample_data(data_16S_filtered)$irrigation <- ordered(sample_data(data_16S_filtered)$irrigation, c("full", "reduced"))
sample_data(data_16S_filtered)$drought_tolerance <- ordered(sample_data(data_16S_filtered)$drought_tolerance, c("HI30", "LO30"))

data_ITS_filtered <- read_csv2phyloseq(otu.file = "ITS-merged-table-filtered.csv", 
                                       taxonomy.file = "ITS-taxonomy.csv", 
                                       metadata.file = "ITS_sample_metadata.csv")
summarize_phyloseq(data_ITS_filtered)

sample_data(data_ITS_filtered)$habitat <- ordered(sample_data(data_ITS_filtered)$habitat, c("soil", "root_endosphere", "rhizosphere"))
sample_data(data_ITS_filtered)$irrigation <- ordered(sample_data(data_ITS_filtered)$irrigation, c("full", "reduced"))
sample_data(data_ITS_filtered)$drought_tolerance <- ordered(sample_data(data_ITS_filtered)$drought_tolerance, c("HI30", "LO30"))

# Rarefaction ------------------------------------------------------------------

# Prune SVs that are not present 5 times in at least 2 samples

data_16S_counts <- data_16S_filtered
data_16S_counts <- prune_samples(sample_sums(data_16S_counts)>=1000, data_16S_counts)
data_16S_filtered <- transform_sample_counts(data_16S_filtered, function(x) x/sum(x))
data_16S_counts <- filter_taxa(data_16S_counts, function(x) sum(x >5) > (0.01058201*length(x)), TRUE)
summarize_phyloseq(data_16S_counts)

sum(colSums(otu_table(data_16S_counts)))
sort(colSums(otu_table(data_16S_counts)))

data_ITS_counts <- data_ITS_filtered
data_ITS_counts <- prune_samples(sample_sums(data_ITS_counts)>=1000, data_ITS_counts)
data_ITS_filtered <- transform_sample_counts(data_ITS_filtered, function(x) x/sum(x))
data_ITS_counts <- filter_taxa(data_ITS_counts, function(x) sum(x>5) > (0.01058201*length(x)), TRUE)
summarize_phyloseq(data_ITS_counts)

sum(colSums(otu_table(data_ITS_counts)))
sort(colSums(otu_table(data_ITS_counts)))

# Rarefaction Curve Function

calculate_rarefaction_curves <- function(psdata, measures, depths, parallel=T) {
  require('plyr') # ldply
  require('reshape2') # melt
  require('doParallel')
  
  # set parallel options if required
  if (parallel) {
    paropts  <- list(.packages=c("phyloseq", "reshape2"))
  } else {
    paropts  <- NULL
  }
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive() && ! parallel, 'text', 'none'), .parallel=parallel, .paropts=paropts)
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

# 16S Rarefaction

detectCores(all.tests=TRUE)

cl <- makeCluster(detectCores(all.tests=TRUE),setup_strategy = "sequential")
registerDoParallel(cl)

rarefaction_curve_data <- calculate_rarefaction_curves(data_16S_counts, c('Observed', 'Shannon'), 
                                                       rep(c(1, 10, 100, 1000, 1:100 * 10000), each = 10))

rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, 
                                        Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))

rarefaction_curve_data_summary_verbose <-  merge(rarefaction_curve_data_summary %>% mutate(Sample = gsub("[.]", "-", Sample)),
                                                 data.frame(sample_data(data_16S_counts)) %>% 
                                                   rownames_to_column(var = "rowname"),
                                                 by.x = 'Sample', by.y = 'rowname')

discrete_palettes <- list(
  RColorBrewer::brewer.pal(3, "Set2")
)
rarefaction_curve_data_summary_verbose$habitat <- factor(rarefaction_curve_data_summary_verbose$habitat, levels = c("root_endosphere", "rhizosphere", "soil"))
curve_16S <- ggplot(rarefaction_curve_data_summary_verbose %>% filter(Measure == "Observed"), 
                    aes(x = Depth, y = Alpha_diversity_mean,  ymin = Alpha_diversity_mean - Alpha_diversity_sd, 
                        ymax = Alpha_diversity_mean + Alpha_diversity_sd, colour = habitat, group = Sample)) + 
  geom_line(alpha=1) +
  geom_pointrange() + 
  scale_x_continuous(trans = "log10", name = "Sequence Depth") +
  ylab("Richness") +
  facet_wrap(~habitat, labeller = labeller(habitat = c("rhizosphere" = "Rhizosphere", 
                                                       "root_endosphere" = "Root Endosphere", 
                                                       "soil" = "Soil")),
             scales = 'free') +
  scale_color_manual("Microbiome Compartment",
                     values = c("#FC8D62", "#8DA0CB", "#66C2A5"),
                     labels = c(bquote(paste("Root Endosphere", sep="")),
                                bquote(paste("Rhizosphere", sep="")),
                                bquote(paste("Soil", sep="")))) +
  theme(legend.text = element_text(hjust = 0)) +
  labs(title = "16S Rarefaction Curves") +
  theme_bw()
curve_16S
# ITS Rarefaction

detectCores(all.tests=TRUE)

cl <- makeCluster(detectCores(all.tests=TRUE),setup_strategy = "sequential")
registerDoParallel(cl)

rarefaction_curve_data <- calculate_rarefaction_curves(data_ITS_counts, c('Observed', 'Shannon'), 
                                                       rep(c(1, 10, 100, 1000, 1:100 * 10000), each = 10))

rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, 
                                        Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))

rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary %>% mutate(Sample = gsub("[.]", "-", Sample)),
                                                data.frame(sample_data(data_ITS_counts)) %>% 
                                                  rownames_to_column(var = "rowname"),
                                                by.x = 'Sample', by.y = 'rowname')


rarefaction_curve_data_summary_verbose$habitat <- factor(rarefaction_curve_data_summary_verbose$habitat, levels = c("root_endosphere", "rhizosphere", "soil"))
curve_ITS <- ggplot(rarefaction_curve_data_summary_verbose %>% filter(Measure == "Observed"), 
                    aes(x = Depth, y = Alpha_diversity_mean,  ymin = Alpha_diversity_mean - Alpha_diversity_sd, 
                        ymax = Alpha_diversity_mean + Alpha_diversity_sd, colour = habitat, group = Sample)) + 
  geom_line(alpha=1) +
  geom_pointrange() + 
  scale_x_continuous(trans = "log10", name = "Sequence Depth") +
  ylab("Richness") +
  facet_wrap(~habitat, labeller = labeller(habitat = c("rhizosphere" = "Rhizosphere", 
                                                       "root_endosphere" = "Root Endosphere", 
                                                       "soil" = "Soil")),
             scales = 'free') +
  scale_color_manual("Microbiome Compartment",
                     values = c("#FC8D62", "#8DA0CB", "#66C2A5"),
                     labels = c(bquote(paste("Root Endosphere", sep="")),
                                bquote(paste("Rhizosphere", sep="")),
                                bquote(paste("Soil", sep="")))) +
  theme(legend.text = element_text(hjust = 0)) +
  labs(title = "ITS Rarefaction Curves") +
  theme_bw()

curve_ITS
# Append both 16S and ITS Plots together 

prow <- plot_grid(curve_16S + theme(legend.position = "none"), curve_ITS + theme(legend.position = "none"), ncol = 1, labels = "AUTO", align = "vh")
leg <- get_legend(curve_16S)
prow <- plot_grid(prow, leg, ncol = 2, rel_widths = c(1, 0.23))
ggsave(filename=paste("rarefaction_curves.tiff", sep=""), plot=prow, width = 9.4, height=6.5, dpi=600)

# Taxonomy --------------------------------------------------------------
# 16S Habitat Overview Plot
# Rarefy data to 1000 reads 
set.seed(01221990)
rare_16S_habitat <- rarefy_even_depth(data_16S_counts, sample.size = 1000)

data_16S_rel <- transform_sample_counts(rare_16S_habitat, function(x) x/sum(x))

phylum_sum_16S <- tapply(taxa_sums(data_16S_rel), 
                         tax_table(data_16S_rel)[, "Phylum"], sum, na.rm=TRUE)
top10Phylum_16S <- names(sort(phylum_sum_16S, TRUE))[1:10]

abundant_16S_phylum <- subset_taxa(data_16S_rel, Phylum %in% top10Phylum_16S)
abundant_16S_phylum_glom <- tax_glom(abundant_16S_phylum, taxrank = "Phylum")
abundant_16S_phylum_melt <- psmelt(abundant_16S_phylum_glom)

sum2 <- abundant_16S_phylum_melt %>% 
  group_by(habitat, Phylum) %>%
  summarise(means = mean(Abundance))

#Everything that's left is considered rare  
rare <- sum2 %>%
  group_by(habitat) %>%
  summarise(means = 1- sum(means)) %>%
  mutate(Phylum = "Rare Taxa")

#concatenate the datasets
sum2 = rbind(sum2, rare)

#order groups
sum2$Phylum <- forcats::fct_relevel(sum2$Phylum, "Rare Taxa", after = Inf)

sum2$habitat<- dplyr::recode(sum2$habitat, soil = "Soil")
sum2$habitat<- dplyr::recode(sum2$habitat, root_endosphere = "Root Endosphere")
sum2$habitat<- dplyr::recode(sum2$habitat, rhizosphere = "Rhizosphere")

#Order habitat
sum2$habitat<- ordered(sum2$habitat, c("Root Endosphere", "Rhizosphere", "Soil"))


p2 <- ggplot(sum2, aes(x = habitat, y = means, fill = Phylum)) +
  geom_bar(position = "stack", stat = "identity", col = "black") +
  theme(legend.position = "right", legend.title = element_blank(), axis.line = element_blank()) +
  scale_fill_manual(values = c("firebrick4","tomato", "darkorange3", "orange1" ,"goldenrod4", "khaki1", "forestgreen",
                               "olivedrab1", "midnightblue", "powderblue", "mediumorchid4")) +
  scale_y_continuous(name = NULL, breaks = NULL) +
  panel_border() +
  theme_bw() +
  ggtitle("") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text = element_text(size = 12, face = "italic"), axis.text.x = element_text(size = 10),
        axis.title = element_blank(),
        title = element_text(size = 12))

# ITS Habitat Overview Plot
# Rarefy data to 1000 reads
set.seed(01221990)
rare_ITS_habitat <- rarefy_even_depth(data_ITS_counts, sample.size = 1000)

data_ITS_rel <- transform_sample_counts(rare_ITS_habitat, function(x) x/sum(x))

phylum_sum_ITS <- tapply(taxa_sums(data_ITS_rel), 
                         tax_table(data_ITS_rel)[, "Phylum"], sum, na.rm=TRUE)
top10Phylum_ITS <- names(sort(phylum_sum_ITS, TRUE))[1:7]

abundant_ITS_phylum <- subset_taxa(data_ITS_rel, Phylum %in% top10Phylum_ITS)
abundant_ITS_phylum_glom <- tax_glom(abundant_ITS_phylum, taxrank = "Phylum")
abundant_ITS_phylum_melt <- psmelt(abundant_ITS_phylum_glom)

sum2 <- abundant_ITS_phylum_melt %>% 
  group_by(habitat, Phylum) %>%
  summarise(means = mean(Abundance))

#Everything that's left is considered rare  
rare <- sum2 %>%
  group_by(habitat) %>%
  summarise(means = 1- sum(means)) %>%
  mutate(Phylum = "Rare Taxa")

#concatenate the datasets
sum2 = rbind(sum2, rare)

#order groups
sum2$Phylum <- dplyr::recode(sum2$Phylum, unidentified = "Unidentified Fungi")
sum2$Phylum <- forcats::fct_relevel(sum2$Phylum, "Rare Taxa", after = Inf)

sum2$habitat<- dplyr::recode(sum2$habitat, soil = "Soil")
sum2$habitat<- dplyr::recode(sum2$habitat, root_endosphere = "Root Endosphere")
sum2$habitat<- dplyr::recode(sum2$habitat, rhizosphere = "Rhizosphere")

#Order habitat
sum2$habitat<- ordered(sum2$habitat, c("Root Endosphere", "Rhizosphere", "Soil"))


p1 <- ggplot(sum2, aes(x = habitat, y = means, fill = Phylum)) +
  geom_bar(position = "stack", stat = "identity", col = "black") +
  theme(legend.position = "right", legend.title = element_blank(), axis.line = element_blank()) +
  scale_fill_manual(values = c("firebrick4","tomato","goldenrod4", "forestgreen",
                               "olivedrab1", "midnightblue", "powderblue", "mediumorchid4")) +
  scale_y_continuous(name = NULL, breaks = NULL) +
  panel_border() +
  theme_bw() +
  ggtitle("") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text = element_text(size = 12, face = "italic"), axis.text.x = element_text(size = 10),
        axis.title = element_blank(),
        title = element_text(size = 12))

# Append phylum plots together
prow <- plot_grid(p2 + theme(legend.position = "right"), p1 + theme(legend.position = "right"), ncol = 2, labels = "AUTO", align = "h")
prow
ggsave(filename=paste("phylum_taxonomy.tiff", sep=""), plot=prow, width=8.7, height=6, dpi=600)

# Alpha Diversity - Habitat ----------------------------------------------------
# 16S Alpha Diversity - Habitat Variability 

OTU_16S <- t(as(rare_16S_habitat@otu_table, "matrix"))
OTU_16S <- as.data.frame(OTU_16S)

rich_16S <- hill_taxa(OTU_16S, q = 0)
rich_16S <- merge(data_16S_counts@sam_data, rich_16S, by = 0)

shan_16S <- hill_taxa(OTU_16S, q = 1)
shan_16S <- merge(data_16S_counts@sam_data, shan_16S, by = 0)

simp_16S <- hill_taxa(OTU_16S, q=2)
simp_16S <- merge(data_16S_counts@sam_data, simp_16S, by=0)

# 16S Alpha Diversity Statistics 
model <- lme(y ~ habitat, random=~1|Row.names, data = rich_16S, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) # Habitat influences species richness (p < 0.05)
summary(glht(model, mcp(habitat="Tukey"))) # (Rhizosphere and Soil Species Richness are not significantly different)

model <- lme(y ~ habitat, random=~1|Row.names, data = shan_16S, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) #Habitat influences shannon diversity (p<0.05)
summary(glht(model, mcp(habitat="Tukey")))  #(Rhizosphere and Soil shannon diversity are not significantly different )

model <- lme(y ~ habitat, random=~1|Row.names, data = simp_16S, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) # Habitat influences simpson's diversity index (p <0.05)
summary(glht(model, mcp(habitat="Tukey"))) # All means for each habitat are significantly different 

# ITS Alpha Diversity - Habitat variability 

OTU_ITS <- t(as(rare_ITS_habitat@otu_table, "matrix"))
OTU_ITS <- as.data.frame(OTU_ITS)

rich_ITS <- hill_taxa(OTU_ITS, q = 0)
rich_ITS <- merge(data_ITS_counts@sam_data, rich_ITS, by = 0)

shan_ITS <- hill_taxa(OTU_ITS, q = 1)
shan_ITS <- merge(data_ITS_counts@sam_data, shan_ITS, by = 0)

simp_ITS <- hill_taxa(OTU_ITS, q=2)
simp_ITS <- merge(data_ITS_counts@sam_data, simp_ITS, by=0)

# ITS Alpha Diversity Statistics

model <- lme(y ~ habitat, random=~1|Row.names, data = rich_ITS, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) #All habitats significantly vary in species richness (p <0.05)
summary(glht(model, mcp(habitat="Tukey"))) 

model <- lme(y ~ habitat, random=~1|Row.names, data = shan_ITS, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) #Habitat influences shannon diversity (p<0.05)
summary(glht(model, mcp(habitat="Tukey")))  #(Rhizosphere and Soil shannon diversity are not significantly different )

model <- lme(y ~ habitat, random=~1|Row.names, data = simp_ITS, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) # Habitat influences simpson's diversity index (p <0.05)
summary(glht(model, mcp(habitat="Tukey"))) # Rhizosphere and Soil simpson diversity means are not significantly different 

# Plot Alpha Diversity
rich_16S$hillr <- "0"
shan_16S$hillr <- "1"
simp_16S$hillr <- "2"
bac_arc<- rbind(rich_16S %>% dplyr::select(habitat, y, Row.names, hillr),
                shan_16S %>% dplyr::select(habitat, y, Row.names, hillr),
                simp_16S %>% dplyr::select(habitat, y, Row.names, hillr))

bac_arc$habitat <- factor(bac_arc$habitat, levels = c("root_endosphere", "rhizosphere", "soil"))
bac_arc %>%
  ggplot(aes(x = habitat, y = y, fill=habitat)) +
  geom_boxplot(alpha=1) + 
  facet_grid(rows = vars(hillr),
             labeller = label_bquote(~phantom()^.(hillr)*italic(D)),
             scales = "free") +
  scale_fill_manual("Plant Habitat",
                    values = c("#FC8D62", "#8DA0CB", "#66C2A5"),
                    labels = c(expression(paste("Root Endosphere", sep="")),
                               expression(paste("Rhizosphere", sep="")),
                               expression(paste("Soil", sep="")))) +
  scale_x_discrete(name = NULL, labels = c("Root Endosphere", "Rhizosphere", "Soil")) +
  scale_y_continuous(name = "Archaeal/Bacterial Alpha Diversity") +
  panel_border() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        legend.title = element_text(size = 12)) -> p_bac_arc


rich_ITS$hillr <- "0"
shan_ITS$hillr <- "1"
simp_ITS$hillr <- "2"

its <- rbind(   rich_ITS %>% dplyr::select(habitat, y, Row.names, hillr),
                shan_ITS %>% dplyr::select(habitat, y, Row.names, hillr),
                simp_ITS %>% dplyr::select(habitat, y, Row.names, hillr))

its$habitat <- factor(its$habitat, levels = c("root_endosphere", "rhizosphere", "soil"))

its %>%
  ggplot(aes(x = habitat, y = y, fill=habitat)) +
  geom_boxplot(alpha=1) +
  facet_grid(rows = vars(hillr),
             labeller = label_bquote(~phantom()^.(hillr)*italic(D)),
             scales = "free") +
  scale_x_discrete(name = NULL, labels = c("Root Endosphere", "Rhizosphere", "Soil")) +
  scale_y_continuous(name = "Fungal Alpha Diversity") +
  scale_fill_manual("Plant Habitat",
                    values = c("#FC8D62", "#8DA0CB", "#66C2A5"),
                    labels = c(expression(paste("Root Endosphere", sep="")),
                               expression(paste("Rhizosphere", sep="")),
                               expression(paste("Soil", sep="")))) +
  panel_border() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        legend.title = element_text(size = 12)) -> p_its
p_its

prow <- plot_grid(p_bac_arc + theme(legend.position = "none"), p_its + theme(legend.position = "none"),
                  ncol = 2, align = "vh", labels ="AUTO")
habitat_alpha <- plot_grid(prow, ncol = 1, rel_widths = c(1))

habitat_alpha
ggsave(filename=paste("alpha_diversity_habitat.tiff",sep=""), plot=habitat_alpha, width=6.5, height=6.5, dpi=600)

# Beta Diversity - Habitat -----------------------------------------------------

# dbRDA Plot Using Jaccard

df_16S <- data.frame(sample_data(rare_16S_habitat))
OTU_16S <- t(as(rare_16S_habitat@otu_table, "matrix"))
OTU_16S <- as.data.frame(OTU_16S)
dist_16S <- vegdist(OTU_16S, "jaccard")

set.seed(1221990)
rda_16S <- dbrda(dist_16S ~ habitat*irrigation*drought_tolerance, df_16S, permutations = perm) 

anova(rda_16S) 
print(anova(rda_16S, by="terms", permutation = 9999))

y <- varpart(dist_16S, ~irrigation, ~drought_tolerance, ~habitat, data = df_16S)
print(y)

# Measure for dispersion
dispersion_16S <- betadisper(dist_16S, df_16S$habitat, type = c("centroid"), bias.adjust = TRUE, sqrt.dist = FALSE, add = FALSE)
anova(dispersion_16S)

set.seed(1221990)
permutest(dispersion_16S, pairwise = TRUE, permutations = 9999)
(dispersion_16S_HSD <- TukeyHSD(dispersion_16S))

# Plot Dispersion
data <- data.frame(Distance_to_centroid=dispersion_16S$distances, Group=dispersion_16S$group)
dispersion_16S$group <- factor(dispersion_16S$group, levels = c("root_endosphere", "rhizosphere", "soil"))
groups <- dispersion_16S$group
level_order <- c("root_endosphere", "rhizosphere", "soil")
disper_16S_plot <- ggplot(data=data, aes(x= factor(Group, level = level_order), y=Distance_to_centroid, colour = groups)) + geom_boxplot(alpha=0.5) + geom_jitter(shape = 16, size =2, position = position_jitter(0.2)) +
  scale_color_manual("Microbiome Compartment",
                     values = c("#FC8D62", "#8DA0CB", "#66C2A5"),
                     labels = c(expression(paste("Root Endosphere", sep="")),
                                expression(paste("Rhizosphere", sep="")),
                                expression(paste("Soil", sep="")))) +
  scale_x_discrete(name = NULL, labels = c("Root Endosphere", "Rhizosphere", "Soil")) +
  scale_y_continuous(name = "Distance to Centroid") +
  ggtitle("16S Dispersion") +
  panel_border() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12))

# Plot Community Composition
plot.data <- merge(summary(rda_16S)$sites, df_16S, by = "row.names") 
plot.data$habitat <- factor(plot.data$habitat, levels = c("root_endosphere", "rhizosphere", "soil"))
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, color = habitat), 
             size = 2, alpha = 1) +
  ggtitle("16S Composition") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_16S)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_16S)$cont)[2,2], 1)), "%)")) +
  theme_bw() +
  scale_color_manual("Microbiome Compartment",
                     values = c("#FC8D62", "#8DA0CB", "#66C2A5"),
                     labels = c(expression(paste("Root Endosphere", sep="")),
                                expression(paste("Rhizosphere", sep="")),
                                expression(paste("Soil", sep="")))) +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        plot.title = element_text(size=12, face="bold"),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Habitat")) -> p_dbrda_16S
p_dbrda_16S
ggsave("dbrda_16S.tiff", plot=p_dbrda_16S, width=5, height=5, dpi=300)

# ITS Distance Evaluation 

# dbRDA Plot Using Jaccard

df_ITS <- data.frame(sample_data(rare_ITS_habitat))
OTU_ITS <- t(as(rare_ITS_habitat@otu_table, "matrix"))
OTU_ITS <- as.data.frame(OTU_ITS)
dist_ITS <- vegdist(OTU_ITS, "jaccard")
set.seed(1221990)
rda_ITS <- dbrda(dist_ITS ~ habitat*irrigation*drought_tolerance, df_ITS)

anova(rda_ITS) # Overall model is significant
anova(rda_ITS, by = "margin", permutations = 9999) # Habitat significantly influences composition
anova(rda_ITS, by="terms", permutation = 9999) # p = 0.001
y <- varpart(dist_ITS, ~habitat, ~irrigation, ~drought_tolerance, data = df_ITS)
print(y)

# Measure for dispersion
dispersion_ITS <- betadisper(dist_ITS, df_ITS$habitat, type = c("centroid"), bias.adjust = TRUE, sqrt.dist = FALSE, add = FALSE)
anova(dispersion_ITS)

set.seed(1221990)
permutest(dispersion_ITS, pairwise = TRUE, permutations = 9999)
(dispersion_ITS_HSD <- TukeyHSD(dispersion_ITS))

# Plot Dispersion
dispersion_ITS$group <- factor(dispersion_ITS$group, levels = c("root_endosphere", "rhizosphere", "soil"))
data_ITS <- data.frame(Distance_to_centroid=dispersion_ITS$distances, Group=dispersion_ITS$group)
groups_ITS <- dispersion_ITS$group
level_order <- c("root_endosphere", "rhizosphere", "soil")
disper_ITS_plot <- ggplot(data=data_ITS, aes(x= factor(Group, level = level_order), y=Distance_to_centroid, color = groups_ITS)) + geom_boxplot(alpha=0.5) + geom_jitter(shape = 16, size =2, position = position_jitter(0.2)) + 
  scale_color_manual("Microbiome Compartment",
                     values = c("#FC8D62", "#8DA0CB", "#66C2A5"),
                     labels = c(expression(paste("Root Endosphere", sep="")),
                                expression(paste("Rhizosphere", sep="")),
                                expression(paste("Soil", sep="")))) +
  scale_x_discrete(name = NULL, labels = c("Root Endosphere", "Rhizosphere", "Soil")) +
  scale_y_continuous(name = "Distance to Centroid") +
  ggtitle("ITS Dispersion") +
  panel_border() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12))

disper_ITS_plot
# Plot Community Composition
plot.data <- merge(summary(rda_ITS)$sites, df_ITS, by = "row.names") 
plot.data$habitat <- factor(plot.data$habitat, levels = c("root_endosphere", "rhizosphere", "soil"))
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, color = habitat), 
             size = 2, alpha = 1) +
  ggtitle("ITS Composition") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_ITS)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_ITS)$cont)[2,2], 1)), "%)")) +
  theme_bw() +
  scale_color_manual("Microbiome Compartment",
                     values = c("#FC8D62", "#8DA0CB", "#66C2A5"),
                     labels = c(expression(paste("Root Endosphere", sep="")),
                                expression(paste("Rhizosphere", sep="")),
                                expression(paste("Soil", sep="")))) +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        plot.title = element_text(size=12, face="bold"),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Habitat")) -> p_dbrda_its
p_dbrda_its

# Append dbRDA and dispersion plots together

leg <- get_legend(p_dbrda_16S)
prow1 <- plot_grid(p_dbrda_16S + theme(legend.position = "none"), disper_16S_plot + theme(legend.position = "none"), labels = c("A","B"), ncol = 1, align = "v")
prow2 <- plot_grid(p_dbrda_its + theme(legend.position = "none"), disper_ITS_plot + theme(legend.position = "none"), labels = c("C","D"), ncol = 1, align = "v")
prow3 <- plot_grid(prow1, prow2, leg, ncol = 3, rel_widths = c(1, 1,0.75))
prow3 
ggsave(filename=paste("beta_diversity_habitat.tiff", sep=""), plot=prow3, width=7.5, height=6.5, dpi=600)
prow

# Core Microbiome ----------------------------------------------

# Bacteria/Archaea Root Endosphere

# subset by habitat
rare_16S_re_core <- subset_samples(rare_16S_habitat, habitat == "root_endosphere")

# Convert to relative abundance
rare_16S_re_core <- microbiome::transform(rare_16S_re_core, "compositional")



# Create a column with appended irrigation.drought tolerance 
rare_16S_re_core@sam_data$irrigation.drought_tolerance <- with(rare_16S_re_core@sam_data, 
                                                               interaction(irrigation, drought_tolerance))

# Make a list of treatments
treatments <- unique(as.character(meta(rare_16S_re_core)$irrigation.drought_tolerance))

# Write a for loop to go through each treatment one by one and combine identified core taxa into a list
list_core_re <- c() 

for (n in treatments){
  core_re_sub <- subset_samples(rare_16S_re_core, irrigation.drought_tolerance == n)
  
  core_m <- core_members(core_re_sub,
                         prevalence = 0.95,
                         detection = 0)
  
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core_re[[n]] <- core_m
}

# Export core taxa
core_re_taxa <- subset(tax_table(rare_16S_re_core), rownames(tax_table(rare_16S_re_core)) %in% c(core_m))
core_re_taxa <- as.data.frame(core_re_taxa)
distinct_family <- core_re_taxa$Family
distinct_family <- unique(distinct_family)
write.csv(distinct_family, file = "16S_core_re.csv")                                                                                                       

#Fungal Root Endosphere

# subset by habitat
rare_ITS_re_core <- subset_samples(rare_ITS_habitat, habitat == "root_endosphere")

# Convert to relative abundance
rare_ITS_re_core <- microbiome::transform(rare_ITS_re_core, "compositional")



# Create a column with appended irrigation.drought tolerance 
rare_ITS_re_core@sam_data$irrigation.drought_tolerance <- with(rare_ITS_re_core@sam_data, 
                                                               interaction(irrigation, drought_tolerance))

# Make a list of treatments
treatments <- unique(as.character(meta(rare_ITS_re_core)$irrigation.drought_tolerance))

# Write a for loop to go through each treatment one by one and combine identified core taxa into a list
list_core_re <- c() 

for (n in treatments){
  core_re_sub <- subset_samples(rare_ITS_re_core, irrigation.drought_tolerance == n)
  
  core_m <- core_members(core_re_sub,
                         prevalence = 0.95,
                         detection = 0)
  
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core_re[[n]] <- core_m
}

# No core taxa were found

# Archaeal/Bacterial Rhizosphere

# subset by habitat
rare_16S_rhizo_core <- subset_samples(rare_16S_habitat, habitat == "rhizosphere")

# Convert to relative abundance
rare_16S_rhizo_core <- microbiome::transform(rare_16S_rhizo_core, "compositional")



# Create a column with appended irrigation.drought tolerance 
rare_16S_rhizo_core@sam_data$irrigation.drought_tolerance <- with(rare_16S_rhizo_core@sam_data, 
                                                                  interaction(irrigation, drought_tolerance))

# Make a list of treatments
treatments <- unique(as.character(meta(rare_16S_rhizo_core)$irrigation.drought_tolerance))

# Write a for loop to go through each treatment one by one and combine identified core taxa into a list
list_core_rhizo <- c() 

for (n in treatments){
  core_rhizo_sub <- subset_samples(rare_16S_rhizo_core, irrigation.drought_tolerance == n)
  
  core_m <- core_members(core_rhizo_sub,
                         prevalence = 0.95,
                         detection = 0)
  
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core_rhizo[[n]] <- core_m
}

# Export core taxa
core_rhizo_taxa <- subset(tax_table(rare_16S_rhizo_core), rownames(tax_table(rare_16S_rhizo_core)) %in% c(core_m))
core_rhizo_taxa <- as.data.frame(core_rhizo_taxa)
distinct_genus <- core_rhizo_taxa$Family
distinct_genus <- unique(distinct_genus)
write.csv(distinct_genus, file = "16S_core_rhizo.csv")                                                                                                       

# Fungal Rhizosphere

# subset by habitat
rare_ITS_rhizo_core <- subset_samples(rare_ITS_habitat, habitat == "rhizosphere")

# Convert to relative abundance
rare_ITS_rhizo_core <- microbiome::transform(rare_ITS_rhizo_core, "compositional")



# Create a column with appended irrigation.drought tolerance 
rare_ITS_rhizo_core@sam_data$irrigation.drought_tolerance <- with(rare_ITS_rhizo_core@sam_data, 
                                                                  interaction(irrigation, drought_tolerance))

# Make a list of treatments
treatments <- unique(as.character(meta(rare_ITS_rhizo_core)$irrigation.drought_tolerance))

# Write a for loop to go through each treatment one by one and combine identified core taxa into a list
list_core_rhizo <- c() 

for (n in treatments){
  core_rhizo_sub <- subset_samples(rare_ITS_rhizo_core, irrigation.drought_tolerance == n)
  
  core_m <- core_members(core_rhizo_sub,
                         prevalence = .95,
                         detection = 0)
  
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core_rhizo[[n]] <- core_m
}

# Export core taxa
core_rhizo_taxa <- subset(tax_table(rare_ITS_rhizo_core), rownames(tax_table(rare_ITS_rhizo_core)) %in% c(core_m))
core_rhizo_taxa <- as.data.frame(core_rhizo_taxa)
distinct_genus <- core_rhizo_taxa$Family
distinct_genus <- unique(distinct_genus)
write.csv(distinct_genus, file = "ITS_core_rhizo.csv")    

# Archaeal/Bacterial Soil
# subset by habitat
rare_16S_soil_core <- subset_samples(rare_16S_habitat, habitat == "soil")

# Convert to relative abundance
rare_16S_soil_core <- microbiome::transform(rare_16S_soil_core, "compositional")



# Create a column with appended irrigation.drought tolerance 
rare_16S_soil_core@sam_data$irrigation.drought_tolerance <- with(rare_16S_soil_core@sam_data, 
                                                                 interaction(irrigation, drought_tolerance))

# Make a list of treatments
treatments <- unique(as.character(meta(rare_16S_soil_core)$irrigation.drought_tolerance))

# Write a for loop to go through each treatment one by one and combine identified core taxa into a list
list_core_soil <- c() 

for (n in treatments){
  core_soil_sub <- subset_samples(rare_16S_soil_core, irrigation.drought_tolerance == n)
  
  core_m <- core_members(core_soil_sub,
                         prevalence = .95,
                         detection = 0)
  
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core_soil[[n]] <- core_m
}

# Export core taxa
core_soil_taxa <- subset(tax_table(rare_16S_soil_core), rownames(tax_table(rare_16S_soil_core)) %in% c(core_m))
core_soil_taxa <- as.data.frame(core_soil_taxa)
distinct_genus <- core_soil_taxa$Family
distinct_genus <- unique(distinct_genus)
write.csv(distinct_genus, file = "16S_core_soil.csv") 

# Fungal Soil
# subset by habitat
rare_ITS_soil_core <- subset_samples(rare_ITS_habitat, habitat == "soil")

# Convert to relative abundance
rare_ITS_soil_core <- microbiome::transform(rare_ITS_soil_core, "compositional")



# Create a column with appended irrigation.drought tolerance 
rare_ITS_soil_core@sam_data$irrigation.drought_tolerance <- with(rare_ITS_soil_core@sam_data, 
                                                                 interaction(irrigation, drought_tolerance))

# Make a list of treatments
treatments <- unique(as.character(meta(rare_ITS_soil_core)$irrigation.drought_tolerance))

# Write a for loop to go through each treatment one by one and combine identified core taxa into a list
list_core_soil <- c() 

for (n in treatments){
  core_soil_sub <- subset_samples(rare_ITS_soil_core, irrigation.drought_tolerance == n)
  
  core_m <- core_members(core_soil_sub,
                         prevalence = 0.95,
                         detection = 0)
  
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core_soil[[n]] <- core_m
}

# Export core taxa
core_soil_taxa <- subset(tax_table(rare_ITS_soil_core), rownames(tax_table(rare_ITS_soil_core)) %in% c(core_m))
core_soil_taxa <- as.data.frame(core_soil_taxa)
distinct_genus <- core_soil_taxa$Genus
distinct_genus <- unique(distinct_genus)
write.csv(distinct_genus, file = "ITS_core_soil.csv")      

# Alpha Diversity - Treatment --------------------------------------------------

# Subset OTU tables by habitat (Soil, Rhizosphere, and Root Endosphere) -- Rarefy each library respectively
data_16S_RE <- subset_samples(data_16S_counts, habitat == "root_endosphere")
data_16S_rhizosphere <- subset_samples(data_16S_counts, habitat == "rhizosphere")
data_16S_soil <- subset_samples(data_16S_counts, habitat == "soil")

data_ITS_RE <- subset_samples(data_ITS_counts, habitat == "root_endosphere")
data_ITS_rhizosphere <- subset_samples(data_ITS_counts, habitat == "rhizosphere")
data_ITS_soil <- subset_samples(data_ITS_counts, habitat == "soil")

set.seed(01221990)
rare_16S_RE <- rarefy_even_depth(data_16S_RE, sample.size = 14000) #14000
ntaxa(rare_16S_RE)

set.seed(01221990)
rare_16S_rhizosphere <- rarefy_even_depth(data_16S_rhizosphere, sample.size = 25000) #25000
ntaxa(rare_16S_rhizosphere)


set.seed(01221990)
rare_16S_soil <- rarefy_even_depth(data_16S_soil, sample.size = 45000) # 45000
ntaxa(rare_16S_soil)


set.seed(01221990)
rare_ITS_RE <- rarefy_even_depth(data_ITS_RE, sample.size = 1000) 
ntaxa(rare_ITS_RE)

set.seed(01221990)
rare_ITS_rhizosphere <- rarefy_even_depth(data_ITS_rhizosphere, sample.size = 1300) #1300
ntaxa(rare_ITS_rhizosphere)


set.seed(01221990)
rare_ITS_soil <- rarefy_even_depth(data_ITS_soil, sample.size = 10000) #10000
ntaxa(rare_ITS_soil)

# 16S Alpha Diversity: Drought tolerance * irrigation interaction

# Soil
OTU_16S_soil <- t(as(rare_16S_soil@otu_table, "matrix"))
OTU_16S_soil <- as.data.frame(OTU_16S_soil)

rich_16S_soil <- hill_taxa(OTU_16S_soil, q = 0)
rich_16S_soil <- merge(data_16S_counts@sam_data, rich_16S_soil, by = 0)

shan_16S_soil <- hill_taxa(OTU_16S_soil, q = 1)
shan_16S_soil <- merge(data_16S_counts@sam_data, shan_16S_soil, by = 0)

simp_16S_soil <- hill_taxa(OTU_16S_soil, q=2)
simp_16S_soil <- merge(data_16S_counts@sam_data, simp_16S_soil, by=0)

# Soil Statistics 
model <- lme(y ~ irrigation + drought_tolerance + irrigation*drought_tolerance, random=~1|Row.names, data = rich_16S_soil, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) 

model <- lme(y ~ irrigation + drought_tolerance + irrigation*drought_tolerance, random=~1|Row.names, data = shan_16S_soil, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) 
model <- lme(y ~ irrigation + drought_tolerance + irrigation*drought_tolerance, random=~1|Row.names, data = simp_16S_soil, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) 

# Rhizosphere

OTU_16S_rhizosphere <- t(as(rare_16S_rhizosphere@otu_table, "matrix"))
OTU_16S_rhizosphere <- as.data.frame(OTU_16S_rhizosphere)

rich_16S_rhizosphere <- hill_taxa(OTU_16S_rhizosphere, q = 0)
rich_16S_rhizosphere <- merge(data_16S_counts@sam_data, rich_16S_rhizosphere, by = 0)

shan_16S_rhizosphere <- hill_taxa(OTU_16S_rhizosphere, q = 1)
shan_16S_rhizosphere <- merge(data_16S_counts@sam_data, shan_16S_rhizosphere, by = 0)

simp_16S_rhizosphere <- hill_taxa(OTU_16S_rhizosphere, q=2)
simp_16S_rhizosphere <- merge(data_16S_counts@sam_data, simp_16S_rhizosphere, by=0)

# Rhizosphere Statistics
model <- lme(y ~ irrigation + drought_tolerance + irrigation*drought_tolerance, random=~1|Row.names, data = rich_16S_rhizosphere, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) 

model <- lme(y ~ irrigation + drought_tolerance + irrigation*drought_tolerance, random=~1|Row.names, data = shan_16S_rhizosphere, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) 

model <- lme(y ~ irrigation + drought_tolerance + irrigation*drought_tolerance, random=~1|Row.names, data = simp_16S_rhizosphere, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) 


# Root Endosphere

OTU_16S_re <- t(as(rare_16S_RE@otu_table, "matrix"))
OTU_16S_re <- as.data.frame(OTU_16S_re)

rich_16S_re <- hill_taxa(OTU_16S_re, q = 0)
rich_16S_re <- merge(data_16S_counts@sam_data, rich_16S_re, by = 0)

shan_16S_re <- hill_taxa(OTU_16S_re, q = 1)
shan_16S_re <- merge(data_16S_counts@sam_data, shan_16S_re, by = 0)

simp_16S_re <- hill_taxa(OTU_16S_re, q=2)
simp_16S_re <- merge(data_16S_counts@sam_data, simp_16S_re, by=0)

# Root Endosphere Statistics 
model <- lme(y ~ irrigation + drought_tolerance + irrigation*drought_tolerance, random=~1|Row.names, data = rich_16S_re, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) 

model <- lme(y ~ irrigation + drought_tolerance + irrigation*drought_tolerance, random=~1|Row.names, data = shan_16S_re, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) 

model <- lme(y ~ irrigation + drought_tolerance + irrigation*drought_tolerance, random=~1|Row.names, data = simp_16S_re, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) 

# ITS Alpha Diversity: Drought tolerance * irrigation interaction

# Calculate hill numbers 

# Soil

OTU_ITS_soil <- t(as(rare_ITS_soil@otu_table, "matrix"))
OTU_ITS_soil <- as.data.frame(OTU_ITS_soil)

rich_ITS_soil <- hill_taxa(OTU_ITS_soil, q = 0)
rich_ITS_soil <- merge(data_ITS_counts@sam_data, rich_ITS_soil, by = 0)

shan_ITS_soil <- hill_taxa(OTU_ITS_soil, q = 1)
shan_ITS_soil <- merge(data_ITS_counts@sam_data, shan_ITS_soil, by = 0)

simp_ITS_soil <- hill_taxa(OTU_ITS_soil, q=2)
simp_ITS_soil <- merge(data_ITS_counts@sam_data, simp_ITS_soil, by=0)

# Soil Statistics 
model <- lme(y ~ irrigation + drought_tolerance + irrigation*drought_tolerance, random=~1|Row.names, data = rich_ITS_soil, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) 

model <- lme(y ~ irrigation + drought_tolerance + irrigation*drought_tolerance, random=~1|Row.names, data = shan_ITS_soil, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) 

model <- lme(y ~ irrigation + drought_tolerance + irrigation*drought_tolerance, random=~1|Row.names, data = simp_ITS_soil, method="REML")
anova.lme(model, type="sequential", adjustSigma = T)

# Rhizosphere

OTU_ITS_rhizosphere <- t(as(rare_ITS_rhizosphere@otu_table, "matrix"))
OTU_ITS_rhizosphere <- as.data.frame(OTU_ITS_rhizosphere)

rich_ITS_rhizosphere <- hill_taxa(OTU_ITS_rhizosphere, q = 0)
rich_ITS_rhizosphere <- merge(data_ITS_counts@sam_data, rich_ITS_rhizosphere, by = 0)

shan_ITS_rhizosphere <- hill_taxa(OTU_ITS_rhizosphere, q = 1)
shan_ITS_rhizosphere <- merge(data_ITS_counts@sam_data, shan_ITS_rhizosphere, by = 0)

simp_ITS_rhizosphere <- hill_taxa(OTU_ITS_rhizosphere, q=2)
simp_ITS_rhizosphere <- merge(data_ITS_counts@sam_data, simp_ITS_rhizosphere, by=0)

# Rhizosphere Statistics 
model <- lme(y ~ irrigation + drought_tolerance + irrigation*drought_tolerance, random=~1|Row.names, data = rich_ITS_rhizosphere, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) # No significant influence on rhizosphere alpha diversity

model <- lme(y ~ irrigation + drought_tolerance + irrigation*drought_tolerance, random=~1|Row.names, data = shan_ITS_rhizosphere, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) # No significant influence on rhizosphere alpha diversity 

model <- lme(y ~ irrigation + drought_tolerance + irrigation*drought_tolerance, random=~1|Row.names, data = simp_ITS_rhizosphere, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) # No significant influence on rhizosphere alpha diversity

# Root Endosphere

OTU_ITS_re <- t(as(rare_ITS_RE@otu_table, "matrix"))
OTU_ITS_re <- as.data.frame(OTU_ITS_re)

rich_ITS_re <- hill_taxa(OTU_ITS_re, q = 0)
rich_ITS_re <- merge(data_ITS_counts@sam_data, rich_ITS_re, by = 0)

shan_ITS_re <- hill_taxa(OTU_ITS_re, q = 1)
shan_ITS_re <- merge(data_ITS_counts@sam_data, shan_ITS_re, by = 0)

simp_ITS_re <- hill_taxa(OTU_ITS_re, q=2)
simp_ITS_re <- merge(data_ITS_counts@sam_data, simp_ITS_re, by=0)

# Root Endosphere 
model <- lme(y ~ irrigation + drought_tolerance + irrigation*drought_tolerance, random=~1|Row.names, data = rich_ITS_re, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) # No significant influence on richness 

model <- lme(y ~ irrigation + drought_tolerance + irrigation*drought_tolerance, random=~1|Row.names, data = shan_ITS_re, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) # No significant influence on shannon diversity 

model <- lme(y ~ irrigation + drought_tolerance + irrigation*drought_tolerance, random=~1|Row.names, data = simp_ITS_re, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) # No significant treatment influence on simpson diversity 

# Plot Shannon Diversity (Main Figure)
shan_16S_re$compartment <- "Root Endosphere"
shan_16S_rhizosphere$compartment <- "Rhizosphere"
shan_16S_soil$compartment <- "Bulk Soil"
shannon_16S <- rbind(shan_16S_re %>% dplyr::select(irrigation, y, drought_tolerance, marker_gene, compartment),
                     shan_16S_rhizosphere %>% dplyr::select(irrigation, y, drought_tolerance, marker_gene, compartment),
                     shan_16S_soil %>% dplyr::select(irrigation, y, drought_tolerance, marker_gene, compartment))

# Shannon Diversity Plot
#16S
shannon_16S %>% 
  ggplot(aes(x = irrigation, y = y, fill = drought_tolerance)) +
  geom_boxplot() +
  facet_grid(rows = vars(compartment), 
             scales = "free") +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) +
  scale_y_continuous(name = bquote(~phantom()^1*italic(D))) +
  scale_fill_manual("Drought Tolerance",
                    values = c("tomato1","steelblue4"),
                    labels = c(expression(paste("Drought Tolerant", sep="")),
                               expression(paste("Drought Susceptible", sep="")))) +
  panel_border() +
  theme_bw() +
  ggtitle(label = "Archaea/Bacterial Shannon Diversity") +
  theme(plot.title = element_text(size = 9.5)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        legend.title = element_text(size = 12)) +
  geom_vline(xintercept = c(1.5), linetype = "dotted") -> p_shannon_16S
p_shannon_16S

#ITS
shan_ITS_re$compartment <- "Root Endosphere"
shan_ITS_rhizosphere$compartment <- "Rhizosphere"
shan_ITS_soil$compartment <- "Bulk Soil"

shannon_ITS <- rbind(shan_ITS_re %>% dplyr::select(irrigation, y, drought_tolerance, marker_gene, compartment),
                     shan_ITS_rhizosphere %>% dplyr::select(irrigation, y, drought_tolerance, marker_gene, compartment),
                     shan_ITS_soil %>% dplyr::select(irrigation, y, drought_tolerance, marker_gene, compartment))

shannon_ITS %>% 
  ggplot(aes(x = irrigation, y = y, fill = drought_tolerance)) +
  geom_boxplot() +
  facet_grid(rows = vars(compartment), 
             scales = "free") +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) +
  scale_y_continuous(name = bquote(~phantom()^1*italic(D))) +
  scale_fill_manual("Drought Tolerance",
                    values = c("tomato1","steelblue4"),
                    labels = c(expression(paste("Drought Tolerant", sep="")),
                               expression(paste("Drought Susceptible", sep="")))) +
  panel_border() +
  theme_bw() +
  ggtitle(label = "Fungal Shannon Diversity") +
  theme(plot.title = element_text(size = 9.5)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        legend.title = element_text(size = 12)) +
  geom_vline(xintercept = c(1.5), linetype = "dotted") -> p_shannon_ITS
p_shannon_ITS

# Append plots together
leg <- get_legend(p_shannon_16S)
prow <- plot_grid(p_shannon_16S + theme(legend.position = "none"), p_shannon_ITS + theme(legend.position = "none"),
                  ncol = 2, labels = "AUTO", align = "vh")
prow <- plot_grid(prow, leg, ncol = 2, rel_widths = c(1, 0.25))
prow
ggsave(filename=paste("Fig_1.tiff",sep=""), plot=prow, width=7.2, height=7, dpi=600)

# Plot Species Richness (Supplementary)
rich_16S_re$compartment <- "Root Endosphere"
rich_16S_rhizosphere$compartment <- "Rhizosphere"
rich_16S_soil$compartment <- "Bulk Soil"
rich_16S <- rbind(rich_16S_re %>% dplyr::select(irrigation, y, drought_tolerance, marker_gene, compartment),
                  rich_16S_rhizosphere %>% dplyr::select(irrigation, y, drought_tolerance, marker_gene, compartment),
                  rich_16S_soil %>% dplyr::select(irrigation, y, drought_tolerance, marker_gene, compartment))

# Richness Diversity Plot
#16S
rich_16S %>% 
  ggplot(aes(x = irrigation, y = y, fill = drought_tolerance)) +
  geom_boxplot() +
  facet_grid(rows = vars(compartment), 
             scales = "free") +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) +
  scale_y_continuous(name = bquote(~phantom()^0*italic(D))) +
  scale_fill_manual("Drought Tolerance",
                    values = c("tomato1","steelblue4"),
                    labels = c(expression(paste("Drought Tolerant", sep="")),
                               expression(paste("Drought Susceptible", sep="")))) +
  panel_border() +
  theme_bw() +
  ggtitle(label = "Archaea/Bacterial Species Richness") +
  theme(plot.title = element_text(size = 9.5)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        legend.title = element_text(size = 12)) +
  geom_vline(xintercept = c(1.5), linetype = "dotted") -> p_rich_16S
p_rich_16S

#ITS
rich_ITS_re$compartment <- "Root Endosphere"
rich_ITS_rhizosphere$compartment <- "Rhizosphere"
rich_ITS_soil$compartment <- "Bulk Soil"

rich_ITS <- rbind(rich_ITS_re %>% dplyr::select(irrigation, y, drought_tolerance, marker_gene, compartment),
                  rich_ITS_rhizosphere %>% dplyr::select(irrigation, y, drought_tolerance, marker_gene, compartment),
                  rich_ITS_soil %>% dplyr::select(irrigation, y, drought_tolerance, marker_gene, compartment))

rich_ITS %>% 
  ggplot(aes(x = irrigation, y = y, fill = drought_tolerance)) +
  geom_boxplot() +
  facet_grid(rows = vars(compartment), 
             scales = "free") +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) +
  scale_y_continuous(name = bquote(~phantom()^0*italic(D))) +
  scale_fill_manual("Drought Tolerance",
                    values = c("tomato1","steelblue4"),
                    labels = c(expression(paste("Drought Tolerant", sep="")),
                               expression(paste("Drought Susceptible", sep="")))) +
  panel_border() +
  theme_bw() +
  ggtitle(label = "Fungal Species Richness") +
  theme(plot.title = element_text(size = 9.5)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        legend.title = element_text(size = 12)) +
  geom_vline(xintercept = c(1.5), linetype = "dotted") -> p_rich_ITS
p_rich_ITS

# Append plots together
leg <- get_legend(p_rich_16S)
prow <- plot_grid(p_rich_16S + theme(legend.position = "none"), p_rich_ITS + theme(legend.position = "none"),
                  ncol = 2, labels = "AUTO", align = "vh")
prow <- plot_grid(prow, leg, ncol = 2, rel_widths = c(1, 0.25))
prow
ggsave(filename=paste("Fig_Supp_richness.tiff",sep=""), plot=prow, width=7.2, height=7, dpi=600)

# Plot Simpson Diversity (Supplementary)
simp_16S_re$compartment <- "Root Endosphere"
simp_16S_rhizosphere$compartment <- "Rhizosphere"
simp_16S_soil$compartment <- "Bulk Soil"
simp_16S <- rbind(simp_16S_re %>% dplyr::select(irrigation, y, drought_tolerance, marker_gene, compartment),
                  simp_16S_rhizosphere %>% dplyr::select(irrigation, y, drought_tolerance, marker_gene, compartment),
                  simp_16S_soil %>% dplyr::select(irrigation, y, drought_tolerance, marker_gene, compartment))

# Simpson Diversity Plot
#16S
simp_16S %>% 
  ggplot(aes(x = irrigation, y = y, fill = drought_tolerance)) +
  geom_boxplot() +
  facet_grid(rows = vars(compartment), 
             scales = "free") +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) +
  scale_y_continuous(name = bquote(~phantom()^2*italic(D))) +
  scale_fill_manual("Drought Tolerance",
                    values = c("tomato1","steelblue4"),
                    labels = c(expression(paste("Drought Tolerant", sep="")),
                               expression(paste("Drought Susceptible", sep="")))) +
  panel_border() +
  theme_bw() +
  ggtitle(label = "Archaea/Bacterial Simpson Diversity") +
  theme(plot.title = element_text(size = 9.5)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        legend.title = element_text(size = 12)) +
  geom_vline(xintercept = c(1.5), linetype = "dotted") -> p_simp_16S
p_simp_16S

#ITS
simp_ITS_re$compartment <- "Root Endosphere"
simp_ITS_rhizosphere$compartment <- "Rhizosphere"
simp_ITS_soil$compartment <- "Bulk Soil"

simp_ITS <- rbind(simp_ITS_re %>% dplyr::select(irrigation, y, drought_tolerance, marker_gene, compartment),
                  simp_ITS_rhizosphere %>% dplyr::select(irrigation, y, drought_tolerance, marker_gene, compartment),
                  simp_ITS_soil %>% dplyr::select(irrigation, y, drought_tolerance, marker_gene, compartment))

simp_ITS %>% 
  ggplot(aes(x = irrigation, y = y, fill = drought_tolerance)) +
  geom_boxplot() +
  facet_grid(rows = vars(compartment), 
             scales = "free") +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) +
  scale_y_continuous(name = bquote(~phantom()^2*italic(D))) +
  scale_fill_manual("Drought Tolerance",
                    values = c("tomato1","steelblue4"),
                    labels = c(expression(paste("Drought Tolerant", sep="")),
                               expression(paste("Drought Susceptible", sep="")))) +
  panel_border() +
  theme_bw() +
  ggtitle(label = "Fungal Simpson Diversity") +
  theme(plot.title = element_text(size = 9.5)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        legend.title = element_text(size = 12)) +
  geom_vline(xintercept = c(1.5), linetype = "dotted") -> p_simp_ITS
p_simp_ITS

# Append plots together
leg <- get_legend(p_simp_16S)
prow <- plot_grid(p_simp_16S + theme(legend.position = "none"), p_simp_ITS + theme(legend.position = "none"),
                  ncol = 2, labels = "AUTO", align = "vh")
prow <- plot_grid(prow, leg, ncol = 2, rel_widths = c(1, 0.25))
prow
ggsave(filename=paste("Fig_Supp_simpson.tiff",sep=""), plot=prow, width=7.2, height=7, dpi=600)

# Beta Diversity - Treatment ---------------------------------------------------

# 16S Library -- Subsetted by Habitat
# Root Endosphere
df_16S_RE <- data.frame(rare_16S_RE@sam_data)
df_16S_RE$irrigation.drought_tolerance <- with(rare_16S_RE@sam_data, interaction(irrigation, drought_tolerance))
OTU_16S <- t(as(rare_16S_RE@otu_table, "matrix"))
OTU_16S <- as.data.frame(OTU_16S)
dist_16S <- vegdist(OTU_16S, "jaccard")
set.seed(1221990)
rda_16S_RE <- dbrda(dist_16S ~ irrigation + drought_tolerance + irrigation*drought_tolerance, df_16S_RE)

# permANOVA
anova(rda_16S_RE) 
print(anova(rda_16S_RE, by="terms", permutation = 9999))
plot(rda_16S_RE)
y <- varpart(dist_16S, ~irrigation, ~drought_tolerance, data = df_16S_RE)
print(y)

#Dispersal
dispersion_16S_RE <- betadisper(dist_16S, df_16S_RE$irrigation.drought_tolerance, type = c("centroid"), bias.adjust = TRUE, sqrt.dist = FALSE, add = FALSE)
anova(dispersion_16S_RE)
set.seed(1221990)
permutest(dispersion_16S_RE, pairwise = TRUE, permutations = 9999)
(dispersion_16S_RE_HSD <- TukeyHSD(dispersion_16S_RE))

# Rhizosphere
df_16S_rhizo <- data.frame(rare_16S_rhizosphere@sam_data)
df_16S_rhizo$irrigation.drought_tolerance <- with(rare_16S_rhizosphere@sam_data, interaction(irrigation, drought_tolerance))
OTU_16S <- t(as(rare_16S_rhizosphere@otu_table, "matrix"))
OTU_16S <- as.data.frame(OTU_16S)
dist_16S <- vegdist(OTU_16S, "jaccard")
set.seed(1221990)
rda_16S_rhizosphere <- dbrda(dist_16S ~ irrigation + drought_tolerance + irrigation*drought_tolerance, df_16S_rhizo)

#permANOVA
anova(rda_16S_rhizosphere) 
print(anova(rda_16S_rhizosphere, by="terms", permutation = 9999))
plot(rda_16S_rhizosphere)
y <- varpart(dist_16S, ~irrigation, ~drought_tolerance, data = df_16S_rhizo)
print(y)

#Dispersal
dispersion_16S_rhizo <- betadisper(dist_16S, df_16S_rhizo$irrigation.drought_tolerance, type = c("centroid"), bias.adjust = TRUE, sqrt.dist = FALSE, add = FALSE)
anova(dispersion_16S_rhizo)
set.seed(1221990)
permutest(dispersion_16S_rhizo, pairwise = TRUE, permutations = 9999)
(dispersion_16S_rhizo_HSD <- TukeyHSD(dispersion_16S_rhizo))

# Bulk Soil
df_16S_soil <- data.frame(rare_16S_soil@sam_data)
df_16S_soil$irrigation.drought_tolerance <- with(rare_16S_soil@sam_data, interaction(irrigation, drought_tolerance))
OTU_16S <- t(as(rare_16S_soil@otu_table, "matrix"))
OTU_16S <- as.data.frame(OTU_16S)
dist_16S <- vegdist(OTU_16S, "jaccard")
set.seed(1221990)
rda_16S_soil <- dbrda(dist_16S ~ irrigation + drought_tolerance + irrigation*drought_tolerance, df_16S_soil)

#permANOVA
anova(rda_16S_soil) 
print(anova(rda_16S_soil, by="terms", permutation = 9999))
plot(rda_16S_soil)
y <- varpart(dist_16S, ~irrigation, ~drought_tolerance, data = df_16S_soil)
print(y)

#Dispersal
dispersion_16S_soil <- betadisper(dist_16S, df_16S_soil$irrigation.drought_tolerance, type = c("centroid"), bias.adjust = TRUE, sqrt.dist = FALSE, add = FALSE)
anova(dispersion_16S_soil)
set.seed(1221990)
permutest(dispersion_16S_soil, pairwise = TRUE, permutations = 9999)
(dispersion_16S_soil_HSD <- TukeyHSD(dispersion_16S_soil))

# Plot Ordinations 
# Root Endosphere
data_RE <- as.data.frame(summary(rda_16S_RE)$biplot) %>% rownames_to_column()
data_RE <- data.frame(lapply(data_RE, function(x) {
  gsub("irrigation.L", "Irrigation", x)}))
data_RE <- data.frame(lapply(data_RE, function(x) {
  gsub("drought_tolerance.L", "Drought Tolerance", x)}))
data_RE <- data.frame(lapply(data_RE, function(x) {
  gsub("Irrigation:Drought Tolerance", "Irrigation*Drought Tolerance", x)}))
data_RE$dbRDA1 <- as.numeric(data_RE$dbRDA1)
data_RE$dbRDA2 <- as.numeric(data_RE$dbRDA2)
data_RE$dbRDA3 <- as.numeric(data_RE$dbRDA3)
data_RE <- data_RE[-c(2), ]

data_RE_arrows <- as.data.frame(summary(rda_16S_RE)$biplot)
data_RE_arrows <- data_RE_arrows[-c(2), ]

plot.data <- merge(summary(rda_16S_RE)$sites, df_16S_RE, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, shape = irrigation, fill = drought_tolerance), 
             size = 3.5, alpha = 1) +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  ggtitle("Archaeal/Bacterial Root Endosphere") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_16S_RE)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_16S_RE)$cont)[2,2], 1)), "%)")) +
  geom_segment(data = data_RE_arrows, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               size = 0.45, color = "black", arrow = arrow()) +
  geom_text(data = data_RE, aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) -> p_ordination_re
p_ordination_re

# Rhizosphere
data_rhizosphere <- as.data.frame(summary(rda_16S_rhizosphere)$biplot) %>% rownames_to_column()
data_rhizosphere <- data.frame(lapply(data_rhizosphere, function(x) {
  gsub("irrigation.L", "Irrigation", x)}))
data_rhizosphere <- data.frame(lapply(data_rhizosphere, function(x) {
  gsub("drought_tolerance.L", "Drought Tolerance", x)}))
data_rhizosphere <- data.frame(lapply(data_rhizosphere, function(x) {
  gsub("Irrigation:Drought Tolerance", "Irrigation*Drought Tolerance", x)}))
data_rhizosphere$dbRDA1 <- as.numeric(data_rhizosphere$dbRDA1)
data_rhizosphere$dbRDA2 <- as.numeric(data_rhizosphere$dbRDA2)
data_rhizosphere$dbRDA3 <- as.numeric(data_rhizosphere$dbRDA3)
data_rhizosphere <- data_rhizosphere[-c(2,3), ]

data_rhizo_arrows <- as.data.frame(summary(rda_16S_rhizosphere)$biplot)
data_rhizo_arrows <- data_rhizo_arrows[-c(2,3), ]

plot.data <- merge(summary(rda_16S_rhizosphere)$sites, df_16S_rhizo, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, shape = irrigation, fill = drought_tolerance), 
             size = 3.5, alpha = 1) +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  ggtitle("Archaeal/Bacterial Rhizosphere") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_16S_rhizosphere)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_16S_rhizosphere)$cont)[2,2], 1)), "%)")) +
  geom_segment(data = data_rhizo_arrows, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               size = 0.45, color = "black", arrow = arrow()) +
  geom_text(data = data_rhizosphere, aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) -> p_ordination_rhizo
p_ordination_rhizo

# Bulk Soil
data_soil <- as.data.frame(summary(rda_16S_soil)$biplot) %>% rownames_to_column()
data_soil <- data.frame(lapply(data_soil, function(x) {
  gsub("irrigation.L", "Irrigation", x)}))
data_soil <- data.frame(lapply(data_soil, function(x) {
  gsub("drought_tolerance.L", "Drought Tolerance", x)}))
data_soil <- data.frame(lapply(data_soil, function(x) {
  gsub("Irrigation:Drought Tolerance", "Irrigation*Drought Tolerance", x)}))
data_soil$dbRDA1 <- as.numeric(data_soil$dbRDA1)
data_soil$dbRDA2 <- as.numeric(data_soil$dbRDA2)
data_soil$dbRDA3 <- as.numeric(data_soil$dbRDA3)
data_soil <- data_soil[-c(2,3), ]

data_soil_arrows <- as.data.frame(summary(rda_16S_soil)$biplot)
data_soil_arrows <- data_soil_arrows[-c(2,3), ]

plot.data <- merge(summary(rda_16S_soil)$sites, df_16S_soil, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, shape = irrigation, fill = drought_tolerance), 
             size = 3.5, alpha = 1) +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  ggtitle("Archaeal/Bacterial Bulk Soil") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_16S_soil)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_16S_soil)$cont)[2,2], 1)), "%)")) +
  geom_segment(data = data_soil_arrows, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               size = 0.45, color = "black", arrow = arrow()) +
  geom_text(data = data_soil, aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) -> p_ordination_soil
p_ordination_soil

# Plot Dispersal
#Root Endosphere
data <- data.frame(Distance_to_centroid=dispersion_16S_RE$distances, Group=dispersion_16S_RE$group)
data$irrigation <- df_16S_RE$irrigation
data$drought_tolerance <- df_16S_RE$drought_tolerance

groups <- dispersion_16S_RE$group
disper_16S_re <- ggplot(data=data, aes(x = irrigation, y=Distance_to_centroid, shape = irrigation, fill = drought_tolerance)) + 
  geom_jitter(size =3.5, position = position_jitter(0.1)) +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) + 
  scale_y_continuous(name = "Distance to Centroid") +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  geom_vline(xintercept = c(1.5), linetype = "dotted") +
  ggtitle("Archaeal/Bacterial Root Endosphere Dispersion") +
  panel_border() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12))
#Rhizosphere
data <- data.frame(Distance_to_centroid=dispersion_16S_rhizo$distances, Group=dispersion_16S_rhizo$group)
data$irrigation <- df_16S_rhizo$irrigation
data$drought_tolerance <- df_16S_rhizo$drought_tolerance

groups <- dispersion_16S_rhizo$group
disper_16S_rhizo <- ggplot(data=data, aes(x = irrigation, y=Distance_to_centroid, shape = irrigation, fill = drought_tolerance)) + 
  geom_jitter(size =3.5, position = position_jitter(0.1)) +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) + 
  scale_y_continuous(name = "Distance to Centroid") +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  geom_vline(xintercept = c(1.5), linetype = "dotted") +
  ggtitle("Archaeal/Bacterial Rhizosphere Dispersion") +
  panel_border() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12))

# Bulk Soil
data <- data.frame(Distance_to_centroid=dispersion_16S_soil$distances, Group=dispersion_16S_soil$group)
data$irrigation <- df_16S_soil$irrigation
data$drought_tolerance <- df_16S_soil$drought_tolerance

groups <- dispersion_16S_soil$group
disper_16S_soil <- ggplot(data=data, aes(x = irrigation, y=Distance_to_centroid, shape = irrigation, fill = drought_tolerance)) + 
  geom_jitter(size =3.5, position = position_jitter(0.1)) +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) + 
  scale_y_continuous(name = "Distance to Centroid") +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  geom_vline(xintercept = c(1.5), linetype = "dotted") +
  ggtitle("Archaeal/Bacterial Bulk Soil Dispersion") +
  panel_border() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12))
# ITS
# Root Endosphere
df_ITS_RE <- data.frame(rare_ITS_RE@sam_data)
df_ITS_RE$irrigation.drought_tolerance <- with(rare_ITS_RE@sam_data, interaction(irrigation, drought_tolerance))
OTU_ITS <- t(as(rare_ITS_RE@otu_table, "matrix"))
OTU_ITS <- as.data.frame(OTU_ITS)
dist_ITS <- vegdist(OTU_ITS, "jaccard")
set.seed(1221990)
rda_ITS_RE <- dbrda(dist_ITS ~ irrigation + drought_tolerance + irrigation*drought_tolerance, df_ITS_RE)

#permANOVA
anova(rda_ITS_RE) 
print(anova(rda_ITS_RE, by="terms", permutation = 9999))
plot(rda_ITS_RE)
y <- varpart(dist_ITS, ~irrigation, ~drought_tolerance, data = df_ITS_RE)
print(y)

#Dispersal
dispersion_ITS_RE <- betadisper(dist_ITS, df_ITS_RE$irrigation.drought_tolerance, type = c("centroid"), bias.adjust = TRUE, sqrt.dist = FALSE, add = FALSE)
anova(dispersion_ITS_RE)
set.seed(1221990)
permutest(dispersion_ITS_RE, pairwise = TRUE, permutations = 9999)
(dispersion_ITS_RE_HSD <- TukeyHSD(dispersion_ITS_RE))

# Rhizosphere
df_ITS_rhizo <- data.frame(rare_ITS_rhizosphere@sam_data)
df_ITS_rhizo$irrigation.drought_tolerance <- with(rare_ITS_rhizosphere@sam_data, interaction(irrigation, drought_tolerance))
OTU_ITS <- t(as(rare_ITS_rhizosphere@otu_table, "matrix"))
OTU_ITS <- as.data.frame(OTU_ITS)
dist_ITS <- vegdist(OTU_ITS, "jaccard")
set.seed(1221990)
rda_ITS_rhizosphere <- dbrda(dist_ITS ~ irrigation + drought_tolerance + irrigation*drought_tolerance, df_ITS_rhizo)

#permANOVA
anova(rda_ITS_rhizosphere) 
print(anova(rda_ITS_rhizosphere, by="terms", permutation = 9999))
plot(rda_ITS_rhizosphere)
y <- varpart(dist_ITS, ~irrigation, ~drought_tolerance, data = df_ITS_rhizo)
print(y)

#Dispersal
dispersion_ITS_rhizo <- betadisper(dist_ITS, df_ITS_rhizo$irrigation.drought_tolerance, type = c("centroid"), bias.adjust = TRUE, sqrt.dist = FALSE, add = FALSE)
anova(dispersion_ITS_rhizo)
set.seed(1221990)
permutest(dispersion_ITS_rhizo, pairwise = TRUE, permutations = 9999)
(dispersion_ITS_rhizo_HSD <- TukeyHSD(dispersion_ITS_rhizo))

# Bulk Soil
df_ITS_soil <- data.frame(rare_ITS_soil@sam_data)
df_ITS_soil$irrigation.drought_tolerance <- with(rare_ITS_soil@sam_data, interaction(irrigation, drought_tolerance))
OTU_ITS <- t(as(rare_ITS_soil@otu_table, "matrix"))
OTU_ITS <- as.data.frame(OTU_ITS)
dist_ITS <- vegdist(OTU_ITS, "jaccard")
set.seed(1221990)
rda_ITS_soil <- dbrda(dist_ITS ~ irrigation + drought_tolerance + irrigation*drought_tolerance, df_ITS_soil)

#permANOVA
anova(rda_ITS_soil) 
print(anova(rda_ITS_soil, by="terms", permutation = 9999))
plot(rda_ITS_soil)
y <- varpart(dist_ITS, ~irrigation, ~drought_tolerance, data = df_ITS_soil)
print(y)

#Dispersal
dispersion_ITS_soil <- betadisper(dist_ITS, df_ITS_soil$irrigation.drought_tolerance, type = c("centroid"), bias.adjust = TRUE, sqrt.dist = FALSE, add = FALSE)
anova(dispersion_ITS_soil)
set.seed(1221990)
permutest(dispersion_ITS_soil, pairwise = TRUE, permutations = 9999)
(dispersion_ITS_soil_HSD <- TukeyHSD(dispersion_ITS_soil))

# Plot Ordinations 
# Root Endosphere
data_RE <- as.data.frame(summary(rda_ITS_RE)$biplot) %>% rownames_to_column()
data_RE <- data.frame(lapply(data_RE, function(x) {
  gsub("irrigation.L", "Irrigation", x)}))
data_RE <- data.frame(lapply(data_RE, function(x) {
  gsub("drought_tolerance.L", "Drought Tolerance", x)}))
data_RE <- data.frame(lapply(data_RE, function(x) {
  gsub("Irrigation:Drought Tolerance", "Irrigation*Drought Tolerance", x)}))
data_RE$dbRDA1 <- as.numeric(data_RE$dbRDA1)
data_RE$dbRDA2 <- as.numeric(data_RE$dbRDA2)
data_RE$dbRDA3 <- as.numeric(data_RE$dbRDA3)
data_RE <- data_RE[-c(1,3), ]

data_RE_arrows <- as.data.frame(summary(rda_ITS_RE)$biplot)
data_RE_arrows <- data_RE_arrows[-c(1,3), ]

plot.data <- merge(summary(rda_ITS_RE)$sites, df_ITS_RE, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, shape = irrigation, fill = drought_tolerance), 
             size = 3.5, alpha = 1) +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  ggtitle("Fungal Root Endosphere") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_ITS_RE)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_ITS_RE)$cont)[2,2], 1)), "%)")) +
  geom_segment(data = data_RE_arrows, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               size = 0.45, color = "black", arrow = arrow()) +
  geom_text(data = data_RE, aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) -> p_ordination_re_ITS
p_ordination_re_ITS

# Rhizosphere
data_rhizosphere <- as.data.frame(summary(rda_ITS_rhizosphere)$biplot) %>% rownames_to_column()
data_rhizosphere <- data.frame(lapply(data_rhizosphere, function(x) {
  gsub("irrigation.L", "Irrigation", x)}))
data_rhizosphere <- data.frame(lapply(data_rhizosphere, function(x) {
  gsub("drought_tolerance.L", "Drought Tolerance", x)}))
data_rhizosphere <- data.frame(lapply(data_rhizosphere, function(x) {
  gsub("Irrigation:Drought Tolerance", "Irrigation*Drought Tolerance", x)}))
data_rhizosphere$dbRDA1 <- as.numeric(data_rhizosphere$dbRDA1)
data_rhizosphere$dbRDA2 <- as.numeric(data_rhizosphere$dbRDA2)
data_rhizosphere$dbRDA3 <- as.numeric(data_rhizosphere$dbRDA3)
data_rhizosphere <- data_rhizosphere[-c(2,3), ]

data_rhizo_arrows <- as.data.frame(summary(rda_ITS_rhizosphere)$biplot)
data_rhizo_arrows <- data_rhizo_arrows[-c(2,3), ]


plot.data <- merge(summary(rda_ITS_rhizosphere)$sites, df_ITS_rhizo, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, shape = irrigation, fill = drought_tolerance), 
             size = 3.5, alpha = 1) +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  ggtitle("Fungal Rhizosphere") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_ITS_rhizosphere)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_ITS_rhizosphere)$cont)[2,2], 1)), "%)")) +
  geom_segment(data = data_rhizo_arrows, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               size = 0.45, color = "black", arrow = arrow()) +
  geom_text(data = data_rhizosphere, aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) -> p_ordination_rhizo_ITS
p_ordination_rhizo_ITS

# Bulk Soil
data_soil <- as.data.frame(summary(rda_ITS_soil)$biplot) %>% rownames_to_column()
data_soil <- data.frame(lapply(data_soil, function(x) {
  gsub("irrigation.L", "Irrigation", x)}))
data_soil <- data.frame(lapply(data_soil, function(x) {
  gsub("drought_tolerance.L", "Drought Tolerance", x)}))
data_soil <- data.frame(lapply(data_soil, function(x) {
  gsub("Irrigation:Drought Tolerance", "Irrigation*Drought Tolerance", x)}))
data_soil$dbRDA1 <- as.numeric(data_soil$dbRDA1)
data_soil$dbRDA2 <- as.numeric(data_soil$dbRDA2)
data_soil$dbRDA3 <- as.numeric(data_soil$dbRDA3)
data_soil <- data_soil[-c(2,3), ]

data_soil_arrows <- as.data.frame(summary(rda_ITS_soil)$biplot)
data_soil_arrows <- data_soil_arrows[-c(2,3), ]

plot.data <- merge(summary(rda_ITS_soil)$sites, df_ITS_soil, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, shape = irrigation, fill = drought_tolerance), 
             size = 3.5, alpha = 1) +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  ggtitle("Fungal Bulk Soil") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_ITS_soil)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_ITS_soil)$cont)[2,2], 1)), "%)")) +
  geom_segment(data = data_soil_arrows, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               size = 0.45, color = "black", arrow = arrow()) +
  geom_text(data = data_soil, aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) -> p_ordination_soil_ITS
p_ordination_soil_ITS

# Plot Dispersal
#Root Endosphere
data <- data.frame(Distance_to_centroid=dispersion_ITS_RE$distances, Group=dispersion_ITS_RE$group)
data$irrigation <- df_ITS_RE$irrigation
data$drought_tolerance <- df_ITS_RE$drought_tolerance

groups <- dispersion_ITS_RE$group
disper_ITS_re <- ggplot(data=data, aes(x = irrigation, y=Distance_to_centroid, shape = irrigation, fill = drought_tolerance)) + 
  geom_jitter(size =3.5, position = position_jitter(0.1)) +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) + 
  scale_y_continuous(name = "Distance to Centroid") +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  geom_vline(xintercept = c(1.5), linetype = "dotted") +
  ggtitle("Fungal Root Endosphere Dispersion") +
  panel_border() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12))

#Rhizosphere
data <- data.frame(Distance_to_centroid=dispersion_ITS_rhizo$distances, Group=dispersion_ITS_rhizo$group)
data$irrigation <- df_ITS_rhizo$irrigation
data$drought_tolerance <- df_ITS_rhizo$drought_tolerance

groups <- dispersion_ITS_rhizo$group
disper_ITS_rhizo <- ggplot(data=data, aes(x = irrigation, y=Distance_to_centroid, shape = irrigation, fill = drought_tolerance)) + 
  geom_jitter(size =3.5, position = position_jitter(0.1)) +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) + 
  scale_y_continuous(name = "Distance to Centroid") +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  geom_vline(xintercept = c(1.5), linetype = "dotted") +
  ggtitle("Fungal Rhizosphere Dispersion") +
  panel_border() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12))

# Bulk Soil
data <- data.frame(Distance_to_centroid=dispersion_ITS_soil$distances, Group=dispersion_ITS_soil$group)
data$irrigation <- df_ITS_soil$irrigation
data$drought_tolerance <- df_ITS_soil$drought_tolerance

groups <- dispersion_ITS_soil$group
disper_ITS_soil <- ggplot(data=data, aes(x = irrigation, y=Distance_to_centroid, shape = irrigation, fill = drought_tolerance)) + 
  geom_jitter(size =3.5, position = position_jitter(0.1)) +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) + 
  scale_y_continuous(name = "Distance to Centroid") +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  geom_vline(xintercept = c(1.5), linetype = "dotted") +
  ggtitle("Fungal Bulk Soil Dispersion") +
  panel_border() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12))

# Append all ordination plots together
prow <- plot_grid(p_ordination_re + theme(legend.position = "none"), p_ordination_re_ITS + theme(legend.position = "none"),
                  p_ordination_rhizo + theme(legend.position = "none"), p_ordination_rhizo_ITS + theme(legend.position = "none"),
                  p_ordination_soil + theme(legend.position="none"), p_ordination_soil_ITS + theme(legend.position = "none"),
                  ncol = 2, labels = "AUTO", align = "vh")
leg <- get_legend(p_ordination_re)
prow <- plot_grid(prow, leg, ncol = 2, rel_widths = c(0.8, 0.15))
prow
ggsave(filename=paste("Fig_2.tiff", sep=""), plot=prow, width=12.5, height=10.5, dpi=600)

# Append all dispersal plots together
prow <- plot_grid(disper_16S_re + theme(legend.position = "none"), disper_ITS_re + theme(legend.position = "none"),
                  disper_16S_rhizo + theme(legend.position = "none"), disper_ITS_rhizo + theme(legend.position = "none"),
                  disper_16S_soil + theme(legend.position="none"), disper_ITS_soil + theme(legend.position = "none"),
                  ncol = 2, labels = "AUTO", align = "vh")
leg <- get_legend(disper_16S_re)
prow <- plot_grid(prow, leg, ncol = 2, rel_widths = c(0.8, 0.15))
prow
ggsave(filename=paste("Supp_Fig_Dispersal.tiff", sep=""), plot=prow, width=12.5, height=10.5, dpi=600)

# Beta Diversity - Genotype ----------------------------------------------------
#16S
# Root Endosphere
df_16S_RE <- data.frame(rare_16S_RE@sam_data)
df_16S_RE$irrigation.drought_tolerance <- with(rare_16S_RE@sam_data, interaction(irrigation, drought_tolerance))
OTU_16S <- t(as(rare_16S_RE@otu_table, "matrix"))
OTU_16S <- as.data.frame(OTU_16S)
dist_16S <- vegdist(OTU_16S, "jaccard")
set.seed(1221990)
rda_16S_RE <- dbrda(dist_16S ~ irrigation + genotype + irrigation*genotype, df_16S_RE)

# permANOVA
anova(rda_16S_RE) 
print(anova(rda_16S_RE, by="terms", permutation = 9999))
plot(rda_16S_RE)
y <- varpart(dist_16S, ~irrigation, ~genotype, data = df_16S_RE)
print(y)

# Rhizosphere
df_16S_rhizo <- data.frame(rare_16S_rhizosphere@sam_data)
df_16S_rhizo$irrigation.drought_tolerance <- with(rare_16S_rhizosphere@sam_data, interaction(irrigation, drought_tolerance))
OTU_16S <- t(as(rare_16S_rhizosphere@otu_table, "matrix"))
OTU_16S <- as.data.frame(OTU_16S)
dist_16S <- vegdist(OTU_16S, "jaccard")
set.seed(1221990)
rda_16S_rhizosphere <- dbrda(dist_16S ~ irrigation + genotype + irrigation*genotype, df_16S_rhizo)

#permANOVA
anova(rda_16S_rhizosphere) 
print(anova(rda_16S_rhizosphere, by="terms", permutation = 9999))
plot(rda_16S_rhizosphere)
y <- varpart(dist_16S, ~irrigation, ~genotype, data = df_16S_rhizo)
print(y)

# Bulk Soil
df_16S_soil <- data.frame(rare_16S_soil@sam_data)
df_16S_soil$irrigation.drought_tolerance <- with(rare_16S_soil@sam_data, interaction(irrigation, drought_tolerance))
OTU_16S <- t(as(rare_16S_soil@otu_table, "matrix"))
OTU_16S <- as.data.frame(OTU_16S)
dist_16S <- vegdist(OTU_16S, "jaccard")
set.seed(1221990)
rda_16S_soil <- dbrda(dist_16S ~ irrigation + genotype + irrigation*genotype, df_16S_soil)

#permANOVA
anova(rda_16S_soil) 
print(anova(rda_16S_soil, by="terms", permutation = 9999))
plot(rda_16S_soil)
y <- varpart(dist_16S, ~irrigation, ~genotype, data = df_16S_soil)
print(y)

# Plot Ordinations 
n = 50
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

# Root Endosphere
data_RE <- as.data.frame(summary(rda_16S_RE)$biplot) %>% rownames_to_column()
data_RE <- data.frame(lapply(data_RE, function(x) {
  gsub("irrigation.L", "Irrigation", x)}))
data_RE$dbRDA1 <- as.numeric(data_RE$dbRDA1)
data_RE$dbRDA2 <- as.numeric(data_RE$dbRDA2)
data_RE$dbRDA3 <- as.numeric(data_RE$dbRDA3)
data_RE <- data_RE[-c(17:31), ]

data_RE_arrows <- as.data.frame(summary(rda_16S_RE)$biplot)
data_RE_arrows <- data_RE_arrows[-c(17:31), ]

plot.data <- merge(summary(rda_16S_RE)$sites, df_16S_RE, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, shape = irrigation, fill = genotype), 
             size = 3.5, alpha = 1) +
  scale_fill_manual(values = col_vector) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  ggtitle("Archaeal/Bacterial Root Endosphere") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_16S_RE)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_16S_RE)$cont)[2,2], 1)), "%)")) +
  geom_segment(data = data_RE_arrows, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               size = 0.45, color = "black", arrow = arrow()) +
  geom_text(data = data_RE, aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Genotype")) -> p_ordination_re
p_ordination_re

# Rhizosphere
data_rhizosphere <- as.data.frame(summary(rda_16S_rhizosphere)$biplot) %>% rownames_to_column()
data_rhizosphere <- data.frame(lapply(data_rhizosphere, function(x) {
  gsub("irrigation.L", "Irrigation", x)}))
data_rhizosphere$dbRDA1 <- as.numeric(data_rhizosphere$dbRDA1)
data_rhizosphere$dbRDA2 <- as.numeric(data_rhizosphere$dbRDA2)
data_rhizosphere$dbRDA3 <- as.numeric(data_rhizosphere$dbRDA3)
data_rhizosphere <- data_rhizosphere[-c(2:31), ]

data_rhizo_arrows <- as.data.frame(summary(rda_16S_rhizosphere)$biplot)
data_rhizo_arrows <- data_rhizo_arrows[-c(2:31), ]

plot.data <- merge(summary(rda_16S_rhizosphere)$sites, df_16S_rhizo, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, shape = irrigation, fill = genotype), 
             size = 3.5, alpha = 1) +
  scale_fill_manual(values = col_vector) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  ggtitle("Archaeal/Bacterial Rhizosphere") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_16S_rhizosphere)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_16S_rhizosphere)$cont)[2,2], 1)), "%)")) +
  geom_segment(data = data_rhizo_arrows, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               size = 0.45, color = "black", arrow = arrow()) +
  geom_text(data = data_rhizosphere, aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) -> p_ordination_rhizo
p_ordination_rhizo

# Bulk Soil
data_soil <- as.data.frame(summary(rda_16S_soil)$biplot) %>% rownames_to_column()
data_soil <- data.frame(lapply(data_soil, function(x) {
  gsub("irrigation.L", "Irrigation", x)}))
data_soil$dbRDA1 <- as.numeric(data_soil$dbRDA1)
data_soil$dbRDA2 <- as.numeric(data_soil$dbRDA2)
data_soil$dbRDA3 <- as.numeric(data_soil$dbRDA3)
data_soil <- data_soil[-c(2:31), ]

data_soil_arrows <- as.data.frame(summary(rda_16S_soil)$biplot)
data_soil_arrows <- data_soil_arrows[-c(2:31), ]

plot.data <- merge(summary(rda_16S_soil)$sites, df_16S_soil, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, shape = irrigation, fill = genotype), 
             size = 3.5, alpha = 1) +
  scale_fill_manual(values = col_vector) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  ggtitle("Archaeal/Bacterial Bulk Soil") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_16S_soil)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_16S_soil)$cont)[2,2], 1)), "%)")) +
  geom_segment(data = data_soil_arrows, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               size = 0.45, color = "black", arrow = arrow()) +
  geom_text(data = data_soil, aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) -> p_ordination_soil
p_ordination_soil

#ITS
# Root Endosphere
df_ITS_RE <- data.frame(rare_ITS_RE@sam_data)
df_ITS_RE$irrigation.drought_tolerance <- with(rare_ITS_RE@sam_data, interaction(irrigation, drought_tolerance))
OTU_ITS <- t(as(rare_ITS_RE@otu_table, "matrix"))
OTU_ITS <- as.data.frame(OTU_ITS)
dist_ITS <- vegdist(OTU_ITS, "jaccard")
set.seed(1221990)
rda_ITS_RE <- dbrda(dist_ITS ~ irrigation + genotype + irrigation*genotype, df_ITS_RE)

# permANOVA
anova(rda_ITS_RE) 
print(anova(rda_ITS_RE, by="terms", permutation = 9999))
plot(rda_ITS_RE)
y <- varpart(dist_ITS, ~irrigation, ~genotype, data = df_ITS_RE)
print(y)

# Rhizosphere
df_ITS_rhizo <- data.frame(rare_ITS_rhizosphere@sam_data)
df_ITS_rhizo$irrigation.drought_tolerance <- with(rare_ITS_rhizosphere@sam_data, interaction(irrigation, drought_tolerance))
OTU_ITS <- t(as(rare_ITS_rhizosphere@otu_table, "matrix"))
OTU_ITS <- as.data.frame(OTU_ITS)
dist_ITS <- vegdist(OTU_ITS, "jaccard")
set.seed(1221990)
rda_ITS_rhizosphere <- dbrda(dist_ITS ~ irrigation + genotype + irrigation*genotype, df_ITS_rhizo)

#permANOVA
anova(rda_ITS_rhizosphere) 
print(anova(rda_ITS_rhizosphere, by="terms", permutation = 9999))
plot(rda_ITS_rhizosphere)
y <- varpart(dist_ITS, ~irrigation, ~genotype, data = df_ITS_rhizo)
print(y)

# Bulk Soil
df_ITS_soil <- data.frame(rare_ITS_soil@sam_data)
df_ITS_soil$irrigation.drought_tolerance <- with(rare_ITS_soil@sam_data, interaction(irrigation, drought_tolerance))
OTU_ITS <- t(as(rare_ITS_soil@otu_table, "matrix"))
OTU_ITS <- as.data.frame(OTU_ITS)
dist_ITS <- vegdist(OTU_ITS, "jaccard")
set.seed(1221990)
rda_ITS_soil <- dbrda(dist_ITS ~ irrigation + genotype + irrigation*genotype, df_ITS_soil)

#permANOVA
anova(rda_ITS_soil) 
print(anova(rda_ITS_soil, by="terms", permutation = 9999))
plot(rda_ITS_soil)
y <- varpart(dist_ITS, ~irrigation, ~genotype, data = df_ITS_soil)
print(y)

# Plot Ordinations 
n = 50
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

# Root Endosphere
data_RE <- as.data.frame(summary(rda_ITS_RE)$biplot) %>% rownames_to_column()
data_RE <- data.frame(lapply(data_RE, function(x) {
  gsub("irrigation.L", "Irrigation", x)}))
data_RE$dbRDA1 <- as.numeric(data_RE$dbRDA1)
data_RE$dbRDA2 <- as.numeric(data_RE$dbRDA2)
data_RE$dbRDA3 <- as.numeric(data_RE$dbRDA3)
data_RE <- data_RE[-c(17:31), ]

data_RE_arrows <- as.data.frame(summary(rda_ITS_RE)$biplot)
data_RE_arrows <- data_RE_arrows[-c(17:31), ]

plot.data <- merge(summary(rda_ITS_RE)$sites, df_ITS_RE, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, shape = irrigation, fill = genotype), 
             size = 3.5, alpha = 1) +
  scale_fill_manual(values = col_vector) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  ggtitle("Fungal Root Endosphere") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_ITS_RE)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_ITS_RE)$cont)[2,2], 1)), "%)")) +
  geom_segment(data = data_RE_arrows, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               size = 0.45, color = "black", arrow = arrow()) +
  geom_text(data = data_RE, aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Genotype")) -> p_ordination_re_ITS
p_ordination_re_ITS

# Rhizosphere
data_rhizosphere <- as.data.frame(summary(rda_ITS_rhizosphere)$biplot) %>% rownames_to_column()
data_rhizosphere <- data.frame(lapply(data_rhizosphere, function(x) {
  gsub("irrigation.L", "Irrigation", x)}))
data_rhizosphere$dbRDA1 <- as.numeric(data_rhizosphere$dbRDA1)
data_rhizosphere$dbRDA2 <- as.numeric(data_rhizosphere$dbRDA2)
data_rhizosphere$dbRDA3 <- as.numeric(data_rhizosphere$dbRDA3)
data_rhizosphere <- data_rhizosphere[-c(2:31), ]

data_rhizo_arrows <- as.data.frame(summary(rda_ITS_rhizosphere)$biplot)
data_rhizo_arrows <- data_rhizo_arrows[-c(2:31), ]

plot.data <- merge(summary(rda_ITS_rhizosphere)$sites, df_ITS_rhizo, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, shape = irrigation, fill = genotype), 
             size = 3.5, alpha = 1) +
  scale_fill_manual(values = col_vector) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  ggtitle("Fungal Rhizosphere") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_ITS_rhizosphere)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_ITS_rhizosphere)$cont)[2,2], 1)), "%)")) +
  geom_segment(data = data_rhizo_arrows, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               size = 0.45, color = "black", arrow = arrow()) +
  geom_text(data = data_rhizosphere, aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) -> p_ordination_rhizo_ITS
p_ordination_rhizo_ITS

# Bulk Soil
data_soil <- as.data.frame(summary(rda_ITS_soil)$biplot) %>% rownames_to_column()
data_soil <- data.frame(lapply(data_soil, function(x) {
  gsub("irrigation.L", "Irrigation", x)}))
data_soil$dbRDA1 <- as.numeric(data_soil$dbRDA1)
data_soil$dbRDA2 <- as.numeric(data_soil$dbRDA2)
data_soil$dbRDA3 <- as.numeric(data_soil$dbRDA3)
data_soil <- data_soil[-c(2:31), ]

data_soil_arrows <- as.data.frame(summary(rda_ITS_soil)$biplot)
data_soil_arrows <- data_soil_arrows[-c(2:31), ]

plot.data <- merge(summary(rda_ITS_soil)$sites, df_ITS_soil, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, shape = irrigation, fill = genotype), 
             size = 3.5, alpha = 1) +
  scale_fill_manual(values = col_vector) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  ggtitle("Fungal Bulk Soil") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_ITS_soil)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_ITS_soil)$cont)[2,2], 1)), "%)")) +
  geom_segment(data = data_soil_arrows, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               size = 0.45, color = "black", arrow = arrow()) +
  geom_text(data = data_soil, aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) -> p_ordination_soil_ITS
p_ordination_soil_ITS

prow <- plot_grid(p_ordination_re + theme(legend.position = "none"), p_ordination_re_ITS + theme(legend.position = "none"),
                  p_ordination_rhizo + theme(legend.position = "none"), p_ordination_rhizo_ITS + theme(legend.position = "none"),
                  p_ordination_soil + theme(legend.position="none"), p_ordination_soil_ITS + theme(legend.position = "none"),
                  ncol = 2, labels = "AUTO", align = "vh")
leg <- get_legend(p_ordination_re)
prow <- plot_grid(prow, leg, ncol = 2, rel_widths = c(0.8, 0.15))
prow
ggsave(filename=paste("Fig_SX.tiff", sep=""), plot=prow, width=12.5, height=10.5, dpi=600)

# Treatment Taxonomy------------------------------------------------------------

# ITS

# Develop overview plot, subset by habitat; filter to most abundant genera
# Phylum 


data_ITS_rel <- transform_sample_counts(data_ITS_counts, function(x) x/sum(x))

abundant_ITS <- subset_taxa(data_ITS_rel, Phylum == "Ascomycota" |
                              Phylum == "Basidiomycota" |
                              Phylum == "Mortierellomycota" |
                              Phylum == "unidentified")

colnames(tax_table(abundant_ITS))[1] <- "group"
tax_table(abundant_ITS)[,"group"] <- tax_table(abundant_ITS)[,"Phylum"]
abundant_ITS_glom <- tax_glom(abundant_ITS, taxrank = "Phylum")
abundant_ITS_melt <- psmelt(abundant_ITS_glom)

sum2 <- abundant_ITS_melt %>% 
  group_by(irrigation, habitat, drought_tolerance, Phylum) %>%
  summarise(means = mean(Abundance))

#Everything that's left is considered rare  
rare <- sum2 %>%
  group_by(irrigation, habitat, drought_tolerance) %>%
  summarise(means = 1- sum(means)) %>%
  mutate(Phylum = "Rare Taxa")

#concatenate the datasets
sum2 = rbind(sum2, rare)

sum2$Phylum <- dplyr::recode(sum2$Phylum, unidentified = "Unidentified")
sum2$habitat<- dplyr::recode(sum2$habitat, soil = "Bulk Soil")
sum2$habitat<- dplyr::recode(sum2$habitat, root_endosphere = "Root Endosphere")
sum2$habitat<- dplyr::recode(sum2$habitat, rhizosphere = "Rhizosphere")

#Order habitat
sum2$habitat<- ordered(sum2$habitat, c("Root Endosphere", "Rhizosphere", "Bulk Soil"))

#order groups
sum2$Phylum <- ordered(sum2$Phylum, c("Ascomycota", "Basidiomycota", "Mortierellomycota", 
                                      "Unidentified", "Rare Taxa"))

rel_plots_ITS2 <- list()
for(i in levels(as.factor(sum2$habitat))){
  sub <- subset(sum2, habitat %in% i)
  p2 <- ggplot(sub, aes(x = drought_tolerance, y = means, fill = Phylum)) +
    geom_bar(position = "stack", stat = "identity", col = "black") +
    facet_wrap(~irrigation,
               labeller = labeller(irrigation = c("full" = "Full", "reduced" = "Reduced"))) +
    theme(legend.position = "right", legend.title = element_blank(), axis.line = element_blank()) +
    scale_fill_manual(values = c("firebrick4", "orange1", "midnightblue", "forestgreen", "mediumorchid4")) +
    scale_x_discrete(labels=c("Tolerant", "Susceptible")) +
    scale_y_continuous(name = NULL, breaks = NULL) +
    panel_border() +
    ggtitle(i) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(strip.text = element_text(size = 12, face = "italic"), axis.text.x = element_text(size = 10),
          axis.title = element_blank(),
          title = element_text(size = 12))
  rel_plots_ITS2[[i]] <- p2 + theme(legend.position = "none")
}

leg <- get_legend(p2 + theme(legend.position = "bottom"))
grid_rel_plots_ITS2 <- plot_grid(plotlist = rel_plots_ITS2, nrow =1, ncol= 3, align = "v")
rel_ITS_2 <- plot_grid(grid_rel_plots_ITS2, leg, ncol = 1, nrow =2, rel_heights = c(1,0.1))
rel_ITS_2
ggsave(filename=paste("phylum_ITS.tiff", sep=""), plot=rel_ITS_2, width=7.5, height=5, dpi=600)

# Class
abundant_ITS <- subset_taxa(data_ITS_rel, Class == "Agaricomycetes" |
                              Class == "Dothideomycetes" |
                              Class == "Eurotiomycetes" |
                              Class == "Leotiomycetes"|
                              Class == "Lecanoromycetes" |
                              Class == "Mortierellomycetes" |
                              Class == "Pezizomycetes" |
                              Class == "Sordariomycetes" |
                              Class == "Tremellomycetes" |
                              Class == "unidentified")

colnames(tax_table(abundant_ITS))[1] <- "group"
tax_table(abundant_ITS)[,"group"] <- tax_table(abundant_ITS)[,"Class"]
abundant_ITS_glom <- tax_glom(abundant_ITS, taxrank = "Class")
abundant_ITS_melt <- psmelt(abundant_ITS_glom)

sum2 <- abundant_ITS_melt %>% 
  group_by(irrigation, drought_tolerance, habitat, Class) %>%
  summarise(means = mean(Abundance))

#Everything that's left is considered rare  
rare <- sum2 %>%
  group_by(irrigation, drought_tolerance, habitat) %>%
  summarise(means = 1- sum(means)) %>%
  mutate(Class = "Rare Taxa")

#concatenate the datasets
sum2 = rbind(sum2, rare)

sum2$Class <- dplyr::recode(sum2$Class, unidentified = "Unidentified")
sum2$habitat<- dplyr::recode(sum2$habitat, soil = "Bulk Soil")
sum2$habitat<- dplyr::recode(sum2$habitat, root_endosphere = "Root Endosphere")
sum2$habitat<- dplyr::recode(sum2$habitat, rhizosphere = "Rhizosphere")

#Order habitat
sum2$habitat<- ordered(sum2$habitat, c("Root Endosphere", "Rhizosphere", "Bulk Soil"))

#order groups
sum2$Class <- ordered(sum2$Class, c("Agaricomycetes", "Dothideomycetes", "Eurotiomycetes", "Lecanoromycetes",
                                    "Leotiomycetes", "Mortierellomycetes", "Pezizomycetes", "Sordariomycetes",
                                    "Tremellomycetes", "Unidentified", "Rare Taxa"))

rel_plots_ITS2 <- list()
for(i in levels(as.factor(sum2$habitat))){
  sub <- subset(sum2, habitat %in% i)
  p2 <- ggplot(sub, aes(x = drought_tolerance, y = means, fill = Class)) +
    geom_bar(position = "stack", stat = "identity", col = "black") +
    facet_wrap(~irrigation,
               labeller = labeller(irrigation = c("full" = "Full", "reduced" = "Reduced"))) +
    theme(legend.position = "right", legend.title = element_blank(), axis.line = element_blank()) +
    scale_fill_manual(values = c("firebrick4","tomato", "darkorange3", "orange1" ,"goldenrod4", "khaki1", "forestgreen",
                                 "olivedrab1", "midnightblue", "powderblue", "mediumorchid4")) +
    scale_x_discrete(labels=c("Tolerant", "Susceptible")) +
    scale_y_continuous(name = NULL, breaks = NULL) +
    panel_border() +
    ggtitle(i) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(strip.text = element_text(size = 12, face = "italic"), axis.text.x = element_text(size = 10),
          axis.title = element_blank(),
          title = element_text(size = 12))
  rel_plots_ITS2[[i]] <- p2 + theme(legend.position = "none")
}

leg <- get_legend(p2 + theme(legend.position = "bottom"))
grid_rel_plots_ITS2 <- plot_grid(plotlist = rel_plots_ITS2, nrow =1, ncol= 3, align = "v")
rel_ITS_2 <- plot_grid(grid_rel_plots_ITS2, leg, ncol = 1, nrow =2, rel_heights = c(1,0.2))
rel_ITS_2
ggsave(filename=paste("class_ITS.tiff", sep=""), plot=rel_ITS_2, width=7.5, height=5, dpi=600)

# Order
order_sum_ITS <- tapply(taxa_sums(data_ITS_rel), 
                        tax_table(data_ITS_rel)[, "Order"], sum, na.rm=TRUE)
top10order_ITS <- names(sort(order_sum_ITS, TRUE))[1:11]

abundant_ITS_order <- subset_taxa(data_ITS_rel, Order %in% top10order_ITS)
abundant_ITS_order_glom <- tax_glom(abundant_ITS_order, taxrank = "Order")
abundant_ITS_order_melt <- psmelt(abundant_ITS_order_glom)

sum2 <- abundant_ITS_order_melt %>% 
  filter(Order != "unidentified") %>%
  group_by(habitat, Order, irrigation, drought_tolerance) %>%
  summarise(means = mean(Abundance))

#Everything that's left is considered rare  
rare <- sum2 %>%
  group_by(habitat, irrigation, drought_tolerance) %>%
  summarise(means = 1- sum(means)) %>%
  mutate(Order = "Other/Unidentified Order")

#concatenate the datasets
sum2 = rbind(sum2, rare)

#order groups
sum2$Order <- forcats::fct_relevel(sum2$Order, "Other/Unidentified Order", after = Inf)

sum2$habitat<- dplyr::recode(sum2$habitat, soil = "Bulk Soil")
sum2$habitat<- dplyr::recode(sum2$habitat, root_endosphere = "Root Endosphere")
sum2$habitat<- dplyr::recode(sum2$habitat, rhizosphere = "Rhizosphere")

#Order habitat
sum2$habitat<- ordered(sum2$habitat, c("Root Endosphere", "Rhizosphere", "Bulk Soil"))

rel_plots_ITS2 <- list()
for(i in levels(as.factor(sum2$habitat))){
  sub <- subset(sum2, habitat %in% i)
  p2 <- ggplot(sub, aes(x = drought_tolerance, y = means, fill = Order)) +
    geom_bar(position = "stack", stat = "identity", col = "black") +
    facet_wrap(~irrigation,
               labeller = labeller(irrigation = c("full" = "Full", "reduced" = "Reduced"))) +
    theme(legend.position = "right", legend.title = element_blank(), axis.line = element_blank()) +
    scale_fill_manual(values = c("firebrick4","tomato", "darkorange3", "orange1" ,"goldenrod4", "khaki1", "forestgreen",
                                 "olivedrab1", "midnightblue", "powderblue", "mediumorchid4")) +
    scale_x_discrete(labels=c("Tolerant", "Susceptible")) +
    scale_y_continuous(name = NULL, breaks = NULL) +
    panel_border() +
    ggtitle(i) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(strip.text = element_text(size = 12, face = "italic"), axis.text.x = element_text(size = 10),
          axis.title = element_blank(),
          title = element_text(size = 12))
  rel_plots_ITS2[[i]] <- p2 + theme(legend.position = "none")
}

leg <- get_legend(p2 + theme(legend.position = "bottom"))
grid_rel_plots_ITS2 <- plot_grid(plotlist = rel_plots_ITS2, nrow =1, ncol= 3, align = "v")
rel_ITS_2 <- plot_grid(grid_rel_plots_ITS2, leg, ncol = 1, nrow =2, rel_heights = c(1,0.2))
rel_ITS_2
ggsave(filename=paste("order_ITS.tiff", sep=""), plot=rel_ITS_2, width=7.5, height=5, dpi=600)

# Family
family_sum_ITS <- tapply(taxa_sums(data_ITS_rel), 
                         tax_table(data_ITS_rel)[, "Family"], sum, na.rm=TRUE)
top10family_ITS <- names(sort(family_sum_ITS, TRUE))[1:11]

abundant_ITS_family <- subset_taxa(data_ITS_rel, Family %in% top10family_ITS)
abundant_ITS_family_glom <- tax_glom(abundant_ITS_family, taxrank = "Family")
abundant_ITS_family_melt <- psmelt(abundant_ITS_family_glom)

sum2 <- abundant_ITS_family_melt %>% 
  filter(Family != "unidentified") %>%
  group_by(habitat, Family, irrigation, drought_tolerance) %>%
  summarise(means = mean(Abundance))

#Everything that's left is considered rare  
rare <- sum2 %>%
  group_by(habitat, irrigation, drought_tolerance) %>%
  summarise(means = 1- sum(means)) %>%
  mutate(Family = "Other/Unidentified Family")

#concatenate the datasets
sum2 = rbind(sum2, rare)

#order groups
sum2$Family <- forcats::fct_relevel(sum2$Family, "Other/Unidentified Family", after = Inf)

sum2$habitat<- dplyr::recode(sum2$habitat, soil = "Bulk Soil")
sum2$habitat<- dplyr::recode(sum2$habitat, root_endosphere = "Root Endosphere")
sum2$habitat<- dplyr::recode(sum2$habitat, rhizosphere = "Rhizosphere")

#Order habitat
sum2$habitat<- ordered(sum2$habitat, c("Root Endosphere", "Rhizosphere", "Bulk Soil"))

rel_plots_ITS2 <- list()
for(i in levels(as.factor(sum2$habitat))){
  sub <- subset(sum2, habitat %in% i)
  p2 <- ggplot(sub, aes(x = drought_tolerance, y = means, fill = Family)) +
    geom_bar(position = "stack", stat = "identity", col = "black") +
    facet_wrap(~irrigation,
               labeller = labeller(irrigation = c("full" = "Full", "reduced" = "Reduced"))) +
    theme(legend.position = "right", legend.title = element_blank(), axis.line = element_blank()) +
    scale_fill_manual(values = c("firebrick4","tomato", "darkorange3", "orange1" ,"goldenrod4", "khaki1", "forestgreen",
                                 "olivedrab1", "midnightblue", "powderblue", "mediumorchid4")) +
    scale_x_discrete(labels=c("Tolerant", "Susceptible")) +
    scale_y_continuous(name = NULL, breaks = NULL) +
    panel_border() +
    ggtitle(i) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(strip.text = element_text(size = 12, face = "italic"), axis.text.x = element_text(size = 10),
          axis.title = element_blank(),
          title = element_text(size = 12))
  rel_plots_ITS2[[i]] <- p2 + theme(legend.position = "none")
}

leg <- get_legend(p2 + theme(legend.position = "bottom"))
grid_rel_plots_ITS2 <- plot_grid(plotlist = rel_plots_ITS2, nrow =1, ncol= 3, align = "v")
rel_ITS_2 <- plot_grid(grid_rel_plots_ITS2, leg, ncol = 1, nrow =2, rel_heights = c(1,0.2))
rel_ITS_2
ggsave(filename=paste("family_ITS.tiff", sep=""), plot=rel_ITS_2, width=7.5, height=5, dpi=600)


# Genus 
genus_sum_ITS <- tapply(taxa_sums(data_ITS_rel), 
                        tax_table(data_ITS_rel)[, "Genus"], sum, na.rm=TRUE)
top10genus_ITS <- names(sort(genus_sum_ITS, TRUE))[1:11]

abundant_ITS_genus <- subset_taxa(data_ITS_rel, Genus %in% top10genus_ITS)
abundant_ITS_genus_glom <- tax_glom(abundant_ITS_genus, taxrank = "Genus")
abundant_ITS_genus_melt <- psmelt(abundant_ITS_genus_glom)

sum2 <- abundant_ITS_genus_melt %>% 
  filter(Genus != "unidentified") %>%
  group_by(habitat, Genus, irrigation, drought_tolerance) %>%
  summarise(means = mean(Abundance))

#Everything that's left is considered rare  
rare <- sum2 %>%
  group_by(habitat, irrigation, drought_tolerance) %>%
  summarise(means = 1- sum(means)) %>%
  mutate(Genus = "Other/Unidentified Genus")

#concatenate the datasets
sum2 = rbind(sum2, rare)

#order groups
sum2$Genus <- forcats::fct_relevel(sum2$Genus, "Other/Unidentified Genus", after = Inf)

sum2$habitat<- dplyr::recode(sum2$habitat, soil = "Bulk Soil")
sum2$habitat<- dplyr::recode(sum2$habitat, root_endosphere = "Root Endosphere")
sum2$habitat<- dplyr::recode(sum2$habitat, rhizosphere = "Rhizosphere")

#Order habitat
sum2$habitat<- ordered(sum2$habitat, c("Root Endosphere", "Rhizosphere", "Bulk Soil"))

rel_plots_ITS2 <- list()
for(i in levels(as.factor(sum2$habitat))){
  sub <- subset(sum2, habitat %in% i)
  p2 <- ggplot(sub, aes(x = drought_tolerance, y = means, fill = Genus)) +
    geom_bar(position = "stack", stat = "identity", col = "black") +
    facet_wrap(~irrigation,
               labeller = labeller(irrigation = c("full" = "Full", "reduced" = "Reduced"))) +
    theme(legend.position = "right", legend.title = element_blank(), axis.line = element_blank()) +
    scale_fill_manual(values = c("firebrick4","tomato", "darkorange3", "orange1" ,"goldenrod4", "khaki1", "forestgreen",
                                 "olivedrab1", "midnightblue", "powderblue", "mediumorchid4")) +
    scale_x_discrete(labels=c("Tolerant", "Susceptible")) +
    scale_y_continuous(name = NULL, breaks = NULL) +
    panel_border() +
    ggtitle(i) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(strip.text = element_text(size = 12, face = "italic"), axis.text.x = element_text(size = 10),
          axis.title = element_blank(),
          title = element_text(size = 12))
  rel_plots_ITS2[[i]] <- p2 + theme(legend.position = "none")
}

leg <- get_legend(p2 + theme(legend.position = "bottom"))
grid_rel_plots_ITS2 <- plot_grid(plotlist = rel_plots_ITS2, nrow =1, ncol= 3, align = "v")
rel_ITS_2 <- plot_grid(grid_rel_plots_ITS2, leg, ncol = 1, nrow =2, rel_heights = c(1,0.2))
rel_ITS_2
ggsave(filename=paste("genus_ITS.tiff", sep=""), plot=rel_ITS_2, width=7.5, height=5, dpi=600)

# 16S 

data_16S_rel <- transform_sample_counts(data_16S_counts, function(x) x/sum(x))

#Seperate Proteos, assign grouping by Class, glom, and melt
proteos <- subset_taxa(data_16S_rel, Class == "Alphaproteobacteria" |
                         Class == "Deltaproteobacteria" |
                         Class == "Gammaproteobacteria")
colnames(tax_table(proteos))[1] <- "group"
tax_table(proteos)[,"group"] <- tax_table(proteos)[,"Class"]
proteos_glom <- tax_glom(proteos, taxrank = "group")
proteos_melt <- psmelt(proteos_glom)

#Do same for Archaea

archaea <- subset_taxa(data_16S_rel, Kingdom == "Archaea")
colnames(tax_table(archaea))[1] <- "group"
archaea_glom <- tax_glom(archaea, taxrank = "group")
archaea_melt <- psmelt(archaea_glom)

abundant_non_proteos <- subset_taxa(data_16S_rel, Phylum == "Acidobacteria" |
                                      Phylum == "Actinobacteria" |
                                      Phylum == "Bacteroidetes" |
                                      Phylum == "Chloroflexi" |
                                      Phylum == "Firmicutes" |
                                      Phylum == "Patescibacteria" |
                                      Phylum == "Gemmatimonadetes" |
                                      Phylum == "Rokubacteria" |
                                      Phylum == "Verrucomicrobia")

colnames(tax_table(abundant_non_proteos))[1] <- "group"
tax_table(abundant_non_proteos)[,"group"] <- tax_table(abundant_non_proteos)[,"Phylum"]
abundant_non_proteos_glom <- tax_glom(abundant_non_proteos, taxrank = "group")
abundant_non_proteos_melt <- psmelt(abundant_non_proteos_glom)

#concatenate the datasets
data_16S_grouping <- rbind(proteos_melt, abundant_non_proteos_melt, archaea_melt)

# Now split by species

sum2 <- data_16S_grouping %>% 
  group_by(irrigation, habitat, group, drought_tolerance) %>%
  summarise(means = mean(Abundance))

rare <- sum2 %>%
  group_by(irrigation, habitat, drought_tolerance) %>%
  summarise(means = 1- sum(means)) %>%
  mutate(group = "Rare Taxa")

sum2 = rbind(sum2, rare)
sum2$habitat<- dplyr::recode(sum2$habitat, soil = "Bulk Soil")
sum2$habitat<- dplyr::recode(sum2$habitat, root_endosphere = "Root Endosphere")
sum2$habitat<- dplyr::recode(sum2$habitat, rhizosphere = "Rhizosphere")

#Order habitat
sum2$habitat<- ordered(sum2$habitat, c("Root Endosphere", "Rhizosphere", "Bulk Soil"))

sum2$group <- ordered(sum2$group, c("Archaea", "Alphaproteobacteria", "Deltaproteobacteria", 
                                    "Gammaproteobacteria", "Acidobacteria",
                                    "Actinobacteria", "Bacteroidetes", "Chloroflexi", "Firmicutes",
                                    "Patescibacteria", "Gemmatimonadetes", "Rokubacteria", "Verrucomicrobia",
                                    "Rare Taxa"))

rel_plots_16S2 <- list()
for(i in levels(as.factor(sum2$habitat))){
  sub <- subset(sum2, habitat %in% i)
  p2 <- ggplot(sub, aes(x = drought_tolerance, y = means, fill = group)) +
    geom_bar(position = "stack", stat = "identity", col = "black") +
    facet_wrap(~irrigation,
               labeller = labeller(irrigation = c("full" = "Full", "reduced" = "Reduced"))) +
    theme(legend.position = "right", legend.title = element_blank(), axis.line = element_blank()) +
    scale_fill_manual(values = c("firebrick4","tomato", "darkorange3", "orange1" ,"goldenrod4", "khaki1", "forestgreen",
                                 "olivedrab1", "midnightblue", "powderblue","hotpink4","hotpink", "mediumvioletred", "mediumorchid4")) +
    scale_x_discrete(labels=c("Tolerant", "Susceptible")) +
    scale_y_continuous(name = NULL, breaks = NULL) +
    panel_border() +
    ggtitle(i) +
    guides(fill=guide_legend(title="Phylum")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(strip.text = element_text(size = 12, face = "italic"), axis.text.x = element_text(size = 10),
          axis.title = element_blank(),
          title = element_text(size = 12))
  rel_plots_16S2[[i]] <- p2 + theme(legend.position = "none")
}

leg <- get_legend(p2 + theme(legend.position = "bottom"))
grid_rel_plots_16S2 <- plot_grid(plotlist = rel_plots_16S2, nrow =1, ncol= 3, align = "v")
rel_16S_2 <- plot_grid(grid_rel_plots_16S2, leg, ncol = 1, nrow =2, rel_heights = c(1,0.2))
rel_16S_2
ggsave(filename=paste("phylum_16S.tiff", sep=""), plot=rel_16S_2, width=7.85, height=5, dpi=600)

# Order
order_sum_16S <- tapply(taxa_sums(data_16S_rel), 
                        tax_table(data_16S_rel)[, "Order"], sum, na.rm=TRUE)
top10order_16S <- names(sort(order_sum_16S, TRUE))[1:10]

abundant_16S_order <- subset_taxa(data_16S_rel, Order %in% top10order_16S)
abundant_16S_order_glom <- tax_glom(abundant_16S_order, taxrank = "Order")
abundant_16S_order_melt <- psmelt(abundant_16S_order_glom)

sum2 <- abundant_16S_order_melt %>% 
  group_by(habitat, Order, irrigation, drought_tolerance) %>%
  summarise(means = mean(Abundance))

#Everything that's left is considered rare  
rare <- sum2 %>%
  group_by(habitat, irrigation, drought_tolerance) %>%
  summarise(means = 1- sum(means)) %>%
  mutate(Order = "Other/Unidentified Order")

#concatenate the datasets
sum2 = rbind(sum2, rare)

#order groups
sum2$Order <- forcats::fct_relevel(sum2$Order, "Other/Unidentified Order", after = Inf)

#order habitat
sum2$habitat<- dplyr::recode(sum2$habitat, soil = "Bulk Soil")
sum2$habitat<- dplyr::recode(sum2$habitat, root_endosphere = "Root Endosphere")
sum2$habitat<- dplyr::recode(sum2$habitat, rhizosphere = "Rhizosphere")
sum2$habitat<- ordered(sum2$habitat, c("Root Endosphere", "Rhizosphere", "Bulk Soil"))

rel_plots_16S2 <- list()
for(i in levels(as.factor(sum2$habitat))){
  sub <- subset(sum2, habitat %in% i)
  p2 <- ggplot(sub, aes(x = drought_tolerance, y = means, fill = Order)) +
    geom_bar(position = "stack", stat = "identity", col = "black") +
    facet_wrap(~irrigation,
               labeller = labeller(irrigation = c("full" = "Full", "reduced" = "Reduced"))) +
    theme(legend.position = "right", legend.title = element_blank(), axis.line = element_blank()) +
    scale_fill_manual(values = c("firebrick4","tomato", "darkorange3", "orange1" ,"goldenrod4", "khaki1", "forestgreen",
                                 "olivedrab1", "midnightblue","powderblue", "mediumorchid4")) +
    scale_x_discrete(labels=c("Tolerant", "Susceptible")) +
    scale_y_continuous(name = NULL, breaks = NULL) +
    panel_border() +
    ggtitle(i) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(strip.text = element_text(size = 12, face = "italic"), axis.text.x = element_text(size = 10),
          axis.title = element_blank(),
          title = element_text(size = 12))
  rel_plots_16S2[[i]] <- p2 + theme(legend.position = "none")
}
leg <- get_legend(p2 + theme(legend.position = "bottom"))
grid_rel_plots_16S2 <- plot_grid(plotlist = rel_plots_16S2, nrow =1, ncol= 3, align = "v")
rel_16S_2 <- plot_grid(grid_rel_plots_16S2, leg, ncol = 1, nrow =2, rel_heights = c(1,0.2))
rel_16S_2
ggsave(filename=paste("order_16S.tiff", sep=""), plot=rel_16S_2, width=7.5, height=5, dpi=600)


# Family
family_sum_16S <- tapply(taxa_sums(data_16S_rel), 
                         tax_table(data_16S_rel)[, "Family"], sum, na.rm=TRUE)
top10family_16S <- names(sort(family_sum_16S, TRUE))[1:11]

abundant_16S_family <- subset_taxa(data_16S_rel, Family %in% top10family_16S)
abundant_16S_family_glom <- tax_glom(abundant_16S_family, taxrank = "Family")
abundant_16S_family_melt <- psmelt(abundant_16S_family_glom)

sum2 <- abundant_16S_family_melt %>% 
  group_by(habitat, Family, irrigation, drought_tolerance) %>%
  summarise(means = mean(Abundance))

sum2$Family <- dplyr::recode(sum2$Family, "uncultured bacterium" = "Other/Unidentified Family")

#Everything that's left is considered rare  
rare <- sum2 %>%
  group_by(habitat, irrigation, drought_tolerance) %>%
  summarise(means = 1- sum(means)) %>%
  mutate(Family = "Other/Unidentified Family")


#concatenate the datasets
sum2 = rbind(sum2, rare)

#family groups
sum2$Family <- forcats::fct_relevel(sum2$Family, "Other/Unidentified Family", after = Inf)

#family habitat
sum2$habitat<- dplyr::recode(sum2$habitat, soil = "Bulk Soil")
sum2$habitat<- dplyr::recode(sum2$habitat, root_endosphere = "Root Endosphere")
sum2$habitat<- dplyr::recode(sum2$habitat, rhizosphere = "Rhizosphere")
sum2$habitat<- ordered(sum2$habitat, c("Root Endosphere", "Rhizosphere", "Bulk Soil"))

rel_plots_16S2 <- list()
for(i in levels(as.factor(sum2$habitat))){
  sub <- subset(sum2, habitat %in% i)
  p2 <- ggplot(sub, aes(x = drought_tolerance, y = means, fill = Family)) +
    geom_bar(position = "stack", stat = "identity", col = "black") +
    facet_wrap(~irrigation,
               labeller = labeller(irrigation = c("full" = "Full", "reduced" = "Reduced"))) +
    theme(legend.position = "right", legend.title = element_blank(), axis.line = element_blank()) +
    scale_fill_manual(values = c("firebrick4","tomato", "darkorange3", "orange1" ,"goldenrod4", "khaki1", "forestgreen",
                                 "olivedrab1", "midnightblue","powderblue", "mediumorchid4")) +
    scale_x_discrete(labels=c("Tolerant", "Susceptible")) +
    scale_y_continuous(name = NULL, breaks = NULL) +
    panel_border() +
    ggtitle(i) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(strip.text = element_text(size = 12, face = "italic"), axis.text.x = element_text(size = 10),
          axis.title = element_blank(),
          title = element_text(size = 12))
  rel_plots_16S2[[i]] <- p2 + theme(legend.position = "none")
}
leg <- get_legend(p2 + theme(legend.position = "bottom"))
grid_rel_plots_16S2 <- plot_grid(plotlist = rel_plots_16S2, nrow =1, ncol= 3, align = "v")
rel_16S_2 <- plot_grid(grid_rel_plots_16S2, leg, ncol = 1, nrow =2, rel_heights = c(1,0.2))
rel_16S_2
ggsave(filename=paste("family_16S.tiff", sep=""), plot=rel_16S_2, width=7.5, height=5, dpi=600)


# Genus
genus_sum_16S <- tapply(taxa_sums(data_16S_rel), 
                        tax_table(data_16S_rel)[, "Genus"], sum, na.rm=TRUE)
top10genus_16S <- names(sort(genus_sum_16S, TRUE))[1:13]

abundant_16S_genus <- subset_taxa(data_16S_rel, Genus %in% top10genus_16S)
abundant_16S_genus_glom <- tax_glom(abundant_16S_genus, taxrank = "Genus")
abundant_16S_genus_melt <- psmelt(abundant_16S_genus_glom)

sum2 <- abundant_16S_genus_melt %>% 
  group_by(habitat, Genus, irrigation, drought_tolerance) %>%
  summarise(means = mean(Abundance))

sum2$Genus <- dplyr::recode(sum2$Genus, "uncultured bacterium" = "Other/Unidentified Genus")
sum2$Genus <- dplyr::recode(sum2$Genus, "uncultured" = "Other/Unidentified Genus")
sum2$Genus <- dplyr::recode(sum2$Genus, "metagenome" = "Other/Unidentified Genus")


#Everything that's left is considered rare  
rare <- sum2 %>%
  group_by(habitat, irrigation, drought_tolerance) %>%
  summarise(means = 1- sum(means)) %>%
  mutate(Genus = "Other/Unidentified Genus")


#concatenate the datasets
sum2 = rbind(sum2, rare)

#genus groups
sum2$Genus <- forcats::fct_relevel(sum2$Genus, "Other/Unidentified Genus", after = Inf)

#genus habitat
sum2$habitat<- dplyr::recode(sum2$habitat, soil = "Bulk Soil")
sum2$habitat<- dplyr::recode(sum2$habitat, root_endosphere = "Root Endosphere")
sum2$habitat<- dplyr::recode(sum2$habitat, rhizosphere = "Rhizosphere")
sum2$habitat<- ordered(sum2$habitat, c("Root Endosphere", "Rhizosphere", "Bulk Soil"))

rel_plots_16S2 <- list()
for(i in levels(as.factor(sum2$habitat))){
  sub <- subset(sum2, habitat %in% i)
  p2 <- ggplot(sub, aes(x = drought_tolerance, y = means, fill = Genus)) +
    geom_bar(position = "stack", stat = "identity", col = "black") +
    facet_wrap(~irrigation,
               labeller = labeller(irrigation = c("full" = "Full", "reduced" = "Reduced"))) +
    theme(legend.position = "right", legend.title = element_blank(), axis.line = element_blank()) +
    scale_fill_manual(values = c("firebrick4","tomato", "darkorange3", "orange1" ,"goldenrod4", "khaki1", "forestgreen",
                                 "olivedrab1", "midnightblue","powderblue", "mediumorchid4")) +
    scale_x_discrete(labels=c("Tolerant", "Susceptible")) +
    scale_y_continuous(name = NULL, breaks = NULL) +
    panel_border() +
    ggtitle(i) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(strip.text = element_text(size = 12, face = "italic"), axis.text.x = element_text(size = 10),
          axis.title = element_blank(),
          title = element_text(size = 12))
  rel_plots_16S2[[i]] <- p2 + theme(legend.position = "none")
}
leg <- get_legend(p2 + theme(legend.position = "bottom"))
grid_rel_plots_16S2 <- plot_grid(plotlist = rel_plots_16S2, nrow =1, ncol= 3, align = "v")
rel_16S_2 <- plot_grid(grid_rel_plots_16S2, leg, ncol = 1, nrow =2, rel_heights = c(1,0.2))
rel_16S_2
ggsave(filename=paste("genus_16S.tiff", sep=""), plot=rel_16S_2, width=7.5, height=5, dpi=600)

# FUNGuild Analysis/Plots ------------------------------------------------------

# Append taxonomic information to the OTU table to run through FUNGuild
ASV_ITS <- data_ITS_counts@otu_table
ASV_ITS <- as.data.frame(ASV_ITS)
taxaID <- rownames(ASV_ITS)
ASV_ITS_t <- t(ASV_ITS)
sampleID <- rownames(ASV_ITS_t)
ASV_ITS_t <- as.data.frame(ASV_ITS_t)

taxa_ITS <- as.data.frame(data_ITS_counts@tax_table) %>% tibble::rownames_to_column(var="SVID")
ASV_ITS <- as.data.frame(ASV_ITS) %>% tibble::rownames_to_column(var="SVID")
otu_table_wtax <- ASV_ITS %>% left_join (taxa_ITS, by="SVID")
write.csv(otu_table_wtax, file="ITS-ASV-table-wtax.csv")

# Evaluate all guild classified OTUs and their % relative abundance across treatment 


FunGuild_rhizo <- read.delim("ITS-ASV-table-wtax.guilds_matched.txt")

funguild_tax_rhizo <- as.data.frame(as(tax_table(rare_ITS_rhizosphere), "matrix"))
funguild_tax_rhizo2 <- merge(funguild_tax_rhizo, FunGuild_rhizo %>% dplyr::select(SVID, Guild), by.x = "row.names", by.y = "SVID", all.x = T)
funguild_tax_rhizo2 <- as(funguild_tax_rhizo2, "matrix")
funguild_tax_rhizo3 <- tax_table(funguild_tax_rhizo2)

taxa_names(funguild_tax_rhizo3) <- funguild_tax_rhizo2[,1]
data_ITS_rhizo_guild <- merge_phyloseq(otu_table(rare_ITS_rhizosphere), funguild_tax_rhizo3, rare_ITS_rhizosphere@sam_data)
rhizo_guild_df <- as.data.frame(data_ITS_rhizo_guild@tax_table) 
rhizo_guild_df$freq <- sapply(strsplit(as.character(rhizo_guild_df$Guild), "-"), length) # Identify the number of FunGuild classifications
rhizo_otu_df <- as.data.frame(data_ITS_rhizo_guild@otu_table)
rhizo_otu_df <- setDT(rhizo_otu_df, keep.rownames = "Row.names")[]
rhizo_merged_df <- merge(rhizo_guild_df, rhizo_otu_df, by = "Row.names")
rhizo_merged_df2 <- rhizo_merged_df[,-1]
rownames(rhizo_merged_df2) <- rhizo_merged_df[,1]
rhizo_merged_df2 <- rhizo_merged_df2[-c(1:10)]
rhizo_merged_df3 <- (rhizo_merged_df2/ rhizo_merged_df2$freq) # divide ASV abundance by the number of FunGuild classifications 
rhizo_merged_df3 <- rhizo_merged_df3[-c(1)]
rhizo_merged_df3 <- as(rhizo_merged_df3, "matrix")
rhizo_merged_df4 <- otu_table(rhizo_merged_df3, taxa_are_rows = T)
data_ITS_rhizo_guild_norm <- merge_phyloseq(otu_table(rhizo_merged_df4), funguild_tax_rhizo3, rare_ITS_rhizosphere@sam_data)

data_ITS_rhizo_guild_norm@sam_data$irrigation.drought_tolerance <- with(data_ITS_rhizo_guild_norm@sam_data, 
                                                                        interaction(irrigation, drought_tolerance))


# Root Endosphere FunGuild Overview
FunGuild_re <- read.delim("ITS-ASV-table-wtax.guilds_matched.txt")
funguild_tax_re <- as.data.frame(as(tax_table(rare_ITS_RE), "matrix"))
funguild_tax_re2 <- merge(funguild_tax_re, FunGuild_re %>% dplyr::select(SVID, Guild), by.x = "row.names", by.y = "SVID", all.x = T)
funguild_tax_re2 <- as(funguild_tax_re2, "matrix")
funguild_tax_re3 <- tax_table(funguild_tax_re2)

taxa_names(funguild_tax_re3) <- funguild_tax_re2[,1]
data_ITS_re_guild <- merge_phyloseq(otu_table(rare_ITS_RE), funguild_tax_re3, rare_ITS_RE@sam_data)
re_guild_df <- as.data.frame(data_ITS_re_guild@tax_table) 
re_guild_df$freq <- sapply(strsplit(as.character(re_guild_df$Guild), "-"), length) # Identify the number of FunGuild classifications
re_otu_df <- as.data.frame(data_ITS_re_guild@otu_table)
re_otu_df <- setDT(re_otu_df, keep.rownames = "Row.names")[]
re_merged_df <- merge(re_guild_df, re_otu_df, by = "Row.names")
re_merged_df2 <- re_merged_df[,-1]
rownames(re_merged_df2) <- re_merged_df[,1]
re_merged_df2 <- re_merged_df2[-c(1:10)]
re_merged_df3 <- (re_merged_df2/ re_merged_df2$freq) # divide ASV abundance by the number of FunGuild classifications 
re_merged_df3 <- re_merged_df3[-c(1)]
re_merged_df3 <- as(re_merged_df3, "matrix")
re_merged_df4 <- otu_table(re_merged_df3, taxa_are_rows = T)
data_ITS_re_guild_norm <- merge_phyloseq(otu_table(re_merged_df4), funguild_tax_re3, rare_ITS_RE@sam_data)

data_ITS_re_guild_norm@sam_data$irrigation.drought_tolerance <- with(data_ITS_re_guild@sam_data, 
                                                                     interaction(irrigation, drought_tolerance))


# Evaluation of Mycorrhizal Fungi in the Root Endosphere and Rhizosphere accross treatment

# Rhizosphere
data_rel_ITS_rhizo_tax <- as.data.frame(data_ITS_rhizo_guild_norm@tax_table) # reformat names in tax table to include all labeled mycorrhizal otus
data_rel_ITS_rhizo_tax$Guild[str_detect(data_rel_ITS_rhizo_tax$Guild, "Ectomycorrhizal") == TRUE] <- "Ectomycorrhizal"
data_rel_ITS_rhizo_tax$Guild[str_detect(data_rel_ITS_rhizo_tax$Guild, "Arbuscular Mycorrhizal") == TRUE] <- "Arbuscular Mycorrhizal"
data_rel_ITS_rhizo_tax <- as(data_rel_ITS_rhizo_tax, "matrix")
data_rel_ITS_rhizo_tax <- tax_table(data_rel_ITS_rhizo_tax)
data_rel_ITS_rhizo <- merge_phyloseq(otu_table(rhizo_merged_df4), data_rel_ITS_rhizo_tax, rare_ITS_rhizosphere@sam_data)
data_rel_ITS_rhizo <- transform_sample_counts(data_ITS_rhizo_guild_norm, function(x) x/sum(x))
Myco <- subset_taxa(data_rel_ITS_rhizo, Guild == "Ectomycorrhizal" | Guild == "Arbuscular Mycorrhizal")
Myco_ecto <- subset_taxa(data_rel_ITS_rhizo, Guild == "Ectomycorrhizal")
Myco_arb <- subset_taxa(data_rel_ITS_rhizo, Guild == "Arbuscular Mycorrhizal")
myco_melt_rhizo <- psmelt(Myco)
myco_melt_ecto <- psmelt(Myco_ecto)
myco_melt_arb <- psmelt(Myco_arb)

# Ectomycorrhizae Abundance
myco_melt_ecto %>%
  group_by(irrigation, drought_tolerance, Guild, Sample) %>%
  summarise(sum = sum(Abundance)*100) %>%
  ungroup(irrigation, drought_tolerance, Guild, Sample) %>%
  group_by(irrigation, drought_tolerance, Guild) %>%
  summarise(mean = mean(sum), se = sd(sum)/sqrt(n())) %>%
  ggplot(aes(x = irrigation, y = mean, fill = drought_tolerance)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), alpha=1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~Guild) +
  scale_y_continuous(name = "Relative Abundance (%)") +
  xlab("Irrigation") +
  scale_x_discrete(labels=c("Full", "Reduced")) +
  theme(axis.title.x = element_blank(),
        legend.position = "right", 
        axis.text.x = element_text(face = "italic"),
        legend.justification = "center",
        legend.title = element_text(size = 12)) +
  ggtitle("Rhizosphere Ectomycorrhizae Abundance") +
  theme_bw() +
  theme(legend.text.align = 0) +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  guides(fill = guide_legend(title = "Drought Tolerance")) -> p_ecto_abundance


# Ectomycorrhizae Taxonomy
Myco_ecto_rel <- psmelt(Myco_ecto)

sum2 <- Myco_ecto_rel %>% 
  group_by(irrigation, drought_tolerance, Genus) %>%
  summarise(sum = sum(Abundance)*100)


ggplot(sum2, aes(x = drought_tolerance, y = sum, fill = Genus)) +
  geom_bar(position = "stack", stat = "identity", col = "black") +
  facet_wrap(~irrigation,
             labeller = labeller(irrigation = c("full" = "Full", "reduced" = "Reduced"))) +
  theme(legend.position = "right", legend.title = element_blank(), axis.line = element_blank()) +
  scale_fill_manual(values = c("firebrick4","tomato", "forestgreen",
                               "olivedrab1", "midnightblue", "powderblue", "mediumorchid4", "red","blue","green","yellow")) +
  scale_x_discrete(labels=c("Tolerant", "Susceptible")) +
  scale_y_continuous(name = NULL, breaks = NULL) +
  panel_border() +
  ggtitle("    Rhizosphere Ectomycorrhizae Taxonomy") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text = element_text(size = 12, face = "italic"), axis.text.x = element_text(size = 10),
        axis.title = element_blank(),
        title = element_text(size = 12)) -> p_ecto_taxonomy

# Append plots together
prow <- plot_grid(p_ecto_abundance + theme(legend.position = "bottom"), 
                  p_ecto_taxonomy + theme(legend.position = "bottom"), 
                  align = "h", axis = "lr", labels = "AUTO", ncol = 2)
prow
ggsave(filename=paste("rhizo_ecto.tiff",sep=""), plot=prow, width=10.5, height=7, dpi=600)

# Arbuscular Mycorrhizal Abundance
myco_melt_arb %>%
  group_by(irrigation, drought_tolerance, Guild, Sample) %>%
  summarise(sum = sum(Abundance)*100) %>%
  ungroup(irrigation, drought_tolerance, Guild, Sample) %>%
  group_by(irrigation, drought_tolerance, Guild) %>%
  summarise(mean = mean(sum), se = sd(sum)/sqrt(n())) %>%
  ggplot(aes(x = irrigation, y = mean, fill = drought_tolerance)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), alpha=1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~Guild) +
  scale_y_continuous(name = "Relative Abundance (%)") +
  xlab("Irrigation") +
  scale_x_discrete(labels=c("Full", "Reduced")) +
  theme(axis.title.x = element_blank(),
        legend.position = "right", 
        axis.text.x = element_text(face = "italic"),
        legend.justification = "center",
        legend.title = element_text(size = 12)) +
  ggtitle("Rhizosphere Arbuscular Mycorrhizal Abundance") +
  theme_bw() +
  theme(legend.text.align = 0) +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  guides(fill = guide_legend(title = "Drought Tolerance")) -> p_arb_abundance

# Arbuscular Mycorrhizal Taxonomy 
Myco_arb_tax_df <- as.data.frame(Myco_arb@tax_table)
Myco_arb_tax_df[is.na(Myco_arb_tax_df)] <- "Unidentified Genus"                  
Myco_arb_tax_df$Genus[str_detect(Myco_arb_tax_df$Genus, "unidentified") == TRUE] <- "Unidentified Genus"
Myco_arb_tax_df <- as(Myco_arb_tax_df, "matrix")
Myco_arb_tax_df <- tax_table(Myco_arb_tax_df)
Myco_arb <- merge_phyloseq(otu_table(Myco_arb), Myco_arb_tax_df, Myco_arb@sam_data)

Myco_arb_rel <- psmelt(Myco_arb)

sum2 <- Myco_arb_rel %>% 
  group_by(irrigation, drought_tolerance, Genus) %>%
  summarise(sum = sum(Abundance)*100)


ggplot(sum2, aes(x = drought_tolerance, y = sum, fill = Genus)) +
  geom_bar(position = "stack", stat = "identity", col = "black") +
  facet_wrap(~irrigation,
             labeller = labeller(irrigation = c("full" = "Full", "reduced" = "Reduced"))) +
  theme(legend.position = "right", legend.title = element_blank(), axis.line = element_blank()) +
  scale_fill_manual(values = c("firebrick4","tomato", "forestgreen",
                               "olivedrab1", "midnightblue", "powderblue", "mediumorchid4", "red","blue","green","yellow")) +
  scale_x_discrete(labels=c("Tolerant", "Susceptible")) +
  scale_y_continuous(name = NULL, breaks = NULL) +
  panel_border() +
  ggtitle("    Rhizosphere Arbuscular Mycorrhizae Taxonomy") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text = element_text(size = 12, face = "italic"), axis.text.x = element_text(size = 10),
        axis.title = element_blank(),
        title = element_text(size = 12)) -> p_arb_taxonomy

# Append Plots Together 
prow <- plot_grid(p_arb_abundance + theme(legend.position = "bottom"), 
                  p_arb_taxonomy + theme(legend.position = "bottom"), 
                  align = "h", axis = "lr", labels = "AUTO", ncol = 2)
prow
ggsave(filename=paste("rhizo_arb.tiff",sep=""), plot=prow, width=10.5, height=7, dpi=600)

# Root Endosphere
data_rel_ITS_re_tax <- as.data.frame(data_ITS_re_guild_norm@tax_table) # reformat names in tax table to include all labeled mycorrhizal otus
data_rel_ITS_re_tax$Guild[str_detect(data_rel_ITS_re_tax$Guild, "Ectomycorrhizal") == TRUE] <- "Ectomycorrhizal"
data_rel_ITS_re_tax$Guild[str_detect(data_rel_ITS_re_tax$Guild, "Arbuscular Mycorrhizal") == TRUE] <- "Arbuscular Mycorrhizal"
data_rel_ITS_re_tax <- as(data_rel_ITS_re_tax, "matrix")
data_rel_ITS_re_tax <- tax_table(data_rel_ITS_re_tax)
data_rel_ITS_re <- merge_phyloseq(otu_table(re_merged_df4), data_rel_ITS_re_tax, rare_ITS_RE@sam_data)
data_rel_ITS_re <- transform_sample_counts(data_ITS_re_guild_norm, function(x) x/sum(x))
Myco <- subset_taxa(data_rel_ITS_re, Guild == "Ectomycorrhizal" | Guild == "Arbuscular Mycorrhizal")
Myco_ecto <- subset_taxa(data_rel_ITS_re, Guild == "Ectomycorrhizal")
Myco_arb <- subset_taxa(data_rel_ITS_re, Guild == "Arbuscular Mycorrhizal")
myco_melt_re <- psmelt(Myco)
myco_melt_ecto <- psmelt(Myco_ecto)
myco_melt_arb <- psmelt(Myco_arb)

# Ectomycorrhizae Abundance
myco_melt_ecto %>%
  group_by(irrigation, drought_tolerance, Guild, Sample) %>%
  summarise(sum = sum(Abundance)*100) %>%
  ungroup(irrigation, drought_tolerance, Guild, Sample) %>%
  group_by(irrigation, drought_tolerance, Guild) %>%
  summarise(mean = mean(sum), se = sd(sum)/sqrt(n())) %>%
  ggplot(aes(x = irrigation, y = mean, fill = drought_tolerance)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), alpha=1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~Guild) +
  scale_y_continuous(name = "Relative Abundance (%)") +
  xlab("Irrigation") +
  scale_x_discrete(labels=c("Full", "Reduced")) +
  theme(axis.title.x = element_blank(),
        legend.position = "right", 
        axis.text.x = element_text(face = "italic"),
        legend.justification = "center",
        legend.title = element_text(size = 12)) +
  ggtitle("Root Endosphere Ectomycorrhizae Abundance") +
  theme_bw() +
  theme(legend.text.align = 0) +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  guides(fill = guide_legend(title = "Drought Tolerance")) -> p_ecto_abundance


# Ectomycorrhizae Taxonomy
Myco_ecto_rel <- psmelt(Myco_ecto)

sum2 <- Myco_ecto_rel %>% 
  group_by(irrigation, drought_tolerance, Genus) %>%
  summarise(sum = sum(Abundance)*100)


ggplot(sum2, aes(x = drought_tolerance, y = sum, fill = Genus)) +
  geom_bar(position = "stack", stat = "identity", col = "black") +
  facet_wrap(~irrigation,
             labeller = labeller(irrigation = c("full" = "Full", "reduced" = "Reduced"))) +
  theme(legend.position = "right", legend.title = element_blank(), axis.line = element_blank()) +
  scale_fill_manual(values = c("firebrick4","tomato", "forestgreen",
                               "olivedrab1", "midnightblue", "powderblue", "mediumorchid4", "red","blue","green","yellow")) +
  scale_x_discrete(labels=c("Tolerant", "Susceptible")) +
  scale_y_continuous(name = NULL, breaks = NULL) +
  panel_border() +
  ggtitle("    Root Endosphere Ectomycorrhizae Taxonomy") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text = element_text(size = 12, face = "italic"), axis.text.x = element_text(size = 10),
        axis.title = element_blank(),
        title = element_text(size = 12)) -> p_ecto_taxonomy

# Append plots together
prow <- plot_grid(p_ecto_abundance + theme(legend.position = "bottom"), 
                  p_ecto_taxonomy + theme(legend.position = "bottom"), 
                  align = "h", axis = "lr", labels = "AUTO", ncol = 2)
prow
ggsave(filename=paste("root_endosphere_ecto.tiff",sep=""), plot=prow, width=9.5, height=6, dpi=600)

# Arbuscular Mycorrhizal Abundance
myco_melt_arb %>%
  group_by(irrigation, drought_tolerance, Guild, Sample) %>%
  summarise(sum = sum(Abundance)*100) %>%
  ungroup(irrigation, drought_tolerance, Guild, Sample) %>%
  group_by(irrigation, drought_tolerance, Guild) %>%
  summarise(mean = mean(sum), se = sd(sum)/sqrt(n())) %>%
  ggplot(aes(x = irrigation, y = mean, fill = drought_tolerance)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), alpha=1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~Guild) +
  scale_y_continuous(name = "Relative Abundance (%)") +
  xlab("Irrigation") +
  scale_x_discrete(labels=c("Full", "Reduced")) +
  theme(axis.title.x = element_blank(),
        legend.position = "right", 
        axis.text.x = element_text(face = "italic"),
        legend.justification = "center",
        legend.title = element_text(size = 12)) +
  ggtitle("Root Endosphere Arbuscular Mycorrhizal Abundance") +
  theme_bw() +
  theme(legend.text.align = 0) +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  guides(fill = guide_legend(title = "Drought Tolerance")) -> p_arb_abundance

# Arbuscular Mycorrhizal Taxonomy 
Myco_arb_tax_df <- as.data.frame(Myco_arb@tax_table)
Myco_arb_tax_df[is.na(Myco_arb_tax_df)] <- "Unidentified Genus"                  
Myco_arb_tax_df$Genus[str_detect(Myco_arb_tax_df$Genus, "unidentified") == TRUE] <- "Unidentified Genus"
Myco_arb_tax_df <- as(Myco_arb_tax_df, "matrix")
Myco_arb_tax_df <- tax_table(Myco_arb_tax_df)
Myco_arb <- merge_phyloseq(otu_table(Myco_arb), Myco_arb_tax_df, Myco_arb@sam_data)

Myco_arb_rel <- psmelt(Myco_arb)

sum2 <- Myco_arb_rel %>% 
  group_by(irrigation, drought_tolerance, Genus) %>%
  summarise(sum = sum(Abundance)*100)


ggplot(sum2, aes(x = drought_tolerance, y = sum, fill = Genus)) +
  geom_bar(position = "stack", stat = "identity", col = "black") +
  facet_wrap(~irrigation,
             labeller = labeller(irrigation = c("full" = "Full", "reduced" = "Reduced"))) +
  theme(legend.position = "right", legend.title = element_blank(), axis.line = element_blank()) +
  scale_fill_manual(values = c("firebrick4","tomato", "forestgreen",
                               "olivedrab1", "midnightblue", "powderblue", "mediumorchid4", "red","blue","green","yellow")) +
  scale_x_discrete(labels=c("Tolerant", "Susceptible")) +
  scale_y_continuous(name = NULL, breaks = NULL) +
  panel_border() +
  ggtitle("   Root Endosphere Arbuscular Mycorrhizae Taxonomy") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text = element_text(size = 12, face = "italic"), axis.text.x = element_text(size = 10),
        axis.title = element_blank(),
        title = element_text(size = 12)) -> p_arb_taxonomy

# Append Plots Together 
prow <- plot_grid(p_arb_abundance + theme(legend.position = "bottom"), 
                  p_arb_taxonomy + theme(legend.position = "bottom"), 
                  align = "h", axis = "lr", labels = "AUTO", ncol = 2)
prow
ggsave(filename=paste("root_endosphere_arb.tiff",sep=""), plot=prow, width=11.5, height=6, dpi=600)


# Statistics 
model <- lme(Abundance ~ irrigation*drought_tolerance, random=~1|Row.names, data = myco_melt_rhizo, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) 

x <- aov(sum ~ irrigation*drought_tolerance, data = myco_melt_rhizo %>%
           group_by(irrigation, drought_tolerance, Guild, Sample) %>%
           summarise(sum = sum(Abundance)*100) %>%
           filter(Guild == "Ectomycorrhizal"))
summary(x) 
TukeyHSD(x)

x <- aov(sum ~ irrigation*drought_tolerance, data = myco_melt_arb %>%
           group_by(irrigation, drought_tolerance, Guild, Sample) %>%
           summarise(sum = sum(Abundance)*100) %>%
           filter(Guild == "Arbuscular Mycorrhizal"))
summary(x) 
TukeyHSD(x)

# Statistics 
model <- lme(Abundance ~ irrigation*drought_tolerance, random=~1|Row.names, data = myco_melt_re, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) 
summary(glht(model, mcp(drought_tolerance="Tukey")))

x <- aov(sum ~ irrigation*drought_tolerance, data = myco_melt_re %>%
           group_by(irrigation, drought_tolerance, Guild, Sample) %>%
           summarise(sum = sum(Abundance)*100) %>%
           filter(Guild == "Ectomycorrhizal"))
summary(x) 
TukeyHSD(x)

x <- aov(sum ~ irrigation*drought_tolerance, data = myco_melt_re %>%
           group_by(irrigation, drought_tolerance, Guild, Sample) %>%
           summarise(sum = sum(Abundance)*100) %>%
           filter(Guild == "Arbuscular Mycorrhizal"))
summary(x) 
TukeyHSD(x)

# Evaluate of Pathogenic Fungi in the Root Endosphere and Rhizosphere accross treatment

# Rhizosphere
data_rel_ITS_rhizo_tax <- as.data.frame(data_ITS_rhizo_guild_norm@tax_table) # reformat names in tax table to include all labeled mycorrhizal otus
data_rel_ITS_rhizo_tax$Guild[str_detect(data_rel_ITS_rhizo_tax$Guild, "Plant Pathogen") == TRUE] <- "Plant Pathogen"
data_rel_ITS_rhizo_tax <- as(data_rel_ITS_rhizo_tax, "matrix")
data_rel_ITS_rhizo_tax <- tax_table(data_rel_ITS_rhizo_tax)
data_rel_ITS_rhizo <- merge_phyloseq(otu_table(rhizo_merged_df4), data_rel_ITS_rhizo_tax, rare_ITS_rhizosphere@sam_data)
data_rel_ITS_rhizo <- transform_sample_counts(data_ITS_rhizo_guild_norm, function(x) x/sum(x))
path_rhizo <- subset_taxa(data_rel_ITS_rhizo, Guild == "Plant Pathogen")
path_rhizo_melt <- psmelt(path_rhizo)

path_rhizo_melt$Guild <- as.character(path_rhizo_melt$Guild)

path_rhizo_melt %>%
  group_by(irrigation, drought_tolerance, Guild) %>%
  summarise(sum = sum(Abundance)) %>%
  mutate(Habitat = "Rhizosphere") -> path_rhizosphere

path_rhizo_melt %>%
  group_by(irrigation, drought_tolerance, Guild, Sample) %>%
  summarise(sum = sum(Abundance)*100) %>%
  ungroup(irrigation, drought_tolerance, Guild, Sample) %>%
  group_by(irrigation, drought_tolerance, Guild) %>%
  summarise(mean = mean(sum), se = sd(sum)/sqrt(n())) %>%
  ggplot(aes(x = irrigation, y = mean, fill = drought_tolerance)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), alpha=1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~Guild, scales = "free") +
  scale_y_continuous(name = "Relative Abundance (%)") +
  xlab("Irrigation") +
  scale_x_discrete(labels=c("Full", "Reduced")) +
  theme(axis.title.x = element_blank(),
        legend.position = "right", 
        axis.text.x = element_text(face = "italic"),
        #legend.justification = "center",
        legend.title = element_text(size = 12)) +
  ggtitle("Rhizosphere Plant Pathogen Abundance") +
  theme_bw() +
  theme(legend.text.align = 0) +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  guides(fill = guide_legend(title = "Drought Tolerance")) -> p_path_rhizo 
p_path_rhizo

# Root Endosphere
data_rel_ITS_re_tax <- as.data.frame(data_ITS_re_guild_norm@tax_table) # reformat names in tax table to include all labeled mycorrhizal otus
data_rel_ITS_re_tax$Guild[str_detect(data_rel_ITS_re_tax$Guild, "Plant Pathogen") == TRUE] <- "Plant Pathogen"
data_rel_ITS_re_tax <- as(data_rel_ITS_re_tax, "matrix")
data_rel_ITS_re_tax <- tax_table(data_rel_ITS_re_tax)
data_rel_ITS_re <- merge_phyloseq(otu_table(re_merged_df4), data_rel_ITS_re_tax, rare_ITS_RE@sam_data)
data_rel_ITS_re <- transform_sample_counts(data_ITS_re_guild_norm, function(x) x/sum(x))
path_re <- subset_taxa(data_rel_ITS_re, Guild == "Plant Pathogen")
path_melt_re <- psmelt(path_re)

path_melt_re$Guild <- as.character(path_melt_re$Guild)

path_melt_re %>%
  group_by(irrigation, drought_tolerance, Guild) %>%
  summarise(sum = sum(Abundance)) %>%
  mutate(Habitat = "Root Endosphere") -> path_re

path_melt_re %>%
  group_by(irrigation, drought_tolerance, Guild, Sample) %>%
  summarise(sum = sum(Abundance)*100) %>%
  ungroup(irrigation, drought_tolerance, Guild, Sample) %>%
  group_by(irrigation, drought_tolerance, Guild) %>%
  summarise(mean = mean(sum), se = sd(sum)/sqrt(n())) %>%
  ggplot(aes(x = irrigation, y = mean, fill = drought_tolerance)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), alpha = 1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~Guild, scales = "free") +
  scale_y_continuous(name = "Relative Abundance (%)") +
  xlab("Irrigation") +
  scale_x_discrete(labels=c("Full", "Reduced")) +
  theme(axis.title.x = element_blank(),
        legend.position = "right", 
        axis.text.x = element_text(face = "italic"),
        #legend.justification = "center",
        legend.title = element_text(size = 12)) +
  ggtitle("Root Endosphere Plant Pathogen Abundance") +
  theme_bw() +
  theme(legend.text.align = 0) +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  guides(fill = guide_legend(title = "Drought Tolerance")) -> p_path_re
p_path_re

prow <- plot_grid(p_path_re + theme(legend.position = "none"), p_path_rhizo + theme(legend.position = "none"), align = "h", axis = "lr", labels = "AUTO", ncol = 2)
leg <- get_legend(p_path_re)
prow <- plot_grid(prow, leg, ncol = 2, rel_widths = c(1.3, 0.23))
prow
ggsave(filename=paste("path_abund_sum.tiff",sep=""), plot=prow, width=10, height=5, dpi=600)

#Statistics
x <- aov(sum ~ irrigation*drought_tolerance, data = path_rhizo_melt %>%
           group_by(irrigation, drought_tolerance, Guild, Sample) %>%
           summarise(sum = sum(Abundance)*100) %>%
           filter(Guild == "Plant Pathogen"))
summary(x) 

x <- aov(sum ~ irrigation*drought_tolerance, data = path_melt_re %>%
           group_by(irrigation, drought_tolerance, Guild, Sample) %>%
           summarise(sum = sum(Abundance)*100) %>% 
           filter(Guild == "Plant Pathogen"))
summary(x) 
TukeyHSD(x)

#FUNGuild Diversity Plots ---------------------------------------------------------------------------
# Rhizosphere - Arbuscular mycorrhizae
Myco_arb_rhizo <- t(as(Myco_arb@otu_table, "matrix"))
Myco_arb_rhizo <- as.data.frame(Myco_arb_rhizo)

rich_ITS <- hill_taxa(Myco_arb_rhizo, q = 0)
rich_ITS <- merge(data_ITS_counts@sam_data, rich_ITS, by = 0)

shan_ITS <- hill_taxa(Myco_arb_rhizo, q = 1)
shan_ITS <- merge(data_ITS_counts@sam_data, shan_ITS, by = 0)

# ITS Alpha Diversity Statistics

model <- lme(y ~ drought_tolerance + irrigation, random=~1|Row.names, data = rich_ITS, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) 

model <- lme(y ~ drought_tolerance + irrigation, random=~1|Row.names, data = shan_ITS, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) 

rich_ITS %>%
  ggplot(aes(x = irrigation, y = y, fill = drought_tolerance)) +
  geom_boxplot() +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) +
  scale_y_continuous(name = "Species Richness") +
  scale_fill_manual("Drought Tolerance",
                    values = c("tomato1","steelblue4"),
                    labels = c(expression(paste("Drought Tolerant", sep="")),
                               expression(paste("Drought Susceptible", sep="")))) +
  panel_border() +
  theme_bw() +
  ggtitle(label = "Arbuscular Mycorrhizae") +
  theme(plot.title = element_text(size = 9.5)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        legend.title = element_text(size = 12)) +
  geom_vline(xintercept = c(1.5), linetype = "dotted") -> p_arb
p_arb

# Rhizosphere - Ectomycorrhizae
Myco_ecto_rhizo <- t(as(Myco_ecto@otu_table, "matrix"))
Myco_ecto_rhizo <- as.data.frame(Myco_ecto_rhizo)

rich_ITS <- hill_taxa(Myco_ecto_rhizo, q = 0)
rich_ITS <- merge(data_ITS_counts@sam_data, rich_ITS, by = 0)

shan_ITS <- hill_taxa(Myco_ecto_rhizo, q = 1)
shan_ITS <- merge(data_ITS_counts@sam_data, shan_ITS, by = 0)

model <- lme(y ~ drought_tolerance + irrigation, random=~1|Row.names, data = rich_ITS, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) 

model <- lme(y ~ drought_tolerance + irrigation, random=~1|Row.names, data = shan_ITS, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) 

rich_ITS %>%
  ggplot(aes(x = irrigation, y = y, fill = drought_tolerance)) +
  geom_boxplot() +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) +
  scale_y_continuous(name = "Species Richness") +
  scale_fill_manual("Drought Tolerance",
                    values = c("tomato1","steelblue4"),
                    labels = c(expression(paste("Drought Tolerant", sep="")),
                               expression(paste("Drought Susceptible", sep="")))) +
  panel_border() +
  theme_bw() +
  ggtitle(label = "Ectomycorrhizae") +
  theme(plot.title = element_text(size = 9.5)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        legend.title = element_text(size = 12)) +
  geom_vline(xintercept = c(1.5), linetype = "dotted") -> p_ecto
p_ecto
leg <- get_legend(p_ecto + theme(legend.position = "right"))
prow <- plot_grid(p_arb + theme(legend.position = "none"), p_ecto + theme(legend.position = "none"), align = "h", ncol = 2)
prow <- plot_grid(prow, leg, ncol = 2, rel_widths = c(1.5, 0.5))
prow
ggsave(filename=paste("Myco_alpha_div.tiff",sep=""), plot=prow, width=8, height=5, dpi=400)

# Calculate the Rarefaction Efficiency Index -----------------------------------
# REI Function is adapted from https://github.com/jcyhong/rarefaction/blob/master/rarefaction_efficiency_index.R)

REI <- function(otu_matrix, rarefied_depth, group, taxa_are_rows=FALSE) {
  # Compute rarefaction efficiency indices based on a given OTU matrix
  # and the rarefied depth.
  #
  # Input(s):
  #  - otu_matrix: a matrix of counts (taxa are columns)
  #  - rarefied_depth: a scalar indicating the rarefied depth
  #  - group: a vector indicating the group memberships
  #  - taxa_are_rows: a boolean, TRUE if each row represents a taxon
  #
  # Returns:
  #  - a vector of REIs
  
  if (taxa_are_rows) {
    otu_matrix <- t(otu_matrix)
  }
  library_size <- rowSums(otu_matrix)
  
  # Drop all observations with library sizes less than rarefied depth.
  otu_matrix <- otu_matrix[library_size >= rarefied_depth, ]
  group <- group[library_size >= rarefied_depth]
  library_size <- library_size[library_size >= rarefied_depth]
  group_lvl <- unique(group)
  n1 <- sum(group == group_lvl[1])
  n2 <- sum(group == group_lvl[2])
  
  var_prop1 <- apply(otu_matrix[group == group_lvl[1], ] / 
                       library_size[group == group_lvl[1]], 
                     2, function(col) {
                       var(col)
                     })
  var_prop_rrf1 <- apply(otu_matrix[group == group_lvl[1], ] / 
                           library_size[group == group_lvl[1]], 
                         2, function(col) {
                           var(col) + 
                             1 / rarefied_depth * mean(col * (1 - col) * (library_size[group == group_lvl[1]] - rarefied_depth) / 
                                                         (library_size[group == group_lvl[1]] - 1))
                         })
  
  var_prop2 <- apply(otu_matrix[group == group_lvl[2], ] / 
                       library_size[group == group_lvl[2]], 
                     2, function(col) {
                       var(col)
                     })
  var_prop_rrf2 <- apply(otu_matrix[group == group_lvl[2], ]  / 
                           library_size[group == group_lvl[2]], 
                         2, function(col) {
                           var(col) + 
                             1 / rarefied_depth * mean(col * (1 - col) * (library_size[group == group_lvl[2]] - rarefied_depth) / 
                                                         (library_size[group == group_lvl[2]] - 1))
                         })
  
  (var_prop1 / n1 + var_prop2 / n2) / 
    (var_prop_rrf1 / n1 + var_prop_rrf2 / n2)
}
# 16S Libaries

# Root Endosphere

Rarefaction <- 14000
RE_metadata <- as.data.frame(data_16S_RE@sam_data)
# Irrigation
bool <- 
  (RE_metadata$irrigation == 'full' |
     RE_metadata$irrigation == 'reduced')
RE_otu_matrix_subset <- otu_table(data_16S_RE)[, bool]
group <- RE_metadata$irrigation[bool]

REIs <- REI(otu_matrix = RE_otu_matrix_subset, rarefied_depth = Rarefaction, group = group,
            taxa_are_rows = TRUE)
mean(REIs, na.rm=TRUE) # 0.6989861

# Drought Tolerance
bool <- 
  (RE_metadata$drought_tolerance == 'HI30' |
     RE_metadata$drought_tolerance == 'LO30')
RE_otu_matrix_subset <- otu_table(data_16S_RE)[, bool]
group <- RE_metadata$drought_tolerance[bool]

REIs <- REI(otu_matrix = RE_otu_matrix_subset, rarefied_depth = Rarefaction, group = group,
            taxa_are_rows = TRUE)
mean(REIs, na.rm=TRUE) # 0.6978382


# Rhizosphere
Rarefaction <- 25000
Rhizo_metadata <- as.data.frame(data_16S_rhizosphere@sam_data)
# Irrigation
bool <- 
  (Rhizo_metadata$irrigation == 'full' |
     Rhizo_metadata$irrigation == 'reduced')
Rhizo_otu_matrix_subset <- otu_table(data_16S_rhizosphere)[, bool]
group <- Rhizo_metadata$irrigation[bool]

REIs <- REI(otu_matrix = Rhizo_otu_matrix_subset, rarefied_depth = Rarefaction, group = group,
            taxa_are_rows = TRUE)
mean(REIs, na.rm=TRUE) #0.8267574

# Drought Tolerance
bool <- 
  (Rhizo_metadata$drought_tolerance == 'HI30' |
     Rhizo_metadata$drought_tolerance == 'LO30')
Rhizo_otu_matrix_subset <- otu_table(data_16S_rhizosphere)[, bool]
group <- Rhizo_metadata$drought_tolerance[bool]

REIs <- REI(otu_matrix = Rhizo_otu_matrix_subset, rarefied_depth = Rarefaction, group = group,
            taxa_are_rows = TRUE)
mean(REIs, na.rm=TRUE) # 0.8274541

# Bulk Soil
Rarefaction <- 45000
Soil_metadata <- as.data.frame(data_16S_soil@sam_data)
# Irrigation
bool <- 
  (Soil_metadata$irrigation == 'full' |
     Soil_metadata$irrigation == 'reduced')
Soil_otu_matrix_subset <- otu_table(data_16S_soil)[, bool]
group <- Soil_metadata$irrigation[bool]

REIs <- REI(otu_matrix = Soil_otu_matrix_subset, rarefied_depth = Rarefaction, group = group,
            taxa_are_rows = TRUE)
mean(REIs, na.rm=TRUE) #0.9249092

# Drought Tolerance
bool <- 
  (Soil_metadata$drought_tolerance == 'HI30' |
     Soil_metadata$drought_tolerance == 'LO30')
Soil_otu_matrix_subset <- otu_table(data_16S_soil)[, bool]
group <- Soil_metadata$drought_tolerance[bool]

REIs <- REI(otu_matrix = Soil_otu_matrix_subset, rarefied_depth = Rarefaction, group = group,
            taxa_are_rows = TRUE)
mean(REIs, na.rm=TRUE) # 0.9252718

# ITS Libraries
# Root Endosphere

Rarefaction <- 1000
RE_metadata <- as.data.frame(data_ITS_RE@sam_data)
# Irrigation
bool <- 
  (RE_metadata$irrigation == 'full' |
     RE_metadata$irrigation == 'reduced')
RE_otu_matrix_subset <- otu_table(data_ITS_RE)[, bool]
group <- RE_metadata$irrigation[bool]

REIs <- REI(otu_matrix = RE_otu_matrix_subset, rarefied_depth = Rarefaction, group = group,
            taxa_are_rows = TRUE)
mean(REIs, na.rm=TRUE) # 0.7044106

# Drought Tolerance
bool <- 
  (RE_metadata$drought_tolerance == 'HI30' |
     RE_metadata$drought_tolerance == 'LO30')
RE_otu_matrix_subset <- otu_table(data_ITS_RE)[, bool]
group <- RE_metadata$drought_tolerance[bool]

REIs <- REI(otu_matrix = RE_otu_matrix_subset, rarefied_depth = Rarefaction, group = group,
            taxa_are_rows = TRUE)
mean(REIs, na.rm=TRUE) # 0.7012424

# Rhizosphere
Rarefaction <- 1300
Rhizo_metadata <- as.data.frame(data_ITS_rhizosphere@sam_data)
# Irrigation
bool <- 
  (Rhizo_metadata$irrigation == 'full' |
     Rhizo_metadata$irrigation == 'reduced')
Rhizo_otu_matrix_subset <- otu_table(data_ITS_rhizosphere)[, bool]
group <- Rhizo_metadata$irrigation[bool]

REIs <- REI(otu_matrix = Rhizo_otu_matrix_subset, rarefied_depth = Rarefaction, group = group,
            taxa_are_rows = TRUE)
mean(REIs, na.rm=TRUE) #0.5293111

# Drought Tolerance
bool <- 
  (Rhizo_metadata$drought_tolerance == 'HI30' |
     Rhizo_metadata$drought_tolerance == 'LO30')
Rhizo_otu_matrix_subset <- otu_table(data_ITS_rhizosphere)[, bool]
group <- Rhizo_metadata$drought_tolerance[bool]

REIs <- REI(otu_matrix = Rhizo_otu_matrix_subset, rarefied_depth = Rarefaction, group = group,
            taxa_are_rows = TRUE)
mean(REIs, na.rm=TRUE) # 0.5307037

# Bulk Soil
Rarefaction <- 10000
Soil_metadata <- as.data.frame(data_ITS_soil@sam_data)
# Irrigation
bool <- 
  (Soil_metadata$irrigation == 'full' |
     Soil_metadata$irrigation == 'reduced')
Soil_otu_matrix_subset <- otu_table(data_ITS_soil)[, bool]
group <- Soil_metadata$irrigation[bool]

REIs <- REI(otu_matrix = Soil_otu_matrix_subset, rarefied_depth = Rarefaction, group = group,
            taxa_are_rows = TRUE)
mean(REIs, na.rm=TRUE) #0.8666217

# Drought Tolerance
bool <- 
  (Soil_metadata$drought_tolerance == 'HI30' |
     Soil_metadata$drought_tolerance == 'LO30')
Soil_otu_matrix_subset <- otu_table(data_ITS_soil)[, bool]
group <- Soil_metadata$drought_tolerance[bool]

REIs <- REI(otu_matrix = Soil_otu_matrix_subset, rarefied_depth = Rarefaction, group = group,
            taxa_are_rows = TRUE)
mean(REIs, na.rm=TRUE) # 0.8668389


# Differential Abundance Analysis ----------------------------------------------
# ANCOM Differential Abundance Analysis 
# 16S Datasets 

# Root Endosphere -- Differentially Abundant ASVs due to Irrigation with Drought Tolerant Covariate

# Pre-Processing 
# Step 1: Data preprocessing
data_RE_genus <- tax_glom(data_16S_RE, "Genus")

meta_data = data_16S_RE@sam_data
meta_data <- as.data.frame(meta_data)
meta_data$names <- rownames(meta_data)

feature_table = data_16S_RE@otu_table
feature_table <- as.data.frame(feature_table)

sample_var = "names" ; group_var = NULL
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM

main_var = "irrigation"; p_adj_method = "BH"; alpha = 0.05
adj_formula = "drought_tolerance"; rand_formula = NULL
t_start = Sys.time()
res_re_irrigation = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
                          alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start

# Append taxanomic information to the ANCOM Output 
re_tax <- as.data.frame(data_16S_RE@tax_table)

res_re_irrigation$out <- merge(res_re_irrigation$out, re_tax %>% rownames_to_column(var = "rowname"),
                               by.x = 'taxa_id', by.y = 'rowname')

write.csv(res_re_irrigation$out, "ancom_re_irrigation_ASV.csv", row.names = F)
# Step 3: Volcano Plot -- WITH Taxanomic annotation added on 

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
res_re_irrigation$fig$data <- merge (res_re_irrigation$fig$data, re_tax %>% rownames_to_column(var = "rowname"),
                                     by.x = 'taxa_id', by.y = 'rowname')
dat_ann = data.frame(x = min(res_re_irrigation$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")

fig_ancom_re_irrigation = res_re_irrigation$fig +  
  geom_point(aes(color = Phylum)) + 
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  scale_fill_brewer(palette="Set1") +
  ggtitle("Root Endosphere - Irrigation") +
  guides(color=guide_legend(title="Phylum")) +
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
fig_ancom_re_irrigation

# Root Endosphere -- Differentially Abundant ASVs due to Drought Tolerance with Irrigation Covariate 

# Step 2: ANCOM

main_var = "drought_tolerance"; p_adj_method = "BH"; alpha = 0.05
adj_formula = "irrigation"; rand_formula = NULL
t_start = Sys.time()
res_re_drought = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
                       alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start

# Append taxanomic information to the ANCOM Output 
re_tax <- as.data.frame(data_16S_RE@tax_table)

res_re_drought$out <- merge(res_re_drought$out, re_tax %>% rownames_to_column(var = "rowname"),
                            by.x = 'taxa_id', by.y = 'rowname')

write.csv(res_re_drought$out, "ancom_re_drought_tolerance_ASV.csv", row.names = F)

# Step 3: Volcano Plot -- WITH Taxanomic annotation added on 

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
res_re_drought$fig$data <- merge (res_re_drought$fig$data, re_tax %>% rownames_to_column(var = "rowname"),
                                  by.x = 'taxa_id', by.y = 'rowname')
dat_ann = data.frame(x = min(res_re_drought$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7 ]")

fig_ancom_re_drought_tolerance = res_re_drought$fig +  
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  geom_point(aes(color = Phylum)) + 
  scale_fill_brewer(palette="Set1") +
  ggtitle("Root Endosphere - Drought Tolerance") +
  guides(color=guide_legend(title="Phylum")) +
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)

fig_ancom_re_drought_tolerance

leg_16S <- get_legend(fig_ancom_re_irrigation + theme(legend.position = "bottom"))
prow <- plot_grid(fig_ancom_re_irrigation + theme(legend.position = "none"), fig_ancom_re_drought_tolerance + theme(legend.position = "none"), align = "h", ncol = 2)
prow <- plot_grid(prow, leg_16S, ncol = 1, rel_heights = c(1, 0.38))
prow
ggsave(filename=paste("16S_ancom_re_ASV.tiff",sep=""), plot=prow, width=10, height=5, dpi=600)

# Rhizosphere -- Differentially Abundant ASVs due to Irrigation with Drought Tolerant Covariate

# Pre-Processing 
# Step 1: Data preprocessing
data_rhizo_genus <- tax_glom(data_16S_rhizosphere, "Genus")

meta_data = data_16S_rhizosphere@sam_data
meta_data <- as.data.frame(meta_data)
meta_data$names <- rownames(meta_data)
feature_table = data_16S_rhizosphere@otu_table
feature_table <- as.data.frame(feature_table)

sample_var = "names" ; group_var = NULL
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM

main_var = "irrigation"; p_adj_method = "BH"; alpha = 0.05
adj_formula = "drought_tolerance"; rand_formula = NULL
t_start = Sys.time()
res_rhizo_irrigation = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
                             alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start

# Append taxanomic information to the ANCOM Output 
re_tax <- as.data.frame(data_16S_rhizosphere@tax_table)

res_rhizo_irrigation$out <- merge(res_rhizo_irrigation$out, re_tax %>% rownames_to_column(var = "rowname"),
                                  by.x = 'taxa_id', by.y = 'rowname')

write.csv(res_rhizo_irrigation$out, "ancom_rhizo_irrigation_ASV.csv", row.names = F)
# Step 3: Volcano Plot -- WITH Taxanomic annotation added on 

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
res_rhizo_irrigation$fig$data <- merge (res_rhizo_irrigation$fig$data, re_tax %>% rownames_to_column(var = "rowname"),
                                        by.x = 'taxa_id', by.y = 'rowname')
dat_ann = data.frame(x = min(res_rhizo_irrigation$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")

fig_ancom_rhizo_irrigation = res_rhizo_irrigation$fig +  
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  geom_point(aes(color = Phylum)) + 
  scale_fill_brewer(palette="Set1") +
  ggtitle("Rhizosphere - Irrigation") +
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
fig_ancom_rhizo_irrigation

# Rhizosphere -- Differentially Abundant ASVs due to Drought Tolerance with Irrigation Covariate 

# Step 2: ANCOM

main_var = "drought_tolerance"; p_adj_method = "BH"; alpha = 0.05
adj_formula = "irrigation"; rand_formula = NULL
t_start = Sys.time()
res_rhizo_drought = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
                          alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start

# Append taxanomic information to the ANCOM Output 
re_tax <- as.data.frame(data_rhizo_genus@tax_table)

res_rhizo_drought$out <- merge(res_rhizo_drought$out, re_tax %>% rownames_to_column(var = "rowname"),
                               by.x = 'taxa_id', by.y = 'rowname')

write.csv(res_rhizo_drought$out, "ancom_rhizo_drought_tolerance_genus.csv", row.names = F)

# Step 3: Volcano Plot -- WITH Taxanomic annotation added on 

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
res_rhizo_drought$fig$data <- merge (res_rhizo_drought$fig$data, re_tax %>% rownames_to_column(var = "rowname"),
                                     by.x = 'taxa_id', by.y = 'rowname')
dat_ann = data.frame(x = min(res_rhizo_drought$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")

fig_ancom_rhizo_drought_tolerance = res_rhizo_drought$fig +  
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  geom_point(aes(color = Phylum)) + 
  scale_fill_brewer(palette="Set1") +
  ggtitle("Rhizosphere - Drought Tolerance") +
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)

fig_ancom_rhizo_drought_tolerance
prow <- plot_grid(fig_ancom_rhizo_irrigation + theme(legend.position = "none"), fig_ancom_rhizo_drought_tolerance + theme(legend.position = "none"), align = "vh", axis = "lr", ncol = 1)
prow <- plot_grid(prow, ncol = 1)


leg_16S <- get_legend(fig_ancom_rhizo_irrigation + theme(legend.position = "bottom"))
prow <- plot_grid(fig_ancom_rhizo_irrigation + theme(legend.position = "none"), fig_ancom_rhizo_drought_tolerance + theme(legend.position = "none"), align = "h", ncol = 2)
prow <- plot_grid(prow, leg_16S, ncol = 1, rel_heights = c(1, 0.50))
prow
ggsave(filename=paste("16S_ancom_rhizo_genus.tiff",sep=""), plot=prow, width=10, height=5, dpi=300)

ggsave(filename=paste("16S_ancom_rhizo.tiff",sep=""), plot=prow, width=4, height=5, dpi=300)

# Soil -- Differentially Abundant ASVs due to Irrigation with Drought Tolerant Covariate 

# Pre-Processing 
# Step 1: Data preprocessing
data_soil_genus <- tax_glom(data_16S_soil, "Genus")

meta_data = data_16S_soil@sam_data
meta_data <- as.data.frame(meta_data)
meta_data$names <- rownames(meta_data)
feature_table = data_soil_genus@otu_table
feature_table <- as.data.frame(feature_table)

sample_var = "names" ; group_var = NULL
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM

main_var = "irrigation"; p_adj_method = "BH"; alpha = 0.05
adj_formula = "drought_tolerance"; rand_formula = NULL
t_start = Sys.time()
res_soil_irrigation = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
                            alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start

# Append taxanomic information to the ANCOM Output 
re_tax <- as.data.frame(data_soil_genus@tax_table)

res_soil_irrigation$out <- merge(res_soil_irrigation$out, re_tax %>% rownames_to_column(var = "rowname"),
                                 by.x = 'taxa_id', by.y = 'rowname')

write.csv(res_soil_irrigation$out, "ancom_soil_irrigation_genus.csv", row.names = F)
# Step 3: Volcano Plot -- WITH Taxanomic annotation added on 

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
res_soil_irrigation$fig$data <- merge (res_soil_irrigation$fig$data, re_tax %>% rownames_to_column(var = "rowname"),
                                       by.x = 'taxa_id', by.y = 'rowname')

dat_ann = data.frame(x = min(res_soil_irrigation$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")

fig_ancom_soil_irrigation = ggplot(data = res_soil_irrigation$fig$data) + aes(x = x, y = y) + 
  labs(x = "CLR mean difference", y = "W statistic") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "top") +
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  geom_point(aes(color = Phylum)) + 
  scale_fill_brewer(palette="Set1") +
  ggtitle("Soil - Irrigation") +
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
fig_ancom_soil_irrigation
ggsave(filename=paste("16S_ancom_soil.tiff",sep=""), plot=fig_ancom_soil_irrigation, width=4, height=5, dpi=600)
# Soil -- Differentially Abundant ASVs due to Drought Tolerance with Irrigation Covariate 

# Step 2: ANCOM

main_var = "drought_tolerance"; p_adj_method = "BH"; alpha = 0.05
adj_formula = "irrigation"; rand_formula = NULL
t_start = Sys.time()
res_soil_drought_tolerance = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
                                   alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start

# Append taxanomic information to the ANCOM Output 
re_tax <- as.data.frame(data_soil_genus@tax_table)

res_soil_drought_tolerance$out <- merge(res_soil_drought_tolerance$out, re_tax %>% rownames_to_column(var = "rowname"),
                                        by.x = 'taxa_id', by.y = 'rowname')

write.csv(res_soil_drought_tolerance$out, "ancom_soil_drought_tolerance_genus.csv", row.names = F)

# Step 3: Volcano Plot -- WITH Taxanomic annotation added on 

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
# Add taxonomic info to the figure 
res_soil_drought_tolerance$fig$data <- merge (res_soil_drought_tolerance$fig$data, re_tax %>% rownames_to_column(var = "rowname"),
                                              by.x = 'taxa_id', by.y = 'rowname')
dat_ann = data.frame(x = min(res_soil_drought_tolerance$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")



fig_ancom_soil_drought_tolerance = res_soil_drought_tolerance$fig +  
  geom_point(aes(color = Phylum)) + 
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  scale_fill_brewer(palette="Set1") +
  ggtitle("Soil - Drought Tolerance") +
  scale_color_discrete(name = "Phylum", drop = FALSE) + 
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)

fig_ancom_soil_drought_tolerance
leg_16S <- get_legend(fig_ancom_soil_irrigation + theme(legend.position = "bottom"))
prow <- plot_grid(fig_ancom_soil_irrigation + theme(legend.position = "none"), fig_ancom_soil_drought_tolerance + theme(legend.position = "none"), align = "hh", ncol = 2)
prow <- plot_grid(prow, leg_16S, ncol = 1, rel_heights = c(1, 0.50))
prow


ggsave(filename=paste("16S_ancom_soil_genus.tiff",sep=""), plot=prow, width=10, height=5, dpi=600)

# ITS Datasets

# Root Endosphere -- Differentially Abundant ASVs due to Irrigation with Drought Tolerant Covariate 


# Pre-Processing 
# Step 1: Data preprocessing
data_RE_genus <- tax_glom(data_ITS_RE, "Genus")

meta_data = data_ITS_RE@sam_data
meta_data <- as.data.frame(meta_data)
meta_data$names <- rownames(meta_data)
feature_table = data_ITS_RE@otu_table
feature_table <- as.data.frame(feature_table)

sample_var = "names" ; group_var = NULL
out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info


# Step 2: ANCOM

main_var = "irrigation"; p_adj_method = "BH"; alpha = 0.05
adj_formula = "drought_tolerance"; rand_formula = NULL
t_start = Sys.time()
res_its_re_irrigation = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
                              alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start

# Append taxanomic information to the ANCOM Output 
re_tax <- as.data.frame(data_ITS_RE@tax_table)

res_its_re_irrigation$out <- merge(res_its_re_irrigation$out, re_tax %>% rownames_to_column(var = "rowname"),
                                   by.x = 'taxa_id', by.y = 'rowname')

write.csv(res_its_re_irrigation$out, "ancom_re_irrigation_ITS_ASV.csv", row.names = F)
# Step 3: Volcano Plot -- WITH Taxanomic annotation added on 

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
res_its_re_irrigation$fig$data <- merge (res_its_re_irrigation$fig$data, re_tax %>% rownames_to_column(var = "rowname"),
                                         by.x = 'taxa_id', by.y = 'rowname')
dat_ann = data.frame(x = min(res_its_re_irrigation$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")

fig_ancom_re_irrigation = res_its_re_irrigation$fig +  
  geom_point(aes(color = Phylum)) + 
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  scale_fill_brewer(palette="Set1") +
  ggtitle("Root Endosphere - Irrigation") +
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
fig_ancom_re_irrigation

# Root Endosphere -- Differentially Abundant ASVs due to Drought Tolerance with Irrigation Covariate 

# Step 2: ANCOM

main_var = "drought_tolerance"; p_adj_method = "BH"; alpha = 0.05
adj_formula = "irrigation"; rand_formula = NULL
t_start = Sys.time()
res_its_re_drought_tolerance = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
                                     alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start

# Append taxanomic information to the ANCOM Output 
re_tax <- as.data.frame(data_ITS_RE@tax_table)

res_its_re_drought_tolerance$out <- merge(res_its_re_drought_tolerance$out, re_tax %>% rownames_to_column(var = "rowname"),
                                          by.x = 'taxa_id', by.y = 'rowname')

write.csv(res_its_re_drought_tolerance$out, "ancom_re_drought_tolerance_ITS_ASV.csv", row.names = F)

# Step 3: Volcano Plot -- WITH Taxanomic annotation added on 

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
res_its_re_drought_tolerance$fig$data <- merge (res_its_re_drought_tolerance$fig$data, re_tax %>% rownames_to_column(var = "rowname"),
                                                by.x = 'taxa_id', by.y = 'rowname')
dat_ann = data.frame(x = min(res_its_re_drought_tolerance$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")

fig_ancom_re_drought_tolerance = res_its_re_drought_tolerance$fig +  
  geom_point(aes(color = Phylum)) + 
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  scale_fill_brewer(palette="Set1") +
  ggtitle("Root Endosphere - Drought Tolerance") +
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)

fig_ancom_re_drought_tolerance
leg_ITS <- get_legend(fig_ancom_re_irrigation + theme(legend.position = "bottom"))
prow <- plot_grid(fig_ancom_re_irrigation + theme(legend.position = "none"), fig_ancom_re_drought_tolerance + theme(legend.position = "none"), align = "h", ncol = 2)
prow <- plot_grid(prow, leg_ITS, ncol = 1, rel_heights = c(1, 0.1))

ggsave(filename=paste("ITS_ancom_re_ASV.tiff",sep=""), plot=prow, width=10, height=5, dpi=600)

# Rhizosphere -- Differentially Abundant ASVs due to Irrigation with Drought Tolerant Covariate 

# Pre-Processing 
# Step 1: Data preprocessing
data_rhizo_genus <- tax_glom(data_ITS_rhizosphere, "Genus")

meta_data = data_ITS_rhizosphere@sam_data
meta_data <- as.data.frame(meta_data)
meta_data$names <- rownames(meta_data)
feature_table = data_ITS_rhizosphere@otu_table
feature_table <- as.data.frame(feature_table)

sample_var = "names" ; group_var = NULL
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM

main_var = "irrigation"; p_adj_method = "BH"; alpha = 0.05
adj_formula = "drought_tolerance"; rand_formula = NULL
t_start = Sys.time()
res_its_rhizo_irrigation = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
                                 alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start

# Append taxanomic information to the ANCOM Output 
re_tax <- as.data.frame(data_ITS_rhizosphere@tax_table)

res_its_rhizo_irrigation$out <- merge(res_its_rhizo_irrigation$out, re_tax %>% rownames_to_column(var = "rowname"),
                                      by.x = 'taxa_id', by.y = 'rowname')

write.csv(res_its_rhizo_irrigation$out, "ancom_rhizo_irrigation_ITS_ASV.csv", row.names = F)
# Step 3: Volcano Plot -- WITH Taxanomic annotation added on 

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
res_its_rhizo_irrigation$fig$data <- merge (res_its_rhizo_irrigation$fig$data, re_tax %>% rownames_to_column(var = "rowname"),
                                            by.x = 'taxa_id', by.y = 'rowname')
dat_ann = data.frame(x = min(res_its_rhizo_irrigation$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")

fig_ancom_rhizo_irrigation = res_its_rhizo_irrigation$fig +  
  geom_point(aes(color = Phylum)) + 
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  scale_fill_brewer(palette="Set1") +
  ggtitle("Rhizosphere - Irrigation") +
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
fig_ancom_rhizo_irrigation

# Rhizosphere -- Differentially Abundant ASVs due to Drought Tolerance with Irrigation Covariate 

# Step 2: ANCOM

main_var = "drought_tolerance"; p_adj_method = "BH"; alpha = 0.05
adj_formula = "irrigation"; rand_formula = NULL
t_start = Sys.time()
res_its_rhizo_drought_tolerance = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
                                        alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start

# Append taxanomic information to the ANCOM Output 
re_tax <- as.data.frame(data_ITS_rhizosphere@tax_table)

res_its_rhizo_drought_tolerance$out <- merge(res_its_rhizo_drought_tolerance$out, re_tax %>% rownames_to_column(var = "rowname"),
                                             by.x = 'taxa_id', by.y = 'rowname')

write.csv(res_its_rhizo_drought_tolerance$out, "ancom_rhizo_drought_tolerance_ITS_ASV.csv", row.names = F)

# Step 3: Volcano Plot -- WITH Taxanomic annotation added on 

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
res_its_rhizo_drought_tolerance$fig$data <- merge (res_its_rhizo_drought_tolerance$fig$data, re_tax %>% rownames_to_column(var = "rowname"),
                                                   by.x = 'taxa_id', by.y = 'rowname')
dat_ann = data.frame(x = min(res_its_re_drought_tolerance$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")

fig_ancom_rhizo_drought_tolerance = res_its_rhizo_drought_tolerance$fig + 
  geom_point(aes(color = Phylum)) + 
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  scale_fill_brewer(palette="Set1") +
  ggtitle("Rhizosphere - Drought Tolerance") +
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)

fig_ancom_rhizo_drought_tolerance
leg_ITS <- get_legend(fig_ancom_rhizo_irrigation + theme(legend.position = "bottom"))
prow <- plot_grid(fig_ancom_rhizo_irrigation + theme(legend.position = "none"), fig_ancom_rhizo_drought_tolerance + theme(legend.position = "none"), align = "h", ncol = 2)
prow <- plot_grid(prow, leg_ITS, ncol = 1, rel_heights = c(1, 0.1))
prow

ggsave(filename=paste("ITS_ancom_rhizo_ASV.tiff",sep=""), plot=prow, width=10, height=5, dpi=600)


# Soil -- Differentially Abundant ASVs due to Irrigation with Drought Tolerant Covariate 

# Pre-Processing 
# Step 1: Data preprocessing
data_soil_genus <- tax_glom(data_ITS_soil, "Genus")

meta_data = data_ITS_soil@sam_data
meta_data <- as.data.frame(meta_data)
meta_data$names <- rownames(meta_data)
feature_table = data_ITS_soil@otu_table
feature_table <- as.data.frame(feature_table)

sample_var = "names" ; group_var = NULL
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM

main_var = "irrigation"; p_adj_method = "BH"; alpha = 0.05
adj_formula = "drought_tolerance"; rand_formula = NULL
t_start = Sys.time()
res_its_soil_irrigation = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
                                alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start

# Append taxanomic information to the ANCOM Output 
re_tax <- as.data.frame(data_ITS_soil@tax_table)

res_its_soil_irrigation$out <- merge(res_its_soil_irrigation$out, re_tax %>% rownames_to_column(var = "rowname"),
                                     by.x = 'taxa_id', by.y = 'rowname')

write.csv(res_its_soil_irrigation$out, "ancom_soil_irrigation_ITS_ASV.csv", row.names = F)
# Step 3: Volcano Plot -- WITH Taxanomic annotation added on 

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
res_its_soil_irrigation$fig$data <- merge (res_its_soil_irrigation$fig$data, re_tax %>% rownames_to_column(var = "rowname"),
                                           by.x = 'taxa_id', by.y = 'rowname')
dat_ann = data.frame(x = min(res_its_soil_irrigation$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")

fig_ancom_soil_irrigation = res_its_soil_irrigation$fig +  
  geom_point(aes(color = Phylum)) +
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  scale_fill_brewer(palette="Set1") +
  ggtitle("Soil - Irrigation") +
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
fig_ancom_soil_irrigation

# Soil -- Differentially Abundant ASVs due to Drought Tolerance with Irrigation Covariate 

# Step 2: ANCOM

main_var = "drought_tolerance"; p_adj_method = "BH"; alpha = 0.05
adj_formula = "irrigation"; rand_formula = NULL
t_start = Sys.time()
res_its_soil_drought_tolerance = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
                                       alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start

# Append taxanomic information to the ANCOM Output 
re_tax <- as.data.frame(data_ITS_soil@tax_table)

res_its_soil_drought_tolerance$out <- merge(res_its_soil_drought_tolerance$out, re_tax %>% rownames_to_column(var = "rowname"),
                                            by.x = 'taxa_id', by.y = 'rowname')

write.csv(res_its_soil_drought_tolerance$out, "ancom_soil_drought_tolerance_ASV_ITS.csv", row.names = F)

# Step 3: Volcano Plot -- WITH Taxanomic annotation added on 

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
res_its_soil_drought_tolerance$fig$data <- merge (res_its_soil_drought_tolerance$fig$data, re_tax %>% rownames_to_column(var = "rowname"),
                                                  by.x = 'taxa_id', by.y = 'rowname')
dat_ann = data.frame(x = min(res_its_soil_drought_tolerance$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")

fig_ancom_soil_drought_tolerance = res_its_soil_drought_tolerance$fig +  
  geom_point(aes(color = Phylum)) +
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  scale_fill_brewer(palette="Set1") +
  ggtitle("Soil - Drought Tolerance") +
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)

fig_ancom_soil_drought_tolerance
leg_ITS <- get_legend(fig_ancom_soil_irrigation + theme(legend.position = "bottom"))
prow <- plot_grid(fig_ancom_soil_irrigation + theme(legend.position = "none"), fig_ancom_soil_drought_tolerance + theme(legend.position = "none"), align = "h", ncol = 2)
prow <- plot_grid(prow, leg_ITS, ncol = 1, rel_heights = c(1, 0.1))
prow

ggsave(filename=paste("ITS_ancom_soil_ASV.tiff",sep=""), plot=prow, width=10, height=5, dpi=600)


# Analyze Differentially Abundant Proteobacteria, Actinobacteria, and Thaumarchaeota ASVs 
# Proteobacteria ------
# Root Endoshere

# Subset ASVs
data_16S_RE_prot <- subset(otu_table(data_16S_RE), rownames(otu_table(data_16S_RE)) %in% c('b3c195907b2347a887d116e237540feb',
                                                                                           'b3bc75b5a73133d32c26ffd8ae90d41d',
                                                                                           '2653f9a6e90d957c951070035e6edeea',
                                                                                           '47a9d10e32ffbceb163e240b65825d7e',
                                                                                           '1c2b8892a0d7f828d72808cd01cfd199',
                                                                                           'bcf0b339e6aacfad9c466f86e0c5971d',
                                                                                           'df71de0e92459aa808364311135bb8fc'))
data_16S_RE_prot <- merge_phyloseq(data_16S_RE_prot, tax_table(data_16S_RE), sample_data(data_16S_RE))


# Transform to CLR Abundance 
data_16S_RE_prot_clr <- transform(
  data_16S_RE_prot,
  transform = "clr",
  target = "OTU",
  shift = 0,
  scale = 1
)

# Export
df <- data_16S_RE_prot_clr@otu_table
df <- as.data.frame(df)
df <- t(df)
df <- merge(data_16S_RE@sam_data, df, by=0)
df <- melt(df)

df <- aggregate(df$value, list(df$irrigation,
                               df$variable), FUN=mean)

write.csv(df,"re_prot_ASVs_CLR_means.csv", row.names = FALSE)

#Rhizosphere
# Subset ASVs
data_16S_rhizo_prot <- subset(otu_table(data_16S_rhizosphere), rownames(otu_table(data_16S_rhizosphere)) %in% c('0ddcd311e02f742e2e0e61ce02cf9c29',
                                                                                                                'aee9f354c80ca7baa872c3da2fe462c2',
                                                                                                                'e1c7d97fe13e9127d225d76a1feb8c78',
                                                                                                                '95d08381352f71ef3bc48445f8ccab02',
                                                                                                                '391c94b85cd374ec9623a32e9f9b4a03',
                                                                                                                'c0d91dde8e97d4e7e3aa705c6ccb18e4',
                                                                                                                '945184b6386c192c0066e0a98a154780',
                                                                                                                '0eb88b722902316b11770922bcdaca7b',
                                                                                                                'ff836ceb31404b3a37e52a4af998742b',
                                                                                                                '91a4ee945ca39dbc22df710686730f5b',
                                                                                                                'daed5a5e05797a803338c4210f1d88c7',
                                                                                                                'df71de0e92459aa808364311135bb8fc',
                                                                                                                '40fa96a2f7fc7b939f5eb0d52d7e5b3c',
                                                                                                                'fb60fb75664eaf6a23ae2c6124fbb005',
                                                                                                                '25cb02fb17426d878afd69a79572b858',
                                                                                                                'fb5c2aa7e7a97d1ba29b2c9405cb6a41',
                                                                                                                'a70b83cbaaf606acfdb78362ffeaa506',
                                                                                                                '3e67e934cc725d230d800c4c980327c8',
                                                                                                                'f6c528c217081be9871de0fc58b89408',
                                                                                                                'f81076e084174a484f7d0c302f94437d',
                                                                                                                '3ebe761bfb1238c87195d431f41bf976',
                                                                                                                'd4cb29080e56a06ef36f99b209f2b225',
                                                                                                                'a0b1f105f827e856c5fdb5f98eaecc09',
                                                                                                                '0e11a8e894a8d9a3241d7cd7915570ac',
                                                                                                                '4bc2853b48494b7dd698fccb74a67f56',
                                                                                                                'd1b7aa5cdc783723fadc6cdc643472f0',
                                                                                                                '2653f9a6e90d957c951070035e6edeea',
                                                                                                                'caa1df79927af52cdcf97e08f266aedd',
                                                                                                                '690e2d312e3f83ab4ef9893d8c5eb3ab',
                                                                                                                '7684072284024ca2b426a6bdd04c45f9',
                                                                                                                '034e2e9539ccb9186ee9b67b0429f29f',
                                                                                                                'c252134efa44f25e744c20e9ade0781e',
                                                                                                                '1bc625bd6d3510d38f0b5b3f460bdc50',
                                                                                                                '7d9564562bb10cea1553ad90158d8545',
                                                                                                                'ec414a9986a44003df0ebd401213b7d0',
                                                                                                                'f66df1898b46126275f1585978e60df0',
                                                                                                                '809b77a90319f0cda021696d7150f1f0',
                                                                                                                '8ff49840b8cdb3c84d65399a75f2f4c3',
                                                                                                                '96799da3b007f0ef0e63b89430bfef93',
                                                                                                                'a34c67c3e970be24a82f8b4c537c08ed',
                                                                                                                'c3b8b45ccbc3b2aac60dda6ef98d2863',
                                                                                                                'f91201adfe8be6ff17177daa8bb13518',
                                                                                                                '1394f8ffb5c637943daae55332628681',
                                                                                                                'bd88cc06608d082fa8df87333f49b4e0',
                                                                                                                'e8bd5a5228cd426c496dcd50f6dad542',
                                                                                                                'bdd15f551778f16db860acaa49aa6713',
                                                                                                                'aa11b7fecbbca304e3b36cce12f5bb66',
                                                                                                                '17e562606dc747dd9f950ce27d3eb95a',
                                                                                                                'f47beff0fc3748fcf424bf33bc490cf2',
                                                                                                                '651f7e388a738bdbd1aae3674ff8aed4',
                                                                                                                'aa441ccaffa9c33d44bfaf008a3cbc22',
                                                                                                                '2f35bdf3cbfaa4cd319a51cc46764656'
))
data_16S_rhizo_prot <- merge_phyloseq(data_16S_rhizo_prot, tax_table(data_16S_rhizosphere), sample_data(data_16S_rhizosphere))

# Transform to CLR Abundance 
data_16S_rhizo_prot_clr <- transform(
  data_16S_rhizo_prot,
  transform = "clr",
  target = "OTU",
  shift = 0,
  scale = 1
)
# Export
df <- data_16S_rhizo_prot_clr@otu_table
df <- as.data.frame(df)
df <- t(df)
df <- merge(data_16S_rhizosphere@sam_data, df, by=0)
df <- melt(df)

df <- aggregate(df$value, list(df$irrigation,
                               df$variable), FUN=mean)

write.csv(df,"rhizo_prot_ASVs_CLR_means.csv", row.names = FALSE)

# Soil
data_16S_soil_prot <- subset(otu_table(data_16S_soil), rownames(otu_table(data_16S_soil)) %in% c('b58912bc0c003512959e9d3f18212929',
                                                                                                 '8d88b330419d842ed44e8e815fca2e6f',
                                                                                                 '43e2703da4799e73d1552863fd34b9bb',
                                                                                                 '0ddcd311e02f742e2e0e61ce02cf9c29',
                                                                                                 '86d6db0b6e15510fa5d27ceae8fd5e0d',
                                                                                                 '7aadb3d89e4c837cb0a4d6907b061139',
                                                                                                 '9ec6c32027e4c9e177a8f06130e4433c',
                                                                                                 'aee9f354c80ca7baa872c3da2fe462c2',
                                                                                                 '43fddf1528d4a98928fd8c3a8ac23bfd',
                                                                                                 '9a7d54b73ad99641f6bbe400eb47fc75',
                                                                                                 '7d05397e621f14fd9f7637cbd89f3c1a',
                                                                                                 'e7042f5ab73df71ce0746260c0aa4c96',
                                                                                                 '2f35bdf3cbfaa4cd319a51cc46764656',
                                                                                                 '0113160d2708229f777dad91dafe6cb1',
                                                                                                 'ec414a9986a44003df0ebd401213b7d0',
                                                                                                 '5ae4f39f43232558d2dff917f7f0600e',
                                                                                                 '53e7c08ea8b58c9459d815cde6476ac1',
                                                                                                 '3271cb570c49ad201fe3f53ac83bbb26',
                                                                                                 '9f33bc1dc3fe95d9e085919a8cb0714a',
                                                                                                 '896083bb019d1c346782d30f7ea49d86',
                                                                                                 'f9c59c858aedd9b435b458497d9d077a',
                                                                                                 '422be06a6540a3f2b64b82538808ab27',
                                                                                                 '95e64ce6249478cb61b1518726b7ce73',
                                                                                                 '23deeb5c3190a63640c425904c2a6cab',
                                                                                                 'd832a45e4d0a008d57d224460d91a738',
                                                                                                 '3d8884f3e0e3b6decd0db35b6c77ca5d',
                                                                                                 'bc2e929954e781759c0730da46614a6f',
                                                                                                 'a0ef3c27be1a25560330752609a6d5e2',
                                                                                                 'afc4ce335cfb4029c93121de9ec74cbe',
                                                                                                 '62dfdd481f414b8b1efc59c28d843fd6',
                                                                                                 '29a36f646e713667bf568f89b85329d9',
                                                                                                 'b3bc75b5a73133d32c26ffd8ae90d41d',
                                                                                                 '3f95aa1b372a657401dcde0038bfb884',
                                                                                                 'd2ee947e12c4809301a982bf5838dd83',
                                                                                                 '7b80b9fe9c9fc11b5de1f5ccbcd1c1c5',
                                                                                                 '23d450899c1481b4f8601c28fc70dcb3',
                                                                                                 'bc93367612afeca1ef58fa6242ea49fa',
                                                                                                 'ab1bd8c6c168a93d080f7ad1332a1e95',
                                                                                                 'd1e00e606bd8a22a6c20622ab4079a05'))


data_16S_soil_prot <- merge_phyloseq(data_16S_soil_prot, tax_table(data_16S_soil), sample_data(data_16S_soil))

# Transform to CLR Abundance 
data_16S_soil_prot_clr <- transform(
  data_16S_soil_prot,
  transform = "clr",
  target = "OTU",
  shift = 0,
  scale = 1
)

# Export
df <- data_16S_soil_prot_clr@otu_table
df <- as.data.frame(df)
df <- t(df)
df <- merge(data_16S_soil@sam_data, df, by=0)
df <- melt(df)

df <- aggregate(df$value, list(df$irrigation,
                               df$variable), FUN=mean)

write.csv(df,"soil_prot_ASVs_CLR_means.csv", row.names = FALSE)

# Thaumarchaeota ------
# Rhizosphere
# Subset ASVs
data_16S_rhizo_thaum <- subset(otu_table(data_16S_rhizosphere), rownames(otu_table(data_16S_rhizosphere)) %in% c('9d7b1963f8036854dd7b9aad0149310f',
                                                                                                                 'c9d59b5dca26cd50adf46d0171c1e087',
                                                                                                                 '0f8ac0d81dd6449a43b34fcb184efd8f'))
data_16S_rhizo_thaum <- merge_phyloseq(data_16S_rhizo_thaum, tax_table(data_16S_rhizosphere), sample_data(data_16S_rhizosphere))

# Transform to CLR Abundance 
data_16S_rhizo_thaum_clr <- transform(
  data_16S_rhizo_thaum,
  transform = "clr",
  target = "OTU",
  shift = 0,
  scale = 1
)
# Export
df <- data_16S_rhizo_thaum_clr@otu_table
df <- as.data.frame(df)
df <- t(df)
df <- merge(data_16S_rhizosphere@sam_data, df, by=0)
df <- melt(df)

df <- aggregate(df$value, list(df$irrigation,
                               df$variable), FUN=mean)

write.csv(df,"rhizo_thaum_ASVs_CLR_means.csv", row.names = FALSE)

#Soil
data_16S_soil_thaum <- subset(otu_table(data_16S_soil), rownames(otu_table(data_16S_soil)) %in% c('2cb07dfb5ce307816d50474901ec63d2',
                                                                                                  '94e44251662880ffc5b82a1aca4ff1d8'))
data_16S_soil_thaum <- merge_phyloseq(data_16S_soil_thaum, tax_table(data_16S_soil), sample_data(data_16S_soil))

# Transform to CLR Abundance 
data_16S_soil_thaum_clr <- transform(
  data_16S_soil_thaum,
  transform = "clr",
  target = "OTU",
  shift = 0,
  scale = 1
)
# Export
df <- data_16S_soil_thaum_clr@otu_table
df <- as.data.frame(df)
df <- t(df)
df <- merge(data_16S_soil@sam_data, df, by=0)
df <- melt(df)

df <- aggregate(df$value, list(df$irrigation,
                               df$variable), FUN=mean)

write.csv(df,"soil_thaum_ASVs_CLR_means.csv", row.names = FALSE)

# Actinobacteria ------
# Root Endosphere
# Subset ASVs
data_16S_RE_actin <- subset(otu_table(data_16S_RE), rownames(otu_table(data_16S_RE)) %in% c('e984a619f22fe030e13c709b0489b52a',
                                                                                            '754c7c6f3b6988d240275e2f44436292'))

data_16S_RE_actin <- merge_phyloseq(data_16S_RE_actin, tax_table(data_16S_RE), sample_data(data_16S_RE))


# Transform to CLR Abundance 
data_16S_RE_actin_clr <- transform(
  data_16S_RE_actin,
  transform = "clr",
  target = "OTU",
  shift = 0,
  scale = 1
)

# Export
df <- data_16S_RE_actin_clr@otu_table
df <- as.data.frame(df)
df <- t(df)
df <- merge(data_16S_RE@sam_data, df, by=0)
df <- melt(df)

df <- aggregate(df$value, list(df$irrigation,
                               df$variable), FUN=mean)

write.csv(df,"re_actin_ASVs_CLR_means.csv", row.names = FALSE)

# Rhizosphere
# Subset ASVs
data_16S_rhizo_actin <- subset(otu_table(data_16S_rhizosphere), rownames(otu_table(data_16S_rhizosphere)) %in% c('572b96fefec73aa11fc6628597853ed7',
                                                                                                                 'c7c40b0f480f3edee216585ec048011d',
                                                                                                                 '2ee683d3ac3574f69eee1a2a1cd89d09',
                                                                                                                 'a1ed615a4ed0b08c534fcc504a4608b2',
                                                                                                                 'b4b365be5234d748dd5812df4ac2c6de',
                                                                                                                 '7df299814b7a37f021b9bc3c3fc40080',
                                                                                                                 'e4e08e332b44c5b0f3c79143fb9c0f63',
                                                                                                                 'ddfaffe2f2f94985e611459236b87647',
                                                                                                                 '1cf8836d9f2282cac057b416706ba561',
                                                                                                                 'cfdf66ad8b27626f14f419a4f6c51b88',
                                                                                                                 'fe170a37c40648f3cd117b40946280fc',
                                                                                                                 'b48a188a24ed126894db7d73bd6617cc',
                                                                                                                 '2e6bc9a35299587b295c308e5f616477',
                                                                                                                 'e07052bb9c0f0d2370537c4b9acbed4e',
                                                                                                                 '50fb93dad525b440dc193bc631872eba',
                                                                                                                 '163076c2d7e10152bb3d0e2f4b28df75',
                                                                                                                 'c4def24d4e9cea12f3f0a8a0ca43afda',
                                                                                                                 '333048cac277d86e040b020b11c3eabf',
                                                                                                                 '875dc0da5f66b1b5644788543a07544f',
                                                                                                                 '7b7f9c073366c233fa85c7a9aee5903c',
                                                                                                                 '13645da5ce1259a647f3bc11b09dcf13',
                                                                                                                 'd51c399452cc3e608552d0aca841b28c',
                                                                                                                 '3a54d4e72af90b9c2fd50d7d1bb0336b',
                                                                                                                 '7034493f95fe775db5f9eeb1aac0cad9',
                                                                                                                 '594061ed4e2ea378773c474cd4803c29',
                                                                                                                 'eea13d94a260c6a4c38e783eb79c68d0',
                                                                                                                 'b9875bae83467e83cbb87f0bb6ed2db0',
                                                                                                                 'd2d78251d7df8d1a835217535550cb59',
                                                                                                                 'e3bbf1bd6c612fbffbbf88f4dc54b805',
                                                                                                                 'd19dbb0e1566e9a2d91adda6f64e0232',
                                                                                                                 'c51c974c591cce75ceff90c9bd817424',
                                                                                                                 '002664222d5c2d5b79c41b33598ea1e4',
                                                                                                                 '52fd95eb476921cf2b9510888ea0cfe8',
                                                                                                                 'd892d40bc79d0c356bd2886ce8d120ee',
                                                                                                                 '4922eb95c3994f40d0784f2cfbf7ede5',
                                                                                                                 '3bed507c2cb7458e3bd9b0ee3692fecb'))
data_16S_rhizo_actin <- merge_phyloseq(data_16S_rhizo_actin, tax_table(data_16S_rhizosphere), sample_data(data_16S_rhizosphere))


# Transform to CLR Abundance 
data_16S_rhizo_actin_clr <- transform(
  data_16S_rhizo_actin,
  transform = "clr",
  target = "OTU",
  shift = 0,
  scale = 1
)

# Export
df <- data_16S_rhizo_actin_clr@otu_table
df <- as.data.frame(df)
df <- t(df)
df <- merge(data_16S_rhizosphere@sam_data, df, by=0)
df <- melt(df)

df <- aggregate(df$value, list(df$irrigation,
                               df$variable), FUN=mean)

write.csv(df,"rhizo_actin_ASVs_CLR_means.csv", row.names = FALSE)

#Bulk Soil
# Subset ASVs
data_16S_soil_actin <- subset(otu_table(data_16S_soil), rownames(otu_table(data_16S_soil)) %in% c('7c9260cb39c10fef45b09ad77d5adb8f',
                                                                                                  'b8f4b792a7569c5f07c33e036d449dc4',
                                                                                                  'fe170a37c40648f3cd117b40946280fc',
                                                                                                  '4f46fb48db4a759699801062a06b93fc',
                                                                                                  '790ae917c0105961dfa23507f17853dd',
                                                                                                  'd19dbb0e1566e9a2d91adda6f64e0232',
                                                                                                  '2ee683d3ac3574f69eee1a2a1cd89d09',
                                                                                                  '50fb93dad525b440dc193bc631872eba',
                                                                                                  'c7c40b0f480f3edee216585ec048011d',
                                                                                                  'c6aac0ac5398c1e9a6503ec5c65f13c1',
                                                                                                  '7286a9b89b5bb4cfd7046f878ea2feab',
                                                                                                  'ef0c7928901deab56b8d8abef337ac5c',
                                                                                                  'c8ad5c63b9707523790d4ae35b5caa54',
                                                                                                  '09de3dde6c465b74bcb8e6f8a478e4e0',
                                                                                                  '3a54d4e72af90b9c2fd50d7d1bb0336b',
                                                                                                  '1c456aadfa289b274ea023f362bbbd30',
                                                                                                  'ee6ee5f00eb1c192da3115ea88c2c0b1',
                                                                                                  '8a089f32ad8b77b8ba5df3d9a8d3ba04',
                                                                                                  '0708edc072e4349909322f96b7b9df4e',
                                                                                                  'c5079e844cab8cd7c10cf5328cfa1c8d',
                                                                                                  '9ea1d464f26978844f7277e91ef53c10',
                                                                                                  '8e7c0224d0b63a3d52e5886b65020a83',
                                                                                                  '94c7a1a05aaac0b19c8cd35e0b430b90'))

data_16S_soil_actin <- merge_phyloseq(data_16S_soil_actin, tax_table(data_16S_soil), sample_data(data_16S_soil))


# Transform to CLR Abundance 
data_16S_soil_actin_clr <- transform(
  data_16S_soil_actin,
  transform = "clr",
  target = "OTU",
  shift = 0,
  scale = 1
)

# Export
df <- data_16S_soil_actin_clr@otu_table
df <- as.data.frame(df)
df <- t(df)
df <- merge(data_16S_soil@sam_data, df, by=0)
df <- melt(df)

df <- aggregate(df$value, list(df$irrigation,
                               df$variable), FUN=mean)

write.csv(df,"soil_actin_ASVs_CLR_means.csv", row.names = FALSE)

# Boardman Genotype Height Data Analysis ---------------------------------------------------------------------
library(readxl)
Genotype_Height_Data <- read_excel("Genotype_Height_Data.xlsx")

my_comparisons <- list(c("HI30", "LO30"))
Genotype_Height_Data %>%
  ggplot(aes(x=`Drought Tolerance`, y=`TPS Corrected Height`, colour = `Drought Tolerance`)) +
  scale_color_manual("Drought Tolerance",
                     values = c("tomato1","steelblue4"),
                     labels = c(expression(paste("HI30", sep="")),
                                expression(paste("LO30", sep="")))) +
  facet_wrap(~Irrigation) +
  geom_boxplot(alpha=0.5) + 
  geom_jitter(shape = 16, size =2, position = position_jitter(0.2)) +
  scale_y_continuous(name = "Average Height") +
  stat_compare_means(comparisons = my_comparisons, method="t.test",
                     label.x = 1.5) +
  theme_bw() +
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_blank()) -> p_Genotype
ggsave("Genotype_Height.tiff", plot = p_Genotype, height = 5, width = 7.5, dpi = 600)




