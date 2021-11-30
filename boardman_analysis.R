# Load Required Packages

#library(ALDEx2)
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
#library(SpiecEasi)
library(tidyverse)
library(vegan)
library(zCompositions)
#source("C:/Users/Brand/OneDrive/Documents/R/win-library/4.1/ANCOM-master/scripts/ancom_v2.1.R")

#Set-up ------------------------------------------------------------------------
# Set Working Directory

setwd("~/Dropbox/Cregger_CBI_drought/data/")

# Read QZA files into dataframe, re-format taxonomic tables, and re-upload them as .csv files

SVs16S <- read_qza("qiime_output/16S-merged-table.qza")
SVs16Stable <- SVs16S$data
write.csv(SVs16Stable, file = "16S-merged-table.csv")
taxonomy16S <- read_qza("qiime_output/16S-taxonomy.qza")
tax16S<-taxonomy16S$data %>% as.tibble() %>%
  mutate(Taxon=gsub("D_[0-9]__", "", Taxon)) %>% 
  separate(Taxon, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))%>%
  mutate(Phylum=replace_na(Phylum,"empty"))
write.csv(tax16S, file = "16S-taxonomy.csv", row.names =F)

SVsITS <- read_qza("qiime_output/ITS-merged-table.qza")
SVsITStable <- SVsITS$data
write.csv(SVsITStable, file = "ITS-merged-table.csv")
taxonomyITS <- read_qza("qiime_output/ITS-taxonomy.qza")
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
  scale_color_manual("Plant Compartment",
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
  scale_color_manual("Plant Compartment",
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
ggsave(filename=paste("rarefaction_curves.tiff", sep=""), plot=prow, width = 8.5, height=6.5, dpi=600)

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
        scale_color_manual("Plant Compartment",
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
  scale_color_manual("Plant Compartment",
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
  scale_color_manual("Plant Compartment",
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
  scale_color_manual("Plant Compartment",
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
prow3 <- plot_grid(prow1, prow2, leg, ncol = 3, rel_widths = c(1, 1,0.65))
prow3 
ggsave(filename=paste("beta_diversity_habitat.tiff", sep=""), plot=prow3, width=6.5, height=6.5, dpi=600)
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
rare_16S_RE <- rarefy_even_depth(data_16S_RE, sample.size = 14000)
ntaxa(rare_16S_RE)

set.seed(01221990)
rare_16S_rhizosphere <- rarefy_even_depth(data_16S_rhizosphere, sample.size = 25000)
ntaxa(rare_16S_rhizosphere)


set.seed(01221990)
rare_16S_soil <- rarefy_even_depth(data_16S_soil, sample.size = 45000)
ntaxa(rare_16S_soil)


set.seed(01221990)
rare_ITS_RE <- rarefy_even_depth(data_ITS_RE, sample.size = 820)
ntaxa(rare_ITS_RE)

set.seed(01221990)
rare_ITS_rhizosphere <- rarefy_even_depth(data_ITS_rhizosphere, sample.size = 1300)
ntaxa(rare_ITS_rhizosphere)


set.seed(01221990)
rare_ITS_soil <- rarefy_even_depth(data_ITS_soil, sample.size = 10000)
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

# Alpha Diversity Plot
# Root Endosphere
rich_16S_re$hilln <- "richness"
shan_16S_re$hilln <- "shannon"
simp_16S_re$hilln <- "simpson"
re <- rbind(rich_16S_re %>% dplyr::select(irrigation, y, drought_tolerance, hilln),
            shan_16S_re %>% dplyr::select(irrigation, y, drought_tolerance, hilln),
            simp_16S_re %>% dplyr::select(irrigation, y, drought_tolerance, hilln))
vnames <- list(
  "richness" = bquote(~phantom()^0*italic(D)), 
  "shannon" = bquote(~phantom()^1*italic(D)), 
  "simpson" = bquote(~phantom()^2*italic(D))
)
vlabeller <- function(variable, value){return(vnames[value])}
re %>%
  ggplot(aes(x = irrigation, y = y, fill = drought_tolerance)) +
  geom_boxplot() +
  facet_grid(rows = vars(hilln), 
             labeller = vlabeller,
             scales = "free") +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) +
  scale_y_continuous(name = "Archaeal/Bacterial Alpha Diversity") +
  scale_fill_manual("Drought Tolerance",
                    values = c("tomato1","steelblue4"),
                    labels = c(expression(paste("Drought Tolerant", sep="")),
                               expression(paste("Drought Susceptible", sep="")))) +
  panel_border() +
  theme_bw() +
  ggtitle(label = "Root Endosphere") +
  theme(plot.title = element_text(size = 9.5)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        legend.title = element_text(size = 12)) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted") -> p_re
p_re

# Rhizosphere 
rich_16S_rhizosphere$hilln <- "richness"
shan_16S_rhizosphere$hilln <- "shannon"
simp_16S_rhizosphere$hilln <- "simpson"
rhizo <- rbind(rich_16S_rhizosphere %>% dplyr::select(irrigation, y, drought_tolerance, hilln),
               shan_16S_rhizosphere %>% dplyr::select(irrigation, y, drought_tolerance, hilln),
               simp_16S_rhizosphere %>% dplyr::select(irrigation, y, drought_tolerance, hilln))

rhizo %>%
  ggplot(aes(x = irrigation, y = y, fill = drought_tolerance)) +
  geom_boxplot() +
  facet_grid(rows = vars(hilln), 
             labeller = vlabeller,
             scales = "free") +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) +
  scale_y_continuous(name = "Archaeal/Bacterial Alpha Diversity") +
  scale_fill_manual("Drought TOlerance",
                    values = c("tomato1", "steelblue4"),
                    labels = c(expression(paste("Drought Tolerant", sep="")),
                               expression(paste("Drought Susceptible", sep="")))) +
  panel_border() +
  theme_bw() +
  ggtitle(label = "Rhizosphere") +
  theme(plot.title = element_text(size = 9.5)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        legend.title = element_text(size = 12)) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted") -> p_rhizo
p_rhizo

# Bulk Soil
rich_16S_soil$hilln <- "richness"
shan_16S_soil$hilln <- "shannon"
simp_16S_soil$hilln <- "simpson"
soil <- rbind(rich_16S_soil %>% dplyr::select(irrigation, y, drought_tolerance, hilln),
              shan_16S_soil %>% dplyr::select(irrigation, y, drought_tolerance, hilln),
              simp_16S_soil %>% dplyr::select(irrigation, y, drought_tolerance, hilln))

soil %>%
  ggplot(aes(x = irrigation, y = y, fill = drought_tolerance)) +
  geom_boxplot() +
  facet_grid(rows = vars(hilln), 
             labeller = vlabeller,
             scales = "free") +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) +
  scale_y_continuous(name = "Archaeal/Bacterial Alpha Diversity") +
  scale_fill_manual("Drought Tolerance",
                    values = c("tomato1", "steelblue4"),
                    labels = c(expression(paste("Drought Tolerant", sep="")),
                               expression(paste("Drought Susceptible", sep="")))) +
  panel_border() +
  theme_bw() +
  ggtitle(label = "Bulk Soil") +
  theme(plot.title = element_text(size = 9.5)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        legend.title = element_text(size = 12)) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted") -> p_soil
p_soil

# Appended 16S Alpha Diversity Plot
leg_16S <- get_legend(p_re)
prow_16S <- plot_grid(p_re + theme(legend.position = "none"), p_rhizo + theme(legend.position = "none"),  p_soil + theme(legend.position = "none"),
                       ncol = 3, labels = "AUTO", align = "vh")
alpha_16S <- plot_grid(prow_16S, leg_16S, ncol = 2, rel_widths = c(1, 0.25))
alpha_16S
ggsave(filename=paste("16S_alpha_diversity_treatments.tiff",sep=""), plot=alpha_16S, width=8.4, height=7, dpi=600)

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

# Root Endosphere

rich_ITS_re$hilln <- "richness"
shan_ITS_re$hilln <- "shannon"
simp_ITS_re$hilln <- "simpson"
re <- rbind(rich_ITS_re %>% dplyr::select(irrigation, y, drought_tolerance, hilln),
            shan_ITS_re %>% dplyr::select(irrigation, y, drought_tolerance, hilln),
            simp_ITS_re %>% dplyr::select(irrigation, y, drought_tolerance, hilln))
vnames <- list(
  "richness" = bquote(~phantom()^0*italic(D)), 
  "shannon" = bquote(~phantom()^1*italic(D)), 
  "simpson" = bquote(~phantom()^2*italic(D))
)
vlabeller <- function(variable, value){return(vnames[value])}
re %>%
  ggplot(aes(x = irrigation, y = y, fill = drought_tolerance)) +
  geom_boxplot() +
  facet_grid(rows = vars(hilln), 
             labeller = vlabeller,
             scales = "free") +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) +
  scale_y_continuous(name = "Fungal Alpha Diversity") +
  scale_fill_manual("Drought Tolerance",
                    values = c("tomato1", "steelblue4"),
                    labels = c(expression(paste("Drought Tolerant", sep="")),
                               expression(paste("Drought Susceptible", sep="")))) +
  panel_border() +
  theme_bw() +
  ggtitle(label = "Root Endosphere") +
  theme(plot.title = element_text(size = 9.5)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        legend.title = element_text(size = 12)) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted") -> p_re
p_re

# Rhizosphere

rich_ITS_rhizosphere$hilln <- "richness"
shan_ITS_rhizosphere$hilln <- "shannon"
simp_ITS_rhizosphere$hilln <- "simpson"
rhizo <- rbind(rich_ITS_rhizosphere %>% dplyr::select(irrigation, y, drought_tolerance, hilln),
               shan_ITS_rhizosphere %>% dplyr::select(irrigation, y, drought_tolerance, hilln),
               simp_ITS_rhizosphere %>% dplyr::select(irrigation, y, drought_tolerance, hilln))
rhizo %>%
  ggplot(aes(x = irrigation, y = y, fill = drought_tolerance)) +
  geom_boxplot() +
  facet_grid(rows = vars(hilln), 
             labeller = vlabeller,
             scales = "free") +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) +
  scale_y_continuous(name = "Fungal Alpha Diversity") +
  scale_fill_manual("Drought Tolerance",
                    values = c("tomato1", "steelblue4"),
                    labels = c(expression(paste("Drought Tolerant", sep="")),
                               expression(paste("Drought Susceptible", sep="")))) +
  panel_border() +
  theme_bw() +
  ggtitle(label = "Rhizosphere") +
  theme(plot.title = element_text(size = 9.5)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        legend.title = element_text(size = 12)) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted") -> p_rhizo
p_rhizo
ggsave(filename="slide_27.tiff", plot=p_rhizo, height=5.5, width=6.5, dpi=300)

# Soil

rich_ITS_soil$hilln <- "richness"
shan_ITS_soil$hilln <- "shannon"
simp_ITS_soil$hilln <- "simpson"
soil <- rbind(rich_ITS_soil %>% dplyr::select(irrigation, y, drought_tolerance, hilln),
               shan_ITS_soil %>% dplyr::select(irrigation, y, drought_tolerance, hilln),
               simp_ITS_soil %>% dplyr::select(irrigation, y, drought_tolerance, hilln))

soil %>%
  ggplot(aes(x = irrigation, y = y, fill = drought_tolerance)) +
  geom_boxplot() +
  facet_grid(rows = vars(hilln), 
             labeller = vlabeller,
             scales = "free") +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) +
  scale_y_continuous(name = "Fungal Alpha Diversity") +
  scale_fill_manual("Drought Tolerance",
                    values = c("tomato1", "steelblue4"),
                    labels = c(expression(paste("Drought Tolerant", sep="")),
                               expression(paste("Drought Susceptible", sep="")))) +
  panel_border() +
  theme_bw() +
  ggtitle(label = "Bulk Soil") +
  theme(plot.title = element_text(size = 9.5)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        legend.title = element_text(size = 12)) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted") -> p_soil
p_soil

#Appended ITS Alpha Diversity Plot
leg_ITS <- get_legend(p_re)
prow_ITS <- plot_grid(p_re + theme(legend.position = "none"), p_rhizo + theme(legend.position = "none"),  p_soil + theme(legend.position = "none"),
                      ncol = 3, labels = "AUTO", align = "vh")
alpha_ITS <- plot_grid(prow_ITS, leg_ITS, ncol = 2, rel_widths = c(1, 0.25))


ggsave(filename=paste("ITS_alpha_diversity_treatments.tiff",sep=""), plot=alpha_ITS, width=8.4, height=7, dpi=600)


# Beta Diversity - Treatment ---------------------------------------------------

# 16S Library -- Subsetted by Habitat


#Root endosphere

#dbRDA Plot for the Root Endosphere
# Create a column with appended irrigation.drought tolerance 
df_16S <- data.frame(rare_16S_RE@sam_data)
df_16S$irrigation.drought_tolerance <- with(rare_16S_RE@sam_data, interaction(irrigation, drought_tolerance))
OTU_16S <- t(as(rare_16S_RE@otu_table, "matrix"))
OTU_16S <- as.data.frame(OTU_16S)
dist_16S <- vegdist(OTU_16S, "jaccard")
set.seed(1221990)
rda_16S <- dbrda(dist_16S ~ irrigation + drought_tolerance + irrigation*drought_tolerance, df_16S)

anova(rda_16S) 
print(anova(rda_16S, by="terms", permutation = 9999))
plot(rda_16S)

y <- varpart(dist_16S, ~irrigation, ~drought_tolerance, data = df_16S)
print(y)

# Measure for dispersion
dispersion_16S <- betadisper(dist_16S, df_16S$irrigation.drought_tolerance, type = c("centroid"), bias.adjust = TRUE, sqrt.dist = FALSE, add = FALSE)
anova(dispersion_16S)

set.seed(1221990)
permutest(dispersion_16S, pairwise = TRUE, permutations = 9999)
(dispersion_16S_HSD <- TukeyHSD(dispersion_16S))

# Plot Dispersion

data <- data.frame(Distance_to_centroid=dispersion_16S$distances, Group=dispersion_16S$group)
data$irrigation <- df_16S$irrigation
data$drought_tolerance <- df_16S$drought_tolerance

groups <- dispersion_16S$group
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
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted") +
  ggtitle("Archaeal/Bacterial Root Endosphere Dispersion") +
  panel_border() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12))




plot.data <- merge(summary(rda_16S)$sites, df_16S, by = "row.names") 
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
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_16S)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_16S)$cont)[2,2], 1)), "%)")) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) -> p_jaccard_re_16S

p_jaccard_re_16S

ggsave(filename="16S_dbrda_jaccard_re.tiff", plot=p_jaccard_re, height=5.5, width=6.5, dpi=300)

#Rhizosphere

#dbRDA Plot for the Rhizosphere
df_16S <- data.frame(rare_16S_rhizosphere@sam_data)
df_16S$irrigation.drought_tolerance <- with(rare_16S_rhizosphere@sam_data, interaction(irrigation, drought_tolerance))
OTU_16S <- t(as(rare_16S_rhizosphere@otu_table, "matrix"))
OTU_16S <- as.data.frame(OTU_16S)
dist_16S <- vegdist(OTU_16S, "jaccard")
set.seed(1221990)
rda_16S <- dbrda(dist_16S ~ irrigation + drought_tolerance + irrigation*drought_tolerance, df_16S)

anova(rda_16S) # Overall model is significant (p=0.001)
print(anova(rda_16S, by="terms", permutation = 9999))
plot(rda_16S)

y <- varpart(dist_16S, ~irrigation, ~drought_tolerance, data = df_16S)
print(y)

# Measure for dispersion
dispersion_16S <- betadisper(dist_16S, df_16S$irrigation.drought_tolerance, type = c("centroid"), bias.adjust = TRUE, sqrt.dist = FALSE, add = FALSE)
anova(dispersion_16S)

set.seed(1221990)
permutest(dispersion_16S, pairwise = TRUE, permutations = 9999)
(dispersion_16S_HSD <- TukeyHSD(dispersion_16S))

# Plot Dispersion

data <- data.frame(Distance_to_centroid=dispersion_16S$distances, Group=dispersion_16S$group)
data$irrigation <- df_16S$irrigation
data$drought_tolerance <- df_16S$drought_tolerance

groups <- dispersion_16S$group
disper_16S_rhizo <- ggplot(data=data, aes(x= irrigation, y=Distance_to_centroid, shape = irrigation, fill = drought_tolerance)) + 
  geom_jitter(size =3.5, position = position_jitter(0.1)) +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) + 
  scale_y_continuous(name = "Distance to Centroid") +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted") +
  ggtitle("Archaeal/Bacterial Rhizosphere Dispersion") +
  panel_border() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        plot.title = element_text(size=12,  face="bold"),
        legend.title = element_text(size = 12))


plot.data <- merge(summary(rda_16S)$sites, df_16S, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, shape = irrigation, fill = drought_tolerance), 
             size = 3.5, alpha = 1) +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression(italic("Full")), 
                                expression(italic("Reduced")))) +
  ggtitle("Archaeal/Bacterial Rhizosphere") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_16S)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_16S)$cont)[2,2], 1)), "%)")) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        plot.title = element_text(size=12,  face="bold"),
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) -> p_jaccard_rhizo_16S

p_jaccard_rhizo_16S

#Bulk Soil

#dbRDA Plot for the Soil

df_16S <- data.frame(rare_16S_soil@sam_data)
df_16S$irrigation.drought_tolerance <- with(rare_16S_soil@sam_data, interaction(irrigation, drought_tolerance))
OTU_16S <- t(as(rare_16S_soil@otu_table, "matrix"))
OTU_16S <- as.data.frame(OTU_16S)
dist_16S <- vegdist(OTU_16S, "jaccard")
set.seed(1221990)
rda_16S <- dbrda(dist_16S ~ irrigation + drought_tolerance + irrigation*drought_tolerance, df_16S)

anova(rda_16S) # Overall model is significant (p=0.001)
print(anova(rda_16S, by="terms", permutation = 9999))
plot(rda_16S)

y <- varpart(dist_16S, ~irrigation, ~drought_tolerance, data = df_16S)
print(y)

# Measure for dispersion
dispersion_16S <- betadisper(dist_16S, df_16S$irrigation.drought_tolerance, type = c("centroid"), bias.adjust = TRUE, sqrt.dist = FALSE, add = FALSE)
anova(dispersion_16S)

set.seed(1221990)
permutest(dispersion_16S, pairwise = TRUE, permutations = 9999)
(dispersion_16S_HSD <- TukeyHSD(dispersion_16S))

# Plot Dispersion

data <- data.frame(Distance_to_centroid=dispersion_16S$distances, Group=dispersion_16S$group)
data$irrigation <- df_16S$irrigation
data$drought_tolerance <- df_16S$drought_tolerance

groups <- dispersion_16S$group
disper_16S_soil <- ggplot(data=data, aes(x= irrigation, y=Distance_to_centroid, shape = irrigation, fill = drought_tolerance)) + 
  geom_jitter(size =3.5, position = position_jitter(0.1)) +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) + 
  scale_y_continuous(name = "Distance to Centroid") +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted") +
  ggtitle("Archaeal/Bacterial Bulk Soil Dispersion") +
  panel_border() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12))




plot.data <- merge(summary(rda_16S)$sites, df_16S, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, shape = irrigation, fill = drought_tolerance), 
             size = 3.5, alpha = 1) +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression (italic("Reduced")))) +
  ggtitle("Archaeal/Bacterial Bulk Soil") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_16S)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_16S)$cont)[2,2], 1)), "%)")) +
  theme_bw() +
  #geom_text(data = as.data.frame(summary(rda_ITS)$biplot) %>% rownames_to_column(), aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        plot.title = element_text(size=12,  face="bold"),
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) -> p_jaccard_soil_16S

p_jaccard_soil_16S

ggsave(filename="16S_dbrda_jaccard_soil.tiff", plot=p_jaccard_soil, height=5.5, width=6.5, dpi=300)

#Append Plots Together
 prow <- plot_grid(p_jaccard_re_16S + theme(legend.position = "none"), disper_16S_re + theme(legend.position = "none"),
                  p_jaccard_rhizo_16S + theme(legend.position = "none"), disper_16S_rhizo + theme(legend.position = "none"),
                  p_jaccard_soil_16S + theme(legend.position="none"), disper_16S_soil + theme(legend.position = "none"),
                  ncol = 2, labels = "AUTO", align = "vh")
leg <- get_legend(p_jaccard_re_16S)
prow <- plot_grid(prow, leg, ncol = 2, rel_widths = c(0.8, 0.23))
ggsave(filename=paste("16S_treatment_beta_diversity.tiff", sep=""), plot=prow, width=7.5, height=7.5, dpi=600)

  # ITS Library -- Subsetted by Habitat 
# Loop through each distance method, save each plot to a list called plist

#Root endosphere
#dbRDA Plot for the Root Endosphere

df_ITS <- data.frame(rare_ITS_RE@sam_data)
df_ITS$irrigation.drought_tolerance <- with(rare_ITS_RE@sam_data, interaction(irrigation, drought_tolerance))
OTU_ITS <- t(as(rare_ITS_RE@otu_table, "matrix"))
OTU_ITS <- as.data.frame(OTU_ITS)
dist_ITS <- vegdist(OTU_ITS, "jaccard")
set.seed(1221990)
rda_ITS <- dbrda(dist_ITS ~ irrigation + drought_tolerance + irrigation*drought_tolerance, df_ITS)

anova(rda_ITS) # Overall model is significant (p=0.036)
anova(rda_ITS, by="terms", permutation = 9999) # Drought tolerance influenced microbiome structure of the root endosphere (p=0.001);

y <- varpart(dist_ITS, ~irrigation, ~drought_tolerance, data = df_ITS)
print(y)

# Measure for dispersion
dispersion_ITS <- betadisper(dist_ITS, df_ITS$irrigation.drought_tolerance, type = c("centroid"), bias.adjust = TRUE, sqrt.dist = FALSE, add = FALSE)
anova(dispersion_ITS)

set.seed(1221990)
permutest(dispersion_ITS, pairwise = TRUE, permutations = 9999)
(dispersion_ITS_HSD <- TukeyHSD(dispersion_ITS))

# Plot Dispersion

data <- data.frame(Distance_to_centroid=dispersion_ITS$distances, Group=dispersion_ITS$group)
data$irrigation <- df_ITS$irrigation
data$drought_tolerance <- df_ITS$drought_tolerance

groups <- dispersion_ITS$group
disper_ITS_re <- ggplot(data=data, aes(x= irrigation, y=Distance_to_centroid, shape = irrigation, fill = drought_tolerance)) + 
  geom_jitter(size =3.5, position = position_jitter(0.1)) +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) + 
  scale_y_continuous(name = "Distance to Centroid") +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted") +
  ggtitle("Root Endosphere") +
  panel_border() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        plot.title = element_text(size=12,  face="bold"),
        legend.title = element_text(size = 12))

plot.data <- merge(summary(rda_ITS)$sites, df_ITS, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, shape = irrigation, fill = drought_tolerance), 
             size = 4, alpha = 1) +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression(italic("Full")), 
                                expression(italic("Reduced")))) +
  ggtitle("Root Endosphere") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_ITS)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_ITS)$cont)[2,2], 1)), "%)")) +
  theme_bw() +
  #geom_text(data = as.data.frame(summary(rda_ITS)$biplot) %>% rownames_to_column(), aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        plot.title = element_text(size=12,  face="bold"),
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) -> p_jaccard_re_ITS

p_jaccard_re_ITS

ggsave(filename="ITS_dbrda_jaccard_re.tiff", plot=p_jaccard_re, height=5.5, width=6.5, dpi=300)

#Rhizosphere
#dbRDA Plot for the Rhizosphere

df_ITS <- data.frame(rare_ITS_rhizosphere@sam_data)
df_ITS$irrigation.drought_tolerance <- with(rare_ITS_rhizosphere@sam_data, interaction(irrigation, drought_tolerance))
OTU_ITS <- t(as(rare_ITS_rhizosphere@otu_table, "matrix"))
OTU_ITS <- as.data.frame(OTU_ITS)
dist_ITS <- vegdist(OTU_ITS, "jaccard")
set.seed(1221990)
rda_ITS <- dbrda(dist_ITS ~ irrigation + drought_tolerance + irrigation*drought_tolerance, df_ITS)

anova(rda_ITS) # Overall model is significant (p=0.005)
anova(rda_ITS, by="terms", permutation = 9999) # Irrigation influenced microbiome structure of the rhizosphere (p=0.004)
plot(rda_ITS)

y <- varpart(dist_ITS, ~irrigation, ~drought_tolerance, data = df_ITS)
print(y)


# Measure for dispersion
dispersion_ITS <- betadisper(dist_ITS, df_ITS$irrigation.drought_tolerance, type = c("centroid"), bias.adjust = TRUE, sqrt.dist = FALSE, add = FALSE)
anova(dispersion_ITS)

set.seed(1221990)
permutest(dispersion_ITS, pairwise = TRUE, permutations = 9999)
(dispersion_ITS_HSD <- TukeyHSD(dispersion_ITS))

# Plot Dispersion

data <- data.frame(Distance_to_centroid=dispersion_ITS$distances, Group=dispersion_ITS$group)
data$irrigation <- df_ITS$irrigation
data$drought_tolerance <- df_ITS$drought_tolerance

groups <- dispersion_ITS$group
disper_ITS_rhizo <- ggplot(data=data, aes(x= irrigation, y=Distance_to_centroid, shape = irrigation, fill = drought_tolerance)) + 
  geom_jitter(size =3.5, position = position_jitter(0.1)) +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) + 
  scale_y_continuous(name = "Distance to Centroid") +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted") +
  ggtitle("Rhizosphere") +
  panel_border() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        plot.title = element_text(size=12,  face="bold"),
        legend.title = element_text(size = 12))
disper_ITS_rhizo

plot.data <- merge(summary(rda_ITS)$sites, df_ITS, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, shape = irrigation, fill = drought_tolerance), 
             size = 4, alpha = 1) +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression(italic("Full")), 
                                expression(italic("Reduced")))) +
  ggtitle("Rhizosphere") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_ITS)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_ITS)$cont)[2,2], 1)), "%)")) +
  theme_bw() +
  #geom_text(data = as.data.frame(summary(rda_ITS)$biplot) %>% rownames_to_column(), aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        plot.title = element_text(size=12,  face="bold"),
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) -> p_jaccard_rhizo_ITS

p_jaccard_rhizo_ITS

ggsave(filename="slide_26.tiff", plot=p_jaccard_rhizo_ITS, height=5.5, width=6.5, dpi=300)

#Bulk Soil
#dbRDA Plot for the Soil

df_ITS <- data.frame(rare_ITS_soil@sam_data)
df_ITS$irrigation.drought_tolerance <- with(rare_ITS_soil@sam_data, interaction(irrigation, drought_tolerance))
OTU_ITS <- t(as(rare_ITS_soil@otu_table, "matrix"))
OTU_ITS <- as.data.frame(OTU_ITS)
dist_ITS <- vegdist(OTU_ITS, "jaccard")
set.seed(1221990)
rda_ITS <- dbrda(dist_ITS ~ irrigation + drought_tolerance + irrigation*drought_tolerance, df_ITS)

anova(rda_ITS)
anova(rda_ITS, by="terms", permutation = 9999) 
plot(rda_ITS)

y <- varpart(dist_ITS, ~irrigation, ~drought_tolerance, data = df_ITS)
print(y)

# Measure for dispersion
dispersion_ITS <- betadisper(dist_ITS, df_ITS$irrigation.drought_tolerance, type = c("centroid"), bias.adjust = TRUE, sqrt.dist = FALSE, add = FALSE)
anova(dispersion_ITS)

set.seed(1221990)
permutest(dispersion_ITS, pairwise = TRUE, permutations = 9999)
(dispersion_ITS_HSD <- TukeyHSD(dispersion_ITS))

# Plot Dispersion

data <- data.frame(Distance_to_centroid=dispersion_ITS$distances, Group=dispersion_ITS$group)
data$irrigation <- df_ITS$irrigation
data$drought_tolerance <- df_ITS$drought_tolerance

groups <- dispersion_ITS$group
disper_ITS_soil <- ggplot(data=data, aes(x= irrigation, y=Distance_to_centroid, shape = irrigation, fill = drought_tolerance)) + 
  geom_jitter(size =3.5, position = position_jitter(0.1)) +
  scale_x_discrete(name = NULL, labels = c("Full", "Reduced")) + 
  scale_y_continuous(name = "Distance to Centroid") +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression("Full"), 
                                expression("Reduced"))) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted") +
  ggtitle("Bulk Soil") +
  panel_border() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        plot.title = element_text(size=12,  face="bold"),
        legend.title = element_text(size = 12))
disper_ITS_soil




plot.data <- merge(summary(rda_ITS)$sites, df_ITS, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, shape = irrigation, fill = drought_tolerance), 
             size = 4, alpha = 1) +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  scale_shape_manual(values = c(21,24), name = "Irrigation",
                     labels = c(expression(italic("Full")), 
                                expression(italic("Reduced")))) +
  ggtitle("Bulk Soil") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_ITS)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_ITS)$cont)[2,2], 1)), "%)")) +
  theme_bw() +
  #geom_text(data = as.data.frame(summary(rda_ITS)$biplot) %>% rownames_to_column(), aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        plot.title = element_text(size=12,  face="bold"),
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(pch = 21), title = "Drought Tolerance")) -> p_jaccard_soil_ITS

p_jaccard_soil_ITS

ggsave(filename="ITS_dbrda_jaccard_soil.tiff", plot=p_jaccard_soil, height=5.5, width=6.5, dpi=300)

# Append All plots Together 
prow <- plot_grid(p_jaccard_re_ITS + theme(legend.position = "none"), disper_ITS_re + theme(legend.position = "none"),
                  p_jaccard_rhizo_ITS + theme(legend.position = "none"), disper_ITS_rhizo + theme(legend.position = "none"),
                  p_jaccard_soil_ITS + theme(legend.position="none"), disper_ITS_soil + theme(legend.position = "none"),
                  ncol = 2, labels = "AUTO", align = "vh")
leg <- get_legend(p_jaccard_re_ITS)
prow <- plot_grid(prow, leg, ncol = 2, rel_widths = c(0.8, 0.23))
ggsave(filename=paste("ITS_treatment_beta_diversity.tiff", sep=""), plot=prow, width=7.5, height=7.5, dpi=600)


# Treatment Taxonomy------------------------------------------------------------

# ITS --------------------------------------------------------------------------

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

# 16S --------------------------------------------------------------------------

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

data_ITS_re_guild@sam_data$irrigation.drought_tolerance <- with(data_ITS_re_guild@sam_data, 
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
myco_melt_rhizo <- psmelt(Myco)

myco_melt_rhizo %>%
  group_by(irrigation, drought_tolerance, Guild) %>%
  summarise(sum = sum(Abundance)) %>%
  mutate(Habitat = "Rhizosphere") -> Myco_Rhizosphere

myco_melt_rhizo %>%
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
  ggtitle("Rhizosphere Mycorrhizal Abundance") +
  theme_bw() +
  theme(legend.text.align = 0) +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  guides(fill = guide_legend(title = "Drought Tolerance")) -> p_rhizo 

p_rhizo

ggsave("slide_28.tiff", plot = p_rhizo, width=7.5, height=6, dpi=300)

# Statistics 
model <- lme(Abundance ~ irrigation*drought_tolerance, random=~1|Row.names, data = myco_melt_rhizo, method="REML")
anova.lme(model, type="sequential", adjustSigma = T) 

x <- aov(sum ~ irrigation*drought_tolerance, data = myco_melt_rhizo %>%
           group_by(irrigation, drought_tolerance, Guild, Sample) %>%
           summarise(sum = sum(Abundance)*100) %>%
           filter(Guild == "Ectomycorrhizal"))
summary(x) 
TukeyHSD(x)

x <- aov(sum ~ irrigation*drought_tolerance, data = myco_melt_rhizo %>%
           group_by(irrigation, drought_tolerance, Guild, Sample) %>%
           summarise(sum = sum(Abundance)*100) %>%
           filter(Guild == "Arbuscular Mycorrhizal"))
summary(x) 
TukeyHSD(x)

# Root Endosphere

data_rel_ITS_re_tax <- as.data.frame(data_ITS_re_guild_norm@tax_table) # reformat names in tax table to include all labeled mycorrhizal otus
data_rel_ITS_re_tax$Guild[str_detect(data_rel_ITS_re_tax$Guild, "Ectomycorrhizal") == TRUE] <- "Ectomycorrhizal"
data_rel_ITS_re_tax$Guild[str_detect(data_rel_ITS_re_tax$Guild, "Arbuscular Mycorrhizal") == TRUE] <- "Arbuscular Mycorrhizal"
data_rel_ITS_re_tax <- as(data_rel_ITS_re_tax, "matrix")
data_rel_ITS_re_tax <- tax_table(data_rel_ITS_re_tax)
data_rel_ITS_re <- merge_phyloseq(otu_table(re_merged_df4), data_rel_ITS_re_tax, rare_ITS_RE@sam_data)
data_rel_ITS_re <- transform_sample_counts(data_ITS_re_guild_norm, function(x) x/sum(x))
Myco <- subset_taxa(data_rel_ITS_re, Guild == "Ectomycorrhizal" | Guild == "Arbuscular Mycorrhizal")
myco_melt_re <- psmelt(Myco)

myco_melt_re %>%
  group_by(irrigation, drought_tolerance, Guild, Sample) %>%
  summarise(sum = sum(Abundance)*100) %>%
  ungroup(irrigation, drought_tolerance, Guild, Sample) %>%
  group_by(irrigation, drought_tolerance, Guild) %>%
  mutate(Habitat = "Root Endosphere") -> Myco_re

myco_melt_re %>%
  group_by(irrigation, drought_tolerance, Guild, Sample) %>%
  summarise(sum = sum(Abundance)*100) %>%
  ungroup(irrigation, drought_tolerance, Guild, Sample) %>%
  group_by(irrigation, drought_tolerance, Guild) %>%
  summarise(mean = mean(sum), se = sd(sum)/sqrt(n())) %>%
  ggplot(aes(x = irrigation, y = mean, fill = drought_tolerance)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), alpha = 1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~Guild) +
  xlab("Irrigation") +
  scale_x_discrete(labels=c("Full", "Reduced")) +
  scale_y_continuous(name = "Relative Abundance (%)") +
  theme(axis.title.x = element_blank(),
        legend.position = "right", 
        axis.text.x = element_text(face = "italic"),
        #legend.justification = "center",
        legend.title = element_text(size = 12)) +
  ggtitle("Root Endosphere Mycorrhizal Abundance") +
  scale_fill_manual(values = c("tomato1", "steelblue4"), name = "Drought Tolerance",
                    labels = c(expression("Drought Tolerant"),
                               expression("Drought Susceptible"))) +
  theme_bw() +
  theme(legend.text.align = 0) +
  guides(fill = guide_legend(title = "Drought Tolerance")) -> p_re 
p_re

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


# Append all plots together
# Compositional Plot 
prow <- plot_grid(p_re + theme(legend.position = "none"), p_rhizo + theme(legend.position = "none"), align = "h", axis = "lr", labels = "AUTO", ncol = 2)
leg <- get_legend(p_re)
prow <- plot_grid(prow, leg, ncol = 2, rel_widths = c(1, 0.25))
prow
ggsave(filename=paste("myco_summary.tiff",sep=""), plot=prow, width=10, height=5, dpi=600)


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


# Differential Abundance Analysis ----------------------------------------------
# ANCOM Differential Abundance Analysis 
# 16S Datasets 

# Root Endosphere -- Differentially Abundant ASVs due to Irrigation with Drought Tolerant Covariate --------------------------

# Pre-Processing 
# Step 1: Data preprocessing

meta_data = rare_16S_RE@sam_data
meta_data <- as.data.frame(meta_data)
meta_data$names <- rownames(meta_data)

feature_table = rare_16S_RE@otu_table
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
re_tax <- as.data.frame(rare_16S_RE_ancom@tax_table)

res_re_irrigation$out <- merge(res_re_irrigation$out, re_tax %>% rownames_to_column(var = "rowname"),
                 by.x = 'taxa_id', by.y = 'rowname')

write.csv(res_re_irrigation$out, "ancom_re_irrigation.csv", row.names = F)
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
re_tax <- as.data.frame(rare_16S_RE_ancom@tax_table)

res_re_drought$out <- merge(res_re_drought$out, re_tax %>% rownames_to_column(var = "rowname"),
                 by.x = 'taxa_id', by.y = 'rowname')

write.csv(res$out, "ancom_re_drought_tolerance.csv", row.names = F)

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
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)

fig_ancom_re_drought_tolerance

leg_16S <- get_legend(fig_ancom_re_irrigation + theme(legend.position = "bottom"))
prow <- plot_grid(fig_ancom_re_irrigation + theme(legend.position = "none"), fig_ancom_re_drought_tolerance + theme(legend.position = "none"), align = "h", ncol = 2)
prow <- plot_grid(prow, leg_16S, ncol = 1, rel_heights = c(1, 0.2))
prow
ggsave(filename=paste("16S_ancom_re.tiff",sep=""), plot=prow, width=10, height=5, dpi=300)

# Rhizosphere -- Differentially Abundant ASVs due to Irrigation with Drought Tolerant Covariate ----------

# Pre-Processing 
# Step 1: Data preprocessing

meta_data = rare_16S_rhizosphere@sam_data
meta_data <- as.data.frame(meta_data)
meta_data$names <- rownames(meta_data)
feature_table = rare_16S_rhizosphere@otu_table
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
re_tax <- as.data.frame(rare_16S_rhizosphere@tax_table)

res_rhizo_irrigation$out <- merge(res_rhizo_irrigation$out, re_tax %>% rownames_to_column(var = "rowname"),
                 by.x = 'taxa_id', by.y = 'rowname')

write.csv(res_rhizo_irrigation$out, "ancom_rhizo_irrigation.csv", row.names = F)
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
re_tax <- as.data.frame(rare_16S_rhizosphere@tax_table)

res_rhizo_drought$out <- merge(res_rhizo_drought$out, re_tax %>% rownames_to_column(var = "rowname"),
                 by.x = 'taxa_id', by.y = 'rowname')

write.csv(res$out, "ancom_rhizo_drought_tolerance.csv", row.names = F)

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
prow <- plot_grid(prow, leg_16S, ncol = 1, rel_heights = c(1, 0.35))
prow
ggsave(filename=paste("16S_ancom_rhizo.tiff",sep=""), plot=prow, width=10, height=5, dpi=300)

ggsave(filename=paste("16S_ancom_rhizo.tiff",sep=""), plot=prow, width=4, height=5, dpi=300)

# Soil -- Differentially Abundant ASVs due to Irrigation with Drought Tolerant Covariate 

# Pre-Processing 
# Step 1: Data preprocessing

meta_data = rare_16S_soil@sam_data
meta_data <- as.data.frame(meta_data)
meta_data$names <- rownames(meta_data)
feature_table = rare_16S_soil@otu_table
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
re_tax <- as.data.frame(rare_16S_soil@tax_table)

res_soil_irrigation$out <- merge(res_soil_irrigation$out, re_tax %>% rownames_to_column(var = "rowname"),
                 by.x = 'taxa_id', by.y = 'rowname')

write.csv(res_soil_irrigation$out, "ancom_soil_irrigation.csv", row.names = F)
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
re_tax <- as.data.frame(rare_16S_soil@tax_table)

res_soil_drought_tolerance$out <- merge(res_soil_drought_tolerance$out, re_tax %>% rownames_to_column(var = "rowname"),
                 by.x = 'taxa_id', by.y = 'rowname')

write.csv(res_soil_drought_tolerance$out, "ancom_soil_drought_tolerance.csv", row.names = F)

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
prow <- plot_grid(prow, leg_16S, ncol = 1, rel_heights = c(1, 0.4))


ggsave(filename=paste("16S_ancom_soil.tiff",sep=""), plot=prow, width=10, height=7, dpi=600)

# ITS Datasets -----------------------------------------------------------

# Root Endosphere -- Differentially Abundant ASVs due to Irrigation with Drought Tolerant Covariate 


# Pre-Processing 
# Step 1: Data preprocessing

meta_data = data_ITS_RE@sam_data
meta_data <- as.data.frame(meta_data)
meta_data$names <- rownames(meta_data)
feature_table = rare_ITS_RE@otu_table
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

write.csv(res_its_re_irrigation$out, "ancom_re_irrigation_ITS.csv", row.names = F)
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

write.csv(res_its_re_drought_tolerance$out, "ancom_re_drought_tolerance_ITS.csv", row.names = F)

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

ggsave(filename=paste("ITS_ancom_re.tiff",sep=""), plot=prow, width=7.5, height=5, dpi=600)

# Rhizosphere -- Differentially Abundant ASVs due to Irrigation with Drought Tolerant Covariate 

# Pre-Processing 
# Step 1: Data preprocessing

meta_data = rare_ITS_rhizosphere@sam_data
meta_data <- as.data.frame(meta_data)
meta_data$names <- rownames(meta_data)
feature_table = rare_ITS_rhizosphere@otu_table
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
re_tax <- as.data.frame(rare_ITS_rhizosphere@tax_table)

res_its_rhizo_irrigation$out <- merge(res_its_rhizo_irrigation$out, re_tax %>% rownames_to_column(var = "rowname"),
                 by.x = 'taxa_id', by.y = 'rowname')

write.csv(res$out, "ancom_rhizo_irrigation_ITS.csv", row.names = F)
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
re_tax <- as.data.frame(rare_ITS_rhizosphere@tax_table)

res_its_rhizo_drought_tolerance$out <- merge(res_its_rhizo_drought_tolerance$out, re_tax %>% rownames_to_column(var = "rowname"),
                 by.x = 'taxa_id', by.y = 'rowname')

write.csv(res_its_rhizo_drought_tolerance$out, "ancom_rhizo_drought_tolerance_ITS.csv", row.names = F)

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

ggsave(filename=paste("ITS_ancom_rhizo.tiff",sep=""), plot=prow, width=7.5, height=5, dpi=600)


# Soil -- Differentially Abundant ASVs due to Irrigation with Drought Tolerant Covariate 

# Pre-Processing 
# Step 1: Data preprocessing

meta_data = rare_ITS_soil@sam_data
meta_data <- as.data.frame(meta_data)
meta_data$names <- rownames(meta_data)
feature_table = rare_ITS_soil@otu_table
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
re_tax <- as.data.frame(rare_ITS_soil@tax_table)

res_its_soil_irrigation$out <- merge(res_its_soil_irrigation$out, re_tax %>% rownames_to_column(var = "rowname"),
                 by.x = 'taxa_id', by.y = 'rowname')

write.csv(res_its_soil_irrigation$out, "ancom_soil_irrigation_ITS.csv", row.names = F)
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
re_tax <- as.data.frame(rare_ITS_soil@tax_table)

res_its_soil_drought_tolerance$out <- merge(res_its_soil_drought_tolerance$out, re_tax %>% rownames_to_column(var = "rowname"),
                 by.x = 'taxa_id', by.y = 'rowname')

write.csv(res_its_soil_drought_tolerance$out, "ancom_soil_drought_tolerance.csv", row.names = F)

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

ggsave(filename=paste("ITS_ancom_soil.tiff",sep=""), plot=prow, width=7.5, height=5, dpi=600)


# Differential Abundance Heatmaps ----------------------------------------------
# Drought Tolerant Taxa in the Root Endosphere
# Subset OTU table to include only significant ASVS
drought_re_asv <- subset(otu_table(data_16S_RE), rownames(otu_table(data_16S_RE)) %in% c('2dd9e191a81b91117b77df5eac7e6048',
                                                                                         '72d325488d710b04688e0d14d0fea79d',
                                                                                         'c757ad670c5d796cb866c941c62c4e33',
                                                                                         'aee9f354c80ca7baa872c3da2fe462c2',
                                                                                         '41300091835845bd9f0d2c9d5efa62e5'))
drought_re_asv <- merge_phyloseq(drought_re_asv, tax_table(data_16S_RE), sample_data(data_16S_RE))
# Plot Heatmap
heatmap1 <- plot_heatmap(drought_re_asv, method = NULL, sample.label = "drought_tolerance", 
                         sample.order = "drought_tolerance", taxa.label = "Family",
                         title = "16S Root Endosphere - Significant Drought Tolerant ASVs")
heatmap1
ggsave("ancom_re_drought_tolerance_heatmap.tiff", plot=heatmap1, width=10, height=5, dpi=300)
# Irrigation Taxa in the Root Endosphere
# Subset OTU table to include only significant ASVS
irrigation_re_asv <- subset(otu_table(data_16S_RE), rownames(otu_table(data_16S_RE)) %in% c('b3c195907b2347a887d116e237540feb',
                                                                                            'bc3676c17839094c4fa8e33268905268',
                                                                                            'e984a619f22fe030e13c709b0489b52a',
                                                                                            '754c7c6f3b6988d240275e2f44436292',
                                                                                            'b3bc75b5a73133d32c26ffd8ae90d41d',
                                                                                            '2653f9a6e90d957c951070035e6edeea',
                                                                                            '47a9d10e32ffbceb163e240b65825d7e',
                                                                                            '2740cf2417c92847cc298cbd71dd1fcd',
                                                                                            '1c2b8892a0d7f828d72808cd01cfd199',
                                                                                            'bcf0b339e6aacfad9c466f86e0c5971d',
                                                                                            '2f6a15af8c9662afab9a6c1f0593eb88',
                                                                                            'df71de0e92459aa808364311135bb8fc'))
irrigation_re_asv <- merge_phyloseq(irrigation_re_asv, tax_table(data_16S_RE), sample_data(data_16S_RE))
# Plot Heatmap
heatmap2 <- plot_heatmap(irrigation_re_asv, method = NULL, sample.label = "irrigation", 
                         sample.order = "irrigation", taxa.label = "Family",
                         title = "16S Root Endosphere - Significant Irrigation ASVs")
heatmap2
ggsave("ancom_re_irrigation_heatmap.tiff", plot=heatmap2, width=10, height=5, dpi=300)

# Drought Tolerant Taxa in the Rhizosphere
# Subset OTU table to include only significant ASVS

drought_rhizo_asv <- subset(otu_table(data_16S_rhizosphere), rownames(otu_table(data_16S_rhizosphere)) %in% c('c993a60d06e12487198bc56fa24c625f',
                                                                                                              '9ac3bb1d7dedb08a012692a6f536b5af'))
drought_rhizo_asv <- merge_phyloseq(drought_rhizo_asv, tax_table(data_16S_rhizosphere), sample_data(data_16S_rhizosphere))

# Plot Heatmap
heatmap3 <- plot_heatmap(drought_rhizo_asv, method = NULL, sample.label = "drought_tolerance", 
                         sample.order = "drought_tolerance", taxa.label = "Family",
                         title = "16S Rhizosphere - Significant Drought Tolerant ASVs")
heatmap3
ggsave("ancom_rhizo_drought_tolerance_heatmap.tiff", plot=heatmap3, width=10, height=5, dpi=300)

# Irrigation Taxa in the Rhizosphere
# Subset OTU table to include only significant ASVS
irrigation_rhizo_asv <- subset(otu_table(data_16S_rhizosphere), rownames(otu_table(data_16S_rhizosphere)) %in% c('0ddcd311e02f742e2e0e61ce02cf9c29',
                                                                                                                 '572b96fefec73aa11fc6628597853ed7',
                                                                                                                 'aee9f354c80ca7baa872c3da2fe462c2',
                                                                                                                 '0ec20bde5b936bac2e59d048604f9358',
                                                                                                                 'e1c7d97fe13e9127d225d76a1feb8c78',
                                                                                                                 '95d08381352f71ef3bc48445f8ccab02',
                                                                                                                 '3125694452d98013214919bea3df321d',
                                                                                                                 'c7c40b0f480f3edee216585ec048011d',
                                                                                                                 '9ebd33169a6be011a108ff2661ae2bc6',
                                                                                                                 '0dcffe423f594ed936b5435a15126443',
                                                                                                                 '391c94b85cd374ec9623a32e9f9b4a03',
                                                                                                                 '945184b6386c192c0066e0a98a154780',
                                                                                                                 'a59a9a9b27bc75e0323284f8ac75d55e',
                                                                                                                 'f9addd8a736ba8f99f715836c4086005',
                                                                                                                 'c0d91dde8e97d4e7e3aa705c6ccb18e4',
                                                                                                                 'de217a040cee6a1b45efadffa736efbf',
                                                                                                                 '9d7b1963f8036854dd7b9aad0149310f',
                                                                                                                 '8597d6b9d76ccfe81ca18d0bda2c6214',
                                                                                                                 'a1ed615a4ed0b08c534fcc504a4608b2',
                                                                                                                 'dda3edc4376cd111ef09921cf1a27da8',
                                                                                                                 'daed5a5e05797a803338c4210f1d88c7',
                                                                                                                 '2ee683d3ac3574f69eee1a2a1cd89d09',
                                                                                                                 'ff836ceb31404b3a37e52a4af998742b',
                                                                                                                 '0eb88b722902316b11770922bcdaca7b',
                                                                                                                 '70fb4117bc0802189a5c54f648668ac9',
                                                                                                                 'b4b365be5234d748dd5812df4ac2c6de',
                                                                                                                 'c9d59b5dca26cd50adf46d0171c1e087',
                                                                                                                 '40fa96a2f7fc7b939f5eb0d52d7e5b3c',
                                                                                                                 'e05c79009fd7877b24e3127bd49986a0',
                                                                                                                 '91a4ee945ca39dbc22df710686730f5b',
                                                                                                                 'fb8ce5c795542b9f42e85c54f344f061',
                                                                                                                 'fb60fb75664eaf6a23ae2c6124fbb005',
                                                                                                                 'b968fc1f8ed4d7fbacac9d1f3a23ffc4',
                                                                                                                 '1cf8836d9f2282cac057b416706ba561',
                                                                                                                 'cfdf66ad8b27626f14f419a4f6c51b88',
                                                                                                                 'e4e08e332b44c5b0f3c79143fb9c0f63',
                                                                                                                 'f6c528c217081be9871de0fc58b89408',
                                                                                                                 'df71de0e92459aa808364311135bb8fc',
                                                                                                                 '9e56d913c8bbdccb6b373fb563f06f2a',
                                                                                                                 '8dad047cb6374c9c0f872f5c36dcfea1',
                                                                                                                 '7df299814b7a37f021b9bc3c3fc40080',
                                                                                                                 'e61a7bc42a30f0ca4f170b5e3f6cad0c',
                                                                                                                 'a70b83cbaaf606acfdb78362ffeaa506',
                                                                                                                 'fb5c2aa7e7a97d1ba29b2c9405cb6a41',
                                                                                                                 'f81076e084174a484f7d0c302f94437d',
                                                                                                                 '3e67e934cc725d230d800c4c980327c8',
                                                                                                                 'ddfaffe2f2f94985e611459236b87647',
                                                                                                                 '25cb02fb17426d878afd69a79572b858',
                                                                                                                 '833d619b70b86d87f185a0c98ac2fe57',
                                                                                                                 '4bc2853b48494b7dd698fccb74a67f56',
                                                                                                                 '3ebe761bfb1238c87195d431f41bf976',
                                                                                                                 '50fb93dad525b440dc193bc631872eba',
                                                                                                                 'c4def24d4e9cea12f3f0a8a0ca43afda',
                                                                                                                 'a0b1f105f827e856c5fdb5f98eaecc09',
                                                                                                                 '2e6bc9a35299587b295c308e5f616477',
                                                                                                                 '41235f78280d0304e45c61651f6237db',
                                                                                                                 '43e1ec8977fc7de1d298cbe9e7ce4fae',
                                                                                                                 '034e2e9539ccb9186ee9b67b0429f29f',
                                                                                                                 '45197ad8ceb6fc7a202b0162f679af5b',
                                                                                                                 'fed4e369086175f7f4bafdf67f1c4dcf',
                                                                                                                 '7684072284024ca2b426a6bdd04c45f9',
                                                                                                                 '163076c2d7e10152bb3d0e2f4b28df75',
                                                                                                                 'd1b7aa5cdc783723fadc6cdc643472f0',
                                                                                                                 '0e11a8e894a8d9a3241d7cd7915570ac',
                                                                                                                 'fe170a37c40648f3cd117b40946280fc',
                                                                                                                 '333048cac277d86e040b020b11c3eabf',
                                                                                                                 'c252134efa44f25e744c20e9ade0781e',
                                                                                                                 'e07052bb9c0f0d2370537c4b9acbed4e',
                                                                                                                 '9fe3cf958c80040f5c4ff05b64278e8a',
                                                                                                                 '96cc99ed28e8cab1001fae99d2ff8c2e',
                                                                                                                 'b48a188a24ed126894db7d73bd6617cc',
                                                                                                                 '17f654dabb5544c03d515eea432c1632',
                                                                                                                 'd554f717d25d5f90590ec252d3c67d0d',
                                                                                                                 '13645da5ce1259a647f3bc11b09dcf13',
                                                                                                                 'ec414a9986a44003df0ebd401213b7d0',
                                                                                                                 '7b7f9c073366c233fa85c7a9aee5903c',
                                                                                                                 'eb6091d870581bf798b9bf7056cee96e',
                                                                                                                 '2bb412e90f3ddf0697b2d449b44a05b8',
                                                                                                                 '594061ed4e2ea378773c474cd4803c29',
                                                                                                                 'caa1df79927af52cdcf97e08f266aedd',
                                                                                                                 '3a54d4e72af90b9c2fd50d7d1bb0336b',
                                                                                                                 '690e2d312e3f83ab4ef9893d8c5eb3ab',
                                                                                                                 '1394f8ffb5c637943daae55332628681',
                                                                                                                 '7d9564562bb10cea1553ad90158d8545',
                                                                                                                 '875dc0da5f66b1b5644788543a07544f',
                                                                                                                 '99861da6145a8efcd3524241496cabe9',
                                                                                                                 'd4cb29080e56a06ef36f99b209f2b225',
                                                                                                                 '29d6eab2f173318980e3f98caf5af437',
                                                                                                                 'd2d78251d7df8d1a835217535550cb59',
                                                                                                                 '3c008d3fda79db47dd7167837a297c48',
                                                                                                                 'd19dbb0e1566e9a2d91adda6f64e0232',
                                                                                                                 '0f8ac0d81dd6449a43b34fcb184efd8f',
                                                                                                                 '7034493f95fe775db5f9eeb1aac0cad9',
                                                                                                                 'd58b6cf2e9af6364330b92566e14d59e',
                                                                                                                 '52fd95eb476921cf2b9510888ea0cfe8',
                                                                                                                 'eea13d94a260c6a4c38e783eb79c68d0',
                                                                                                                 '809b77a90319f0cda021696d7150f1f0',
                                                                                                                 '334da1236bd220d991a0acc859908b51',
                                                                                                                 'c3b8b45ccbc3b2aac60dda6ef98d2863',
                                                                                                                 '002664222d5c2d5b79c41b33598ea1e4',
                                                                                                                 '802e0feee197426ddbe8feb3f3260cd3',
                                                                                                                 'f47beff0fc3748fcf424bf33bc490cf2',
                                                                                                                 '63321f408932c626da9ed1f1e73561e9',
                                                                                                                 '8ff49840b8cdb3c84d65399a75f2f4c3',
                                                                                                                 'e3bbf1bd6c612fbffbbf88f4dc54b805',
                                                                                                                 'b9875bae83467e83cbb87f0bb6ed2db0',
                                                                                                                 '2653f9a6e90d957c951070035e6edeea',
                                                                                                                 '17e562606dc747dd9f950ce27d3eb95a',
                                                                                                                 'd51c399452cc3e608552d0aca841b28c',
                                                                                                                 'a7b330d77a6e4aab2ff36ba6314ac038',
                                                                                                                 'a34c67c3e970be24a82f8b4c537c08ed',
                                                                                                                 'bd88cc06608d082fa8df87333f49b4e0',
                                                                                                                 'f91201adfe8be6ff17177daa8bb13518',
                                                                                                                 '816ebb72c57410cfcc6c9c5db3a37129',
                                                                                                                 'e8bd5a5228cd426c496dcd50f6dad542',
                                                                                                                 '3bed507c2cb7458e3bd9b0ee3692fecb',
                                                                                                                 'bdd15f551778f16db860acaa49aa6713',
                                                                                                                 '4922eb95c3994f40d0784f2cfbf7ede5',
                                                                                                                 '3706633740d9189bb93aa193fb22ba60',
                                                                                                                 'c4b069153ae53489ea1cbfb2de1cfb1d',
                                                                                                                 '96799da3b007f0ef0e63b89430bfef93',
                                                                                                                 'd0c8bf44946dc8315c6b9629eb472477',
                                                                                                                 '614ed31805095132a7610c0e4ecbc205',
                                                                                                                 '2fac845da336657d55f1a9f350c6bc09',
                                                                                                                 '0a1967051180fbda882271ec5a9df0d9',
                                                                                                                 'd892d40bc79d0c356bd2886ce8d120ee',
                                                                                                                 'fa96cbd7b2df953f48e159b78c51527d',
                                                                                                                 'f66df1898b46126275f1585978e60df0',
                                                                                                                 '37db27274d1534b5ff701cfe737c1803',
                                                                                                                 '1bc625bd6d3510d38f0b5b3f460bdc50',
                                                                                                                 '91878d7c16811333ee5db6102c3be33e',
                                                                                                                 '2f35bdf3cbfaa4cd319a51cc46764656',
                                                                                                                 '6376ea6dab7d0fde3cd66f53b57e1484',
                                                                                                                 'aa11b7fecbbca304e3b36cce12f5bb66',
                                                                                                                 'ac68c7c702432c26bc13677cb4c1ea8f',
                                                                                                                 '0b5f0a192b92c230e826d1ceb10f9d39',
                                                                                                                 '7cf569dfce221cef78261d24e0ee9aba',
                                                                                                                 '1e7808d3462d49317c9bcbe880a8cbb3'))

irrigation_rhizo_asv <- merge_phyloseq(irrigation_rhizo_asv, tax_table(data_16S_rhizosphere), sample_data(data_16S_rhizosphere))
# Plot Heatmap
heatmap4 <- plot_heatmap(irrigation_rhizo_asv, method = NULL, sample.label = "irrigation", 
                         sample.order = "irrigation", taxa.label = "Phylum", taxa.order = "Phylum", 
                         title = "16S Rhizosphere - Significant Irrigation ASVs")
heatmap4
ggsave("ancom_rhizo_irrigation_heatmap.tiff", plot=heatmap4, width=10, height=7.5, dpi=600)

# ITS Heatmaps
# Rhizosphere 
# Irrigation 
irrigation_rhizo_asv <- subset(otu_table(data_ITS_rhizosphere), rownames(otu_table(data_ITS_rhizosphere)) %in% c('ddfb725268acf837fad0ee338ff65c48',
                                                                                                                 '0ae3b4f6852b6ff563f44369fb98ad5b',
                                                                                                                 '6ef2e36136cf7abb4b761ea5d2a09cde',
                                                                                                                 'e93a30988c5ec5ebfb7fb7f63173a46b',
                                                                                                                 '2806475cca11a74a5c4a7d25237ddb82',
                                                                                                                 '03fde61dc4c88d74269ba23df5bd6228',
                                                                                                                 'd5f74c08d012c17c4eea05a0bf0e82f8',
                                                                                                                 'e9211240d9c06829f6e32a3d70c158a0',
                                                                                                                 'f1b6d8b651c97d9c236773b5455c75f2',
                                                                                                                 '90b66fdf3dbc01ffe01ff138d0d137b8',
                                                                                                                 '270ab69cdcce97bb7a5c464fbde7bfaf',
                                                                                                                 '8f797c8b6e732c0f7bf9f39dffa0b7bd',
                                                                                                                 '146cc8e405249ba0e1cf5ab510e13934',
                                                                                                                 '40c07191656c78fb28cc3167802945a6',
                                                                                                                 '7254eb7970555c4265513095417bd51c',
                                                                                                                 'f4d6675d1feb2de01ddccb53ecdcb556',
                                                                                                                 '758895697a468b41e7a731e4e8c10915',
                                                                                                                 '4411be6aa5f6f72703173d7def1eb076',
                                                                                                                 '91f57df75fcdea17350c8d42d7ed4b4a',
                                                                                                                 '378b3b7b4bc00ee26fdd851f0ee00849',
                                                                                                                 '3b5b00c96015fc6801eb8a75ce491c15',
                                                                                                                 'dbad2bc744548baef8805553e9e7c829',
                                                                                                                 '52f19787a7ce17c334d1b4dac36ad475',
                                                                                                                 '0182952ed15b504a960fcf13ecd856d4',
                                                                                                                 'ae14245b252cfdfbf0cffd0fbf6b6abd',
                                                                                                                 '84040ce2d7a10c0da47188d46f7cc558'))

irrigation_rhizo_asv <- merge_phyloseq(irrigation_rhizo_asv, tax_table(data_ITS_rhizosphere), sample_data(data_ITS_rhizosphere))
# Plot Heatmap
heatmap4 <- plot_heatmap(irrigation_rhizo_asv, method = NULL, sample.label = "irrigation", 
                         sample.order = "irrigation", taxa.label = "Family", taxa.order = "Family", 
                         title = "ITS Rhizosphere - Significant Irrigation ASVs")
heatmap4
ggsave("ancom_rhizo_irrigation_heatmap_ITS.tiff", plot=heatmap4, width=10, height=7.5, dpi=600)

# Drought Tolerance 
drought_rhizo_asv <- subset(otu_table(data_ITS_rhizosphere), rownames(otu_table(data_ITS_rhizosphere)) %in% c('894a0203a42894a47bdde5f2c94e62e8',
                                                                                                              '74c4e4f306d9988a8d4e4e08a2ea0670'))
drought_rhizo_asv <- merge_phyloseq(drought_rhizo_asv, tax_table(data_ITS_rhizosphere), sample_data(data_ITS_rhizosphere))

# Plot Heatmap
heatmap3 <- plot_heatmap(drought_rhizo_asv, method = NULL, sample.label = "drought_tolerance", 
                         sample.order = "drought_tolerance", taxa.label = "Family",
                         title = "ITS Rhizosphere - Significant Drought Tolerant ASVs")
heatmap3
ggsave("ancom_rhizo_drought_tolerance_heatmap_ITS.tiff", plot=heatmap3, width=10, height=5, dpi=300)

# Network Analysis ---------------------------------------------------
# Root Endosphere Network --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
output.path <- "C:/Users/Brand/Dropbox/Cregger_CBI_drought/results/treatment_analysis/network/"

# Increase filtering parameters to improve network resolution
data_16S_root_endosphere <- subset_samples(data_16S_filtered, habitat=="root_endosphere")%>%
  filter_taxa(function(x) sum(x > 5) > (0.2*length(x)), TRUE)
data_16S_re_reduced <- subset_samples(data_16S_root_endosphere, irrigation=="reduced")
data_16S_re_reduced_l030 <- subset_samples(data_16S_re_reduced, drought_tolerance =="LO30")
data_16S_re_reduced_hi30 <- subset_samples(data_16S_re_reduced, drought_tolerance == "HI30")
data_16S_re_full <- subset_samples(data_16S_root_endosphere, irrigation =="full")
data_16S_re_full_l030 <- subset_samples(data_16S_re_full, drought_tolerance == "LO30")
data_16S_re_full_hi30 <- subset_samples(data_16S_re_full, drought_tolerance == "HI30")

data_ITS_root_endosphere <-subset_samples(data_ITS_filtered, habitat=="root_endosphere")%>%
  filter_taxa(function(x) sum(x > 5) > (0.2*length(x)), TRUE)
data_ITS_re_reduced <- subset_samples(data_ITS_root_endosphere, irrigation =="reduced")
data_ITS_re_reduced_l030 <- subset_samples(data_ITS_re_reduced, drought_tolerance =="LO30")
data_ITS_re_reduced_hi30 <- subset_samples(data_ITS_re_reduced, drought_tolerance == "HI30")
data_ITS_re_full <- subset_samples(data_ITS_root_endosphere, irrigation =="full")
data_ITS_re_full_l030 <- subset_samples(data_ITS_re_full, drought_tolerance == "LO30")
data_ITS_re_full_hi30 <- subset_samples(data_ITS_re_full, drought_tolerance == "HI30")

# Generate subnetworks across each treatment
data_ITS_re_reduced_l030 <- subset_samples(data_ITS_re_reduced_l030, data_ITS_re_reduced_l030@sam_data@row.names!="S-GW9583-3-5-37-rootendo")

se.cor_RE_reduced_l030 <- spiec.easi(list(data_ITS_re_reduced_l030, data_16S_re_reduced_l030), 
                                     method='mb', 
                                     nlambda=100, 
                                     lambda.min.ratio=1e-1, 
                                     pulsar.params = list(thresh = 0.01))
data_16S_re_reduced_hi30 <- subset_samples(data_16S_re_reduced_hi30, data_16S_re_reduced_hi30@sam_data@row.names!="S-56-4-25-9-rootendo")
se.cor_RE_reduced_hi30 <- spiec.easi(list(data_ITS_re_reduced_hi30, data_16S_re_reduced_hi30),
                                     method='mb',
                                     nlambda=100,
                                     lambda.min.ratio=1e-1,
                                     pulsar.params = list(thresh= 0.01))

se.cor_RE_full_l030 <- spiec.easi(list(data_ITS_re_full_l030, data_16S_re_full_l030),
                                  method='mb',
                                  nlambda=100,
                                  lambda.min.ratio=1e-1,
                                  pulsar.params = list(thresh= 0.01))

se.cor_RE_full_hi30 <- spiec.easi(list(data_ITS_re_full_hi30, data_16S_re_full_hi30),
                                  method='mb',
                                  nlambda=100,
                                  lambda.min.ratio=1e-1,
                                  pulsar.params= list(thresh= 0.01))

physeq_combo_reduced_l030 <-merge_phyloseq(data_ITS_re_reduced_l030, data_16S_re_reduced_l030)
ig_reduced_l030 <- adj2igraph(getRefit(se.cor_RE_reduced_l030),vertex.attr=list(name=c(taxa_names(data_ITS_re_reduced_l030),taxa_names(data_16S_re_reduced_l030))))
write.table(physeq_combo_reduced_l030@tax_table,file=file.path(output.path,"reduced_l030_tax.txt"),sep="\t", quote=FALSE)

physeq_combo_reduced_hi30 <-merge_phyloseq(data_ITS_re_reduced_hi30, data_16S_re_reduced_hi30)
ig_reduced_hi30 <- adj2igraph(getRefit(se.cor_RE_reduced_hi30),vertex.attr=list(name=c(taxa_names(data_ITS_re_reduced_hi30),taxa_names(data_16S_re_reduced_hi30))))
write.table(physeq_combo_reduced_hi30@tax_table,file=file.path(output.path,"reduced_hi30_tax.txt"),sep="\t", quote=FALSE)

physeq_combo_full_l030 <-merge_phyloseq(data_ITS_re_full_l030, data_16S_re_full_l030)
ig_full_l030 <- adj2igraph(getRefit(se.cor_RE_full_l030),vertex.attr=list(name=c(taxa_names(data_ITS_re_full_l030),taxa_names(data_16S_re_full_l030))))
write.table(physeq_combo_full_l030@tax_table,file=file.path(output.path,"full_l030_tax.txt"),sep="\t", quote=FALSE)

physeq_combo_full_hi30 <-merge_phyloseq(data_ITS_re_full_hi30, data_16S_re_full_hi30)
ig_full_hi30 <- adj2igraph(getRefit(se.cor_RE_full_hi30),vertex.attr=list(name=c(taxa_names(data_ITS_re_full_hi30),taxa_names(data_16S_re_full_hi30))))
write.table(physeq_combo_full_hi30@tax_table,file=file.path(output.path,"full_hi30_tax.txt"),sep="\t", quote=FALSE)

library(RColorBrewer)
library(pals)
#make vectors to plot all shapes and colors the same
Kingdom_reduced_l030 <-unique(tax_table(physeq_combo_reduced_l030)[,"Kingdom"])
Kingdom_reduced_hi30 <-unique(tax_table(physeq_combo_reduced_hi30)[,"Kingdom"])
Kingdom_full_l030 <-unique(tax_table(physeq_combo_full_l030)[,"Kingdom"])
Kingdom_full_hi30 <-unique(tax_table(physeq_combo_full_hi30)[,"Kingdom"])
KingdomList<-do.call("rbind", list(Kingdom_reduced_l030,Kingdom_reduced_hi30,Kingdom_full_l030,Kingdom_full_hi30))%>%unique()
KingdomShape = c(16,17,15)
names(KingdomShape) = KingdomList


getPalette = colorRampPalette(brewer.pal(9, "Set1"))
Phylum_reduced_l030 <-unique(tax_table(physeq_combo_reduced_l030)[,"Phylum"])
Phylum_reduced_hi30 <-unique(tax_table(physeq_combo_reduced_hi30)[,"Phylum"])
Phylum_full_l030 <-unique(tax_table(physeq_combo_full_l030)[,"Phylum"])
Phylum_full_hi30 <-unique(tax_table(physeq_combo_full_hi30)[,"Phylum"])
PhylumList<-do.call("rbind", list(Phylum_reduced_l030, Phylum_reduced_hi30, Phylum_full_l030, Phylum_full_hi30))%>%unique()

PhylumPalette = getPalette(length(PhylumList))
names(PhylumPalette) = PhylumList

#plot networks
reduced_l030<-plot_network(ig_reduced_l030,physeq_combo_reduced_l030,type="taxa",color="Phylum",shape="Kingdom",label=NULL,point_size = 1)+scale_shape_manual(values=KingdomShape,drop=FALSE)+scale_color_manual(values=PhylumPalette,drop=FALSE)
reduced_hi30<-plot_network(ig_reduced_hi30,physeq_combo_reduced_hi30,type="taxa",color="Phylum",shape="Kingdom",label=NULL,point_size=1)+scale_shape_manual(values=KingdomShape)+scale_color_manual(values=PhylumPalette)
full_l030<-plot_network(ig_full_l030,physeq_combo_full_l030,type="taxa",color="Phylum",shape="Kingdom",label=NULL,point_size=1)+scale_shape_manual(values=KingdomShape)+scale_color_manual(values=PhylumPalette)
full_hi30<-plot_network(ig_full_hi30,physeq_combo_full_hi30,type="taxa",color="Phylum",shape="Kingdom",label=NULL,point_size=1)+scale_shape_manual(values=KingdomShape)+scale_color_manual(values=PhylumPalette)


re_network<-ggarrange(reduced_l030,reduced_hi30,full_l030,full_hi30,common.legend = TRUE,legend="bottom", labels=c("A","B","C","D"))
ggsave("network_re.tiff",plot=re_network, width=11,height=7, device="tiff",dpi=600)

# Export Network Information
taxa.RE_reduced_l030 <-data.frame(tax_table(physeq_combo_reduced_l030))
betaMat_RE_reduced_l030=as.matrix(symBeta(getOptBeta(se.cor_RE_reduced_l030)))
total_re_reduced_l030 <- length(betaMat_RE_reduced_l030[betaMat_RE_reduced_l030!=0])/2 

taxa.RE_reduced_hi30 <-data.frame(tax_table(physeq_combo_reduced_hi30))
betaMat_RE_reduced_hi30=as.matrix(symBeta(getOptBeta(se.cor_RE_reduced_hi30)))
total_re_reduced_hi30 <- length(betaMat_RE_reduced_hi30[betaMat_RE_reduced_hi30!=0])/2 

taxa.RE_full_l030 <-data.frame(tax_table(physeq_combo_full_l030))
betaMat_RE_full_l030=as.matrix(symBeta(getOptBeta(se.cor_RE_full_l030)))
total_re_full_l030 <- length(betaMat_RE_full_l030[betaMat_RE_full_l030!=0])/2 

taxa.RE_full_hi30 <-data.frame(tax_table(physeq_combo_full_hi30))
betaMat_RE_full_hi30 =as.matrix(symBeta(getOptBeta(se.cor_RE_full_hi30)))

wtRE_reduced_l030 <- symBeta(getOptBeta(se.cor_RE_reduced_l030), mode="maxabs")
weightRE_reduced_l030 <-Matrix::summary(t(wtRE_reduced_l030))[,3]
igRE_reduced_l030_wt <-adj2igraph(getRefit(se.cor_RE_reduced_l030),edge.attr=list(weight=weightRE_reduced_l030),vertex.attr=list(name=c(taxa_names(data_16S_re_reduced_l030),taxa_names(data_ITS_re_reduced_l030))))
write.graph(igRE_reduced_l030_wt,file=file.path("spieceasi_reduced_l030.txt"),format="ncol") 

wtRE_reduced_hi30 <- symBeta(getOptBeta(se.cor_RE_reduced_hi30), mode="maxabs")
weightRE_reduced_hi30 <-Matrix::summary(t(wtRE_reduced_hi30))[,3]
igRE_reduced_hi30_wt <-adj2igraph(getRefit(se.cor_RE_reduced_hi30),edge.attr=list(weight=weightRE_reduced_hi30),vertex.attr=list(name=c(taxa_names(data_16S_re_reduced_hi30),taxa_names(data_ITS_re_reduced_hi30))))
write.graph(igRE_reduced_hi30_wt,file=file.path("spieceasi_reduced_hi30.txt"),format="ncol") 

wtRE_full_l030 <- symBeta(getOptBeta(se.cor_RE_full_l030), mode="maxabs")
weightRE_full_l030 <-Matrix::summary(t(wtRE_full_l030))[,3]
igRE_full_l030_wt <-adj2igraph(getRefit(se.cor_RE_full_l030),edge.attr=list(weight=weightRE_full_l030),vertex.attr=list(name=c(taxa_names(data_16S_re_full_l030),taxa_names(data_ITS_re_full_l030))))
write.graph(igRE_full_l030_wt,file=file.path("spieceasi_full_l030.txt"),format="ncol") 

wtRE_full_hi30 <- symBeta(getOptBeta(se.cor_RE_full_hi30), mode="maxabs")
weightRE_full_hi30 <-Matrix::summary(t(wtRE_full_hi30))[,3]
igRE_full_hi30_wt <-adj2igraph(getRefit(se.cor_RE_full_hi30),edge.attr=list(weight=weightRE_full_hi30),vertex.attr=list(name=c(taxa_names(data_16S_re_full_hi30),taxa_names(data_ITS_re_full_hi30))))
write.graph(igRE_full_hi30_wt,file=file.path("spieceasi_full_hi30.txt"),format="ncol") 

# Rhizosphere Network --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
data_16S_rhizosphere <- subset_samples(data_16S_filtered, habitat=="rhizosphere")%>%
  filter_taxa(function(x) sum(x > 5) > (0.4*length(x)), TRUE)
data_16S_rhizo_reduced <- subset_samples(data_16S_rhizosphere, irrigation=="reduced")
data_16S_rhizo_reduced_l030 <- subset_samples(data_16S_rhizo_reduced, drought_tolerance =="LO30")
data_16S_rhizo_reduced_hi30 <- subset_samples(data_16S_rhizo_reduced, drought_tolerance == "HI30")
data_16S_rhizo_full <- subset_samples(data_16S_rhizosphere, irrigation =="full")
data_16S_rhizo_full_l030 <- subset_samples(data_16S_rhizo_full, drought_tolerance == "LO30")
data_16S_rhizo_full_hi30 <- subset_samples(data_16S_rhizo_full, drought_tolerance == "HI30")

data_ITS_rhizosphere <-subset_samples(data_ITS_filtered, habitat=="rhizosphere")%>%
  filter_taxa(function(x) sum(x > 5) > (0.4*length(x)), TRUE)
data_ITS_rhizo_reduced <- subset_samples(data_ITS_rhizosphere, irrigation=="reduced")
data_ITS_rhizo_reduced_l030 <- subset_samples(data_ITS_rhizo_reduced, drought_tolerance =="LO30")
data_ITS_rhizo_reduced_hi30 <- subset_samples(data_ITS_rhizo_reduced, drought_tolerance == "HI30")
data_ITS_rhizo_full <- subset_samples(data_ITS_rhizosphere, irrigation =="full")
data_ITS_rhizo_full_l030 <- subset_samples(data_ITS_rhizo_full, drought_tolerance == "LO30")
data_ITS_rhizo_full_hi30 <- subset_samples(data_ITS_rhizo_full, drought_tolerance == "HI30")

#remove sample not present in the ITS for reduced x LO30: 1024-4-17-37
data_ITS_rhizo_reduced_l030 <- subset_samples(data_ITS_rhizo_reduced_l030, data_ITS_rhizo_reduced_l030@sam_data@row.names!="1024-4-17-37")


se.cor_rhizo_reduced_l030 <- spiec.easi(list(data_ITS_rhizo_reduced_l030, data_16S_rhizo_reduced_l030), 
                                        method='mb', 
                                        nlambda=100, 
                                        lambda.min.ratio=1e-1, 
                                        pulsar.params = list(thresh = 0.01))

se.cor_rhizo_reduced_hi30 <- spiec.easi(list(data_ITS_rhizo_reduced_hi30, data_16S_rhizo_reduced_hi30),
                                        method='mb',
                                        nlambda=100,
                                        lambda.min.ratio=1e-1,
                                        pulsar.params = list(thresh= 0.01))

se.cor_rhizo_full_l030 <- spiec.easi(list(data_ITS_rhizo_full_l030, data_16S_rhizo_full_l030),
                                     method='mb',
                                     nlambda=100,
                                     lambda.min.ratio=1e-1,
                                     pulsar.params = list(thresh= 0.01))

#remove sample not present in the ITS for reduced x LO30: 1024-4-17-37
data_16S_rhizo_full_hi30 <- subset_samples(data_16S_rhizo_full_hi30, data_16S_rhizo_full_hi30@sam_data@row.names!="GW9589-2-11-71")
se.cor_rhizo_full_hi30 <- spiec.easi(list(data_16S_rhizo_full_hi30, data_16S_rhizo_full_hi30),
                                     method='mb',
                                     nlambda=100,
                                     lambda.min.ratio=5e-2,
                                     pulsar.params = list(thresh= 0.01))

 library(RColorBrewer)
library(pals)
#make vectors to plot all shapes and colors the same
Kingdom_reduced_l030 <-unique(tax_table(physeq_combo_reduced_l030)[,"Kingdom"])
Kingdom_reduced_hi30 <-unique(tax_table(physeq_combo_reduced_hi30)[,"Kingdom"])
Kingdom_full_l030 <-unique(tax_table(physeq_combo_full_l030)[,"Kingdom"])
Kingdom_full_hi30 <-unique(tax_table(physeq_combo_full_hi30)[,"Kingdom"])
KingdomList<-do.call("rbind", list(Kingdom_reduced_l030, Kingdom_reduced_hi30, Kingdom_full_l030, Kingdom_full_hi30))%>%unique()
KingdomShape = c(16,17,15)
names(KingdomShape) = KingdomList


getPalette = colorRampPalette(brewer.pal(9, "Set1"))
Phylum_reduced_l030 <-unique(tax_table(physeq_combo_reduced_l030)[,"Phylum"])
Phylum_reduced_hi30 <-unique(tax_table(physeq_combo_reduced_hi30)[,"Phylum"])
Phylum_full_l030 <-unique(tax_table(physeq_combo_full_l030)[,"Phylum"])
Phylum_full_hi30 <-unique(tax_table(physeq_combo_full_hi30)[,"Phylum"])
PhylumList<-do.call("rbind", list(Phylum_reduced_l030, Phylum_reduced_hi30, Phylum_full_l030, Phylum_full_hi30))%>%unique()

PhylumPalette = getPalette(length(PhylumList))
names(PhylumPalette) = PhylumList

#plot networks
reduced_l030<-plot_network(ig_reduced_l030,physeq_combo_reduced_l030,type="taxa",color="Phylum",shape="Kingdom",label=NULL,point_size = 1)+scale_shape_manual(values=KingdomShape,drop=FALSE)+scale_color_manual(values=PhylumPalette,drop=FALSE)
reduced_hi30<-plot_network(ig_reduced_hi30,physeq_combo_reduced_hi30,type="taxa",color="Phylum",shape="Kingdom",label=NULL,point_size=1)+scale_shape_manual(values=KingdomShape)+scale_color_manual(values=PhylumPalette)
full_l030<-plot_network(ig_full_l030,physeq_combo_full_l030,type="taxa",color="Phylum",shape="Kingdom",label=NULL,point_size=1)+scale_shape_manual(values=KingdomShape)+scale_color_manual(values=PhylumPalette)
full_hi30<-plot_network(ig_full_hi30,physeq_combo_full_hi30,type="taxa",color="Phylum",shape="Kingdom",label=NULL,point_size=1)+scale_shape_manual(values=KingdomShape)+scale_color_manual(values=PhylumPalette)

reduced_l030

rhizo_network<-ggarrange(reduced_l030,reduced_hi30,full_l030,full_hi30,common.legend = TRUE,legend="bottom", labels=c("A","B","C","D"))
ggsave("network_rhizo.tiff",plot=rhizo_network, width=11,height=7, device="tiff",dpi=600)

# Export Network Information
taxa.rhizo_reduced_l030 <-data.frame(tax_table(physeq_combo_reduced_l030))
betaMat_rhizo_reduced_l030=as.matrix(symBeta(getOptBeta(se.cor_rhizo_reduced_l030)))
total_rhizo_reduced_l030 <- length(betaMat_rhizo_reduced_l030[betaMat_rhizo_reduced_l030!=0])/2 

taxa.rhizo_reduced_hi30 <-data.frame(tax_table(physeq_combo_reduced_hi30))
betaMat_rhizo_reduced_hi30=as.matrix(symBeta(getOptBeta(se.cor_rhizo_reduced_hi30)))
total_rhizo_reduced_hi30 <- length(betaMat_rhizo_reduced_hi30[betaMat_rhizo_reduced_hi30!=0])/2 

taxa.rhizo_full_l030 <-data.frame(tax_table(physeq_combo_full_l030))
betaMat_rhizo_full_l030=as.matrix(symBeta(getOptBeta(se.cor_rhizo_full_l030)))
total_rhizo_full_l030 <- length(betaMat_rhizo_full_l030[betaMat_rhizo_full_l030!=0])/2 

taxa.rhizo_full_hi30 <-data.frame(tax_table(physeq_combo_full_hi30))
betaMat_rhizo_full_hi30 =as.matrix(symBeta(getOptBeta(se.cor_rhizo_full_hi30)))

wtrhizo_reduced_l030 <- symBeta(getOptBeta(se.cor_rhizo_reduced_l030), mode="maxabs")
weightrhizo_reduced_l030 <-Matrix::summary(t(wtrhizo_reduced_l030))[,3]
igrhizo_reduced_l030_wt <-adj2igraph(getRefit(se.cor_rhizo_reduced_l030),edge.attr=list(weight=weightrhizo_reduced_l030),vertex.attr=list(name=c(taxa_names(data_16S_rhizo_reduced_l030),taxa_names(data_ITS_rhizo_reduced_l030))))
write.graph(igrhizo_reduced_l030_wt,file=file.path("spieceasi_reduced_l030.txt"),format="ncol") 

wtrhizo_reduced_hi30 <- symBeta(getOptBeta(se.cor_rhizo_reduced_hi30), mode="maxabs")
weightrhizo_reduced_hi30 <-Matrix::summary(t(wtrhizo_reduced_hi30))[,3]
igrhizo_reduced_hi30_wt <-adj2igraph(getRefit(se.cor_rhizo_reduced_hi30),edge.attr=list(weight=weightrhizo_reduced_hi30),vertex.attr=list(name=c(taxa_names(data_16S_rhizo_reduced_hi30),taxa_names(data_ITS_rhizo_reduced_hi30))))
write.graph(igrhizo_reduced_hi30_wt,file=file.path("spieceasi_reduced_hi30.txt"),format="ncol") 

wtrhizo_full_l030 <- symBeta(getOptBeta(se.cor_rhizo_full_l030), mode="maxabs")
weightrhizo_full_l030 <-Matrix::summary(t(wtrhizo_full_l030))[,3]
igrhizo_full_l030_wt <-adj2igraph(getRefit(se.cor_rhizo_full_l030),edge.attr=list(weight=weightrhizo_full_l030),vertex.attr=list(name=c(taxa_names(data_16S_rhizo_full_l030),taxa_names(data_ITS_rhizo_full_l030))))
write.graph(igrhizo_full_l030_wt,file=file.path("spieceasi_full_l030.txt"),format="ncol") 

wtrhizo_full_hi30 <- symBeta(getOptBeta(se.cor_rhizo_full_hi30), mode="maxabs")
weightrhizo_full_hi30 <-Matrix::summary(t(wtrhizo_full_hi30))[,3]
igrhizo_full_hi30_wt <-adj2igraph(getRefit(se.cor_rhizo_full_hi30),edge.attr=list(weight=weightrhizo_full_hi30),vertex.attr=list(name=c(taxa_names(data_16S_rhizo_full_hi30),taxa_names(data_ITS_rhizo_full_hi30))))
write.graph(igrhizo_full_hi30_wt,file=file.path("spieceasi_full_hi30.txt"),format="ncol") 
