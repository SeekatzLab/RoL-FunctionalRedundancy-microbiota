# For Reduced carbohydrate complexity and diversity 
# alter gut microbial structure independent of total intake Flores & Seekatz (2025)

# Figure 3

library(readxl)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(ggpubr) # ggarange for plots
library(ggbreak) # customize ggplot
library(ggtext)
library(vegan) # Ordination methods, diversity analysis and other functions for community and vegetation ecologists.
library(lme4) # for stats
library(lmerTest)
library(car)
library(patchwork) # customize multiple plots into 1 figure
library(FSA) # for stats, specifically dunn's post hoc test
library(ggsignif) # add significance brackets
library(Maaslin2)
library(broom) #stats
library(pheatmap) # pretty heatmaps
library(circlize)
library(writexl) # write excel files
library(rstatix) # for stat tests
library(ggbeeswarm)

##########################################################################
Figure 3 Relab, Differential Abundance & Heat Tree
Note: I am working with mothur data here 
Files needed: metadata, final.asv.ASV.cons.taxonomy, final.asv.ASV.subsample.shared [This file has cols: 
                                                          label, Group, numOtus, ASV0001, etc]
#########################################################################
# Read in feature table
############################################################################################################

# Filtered table done for Figure 2

cec_relabund <- inner_join(meta.cecal, filtered_asv_counts, by="group") %>%
  inner_join(., taxonomy_asv, by="otu") %>% # ., takes data coming through the pipeline
  group_by(group) %>%
  dplyr::mutate(sum_count = sum(count)) %>% 
  dplyr::mutate(rel_abund = count / sum_count) 

cec_relabund.long <- inner_join(meta.cecal, filtered_asv_counts, by="group") %>%
  inner_join(., taxonomy_asv, by="otu") %>% # ., takes data coming through the pipeline
  group_by(group) %>%
  dplyr::mutate(sum_count = sum(count)) %>% 
  dplyr::mutate(rel_abund = count / sum_count) %>%
  ungroup %>%
  pivot_longer(cols=c("kingdom", "phylum", "class", "order", "family", "genus", "otu"),
               names_to = "level",
               values_to = "taxon") %>%
  dplyr::mutate(tx_group = factor(tx_group, 
                                  levels=c("LCD", "MCD", "HCD"))) 

# Can run the following to check relab was calculated correctly
cec_relabund %>% 
  dplyr::group_by(group) %>% 
  dplyr::summarise(sum_relab = sum(rel_abund))

unique(cec_relabund$phylum)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Color Scheme %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Make a new df with just genus and phylum from meta.OTU2 data
genus <- unique(cec_relabund[c("genus", "phylum")])
# From this df it looks like we have a grand total of 72 unique genera

# Find the unique phyla in this data set
unique(genus$phylum)
# And we find that we have 7 unique phyla...same as above
# Check to see the total abundance of all phyla across samples
cec_relabund %>%
  group_by(phylum) %>%
  summarise(sum = sum(count)) %>% arrange(desc(sum))

# lets re-order based on phylum 
genus$phylum <- factor(genus$phylum, levels = c("Firmicutes","Bacteroidetes","Verrucomicrobia",
                                                "Actinobacteria", "Tenericutes", 
                                                "Bacteria_unclassified", "Proteobacteria")) 

genus <- genus[order(genus$phylum),]
genus <- group_by(genus, phylum)

# if you need to check how many of each group 
table(genus$phylum)

# define colors 
firm<-colorRampPalette(c("blue","dodgerblue4","dodgerblue1","deepskyblue4","deepskyblue1","skyblue3","skyblue","mediumblue","steelblue4","steelblue1","royalblue4","royalblue1","slateblue4","lightskyblue","lightskyblue4","cornflowerblue"))(n=50)
bac<-c("darkolivegreen","darkolivegreen1","darkolivegreen2","darkolivegreen3","darkolivegreen4","green3","lightgreen","seagreen", "lightgreen")
pro<-c("gold", "gold1","gold2","gold3","gold4")
actino<-c("tan","salmon1","tan1","salmon2","tan2")
verruco<-c("hotpink")
tener<-c("purple3")
unclass<-c("grey")

Color <- c(firm, bac, verruco, actino,  tener, unclass, pro)
df <- as.data.frame(Color)
df <- cbind(genus, df)

# Count unique genus for the phylum present in our sample
aggregate(data = genus,                # Applying aggregate
          genus ~ phylum,
          function(genus) length(unique(genus))) %>% print() # Print counts

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of Color Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####################################### 
Figure: 3A Relative Abundance Bar Graph
#######################################
cec_relabund %>% 
  group_by(phylum) %>%
  summarise(sum_relab = sum(rel_abund)) %>%
  arrange(desc(sum_relab)) %>% print()

# Plot cecal data
cec_relabund %>%
  mutate(tx_group = ifelse(tx_group %in% c("LCD"),
                           paste0("x", tx_group),
                           tx_group)) %>%
  ggplot(., aes(x = sampleid, y = rel_abund, fill=factor(genus),
                group = factor(phylum, 
                               levels = c("Bacteria_unclassified","Proteobacteria","Tenericutes",
                                          "Actinobacteria","Verrucomicrobia","Bacteroidetes",
                                          "Firmicutes"
                               )))) +
  #facet_wrap('tx_group', scales = "free_x") +
  facet_grid(~ tx_group + sex, scales = "free_x", space = "free_x") +
  geom_col() +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(labels = scales::percent) +
  ylab('Relative Abundance %') +
  xlab('') +
  scale_fill_manual(values = setNames(df$Color, df$genus)) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none", axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside") 

####################################################
Figure 3B Differential abundance heatmap (Maaslin2)
Files needed: metadata, final.asv.ASV.cons.taxonomy, final.asv.ASV.subsample.shared [This file has cols: 
                                                          label, Group, numOtus, ASV0001, etc]
###################################################

# Format metadata, asv_counts, and taxa into appropriate df to use Maaslin
# For the version we are using we need df that has cols ID, ASV_Bac1, ASV_Bac2, etc

Maaslin_tax <- inner_join(meta.cecal, filtered_asv_counts, by="group") %>%
  inner_join(., taxonomy_asv, by="otu") %>% # ., takes data coming through the pipeline
  unite(taxa, c("otu", "genus")) %>%
  group_by(group) %>%
  dplyr::mutate(sum_count = sum(count)) %>% 
  dplyr::mutate(rel_abund = count / sum_count) %>%
  ungroup() %>%
  select(group, taxa, rel_abund) %>%
  dplyr::rename(ID = group) %>%
  pivot_wider(names_from = taxa, values_from = rel_abund)
#write_tsv(Maaslin_tax, "/path/to/Maaslin_tax.tsv")

Masslin_meta <- meta.cecal %>% dplyr::rename(ID = group)
#write_tsv(Masslin_meta, "/path/to/Maaslin_meta.tsv")

input_data = "/path/to/Maaslin_tax.tsv" # The abundance table file
input_data
input_metadata = "/path/to/Maaslin_meta.tsv" # The metadata table file
input_metadata

#Saving inputs as data frames
df_input_data = read.table(file = input_data, header = TRUE, sep = "\t",
                           row.names = 1,
                           stringsAsFactors = FALSE)
df_input_data[1:5, 1:5]
df_input_metadata = read.table(file = input_metadata, header = TRUE, sep = "\t",
                               row.names = 1,
                               stringsAsFactors = FALSE)
df_input_metadata[1:5, ]

df_input_metadata$tx_modified = factor(df_input_metadata$tx_group,
                                       levels = c("HCD", "MCD", "LCD"))

fit_data = Maaslin2(
  input_data = df_input_data, 
  input_metadata = df_input_metadata, 
  output = "dFR_Masoutput", 
  fixed_effects = c("tx_modified"),
  reference = c("tx_modified,HCD")) #  adding a space between the variable and level might result in the wrong reference level being used

fit_data2 = Maaslin2(
  input_data = df_input_data, 
  input_metadata = df_input_metadata, 
  output = "dFR_Masoutput_min_abund", 
  fixed_effects = c("tx_modified"),
  min_abundance = 0.001)

# Read in significant results from Masslin output 
# For this we will use output from results from our min abundance run and join with rel_abund data

maas_sig_results <- read_tsv("/path/to/significant_results.tsv")

sig_sig_tax <- Maaslin_tax %>% column_to_rownames("ID" ) %>%
  t(.) %>% as.data.frame() %>% 
  rownames_to_column("feature") %>% 
  inner_join(., maas_sig_results, by = "feature") %>%
  filter(qval < 0.005) 

# Remove the treatment_group column before reshaping
maas.sig.rel.abund <- sig_sig_tax %>% select(feature, starts_with("cecal")) %>%
  pivot_longer(
    cols = starts_with("cecal"),  # Adjust if seqIDs have different naming pattern
    names_to = "seqID",
    values_to = "rel_abund"
  ) %>% mutate(
    tx_group = case_when(
      grepl("^cecal_901_", seqID) ~ "HCD",   # Matches any code starting with "901-"
      grepl("^cecal_902_", seqID) ~ "HCD",   
      grepl("^cecal_903_", seqID) ~ "MCD",   # Example for other patterns
      grepl("^cecal_904_", seqID) ~ "MCD",   
      grepl("^cecal_905_", seqID) ~ "LCD",
      grepl("^cecal_906_", seqID) ~ "LCD",
      grepl("^cecal_907_", seqID) ~ "HCD",
      grepl("^cecal_908_", seqID) ~ "HCD",
      grepl("^cecal_909_", seqID) ~ "HCD",
      grepl("^cecal_910_", seqID) ~ "HCD",
      grepl("^cecal_911_", seqID) ~ "MCD",
      grepl("^cecal_912_", seqID) ~ "MCD",
      grepl("^cecal_913_", seqID) ~ "MCD",
      grepl("^cecal_914_", seqID) ~ "MCD",
      grepl("^cecal_915_", seqID) ~ "LCD",
      grepl("^cecal_916_", seqID) ~ "LCD",
      grepl("^cecal_917_", seqID) ~ "LCD",
      grepl("^cecal_918_", seqID) ~ "LCD",
      TRUE ~ NA_character_            # Default case if no conditions match
    )
  )

sig.feat <- maas.sig.rel.abund %>% 
  group_by(tx_group, feature) %>%
  mutate(sum_relabund = sum(rel_abund)) %>%
  filter(sum_relabund > 0.1)

sig.feat <- as.data.frame(sig.feat$feature) %>%
  rename(feature = 1) %>%
  distinct()

# Plot 3B

# with pheatmap to show phyla sig.diffs belong to
ASV_logRA <- inner_join(sig.feat, maas.sig.rel.abund, by = "feature", relationship = "many-to-many") %>%
  mutate(tx_group = factor(tx_group, 
                           levels=c("HCD", "MCD", "LCD"))) %>% 
  mutate(relab_plus1 = log2(rel_abund + 1)) %>%
  select(-rel_abund) %>%
  group_by(tx_group, feature) %>%
  mutate(mean_logRA = mean(relab_plus1)) %>%
  ungroup() %>%
  select(tx_group, feature, mean_logRA) %>%
  distinct() %>%
  pivot_wider(names_from = "feature",
              values_from = "mean_logRA") 

# Convert tibble to matrix with first row as colnames and first column as rownames
ASV_logRA_mat <- ASV_logRA %>%
  column_to_rownames(var = "tx_group") %>%  # Set the first column as rownames
  as.matrix()                  # Convert to matrix

# Create a vector of column names from ASV_logRA_mat
asv_columns <- colnames(ASV_logRA_mat)

# Create matrix that has phylum as col_name and rows as ASVs
ann_tax <- taxonomy_asv %>% mutate(feature = str_c(otu, genus, sep = "_")) %>%
  select(feature, phylum) %>%
  distinct() %>%
  filter(feature %in% asv_columns) %>%
  column_to_rownames(var = "feature")

# Ensure the 'phylum' column in 'ann_tax' is a factor with the correct order
ann_tax$phylum <- factor(ann_tax$phylum, levels = c("Firmicutes", 
                                                    "Bacteroidetes",
                                                    "Tenericutes",
                                                    "Actinobacteria",
                                                    "Verrucomicrobia",
                                                    "Proteobacteria",
                                                    "Bacteria_unclassified"))
# Reorder the entire ann_tax data frame by the 'phylum' column (as a data frame)
ann_tax <- ann_tax[order(ann_tax$phylum), , drop = FALSE]

# Now reorder the columns of ASV_logRA_mat to match the order of ann_tax
ASV_logRA_mat_reordered <- ASV_logRA_mat[, rownames(ann_tax)]
# remove _ in tax names, this is not matching ann_tax now so...
#colnames(ASV_logRA_mat_reordered) <- gsub("_", " ", colnames(ASV_logRA_mat_reordered))

ann_colors <- list(
  phylum = c("Firmicutes" = "darkblue",
             "Bacteroidetes" = "darkgreen",
             "Verrucomicrobia" = "hotpink",
             "Actinobacteria" = "salmon2",
             "Proteobacteria" = "gold",
             "Tenericutes" = "purple3",
             "Bacteria_unclassified" = "grey"))

# Define a smoother color palette with more colors
pheat_palette <- colorRampPalette(c("white", "#E6E6E6", "#CCCCCC", "black"))(20)

# Create custom breaks for smoother transitions
# Create 20 breaks, manually
pheat_breaks <- c(0.0, 0.00002, 0.0001, 0.0002, 0.0005, 0.0006,
                  0.0008, 0.002, 0.004, 0.008, 0.01, 0.02,
                  0.03, 0.04, 0.05, 0.06, 0.08, 0.15,
                  0.2, 0.3)

heat_plot <- pheatmap(ASV_logRA_mat_reordered,
                      cluster_cols=F, 
                      cluster_rows=F,
                      color = pheat_palette,
                      border_color="grey20",
                      #scale = "column", # 3 options: none, row, and column
                      cellheight = 20,
                      cellwidth = 20,
                      fontsize = 8,
                      annotation_col = ann_tax,
                      annotation_colors = ann_colors,
                      legend = TRUE, 
                      annotation_legend = FALSE,
                      breaks = pheat_breaks)

pdf(file = "/path/to/dFR_Fig3B_DA_fixed.pdf",   # The directory you want to save the file in
    width = 8.5, # The width of the plot in inches
    height = 5) # The height of the plot in inches
heat_plot
dev.off()

###############################################################################
Figure 3C
Heat Tree
Files needed: 
###############################################################################
#Archived on CRAN...except for conditionz 
install.packages("https://cran.r-project.org/src/contrib/Archive/bold/bold_1.3.0.tar.gz", repos = NULL, type = "source")
install.packages("https://cran.r-project.org/src/contrib/conditionz_0.1.0.tar.gz", repos = NULL, type = "source")
install.packages("https://cran.r-project.org/src/contrib/Archive/taxize/taxize_0.9.99.tar.gz", repos = NULL, type = "source")
install.packages("https://cran.r-project.org/src/contrib/Archive/metacoder/metacoder_0.3.7.tar.gz", repos = NULL, type = "source")

library(metacoder)

# Read in metadata
metadata <- read_excel("~/path/to/dFR_illumina16s_meta_metadata.xlsx", na="NA") %>%
  mutate(seqID = str_replace_all(seqID, "-", "_")) %>%
  mutate(SampleID = str_replace_all(SampleID, "-", "_"))

# Filter out cecal meta data
# rename(new1=old1)
meta.cecal <- metadata %>% filter(str_detect(seqID, 'cecal')) %>%
  rename_all(tolower) %>% dplyr::rename(group = seqid) 

# Read in feature table
############################################################################################################
asv_counts.long <- read_tsv("~/path/to/final.asv.ASV.subsample.shared") %>% 
  select(Group, starts_with("ASV")) %>% 
  pivot_longer(-Group, names_to="otu", values_to="count") %>%
  rename_all(tolower) 

asv_counts <- read_tsv("/path/to/final.asv.ASV.subsample.shared") %>% 
  select(Group, starts_with("ASV")) %>% 
  filter(str_detect(Group, 'cecal')) %>%
  rename(group = Group) %>% column_to_rownames("group") %>%
  t(.) %>% as.data.frame() %>%
  rownames_to_column("otu")

#********************** Filter out 10 reads or less *****************
# Step 1: Filtering out ASVs with less than 10 reads
filtered_asv_df <- asv_counts.long %>%
  group_by(otu) %>%
  summarise(total_count = sum(count)) %>%
  filter(total_count > 9)
# Step 2: Extract the taxa that meet the condition
filtered_taxa <- filtered_asv_df$otu
# Step 3: Subset the OTU table to keep only the filtered taxa, & fix the wrong sample name
filtered_asv_counts <- asv_counts %>%
  filter(otu %in% filtered_taxa) #%>%
#mutate(across(everything(), ~str_replace_all(., "d57_905_1_2", "d22_905_1"))) %>%

# Read in taxonomy
############################################################################################################
taxonomy_asv <- read_tsv("~/path/to/final.asv.ASV.cons.taxonomy") %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy,"\\(\\d+\\)","")) %>% #pattern i want to find is something that starts with ( then we want to match a digit (d) + = match 1 or more digits
  mutate(taxonomy = str_replace(taxonomy,";$","")) %>% # $ match the current character at the end of the string
  separate(taxonomy, 
           into=c("kingdom", "phylum", "class", "order", "family", "genus"), 
           sep=";") %>%
  mutate(
    phylum = paste0("p__", phylum),
    class = paste0("c__", class),
    order = paste0("o__", order),
    family = paste0("f__", family), 
    genus = paste0("g__", genus) 
  ) %>%
  mutate(root = "r__Root") %>% select(-kingdom) %>% 
  select(otu, last_col(), everything()) %>%
  unite("lineage", root, phylum, class, order, family, genus, sep = ";")

# Join taxa and asv_counts together
############################################################################################################
dfr_asvs <- inner_join(filtered_asv_counts, taxonomy_asv, "otu") %>%
  select(otu, last_col(), everything())

obj <- parse_tax_data(dfr_asvs,
                      class_cols = "lineage", # the column that contains taxonomic information
                      class_sep = ";", # The character used to separate taxa in the classification
                      class_regex = "^(.+)__(.+)$", # Regex identifying where the data for each taxon is
                      class_key = c(tax_rank = "info", # A key describing each regex capture group
                                    tax_name = "taxon_name"))
print(obj)

################################
Accounting for uneven sampling
################################
# So far we’ve used raw counts, but people typically work with rarefied counts 
# or proportions to try to avoid the possibility of sampling depth biasing the results. 
# Here we use the function calc_obs_props to divide each sample’s counts by 
# the total number of counts observed for each sample, resulting in a proportion.
obj$data$tax_data <- calc_obs_props(obj, "tax_data")
print(obj)

#################################
Getting per taxon information
################################
# Currently, we have values for the abundance of each OTU, not each taxon. 
# To get information on the taxa, we can sum the abundance per-taxon and 
# add the results to the taxmap object in a new table:
obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data",
                                       cols = meta.cecal$sampleid)
# We can also easily calculate the number of samples that have reads for each taxon:
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = meta.cecal$tx_group, cols = meta.cecal$sampleid)
######################################################
Plotting taxonomic data
######################################################

# Now that we have per-taxon information (The tax_abund and tax_occ tables), 
# we can plot the information using heat trees. Heat trees are what we call 
# taxonomic trees in which the size and color of tree parts correspond to 
# some statistic of interest. The code below plots the number of “HCD” 
# samples that have reads for each taxon as the size of each taxon. 
# It also plots the number of OTUs assigned to each taxon in the overall dataset as color.

set.seed(8597) # This makes the plot appear the same each time it is run 

hcd.tree <- heat_tree(obj, 
                      node_label = taxon_names,
                      node_size = n_obs,
                      node_color = HCD, 
                      node_color_range = c("#e3d9d7","#a5a6ad","#858793", "#666878"),
                      node_size_axis_label = "ASV count",
                      node_color_axis_label = "Samples with reads",
                      layout = "davidson-harel", # The primary layout algorithm
                      initial_layout = "reingold-tilford") #, # The layout algorithm that initializes node locations
#output_file = "differential_heat_tree.pdf") 
# This should check out because our filtered_taxa contain 1,399 ASVs
#dev.off()






