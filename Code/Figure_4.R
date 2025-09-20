# For Reduced carbohydrate complexity and diversity 
# alter gut microbial structure independent of total intake Flores & Seekatz (2025)

# Figure 4

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

############################################################################################
Figure 4
Taxonomic Redundancy:
Files needed: metadata, final.asv.ASV.subsample.shared
############################################################################################
####################################### 
Figure: 4A Box plots of ASVs/genus
Tax redundancy @ each taxa level phyla, class...genus
Note: This is just for cecal samples...
#######################################

# Create df for phyla level redundancy
unique(cec_relabund$genus)

# Define a function that handles the process for each taxonomic level
add_missing_taxa <- function(cec_relabund.long, level) {
  # Filter for the specific level and remove rows with count 0
  filtered_df <- cec_relabund.long %>%
    filter(level == !!level, count != 0) %>%
    group_by(group, tx_group, taxon) %>%
    tally(name = "count")
  
  # Get the list of all distinct taxa for the current level
  all_taxa <- cec_relabund.long %>%
    filter(level == !!level) %>%
    distinct(taxon) %>%
    pull(taxon)
  
  # Create all combinations of group, tx_group, and taxon
  complete_df <- filtered_df %>%
    distinct(group, tx_group) %>%
    crossing(taxon = all_taxa) %>%
    left_join(filtered_df, by = c("group", "tx_group", "taxon")) %>%
    mutate(count = ifelse(is.na(count), 0, count))
  
  return(complete_df)
}

phy.redund.combo <- add_missing_taxa(cec_relabund.long, "phylum")
class.redund.combo <- add_missing_taxa(cec_relabund.long, "class")
ord.redund.combo <- add_missing_taxa(cec_relabund.long, "order")
fam.redund.combo <- add_missing_taxa(cec_relabund.long, "family")
gen.redund.combo <- add_missing_taxa(cec_relabund.long, "genus")

# Save redundancy dfs
# Create a list of your data frames
dFR_redund_list <- list(
  "Phylum" = phy.redund.combo,
  "Class" = class.redund.combo,
  "Order" = ord.redund.combo,
  "Family" = fam.redund.combo,
  "Genus" = gen.redund.combo
)

# Write to an Excel file with each data frame in a different sheet
#write_xlsx(dFR_redund_list, path = "dFR_tax_redundant_data.xlsx")

# Graph box plots for redundancy at each level
phy.redund <- read_excel("~/path/to/dFR_tax_redundant_data.xlsx",
                         sheet = "Phylum") 

# Step 1: Summarize the data to get the median (or mean) count per taxon
taxon_order <- phy.redund %>%
  group_by(taxon) %>%
  summarize(median_count = median(count + 1)) %>%
  arrange(desc(median_count)) %>%
  pull(taxon)

# Step 2: Reorder the taxon factor based on the summary statistic (median_count)
phy.redund <- phy.redund %>%
  mutate(taxon = factor(taxon, levels = taxon_order))

ggplot(phy.redund, aes(x=tx_group, y=log(count+1), fill=tx_group)) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, alpha = 0.5, fill = NA) +
  geom_jitter(show.legend = FALSE, width = 0.25, shape = 21, color = "black") +
  # stat_summary(fun = mean, show.legend = FALSE, geom = "crossbar",
  #              color = "black", width=0.6, size=0.5) +
  labs(x=NULL, 
       y="log(counts + 1)") +
  scale_x_discrete(limits=c("HCD","MCD","LCD")) +
  scale_fill_manual(name=NULL,
                    values=c("HCD"="orange", "MCD"="green", "LCD"="purple")) +
  facet_wrap("taxon") +
  theme_bw() +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        strip.text = element_text(color = "black"), 
        aspect.ratio = 1) 

# Define the function to plot redundancy at each taxonomic level
tax_redund_boxplot <- function(taxon_level) {
  
  # Read the data from the appropriate sheet (based on taxon level)
  tr.df <- read_excel("~/Desktop/thesis analysis dFR/dFR/dFR_tax_redundant_data.xlsx",
                      sheet = taxon_level)
  
  # Filter out unwanted taxa (Alphaproteobacteria and Bacteroidetes_unclassified)
  # These have 0 counts in cecal data
  if(taxon_level == "Class") {
    tr.df <- tr.df %>%
      filter(!taxon %in% c("Alphaproteobacteria", "Bacteroidetes_unclassified"))
  }
  
  if(taxon_level == "Order") {
    tr.df <- tr.df %>%
      filter(!taxon %in% c("Rickettsiales", "Bacteroidetes_unclassified"))
  }
  
  if(taxon_level == "Family") {
    tr.df <- tr.df %>%
      filter(!taxon %in% c("Rickettsiaceae", "Bacteroidetes_unclassified"))
  }
  
  if(taxon_level == "Genus") {
    tr.df <- tr.df %>%
      filter(!taxon %in% c("Bacteroidetes_unclassified", 
                           "Clostridium_XlVb",
                           "Faecalibacterium",
                           "Hungatella",
                           "Lacrimispora",
                           "Monoglobus",
                           "Paludicola",
                           "Rickettsia",
                           "Roseburia",
                           "Ruminococcus"))
  }
  
  # Step 1: Summarize the data to get the median count per taxon
  taxon_order <- tr.df %>%
    group_by(taxon) %>%
    summarize(median_count = median(count + 1)) %>%
    arrange(desc(median_count)) %>%
    pull(taxon)
  
  # Step 2: Reorder the taxon factor based on the summary statistic (median_count)
  tr.df <- tr.df %>%
    mutate(taxon = factor(taxon, levels = taxon_order))
  
  # Step 3: Create the box plot
  plot <- ggplot(tr.df, aes(x=tx_group, y=log2(count+1), fill=tx_group)) +
    geom_boxplot(show.legend = FALSE, outlier.shape = NA, alpha = 0.5, fill = NA) +
    geom_jitter(show.legend = FALSE, width = 0.25, shape = 21, color = "black") +
    labs(x=NULL, 
         y="log2 (counts + 1)") +
    scale_x_discrete(limits=c("HCD", "MCD", "LCD")) +
    scale_fill_manual(name=NULL,
                      values=c("HCD"="orange", "MCD"="green", "LCD"="purple")) +
    facet_wrap("taxon") +
    theme_bw() +
    theme(strip.placement = "outside",
          strip.background = element_blank(),
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          strip.text = element_text(color = "black"),
          aspect.ratio = 1)
  
  # Return the plot
  return(plot)
}

# Now call the function for each level (Phylum, Class, Order, etc.)
plot_phylum <- tax_redund_boxplot("Phylum")
plot_class <- tax_redund_boxplot("Class")
plot_order <- tax_redund_boxplot("Order")
plot_family <- tax_redund_boxplot("Family")
plot_genus <- tax_redund_boxplot("Genus") 

# Print the plots
print(plot_phylum)
#ggsave("phylum.tax_redund_boxplot_log2.pdf", height = 8, width = 6)
print(plot_class)
#ggsave("class.tax_redund_boxplot.pdf", height = 8, width = 6)
print(plot_order)
#ggsave("order.tax_redund_boxplot.pdf", height = 8, width = 6)
print(plot_family)
#ggsave("family.tax_redund_boxplot.pdf", height = 7, width = 7)
print(plot_genus)
#ggsave("genus.tax_redund_boxplot.pdf", height = 7, width = 7)

####################################################################
Statistics
####################################################################
phy.kw <- phy.redund %>%
  group_by(taxon) %>% 
  do({
    # perform Kruskal-Wallis test
    kw_result <- kruskal.test(count ~ tx_group, data = .)
    
    # Initialize result df 
    data.frame(
      taxon = unique(.$taxon),
      kw_statistic = kw_result$statistic,
      kw_p_value = kw_result$p.value,
      df = kw_result$parameter
    )
  }) 
print(phy.kw)

# Perform pairwise comparisons if Kruskal-Wallis is significant
# Filter for every KW sig tax
pw_filt <- phy.redund %>%
  filter(taxon == "Tenericutes") 
pairwise.wilcox.test(pw_filt$count, pw_filt$tx_group, p.adjust.method = "BH")

class.redund <- read_excel("~/path/to/dFR_tax_redundant_data.xlsx",
                           sheet = "Class") 
class.kw <- class.redund %>%
  group_by(taxon) %>% 
  do({
    # perform Kruskal-Wallis test
    kw_result <- kruskal.test(count ~ tx_group, data = .)
    
    # Initialize result df 
    data.frame(
      taxon = unique(.$taxon),
      kw_statistic = kw_result$statistic,
      kw_p_value = kw_result$p.value,
      df = kw_result$parameter
    )
  }) 
print(class.kw)

# Perform pairwise comparisons if Kruskal-Wallis is significant
# Filter for every KW sig tax
pw_filt <- class.redund %>%
  filter(taxon == "Mollicutes") 
pairwise.wilcox.test(pw_filt$count, pw_filt$tx_group, p.adjust.method = "BH")

order.redund <- read_excel("~/path/to/dFR_tax_redundant_data.xlsx",
                           sheet = "Order") 
order.kw <- order.redund %>%
  group_by(taxon) %>% 
  do({
    # perform Kruskal-Wallis test
    kw_result <- kruskal.test(count ~ tx_group, data = .)
    
    # Initialize result df 
    data.frame(
      taxon = unique(.$taxon),
      kw_statistic = kw_result$statistic,
      kw_p_value = kw_result$p.value,
      df = kw_result$parameter
    )
  }) 
print(order.kw)

# Perform pairwise comparisons if Kruskal-Wallis is significant
# Filter for every KW sig tax
pw_filt <- order.redund %>%
  filter(taxon == "Coriobacteriales") 
pairwise.wilcox.test(pw_filt$count, pw_filt$tx_group, p.adjust.method = "BH")

family.redund <- read_excel("~/path/to/dFR_tax_redundant_data.xlsx",
                            sheet = "Family") 
family.kw <- family.redund %>%
  group_by(taxon) %>% 
  do({
    # perform Kruskal-Wallis test
    kw_result <- kruskal.test(count ~ tx_group, data = .)
    
    # Initialize result df 
    data.frame(
      taxon = unique(.$taxon),
      kw_statistic = kw_result$statistic,
      kw_p_value = kw_result$p.value,
      df = kw_result$parameter
    )
  }) 
print(family.kw)

# Perform pairwise comparisons if Kruskal-Wallis is significant
# Filter for every KW sig tax
pw_filt <- family.redund %>%
  filter(taxon == "Muribaculaceae") 
pairwise.wilcox.test(pw_filt$count, pw_filt$tx_group, p.adjust.method = "BH")

#######################################################
genus.redund <- read_excel("~/path/to/dFR_tax_redundant_data.xlsx",
                           sheet = "Genus") 
genus.kw <- genus.redund %>%
  group_by(taxon) %>% 
  do({
    # perform Kruskal-Wallis test
    kw_result <- kruskal.test(count ~ tx_group, data = .)
    
    # Initialize result df 
    data.frame(
      taxon = unique(.$taxon),
      kw_statistic = kw_result$statistic,
      kw_p_value = kw_result$p.value,
      df = kw_result$parameter
    )
  }) 
print(genus.kw)

# Perform pairwise comparisons if Kruskal-Wallis is significant
# Filter for every KW sig tax
pw_filt <- genus.redund %>%
  filter(taxon == "Staphylococcaceae_unclassified") 
pairwise.wilcox.test(pw_filt$count, pw_filt$tx_group, p.adjust.method = "BH")

####################################### 
Figure 4B: Plots of taxonomic redundancy
#######################################

# Step 1: Create a summary of mean values by taxon
top_taxa <- genus.redund %>%
  filter(!taxon %in% c("Bacteroidetes_unclassified", 
                       "Clostridium_XlVb",
                       "Faecalibacterium",
                       "Hungatella",
                       "Lacrimispora",
                       "Monoglobus",
                       "Paludicola",
                       "Rickettsia",
                       "Roseburia",
                       "Ruminococcus")) %>%
  mutate(log2_count = log2(count + 1)) %>%
  group_by(taxon) %>%
  summarise(mean_log2_count = mean(log2_count)) %>%
  arrange(desc(mean_log2_count)) %>%
  slice_head(n = 31)  # Top 31 unique taxa

# Step 1: Create a summary of mean values by taxon
bottom_taxa <- genus.redund %>%
  filter(!taxon %in% c("Bacteroidetes_unclassified", 
                       "Clostridium_XlVb",
                       "Faecalibacterium",
                       "Hungatella",
                       "Lacrimispora",
                       "Monoglobus",
                       "Paludicola",
                       "Rickettsia",
                       "Roseburia",
                       "Ruminococcus")) %>%
  mutate(log2_count = log2(count + 1)) %>%
  group_by(taxon) %>%
  summarise(mean_log2_count = mean(log2_count)) %>%
  arrange(desc(mean_log2_count)) %>%
  slice_tail(n = 31)  # Bottom 31 unique taxa

# Step 2: Use those top taxa to filter main data for plot

genus.redund %>%
  filter(taxon %in% top_taxa$taxon) %>%
  mutate(log2_count = log2(count + 1)) %>%  # Log transform counts
  group_by(taxon, tx_group) %>%
  mutate(mean_log2_count = mean(log2_count)) %>%  # Calculate the mean per taxon
  ungroup() %>%
  mutate(tx_group = factor(tx_group, levels = c("HCD", "MCD", "LCD"))) %>%
  mutate(taxon = reorder(taxon, -mean_log2_count)) %>%  # Reorder taxa based on mean, put - sign in front if you want big-small
  ggplot(aes(x = taxon, y = log2(count + 1), color = tx_group)) +
  
  geom_point(
    aes(x = taxon, y = log2(count + 1), 
        color = tx_group),
    position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8),
    alpha = 0.6,
    shape = 16,
    size = 2,
    show.legend = FALSE
  ) +
  # Use stat_summary to calculate and plot mean and error bars
  stat_summary(fun.data = "mean_sdl", aes(x = taxon, y = log2(count + 1)), 
               fun.args = list(mult = 1), 
               position = position_dodge(width = 0.8),
               geom = "errorbar", width = 0.4, show.legend = FALSE) +  # Min-Max error bars
  
  stat_summary(fun = "mean", aes(x = taxon, y = log2(count + 1), fill = tx_group), 
               geom = "point", size = 3,
               alpha = 0.5,
               shape = 21, 
               color = "red",
               show.legend = FALSE,
               position = position_dodge(width = 0.8)) +  # Mean points
  
  labs(x = "", y = "log2 (count +1)") +  # Axis labels
  scale_x_discrete(labels = function(x) gsub("unclassified", "`", gsub("_", " ", x))) +  # Clean up taxon names
  scale_fill_manual(name = NULL, values = c("HCD" = "orange", "MCD" = "green", "LCD" = "purple")) +  # Color for tx_group
  scale_color_manual(name = NULL, values = c("HCD" = "orange", "MCD" = "green", "LCD" = "purple")) +  # Color for tx_group
  
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", angle = 90, hjust = 1),  # Rotate x-axis labels for readability
        axis.text.y = element_text(color = "black"),
        strip.text = element_text(color = "black")) 
# + coor_flip() if you want tall plot
#ggsave("top31_genus_redund_cat_plot.pdf", width = 15, height = 5)

# PLOT  bottom 31 genus 
genus.redund %>%
  filter(taxon %in% bottom_taxa$taxon) %>%
  mutate(log2_count = log2(count + 1)) %>%  # Log transform counts
  group_by(taxon, tx_group) %>%
  mutate(mean_log2_count = mean(log2_count)) %>%  # Calculate the mean per taxon
  ungroup() %>%
  mutate(tx_group = factor(tx_group, levels = c("HCD", "MCD", "LCD"))) %>%
  mutate(taxon = reorder(taxon, -mean_log2_count)) %>%  # Reorder taxa based on mean, put - sign in front if you want big-small
  ggplot(aes(x = taxon, y = log2(count + 1), color = tx_group)) +
  
  geom_point(
    aes(x = taxon, y = log2(count + 1), 
        color = tx_group),
    position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8),
    alpha = 0.6,
    shape = 16,
    size = 2,
    show.legend = FALSE
  ) +
  
  # Use stat_summary to calculate and plot mean and error bars
  stat_summary(fun.data = "mean_sdl", aes(x = taxon, y = log2(count + 1)), 
               fun.args = list(mult = 1), 
               position = position_dodge(width = 0.8),
               geom = "errorbar", width = 0.4, show.legend = FALSE) +  # Min-Max error bars
  
  stat_summary(fun = "mean", aes(x = taxon, y = log2(count + 1), fill = tx_group), 
               geom = "point", size = 3,
               alpha = 0.5,
               shape = 21, 
               color = "red",
               show.legend = FALSE,
               position = position_dodge(width = 0.8)) +  # Mean points
  
  labs(x = "", y = "log2 (count +1)") +  # Axis labels
  scale_x_discrete(labels = function(x) gsub("unclassified", "`", gsub("_", " ", x))) +  # Clean up taxon names
  scale_fill_manual(name = NULL, values = c("HCD" = "orange", "MCD" = "green", "LCD" = "purple")) +  # Color for tx_group
  scale_color_manual(name = NULL, values = c("HCD" = "orange", "MCD" = "green", "LCD" = "purple")) +  # Color for tx_group
  
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", angle = 90, hjust = 1),  # Rotate x-axis labels for readability
        axis.text.y = element_text(color = "black"),
        strip.text = element_text(color = "black")) 
# + coor_flip() if you want tall plot
#ggsave("bottom31_genus_redund_cat_plot.pdf", width = 15, height = 5)





