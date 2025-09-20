# For Reduced carbohydrate complexity and diversity 
# alter gut microbial structure independent of total intake Flores & Seekatz (2025)

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
Figure 2
Taxonomic Redundancy:
Files needed: metadata, final.asv.ASV.subsample.shared, final.asv.ASV.cons.taxonomy
############################################################################################

# Read in metadata

metadata <- read_excel("~/path/to/dFR_illumina16s_meta_metadata.xlsx", na="NA") %>%
  mutate(seqID = str_replace_all(seqID, "-", "_")) %>%
  mutate(SampleID = str_replace_all(SampleID, "-", "_"))

# Filter out cecal meta data
# rename(new1=old1)
meta.cecal <- metadata %>% filter(str_detect(seqID, 'cecal')) %>%
  rename_all(tolower) %>% dplyr::rename(group = seqid)

# Now fecal data
meta.fecal <- metadata %>% filter(!str_detect(seqID, 'cecal')) %>%
  rename_all(tolower) %>% dplyr::rename(group = seqid) %>% filter(!is.na(tx_group))

# Filter out rows with NAs
meta.nona <- metadata %>% rename_all(tolower) %>% rename(group = seqid) %>% filter(!is.na(tx_group))

# Read in feature table
############################################################################################################
asv_counts <- read_tsv("~/path/to/final.asv.ASV.subsample.shared") %>% 
  select(Group, starts_with("ASV")) %>%
  pivot_longer(-Group, names_to="otu", values_to="count") %>%
  rename_all(tolower) 

length(unique(asv_counts$otu))

#********************** Filter out 10 reads or less *****************
# Step 1: Filtering out ASVs with less than 10 reads
filtered_asv_df <- asv_counts %>%
  group_by(otu) %>%
  summarise(total_count = sum(count)) %>%
  filter(total_count > 9)
# Step 2: Extract the taxa that meet the condition
filtered_taxa <- filtered_asv_df$otu
# Step 3: Subset the OTU table to keep only the filtered taxa, & fix the wrong sample name
filtered_asv_counts <- asv_counts %>%
  filter(otu %in% filtered_taxa) #%>%
#mutate(across(everything(), ~str_replace_all(., "d57_905_1_2", "d22_905_1"))) 

# Read in taxonomy
# As a Note: d57_905_1_2 should be d57_905_1 wrong in group col correct in sampleid col
############################################################################################################
taxonomy_asv <- read_tsv("~/path/to/final.asv.ASV.cons.taxonomy") %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy,"\\(\\d+\\)","")) %>% #pattern i want to find is something that starts with ( then we want to match a digit (d) + = match 1 or more digits
  mutate(taxonomy = str_replace(taxonomy,";$","")) %>% # $ match the current character at the end of the string
  separate(taxonomy, 
           into=c("kingdom", "phylum", "class", "order", "family", "genus"), 
           sep=";")

asv_relabund <- inner_join(meta.nona, filtered_asv_counts, by="group") %>%
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

#*************************************************************************************************************************
Figure 2A & B  ALPHA DIVERSITY 
#*************************************************************************************************************************
# Calculate alpha diversity stats using the Vegan package from the filtered asv table we made

alpha_div_veg <- filtered_asv_counts %>% 
  group_by(group) %>% 
  summarise(sobs = specnumber(count),
            shannon = diversity(count, index = "shannon"),
            simpson = diversity(count, index = "simpson"),
            invsimpson = diversity(count, index = "invsimpson"),
            my_simpson = my_simpson(count),
            my_invsimpson = 1/my_simpson,
            abundance = sum(count))

# join metadata a diversity metric df together - Cecal data
joined_ASVdiv <- inner_join(meta.cecal, alpha_div_veg, by = 'group') %>%
  mutate(day = case_when( # Trying to hack the endpoint data so I can append to longitudinal graph...note: all endpoint data @ d59!
    tx_group == "HCD" ~ 62,  # If treatment_group is "HCD", set day to 59
    tx_group == "MCD" ~ 67,  # If treatment_group is "MCD", set day to 62
    tx_group == "LCD" ~ 72, # If treatment_group is "LCD", set day to 65
    TRUE ~ day  # Keep the existing value of 'day' if neither condition is met
  ))

ASVcecal.richness <- ggplot(joined_ASVdiv, aes(x=tx_group, y=sobs)) +
  geom_jitter(aes(color=tx_group), show.legend = FALSE) +
  geom_boxplot(aes(fill=tx_group, alpha = 0.5), show.legend = FALSE, width=0.5) +
  scale_fill_manual(values=c("HCD"="orange", "MCD"="green", "LCD"="purple")) +
  scale_color_manual(values=c("HCD"="orange", "MCD"="green", "LCD"="purple")) +
  scale_x_discrete(limits=c("HCD","MCD","LCD")) +
  scale_y_continuous(breaks=c(125,150,175,200,225)) +
  labs(y = "Richness", x="") +
  theme_classic() 

# Can also do Shannon and Simpson
ASVcecal.invsimpson <- ggplot(joined_ASVdiv, aes(x=tx_group, y=invsimpson)) + 
  geom_jitter(aes(color=tx_group), show.legend = FALSE) +
  geom_boxplot(aes(fill=tx_group, alpha = 0.5), show.legend = FALSE, width=0.5, outlier.shape = NA) +
  scale_fill_manual(values=c("HCD"="orange", "MCD"="green", "LCD"="purple")) +
  scale_color_manual(values=c("HCD"="orange", "MCD"="green", "LCD"="purple")) +
  scale_x_discrete(limits=c("HCD","MCD","LCD")) +
  labs(y = "Inverse Simpson", x="") +
  theme_classic() 

# join metadata a diversity metric df together - Fecal data
joined_fecal.ASVdiv <- inner_join(meta.fecal, alpha_div_veg, by = 'group')

invsimp_longWcecal <- ggplot(joined_fecal.ASVdiv, aes(x = day, y = invsimpson, color = tx_group, fill = tx_group)) +
  geom_jitter(width=0.25, alpha=0.4, size = 1, na.rm=TRUE) +
  scale_fill_manual(values=c("HCD"="orange", "MCD"="green", "LCD"="purple")) +
  scale_color_manual(values=c("HCD"="orange", "MCD"="green", "LCD"="purple")) +
  stat_summary(fun.data=median_hilow, na.rm= TRUE, geom="line", linewidth=1,
               fun.args = list(conf.int=0.50)) +
  stat_summary(fun.data = median_hilow, geom = "ribbon", 
               fun.args = list(conf.int = 0.25), alpha = 0.2, color = NA) +
  scale_y_continuous(breaks=c(5,10,15,20,25,30)) +
  scale_x_continuous(limits=c(0, 74), breaks=c(0, 1, 7, 14, 21, 28, 35, 42, 50, 57)) +
  labs(x="day", y="Inverse Simpson") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  theme_classic() +
  geom_vline(xintercept = 59.5, linetype = "solid", color = "black", size = 0.5) +
  geom_boxplot(data = joined_ASVdiv, aes(fill=tx_group, alpha = 0.5), show.legend = FALSE, width=3) +
  geom_jitter(data = joined_ASVdiv, width=0.25, alpha=0.4, size = 1, na.rm=TRUE) +
  theme(legend.position = "none")

#*************************************************************************************************************************
Figure 2 Figure 2A & B  Statistics
#*************************************************************************************************************************

# Step 1: Pivot the data longer for all alpha diversity metrics
alpha_metrics <- c("sobs", "shannon", "simpson", "invsimpson", "abundance")

# Need to factor tx_group so that x axis will be in correct order
order_x <- c("HCD", "MCD", "LCD")

joined_ASV_div_long <- joined_ASVdiv %>%
  select(-c("my_simpson", "my_invsimpson")) %>%
  pivot_longer(
    cols = all_of(alpha_metrics),
    names_to = "alpha_metric",
    values_to = "value") %>%
  mutate(tx_group = factor(tx_group, levels = order_x))

# Step 2: Summarize or run stats by tx_group for each alpha metric
# Example: Summary stats
adiv_summary_stats <- joined_ASV_div_long %>%
  group_by(alpha_metric, tx_group) %>%
  summarise(
    n = n(),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    IQR = IQR(value, na.rm = TRUE),
    .groups = "drop"
  )

# Step 3: Run ANOVA or Kruskal-Wallis for each alpha metric
# Example with Kruskal-Wallis (non-parametric test)
adiv_kruskal_results <- joined_ASV_div_long %>%
  group_by(alpha_metric) %>%
  kruskal_test(value ~ tx_group)

#*************************************************************************************************************************
Figure 2C, D and E  BETA DIVERSITY 
#*************************************************************************************************************************

filtered_asv_mat <- filtered_asv_counts %>%
  pivot_wider(
    names_from = otu,  # Column names will be from the 'asv' column
    values_from = count  # Values will be from the 'count' column
  ) %>%
  column_to_rownames("group")

filt_asv_tib <- filtered_asv_mat %>% rownames_to_column("group")

bc_dist = vegdist(filtered_asv_mat, method = "bray")
PCOA = cmdscale(bc_dist, eig = TRUE, k = 2)
eigenvectors = PCOA$eig / sum(PCOA$eig)
eigenvectors

bc_all <- PCOA[["points"]] 
bc_all <- as.data.frame(bc_all) %>% 
  rename(PC1 = V1, PC2 = V2) 
bc_all <- bc_all %>% rownames_to_column("group")

beta_join <- inner_join(x = meta.nona, y = bc_all, by = 'group')

# make new column 
beta_join$new_group <- beta_join$tx_group
# populate new_group column with 'Pre' for those samples that were 'Pre"
beta_join$new_group[beta_join$tx_status %in% c("Pre")] <- "pre"

# Plot PCoA 
# Figure 2C
ggplot(beta_join, aes(x = PC1, y = PC2, fill = new_group)) +
  geom_point(aes(shape = type, stroke = 1), size = 3) + 
  scale_shape_manual(values=c(21, 24)) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") + 
  # to find the percentages for your axis you will need to look at your the eigenvalues in PCOA list 
  # Be sure to double check that you are looking at the right file!
  xlab("PC1 (15.97%)") +
  ylab("PC2 (11.60%)") +
  scale_fill_manual(values = c("HCD" = "orange", 'MCD' = 'green', 'LCD' = 'purple', 'pre' = 'grey')) #+

# Now lets see if we can see what species are driving the distinct clustering we see in our PCoA plot above
# First let's look at all of the otus (fecal and cecal), then we can be a bit more selective
# All fecal + cecal
PCoA_shared <- inner_join(filtered_asv_counts, bc_all, by = "group")

cor_x <- PCoA_shared %>%
  nest(data = - otu) %>%
  mutate(cor_x = map(data, 
                     ~cor.test(.x$count, 
                               .x$PC1, method = "spearman",
                               exact = FALSE) %>% tidy())) %>% #map to apply cor.test to each df
  unnest(cor_x) %>%
  select(otu, estimate, p.value)

cor_y <- PCoA_shared %>%
  nest(data = - otu) %>%
  mutate(cor_y = map(data, 
                     ~cor.test(.x$count, 
                               .x$PC2, method = "spearman",
                               exact = FALSE) %>% tidy())) %>% #map to apply cor.test to each df
  unnest(cor_y) %>%
  select(otu, estimate, p.value)

cor.axes.all <- inner_join(cor_x, cor_y, by = "otu")
#*********************************************************
#*Make vegdist from just cecal samples
cecal_asv_mat <- filtered_asv_counts %>% filter(str_detect(group, 'cecal')) %>%
  pivot_wider(
    names_from = otu,  # Column names will be from the 'asv' column
    values_from = count  # Values will be from the 'count' column
  ) %>%
  column_to_rownames("group")

bc_dist.cecal = vegdist(cecal_asv_mat, method = "bray")
PCOA.cecal = cmdscale(bc_dist.cecal, eig = TRUE, k = 2)
eigenvectors.cecal = PCOA.cecal$eig / sum(PCOA.cecal$eig)

bc_cecal <- PCOA.cecal[["points"]] 
bc_cecal <- as.data.frame(bc_cecal) %>% 
  rename(PC1 = V1, PC2 = V2) 
bc_cecal <- bc_cecal %>% rownames_to_column("group")

beta.cec_join <- inner_join(x = meta.nona, y = bc_cecal, by = 'group')
# Filter asv_table for just cecal samples
filt.cec_asv_table <- filtered_asv_counts %>% filter(str_detect(group, "cecal"))
PCoA.cec_shared <- inner_join(filt.cec_asv_table, bc_cecal, by = "group")

# Using cor.test to see what asvs are driving our cecal clustering
cor_x <- PCoA.cec_shared %>%
  nest(data = - otu) %>%
  mutate(cor_x = map(data, 
                     ~cor.test(.x$count, 
                               .x$PC1, method = "spearman",
                               exact = FALSE) %>% tidy())) %>% #map to apply cor.test to each df
  unnest(cor_x) %>%
  select(otu, estimate, p.value)

cor_y <- PCoA.cec_shared %>%
  nest(data = - otu) %>%
  mutate(cor_y = map(data, 
                     ~cor.test(.x$count, 
                               .x$PC2, method = "spearman",
                               exact = FALSE) %>% tidy())) %>% #map to apply cor.test to each df
  unnest(cor_y) %>%
  select(otu, estimate, p.value)

cor.axes.cecal <- inner_join(cor_x, cor_y, by = "otu")

high_corr_cecal <- cor.axes.cecal %>% 
  filter(abs(estimate.x) > 0.8 | abs(estimate.y) > 0.8) %>%
  filter(p.value.x < 0.005 | p.value.y < 0.005)

# Let's now do the B diversity overtime (aka fecal data)

filtered.fec_asv_mat <- filtered_asv_counts %>%
  filter(!str_detect(group, 'cecal')) %>%
  pivot_wider(
    names_from = otu,  # Column names will be from the 'asv' column
    values_from = count  # Values will be from the 'count' column
  ) %>%
  column_to_rownames("group")


bray.dist.fecal <- avgdist(filtered.fec_asv_mat, sample = 9000, dmethod = "bray")
bd.fecal.mat <- as.matrix(bray.dist.fecal)
bd.fecal.df <- as.data.frame(bd.fecal.mat) %>% rownames_to_column("key")

bray.dist.fecal.to.pre <- bd.fecal.df %>% pivot_longer(cols = -key,
                                                       names_to = "group",
                                                       values_to = "bray_value") %>% distinct() %>%
  mutate(group = str_replace(group, "d57_905_1_2", "d22_905_1")) %>%
  mutate(key = str_replace(key, "d57_905_1_2", "d22_905_1")) %>%
  filter(str_detect(key, 'd')) %>% 
  filter(str_detect(group, 'd')) %>% 
  filter(bray_value > 0) %>% # these should all be comparisons between same sample
  separate(key, into = c("key_day", "key_id"), sep = "_", extra = "merge", remove = FALSE) %>%
  separate(group, into = c("group_day", "group_id"), sep = "_", extra = "merge", remove = FALSE) %>%
  mutate(key_tx = case_when(
    grepl("^901", key_id) | grepl("^902", key_id) | grepl("^907", key_id) |
      grepl("^908", key_id) | grepl("^909", key_id) | grepl("^910", key_id)~ "HCD",  # Condition for "901*" or "902*" etc
    grepl("^903", key_id) | grepl("^904", key_id) | grepl("^911", key_id) |
      grepl("^912", key_id) | grepl("^913", key_id) | grepl("^914", key_id)~ "MCD",  # Condition for "903*" or "904*" etc
    TRUE ~ "xLCD"
  )) %>%
  mutate(group_tx = case_when(
    grepl("^901", group_id) | grepl("^902", group_id) | grepl("^907", group_id) |
      grepl("^908", group_id) | grepl("^909", group_id) | grepl("^910", group_id)~ "HCD",  # Condition for "901*" or "902*" etc
    grepl("^903", group_id) | grepl("^904", group_id) | grepl("^911", group_id) |
      grepl("^912", group_id) | grepl("^913", group_id) | grepl("^914", group_id)~ "MCD",  # Condition for "903*" or "904*" etc
    TRUE ~ "xLCD"
  )) %>%
  filter(key_tx == group_tx)

# Check if 'string' is found anywhere in the data frame
found_str <- apply(test, 2, function(x) any(grepl("d22_905_1", x)))

# Print TRUE if 'apple' is found in any column, otherwise FALSE
if (any(found_str)) {
  print(TRUE)
} else {
  print(FALSE)
}

# Filter by all day comparisons....from each day to d29
bray_dist_d0 <- bray.dist.fecal.to.pre %>% unite(cat_day, key_day, group_day, sep = "-", remove = FALSE) %>%
  filter(str_detect(cat_day, "d0-d29")) %>%
  mutate(combined = paste(pmin(key_id, group_id), pmax(key_id, group_id), sep = "_")) %>%
  distinct(combined, .keep_all = TRUE) %>%
  select(-combined)

bray_dist_d1 <- bray.dist.fecal.to.pre %>% unite(cat_day, key_day, group_day, sep = "-", remove = FALSE) %>%
  filter(str_detect(cat_day, "d1-d29")) %>%
  mutate(combined = paste(pmin(key_id, group_id), pmax(key_id, group_id), sep = "_")) %>%
  distinct(combined, .keep_all = TRUE) %>%
  select(-combined) 

bray_dist_d8 <- bray.dist.fecal.to.pre %>% unite(cat_day, key_day, group_day, sep = "-", remove = FALSE) %>%
  filter(str_detect(cat_day, "d8-d29")) %>%
  mutate(combined = paste(pmin(key_id, group_id), pmax(key_id, group_id), sep = "_")) %>%
  distinct(combined, .keep_all = TRUE) %>%
  select(-combined)

bray_dist_d15 <- bray.dist.fecal.to.pre %>% unite(cat_day, key_day, group_day, sep = "-", remove = FALSE) %>%
  filter(str_detect(cat_day, "d15-d29")) %>%
  mutate(combined = paste(pmin(key_id, group_id), pmax(key_id, group_id), sep = "_")) %>%
  distinct(combined, .keep_all = TRUE) %>%
  select(-combined)

bray_dist_d22 <- bray.dist.fecal.to.pre %>% unite(cat_day, key_day, group_day, sep = "-", remove = FALSE) %>%
  filter(str_detect(cat_day, "d22-d29")) %>%
  mutate(combined = paste(pmin(key_id, group_id), pmax(key_id, group_id), sep = "_")) %>%
  distinct(combined, .keep_all = TRUE) %>%
  select(-combined)

# Post diet change
bray_dist_d29 <- bray.dist.fecal.to.pre %>% unite(cat_day, key_day, group_day, sep = "-", remove = FALSE) %>%
  filter(str_detect(cat_day, "d29-d29")) %>%
  mutate(combined = paste(pmin(key_id, group_id), pmax(key_id, group_id), sep = "_")) %>%
  distinct(combined, .keep_all = TRUE) %>%
  select(-combined)

bray_dist_d36 <- bray.dist.fecal.to.pre %>% unite(cat_day, key_day, group_day, sep = "-", remove = FALSE) %>%
  filter(str_detect(cat_day, "d29-d36")) %>%
  mutate(combined = paste(pmin(key_id, group_id), pmax(key_id, group_id), sep = "_")) %>%
  distinct(combined, .keep_all = TRUE) %>%
  select(-combined)

bray_dist_d43 <- bray.dist.fecal.to.pre %>% unite(cat_day, key_day, group_day, sep = "-", remove = FALSE) %>%
  filter(str_detect(cat_day, "d29-d43")) %>%
  mutate(combined = paste(pmin(key_id, group_id), pmax(key_id, group_id), sep = "_")) %>%
  distinct(combined, .keep_all = TRUE) %>%
  select(-combined)

bray_dist_d50 <- bray.dist.fecal.to.pre %>% unite(cat_day, key_day, group_day, sep = "-", remove = FALSE) %>%
  filter(str_detect(cat_day, "d29-d50")) %>%
  mutate(combined = paste(pmin(key_id, group_id), pmax(key_id, group_id), sep = "_")) %>%
  distinct(combined, .keep_all = TRUE) %>%
  select(-combined)

bray_dist_d57 <- bray.dist.fecal.to.pre %>% unite(cat_day, key_day, group_day, sep = "-", remove = FALSE) %>%
  filter(str_detect(cat_day, "d29-d57")) %>%
  mutate(combined = paste(pmin(key_id, group_id), pmax(key_id, group_id), sep = "_")) %>%
  distinct(combined, .keep_all = TRUE) %>%
  select(-combined)

# Combine all five data frames
combined_bray <- rbind(bray_dist_d0, bray_dist_d1, bray_dist_d8, bray_dist_d15, bray_dist_d22,
                       bray_dist_d29, bray_dist_d36, bray_dist_d43, bray_dist_d50, bray_dist_d57) %>% 
  group_by(cat_day, key_tx) %>%
  mutate(
    mean_bray = mean(bray_value),
    sd_bray = sd(bray_value),
    se_bray = sd(bray_value) / sqrt(n()), 
    median = median(bray_value),
    Q1 = quantile(bray_value, 0.25),
    Q3 = quantile(bray_value, 0.75)
  )
# Make a new day column to plot continuous x-axis 
combined_bray1 <- combined_bray %>%
  mutate(new_day = case_when(
    cat_day == "d0-d29" ~ 0,
    cat_day == "d1-d29" ~ 1,
    cat_day == "d8-d29" ~ 8,
    cat_day == "d15-d29" ~ 15,
    cat_day == "d22-d29" ~ 22,
    cat_day == "d29-d29" ~ 29,
    cat_day == "d29-d36" ~ 36,
    cat_day == "d29-d43" ~ 43,
    cat_day == "d29-d50" ~ 50,
    cat_day == "d29-d57" ~ 57,
    TRUE ~ NA_real_  # Any other cases that don't match, set to NA (optional)
  )) 

# Figure 2E
# Plot your data
ggplot(combined_bray1, aes(x=new_day, y=bray_value, color = key_tx, fill = key_tx)) +
  geom_jitter(width=0.25, alpha=0.4, size = 1, na.rm=TRUE) +
  scale_color_manual(values=c("HCD"="orange", "MCD"="green", "xLCD"="purple")) +
  stat_summary(fun.data=median_hilow, na.rm= TRUE, geom="line", linewidth=1,
               fun.args = list(conf.int=0.50))  +
  geom_ribbon(aes(ymin = Q1, ymax = Q3, group=key_tx), alpha=0.2, color = NA) +
  scale_fill_manual(values=c("HCD"="orange", "MCD"="green", "xLCD"="purple")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,1)) +
  scale_y_break()
scale_x_continuous(breaks = c(0, 1, 7, 14, 21, 29, 35, 42, 50, 57)) +
  geom_vline(xintercept = 29, linetype = "solid", color = "red", size = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(x = "day", y = "Mean distsance to\n pre diet change community") # add \n to add a line break in the label

# Make dist matrices for each comparison we want to do
# Figure 2D
filtered.cecal_asv_mat <- filtered_asv_counts %>%
  filter(str_detect(group, 'cecal')) %>%
  pivot_wider(
    names_from = otu,  # Column names will be from the 'asv' column
    values_from = count  # Values will be from the 'count' column
  ) %>%
  column_to_rownames("group")

# HCD:MCD first
HMcecal.mat <- filtered_asv_counts %>%
  filter(str_detect(group, 'cecal')) %>%
  filter(str_detect(group, "901|907|908|903|911|912")) %>%
  pivot_wider(
    names_from = otu,  # Column names will be from the 'asv' column
    values_from = count  # Values will be from the 'count' column
  ) %>%
  column_to_rownames("group")

bray.dist.HMcecal <- avgdist(HMcecal.mat, sample = 9000, dmethod = "bray")
bd.HMcecal.mat <- as.matrix(bray.dist.HMcecal)
bd.HMcecal.df <- as.data.frame(bd.HMcecal.mat) %>% rownames_to_column("key") %>%
  filter(str_detect(key, "901|907|908")) %>%
  select(matches("key|903|911|912"))

HM_cecal.df <- bd.HMcecal.df %>% 
  pivot_longer(-key, names_to="group", values_to="bray.dist") %>%
  mutate(comparison = "HCD:MCD")

# HCD:LCD second
HLcecal.mat <- filtered_asv_counts %>%
  filter(str_detect(group, 'cecal')) %>%
  filter(str_detect(group, "901|907|908|905|915|916")) %>%
  pivot_wider(
    names_from = otu,  # Column names will be from the 'asv' column
    values_from = count  # Values will be from the 'count' column
  ) %>%
  column_to_rownames("group")

bray.dist.HLcecal <- avgdist(HLcecal.mat, sample = 9000, dmethod = "bray")
bd.HLcecal.mat <- as.matrix(bray.dist.HLcecal)
bd.HLcecal.df <- as.data.frame(bd.HLcecal.mat) %>% rownames_to_column("key") %>%
  filter(str_detect(key, "901|907|908")) %>%
  select(matches("key|905|915|916"))

HL_cecal.df <- bd.HLcecal.df %>% 
  pivot_longer(-key, names_to="group", values_to="bray.dist") %>%
  mutate(comparison = "HCD:LCD")

# MCD:LCD second
MLcecal.mat <- filtered_asv_counts %>%
  filter(str_detect(group, 'cecal')) %>%
  filter(str_detect(group, "903|911|912|905|915|916")) %>%
  pivot_wider(
    names_from = otu,  # Column names will be from the 'asv' column
    values_from = count  # Values will be from the 'count' column
  ) %>%
  column_to_rownames("group")

bray.dist.MLcecal <- avgdist(MLcecal.mat, sample = 9000, dmethod = "bray")
bd.MLcecal.mat <- as.matrix(bray.dist.MLcecal)
bd.MLcecal.df <- as.data.frame(bd.MLcecal.mat) %>% rownames_to_column("key") %>%
  filter(str_detect(key, "903|911|912")) %>%
  select(matches("key|905|915|916"))

ML_cecal.df <- bd.MLcecal.df %>% 
  pivot_longer(-key, names_to="group", values_to="bray.dist") %>%
  mutate(comparison = "MCD:LCD")

comb.cecal.bray <- rbind(HM_cecal.df, HL_cecal.df, ML_cecal.df)

ggplot(comb.cecal.bray, aes(x=comparison, y=bray.dist)) +
  geom_jitter(width=0.25, color="grey", alpha=0.6) +
  stat_summary(fun.data=median_hilow, color="red", size=1,
               fun.args = list(conf.int=0.50)) +
  labs(x=NULL, y="Bray-Curtis distances") +
  scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.2)) +
  scale_y_break(c(0.05, 0.35), space = 0.3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## stats:
library(FSA)
# inputs:
comb.cecal.bray %>% subset(., comparison %in% c("HCD:MCD", "HCD:LCD", "MCD:LCD")) %>%
  dunnTest(bray.dist ~ comparison, data = ., method="bh")






