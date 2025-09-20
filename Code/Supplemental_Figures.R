# For Reduced carbohydrate complexity and diversity 
# alter gut microbial structure independent of total intake Flores & Seekatz (2025)

# Supplemental Figures

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

#######################################################################
Supplemental Figures
#######################################################################
####################################### 
Supplemental S1
Description: Caloric intake and DEXA parameters
Files needed: food_metadata, dexa_metadata
#######################################
# Figure S1A-D
Lets look at weekly food consumption (g and kcal), & total eaten (kcal)
# kcal/gram from diet breakdown sent by Research Diets for all mice
food <- read_excel("~/path/to/dFR_all_g_intake.xlsx") %>%
  mutate(sex = if_else(cage %in% c(901, 903, 905, 907, 908, 911, 912, 915, 916), "F", "M")) %>%
  mutate(kcal_consumed = case_when(
    tx_group == 'HCD' ~ 3.5 * g_food,
    tx_group == 'MCD' ~ 3.53 * g_food,
    tx_group == 'LCD' ~ 3.55 * g_food
  ),
  #Calculate per mouse consumption (divide by 3, since there are 3 mice per cage)
  g_food_per_mouse = g_food / n_mice,
  kcal_per_mouse = kcal_consumed / n_mice
  ) %>%
  #Calculate average per treatment group, week
  group_by(tx_group, weeks) %>%
  mutate(
    avg_g_food_per_mouse = mean(g_food_per_mouse, na.rm = TRUE),
    sd_g_food_per_mouse = sd(g_food_per_mouse),
    se_g_food_per_mouse = sd_g_food_per_mouse / sqrt(n()),
    avg_kcal_per_mouse = mean(kcal_per_mouse, na.rm = TRUE),
    sd_kcal_per_mouse = sd(kcal_per_mouse),
    se_kcal_per_mouse = sd_kcal_per_mouse / sqrt(n())
  ) %>%
  mutate(
    # Calculate 95% Confidence Interval (CI) based on standard error (1.96 * SE)
    g_ci_lower = avg_g_food_per_mouse - 1.96 * se_g_food_per_mouse,
    g_ci_upper = avg_g_food_per_mouse + 1.96 * se_g_food_per_mouse,
    kcal_ci_lower = avg_kcal_per_mouse - 1.96 * se_kcal_per_mouse,
    kcal_ci_upper = avg_kcal_per_mouse + 1.96 * se_kcal_per_mouse
  ) %>% ungroup()

food_summary <- food %>%
  select(tx_group, weeks, sex, cage, g_food_per_mouse, kcal_per_mouse) %>%
  group_by(tx_group, weeks, sex) %>%
  mutate(
    avg_g_food_per_mouse = mean(g_food_per_mouse, na.rm = TRUE),
    sd_g_food_per_mouse = sd(g_food_per_mouse),
    se_g_food_per_mouse = sd_g_food_per_mouse / sqrt(n()),
    avg_kcal_per_mouse = mean(kcal_per_mouse, na.rm = TRUE),
    sd_kcal_per_mouse = sd(kcal_per_mouse),
    se_kcal_per_mouse = sd_kcal_per_mouse / sqrt(n()),
    g_ci_lower = avg_g_food_per_mouse - 1.96 * se_g_food_per_mouse,
    g_ci_upper = avg_g_food_per_mouse + 1.96 * se_g_food_per_mouse,
    kcal_ci_lower = avg_kcal_per_mouse - 1.96 * se_kcal_per_mouse,
    kcal_ci_upper = avg_kcal_per_mouse + 1.96 * se_kcal_per_mouse
  ) %>%
  ungroup() %>%
  group_by(tx_group, cage) %>%
  # remember this is an avg per mouse per cage since individual mouse food consuption was not tracked
  mutate(total_kcal_per_mouse = sum(kcal_per_mouse, na.rm = TRUE),
         total_g_per_mouse = sum(g_food_per_mouse, na.rm = TRUE)) %>%
  ungroup()

# Plot the kcal consumed per mouse over time by treatment group
p1 <- ggplot(food_summary, aes(x = weeks, y = kcal_per_mouse, color = tx_group)) +
  # Jitter points for individual data
  geom_jitter(aes(color = tx_group), width = 0.1, alpha = 0.6) +
  
  # Compute the mean for each group and plot it as a line
  stat_summary(
    fun = "mean",  # Compute the mean for each group
    geom = "line", # Add a line for mean
    aes(group = tx_group),
    size = 1.2,
    show.legend = FALSE
  ) +
  # Add a black line over the colored treatment group lines
  stat_summary(
    fun = "mean",  # Compute the mean for each group again
    geom = "line", # Add the black line in the center
    aes(group = tx_group),
    color = "black", # Set color to black
    size = 0.5, # Smaller size for the black line
    show.legend = FALSE
  ) +
  # Use pre-calculated confidence intervals for the ribbons
  geom_ribbon(
    aes(ymin = kcal_ci_lower, ymax = kcal_ci_upper, fill = tx_group), alpha = 0.2, color = NA,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("HCD" = "orange", "MCD" = "green", "LCD" = "purple")) +
  scale_fill_manual(values = c("HCD" = "orange", "MCD" = "green", "LCD" = "purple")) + 
  labs(
    x = "Weeks",
    y = "Weekly Consumption / Mouse (kcal)",
    color = "Treatment Group",
    fill = "Treatment Group"
  ) +
  facet_wrap(~sex, ncol = 1) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        panel.grid.major = element_blank(),   # Removes major grid lines
        panel.grid.minor = element_blank(),    # Removes minor grid lines
        legend.title = element_blank(),
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

# Total eaten
p2 <- food_summary %>%
  filter(weeks == 8) %>%
  ggplot(., aes(x = factor(tx_group, levels = c("HCD", "MCD", "LCD")), y = total_kcal_per_mouse, color = tx_group)) +
  geom_jitter() +
  geom_boxplot(fill = NA, show.legend = FALSE) +
  facet_wrap(~sex, ncol =1) +
  scale_color_manual(values = c("HCD" = "orange", "MCD" = "green", "LCD" = "purple")) +
  labs(x = "", y = "Total consumption / mouse (kcal)") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        panel.grid.major = element_blank(),   # Removes major grid lines
        panel.grid.minor = element_blank(),    # Removes minor grid lines
        legend.title = element_blank(),
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

(p1 | p2) + 
  plot_layout(widths = c(4, 0.5), heights = c(4, 0.01))k(),
        strip.placement = "outside")

# Figure S1E-J
# Now to graph the DXA scan data for experiment 1 mice
dxa <- read_excel("~/path/to/DXA_Results.xlsx")
dxa_long <- dxa %>%
  pivot_longer(cols = c(BMC, BMD, Fat_Mass,
                        Lean_BMC, Fat_Percent,
                        Total_Mass),
               names_to = "measure",
               values_to = "value")

ggplot(dxa_long, aes(x = factor(diet, levels = c("HCD", "MCD", "LCD")), y = value)) +
  geom_jitter() +
  geom_boxplot(fill = NA) +
  labs(x = "", y = "") +
  facet_wrap(~ factor(measure,
                      levels = c("BMC", "BMD", "Lean_BMC",
                                 "Fat_Percent", "Fat_Mass", "Total_Mass")),
             scales = "free_y") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        panel.grid.major = element_blank(),   # Removes major grid lines
        panel.grid.minor = element_blank(),    # Removes minor grid lines
        legend.title = element_blank(),
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none") +
  # Add significance brackets using ggsignif
  geom_signif(comparisons = list(c("HCD", "MCD"), c("MCD", "LCD"), c("HCD", "LCD")), 
              map_signif_level = TRUE,  # Automatically map significance levels to the brackets
              test = "kruskal.test"    # You can change this to your statistical test (e.g., "t.test")
              #textsize = 5
              #y_position = c(2.5, 3, 3.5)# Size of the p-value text
  )  # Adjust vertical position of the brackets

####################################### 
Supplemental S2
Description: Transit time and intestinal permeability
Files needed: transit_time_metadata, ip_metadata
#######################################

# Figure S2A
gi_tt <- read_excel("~/path/to/GI_transit_time.xlsx", sheet = 2) %>%
  mutate(hours = minutes/60)

ggplot(gi_tt, aes(x = factor(tx_group, levels = c("HCD", "MCD", "LCD")), y = hours)) +
  geom_jitter() +
  geom_boxplot(fill = NA, outlier.shape = NA) +
  #facet_wrap(~sex) +
  labs(x = "", y = "Hours to first sign of dye") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        panel.grid.major = element_blank(),   # Removes major grid lines
        panel.grid.minor = element_blank(),    # Removes minor grid lines
        legend.title = element_blank(),
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none") +
  annotate("text", 
           x = Inf, y = Inf, 
           label = paste("ANOVA p =", format.pval(summary(anova_gi.tt)[[1]]$"Pr(>F)"[1], digits = 3)),
           hjust = 1.1, vjust = 1.4,
           size = 3)

# Figure S2B-C

ip_metadata <- read_excel("/path/to/20220616_FITC_RhodB.xlsx",
                          sheet = 'metadata',
                          na = c("NA", "na", ""))

# Make standard curves by averaging replicates per dye and concentration
ip_metadata %>%
  filter(type == "standard") %>%
  group_by(dye, concentration_ulmL) %>%
  ggplot(., aes(x = concentration_ulmL, y = normalized_fluorescence)) +
  geom_point() +
  stat_smooth(method = "lm", se = TRUE, color = "blue", formula = y ~ x) +
  facet_wrap(~ dye, scales = "free_y") +  # free_y allows each dye to have its own y-axis scale
  theme_minimal()

ip_metadata %>%
  filter(type == "standard") %>%
  group_by(dye) %>%
  do({
    fit = lm(normalized_fluorescence ~ concentration_ulmL, data = .)
    tidy_fit = broom::tidy(fit)
    glance_fit = broom::glance(fit)
    data.frame(
      intercept = tidy_fit$estimate[1],
      slope = tidy_fit$estimate[2],
      r_squared = glance_fit$r.squared
    )
  })

ip_metadata %>%
  filter(type == "standard") %>%
  group_by(dye) %>%
  do(
    model = lm(normalized_fluorescence ~ concentration_ulmL, data = .)
  )

ip_samples <- read_excel("/path/to/20220616_FITC_RhodB.xlsx",
                         sheet = 'samples',
                         na = c("NA", "na", ""))

ggplot(ip_samples, aes(x = treatment_group, y = predicted_concentration)) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2)) +
  facet_wrap(~ dye, scales = "free_y") +
  labs(x=NULL, 
       y="Serum 4-kDa Fitc Dextran (ug/mL)") +
  scale_x_discrete(limits=c("HCD", "MCD", "LCD")) + 
  theme_bw() +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        strip.text = element_text(color = "black"),
        aspect.ratio = 1)

####################################### 
Supplemental S3
Description: Alpha diversity and firmicutes/bacteroides ratio
Files needed: metadata, final.asv.ASV.subsample.shared, final.asv.ASV.cons.taxonomy
#######################################

p_supp_cecal.alphadiv <- joined_ASV_div_long %>%
  filter(alpha_metric != "invsimpson") %>%  # remove invsimpson
  ggplot(aes(x = tx_group, y = value)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~ alpha_metric, scales = "free_y") +
  stat_pvalue_manual(
    adiv_pairwise_results,
    label = "p.adj.signif",  # uses symbols: *, **, ***
    tip.length = 0.01,
    hide.ns = TRUE
  ) +
  theme_bw() +
  labs(
    x = "",
    y = "Alpha Diversity Value",
    title = ""
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.grid.major = element_blank(),  # remove major grid lines
    panel.grid.minor = element_blank(),   # remove minor grid lines
    panel.border = element_rect(color = "black", size = 1)  # thicker border
  )

####################################### 
Supplemental S4-6
Description: Tax redundancy box plots: class 
Files needed: tax_redundant_data
#######################################

# Now create function to do ^^ at every tax level
# Define the function to plot redundancy at each taxonomic level
tax_redund_boxplot <- function(taxon_level) {
  
  # Read the data from the appropriate sheet (based on taxon level)
  tr.df <- read_excel("~/path/to/dFR_tax_redundant_data.xlsx",
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