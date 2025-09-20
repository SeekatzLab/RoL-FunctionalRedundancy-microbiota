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
Figure 1
Taxonomic Redundancy:
Files needed: weights_metadata, food_metadata
############################################################################################

# Read in data
wgts <- read_delim("~/path/to/dFR_merged_weights.txt")

# select only d0 or d29 weights
d0 <- wgts[wgts$day=="0", c("specimen", "weights")]
d29 <- wgts[wgts$day=="29", c("specimen", "weights")]

# rename 'weights' column:
colnames(d0)[2] <- "d0.weight"
colnames(d29)[2] <- "d29.weight"

# merge with reg file:
wgts2 <- merge(wgts, d0, by.x="specimen", by.y="specimen") %>%
  merge(., d29, by ="specimen")

# new percent:
wgts2$percent_d0 <- (wgts2$weights / wgts2$d0.weight )*100
wgts2$percent_d29 <- (wgts2$weights / wgts2$d29.weight )*100

# look at new data:
wgts2[, c("specimen", "percent_weight_after_diet_change", "percent_d29")]	#this is calculated for diet change day
wgts2[, c("specimen", "percent_weights", "percent_d0")]						#this is calculated for start date

# read out data:
wgts2[order(wgts2$specimen, wgts2$day), c("specimen", "day", "percent_weight_after_diet_change", "percent_d29", "percent_weights", "percent_d0")]

### let's graph:
# renaming a few things to keep it the same:
colnames(wgts2)[4] <- "group"
levels(as.factor(wgts2$group))
#[1] "HCD" "LCD" "MCD"
# reorder to your liking:
wgts2$group <- factor(wgts2$group, levels=c('HCD', 'MCD', 'LCD'))
dcol<-c("orange", "green", "purple")

data<-wgts2
data$group<-as.factor(data$group)
data <- data[order(data$group, data$day),]

#*************************************************************************************************************************
Figure 1C Weight change
#*************************************************************************************************************************

# Using ggplot and tidyverse
longiPlot_ggplot <- function(n) {
  # Grouping and summarizing the data using dplyr
  means <- data %>%
    group_by(group, day) %>%
    summarise(mean_value = mean(.data[[n]], na.rm = TRUE)) %>%
    ungroup()
  
  # Creating jittered plot with ggplot2
  ggplot(data, aes(x = day, y = .data[[n]], color = group)) +
    geom_rect(aes(xmin = 0, xmax = 29, ymin = -Inf, ymax = Inf), 
              fill = "antiquewhite", alpha = 0.2, show.legend = FALSE,
              color = NA) +
    geom_rect(aes(xmin = 29, xmax = 57, ymin = -Inf, ymax = Inf), 
              fill = "azure", alpha = 0.2, show.legend = FALSE,
              color = NA) +
    geom_vline(xintercept = 29, color = "red", alpha = 0.5, size = 1) +
    geom_line(data = means, aes(x = day, y = mean_value, group = group), size = 2, show.legend = FALSE) +
    geom_line(data = means, aes(x = day, y = mean_value, group = group), color = "black", size = 0.5, show.legend = FALSE) +
    geom_jitter(alpha = 0.5, size = 1) +
    scale_color_manual(values = dcol) +
    scale_x_continuous(breaks = unique(data$day)) +
    labs(x = "Day", y = names(data)[n]) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 10, hjust = 1),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12),
          axis.ticks.x = element_line(size = 0.5),
          legend.position = c(0.94, 0.01),
          legend.justification = c("right", "bottom"),
          legend.margin = margin(t = 0, r = 6, b = 0, l = 0),
          legend.title = element_blank(),
          legend.box.background = element_rect(color = "black", size = 0.5),
          panel.grid.major = element_blank(),   # Removes major grid lines
          panel.grid.minor = element_blank()    # Removes minor grid lines
    ) +
    coord_fixed(ratio = 0.3) # fix aspect ratio
}
p1 <- longiPlot_ggplot("percent_d29") +
  ylab("percent_d29")

# Stats
# Fit a mixed model with treatment, time, and their interaction as fixed effects, and subject (mouse) as a random effect
model <- lmer(weights ~ 
                # fixed effects
                group * day + 
                # random effects
                (1|specimen), data = data
)
summary(model)

#*************************************************************************************************************************
Figure 1D Weight gain
#*************************************************************************************************************************
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

# Total eaten
food_summary %>%
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


