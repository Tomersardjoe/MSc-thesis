#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(rstatix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript confusion_plots_all.R dup_summary.tsv out_dir")
}

dup_summary_path <- args[1]
out_dir          <- args[2]

dup_summary <- readr::read_tsv(dup_summary_path, show_col_types = FALSE)

# Harmonise tool names
method_colors <- c(
  Coinfinder = "#D95F5F",
  Goldfinder = "#FFD966",
  PanForest  = "#5F9ED1"
)
set_order <- c("Coinfinder","Goldfinder","PanForest")

# Long format for precision/recall/F1
dup_long <- dup_summary %>%
  select(tool, pangenome_id, category, fluidity, openness,
         matches("precision|recall|f1", ignore.case = TRUE)) %>%
  pivot_longer(cols = matches("precision|recall|f1", ignore.case = TRUE),
               names_to = "Metric", values_to = "Value") %>%
  mutate(
    tool = case_when(
      tool == "coinfinder" ~ "Coinfinder",
      tool == "goldfinder" ~ "Goldfinder",
      tool == "panforest"  ~ "PanForest",
      TRUE ~ tool
    ),
    tool = factor(tool, levels = set_order),
    # Make categories and metrics capitalized
    category = recode(category,
                      "closed"   = "Closed",
                      "moderate" = "Moderate",
                      "open"     = "Open",
                      "wildcard" = "Wildcard"),
    Metric = recode(Metric,
                    "precision" = "Precision",
                    "recall"    = "Recall",
                    "f1"        = "F1"),
    Metric = factor(Metric, levels = c("Precision","Recall","F1"))
  )

# Statistical tests

# Kruskal–Wallis: within each tool, do categories differ?
kw_category <- dup_long %>%
  group_by(tool, Metric) %>%
  kruskal_test(Value ~ category)

print(kw_category)

# Kruskal–Wallis: within each category, do tools differ?
kw_tool <- dup_long %>%
  group_by(category, Metric) %>%
  kruskal_test(Value ~ tool)

print(kw_tool)

# Post-hoc pairwise Wilcoxon tests

# Only run post-hoc where omnibus p < 0.05
sig_kw_category <- kw_category %>% filter(p < 0.05)
sig_kw_tool     <- kw_tool %>% filter(p < 0.05)

# Pairwise across categories (per tool & metric)
pairwise_category <- dup_long %>%
  semi_join(sig_kw_category, by = c("tool","Metric")) %>%
  group_by(tool, Metric) %>%
  pairwise_wilcox_test(Value ~ category, p.adjust.method = "BH")

# Pairwise across tools (per category & metric)
pairwise_tool <- dup_long %>%
  semi_join(sig_kw_tool, by = c("category","Metric")) %>%
  group_by(category, Metric) %>%
  pairwise_wilcox_test(Value ~ tool, p.adjust.method = "BH")

# Save post-hoc results
write_tsv(pairwise_category, file.path(out_dir, "pairwise_category_results.tsv"))
write_tsv(pairwise_tool,     file.path(out_dir, "pairwise_tool_results.tsv"))

# Aggregate across pangenomes
dup_long_fluidity <- dup_long %>%
  group_by(tool, Metric, fluidity = round(fluidity, 2)) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

dup_long_openness <- dup_long %>%
  group_by(tool, Metric, openness = round(openness, 2)) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

# Fluidity plot
p_confusion_fluidity <- ggplot(dup_long_fluidity,
                               aes(x = factor(fluidity), y = Value, fill = tool)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, colour = "black") +
  facet_wrap(~Metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(name = "Method", values = method_colors,
                    limits = set_order, breaks = set_order, drop = FALSE) +
  theme_minimal() +
  theme(legend.key = element_rect(fill = "white", colour = NA),
        legend.position = "bottom")

ggsave(file.path(out_dir, "dup_confusion_fluidity.png"),
       plot = p_confusion_fluidity, width = 18, height = 5, dpi = 300, bg = "white")

# Openness plot
p_confusion_open <- ggplot(dup_long_openness,
                           aes(x = factor(openness), y = Value, fill = tool)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, colour = "black") +
  facet_wrap(~Metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(name = "Method", values = method_colors,
                    limits = set_order, breaks = set_order, drop = FALSE) +
  theme_minimal() +
  theme(legend.key = element_rect(fill = "white", colour = NA),
        legend.position = "bottom")

ggsave(file.path(out_dir, "dup_confusion_openness.png"),
       plot = p_confusion_open, width = 18, height = 5, dpi = 300, bg = "white")
       
# Category-based plots
dup_long_cat <- dup_long %>%
  select(tool, pangenome_id, category, Metric, Value)

dup_summary_cat <- dup_long_cat %>%
  group_by(tool, Metric, category) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

# Excluding wildcard/NA
dup_summary_cat_no_wild <- dup_summary_cat %>%
  filter(!is.na(category), category != "Wildcard") %>%
  mutate(category = factor(category, levels = c("Closed","Moderate","Open")))

p_confusion_cat_no_wild <- ggplot(dup_summary_cat_no_wild,
                                  aes(x = category, y = Value, fill = tool)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, colour = "black") +
  facet_wrap(~Metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(name = "Method", values = method_colors,
                    limits = set_order, breaks = set_order, drop = FALSE) +
  labs(x = "Pangenome category", y = "Score",
       title = "Precision, Recall, and F1") +
  theme_minimal() +
  theme(legend.key = element_rect(fill = "white", colour = NA),
        legend.position = "bottom")

ggsave(file.path(out_dir, "dup_confusion_categories_no_wildcard.png"),
       plot = p_confusion_cat_no_wild, width = 12, height = 5, dpi = 300, bg = "white")

# Including wildcard (map NA and 'wildcard' -> 'Wildcard')
dup_summary_cat_wild <- dup_summary_cat %>%
  mutate(category = case_when(
           is.na(category)        ~ "Wildcard",
           category == "Wildcard" ~ "Wildcard",
           TRUE                   ~ category
         ),
         category = factor(category, levels = c("Closed","Moderate","Open","Wildcard")))

p_confusion_cat_wild <- ggplot(dup_summary_cat_wild,
                               aes(x = category, y = Value, fill = tool)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, colour = "black") +
  facet_wrap(~Metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(name = "Method", values = method_colors,
                    limits = set_order, breaks = set_order, drop = FALSE) +
  labs(x = "Pangenome category", y = "Score",
       title = "Precision, Recall, and F1 (including wildcard pangenome)") +
  theme_minimal() +
  theme(legend.key = element_rect(fill = "white", colour = NA),
        legend.position = "bottom")

ggsave(file.path(out_dir, "dup_confusion_categories_wildcard.png"),
       plot = p_confusion_cat_wild, width = 12, height = 5, dpi = 300, bg = "white")

# Define outliers to exclude
outlier_ids <- c(1719, 35790, 783, 83554)

# Filter out outliers
dup_long_no_outliers <- dup_long %>%
  filter(!pangenome_id %in% outlier_ids)

# Aggregate across pangenomes without outliers
dup_long_fluidity_no_outliers <- dup_long_no_outliers %>%
  group_by(tool, Metric, fluidity = round(fluidity, 2)) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

dup_long_openness_no_outliers <- dup_long_no_outliers %>%
  group_by(tool, Metric, openness = round(openness, 2)) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

dup_summary_cat_no_outliers <- dup_long_no_outliers %>%
  group_by(tool, Metric, category) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

# Fluidity plot without outliers
p_confusion_fluidity_no_outliers <- ggplot(dup_long_fluidity_no_outliers,
                               aes(x = factor(fluidity), y = Value, fill = tool)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, colour = "black") +
  facet_wrap(~Metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(name = "Method", values = method_colors,
                    limits = set_order, breaks = set_order, drop = FALSE) +
  theme_minimal() +
  theme(legend.key = element_rect(fill = "white", colour = NA),
        legend.position = "bottom")

ggsave(file.path(out_dir, "dup_confusion_fluidity_no_outliers.png"),
       plot = p_confusion_fluidity_no_outliers, width = 18, height = 5, dpi = 300, bg = "white")

# Openness plot without outliers
p_confusion_open_no_outliers <- ggplot(dup_long_openness_no_outliers,
                           aes(x = factor(openness), y = Value, fill = tool)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, colour = "black") +
  facet_wrap(~Metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(name = "Method", values = method_colors,
                    limits = set_order, breaks = set_order, drop = FALSE) +
  theme_minimal() +
  theme(legend.key = element_rect(fill = "white", colour = NA),
        legend.position = "bottom")

ggsave(file.path(out_dir, "dup_confusion_openness_no_outliers.png"),
       plot = p_confusion_open_no_outliers, width = 18, height = 5, dpi = 300, bg = "white")

# Category-based plots without outliers
dup_summary_cat_no_wild_no_outliers <- dup_summary_cat_no_outliers %>%
  filter(!is.na(category), category != "Wildcard") %>%
  mutate(category = factor(category, levels = c("Closed","Moderate","Open")))

p_confusion_cat_no_wild_no_outliers <- ggplot(dup_summary_cat_no_wild_no_outliers,
                                  aes(x = category, y = Value, fill = tool)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, colour = "black") +
  facet_wrap(~Metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(name = "Method", values = method_colors,
                    limits = set_order, breaks = set_order, drop = FALSE) +
  labs(x = "Pangenome category", y = "Score",
       title = "Precision, Recall, and F1 (excluding outliers)") +
  theme_minimal() +
  theme(legend.key = element_rect(fill = "white", colour = NA),
        legend.position = "bottom")

ggsave(file.path(out_dir, "dup_confusion_categories_no_wildcard_no_outliers.png"),
       plot = p_confusion_cat_no_wild_no_outliers, width = 12, height = 5, dpi = 300, bg = "white")
       
       # Including wildcard (map NA and 'Wildcard' -> 'Wildcard') for no-outliers
dup_summary_cat_wild_no_outliers <- dup_summary_cat_no_outliers %>%
  mutate(category = case_when(
           is.na(category)        ~ "Wildcard",
           category == "Wildcard" ~ "Wildcard",
           TRUE                   ~ category
         ),
         category = factor(category,
                           levels = c("Closed","Moderate","Open","Wildcard")))

p_confusion_cat_wild_no_outliers <- ggplot(dup_summary_cat_wild_no_outliers,
                               aes(x = category, y = Value, fill = tool)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, colour = "black") +
  facet_wrap(~Metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(name = "Method", values = method_colors,
                    limits = set_order, breaks = set_order, drop = FALSE) +
  labs(x = "Pangenome category", y = "Score",
       title = "Precision, Recall, and F1 (including Wildcard, no outliers)") +
  theme_minimal() +
  theme(legend.key = element_rect(fill = "white", colour = NA),
        legend.position = "bottom")

ggsave(file.path(out_dir, "dup_confusion_categories_wildcard_no_outliers.png"),
       plot = p_confusion_cat_wild_no_outliers, width = 12, height = 5, dpi = 300, bg = "white")