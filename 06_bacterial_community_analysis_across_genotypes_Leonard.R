# Rhizosphere and endosphere bacterial community analysis across the soybean genotypes 
# under no fertilizer condition at Leonard, the iron-deficient site over the three years (2021-2023)

# Load required packages
library(phyloseq)
library(vegan)
library(pairwiseAdonis)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(rstatix)
library(dplyr)
library(forcats)
library(openxlsx)
library(writexl)

# Load field phyloseq objects
physeq <- readRDS("Condo7_PS.all_field_data(2021_to_2023)")
ps.rare <- readRDS("ps.rare.all.rds")

# Output directory
out_dir <- "Microbial_community_analysis_no_fertilizer_Leonard_across_genotypes_years"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Shannon diversity analysis of rhizosphere and endosphere bacterial communities between
# tolerant and sensitive genotypes across years at Leonard (2021-2023).

# Subset all years and compartments
IDCY1.all <- subset_samples(ps.rare, Year == "Y2021")
IDCY2.all <- subset_samples(ps.rare, Year == "Y2022")
IDCY3.all <- subset_samples(ps.rare, Year == "Y2023")

Leonard.IDCY1 <- subset_samples(IDCY1.all, Location == "Leonard")
Leonard.IDCY2 <- subset_samples(IDCY2.all, Location == "Leonard")
Leonard.IDCY3 <- subset_samples(IDCY3.all, Location == "Leonard")

Leonard.IDCY1.Rhizo.none <- subset_samples(Leonard.IDCY1, Compartments == "Rhizosphere")
Leonard.IDCY1.Endo.none <- subset_samples(Leonard.IDCY1, Compartments == "Endosphere")

Leonard.IDCY2.Rhizo.none <- subset_samples(Leonard.IDCY2, Compartments == "Rhizosphere" & Fertilizer == "None")
Leonard.IDCY2.Endo.none <- subset_samples(Leonard.IDCY2, Compartments == "Endosphere" & Fertilizer == "None")

Leonard.IDCY3.Rhizo.none <- subset_samples(Leonard.IDCY3, Compartments == "Rhizosphere" & Fertilizer == "None")
Leonard.IDCY3.Endo.none <- subset_samples(Leonard.IDCY3, Compartments == "Endosphere" & Fertilizer == "None")

# Genotype order in original metadata
geno_order <- c("Sensitive", "Resistant")

# Genotype order displayed in figures and tables
geno_display_order <- c("Sensitive", "Tolerant")

# Point colors
geno_point_colors <- c(
  "Sensitive" = "black",
  "Tolerant" = "black"
)

# Shannon diversity across genotypes within plant compartments - Leonard 2021
# Fertilizer = None

# Rhizosphere 2021

ps <- Leonard.IDCY1.Rhizo.none
ps <- subset_samples(ps, Genotype %in% geno_order)
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

sample_data(ps)$Genotype <- fct_relevel(factor(sample_data(ps)$Genotype), geno_order)

# Calculate Shannon diversity
alpha_df <- estimate_richness(ps, measures = "Shannon")
alpha_df$SampleID <- rownames(alpha_df)
meta <- as(sample_data(ps), "data.frame")
meta$SampleID <- rownames(meta)
alpha_df <- left_join(alpha_df, meta, by = "SampleID")

# Display Resistant as Tolerant
alpha_df$Genotype <- fct_recode(factor(alpha_df$Genotype), "Tolerant" = "Resistant")
alpha_df$Genotype <- fct_relevel(alpha_df$Genotype, geno_display_order)

# Shapiro-Wilk normality test
sw5a_results <- do.call(rbind, lapply(geno_display_order, function(g) {
  sub <- alpha_df[alpha_df$Genotype == g, ]
  sw <- shapiro.test(sub$Shannon)
  data.frame(
    Compartment = "Rhizosphere",
    Year = "2021",
    Genotype = g,
    n = nrow(sub),
    W = round(sw$statistic, 4),
    p = round(sw$p.value, 4),
    Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
}))

print(sw5a_results)

# Normality met: one-way ANOVA + Tukey HSD
# One-way ANOVA
aov5a <- aov(Shannon ~ Genotype, data = alpha_df)
print(summary(aov5a))

# Tukey HSD post-hoc test
tukey5a_all <- alpha_df %>%
  tukey_hsd(Shannon ~ Genotype) %>%
  add_significance("p.adj")

tukey5a <- tukey5a_all %>%
  filter(p.adj.signif != "ns")

print(as.data.frame(tukey5a_all))

# Position significance brackets
y_max5a <- max(alpha_df$Shannon, na.rm = TRUE)

tukey5a <- tukey5a %>%
  arrange(group1, group2) %>%
  mutate(y.position = y_max5a + seq(0.25, by = 0.55, length.out = n()))

y_upper5a <- if (nrow(tukey5a) > 0) {
  max(tukey5a$y.position) + 0.5
} else {
  y_max5a + 1
}

# Save statistics
aov5a_export <- data.frame(
  test = "One-way ANOVA",
  group1 = NA_character_,
  group2 = NA_character_,
  statistic = summary(aov5a)[[1]]$`F value`[1],
  p = summary(aov5a)[[1]]$`Pr(>F)`[1],
  p.adj = NA_real_,
  p.adj.signif = NA_character_
)

tukey5a_export <- as.data.frame(tukey5a_all) %>%
  mutate(test = "Tukey HSD", statistic = NA_real_, p = NA_real_) %>%
  select(test, group1, group2, statistic, p, p.adj, p.adj.signif)

write_xlsx(bind_rows(aov5a_export, tukey5a_export),
           file.path(out_dir, "Alpha_diversity_statistics_no_fertilizer_rhizosphere_Leonard_2021.xlsx"))

# Figure
p5a <- ggplot(alpha_df, aes(x = Genotype, y = Shannon, fill = Genotype)) +
  geom_boxplot(alpha = 1, outlier.shape = NA, width = 0.6, colour = "black", size = 1.5) +
  geom_jitter(aes(color = Genotype), width = 0.15, height = 0, size = 10, alpha = 0.8) +
  {if (nrow(tukey5a) > 0) stat_pvalue_manual(tukey5a, label = "p.adj.signif", y.position = "y.position", tip.length = 0.02, bracket.size = 1.2, size = 30, fontface = "bold", color = "black", inherit.aes = FALSE)} +
  scale_x_discrete(limits = geno_display_order, drop = FALSE) +
  scale_y_continuous(breaks = seq(0, 6, by = 1), limits = c(0, y_upper5a), expand = expansion(mult = c(0.02, 0.02))) +
  scale_fill_manual(values = c("Sensitive" = "darkgreen", "Tolerant" = "darkgreen"), limits = geno_display_order, drop = FALSE) +
  scale_color_manual(values = geno_point_colors, limits = geno_display_order, drop = FALSE) +
  labs(x = "", y = "Shannon diversity") +
  facet_wrap(~ "Rhizosphere 2021") +
  theme_cowplot() +
  theme(axis.line = element_line(linewidth = 2), legend.position = "none",
        axis.title = element_text(size = 60, face = "bold"), axis.text.x = element_text(size = 60, angle = 0, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 60, face = "bold"), strip.text = element_text(size = 60, face = "bold"),
        strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 1.5), panel.border = element_rect(fill = NA))

ggsave(filename = file.path(out_dir, "Alpha_diversity_no_fertilizer_rhizosphere_Leonard_2021.png"),
       plot = p5a, device = "png", width = 15, height = 12, units = "in", dpi = 600, bg = "white")



# Endosphere 2021
ps <- Leonard.IDCY1.Endo.none
ps <- subset_samples(ps, Genotype %in% geno_order)
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

sample_data(ps)$Genotype <- fct_relevel(factor(sample_data(ps)$Genotype), geno_order)

# Calculate Shannon diversity
alpha_df <- estimate_richness(ps, measures = "Shannon")
alpha_df$SampleID <- rownames(alpha_df)
meta <- as(sample_data(ps), "data.frame")
meta$SampleID <- rownames(meta)
alpha_df <- left_join(alpha_df, meta, by = "SampleID")

# Display Resistant as Tolerant
alpha_df$Genotype <- fct_recode(factor(alpha_df$Genotype), "Tolerant" = "Resistant")
alpha_df$Genotype <- fct_relevel(alpha_df$Genotype, geno_display_order)

# Shapiro-Wilk normality test
sw5b_results <- do.call(rbind, lapply(geno_display_order, function(g) {
  sub <- alpha_df[alpha_df$Genotype == g, ]
  sw <- shapiro.test(sub$Shannon)
  data.frame(
    Compartment = "Endosphere",
    Year = "2021",
    Genotype = g,
    n = nrow(sub),
    W = round(sw$statistic, 4),
    p = round(sw$p.value, 4),
    Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
}))

print(sw5b_results)

# Normality not met: Kruskal-Wallis + pairwise Wilcoxon with BH correction
# Kruskal-Wallis test
kw5b <- alpha_df %>%
  kruskal_test(Shannon ~ Genotype)

print(as.data.frame(kw5b))

# Pairwise Wilcoxon test with BH correction
pwc5b_all <- alpha_df %>%
  wilcox_test(Shannon ~ Genotype, p.adjust.method = "BH", paired = FALSE) %>%
  add_significance("p.adj")

pwc5b <- pwc5b_all %>%
  filter(p.adj.signif != "ns")

print(as.data.frame(pwc5b_all))

# Position significance brackets
y_max5b <- max(alpha_df$Shannon, na.rm = TRUE)

pwc5b <- pwc5b %>%
  arrange(group1, group2) %>%
  mutate(y.position = y_max5b + seq(0.25, by = 0.55, length.out = n()))

y_upper5b <- if (nrow(pwc5b) > 0) {
  max(pwc5b$y.position) + 0.5
} else {
  y_max5b + 1
}

# Save statistics
kw5b_export <- as.data.frame(kw5b) %>%
  mutate(
    test = "Kruskal-Wallis",
    group1 = NA_character_,
    group2 = NA_character_,
    p.adj = NA_real_,
    p.adj.signif = NA_character_
  ) %>%
  select(test, group1, group2, statistic, p, p.adj, p.adj.signif)

pwc5b_export <- as.data.frame(pwc5b_all) %>%
  mutate(test = "Pairwise Wilcoxon (BH adjusted)") %>%
  select(test, group1, group2, statistic, p, p.adj, p.adj.signif)

write_xlsx(bind_rows(kw5b_export, pwc5b_export),
           file.path(out_dir, "Alpha_diversity_statistics_no_fertilizer_endosphere_Leonard_2021.xlsx"))

# Figure
p5b <- ggplot(alpha_df, aes(x = Genotype, y = Shannon, fill = Genotype)) +
  geom_boxplot(alpha = 1, outlier.shape = NA, width = 0.6, colour = "black", size = 1.5) +
  geom_jitter(aes(color = Genotype), width = 0.15, height = 0, size = 10, alpha = 0.8) +
  {if (nrow(pwc5b) > 0) stat_pvalue_manual(pwc5b, label = "p.adj.signif", y.position = "y.position", tip.length = 0.02, bracket.size = 1.2, size = 30, fontface = "bold", color = "black", inherit.aes = FALSE)} +
  scale_x_discrete(limits = geno_display_order, drop = FALSE) +
  scale_y_continuous(breaks = seq(0, 8, by = 1), limits = c(0, y_upper5b), expand = expansion(mult = c(0.02, 0.02))) +
  scale_fill_manual(values = c("Sensitive" = "darkgreen", "Tolerant" = "darkgreen"), limits = geno_display_order, drop = FALSE) +
  scale_color_manual(values = geno_point_colors, limits = geno_display_order, drop = FALSE) +
  labs(x = "", y = "Shannon diversity") +
  facet_wrap(~ "Endosphere 2021") +
  theme_cowplot() +
  theme(axis.line = element_line(linewidth = 2), legend.position = "none",
        axis.title = element_text(size = 60, face = "bold"), axis.text.x = element_text(size = 60, angle = 0, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 60, face = "bold"), strip.text = element_text(size = 60, face = "bold"),
        strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 1.5), panel.border = element_rect(fill = NA))

ggsave(filename = file.path(out_dir, "Alpha_diversity_no_fertilizer_endosphere_Leonard_2021.png"),
       plot = p5b, device = "png", width = 15, height = 12, units = "in", dpi = 600, bg = "white")

# Shannon diversity across genotypes within plant compartments  - Leonard 2022
# Fertilizer = None

# Rhizosphere 2022
ps <- Leonard.IDCY2.Rhizo.none
ps <- subset_samples(ps, Genotype %in% geno_order)
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

sample_data(ps)$Genotype <- fct_relevel(factor(sample_data(ps)$Genotype), geno_order)

# Calculate Shannon diversity
alpha_df <- estimate_richness(ps, measures = "Shannon")
alpha_df$SampleID <- rownames(alpha_df)
meta <- as(sample_data(ps), "data.frame")
meta$SampleID <- rownames(meta)
alpha_df <- left_join(alpha_df, meta, by = "SampleID")

# Display Resistant as Tolerant
alpha_df$Genotype <- fct_recode(factor(alpha_df$Genotype), "Tolerant" = "Resistant")
alpha_df$Genotype <- fct_relevel(alpha_df$Genotype, geno_display_order)

# Shapiro-Wilk normality test
sw6a_results <- do.call(rbind, lapply(geno_display_order, function(g) {
  sub <- alpha_df[alpha_df$Genotype == g, ]
  sw <- shapiro.test(sub$Shannon)
  data.frame(
    Compartment = "Rhizosphere",
    Year = "2022",
    Genotype = g,
    n = nrow(sub),
    W = round(sw$statistic, 4),
    p = round(sw$p.value, 4),
    Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
}))

print(sw6a_results)

# Normality met: one-way ANOVA + Tukey HSD
# One-way ANOVA
aov6a <- aov(Shannon ~ Genotype, data = alpha_df)
print(summary(aov6a))

# Tukey HSD post-hoc test
tukey6a_all <- alpha_df %>%
  tukey_hsd(Shannon ~ Genotype) %>%
  add_significance("p.adj")

tukey6a <- tukey6a_all %>%
  filter(p.adj.signif != "ns")

print(as.data.frame(tukey6a_all))

# Position significance brackets
y_max6a <- max(alpha_df$Shannon, na.rm = TRUE)

tukey6a <- tukey6a %>%
  arrange(group1, group2) %>%
  mutate(y.position = y_max6a + seq(0.25, by = 0.55, length.out = n()))

y_upper6a <- if (nrow(tukey6a) > 0) {
  max(tukey6a$y.position) + 0.5
} else {
  y_max6a + 1
}

# Save statistics
aov6a_export <- data.frame(
  test = "One-way ANOVA",
  group1 = NA_character_,
  group2 = NA_character_,
  statistic = summary(aov6a)[[1]]$`F value`[1],
  p = summary(aov6a)[[1]]$`Pr(>F)`[1],
  p.adj = NA_real_,
  p.adj.signif = NA_character_
)

tukey6a_export <- as.data.frame(tukey6a_all) %>%
  mutate(test = "Tukey HSD", statistic = NA_real_, p = NA_real_) %>%
  select(test, group1, group2, statistic, p, p.adj, p.adj.signif)

write_xlsx(bind_rows(aov6a_export, tukey6a_export),
           file.path(out_dir, "Alpha_diversity_statistics_no_fertilizer_rhizosphere_Leonard_2022.xlsx"))

# Figure
p6a <- ggplot(alpha_df, aes(x = Genotype, y = Shannon, fill = Genotype)) +
  geom_boxplot(alpha = 1, outlier.shape = NA, width = 0.6, colour = "black", size = 1.5) +
  geom_jitter(aes(color = Genotype), width = 0.15, height = 0, size = 10, alpha = 0.8) +
  {if (nrow(tukey6a) > 0) stat_pvalue_manual(tukey6a, label = "p.adj.signif", y.position = "y.position", tip.length = 0.02, bracket.size = 1.2, size = 30, fontface = "bold", color = "black", inherit.aes = FALSE)} +
  scale_x_discrete(limits = geno_display_order, drop = FALSE) +
  scale_y_continuous(breaks = seq(0, 6, by = 1), limits = c(0, y_upper6a), expand = expansion(mult = c(0.02, 0.02))) +
  scale_fill_manual(values = c("Sensitive" = "#6D712E", "Tolerant" = "#6D712E"), limits = geno_display_order, drop = FALSE) +
  scale_color_manual(values = geno_point_colors, limits = geno_display_order, drop = FALSE) +
  labs(x = "", y = "Shannon diversity") +
  facet_wrap(~ "Rhizosphere 2022") +
  theme_cowplot() +
  theme(axis.line = element_line(linewidth = 2), legend.position = "none",
        axis.title = element_text(size = 60, face = "bold"), axis.text.x = element_text(size = 60, angle = 0, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 60, face = "bold"), strip.text = element_text(size = 60, face = "bold"),
        strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 1.5), panel.border = element_rect(fill = NA))

ggsave(filename = file.path(out_dir, "Alpha_diversity_no_fertilizer_rhizosphere_Leonard_2022.png"),
       plot = p6a, device = "png", width = 15, height = 12, units = "in", dpi = 600, bg = "white")


# Endosphere 2022
ps <- Leonard.IDCY2.Endo.none
ps <- subset_samples(ps, Genotype %in% geno_order)
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

sample_data(ps)$Genotype <- fct_relevel(factor(sample_data(ps)$Genotype), geno_order)

# Calculate Shannon diversity
alpha_df <- estimate_richness(ps, measures = "Shannon")
alpha_df$SampleID <- rownames(alpha_df)
meta <- as(sample_data(ps), "data.frame")
meta$SampleID <- rownames(meta)
alpha_df <- left_join(alpha_df, meta, by = "SampleID")

# Display Resistant as Tolerant
alpha_df$Genotype <- fct_recode(factor(alpha_df$Genotype), "Tolerant" = "Resistant")
alpha_df$Genotype <- fct_relevel(alpha_df$Genotype, geno_display_order)

# Shapiro-Wilk normality test
sw6b_results <- do.call(rbind, lapply(geno_display_order, function(g) {
  sub <- alpha_df[alpha_df$Genotype == g, ]
  sw <- shapiro.test(sub$Shannon)
  data.frame(
    Compartment = "Endosphere",
    Year = "2022",
    Genotype = g,
    n = nrow(sub),
    W = round(sw$statistic, 4),
    p = round(sw$p.value, 4),
    Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
}))

print(sw6b_results)


# Normality met: one-way ANOVA + Tukey HSD
# One-way ANOVA
aov6b <- aov(Shannon ~ Genotype, data = alpha_df)
print(summary(aov6b))

# Tukey HSD post-hoc test
tukey6b_all <- alpha_df %>%
  tukey_hsd(Shannon ~ Genotype) %>%
  add_significance("p.adj")

tukey6b <- tukey6b_all %>%
  filter(p.adj.signif != "ns")

print(as.data.frame(tukey6b_all))

# Position significance brackets
y_max6b <- max(alpha_df$Shannon, na.rm = TRUE)

tukey6b <- tukey6b %>%
  arrange(group1, group2) %>%
  mutate(y.position = y_max6b + seq(0.25, by = 0.55, length.out = n()))

y_upper6b <- if (nrow(tukey6b) > 0) {
  max(tukey6b$y.position) + 0.5
} else {
  y_max6b + 1
}

# Save statistics
aov6b_export <- data.frame(
  test = "One-way ANOVA",
  group1 = NA_character_,
  group2 = NA_character_,
  statistic = summary(aov6b)[[1]]$`F value`[1],
  p = summary(aov6b)[[1]]$`Pr(>F)`[1],
  p.adj = NA_real_,
  p.adj.signif = NA_character_
)

tukey6b_export <- as.data.frame(tukey6b_all) %>%
  mutate(test = "Tukey HSD", statistic = NA_real_, p = NA_real_) %>%
  select(test, group1, group2, statistic, p, p.adj, p.adj.signif)

write_xlsx(bind_rows(aov6b_export, tukey6b_export),
           file.path(out_dir, "Alpha_diversity_statistics_no_fertilizer_endosphere_Leonard_2022.xlsx"))

# Figure
p6b <- ggplot(alpha_df, aes(x = Genotype, y = Shannon, fill = Genotype)) +
  geom_boxplot(alpha = 1, outlier.shape = NA, width = 0.6, colour = "black", size = 1.5) +
  geom_jitter(aes(color = Genotype), width = 0.15, height = 0, size = 10, alpha = 0.8) +
  {if (nrow(tukey6b) > 0) stat_pvalue_manual(tukey6b, label = "p.adj.signif", y.position = "y.position", tip.length = 0.02, bracket.size = 1.2, size = 30, fontface = "bold", color = "black", inherit.aes = FALSE)} +
  scale_x_discrete(limits = geno_display_order, drop = FALSE) +
  scale_y_continuous(breaks = seq(0, 6, by = 1), limits = c(0, y_upper6b), expand = expansion(mult = c(0.02, 0.02))) +
  scale_fill_manual(values = c("Sensitive" = "#6D712E", "Tolerant" = "#6D712E"), limits = geno_display_order, drop = FALSE) +
  scale_color_manual(values = geno_point_colors, limits = geno_display_order, drop = FALSE) +
  labs(x = "", y = "Shannon diversity") +
  facet_wrap(~ "Endosphere 2022") +
  theme_cowplot() +
  theme(axis.line = element_line(linewidth = 2), legend.position = "none",
        axis.title = element_text(size = 60, face = "bold"), axis.text.x = element_text(size = 60, angle = 0, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 60, face = "bold"), strip.text = element_text(size = 60, face = "bold"),
        strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 1.5), panel.border = element_rect(fill = NA))

ggsave(filename = file.path(out_dir, "Alpha_diversity_no_fertilizer_endosphere_Leonard_2022.png"),
       plot = p6b, device = "png", width = 15, height = 12, units = "in", dpi = 600, bg = "white")

# Shannon diversity across genotypes within plant compartments - Leonard 2023
# Fertilizer = None

# Rhizosphere 2023
ps <- Leonard.IDCY3.Rhizo.none
ps <- subset_samples(ps, Genotype %in% geno_order)
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

sample_data(ps)$Genotype <- fct_relevel(factor(sample_data(ps)$Genotype), geno_order)

# Calculate Shannon diversity
alpha_df <- estimate_richness(ps, measures = "Shannon")
alpha_df$SampleID <- rownames(alpha_df)
meta <- as(sample_data(ps), "data.frame")
meta$SampleID <- rownames(meta)
alpha_df <- left_join(alpha_df, meta, by = "SampleID")

# Display Resistant as Tolerant
alpha_df$Genotype <- fct_recode(factor(alpha_df$Genotype), "Tolerant" = "Resistant")
alpha_df$Genotype <- fct_relevel(alpha_df$Genotype, geno_display_order)

# Shapiro-Wilk normality test
sw7a_results <- do.call(rbind, lapply(geno_display_order, function(g) {
  sub <- alpha_df[alpha_df$Genotype == g, ]
  sw <- shapiro.test(sub$Shannon)
  data.frame(
    Compartment = "Rhizosphere",
    Year = "2023",
    Genotype = g,
    n = nrow(sub),
    W = round(sw$statistic, 4),
    p = round(sw$p.value, 4),
    Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
}))

print(sw7a_results)

# Normality met: one-way ANOVA + Tukey HSD
# One-way ANOVA
aov7a <- aov(Shannon ~ Genotype, data = alpha_df)
print(summary(aov7a))

# Tukey HSD post-hoc test
tukey7a_all <- alpha_df %>%
  tukey_hsd(Shannon ~ Genotype) %>%
  add_significance("p.adj")

tukey7a <- tukey7a_all %>%
  filter(p.adj.signif != "ns")

print(as.data.frame(tukey7a_all))

# Position significance brackets
y_max7a <- max(alpha_df$Shannon, na.rm = TRUE)

tukey7a <- tukey7a %>%
  arrange(group1, group2) %>%
  mutate(y.position = y_max7a + seq(0.25, by = 0.55, length.out = n()))

y_upper7a <- if (nrow(tukey7a) > 0) {
  max(tukey7a$y.position) + 0.5
} else {
  y_max7a + 1
}

# Save statistics
aov7a_export <- data.frame(
  test = "One-way ANOVA",
  group1 = NA_character_,
  group2 = NA_character_,
  statistic = summary(aov7a)[[1]]$`F value`[1],
  p = summary(aov7a)[[1]]$`Pr(>F)`[1],
  p.adj = NA_real_,
  p.adj.signif = NA_character_
)

tukey7a_export <- as.data.frame(tukey7a_all) %>%
  mutate(test = "Tukey HSD", statistic = NA_real_, p = NA_real_) %>%
  select(test, group1, group2, statistic, p, p.adj, p.adj.signif)

write_xlsx(bind_rows(aov7a_export, tukey7a_export),
           file.path(out_dir, "Alpha_diversity_statistics_no_fertilizer_rhizosphere_Leonard_2023.xlsx"))

# Figure
p7a <- ggplot(alpha_df, aes(x = Genotype, y = Shannon, fill = Genotype)) +
  geom_boxplot(alpha = 1, outlier.shape = NA, width = 0.6, colour = "black", size = 1.5) +
  geom_jitter(aes(color = Genotype), width = 0.15, height = 0, size = 10, alpha = 0.8) +
  {if (nrow(tukey7a) > 0) stat_pvalue_manual(tukey7a, label = "p.adj.signif", y.position = "y.position", tip.length = 0.02, bracket.size = 1.2, size = 30, fontface = "bold", color = "black", inherit.aes = FALSE)} +
  scale_x_discrete(limits = geno_display_order, drop = FALSE) +
  scale_y_continuous(breaks = seq(0, 6, by = 1), limits = c(0, y_upper7a), expand = expansion(mult = c(0.02, 0.02))) +
  scale_fill_manual(values = c("Sensitive" = "#7EBD01", "Tolerant" = "#7EBD01"), limits = geno_display_order, drop = FALSE) +
  scale_color_manual(values = geno_point_colors, limits = geno_display_order, drop = FALSE) +
  labs(x = "", y = "Shannon diversity") +
  facet_wrap(~ "Rhizosphere 2023") +
  theme_cowplot() +
  theme(axis.line = element_line(linewidth = 2), legend.position = "none",
        axis.title = element_text(size = 60, face = "bold"), axis.text.x = element_text(size = 60, angle = 0, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 60, face = "bold"), strip.text = element_text(size = 60, face = "bold"),
        strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 1.5), panel.border = element_rect(fill = NA))

ggsave(filename = file.path(out_dir, "Alpha_diversity_no_fertilizer_rhizosphere_Leonard_2023.png"),
       plot = p7a, device = "png", width = 15, height = 12, units = "in", dpi = 600, bg = "white")



# Endosphere 2023
ps <- Leonard.IDCY3.Endo.none
ps <- subset_samples(ps, Genotype %in% geno_order)
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

sample_data(ps)$Genotype <- fct_relevel(factor(sample_data(ps)$Genotype), geno_order)

# Calculate Shannon diversity
alpha_df <- estimate_richness(ps, measures = "Shannon")
alpha_df$SampleID <- rownames(alpha_df)
meta <- as(sample_data(ps), "data.frame")
meta$SampleID <- rownames(meta)
alpha_df <- left_join(alpha_df, meta, by = "SampleID")

# Display Resistant as Tolerant
alpha_df$Genotype <- fct_recode(factor(alpha_df$Genotype), "Tolerant" = "Resistant")
alpha_df$Genotype <- fct_relevel(alpha_df$Genotype, geno_display_order)

# Shapiro-Wilk normality test
sw7b_results <- do.call(rbind, lapply(geno_display_order, function(g) {
  sub <- alpha_df[alpha_df$Genotype == g, ]
  sw <- shapiro.test(sub$Shannon)
  data.frame(
    Compartment = "Endosphere",
    Year = "2023",
    Genotype = g,
    n = nrow(sub),
    W = round(sw$statistic, 4),
    p = round(sw$p.value, 4),
    Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
}))

print(sw7b_results)


# Normality met: one-way ANOVA + Tukey HSD
# One-way ANOVA
aov7b <- aov(Shannon ~ Genotype, data = alpha_df)
print(summary(aov7b))

# Tukey HSD post-hoc test
tukey7b_all <- alpha_df %>%
  tukey_hsd(Shannon ~ Genotype) %>%
  add_significance("p.adj")

tukey7b <- tukey7b_all %>%
  filter(p.adj.signif != "ns")

print(as.data.frame(tukey7b_all))

# Position significance brackets
y_max7b <- max(alpha_df$Shannon, na.rm = TRUE)

tukey7b <- tukey7b %>%
  arrange(group1, group2) %>%
  mutate(y.position = y_max7b + seq(0.25, by = 0.55, length.out = n()))

y_upper7b <- if (nrow(tukey7b) > 0) {
  max(tukey7b$y.position) + 0.5
} else {
  y_max7b + 1
}

# Save statistics
aov7b_export <- data.frame(
  test = "One-way ANOVA",
  group1 = NA_character_,
  group2 = NA_character_,
  statistic = summary(aov7b)[[1]]$`F value`[1],
  p = summary(aov7b)[[1]]$`Pr(>F)`[1],
  p.adj = NA_real_,
  p.adj.signif = NA_character_
)

tukey7b_export <- as.data.frame(tukey7b_all) %>%
  mutate(test = "Tukey HSD", statistic = NA_real_, p = NA_real_) %>%
  select(test, group1, group2, statistic, p, p.adj, p.adj.signif)

write_xlsx(bind_rows(aov7b_export, tukey7b_export),
           file.path(out_dir, "Alpha_diversity_statistics_no_fertilizer_endosphere_Leonard_2023.xlsx"))

# Figure
p7b <- ggplot(alpha_df, aes(x = Genotype, y = Shannon, fill = Genotype)) +
  geom_boxplot(alpha = 1, outlier.shape = NA, width = 0.6, colour = "black", size = 1.5) +
  geom_jitter(aes(color = Genotype), width = 0.15, height = 0, size = 10, alpha = 0.8) +
  {if (nrow(tukey7b) > 0) stat_pvalue_manual(tukey7b, label = "p.adj.signif", y.position = "y.position", tip.length = 0.02, bracket.size = 1.2, size = 30, fontface = "bold", color = "black", inherit.aes = FALSE)} +
  scale_x_discrete(limits = geno_display_order, drop = FALSE) +
  scale_y_continuous(breaks = seq(0, 6, by = 1), limits = c(0, y_upper7b), expand = expansion(mult = c(0.02, 0.02))) +
  scale_fill_manual(values = c("Sensitive" = "#7EBD01", "Tolerant" = "#7EBD01"), limits = geno_display_order, drop = FALSE) +
  scale_color_manual(values = geno_point_colors, limits = geno_display_order, drop = FALSE) +
  labs(x = "", y = "Shannon diversity") +
  facet_wrap(~ "Endosphere 2023") +
  theme_cowplot() +
  theme(axis.line = element_line(linewidth = 2), legend.position = "none",
        axis.title = element_text(size = 60, face = "bold"), axis.text.x = element_text(size = 60, angle = 0, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 60, face = "bold"), strip.text = element_text(size = 60, face = "bold"),
        strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 1.5), panel.border = element_rect(fill = NA))

ggsave(filename = file.path(out_dir, "Alpha_diversity_no_fertilizer_endosphere_Leonard_2023.png"),
       plot = p7b, device = "png", width = 15, height = 12, units = "in", dpi = 600, bg = "white")






# Non-metric Multidimensional Scaling or NMDS analysis using Bray-Curtis distance of rhizosphere bacterial communities 
#across genotypes at Leonard over the years (2021-2023)


# Set seed for reproducibility
set.seed(777)

# Subset Leonard rhizosphere no-fertilizer samples across the years
Leonard <- subset_samples(ps.rare, Location %in% c("Leonard"))
Leonard.Rhizo <- subset_samples(Leonard, Compartments %in% c("Rhizosphere"))
Leonard.Rhizo.none <- subset_samples(Leonard.Rhizo, Fertilizer %in% c("None"))

# Convert to relative abundance
ps.prop.Leonard.Rhizo.none <- transform_sample_counts(Leonard.Rhizo.none, function(otu) otu / sum(otu))


# Bray-Curtis NMDS
ord.nmds.bray.Leonard.Rhizo.none <- ordinate(
  ps.prop.Leonard.Rhizo.none,
  method = "NMDS",
  distance = "bray"
)

# Year colors
custom_colors <- c(
  "Y2021" = "darkgreen",
  "Y2022" = "#6D712E",
  "Y2023" = "#7EBD01"
)

# Plot settings
POINT_SIZE <- 45
POINT_STROKE <- 3.5
JITTER_W <- 0.03
JITTER_H <- 0.03
TITLE_SIZE <- 90
AXIS_TITLE <- 90
AXIS_TEXT <- 90
FIG_W <- 35
FIG_H <- 20
FIG_DPI <- 600

# NMDS coordinates
ord_df <- plot_ordination(
  ps.prop.Leonard.Rhizo.none,
  ord.nmds.bray.Leonard.Rhizo.none,
  color = "Year",
  shape = "Genotype",
  justDF = TRUE
)

# NMDS plot
Plot.Leonard.Rhizo.none <- ggplot(
  ord_df,
  aes(x = NMDS1, y = NMDS2, color = Year, shape = Genotype)
) +
  theme_cowplot() +
  geom_point(
    size = POINT_SIZE,
    stroke = POINT_STROKE,
    position = position_jitter(width = JITTER_W, height = JITTER_H)
  ) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = c(15, 19, 17)) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme(
    legend.position = "right",
    plot.title = element_text(size = TITLE_SIZE, face = "bold", hjust = 0.5),
    axis.title = element_text(size = AXIS_TITLE, face = "bold"),
    axis.text = element_text(size = AXIS_TEXT, face = "bold"),
    axis.line.x = element_line(linewidth = 2, colour = "black"),
    axis.line.y = element_line(linewidth = 2, colour = "black")
  )

ggsave(
  filename = file.path(out_dir, "NMDS_Leonard_rhizosphere_across_genotype_years_no_fertilizer.png"),
  plot = Plot.Leonard.Rhizo.none,
  width = FIG_W,
  height = FIG_H,
  units = "in",
  dpi = FIG_DPI,
  bg = "white"
)



# PERMANOVA and pairwise PERMANOVA comparisons of rhizosphere bacterial communities across genotypes at Leonard over the years

years <- c("Y2021", "Y2022", "Y2023")
dispersion_list <- list()
permanova_list <- list()
pairwise_list <- list()

for (yr in years) {
  
ps_yr <- subset_samples(Leonard.Rhizo.none, Year %in% yr)
  
bray_yr <- phyloseq::distance(ps_yr, method = "bray")
meta_yr <- data.frame(sample_data(ps_yr))
  

# Check dispersion across genotypes
set.seed(777)
disp_yr <- betadisper(bray_yr, meta_yr$Genotype)
dispersion_list[[yr]] <- anova(disp_yr)
  
# PERMANOVA across genotypes
set.seed(777)
permanova_list[[yr]] <- adonis2(bray_yr ~ Genotype, data = meta_yr, permutations = 999)
  
# Pairwise PERMANOVA across genotypes
set.seed(777)
pairwise_list[[yr]] <- pairwise.adonis2(bray_yr ~ Genotype, data = meta_yr, permutations = 999, p.adjust.m = "bonferroni")
}

# Save PERMANOVA and pairwise PERMANOVA results
pairwise_df <- do.call(rbind, lapply(names(pairwise_list), function(nm) {
  x <- pairwise_list[[nm]]
  do.call(rbind, lapply(names(x), function(pair) {
    if (!is.data.frame(x[[pair]])) return(NULL)
    cbind(Year = nm, Comparison = pair, x[[pair]])
  }))
}))

results <- c(
  list(Genotype_pairwise_PERMANOVA = pairwise_df),
  permanova_list
)

write.xlsx(
  results,
  file.path(out_dir, "PERMANOVA_and_pairwise_PERMANOVA_Leonard_rhizosphere_across_genotype_no_fertilizer.xlsx"),
  overwrite = TRUE
)


# NMDS analysis using Bray-Curtis distance of endosphere bacterial communities 
#across genotypes at Leonard over the  years (2021-2023).

# Set seed for reproducibility
set.seed(777)

# Subset Leonard endosphere no-fertilizer samples across the years
Leonard <- subset_samples(ps.rare, Location %in% c("Leonard"))
Leonard.Endo <- subset_samples(Leonard, Compartments %in% c("Endosphere"))
Leonard.Endo.None <- subset_samples(Leonard.Endo, Fertilizer %in% c("None"))

# Convert to relative abundance
ps.prop.Leonard.Endo.none <- transform_sample_counts(Leonard.Endo.None, function(otu) otu / sum(otu))

# Bray-Curtis NMDS
ord.nmds.bray.Leonard.Endo.none <- ordinate(
  ps.prop.Leonard.Endo.none,
  method = "NMDS",
  distance = "bray"
)

# Year colors
custom_colors <- c(
  "Y2021" = "darkgreen",
  "Y2022" = "#6D712E",
  "Y2023" = "#7EBD01"
)

# Plot settings
POINT_SIZE <- 45
POINT_STROKE <- 3.5
JITTER_W <- 0.03
JITTER_H <- 0.03
TITLE_SIZE <- 90
AXIS_TITLE <- 90
AXIS_TEXT <- 90
FIG_W <- 35
FIG_H <- 20
FIG_DPI <- 600

# NMDS coordinates
ord_df <- plot_ordination(
  ps.prop.Leonard.Endo.none,
  ord.nmds.bray.Leonard.Endo.none,
  color = "Year",
  shape = "Genotype",
  justDF = TRUE
)

# NMDS plot
Plot.Leonard.Endo.none <- ggplot(
  ord_df,
  aes(x = NMDS1, y = NMDS2, color = Year, shape = Genotype)
) +
  theme_cowplot() +
  geom_point(
    size = POINT_SIZE,
    stroke = POINT_STROKE,
    position = position_jitter(width = JITTER_W, height = JITTER_H)
  ) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = c(15, 19)) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme(
    legend.position = "right",
    plot.title = element_text(size = TITLE_SIZE, face = "bold", hjust = 0.5),
    axis.title = element_text(size = AXIS_TITLE, face = "bold"),
    axis.text = element_text(size = AXIS_TEXT, face = "bold"),
    axis.line.x = element_line(linewidth = 2, colour = "black"),
    axis.line.y = element_line(linewidth = 2, colour = "black")
  )

ggsave(
  filename = file.path(out_dir, "NMDS_Leonard_endosphere_across_genotype_no_fertilizer.png"),
  plot = Plot.Leonard.Endo.none,
  width = FIG_W,
  height = FIG_H,
  units = "in",
  dpi = FIG_DPI,
  bg = "white"
)



# PERMANOVA and pairwise PERMANOVA comparisons of endosphere bacterial communities 
#across genotypes at Leonard over the years

years <- c("Y2021", "Y2022", "Y2023")
dispersion_list <- list()
permanova_list <- list()
pairwise_list <- list()

for (yr in years) {
  
ps_yr <- subset_samples(Leonard.Endo.None, Year %in% yr)
  
bray_yr <- phyloseq::distance(ps_yr, method = "bray")
meta_yr <- data.frame(sample_data(ps_yr))
  
# Check dispersion across genotypes
set.seed(777)
disp_yr <- betadisper(bray_yr, meta_yr$Genotype)
dispersion_list[[yr]] <- anova(disp_yr)
  
# PERMANOVA across genotypes
set.seed(777)
permanova_list[[yr]] <- adonis2(bray_yr ~ Genotype, data = meta_yr, permutations = 999)
  
# Pairwise PERMANOVA across genotypes
set.seed(777)
pairwise_list[[yr]] <- pairwise.adonis2(bray_yr ~ Genotype, data = meta_yr, permutations = 999, p.adjust.m = "bonferroni")
}

# Save PERMANOVA and pairwise PERMANOVA results
pairwise_df <- do.call(rbind, lapply(names(pairwise_list), function(nm) {
  x <- pairwise_list[[nm]]
  do.call(rbind, lapply(names(x), function(pair) {
    if (!is.data.frame(x[[pair]])) return(NULL)
    cbind(Year = nm, Comparison = pair, x[[pair]])
  }))
}))

results <- c(
  list(Genotype_pairwise_PERMANOVA = pairwise_df),
  permanova_list
)

write.xlsx(
  results,
  file.path(out_dir, "PERMANOVA_and_pairwise_PERMANOVA_Leonard_endosphere_across_genotype_no_fertilizer.xlsx"),
  overwrite = TRUE
)
