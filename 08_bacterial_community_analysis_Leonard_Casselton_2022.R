# Soybean bacterial community analysis across different fertilizer treatments 
# at Leonard and Casselton in 2022

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
library(DESeq2)
library(writexl)
library(ggh4x)
library(grid)
library(tidyr)

# Load field phyloseq objects
physeq <- readRDS("Condo7_PS.all_field_data(2021_to_2023)")
ps.rare <- readRDS("ps.rare.all.rds")


# Output directory
out_dir <- "bacterial_community_analysis_Leonard_Casselton_2022"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# Fertilizer order and colors
fert_order <- c("No Fertilizer", "Medium Fertilizer", "High Fertilizer")

fert_box_colors <- c(
  "No Fertilizer" = "#8B1C62",
  "Medium Fertilizer" = "#00008B",
  "High Fertilizer" = "#009E73"
)

fert_point_colors <- c(
  "No Fertilizer" = "black",
  "Medium Fertilizer" = "black",
  "High Fertilizer" = "black"
)



# Shannon diversity analysis across fertilizer treatments within each plant compartment at Leonard in 2022
# Subset by plant compartments across the fertilizer levels

IDCY2.all <- subset_samples(ps.rare, Year == "Y2022")
Leonard.IDCY2 <- subset_samples(IDCY2.all, Location == "Leonard")
Leonard.IDCY2.Bulk.fert  <- subset_samples(Leonard.IDCY2, Compartments == "Bulk")
Leonard.IDCY2.Rhizo.fert <- subset_samples(Leonard.IDCY2, Compartments == "Rhizosphere")
Leonard.IDCY2.Endo.fert  <- subset_samples(Leonard.IDCY2, Compartments == "Endosphere")


# Bulk soil
ps <- Leonard.IDCY2.Bulk.fert
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

f <- sample_data(ps)$Fertilizer |> as.character() |> trimws()
f <- dplyr::recode(f,
                   "None"   = "No Fertilizer",
                   "Medium" = "Medium Fertilizer",
                   "High"   = "High Fertilizer"
)
sample_data(ps)$Fertilizer <- fct_relevel(factor(f), fert_order)


# Calculate Shannon diversity
alpha_df <- estimate_richness(ps, measures = "Shannon")
alpha_df$SampleID <- rownames(alpha_df)
meta <- as(sample_data(ps), "data.frame")
meta$SampleID <- rownames(meta)
alpha_df <- left_join(alpha_df, meta, by = "SampleID")
alpha_df$Fertilizer <- fct_relevel(
  factor(trimws(as.character(alpha_df$Fertilizer))),
  fert_order
)


# Shapiro-Wilk normality test
sw_bulk <- do.call(rbind, lapply(fert_order, function(g) {
  sub <- alpha_df[alpha_df$Fertilizer == g, ]
  sw  <- shapiro.test(sub$Shannon)
  data.frame(
    Compartment = "Bulk",
    Fertilizer  = g,
    n           = nrow(sub),
    W           = round(sw$statistic, 4),
    p           = round(sw$p.value, 4),
    Normal      = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
}))

print(sw_bulk)

# Normality met: one-way ANOVA + Tukey HSD
# One-way ANOVA
aov_bulk <- aov(Shannon ~ Fertilizer, data = alpha_df)

print(summary(aov_bulk))


# Tukey HSD post-hoc test
tukey_bulk_all <- alpha_df %>%
  tukey_hsd(Shannon ~ Fertilizer) %>%
  add_significance("p.adj")

tukey_bulk <- tukey_bulk_all %>%
  filter(p.adj.signif != "ns")

print(as.data.frame(tukey_bulk_all))


# Position significance brackets
y_max_bulk <- max(alpha_df$Shannon, na.rm = TRUE)

tukey_bulk <- tukey_bulk %>%
  arrange(group1, group2) %>%
  mutate(
    y.position = y_max_bulk + seq(0.25, by = 0.65, length.out = n())
  )

y_upper_bulk <- if (nrow(tukey_bulk) > 0) {
  max(tukey_bulk$y.position) + 0.5
} else {
  y_max_bulk + 1
}


# Save statistics
aov_bulk_export <- data.frame(
  test         = "One-way ANOVA",
  group1       = NA_character_,
  group2       = NA_character_,
  statistic    = summary(aov_bulk)[[1]]$`F value`[1],
  p            = summary(aov_bulk)[[1]]$`Pr(>F)`[1],
  p.adj        = NA_real_,
  p.adj.signif = NA_character_
)

tukey_bulk_export <- as.data.frame(tukey_bulk_all) %>%
  mutate(test = "Tukey HSD", statistic = NA_real_, p = NA_real_) %>%
  select(test, group1, group2, statistic, p, p.adj, p.adj.signif)

write_xlsx(
  bind_rows(aov_bulk_export, tukey_bulk_export),
  file.path(out_dir, "Alpha_diversity_statistics_bulk_fertilizerlevels_Leonard_2022.xlsx")
)


# Figure
p_bulk <- ggplot(
  alpha_df,
  aes(x = Fertilizer, y = Shannon, fill = Fertilizer)
) +
  geom_boxplot(
    alpha = 1, outlier.shape = NA,
    width = 0.6, colour = "black", size = 1.5
  ) +
  geom_jitter(
    aes(color = Fertilizer),
    width = 0.15, height = 0, size = 10, alpha = 0.8
  ) +
  {
    if (nrow(tukey_bulk) > 0)
      stat_pvalue_manual(
        tukey_bulk,
        label        = "p.adj.signif",
        y.position   = "y.position",
        tip.length   = 0.02,
        bracket.size = 1.2,
        size         = 30,
        fontface     = "bold",
        color        = "black",
        inherit.aes  = FALSE
      )
  } +
  scale_x_discrete(limits = fert_order, drop = FALSE) +
  scale_y_continuous(
    breaks = seq(0, 10, by = 2),
    limits = c(0, y_upper_bulk),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_fill_manual(values  = fert_box_colors,   limits = fert_order, drop = FALSE) +
  scale_color_manual(values = fert_point_colors, limits = fert_order, drop = FALSE) +
  labs(x = "", y = "Shannon diversity") +
  facet_wrap(~ "Bulk soil") +
  theme_cowplot() +
  theme(
    axis.line        = element_line(linewidth = 2),
    legend.position  = "right",
    axis.title       = element_text(size = 60, face = "bold"),
    axis.text.x      = element_text(size = 60, angle = 45, hjust = 1, face = "bold"),
    axis.text.y      = element_text(size = 60, face = "bold"),
    strip.text       = element_text(size = 60, face = "bold"),
    strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 1.5),
    panel.border     = element_rect(fill = NA)
  )

ggsave(filename = file.path(out_dir, "Alpha_diversity_fertilizerlevels_bulk_Leonard_2022.png"),
  plot = p_bulk, device = "png", width = 15, height = 15, units = "in", dpi = 600, bg = "white"
)


# Rhizosphere
ps <- Leonard.IDCY2.Rhizo.fert
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

f <- sample_data(ps)$Fertilizer |> as.character() |> trimws()
f <- dplyr::recode(f,
                   "None"   = "No Fertilizer",
                   "Medium" = "Medium Fertilizer",
                   "High"   = "High Fertilizer"
)
sample_data(ps)$Fertilizer <- fct_relevel(factor(f), fert_order)


# Calculate Shannon diversity
alpha_df <- estimate_richness(ps, measures = "Shannon")
alpha_df$SampleID <- rownames(alpha_df)
meta <- as(sample_data(ps), "data.frame")
meta$SampleID <- rownames(meta)
alpha_df <- left_join(alpha_df, meta, by = "SampleID")
alpha_df$Fertilizer <- fct_relevel(
  factor(trimws(as.character(alpha_df$Fertilizer))),
  fert_order
)


# Shapiro-Wilk normality test
sw_rhiz <- do.call(rbind, lapply(fert_order, function(g) {
  sub <- alpha_df[alpha_df$Fertilizer == g, ]
  sw  <- shapiro.test(sub$Shannon)
  data.frame(
    Compartment = "Rhizosphere",
    Fertilizer  = g,
    n           = nrow(sub),
    W           = round(sw$statistic, 4),
    p           = round(sw$p.value, 4),
    Normal      = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
}))

print(sw_rhiz)

# Normality met: one-way ANOVA + Tukey HSD
# One-way ANOVA
aov_rhiz <- aov(Shannon ~ Fertilizer, data = alpha_df)

print(summary(aov_rhiz))


# Tukey HSD post-hoc test
tukey_rhiz_all <- alpha_df %>%
  tukey_hsd(Shannon ~ Fertilizer) %>%
  add_significance("p.adj")

tukey_rhiz <- tukey_rhiz_all %>%
  filter(p.adj.signif != "ns") %>%
  mutate(p.adj.signif = ifelse(p.adj.signif == "****", "***", p.adj.signif))




# Position significance brackets
y_max_rhiz <- max(alpha_df$Shannon, na.rm = TRUE)

tukey_rhiz <- tukey_rhiz %>%
  arrange(group1, group2) %>%
  mutate(
    y.position = y_max_rhiz + seq(0.25, by = 0.65, length.out = n())
  )

y_upper_rhiz <- if (nrow(tukey_rhiz) > 0) {
  max(tukey_rhiz$y.position) + 0.5
} else {
  y_max_rhiz + 1
}


# Save statistics
aov_rhiz_export <- data.frame(
  test         = "One-way ANOVA",
  group1       = NA_character_,
  group2       = NA_character_,
  statistic    = summary(aov_rhiz)[[1]]$`F value`[1],
  p            = summary(aov_rhiz)[[1]]$`Pr(>F)`[1],
  p.adj        = NA_real_,
  p.adj.signif = NA_character_
)

tukey_rhiz_export <- as.data.frame(tukey_rhiz_all) %>%
  mutate(test = "Tukey HSD", statistic = NA_real_, p = NA_real_) %>%
  select(test, group1, group2, statistic, p, p.adj, p.adj.signif)

write_xlsx(
  bind_rows(aov_rhiz_export, tukey_rhiz_export),
  file.path(out_dir, "Alpha_diversity_statistics_rhizosphere_fertilizerlevels_Leonard_2022.xlsx")
)


# Figure
p_rhiz <- ggplot(
  alpha_df,
  aes(x = Fertilizer, y = Shannon, fill = Fertilizer)
) +
  geom_boxplot(
    alpha = 1, outlier.shape = NA,
    width = 0.6, colour = "black", size = 1.5
  ) +
  geom_jitter(
    aes(color = Fertilizer),
    width = 0.15, height = 0, size = 10, alpha = 0.8
  ) +
  {
    if (nrow(tukey_rhiz) > 0)
      stat_pvalue_manual(
        tukey_rhiz,
        label        = "p.adj.signif",
        y.position   = "y.position",
        tip.length   = 0.02,
        bracket.size = 1.2,
        size         = 30,
        fontface     = "bold",
        color        = "black",
        inherit.aes  = FALSE
      )
  } +
  scale_x_discrete(limits = fert_order, drop = FALSE) +
  scale_y_continuous(
    breaks = seq(0, 10, by = 2),
    limits = c(0, y_upper_rhiz),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_fill_manual(values  = fert_box_colors,   limits = fert_order, drop = FALSE) +
  scale_color_manual(values = fert_point_colors, limits = fert_order, drop = FALSE) +
  labs(x = "", y = "Shannon diversity") +
  facet_wrap(~ "Rhizosphere") +
  theme_cowplot() +
  theme(
    axis.line        = element_line(linewidth = 2),
    legend.position  = "right",
    axis.title       = element_text(size = 60, face = "bold"),
    axis.text.x      = element_text(size = 60, angle = 45, hjust = 1, face = "bold"),
    axis.text.y      = element_text(size = 60, face = "bold"),
    strip.text       = element_text(size = 60, face = "bold"),
    strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 1.5),
    panel.border     = element_rect(fill = NA)
  )

ggsave(
  filename = file.path(out_dir, "Alpha_diversity_fertilizerlevels_rhizosphere_Leonard_2022.png"),
  plot = p_rhiz, device = "png", width = 15, height = 15, units = "in", dpi = 600, bg = "white"
)


# Endosphere
ps <- Leonard.IDCY2.Endo.fert
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

f <- sample_data(ps)$Fertilizer |> as.character() |> trimws()
f <- dplyr::recode(f,
                   "None"   = "No Fertilizer",
                   "Medium" = "Medium Fertilizer",
                   "High"   = "High Fertilizer"
)
sample_data(ps)$Fertilizer <- fct_relevel(factor(f), fert_order)


# Calculate Shannon diversity
alpha_df <- estimate_richness(ps, measures = "Shannon")
alpha_df$SampleID <- rownames(alpha_df)
meta <- as(sample_data(ps), "data.frame")
meta$SampleID <- rownames(meta)
alpha_df <- left_join(alpha_df, meta, by = "SampleID")
alpha_df$Fertilizer <- fct_relevel(
  factor(trimws(as.character(alpha_df$Fertilizer))),
  fert_order
)


# Shapiro-Wilk normality test
sw_endo <- do.call(rbind, lapply(fert_order, function(g) {
  sub <- alpha_df[alpha_df$Fertilizer == g, ]
  sw  <- shapiro.test(sub$Shannon)
  data.frame(
    Compartment = "Endosphere",
    Fertilizer  = g,
    n           = nrow(sub),
    W           = round(sw$statistic, 4),
    p           = round(sw$p.value, 4),
    Normal      = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
}))

print(sw_endo)

## Normality met: one-way ANOVA + Tukey HSD
# One-way ANOVA
aov_endo <- aov(Shannon ~ Fertilizer, data = alpha_df)

print(summary(aov_endo))


# Tukey HSD post-hoc test
tukey_endo_all <- alpha_df %>%
  tukey_hsd(Shannon ~ Fertilizer) %>%
  add_significance("p.adj")

tukey_endo <- tukey_endo_all %>%
  filter(p.adj.signif != "ns")

print(as.data.frame(tukey_endo_all))


# Position significance brackets
y_max_endo <- max(alpha_df$Shannon, na.rm = TRUE)

tukey_endo <- tukey_endo %>%
  arrange(group1, group2) %>%
  mutate(
    y.position = y_max_endo + seq(0.25, by = 0.65, length.out = n())
  )

y_upper_endo <- if (nrow(tukey_endo) > 0) {
  max(tukey_endo$y.position) + 0.5
} else {
  y_max_endo + 1
}


# Save statistics
aov_endo_export <- data.frame(
  test         = "One-way ANOVA",
  group1       = NA_character_,
  group2       = NA_character_,
  statistic    = summary(aov_endo)[[1]]$`F value`[1],
  p            = summary(aov_endo)[[1]]$`Pr(>F)`[1],
  p.adj        = NA_real_,
  p.adj.signif = NA_character_
)

tukey_endo_export <- as.data.frame(tukey_endo_all) %>%
  mutate(test = "Tukey HSD", statistic = NA_real_, p = NA_real_) %>%
  select(test, group1, group2, statistic, p, p.adj, p.adj.signif)

write_xlsx(
  bind_rows(aov_endo_export, tukey_endo_export),
  file.path(out_dir, "Alpha_diversity_statistics_endosphere_fertilizerlevels_Leonard_2022.xlsx")
)


# Figure
p_endo <- ggplot(
  alpha_df,
  aes(x = Fertilizer, y = Shannon, fill = Fertilizer)
) +
  geom_boxplot(
    alpha = 1, outlier.shape = NA,
    width = 0.6, colour = "black", size = 1.5
  ) +
  geom_jitter(
    aes(color = Fertilizer),
    width = 0.15, height = 0, size = 10, alpha = 0.8
  ) +
  {
    if (nrow(tukey_endo) > 0)
      stat_pvalue_manual(
        tukey_endo,
        label        = "p.adj.signif",
        y.position   = "y.position",
        tip.length   = 0.02,
        bracket.size = 1.2,
        size         = 30,
        fontface     = "bold",
        color        = "black",
        inherit.aes  = FALSE
      )
  } +
  scale_x_discrete(limits = fert_order, drop = FALSE) +
  scale_y_continuous(
    breaks = seq(0, 10, by = 2),
    limits = c(0, y_upper_endo),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_fill_manual(values  = fert_box_colors,   limits = fert_order, drop = FALSE) +
  scale_color_manual(values = fert_point_colors, limits = fert_order, drop = FALSE) +
  labs(x = "", y = "Shannon diversity") +
  facet_wrap(~ "Endosphere") +
  theme_cowplot() +
  theme(
    axis.line        = element_line(linewidth = 2),
    legend.position  = "right",
    axis.title       = element_text(size = 60, face = "bold"),
    axis.text.x      = element_text(size = 60, angle = 45, hjust = 1, face = "bold"),
    axis.text.y      = element_text(size = 60, face = "bold"),
    strip.text       = element_text(size = 60, face = "bold"),
    strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 1.5),
    panel.border     = element_rect(fill = NA)
  )

ggsave(filename = file.path(out_dir, "Alpha_diversity_fertilizerlevels_endosphere_Leonard_2022.png"),
  plot = p_endo, device = "png", width = 15, height = 15, units = "in", dpi = 600, bg = "white"
)





# Shannon diversity analysis across fertilizer treatments within each plant compartment at Casselton in 2022
# Subset by plant compartments across the fertilizer levels

IDCY2.all <- subset_samples(ps.rare, Year == "Y2022")
Casselton.IDCY2 <- subset_samples(IDCY2.all, Location == "Casselton")
Casselton.IDCY2.Bulk.fert  <- subset_samples(Casselton.IDCY2, Compartments == "Bulk")
Casselton.IDCY2.Rhizo.fert <- subset_samples(Casselton.IDCY2, Compartments == "Rhizosphere")
Casselton.IDCY2.Endo.fert  <- subset_samples(Casselton.IDCY2, Compartments == "Endosphere")


# Bulk soil
ps <- Casselton.IDCY2.Bulk.fert
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

f <- sample_data(ps)$Fertilizer |> as.character() |> trimws()
f <- dplyr::recode(f,
                   "None"   = "No Fertilizer",
                   "Medium" = "Medium Fertilizer",
                   "High"   = "High Fertilizer"
)
sample_data(ps)$Fertilizer <- fct_relevel(factor(f), fert_order)

alpha_df <- estimate_richness(ps, measures = "Shannon")
alpha_df$SampleID <- rownames(alpha_df)
meta <- as(sample_data(ps), "data.frame")
meta$SampleID <- rownames(meta)
alpha_df <- left_join(alpha_df, meta, by = "SampleID")
alpha_df$Fertilizer <- fct_relevel(
  factor(trimws(as.character(alpha_df$Fertilizer))),
  fert_order
)

# Shapiro-Wilk normality test
sw_bulk <- do.call(rbind, lapply(fert_order, function(g) {
  sub <- alpha_df[alpha_df$Fertilizer == g, ]
  sw  <- shapiro.test(sub$Shannon)
  data.frame(
    Location = "Casselton", Compartment = "Bulk", Fertilizer = g,
    n = nrow(sub), W = round(sw$statistic, 4), p = round(sw$p.value, 4),
    Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
}))

print(sw_bulk)

# Normality met - one-way ANOVA + Tukey HSD
## One-way ANOVA
aov_bulk <- aov(Shannon ~ Fertilizer, data = alpha_df)

print(summary(aov_bulk))

# Tukey HSD post-hoc test
tukey_bulk_all <- alpha_df %>%
  tukey_hsd(Shannon ~ Fertilizer) %>%
  add_significance("p.adj")

print(as.data.frame(tukey_bulk_all))

aov_bulk_export <- data.frame(
  Location     = "Casselton",
  Compartment  = "Bulk",
  test         = "One-way ANOVA",
  group1       = NA_character_,
  group2       = NA_character_,
  statistic    = summary(aov_bulk)[[1]]$`F value`[1],
  p            = summary(aov_bulk)[[1]]$`Pr(>F)`[1],
  p.adj        = NA_real_,
  p.adj.signif = NA_character_
)

tukey_bulk_export <- as.data.frame(tukey_bulk_all) %>%
  mutate(
    Location    = "Casselton",
    Compartment = "Bulk",
    test        = "Tukey HSD",
    statistic   = NA_real_,
    p           = NA_real_
  ) %>%
  select(Location, Compartment, test, group1, group2, statistic, p, p.adj, p.adj.signif)


# Rhizosphere
ps <- Casselton.IDCY2.Rhizo.fert
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

f <- sample_data(ps)$Fertilizer |> as.character() |> trimws()
f <- dplyr::recode(f,
                   "None"   = "No Fertilizer",
                   "Medium" = "Medium Fertilizer",
                   "High"   = "High Fertilizer"
)
sample_data(ps)$Fertilizer <- fct_relevel(factor(f), fert_order)

alpha_df <- estimate_richness(ps, measures = "Shannon")
alpha_df$SampleID <- rownames(alpha_df)
meta <- as(sample_data(ps), "data.frame")
meta$SampleID <- rownames(meta)
alpha_df <- left_join(alpha_df, meta, by = "SampleID")
alpha_df$Fertilizer <- fct_relevel(
  factor(trimws(as.character(alpha_df$Fertilizer))),
  fert_order
)

# Shapiro-Wilk normality test
sw_rhiz <- do.call(rbind, lapply(fert_order, function(g) {
  sub <- alpha_df[alpha_df$Fertilizer == g, ]
  sw  <- shapiro.test(sub$Shannon)
  data.frame(
    Location = "Casselton", Compartment = "Rhizosphere", Fertilizer = g,
    n = nrow(sub), W = round(sw$statistic, 4), p = round(sw$p.value, 4),
    Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
}))

print(sw_rhiz)

# Normality met - one-way ANOVA + Tukey HSD
## One-way ANOVA
aov_rhiz <- aov(Shannon ~ Fertilizer, data = alpha_df)

print(summary(aov_rhiz))

# Tukey HSD post-hoc test
tukey_rhiz_all <- alpha_df %>%
  tukey_hsd(Shannon ~ Fertilizer) %>%
  add_significance("p.adj")

print(as.data.frame(tukey_rhiz_all))

aov_rhiz_export <- data.frame(
  Location     = "Casselton",
  Compartment  = "Rhizosphere",
  test         = "One-way ANOVA",
  group1       = NA_character_,
  group2       = NA_character_,
  statistic    = summary(aov_rhiz)[[1]]$`F value`[1],
  p            = summary(aov_rhiz)[[1]]$`Pr(>F)`[1],
  p.adj        = NA_real_,
  p.adj.signif = NA_character_
)

tukey_rhiz_export <- as.data.frame(tukey_rhiz_all) %>%
  mutate(
    Location    = "Casselton",
    Compartment = "Rhizosphere",
    test        = "Tukey HSD",
    statistic   = NA_real_,
    p           = NA_real_
  ) %>%
  select(Location, Compartment, test, group1, group2, statistic, p, p.adj, p.adj.signif)


# Endosphere
ps <- Casselton.IDCY2.Endo.fert
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

f <- sample_data(ps)$Fertilizer |> as.character() |> trimws()
f <- dplyr::recode(f,
                   "None"   = "No Fertilizer",
                   "Medium" = "Medium Fertilizer",
                   "High"   = "High Fertilizer"
)
sample_data(ps)$Fertilizer <- fct_relevel(factor(f), fert_order)

alpha_df <- estimate_richness(ps, measures = "Shannon")
alpha_df$SampleID <- rownames(alpha_df)
meta <- as(sample_data(ps), "data.frame")
meta$SampleID <- rownames(meta)
alpha_df <- left_join(alpha_df, meta, by = "SampleID")
alpha_df$Fertilizer <- fct_relevel(
  factor(trimws(as.character(alpha_df$Fertilizer))),
  fert_order
)

# Shapiro-Wilk normality test
sw_endo <- do.call(rbind, lapply(fert_order, function(g) {
  sub <- alpha_df[alpha_df$Fertilizer == g, ]
  sw  <- shapiro.test(sub$Shannon)
  data.frame(
    Location = "Casselton", Compartment = "Endosphere", Fertilizer = g,
    n = nrow(sub), W = round(sw$statistic, 4), p = round(sw$p.value, 4),
    Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
}))

print(sw_endo)

# Normality not met - Kruskal-Wallis + pairwise Wilcoxon with BH correction
#Kruskal-Wallis test
kw_endo <- alpha_df %>%
  kruskal_test(Shannon ~ Fertilizer)

print(as.data.frame(kw_endo))

# Pairwise Wilcoxon with BH correction
pwc_endo_all <- alpha_df %>%
  wilcox_test(Shannon ~ Fertilizer, p.adjust.method = "BH", paired = FALSE) %>%
  add_significance("p.adj")

print(as.data.frame(pwc_endo_all))

kw_endo_export <- as.data.frame(kw_endo) %>%
  mutate(
    Location     = "Casselton",
    Compartment  = "Endosphere",
    test         = "Kruskal-Wallis",
    group1       = NA_character_,
    group2       = NA_character_,
    p.adj        = NA_real_,
    p.adj.signif = NA_character_
  ) %>%
  select(Location, Compartment, test, group1, group2, statistic, p, p.adj, p.adj.signif)

pwc_endo_export <- as.data.frame(pwc_endo_all) %>%
  mutate(
    Location    = "Casselton",
    Compartment = "Endosphere",
    test        = "Pairwise Wilcoxon (BH adjusted)"
  ) %>%
  select(Location, Compartment, test, group1, group2, statistic, p, p.adj, p.adj.signif)


write_xlsx(bind_rows(aov_bulk_export, tukey_bulk_export, aov_rhiz_export, 
                     tukey_rhiz_export, kw_endo_export, pwc_endo_export), 
           file.path(out_dir, "Alpha_diversity_statistics_allcompartments_fertilizerlevels_Casselton_2022.xlsx"))




# Non-metric Multidimensional Scaling or NMDS analysis using Bray-Curtis distance of bacterial communities 
# across fertilizer treatments at Leonard in 2022

# Set seed for reproducibility
set.seed(777)

# Convert to relative abundance
ps.prop.Leo.IDCY2.all <- transform_sample_counts(Leonard.IDCY2, function(otu) otu / sum(otu))

# Bray-Curtis NMDS
ord.nmds.Leo.bray.IDCY2.all <- ordinate(
  ps.prop.Leo.IDCY2.all,
  method = "NMDS",
  distance = "bray"
)

# Fertilizer colors
custom_colors <- c(
  "None" = "#8B1C62",
  "Medium" = "#00008B",
  "High" = "#009E73"
)

# Compartment shapes
comp_shapes <- c(
  "Bulk" = 15,
  "Rhizosphere" = 19,
  "Endosphere" = 17
)

# Plot settings
POINT_SIZE <- 45
POINT_STROKE <- 3.5
JITTER_W <- 0.03
JITTER_H <- 0.03
TITLE_SIZE <- 130
AXIS_TITLE <- 130
AXIS_TEXT <- 140
FIG_W <- 38
FIG_H <- 25
FIG_DPI <- 600

# NMDS coordinates
ord_df <- plot_ordination(
  ps.prop.Leo.IDCY2.all,
  ord.nmds.Leo.bray.IDCY2.all,
  color = "Fertilizer",
  shape = "Compartments",
  justDF = TRUE
)


# NMDS plot
leonard_IDCY2 <- ggplot(
  ord_df,
  aes(x = NMDS1, y = NMDS2, color = Fertilizer, shape = Compartments)
) +
  theme_cowplot() +
  geom_point(
    size = POINT_SIZE,
    stroke = POINT_STROKE,
    position = position_jitter(width = JITTER_W, height = JITTER_H)
  ) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = comp_shapes) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme(
    axis.line = element_line(linewidth = 2),
    legend.position = "right",
    plot.title = element_text(size = TITLE_SIZE, face = "bold", hjust = 0.5),
    axis.title = element_text(size = AXIS_TITLE, face = "bold"),
    axis.text = element_text(size = AXIS_TEXT, face = "bold")
  )

ggsave(
  filename = file.path(out_dir, "NMDS_soybean_microbiome_across_fertilizerlevels_Leonard_2022.png"),
  plot = leonard_IDCY2,
  width = FIG_W,
  height = FIG_H,
  units = "in",
  dpi = FIG_DPI,
  bg = "white"
)



# PERMANOVA and pairwise PERMANOVA comparisons across fertilizer treatments within each compartment
# at Leonard and Casselton in 2022

sites <- c("Leonard", "Casselton")
compartments <- c("Bulk", "Rhizosphere", "Endosphere")

dispersion_list <- list()
permanova_list <- list()
pairwise_list <- list()

for (site in sites) {
  for (comp in compartments) {
    
ps_sub <- subset_samples(ps.rare, Location %in% site & Compartments %in% comp & Year %in% "Y2022")
 if (nsamples(ps_sub) < 4 || length(unique(sample_data(ps_sub)$Fertilizer)) < 2) next
    
bray_sub <- phyloseq::distance(ps_sub, method = "bray")
meta_sub <- data.frame(sample_data(ps_sub))
name <- paste(site, comp, sep = "_")
    
    
    
# Check dispersion across fertilizer treatments
set.seed(777)
disp_sub <- betadisper(bray_sub, meta_sub$Fertilizer)
dispersion_list[[name]] <- anova(disp_sub)
    

# PERMANOVA across fertilizer treatments
set.seed(777)
permanova_list[[name]] <- adonis2(bray_sub ~ Fertilizer, data = meta_sub, permutations = 999)
 
   
# Pairwise PERMANOVA across fertilizer treatments
set.seed(777)
pairwise_list[[name]] <- pairwise.adonis2(bray_sub ~ Fertilizer, data = meta_sub, permutations = 999, p.adjust.m = "bonferroni")
  }
}

# Save PERMANOVA and pairwise PERMANOVA results
pairwise_df <- do.call(rbind, lapply(names(pairwise_list), function(nm) {
  x <- pairwise_list[[nm]]
  do.call(rbind, lapply(names(x), function(pair) {
    if (!is.data.frame(x[[pair]])) return(NULL)
    cbind(Group = nm, Comparison = pair, x[[pair]])
  }))
}))

results <- c(
  list(Fertilizer_pairwise_PERMANOVA = pairwise_df),
  permanova_list
)

write.xlsx(results, file.path(out_dir,
      "PERMANOVA_and_pairwise_PERMANOVA_fertilizer_Leonard_Casselton_2022.xlsx"), overwrite = TRUE)







# Stacked barplot showing relative abundance of rhizosphere bacterial genera (top 40) across the fertilizer treatments 
#and soybean genotypes at Leonard in 2022


# Subset Leonard rhizosphere samples from 2022
IDCY2.all <- subset_samples(physeq, Year == "Y2022")
Leonard.IDCY2 <- subset_samples(IDCY2.all, Location == "Leonard")
Leonard.IDCY2.Rhizo <- subset_samples(Leonard.IDCY2, Compartments == "Rhizosphere")

# Keep taxa with genus-level classification
ps.genus <- subset_taxa(Leonard.IDCY2.Rhizo, !is.na(Genus) & Genus != "")

# Agglomerate taxa at the genus level
ps.genus <- tax_glom(ps.genus, taxrank = "Genus", NArm = FALSE)

# Keep the 40 most abundant genera
top40 <- names(sort(taxa_sums(ps.genus), decreasing = TRUE)[1:40])
ps.genus <- prune_taxa(top40, ps.genus)

# Convert counts to relative abundance
ps.rel_amp <- transform_sample_counts(ps.genus, function(x) x / sum(x))
# Sample names

if (!"Sample.names" %in% colnames(sample_data(ps.rel_amp))) {
  sample_data(ps.rel_amp)$Sample.names <- sample_names(ps.rel_amp)
}

# Fertilizer and genotype facets
sample_data(ps.rel_amp)$FertilizerFacet <- factor(
  c(
    "None" = "No Fertilizer",
    "Medium" = "Medium Fertilizer",
    "High" = "High Fertilizer"
  )[as.character(sample_data(ps.rel_amp)$Fertilizer)],
  levels = c("No Fertilizer", "Medium Fertilizer", "High Fertilizer")
)

sample_data(ps.rel_amp)$GenotypeFacet <- factor(
  as.character(sample_data(ps.rel_amp)$Genotype),
  levels = c("Resistant", "Sensitive")
)

# Color palette
my_colors <- c(
  "#E53935", "#3949AB", "#8E24AA", "#1E88E5", "#00ACC1", "#00897B", "#43A047", "#C0CA33", "#FDD835",
  "#FB8C00", "#00BFA5", "#6D4C41", "#FF8F00", "#303F9F", "#D81B60", "#5E35B1", "#546E7A", "#039BE5",
  "#00B8D4", "#F4511E", "#7CB342", "#CDDC39", "#FFEB3B", "#FFB300", "#FF7043", "#8D6E63", "#9E9D24",
  "#BDBDBD", "#AD1457", "#4527A0", "#283593", "#0277BD", "#00838F", "#00695C", "#2E7D32", "#78909C",
  "#F9A825", "#4E342E", "#D84315", "#757575"
)

# Stacked relative abundance plot
p <- plot_bar(ps.rel_amp, x = "Sample.names", fill = "Genus") +
  ggh4x::facet_nested(
    ~ FertilizerFacet + GenotypeFacet,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(values = my_colors) +
  labs(x = "Sample", y = "Relative abundance") +
  theme_bw(base_size = 26) +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 44, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_line(),
    axis.text.y = element_text(size = 44),
    axis.title = element_text(size = 44, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(0, "pt"),
    ggh4x.facet.nestline = element_line(linewidth = 2)
  )

#Save plot
ggsave(
  file.path(out_dir, "Stacked_barplot_top40_genus_across_fertilizerlevels_Leonard_rhizosphere_2022.png"),
  p, width = 24, height = 12, dpi = 600, bg = "white"
)



# Save genus legend separately

p <- plot_bar(ps.rel_amp, x = "Sample.names", fill = "Genus") +
  scale_fill_manual(values = my_colors) +
  theme_bw(base_size = 26)

get_only_legend <- function(myplot) {
  tmp <- ggplot_gtable(ggplot_build(myplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

legend_only <- get_only_legend(p)

ggsave(
  file.path(out_dir, "Legend_top40_genus_across_fertilizerlevels_Leonard_rhizosphere_2022.png"),
  plot = legend_only, width = 8, height = 12, units = "in", dpi = 600, bg = "white"
)




# Stacked barplot showing relative abundance of endosphere bacterial genera (top 40) across the fertilizer treatments 
#and soybean genotypes at Leonard in 2022

# Subset Leonard endosphere samples from 2022
IDCY2.all <- subset_samples(physeq, Year == "Y2022")
Leonard.IDCY2 <- subset_samples(IDCY2.all, Location == "Leonard")
Leonard.IDCY2.Endo <- subset_samples(Leonard.IDCY2, Compartments == "Endosphere")

# Keep taxa with genus-level classification
ps.genus <- subset_taxa(Leonard.IDCY2.Endo, !is.na(Genus) & Genus != "")

# Agglomerate taxa at the genus level
ps.genus <- tax_glom(ps.genus, taxrank = "Genus", NArm = FALSE)

# Keep the 40 most abundant genera
top40 <- names(sort(taxa_sums(ps.genus), decreasing = TRUE)[1:40])
ps.genus <- prune_taxa(top40, ps.genus)

# Convert counts to relative abundance
ps.rel_amp <- transform_sample_counts(ps.genus, function(x) x / sum(x))

# Ensure Sample.names exists
if (!"Sample.names" %in% colnames(sample_data(ps.rel_amp))) {
  sample_data(ps.rel_amp)$Sample.names <- sample_names(ps.rel_amp)
}

# Fertilizer and genotype facets
sample_data(ps.rel_amp)$FertilizerFacet <- factor(
  c(
    "None" = "No Fertilizer",
    "Medium" = "Medium Fertilizer",
    "High" = "High Fertilizer"
  )[as.character(sample_data(ps.rel_amp)$Fertilizer)],
  levels = c("No Fertilizer", "Medium Fertilizer", "High Fertilizer")
)

sample_data(ps.rel_amp)$GenotypeFacet <- factor(
  as.character(sample_data(ps.rel_amp)$Genotype),
  levels = c("Resistant", "Sensitive")
)

# 40-color genus palette
my_colors <- c(
  "#E53935", "#3949AB", "#8E24AA", "#1E88E5", "#00ACC1", "#00897B", "#43A047", "#C0CA33", "#FDD835",
  "#FB8C00", "#00BFA5", "#6D4C41", "#FF8F00", "#303F9F", "#D81B60", "#5E35B1", "#546E7A", "#039BE5",
  "#00B8D4", "#F4511E", "#7CB342", "#CDDC39", "#FFEB3B", "#FFB300", "#FF7043", "#8D6E63", "#9E9D24",
  "#BDBDBD", "#AD1457", "#4527A0", "#283593", "#0277BD", "#00838F", "#00695C", "#2E7D32", "#78909C",
  "#F9A825", "#4E342E", "#D84315", "#757575"
)

# Stacked relative abundance plot
p <- plot_bar(ps.rel_amp, x = "Sample.names", fill = "Genus") +
  ggh4x::facet_nested(~ FertilizerFacet + GenotypeFacet, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = my_colors) +
  labs(x = "", y = "Relative abundance") +
  theme_bw(base_size = 26) +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 44, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 44),
    axis.title = element_text(size = 44, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = grid::unit(0, "pt"),
    ggh4x.facet.nestline = element_line(linewidth = 2)
  )


# Save plot
ggsave(
  file.path(out_dir, "Stacked_barplot_top40_genus_across_fertilizerlevels_Leonard_endosphere_2022.png"),
  plot = p, width = 24, height = 12, units = "in", dpi = 600, bg = "white"
)



# Save genus legend separately
p_legend <- plot_bar(ps.rel_amp, x = "Sample.names", fill = "Genus") +
  scale_fill_manual(values = my_colors) +
  theme_bw(base_size = 26)

get_only_legend <- function(myplot) {
  tmp <- ggplot_gtable(ggplot_build(myplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

legend_only <- get_only_legend(p_legend)

ggsave(
  file.path(out_dir, "Legend_top40_genus_across_fertilizerlevels_Leonard_endosphere_2022.png"),
  plot = legend_only, width = 8, height = 12, units = "in", dpi = 600, bg = "white"
)




# Differntial abundance analysis by DESeq2 of rhizosphere bacterial genera between the no fertilizer and high fertilizer treatments 
# at Leonard in 2022


# Subset Leonard rhizosphere samples from 2022
IDCY2.all <- subset_samples(physeq, Year == "Y2022")
Leonard.IDCY2 <- subset_samples(IDCY2.all, Location == "Leonard")
Leonard.IDCY2.Rhizo <- subset_samples(Leonard.IDCY2, Compartments == "Rhizosphere")

# Keep None and High fertilizer treatments
sample_data(Leonard.IDCY2.Rhizo)$Fertilizer <- trimws(as.character(sample_data(Leonard.IDCY2.Rhizo)$Fertilizer))
sample_data(Leonard.IDCY2.Rhizo)$Fertilizer <- ifelse(sample_data(Leonard.IDCY2.Rhizo)$Fertilizer %in% c("None", "High"),
                                                      sample_data(Leonard.IDCY2.Rhizo)$Fertilizer, NA)

sample_data(Leonard.IDCY2.Rhizo)$Fertilizer <- factor(
  sample_data(Leonard.IDCY2.Rhizo)$Fertilizer,
  levels = c("None", "High")
)

Leonard.IDCY2.Rhizo <- prune_samples(
  !is.na(sample_data(Leonard.IDCY2.Rhizo)$Fertilizer),
  Leonard.IDCY2.Rhizo
)

# Agglomerate taxa at the genus level
ps_genus <- tax_glom(Leonard.IDCY2.Rhizo, taxrank = "Genus", NArm = TRUE)

# Keep taxa with valid genus names
tx <- as.data.frame(tax_table(ps_genus))
ps_genus <- prune_taxa(!is.na(tx$Genus) & nzchar(tx$Genus), ps_genus)

tx <- as.data.frame(tax_table(ps_genus))
tx$feature <- rownames(tx)

# DESeq2 analysis
dds <- phyloseq_to_deseq2(ps_genus, ~ Fertilizer)
dds <- estimateSizeFactors(dds, type = "poscounts")
dds$Fertilizer <- relevel(dds$Fertilizer, ref = "None")
dds <- DESeq(dds)

alpha <- 0.05
res_H <- results(dds, contrast = c("Fertilizer", "High", "None"), alpha = alpha)

# Keep significant genera
df_H <- as.data.frame(res_H)
df_H$feature <- rownames(df_H)

df_H <- df_H %>%
  left_join(tx[, c("feature", "Genus")], by = "feature") %>%
  filter(!is.na(padj), padj < alpha, !is.na(Genus), nzchar(Genus))

# Order genera by log2 fold change
genus_levels <- df_H %>%
  distinct(Genus, log2FoldChange) %>%
  arrange(log2FoldChange) %>%
  pull(Genus)

# Prepare plotting data
plot_long <- df_H %>%
  transmute(
    Genus = factor(Genus, levels = genus_levels),
    Enriched = ifelse(log2FoldChange >= 0, "High Fertilizer", "No Fertilizer"),
    log2FC = log2FoldChange,
    log10_baseMean = log10(baseMean + 1)
  ) %>%
  pivot_longer(cols = c(log2FC, log10_baseMean), names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("log2FC", "log10_baseMean"),
      labels = c("High fertilizer vs No fertilizer", "log10(BaseMean + 1)")
    ),
    FillGroup = ifelse(Metric == "log10(BaseMean + 1)", "BaseMean", Enriched)
  )

# Plot significant genera
p_merged <- ggplot(plot_long, aes(x = Genus, y = Value, fill = FillGroup)) +
  geom_col(width = 0.82) +
  coord_flip(clip = "off") +
  facet_grid(. ~ Metric, scales = "free_x", space = "free_x") +
  geom_hline(yintercept = 0, linewidth = 0.6, color = "grey40") +
  scale_fill_manual(
    values = c(
      "No Fertilizer" = "#8B1C62",
      "High Fertilizer" = "#009E73",
      "BaseMean" = "orange"
    ),
    name = NULL
  ) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_cowplot(font_size = 90) +
  theme(
    legend.position = "right",
    plot.title = element_blank(),
    strip.text = element_text(size = 60, face = "bold"),
    axis.text.y = element_text(size = 80, face = "bold", margin = margin(r = 8)),
    axis.text.x = element_text(size = 90, face = "bold"),
    axis.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 70),
    plot.margin = margin(t = 12, r = 24, b = 12, l = 12)
  )


# Save figure
ggsave(
  file.path(out_dir, "DESeq2_log2FC_and_log10BaseMean_high_vs_no_fertilizer_rhizosphere_Leonard_2022.png"),
  plot = p_merged, width = 35, height = 45, units = "in", dpi = 300, bg = "white"
)




# Differential abundance analysis by DESeq2 of endosphere bacterial genera between the no fertilizer and high fertilizer treatments
# at Leonard in 2022

# Subset Leonard endosphere samples from 2022
IDCY2.all <- subset_samples(physeq, Year == "Y2022")
Leonard.IDCY2 <- subset_samples(IDCY2.all, Location == "Leonard")
Leonard.IDCY2.Endo <- subset_samples(Leonard.IDCY2, Compartments == "Endosphere")

# Keep None and High fertilizer treatments
sample_data(Leonard.IDCY2.Endo)$Fertilizer <- trimws(as.character(sample_data(Leonard.IDCY2.Endo)$Fertilizer))
sample_data(Leonard.IDCY2.Endo)$Fertilizer <- ifelse(sample_data(Leonard.IDCY2.Endo)$Fertilizer %in% c("None", "High"),
                                                     sample_data(Leonard.IDCY2.Endo)$Fertilizer, NA)

sample_data(Leonard.IDCY2.Endo)$Fertilizer <- factor(
  sample_data(Leonard.IDCY2.Endo)$Fertilizer,
  levels = c("None", "High")
)

Leonard.IDCY2.Endo <- prune_samples(
  !is.na(sample_data(Leonard.IDCY2.Endo)$Fertilizer),
  Leonard.IDCY2.Endo
)

# Agglomerate taxa at the genus level
ps_genus <- tax_glom(Leonard.IDCY2.Endo, taxrank = "Genus", NArm = TRUE)

# Keep taxa with valid genus names
tx <- as.data.frame(tax_table(ps_genus))
ps_genus <- prune_taxa(!is.na(tx$Genus) & nzchar(tx$Genus), ps_genus)

tx <- as.data.frame(tax_table(ps_genus))
tx$feature <- rownames(tx)

# DESeq2 analysis
dds <- phyloseq_to_deseq2(ps_genus, ~ Fertilizer)
dds <- estimateSizeFactors(dds, type = "poscounts")
dds$Fertilizer <- relevel(dds$Fertilizer, ref = "None")
dds <- DESeq(dds)

alpha <- 0.05
res_H <- results(dds, contrast = c("Fertilizer", "High", "None"), alpha = alpha)

# Keep significant genera
df_H <- as.data.frame(res_H)
df_H$feature <- rownames(df_H)

df_H <- df_H %>%
  left_join(tx[, c("feature", "Genus")], by = "feature") %>%
  filter(!is.na(padj), padj < alpha, !is.na(Genus), nzchar(Genus))

# Order genera by log2 fold change
genus_levels <- df_H %>%
  distinct(Genus, log2FoldChange) %>%
  arrange(log2FoldChange) %>%
  pull(Genus)

# Prepare plotting data
plot_long <- df_H %>%
  transmute(
    Genus = factor(Genus, levels = genus_levels),
    Enriched = ifelse(log2FoldChange >= 0, "High Fertilizer", "No Fertilizer"),
    log2FC = log2FoldChange,
    log10_baseMean = log10(baseMean + 1)
  ) %>%
  pivot_longer(cols = c(log2FC, log10_baseMean), names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("log2FC", "log10_baseMean"),
      labels = c("High fertilizer vs No fertilizer", "log10(BaseMean + 1)")
    ),
    FillGroup = ifelse(Metric == "log10(BaseMean + 1)", "BaseMean", Enriched)
  )

# Plot significant genera
p_merged <- ggplot(plot_long, aes(x = Genus, y = Value, fill = FillGroup)) +
  geom_col(width = 0.82) +
  coord_flip(clip = "off") +
  facet_grid(. ~ Metric, scales = "free_x", space = "free_x") +
  geom_hline(yintercept = 0, linewidth = 0.6, color = "grey40") +
  scale_fill_manual(
    values = c(
      "No Fertilizer" = "#8B1C62",
      "High Fertilizer" = "#009E73",
      "BaseMean" = "orange"
    ),
    name = NULL
  ) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_cowplot(font_size = 90) +
  theme(
    legend.position = "right",
    plot.title = element_blank(),
    strip.text = element_text(size = 60, face = "bold"),
    axis.text.y = element_text(size = 80, face = "bold", margin = margin(r = 8)),
    axis.text.x = element_text(size = 90, face = "bold"),
    axis.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 70),
    plot.margin = margin(t = 12, r = 24, b = 12, l = 12)
  )


# Save figure
ggsave(
  file.path(out_dir, "DESeq2_log2FC_and_log10BaseMean_high_vs_no_fertilizer_endosphere_Leonard_2022.png"),
  plot = p_merged, width = 35, height = 45, units = "in", dpi = 300, bg = "white"
)
