# Soybean microbial community analyses across compartments, sites, and years
# Includes Shannon diversity, NMDS, PERMANOVA, and dbRDA analyses

# Load required packages
library(phyloseq)
library(vegan)
library(pairwiseAdonis)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(rstatix)
library(dplyr)
library(openxlsx)
library(writexl)
library(forcats)
library(ggpubr)
library(scales)

# Load field phyloseq objects
physeq <- readRDS("Condo7_PS.all_field_data(2021_to_2023)")
ps.rare <- readRDS("ps.rare.all.rds")

# Output directory
out_dir <- "Soybean_microbiome_across_compartments_sites_years"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Location order
newSorder <- c(
  "Prosper",
  "Casselton",
  "Colfax",
  "Leonard"
)


# Plot colors
box_colors <- c(
  "Prosper"   = "orange",
  "Casselton" = "red",
  "Colfax"    = "blue",
  "Leonard"   = "darkgreen"
)

point_colors <- c(
  "Prosper"   = "#7f4f00",
  "Casselton" = "#7f0000",
  "Colfax"    = "#00007f",
  "Leonard"   = "#003300"
)


#Comparisons of bacterial alpha diversity (Shannon index) across bulk soil, rhizosphere and endosphere 
#communities pooled across four sites and three years (2021-2023)

# 1. Shannon diversity across compartments
# Fertilizer = None; locations and years pooled

ps.none <- subset_samples(
  ps.rare,
  Fertilizer == "None"
)

ps.none <- prune_samples(
  sample_sums(ps.none) > 0,
  ps.none
)

ps.none <- prune_taxa(
  taxa_sums(ps.none) > 0,
  ps.none
)


# Calculate Shannon diversity
alpha <- estimate_richness(
  ps.none,
  measures = "Shannon"
)

alpha$SampleID <- rownames(alpha)

meta <- as(
  sample_data(ps.none),
  "data.frame"
)

meta$SampleID <- rownames(meta)

df <- left_join(
  alpha,
  meta,
  by = "SampleID"
)

df$Compartments <- factor(
  trimws(as.character(df$Compartments)),
  levels = c(
    "Bulk",
    "Rhizosphere",
    "Endosphere"
  )
)


# Shapiro-Wilk normality test
sw1_results <- do.call(
  rbind,
  lapply(
    c("Bulk", "Rhizosphere", "Endosphere"),
    function(comp) {
      sub <- df[df$Compartments == comp, ]
      sw  <- shapiro.test(sub$Shannon)
      data.frame(
        Compartment = comp,
        n      = nrow(sub),
        W      = round(sw$statistic, 4),
        p      = round(sw$p.value, 4),
        Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
      )
    }
  )
)

print(sw1_results)

# Rhizosphere failed normality; Kruskal-Wallis and
# pairwise Wilcoxon with BH correction

# Kruskal-Wallis test
kw1 <- df %>%
  kruskal_test(Shannon ~ Compartments)

print(as.data.frame(kw1))


# Pairwise Wilcoxon test with BH correction
pwc1_all <- df %>%
  wilcox_test(
    Shannon ~ Compartments,
    p.adjust.method = "BH",
    paired = FALSE
  ) %>%
  add_significance("p.adj")

pwc1 <- pwc1_all %>%
  filter(p.adj.signif != "ns")

print(as.data.frame(pwc1_all))


# Position significance brackets
y_max1 <- max(df$Shannon, na.rm = TRUE)

pwc1 <- pwc1 %>%
  arrange(group1, group2) %>%
  mutate(
    y.position = y_max1 + seq(0.25, by = 0.65, length.out = n())
  )

y_upper1 <- if (nrow(pwc1) > 0) {
  max(pwc1$y.position) + 0.5
} else {
  y_max1 + 1
}


# Save statistics
kw1_export <- as.data.frame(kw1) %>%
  mutate(
    test         = "Kruskal-Wallis",
    group1       = NA_character_,
    group2       = NA_character_,
    p.adj        = NA_real_,
    p.adj.signif = NA_character_
  ) %>%
  select(test, group1, group2, statistic, p, p.adj, p.adj.signif)

pwc1_export <- as.data.frame(pwc1_all) %>%
  mutate(test = "Pairwise Wilcoxon (BH adjusted)") %>%
  select(test, group1, group2, statistic, p, p.adj, p.adj.signif)

write_xlsx(
  bind_rows(kw1_export, pwc1_export),
  file.path(out_dir, "Alpha_diversity_statistics_fertilizerNone_compartments_pooled_by_site_years.xlsx")
)


# Figure
p1 <- ggplot(
  df,
  aes(x = Compartments, y = Shannon, fill = Compartments)
) +
  geom_boxplot(
    alpha = 1, outlier.shape = NA,
    width = 0.6, colour = "black", size = 3.5
  ) +
  geom_jitter(
    aes(color = Compartments),
    width = 0.15, height = 0, size = 16, alpha = 1
  ) +
  {
    if (nrow(pwc1) > 0)
      stat_pvalue_manual(
        pwc1,
        label        = "p.adj.signif",
        y.position   = "y.position",
        tip.length   = 0.02,
        bracket.size = 6,
        size         = 70,
        fontface     = "bold",
        color        = "black",
        inherit.aes  = FALSE
      )
  } +
  scale_y_continuous(
    breaks = seq(0, 10, by = 2),
    limits = c(0, y_upper1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(x = "", y = "Shannon diversity") +
  facet_wrap(~ "Compartments") +
  theme_cowplot() +
  theme(
    axis.line        = element_line(linewidth = 2),
    legend.position  = "right",
    axis.title       = element_text(size = 120, face = "bold"),
    axis.text.x      = element_text(size = 120, angle = 45, hjust = 1, face = "bold"),
    axis.text.y      = element_text(size = 120, face = "bold"),
    strip.text       = element_text(size = 120, face = "bold"),
    strip.background = element_rect(fill = "grey90", colour = NA),
    panel.border     = element_rect()
  )

ggsave(
  filename = file.path(out_dir, "Alpha_diversity_fertilizerNone_compartments_pooled_by_site_years.png"),
  plot = p1, device = "png", width = 24, height = 24, units = "in", dpi = 600, bg = "white"
)


#Comparison of bacterial alpha diversity (Shannon index) between four 
#locations within each plant compartment, pooled across three years (2021-2023)

# 2. Shannon diversity across locations - Bulk soil
#    Fertilizer = None; years pooled

Bulk <- subset_samples(
  ps.rare,
  Compartments == "Bulk" & Fertilizer == "None"
)

Bulk <- prune_samples(sample_sums(Bulk) > 0, Bulk)
Bulk <- prune_taxa(taxa_sums(Bulk) > 0, Bulk)

f <- sample_data(Bulk)$Location |> as.character() |> trimws()

sample_data(Bulk)$Location <- fct_relevel(factor(f), newSorder)


# Calculate Shannon diversity
alpha_bulk <- estimate_richness(Bulk, measures = "Shannon")
alpha_bulk$SampleID <- rownames(alpha_bulk)

meta_bulk <- as(sample_data(Bulk), "data.frame")
meta_bulk$SampleID <- rownames(meta_bulk)

alpha_bulk <- left_join(alpha_bulk, meta_bulk, by = "SampleID")

alpha_bulk$Location <- fct_relevel(
  factor(trimws(as.character(alpha_bulk$Location))),
  newSorder
)


# Shapiro-Wilk normality test
sw2_results <- do.call(
  rbind,
  lapply(newSorder, function(loc) {
    sub <- alpha_bulk[alpha_bulk$Location == loc, ]
    sw  <- shapiro.test(sub$Shannon)
    data.frame(
      Compartment = "Bulk",
      Location = loc,
      n      = nrow(sub),
      W      = round(sw$statistic, 4),
      p      = round(sw$p.value, 4),
      Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
    )
  })
)

print(sw2_results)

# Normality met; one-way ANOVA and Tukey HSD 
# One-way ANOVA
aov2 <- aov(Shannon ~ Location, data = alpha_bulk)

print(summary(aov2))


# Tukey HSD post-hoc test
tukey2_all <- alpha_bulk %>%
  tukey_hsd(Shannon ~ Location) %>%
  add_significance("p.adj")

tukey2 <- tukey2_all %>%
  filter(p.adj.signif != "ns")

print(as.data.frame(tukey2_all))


# Position significance brackets
y_max2 <- max(alpha_bulk$Shannon, na.rm = TRUE)

tukey2 <- tukey2 %>%
  arrange(group1, group2) %>%
  mutate(
    y.position = y_max2 + seq(0.25, by = 0.65, length.out = n())
  )

y_upper2 <- if (nrow(tukey2) > 0) {
  max(tukey2$y.position) + 0.5
} else {
  y_max2 + 1
}


# Save statistics
aov2_export <- data.frame(
  test         = "One-way ANOVA",
  group1       = NA_character_,
  group2       = NA_character_,
  statistic    = summary(aov2)[[1]]$`F value`[1],
  p            = summary(aov2)[[1]]$`Pr(>F)`[1],
  p.adj        = NA_real_,
  p.adj.signif = NA_character_
)

tukey2_export <- as.data.frame(tukey2_all) %>%
  mutate(test = "Tukey HSD", statistic = NA_real_, p = NA_real_) %>%
  select(test, group1, group2, statistic, p, p.adj, p.adj.signif)

write_xlsx(
  bind_rows(aov2_export, tukey2_export),
  file.path(out_dir, "Alpha_diversity_statistics_no_fertilizer_bulk_pooled_by_years.xlsx")
)


# Figure
p2 <- ggplot(
  alpha_bulk,
  aes(x = Location, y = Shannon, fill = Location)
) +
  geom_boxplot(
    alpha = 1, outlier.shape = NA,
    width = 0.6, colour = "black", size = 3.5
  ) +
  geom_jitter(
    aes(color = Location),
    width = 0.15, height = 0, size = 16, alpha = 1
  ) +
  {
    if (nrow(tukey2) > 0)
      stat_pvalue_manual(
        tukey2,
        label        = "p.adj.signif",
        y.position   = "y.position",
        tip.length   = 0.02,
        bracket.size = 6,
        size         = 70,
        fontface     = "bold",
        color        = "black",
        inherit.aes  = FALSE
      )
  } +
  scale_x_discrete(limits = newSorder, drop = FALSE) +
  scale_y_continuous(
    breaks = seq(0, 10, by = 2),
    limits = c(0, y_upper2),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_fill_manual(values  = box_colors,   limits = newSorder, drop = FALSE) +
  scale_color_manual(values = point_colors, limits = newSorder, drop = FALSE) +
  labs(x = "", y = "Shannon diversity") +
  facet_wrap(~ "Bulk soil") +
  theme_cowplot() +
  theme(
    axis.line        = element_line(linewidth = 2),
    legend.position  = "right",
    axis.title       = element_text(size = 120, face = "bold"),
    axis.text.x      = element_text(size = 120, angle = 45, hjust = 1, face = "bold"),
    axis.text.y      = element_text(size = 120, face = "bold"),
    strip.text       = element_text(size = 120, face = "bold"),
    strip.background = element_rect(fill = "grey90", colour = NA),
    panel.border     = element_rect()
  )

ggsave(
  filename = file.path(out_dir, "Alpha_diversity_no_fertilizer_bulk_pooled_by_years.png"),
  plot = p2, device = "png", width = 24, height = 24, units = "in", dpi = 600, bg = "white"
)



# 3. Shannon diversity across locations - Rhizosphere
#    Fertilizer = None; years pooled

Rhizosphere <- subset_samples(
  ps.rare,
  Compartments == "Rhizosphere" & Fertilizer == "None"
)

Rhizosphere <- prune_samples(sample_sums(Rhizosphere) > 0, Rhizosphere)
Rhizosphere <- prune_taxa(taxa_sums(Rhizosphere) > 0, Rhizosphere)

f <- sample_data(Rhizosphere)$Location |> as.character() |> trimws()

sample_data(Rhizosphere)$Location <- fct_relevel(factor(f), newSorder)


# Calculate Shannon diversity
alpha_rhiz <- estimate_richness(Rhizosphere, measures = "Shannon")
alpha_rhiz$SampleID <- rownames(alpha_rhiz)

meta_rhiz <- as(sample_data(Rhizosphere), "data.frame")
meta_rhiz$SampleID <- rownames(meta_rhiz)

alpha_rhiz <- left_join(alpha_rhiz, meta_rhiz, by = "SampleID")

alpha_rhiz$Location <- fct_relevel(
  factor(trimws(as.character(alpha_rhiz$Location))),
  newSorder
)


# Shapiro-Wilk normality test
sw3_results <- do.call(
  rbind,
  lapply(newSorder, function(loc) {
    sub <- alpha_rhiz[alpha_rhiz$Location == loc, ]
    sw  <- shapiro.test(sub$Shannon)
    data.frame(
      Compartment = "Rhizosphere",
      Location = loc,
      n      = nrow(sub),
      W      = round(sw$statistic, 4),
      p      = round(sw$p.value, 4),
      Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
    )
  })
)

print(sw3_results)


# Normality not met; Kruskal-Wallis and
# pairwise Wilcoxon with BH correction 


# Kruskal-Wallis test
kw3 <- alpha_rhiz %>%
  kruskal_test(Shannon ~ Location)

print(as.data.frame(kw3))


# Pairwise Wilcoxon test with BH correction
pwc3_all <- alpha_rhiz %>%
  wilcox_test(
    Shannon ~ Location,
    p.adjust.method = "BH",
    paired = FALSE
  ) %>%
  add_significance("p.adj")

pwc3 <- pwc3_all %>%
  filter(p.adj.signif != "ns")

print(as.data.frame(pwc3_all))


# Position significance brackets
y_max3 <- max(alpha_rhiz$Shannon, na.rm = TRUE)

pwc3 <- pwc3 %>%
  arrange(group1, group2) %>%
  mutate(
    y.position = y_max3 + seq(0.25, by = 0.65, length.out = n())
  )

y_upper3 <- if (nrow(pwc3) > 0) {
  max(pwc3$y.position) + 0.5
} else {
  y_max3 + 1
}


# Save statistics
kw3_export <- as.data.frame(kw3) %>%
  mutate(
    test         = "Kruskal-Wallis",
    group1       = NA_character_,
    group2       = NA_character_,
    p.adj        = NA_real_,
    p.adj.signif = NA_character_
  ) %>%
  select(test, group1, group2, statistic, p, p.adj, p.adj.signif)

pwc3_export <- as.data.frame(pwc3_all) %>%
  mutate(test = "Pairwise Wilcoxon (BH adjusted)") %>%
  select(test, group1, group2, statistic, p, p.adj, p.adj.signif)

write_xlsx(
  bind_rows(kw3_export, pwc3_export),
  file.path(out_dir, "Alpha_diversity_statistics_no_fertilizer_rhizosphere_pooled_by_years.xlsx")
)


# Figure
p3 <- ggplot(
  alpha_rhiz,
  aes(x = Location, y = Shannon, fill = Location)
) +
  geom_boxplot(
    alpha = 1, outlier.shape = NA,
    width = 0.6, colour = "black", size = 3.5
  ) +
  geom_jitter(
    aes(color = Location),
    width = 0.15, height = 0, size = 16, alpha = 1
  ) +
  {
    if (nrow(pwc3) > 0)
      stat_pvalue_manual(
        pwc3,
        label        = "p.adj.signif",
        y.position   = "y.position",
        tip.length   = 0.02,
        bracket.size = 6,
        size         = 70,
        fontface     = "bold",
        color        = "black",
        inherit.aes  = FALSE
      )
  } +
  scale_x_discrete(limits = newSorder, drop = FALSE) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 5),
    limits = c(0, y_upper3),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_fill_manual(values  = box_colors,   limits = newSorder, drop = FALSE) +
  scale_color_manual(values = point_colors, limits = newSorder, drop = FALSE) +
  labs(x = "", y = "Shannon diversity") +
  facet_wrap(~ "Rhizosphere") +
  theme_cowplot() +
  theme(
    axis.line        = element_line(linewidth = 2),
    legend.position  = "right",
    axis.title       = element_text(size = 120, face = "bold"),
    axis.text.x      = element_text(size = 120, angle = 45, hjust = 1, face = "bold"),
    axis.text.y      = element_text(size = 120, face = "bold"),
    strip.text       = element_text(size = 120, face = "bold"),
    strip.background = element_rect(fill = "grey90", colour = NA),
    panel.border     = element_rect()
  )

ggsave(
  filename = file.path(out_dir, "Alpha_diversity_no_fertilizer_rhizosphere_pooled_by_years.png"),
  plot = p3, device = "png", width = 24, height = 24, units = "in", dpi = 600, bg = "white"
)



# 4. Shannon diversity across locations - Endosphere
#   Fertilizer = None; years pooled

Endosphere <- subset_samples(
  ps.rare,
  Compartments == "Endosphere" & Fertilizer == "None"
)

Endosphere <- prune_samples(sample_sums(Endosphere) > 0, Endosphere)
Endosphere <- prune_taxa(taxa_sums(Endosphere) > 0, Endosphere)

f <- sample_data(Endosphere)$Location |> as.character() |> trimws()

sample_data(Endosphere)$Location <- fct_relevel(factor(f), newSorder)


# Calculate Shannon diversity
alpha_endo <- estimate_richness(Endosphere, measures = "Shannon")
alpha_endo$SampleID <- rownames(alpha_endo)

meta_endo <- as(sample_data(Endosphere), "data.frame")
meta_endo$SampleID <- rownames(meta_endo)

alpha_endo <- left_join(alpha_endo, meta_endo, by = "SampleID")

alpha_endo$Location <- fct_relevel(
  factor(trimws(as.character(alpha_endo$Location))),
  newSorder
)


# Shapiro-Wilk normality test
sw4_results <- do.call(
  rbind,
  lapply(newSorder, function(loc) {
    sub <- alpha_endo[alpha_endo$Location == loc, ]
    sw  <- shapiro.test(sub$Shannon)
    data.frame(
      Compartment = "Endosphere",
      Location = loc,
      n      = nrow(sub),
      W      = round(sw$statistic, 4),
      p      = round(sw$p.value, 4),
      Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
    )
  })
)

print(sw4_results)


# Normality met; one-way ANOVA and Tukey HSD 
# One-way ANOVA
aov4 <- aov(Shannon ~ Location, data = alpha_endo)
print(summary(aov4))


# Tukey HSD post-hoc test
tukey4_all <- alpha_endo %>%
  tukey_hsd(Shannon ~ Location) %>%
  add_significance("p.adj")

tukey4 <- tukey4_all %>%
  filter(p.adj.signif != "ns")

print(as.data.frame(tukey4_all))


# Position significance brackets
y_max4 <- max(alpha_endo$Shannon, na.rm = TRUE)

tukey4 <- tukey4 %>%
  arrange(group1, group2) %>%
  mutate(
    y.position = y_max4 + seq(0.25, by = 0.55, length.out = n())
  )

y_upper4 <- if (nrow(tukey4) > 0) {
  max(tukey4$y.position) + 0.5
} else {
  y_max4 + 1
}


# Save statistics
aov4_export <- data.frame(
  test         = "One-way ANOVA",
  group1       = NA_character_,
  group2       = NA_character_,
  statistic    = summary(aov4)[[1]]$`F value`[1],
  p            = summary(aov4)[[1]]$`Pr(>F)`[1],
  p.adj        = NA_real_,
  p.adj.signif = NA_character_
)

tukey4_export <- as.data.frame(tukey4_all) %>%
  mutate(test = "Tukey HSD", statistic = NA_real_, p = NA_real_) %>%
  select(test, group1, group2, statistic, p, p.adj, p.adj.signif)

write_xlsx(
  bind_rows(aov4_export, tukey4_export),
  file.path(out_dir, "Alpha_diversity_statistics_no_fertilizer_endosphere_pooled_by_years.xlsx")
)


# Figure
p4 <- ggplot(
  alpha_endo,
  aes(x = Location, y = Shannon, fill = Location)
) +
  geom_boxplot(
    alpha = 1, outlier.shape = NA,
    width = 0.6, colour = "black", size = 3.5
  ) +
  geom_jitter(
    aes(color = Location),
    width = 0.15, height = 0, size = 16, alpha = 1
  ) +
  {
    if (nrow(tukey4) > 0)
      stat_pvalue_manual(
        tukey4,
        label        = "p.adj.signif",
        y.position   = "y.position",
        tip.length   = 0.02,
        bracket.size = 6,
        size         = 70,
        fontface     = "bold",
        color        = "black",
        inherit.aes  = FALSE
      )
  } +
  scale_x_discrete(limits = newSorder, drop = FALSE) +
  scale_y_continuous(
    breaks = seq(0, 10, by = 2),
    limits = c(0, y_upper4),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_fill_manual(values  = box_colors,   limits = newSorder, drop = FALSE) +
  scale_color_manual(values = point_colors, limits = newSorder, drop = FALSE) +
  labs(x = "", y = "Shannon diversity") +
  facet_wrap(~ "Endosphere") +
  theme_cowplot() +
  theme(
    axis.line        = element_line(linewidth = 2),
    legend.position  = "right",
    axis.title       = element_text(size = 120, face = "bold"),
    axis.text.x      = element_text(size = 120, angle = 45, hjust = 1, face = "bold"),
    axis.text.y      = element_text(size = 120, face = "bold"),
    strip.text       = element_text(size = 120, face = "bold"),
    strip.background = element_rect(fill = "grey90", colour = NA),
    panel.border     = element_rect()
  )

ggsave(
  filename = file.path(out_dir, "Alpha_diversity_no_fertilizer_endosphere_pooled_by_years.png"),
  plot = p4, device = "png", width = 24, height = 24, units = "in", dpi = 600, bg = "white"
)



# Non-metric Multidimensional Scaling or NMDS analysis using Bray-Curtis distance of bacterial communities across the plant compartments 
#and the sites within each compartment over the years (2021-2023)

# Set seed for reproducibility
set.seed(777)

# Subset no-fertilizer samples
None <- subset_samples(ps.rare, Fertilizer %in% c("None"))

# Convert to relative abundance
ps.prop.none <- transform_sample_counts(None, function(otu) otu / sum(otu))

# Bray-Curtis NMDS
ord.nmds.bray.none <- ordinate(
  ps.prop.none,
  method = "NMDS",
  distance = "bray"
)

# Location colors
custom_colors <- c(
  "Casselton" = "red",
  "Colfax" = "blue",
  "Leonard" = "darkgreen",
  "Prosper" = "orange")

# Plot settings
POINT_SIZE <- 34
POINT_STROKE <- 3.5
JITTER_W <- 0.03
JITTER_H <- 0.03
TITLE_SIZE <- 160
AXIS_TITLE <- 160
AXIS_TEXT <- 160
FIG_W <- 48
FIG_H <- 32
FIG_DPI <- 600

# NMDS coordinates
ord_df <- plot_ordination(
  ps.prop.none,
  ord.nmds.bray.none,
  color = "Location",
  shape = "Samples",
  justDF = TRUE
)

# NMDS plot
location.none_pub_noleg <- ggplot(
  ord_df,
  aes(x = NMDS1, y = NMDS2, color = Location, shape = Samples)
) +
  theme_cowplot() +
  geom_point(
    size = POINT_SIZE,
    stroke = POINT_STROKE,
    position = position_jitter(width = JITTER_W, height = JITTER_H)
  ) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = c(8, 3, 4, 22, 1, 2, 15, 19, 17)) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme(
    axis.line = element_line(linewidth = 2),
    legend.position = "right",
    plot.title = element_text(size = TITLE_SIZE, face = "bold", hjust = 0.5),
    axis.title = element_text(size = AXIS_TITLE, face = "bold"),
    axis.text = element_text(size = AXIS_TEXT, face = "bold")
  )

ggsave(
  filename = file.path(out_dir, "NMDS_soybean_microbiome_across_sites_years.png"),
  plot = location.none_pub_noleg,
  device = "png",
  width = FIG_W,
  height = FIG_H,
  units = "in",
  dpi = FIG_DPI,
  limitsize = FALSE,
  bg = "white"
)



# PERMANOVA and pairwise PERMANOVA comparisons of bacterial communities across the plant compartments 
#and the sites within each compartment over the years (2021-2023)


# PERMANOVA across compartments, pooled by years
bray_none <- phyloseq::distance(ps.prop.none, method = "bray")
meta_none <- data.frame(sample_data(ps.prop.none))
meta_none$Compartments <- factor(trimws(as.character(meta_none$Compartments)), levels = c("Bulk", "Rhizosphere", "Endosphere"))

# Check dispersion across compartments
set.seed(777)
disp_comp <- betadisper(bray_none, meta_none$Compartments)
disp_comp_anova <- anova(disp_comp)
print(disp_comp_anova)

# PERMANOVA across compartments
set.seed(777)
perm_compartment <- adonis2(bray_none ~ Compartments, data = meta_none, permutations = 999)
print(perm_compartment)

# Pairwise PERMANOVA across compartments
set.seed(777)
pair_compartment <- pairwise.adonis2(bray_none ~ Compartments, data = meta_none, permutations = 999, p.adjust.m = "bonferroni")
print(pair_compartment)

# PERMANOVA across locations within each compartment and year
compartments <- c("Bulk", "Rhizosphere", "Endosphere")
years <- c("Y2021", "Y2022", "Y2023")

dispersion_list <- list()
permanova_list <- list()
pairwise_list <- list()

for (comp in compartments) {
  for (yr in years) {
    
ps_sub <- subset_samples(None, Compartments %in% comp & Year %in% yr)
if (nsamples(ps_sub) < 4 || length(unique(sample_data(ps_sub)$Location)) < 2) next
    
bray_sub <- phyloseq::distance(ps_sub, method = "bray")
meta_sub <- data.frame(sample_data(ps_sub))
name <- paste(comp, yr, sep = "_")
    
# Check dispersion across locations
set.seed(777)
disp_sub <- betadisper(bray_sub, meta_sub$Location)
dispersion_list[[name]] <- anova(disp_sub)    
    
# PERMANOVA across locations
set.seed(777)
permanova_list[[name]] <- adonis2(bray_sub ~ Location, data = meta_sub, permutations = 999)
    

# Pairwise PERMANOVA across locations
set.seed(777)
pairwise_list[[name]] <- pairwise.adonis2(bray_sub ~ Location, data = meta_sub, permutations = 999, p.adjust.m = "bonferroni")
  }
}

# Save PERMANOVA and pairwise PERMANOVA results
pair_compartment_df <- do.call(rbind, lapply(names(pair_compartment), function(nm) {
  x <- pair_compartment[[nm]]
  if (!is.data.frame(x)) return(NULL)
  cbind(Comparison = nm, x)
}))

pairwise_df <- do.call(rbind, lapply(names(pairwise_list), function(nm) {
  x <- pairwise_list[[nm]]
  do.call(rbind, lapply(names(x), function(pair) {
    if (!is.data.frame(x[[pair]])) return(NULL)
    cbind(Group = nm, Comparison = pair, x[[pair]])
  }))
}))

results <- c(
  list(Compartment_PERMANOVA = as.data.frame(perm_compartment),
       Compartment_pairwise_PERMANOVA = pair_compartment_df,
       Location_pairwise_PERMANOVA = pairwise_df),
  permanova_list
)

write.xlsx(results, file.path(out_dir, 
    "PERMANOVA_and_pairwise_PERMANOVA_soybean_microbiome_across_sites_years.xlsx"), overwrite = TRUE)




# distance-based Redundancy Analysis or dbRDA analysis of soil variables associated with soybean rhizosphere and endosphere microbial
# community variations across the sites and years

# dbRDA Function
run_dbrda <- function(ps.rare, compartment_name) {
  
# Extract tables
otu_c <- as.data.frame(otu_table(ps.rare))
meta_c <- as.data.frame(sample_data(ps.rare))
names(meta_c) <- make.names(names(meta_c))
  
# Remove zero-sum samples
keep <- rowSums(otu_c) > 0
otu_c <- otu_c[keep, ]
meta_c <- meta_c[keep, ]
  
# Environmental variables
env <- data.frame(
CCE = as.numeric(meta_c$CCE),
IDC_scores = as.numeric(meta_c$IDC.scores),
pH = as.numeric(meta_c$pH),
NO3N = as.numeric(meta_c$NO3N),
P = as.numeric(meta_c$P),
K = as.numeric(meta_c$K),
OM = as.numeric(meta_c$OM),
SO4 = as.numeric(meta_c$SO4),
Na = as.numeric(meta_c$Na),
NH4 = as.numeric(meta_c$NH4),
row.names = rownames(meta_c))


# Remove samples with missing environmental data
complete <- complete.cases(env)
cat("Samples retained after removing NAs:", sum(complete), "/", nrow(env), "\n")
env_scaled <- as.data.frame(scale(env[complete, ]))
otu_cc <- otu_c[complete, ]
meta_cc <- meta_c[complete, ]
  
# Bray-Curtis and dbRDA
bc <- vegdist(otu_cc, method = "bray")
mod <- dbrda(bc ~ ., data = env_scaled)
  
print(summary(mod))
  
# Permutation test
sig_table <- as.data.frame(anova(mod, by = "margin", permutations = 999))
sig_table$Variable <- rownames(sig_table)
write_xlsx(sig_table, file.path(out_dir, paste0("dbrda_sig_", compartment_name, ".xlsx")))
  
print(sig_table)
  
# Axis variance explained
eig <- eigenvals(mod)
eig_pos <- eig[eig > 0]
ax1_pct <- round(eig_pos[1] / sum(eig_pos) * 100, 1)
ax2_pct <- round(eig_pos[2] / sum(eig_pos) * 100, 1)
  
# Scores
site_sc <- as.data.frame(scores(mod, display = "sites"))
site_sc$Location <- meta_cc$Location
site_sc$Samples <- meta_cc$Samples
bp_sc <- as.data.frame(scores(mod, display = "bp"))
scale_factor <- min(diff(range(site_sc$dbRDA1)), diff(range(site_sc$dbRDA2))) * 0.55
bp_plot <- bp_sc * scale_factor
bp_plot$label <- rownames(bp_plot)
  
# Plot limits
x_range <- range(c(site_sc$dbRDA1, bp_plot$dbRDA1))
y_range <- range(c(site_sc$dbRDA2, bp_plot$dbRDA2))
x_pad <- diff(x_range) * 0.50
y_pad <- diff(y_range) * 0.50
  
# Plot settings
shape_vals <- c(
"Rhizosphere_2021" = 15, "Rhizosphere_2022" = 19, "Rhizosphere_2023" = 17,
"Endosphere_2021" = 0, "Endosphere_2022" = 1, "Endosphere_2023" = 2,
"Bulk_2021" = 3, "Bulk_2022" = 4, "Bulk_2023" = 8)

custom_colors <- c(
    "Casselton" = "red",
    "Colfax" = "blue",
    "Leonard" = "darkgreen",
    "Prosper" = "orange")
  
# Plot
set.seed(42)
      p <- ggplot() +
      geom_segment(
      data = bp_plot,
      aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
      arrow = arrow(length = unit(0.4, "cm"), type = "closed"),
      color = "black",
      linewidth = 1.0
    ) +
    geom_jitter(
      data = site_sc,
      aes(x = dbRDA1, y = dbRDA2, color = Location, shape = Samples),
      width = 0.6,
      height = 0.6,
      size = 7,
      stroke = 1.5,
      alpha = 0.80
    ) +
    geom_text_repel(
      data = bp_plot,
      aes(x = dbRDA1, y = dbRDA2, label = label),
      size = 12,
      fontface = "bold",
      color = "black",
      bg.color = "white",
      bg.r = 0.25,
      box.padding = 2.0,
      point.padding = 1.0,
      segment.color = "grey30",
      segment.size = 0.8,
      segment.linetype = "solid",
      max.overlaps = Inf,
      force = 15,
      force_pull = 0.2,
      min.segment.length = 0,
      xlim = c(x_range[1] - x_pad, x_range[2] + x_pad),
      ylim = c(y_range[1] - y_pad, y_range[2] + y_pad),
      seed = 42
    ) +
    scale_color_manual(values = custom_colors, name = "Location") +
    scale_shape_manual(values = shape_vals, name = "Sample x Year") +
    coord_cartesian(
      xlim = c(x_range[1] - x_pad, x_range[2] + x_pad),
      ylim = c(y_range[1] - y_pad, y_range[2] + y_pad)
    ) +
    labs(
      x = paste0("dbRDA1 (", ax1_pct, "%)"),
      y = paste0("dbRDA2 (", ax2_pct, "%)"),
      title = paste0(compartment_name, " - dbRDA")
    ) +
    theme_cowplot(font_size = 16) +
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 15, face = "bold"),
      legend.position = "right"
    )
  
  
# Save plot
fname_png <- file.path(out_dir, paste0("dbrda_", compartment_name, "_final.png"))
ggsave(fname_png, plot = p, width = 12, height = 9, units = "in", dpi = 600, bg = "white")
return(list(model = mod, plot = p, sig = sig_table))}


# Run dbRDA
ps.rhizo <- subset_samples(ps.rare, Compartments == "Rhizosphere" & Fertilizer == "None")
res_rhizo <- run_dbrda(ps.rhizo, "Rhizosphere")
ps.endo <- subset_samples(ps.rare, Compartments == "Endosphere" & Fertilizer == "None")
res_endo <- run_dbrda(ps.endo, "Endosphere")

