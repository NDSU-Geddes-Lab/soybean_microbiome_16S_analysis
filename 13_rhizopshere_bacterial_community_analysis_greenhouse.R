# Rhizosphere bacterial community analysis across different fertilizer treatments 
# in greenhouse experiment

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



# Load greenhouse phyloseq objects
greenhouse.physeq <- readRDS("Final.ps.GH.Rhizo.2021.rds")
greenhouse.ps.rare <- readRDS("Final.ps.GH.Rhizo.2021.rarefied.rds")


#output directory
out_dir <- "Rhizosphere_bacterial_community_analysis_greenhouse"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# Fertilizer order and colors
fert_order <- c("No Fertilizer", "Medium Fertilizer", "High Fertilizer")

box_colors <- c(
  "No Fertilizer" = "#8B1C62",
  "Medium Fertilizer" = "#00008B",
  "High Fertilizer" = "#009E73"
)

jitter_colors <- c(
  "No Fertilizer" = "black",
  "Medium Fertilizer" = "black",
  "High Fertilizer" = "black"
)


# Subset rhizosphere samples from greenhouse experiment
ps <- subset_samples(greenhouse.ps.rare, Location %in% "Jay")
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

f <- sample_data(ps)$Fertilizer |> as.character() |> trimws()
f <- dplyr::recode(f,
                   "None"   = "No Fertilizer",
                   "Medium" = "Medium Fertilizer",
                   "High"   = "High Fertilizer"
)
sample_data(ps)$Fertilizer <- fct_relevel(factor(f), fert_order)
ls()

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
sw_results <- do.call(rbind, lapply(fert_order, function(g) {
  sub <- alpha_df[alpha_df$Fertilizer == g, ]
  sw  <- shapiro.test(sub$Shannon)
  data.frame(
    Fertilizer = g,
    n          = nrow(sub),
    W          = round(sw$statistic, 4),
    p          = round(sw$p.value, 4),
    Normal     = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
}))

print(sw_results)


# Normality met - one-way ANOVA + Tukey HSD
## One-way ANOVA
aov_result <- aov(Shannon ~ Fertilizer, data = alpha_df)

print(summary(aov_result))


# Tukey HSD post-hoc test
tukey_all <- alpha_df %>%
  tukey_hsd(Shannon ~ Fertilizer) %>%
  add_significance("p.adj")

tukey <- tukey_all %>%
  filter(p.adj.signif != "ns") %>%
  mutate(p.adj.signif = ifelse(p.adj.signif == "****", "***", p.adj.signif))

print(as.data.frame(tukey_all))


# Position significance brackets
y_max <- max(alpha_df$Shannon, na.rm = TRUE)

tukey <- tukey %>%
  arrange(group1, group2) %>%
  mutate(
    y.position = y_max + seq(0.25, by = 0.55, length.out = n())
  )

y_upper <- if (nrow(tukey) > 0) {
  max(tukey$y.position) + 0.5
} else {
  y_max + 1
}


# Save statistics
aov_export <- data.frame(
  test         = "One-way ANOVA",
  group1       = NA_character_,
  group2       = NA_character_,
  statistic    = summary(aov_result)[[1]]$`F value`[1],
  p            = summary(aov_result)[[1]]$`Pr(>F)`[1],
  p.adj        = NA_real_,
  p.adj.signif = NA_character_
)

tukey_export <- as.data.frame(tukey_all) %>%
  mutate(test = "Tukey HSD", statistic = NA_real_, p = NA_real_) %>%
  select(test, group1, group2, statistic, p, p.adj, p.adj.signif)

write_xlsx(
  bind_rows(aov_export, tukey_export),
  file.path(out_dir, "Alpha_diversity_statistics_fertilizerlevels_rhizosphere_greenhouse.xlsx")
)



# Figure
richness.rhizo <- ggplot(
  alpha_df,
  aes(x = Fertilizer, y = Shannon, fill = Fertilizer)
) +
  geom_boxplot(
    alpha = 1, outlier.shape = NA,
    width = 0.6, colour = "black", size = 3.5
  ) +
  geom_jitter(
    aes(color = Fertilizer),
    width = 0.15, height = 0,
    size = 16, alpha = 1
  ) +
  {
    if (nrow(tukey) > 0)
      stat_pvalue_manual(
        tukey,
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
  scale_x_discrete(limits = fert_order, drop = FALSE) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 5),
    limits = c(0, y_upper),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_fill_manual(values  = box_colors,    limits = fert_order, drop = FALSE) +
  scale_color_manual(values = jitter_colors, limits = fert_order, drop = FALSE) +
  labs(x = "", y = "Shannon diversity") +
  facet_wrap(~ "Rhizosphere") +
  theme_cowplot() +
  theme(
    axis.line        = element_line(linewidth = 2),
    legend.position  = "right",
    legend.text      = element_text(size = 160, face = "bold"),
    legend.title     = element_text(size = 160, face = "bold"),
    axis.title       = element_text(size = 160, face = "bold"),
    axis.text.x      = element_text(size = 160, angle = 45, hjust = 1, face = "bold"),
    axis.text.y      = element_text(size = 160, face = "bold"),
    strip.text       = element_text(size = 160, face = "bold"),
    strip.background = element_rect(fill = "grey90", colour = NA),
    panel.border     = element_rect()
  )

ggsave(
  filename = file.path(out_dir, "Alpha_diversity_fertilizerlevels_rhizosphere_greenhouse.png"),
  plot     = richness.rhizo,
  device   = "png",
  width    = 32,
  height   = 34,
  units    = "in",
  dpi      = 600,
  bg       = "white"
)




# Non-metric Multidimensional Scaling or NMDS analysis using Bray-Curtis distance of soybean rhizosphere bacterial communities 
# across fertilizer treatments in the greenhouse experiment

# Set seed for reproducibility
set.seed(777)

# Convert to relative abundance
ps.prop.greenhouse <- transform_sample_counts(ps, function(otu) otu / sum(otu))

# Bray-Curtis NMDS
ord.nmds.greenhouse <- ordinate(
  ps.prop.greenhouse,
  method = "NMDS",
  distance = "bray"
)

# Fertilizer colors
custom_colors <- c(
  "No Fertilizer" = "#8B1C62",
  "Medium Fertilizer" = "#00008B",
  "High Fertilizer" = "#009E73"
)

# Plot settings
POINT_SIZE <- 45
POINT_STROKE <- 3.5
JITTER_W <- 0.03
JITTER_H <- 0.03
AXIS_TITLE <- 130
AXIS_TEXT <- 140
FIG_W <- 38
FIG_H <- 25
FIG_DPI <- 600

# NMDS coordinates
ord_df <- plot_ordination(
  ps.prop.greenhouse,
  ord.nmds.greenhouse,
  color = "Fertilizer",
  justDF = TRUE
)

# NMDS plot
greenhouse_nmds <- ggplot(
  ord_df,
  aes(x = NMDS1, y = NMDS2, color = Fertilizer)
) +
  theme_cowplot() +
  geom_point(
    size = POINT_SIZE,
    stroke = POINT_STROKE,
    position = position_jitter(width = JITTER_W, height = JITTER_H)
  ) +
  scale_color_manual(values = custom_colors, limits = fert_order, drop = FALSE) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme(
    axis.line = element_line(linewidth = 2),
    legend.position = "right",
    axis.title = element_text(size = AXIS_TITLE, face = "bold"),
    axis.text = element_text(size = AXIS_TEXT, face = "bold")
  )

ggsave(
  filename = file.path(out_dir, "NMDS_rhizosphere_microbiome_across_fertilizerlevels_greenhouse.png"),
  plot = greenhouse_nmds,
  width = FIG_W,
  height = FIG_H,
  units = "in",
  dpi = FIG_DPI,
  bg = "white"
)

# PERMANOVA and pairwise PERMANOVA comparisons of soybean rhizosphere bacterial communities across fertilizer treatments
# in the greenhouse experiment

bray_greenhouse <- phyloseq::distance(ps, method = "bray")
meta_greenhouse <- data.frame(sample_data(ps))

# Check dispersion across fertilizer treatments
set.seed(777)
disp_greenhouse <- betadisper(bray_greenhouse, meta_greenhouse$Fertilizer)
dispersion_greenhouse <- anova(disp_greenhouse)

print(dispersion_greenhouse)

# PERMANOVA across fertilizer treatments
set.seed(777)
permanova_greenhouse <- adonis2(
  bray_greenhouse ~ Fertilizer,
  data = meta_greenhouse,
  permutations = 999
)

print(permanova_greenhouse)

# Pairwise PERMANOVA across fertilizer treatments
set.seed(777)
pairwise_greenhouse <- pairwise.adonis2(
  bray_greenhouse ~ Fertilizer,
  data = meta_greenhouse,
  permutations = 999,
  p.adjust.m = "bonferroni"
)

print(pairwise_greenhouse)

# Save PERMANOVA and pairwise PERMANOVA results
pairwise_df <- do.call(rbind, lapply(names(pairwise_greenhouse), function(pair) {
  x <- pairwise_greenhouse[[pair]]
  if (!is.data.frame(x)) return(NULL)
  cbind(Comparison = pair, x)
}))

results <- c(
  list(Fertilizer_pairwise_PERMANOVA = pairwise_df),
  list(Fertilizer_PERMANOVA = permanova_greenhouse)
)

write.xlsx(
  results,
  file.path(out_dir, "PERMANOVA_and_pairwise_PERMANOVA_rhizosphere_microbiome_across_fertilizerlevels_greenhouse.xlsx"),
  overwrite = TRUE
)




# Stacked barplot showing relative abundance of rhizosphere bacterial genera (top 40) across the fertilizer treatments 
#and soybean genotypes in the greenhosue experiment


# Subset rhizosphere samples from the greenhouse experiment
Jay.rhizo <- subset_samples(greenhouse.physeq, Location == "Jay")
Jay.rhizo <- prune_samples(sample_sums(Jay.rhizo) > 0, Jay.rhizo)
Jay.rhizo <- prune_taxa(taxa_sums(Jay.rhizo) > 0, Jay.rhizo)

# Keep taxa with genus-level classification
ps.genus <- subset_taxa(Jay.rhizo, !is.na(Genus) & Genus != "")

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

# Save plot
ggsave(file.path(out_dir, "Stacked_barplot_top40_genus_across_fertilizerlevels_rhizosphere_greenhouse.png"),
       p, width = 24, height = 12, dpi = 600)

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

ggsave(file.path(out_dir, "Legend_top40_genus_across_fertilizerlevels_rhizosphere_greenhouse.png"),
       plot = legend_only, width = 8, height = 12, units = "in", dpi = 600, bg = "white")





# Differential abundance analysis by DESeq2 of rhizosphere bacterial genera between the no fertilizer and high fertilizer treatments
# in the greenhouse experiment

# # Subset rhizosphere samples from the greenhouse experiment
Jay.rhizo <- subset_samples(greenhouse.physeq, Location == "Jay")
Jay.rhizo <- prune_samples(sample_sums(Jay.rhizo) > 0, Jay.rhizo)
Jay.rhizo <- prune_taxa(taxa_sums(Jay.rhizo) > 0, Jay.rhizo)

# Keep None and High fertilizer treatments
sample_data(Jay.rhizo)$Fertilizer <- trimws(as.character(sample_data(Jay.rhizo)$Fertilizer))
sample_data(Jay.rhizo)$Fertilizer <- ifelse(sample_data(Jay.rhizo)$Fertilizer %in% c("None", "High"),
                                            sample_data(Jay.rhizo)$Fertilizer, NA)

sample_data(Jay.rhizo)$Fertilizer <- factor(
  sample_data(Jay.rhizo)$Fertilizer,
  levels = c("None", "High")
)

Jay.rhizo <- prune_samples(
  !is.na(sample_data(Jay.rhizo)$Fertilizer),
  Jay.rhizo
)

# Agglomerate taxa at the genus level
ps_genus <- tax_glom(Jay.rhizo, taxrank = "Genus", NArm = TRUE)

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
  file.path(out_dir, "DESeq2_log2FC_and_log10BaseMean_high_vs_no_fertilizer_rhizosphere_greenhouse.png"),
  plot = p_merged, width = 35, height = 45, units = "in", dpi = 300, bg = "white"
)





