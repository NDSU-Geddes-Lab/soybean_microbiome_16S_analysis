# Soybean bacterial community analysis across different fertilizer treatments 
# at Leonard and Casselton in 2023 

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
out_dir <- "bacterial_community_analysis_Leonard_Casselton_2023"
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


# Shannon diversity analysis across fertilizer treatments within each plant compartment at Leonard and Casselton in 2023
# Subset by plant compartments across the fertilizer levels

IDCY3.all <- subset_samples(ps.rare, Year == "Y2023")
Leonard.IDCY3 <- subset_samples(IDCY3.all, Location == "Leonard")
Leonard.IDCY3.Bulk.fert  <- subset_samples(Leonard.IDCY3, Compartments == "Bulk")
Leonard.IDCY3.Rhizo.fert <- subset_samples(Leonard.IDCY3, Compartments == "Rhizosphere")
Leonard.IDCY3.Endo.fert  <- subset_samples(Leonard.IDCY3, Compartments == "Endosphere")

Casselton.IDCY3 <- subset_samples(IDCY3.all, Location == "Casselton")
Casselton.IDCY3.Bulk.fert  <- subset_samples(Casselton.IDCY3, Compartments == "Bulk")
Casselton.IDCY3.Rhizo.fert <- subset_samples(Casselton.IDCY3, Compartments == "Rhizosphere")
Casselton.IDCY3.Endo.fert  <- subset_samples(Casselton.IDCY3, Compartments == "Endosphere")


# Casselton 2023 - Bulk soil
ps <- Casselton.IDCY3.Bulk.fert
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
sw_cass_bulk <- do.call(rbind, lapply(fert_order, function(g) {
  sub <- alpha_df[alpha_df$Fertilizer == g, ]
  sw  <- shapiro.test(sub$Shannon)
  data.frame(
    Location = "Casselton", Compartment = "Bulk", Fertilizer = g,
    n = nrow(sub), W = round(sw$statistic, 4), p = round(sw$p.value, 4),
    Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
}))

print(sw_cass_bulk)

# Normality met - one-way ANOVA + Tukey HSD
#one-way ANOVA
aov_cass_bulk <- aov(Shannon ~ Fertilizer, data = alpha_df)

print(summary(aov_cass_bulk))

# Tukey HSD post-hoc test
tukey_cass_bulk_all <- alpha_df %>%
  tukey_hsd(Shannon ~ Fertilizer) %>%
  add_significance("p.adj")

print(as.data.frame(tukey_cass_bulk_all))

aov_cass_bulk_export <- data.frame(
  Location     = "Casselton",
  Compartment  = "Bulk",
  test         = "One-way ANOVA",
  group1       = NA_character_,
  group2       = NA_character_,
  statistic    = summary(aov_cass_bulk)[[1]]$`F value`[1],
  p            = summary(aov_cass_bulk)[[1]]$`Pr(>F)`[1],
  p.adj        = NA_real_,
  p.adj.signif = NA_character_
)

tukey_cass_bulk_export <- as.data.frame(tukey_cass_bulk_all) %>%
  mutate(
    Location    = "Casselton",
    Compartment = "Bulk",
    test        = "Tukey HSD",
    statistic   = NA_real_,
    p           = NA_real_
  ) %>%
  select(Location, Compartment, test, group1, group2, statistic, p, p.adj, p.adj.signif)


# Casselton 2023 - Rhizosphere
ps <- Casselton.IDCY3.Rhizo.fert
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
sw_cass_rhiz <- do.call(rbind, lapply(fert_order, function(g) {
  sub <- alpha_df[alpha_df$Fertilizer == g, ]
  sw  <- shapiro.test(sub$Shannon)
  data.frame(
    Location = "Casselton", Compartment = "Rhizosphere", Fertilizer = g,
    n = nrow(sub), W = round(sw$statistic, 4), p = round(sw$p.value, 4),
    Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
}))

print(sw_cass_rhiz)

# Normality not met - Kruskal-Wallis + pairwise Wilcoxon with BH correction
# Kruskal-Wallis test
kw_cass_rhiz <- alpha_df %>%
  kruskal_test(Shannon ~ Fertilizer)

print(as.data.frame(kw_cass_rhiz))

# Pairwise Wilcoxon with BH correction
pwc_cass_rhiz_all <- alpha_df %>%
  wilcox_test(Shannon ~ Fertilizer, p.adjust.method = "BH", paired = FALSE) %>%
  add_significance("p.adj")

print(as.data.frame(pwc_cass_rhiz_all))

kw_cass_rhiz_export <- as.data.frame(kw_cass_rhiz) %>%
  mutate(
    Location     = "Casselton",
    Compartment  = "Rhizosphere",
    test         = "Kruskal-Wallis",
    group1       = NA_character_,
    group2       = NA_character_,
    p.adj        = NA_real_,
    p.adj.signif = NA_character_
  ) %>%
  select(Location, Compartment, test, group1, group2, statistic, p, p.adj, p.adj.signif)

pwc_cass_rhiz_export <- as.data.frame(pwc_cass_rhiz_all) %>%
  mutate(
    Location    = "Casselton",
    Compartment = "Rhizosphere",
    test        = "Pairwise Wilcoxon (BH adjusted)"
  ) %>%
  select(Location, Compartment, test, group1, group2, statistic, p, p.adj, p.adj.signif)


# Casselton 2023 - Endosphere
ps <- Casselton.IDCY3.Endo.fert
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
sw_cass_endo <- do.call(rbind, lapply(fert_order, function(g) {
  sub <- alpha_df[alpha_df$Fertilizer == g, ]
  sw  <- shapiro.test(sub$Shannon)
  data.frame(
    Location = "Casselton", Compartment = "Endosphere", Fertilizer = g,
    n = nrow(sub), W = round(sw$statistic, 4), p = round(sw$p.value, 4),
    Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
}))

print(sw_cass_endo)

# Normality met - one-way ANOVA + Tukey HSD
# one-way ANOVA
aov_cass_endo <- aov(Shannon ~ Fertilizer, data = alpha_df)

print(summary(aov_cass_endo))

# Tukey HSD post-hoc test
tukey_cass_endo_all <- alpha_df %>%
  tukey_hsd(Shannon ~ Fertilizer) %>%
  add_significance("p.adj")

print(as.data.frame(tukey_cass_endo_all))

aov_cass_endo_export <- data.frame(
  Location     = "Casselton",
  Compartment  = "Endosphere",
  test         = "One-way ANOVA",
  group1       = NA_character_,
  group2       = NA_character_,
  statistic    = summary(aov_cass_endo)[[1]]$`F value`[1],
  p            = summary(aov_cass_endo)[[1]]$`Pr(>F)`[1],
  p.adj        = NA_real_,
  p.adj.signif = NA_character_
)

tukey_cass_endo_export <- as.data.frame(tukey_cass_endo_all) %>%
  mutate(
    Location    = "Casselton",
    Compartment = "Endosphere",
    test        = "Tukey HSD",
    statistic   = NA_real_,
    p           = NA_real_
  ) %>%
  select(Location, Compartment, test, group1, group2, statistic, p, p.adj, p.adj.signif)


# Leonard 2023 - Bulk soil
ps <- Leonard.IDCY3.Bulk.fert
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
sw_leon_bulk <- do.call(rbind, lapply(fert_order, function(g) {
  sub <- alpha_df[alpha_df$Fertilizer == g, ]
  sw  <- shapiro.test(sub$Shannon)
  data.frame(
    Location = "Leonard", Compartment = "Bulk", Fertilizer = g,
    n = nrow(sub), W = round(sw$statistic, 4), p = round(sw$p.value, 4),
    Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
}))

print(sw_leon_bulk)

# Normality met - one-way ANOVA + Tukey HSD
# one-way ANOVA
aov_leon_bulk <- aov(Shannon ~ Fertilizer, data = alpha_df)

print(summary(aov_leon_bulk))

# Tukey HSD post-hoc test
tukey_leon_bulk_all <- alpha_df %>%
  tukey_hsd(Shannon ~ Fertilizer) %>%
  add_significance("p.adj")

print(as.data.frame(tukey_leon_bulk_all))

aov_leon_bulk_export <- data.frame(
  Location     = "Leonard",
  Compartment  = "Bulk",
  test         = "One-way ANOVA",
  group1       = NA_character_,
  group2       = NA_character_,
  statistic    = summary(aov_leon_bulk)[[1]]$`F value`[1],
  p            = summary(aov_leon_bulk)[[1]]$`Pr(>F)`[1],
  p.adj        = NA_real_,
  p.adj.signif = NA_character_
)

tukey_leon_bulk_export <- as.data.frame(tukey_leon_bulk_all) %>%
  mutate(
    Location    = "Leonard",
    Compartment = "Bulk",
    test        = "Tukey HSD",
    statistic   = NA_real_,
    p           = NA_real_
  ) %>%
  select(Location, Compartment, test, group1, group2, statistic, p, p.adj, p.adj.signif)


# Leonard 2023 - Rhizosphere
ps <- Leonard.IDCY3.Rhizo.fert
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
sw_leon_rhiz <- do.call(rbind, lapply(fert_order, function(g) {
  sub <- alpha_df[alpha_df$Fertilizer == g, ]
  sw  <- shapiro.test(sub$Shannon)
  data.frame(
    Location = "Leonard", Compartment = "Rhizosphere", Fertilizer = g,
    n = nrow(sub), W = round(sw$statistic, 4), p = round(sw$p.value, 4),
    Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
}))

print(sw_leon_rhiz)

# Normality met - one-way ANOVA + Tukey HSD
# one-way ANOVA
aov_leon_rhiz <- aov(Shannon ~ Fertilizer, data = alpha_df)

print(summary(aov_leon_rhiz))

# Tukey HSD post-hoc test
tukey_leon_rhiz_all <- alpha_df %>%
  tukey_hsd(Shannon ~ Fertilizer) %>%
  add_significance("p.adj")

print(as.data.frame(tukey_leon_rhiz_all))

aov_leon_rhiz_export <- data.frame(
  Location     = "Leonard",
  Compartment  = "Rhizosphere",
  test         = "One-way ANOVA",
  group1       = NA_character_,
  group2       = NA_character_,
  statistic    = summary(aov_leon_rhiz)[[1]]$`F value`[1],
  p            = summary(aov_leon_rhiz)[[1]]$`Pr(>F)`[1],
  p.adj        = NA_real_,
  p.adj.signif = NA_character_
)

tukey_leon_rhiz_export <- as.data.frame(tukey_leon_rhiz_all) %>%
  mutate(
    Location    = "Leonard",
    Compartment = "Rhizosphere",
    test        = "Tukey HSD",
    statistic   = NA_real_,
    p           = NA_real_
  ) %>%
  select(Location, Compartment, test, group1, group2, statistic, p, p.adj, p.adj.signif)


# Leonard 2023 - Endosphere
ps <- Leonard.IDCY3.Endo.fert
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
sw_leon_endo <- do.call(rbind, lapply(fert_order, function(g) {
  sub <- alpha_df[alpha_df$Fertilizer == g, ]
  sw  <- shapiro.test(sub$Shannon)
  data.frame(
    Location = "Leonard", Compartment = "Endosphere", Fertilizer = g,
    n = nrow(sub), W = round(sw$statistic, 4), p = round(sw$p.value, 4),
    Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
}))

print(sw_leon_endo)

# Normality met - one-way ANOVA + Tukey HSD
# one-way ANOVA
aov_leon_endo <- aov(Shannon ~ Fertilizer, data = alpha_df)

print(summary(aov_leon_endo))

# Tukey HSD post-hoc test
tukey_leon_endo_all <- alpha_df %>%
  tukey_hsd(Shannon ~ Fertilizer) %>%
  add_significance("p.adj")

print(as.data.frame(tukey_leon_endo_all))

aov_leon_endo_export <- data.frame(
  Location     = "Leonard",
  Compartment  = "Endosphere",
  test         = "One-way ANOVA",
  group1       = NA_character_,
  group2       = NA_character_,
  statistic    = summary(aov_leon_endo)[[1]]$`F value`[1],
  p            = summary(aov_leon_endo)[[1]]$`Pr(>F)`[1],
  p.adj        = NA_real_,
  p.adj.signif = NA_character_
)

tukey_leon_endo_export <- as.data.frame(tukey_leon_endo_all) %>%
  mutate(
    Location    = "Leonard",
    Compartment = "Endosphere",
    test        = "Tukey HSD",
    statistic   = NA_real_,
    p           = NA_real_
  ) %>%
  select(Location, Compartment, test, group1, group2, statistic, p, p.adj, p.adj.signif)


# Save combined table
write_xlsx(
  bind_rows(
    aov_cass_bulk_export,
    tukey_cass_bulk_export,
    kw_cass_rhiz_export,
    pwc_cass_rhiz_export,
    aov_cass_endo_export,
    tukey_cass_endo_export,
    aov_leon_bulk_export,
    tukey_leon_bulk_export,
    aov_leon_rhiz_export,
    tukey_leon_rhiz_export,
    aov_leon_endo_export,
    tukey_leon_endo_export
  ),
  file.path(
    out_dir,
    "Alpha_diversity_statistics_allcompartments_fertilizerlevels_Casselton_Leonard_2023.xlsx"
  )
)


#Beta diversity analysis
# PERMANOVA and pairwise PERMANOVA comparisons across fertilizer treatments within each compartment
# at Leonard and Casselton in 2023

sites <- c("Leonard", "Casselton")
compartments <- c("Bulk", "Rhizosphere", "Endosphere")

dispersion_list <- list()
permanova_list <- list()
pairwise_list <- list()

for (site in sites) {
  for (comp in compartments) {
    
ps_sub <- subset_samples(ps.rare, Location %in% site & Compartments %in% comp & Year %in% "Y2023")
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
       "PERMANOVA_and_pairwise_PERMANOVA_fertilizer_Leonard_Casselton_2023.xlsx"), overwrite = TRUE)






# Stacked barplot showing relative abundance of rhizosphere bacterial genera (top 40) across the fertilizer treatments 
#and soybean genotypes at Leonard in 2023


# Subset Leonard rhizosphere samples from 2023
IDCY3.all <- subset_samples(physeq, Year == "Y2023")
Leonard.IDCY3 <- subset_samples(IDCY3.all, Location == "Leonard")
Leonard.IDCY3.Rhizo <- subset_samples(Leonard.IDCY3, Compartments == "Rhizosphere")

# Keep taxa with genus-level classification
ps.genus <- subset_taxa(Leonard.IDCY3.Rhizo, !is.na(Genus) & Genus != "")

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
  "#00B8D4", "#F4511E", "#7CB342", "#CDDC39", "#FFEB3B", "#FFB300", "#FF7043", "#8D6E63", "#BDBDBD",
  "#AD1457", "#9E9D24", "#4527A0", "#283593", "#0277BD", "#00838F", "#00695C", "#2E7D32", "#78909C",
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
ggsave(file.path(out_dir, "Stacked_barplot_top40_genus_across_fertilizerlevels_Leonard_rhizosphere_2023.png"),
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

ggsave(file.path(out_dir, "Legend_top40_genus_across_fertilizerlevels_Leonard_rhizosphere_2023.png"),
       plot = legend_only, width = 8, height = 12, units = "in", dpi = 600, bg = "white")





# Stacked barplot showing relative abundance of endosphere bacterial genera (top 40) across the fertilizer treatments 
#and soybean genotypes at Leonard in 2023

# Subset Leonard endosphere samples from 2023
IDCY3.all <- subset_samples(physeq, Year == "Y2023")
Leonard.IDCY3 <- subset_samples(IDCY3.all, Location == "Leonard")
Leonard.IDCY3.Endo <- subset_samples(Leonard.IDCY3, Compartments == "Endosphere")

# Keep taxa with genus-level classification
ps.genus <- subset_taxa(Leonard.IDCY3.Endo, !is.na(Genus) & Genus != "")

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
  "#BDBDBD", "#00695C", "#4527A0", "#283593", "#0277BD", "#AD1457", "#00838F", "#2E7D32", "#78909C",
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
ggsave(file.path(out_dir, "Stacked_barplot_top40_genus_across_fertilizerlevels_Leonard_endosphere_2023.png"), 
       plot = p, width = 24, height = 12, units = "in", dpi = 600, bg = "white")




# Save genus legend separately
p_legend <- plot_bar(ps.rel_amp, x = "Sample.names", fill = "Genus") +
  scale_fill_manual(values = my_colors) +
  theme_bw(base_size = 26) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 26, face = "bold")
  )

# Extract legend only
legend_only <- cowplot::get_legend(p_legend)


ggsave(file.path(out_dir, "Legend_top40_genus_across_fertilizerlevels_Leonard_endosphere_2023.png"), 
       plot = legend_only, width = 18, height = 12, units = "in", dpi = 600, bg = "white")


