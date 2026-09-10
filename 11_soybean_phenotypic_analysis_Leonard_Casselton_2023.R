# Phenotypic responses across fertilizer treatments at Leonard and Casselton in 2023

# Load required packages
library(tidyverse)
library(rstatix)
library(ggpubr)

# Output directory
outdir <- "Leonard_Casselton_soybean_phenotype_2023"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Fertilizer order and colors
fert_order <- c("No fertilizer", "Medium fertilizer", "High fertilizer")

fert_colors <- c(
  "No fertilizer" = "#8B1C62",
  "Medium fertilizer" = "#00008B",
  "High fertilizer" = "#009E73"
)

# Phenotype data
datasets <- list(
  
# Leonard
  
  "Leonard_ND18-17323_IDC" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    3, 3, 2,
    3.5, 3.5, 1.5,
    3, 2, 2.5,
    2.5, 1.5, 1.5
  ),
  
  "Leonard_ND18-17323_SPAD" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    20.5, 26.5, 30.6,
    16.2, 18.6, 33.8,
    21, 25.2, 28.5,
    19.5, 28, 38.8
  ),
  
  "Leonard_ND18-17323_SDW" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    9063.4, 6764.5, 15766.2,
    2358.8, 7060.3, 12723.3,
    4706.3, 11496.8, 16871.2,
    8782.3, 10166.9, 18136.3
  ),
  
  "Leonard_ND17009GT_IDC" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    5, 3, 1.5,
    4.5, 3, 2.5,
    4, 3.5, 1.5,
    5, 3.5, 2
  ),
  
  "Leonard_ND17009GT_SPAD" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    6.5, 24.2, 39.7,
    11.7, 22.5, 24.1,
    12, 21.8, 46.2,
    14.8, 19.4, 27.7
  ),
  
  "Leonard_ND17009GT_SDW" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    2594.1, 6968.1, 16513.6,
    1870.2, 11593.4, 12677.9,
    1903.2, 14497.3, 15889.2,
    1090.8, 13081.4, 11266.1
  ),
  
  "Leonard_A11_IDC" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    1, 1.5, 1,
    1, 1.5, 1.5,
    1.5, 1.5, 1,
    1.5, 1.5, 1.5
  ),
  
  "Leonard_A11_SPAD" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    34.2, 32.6, 40.9,
    38.6, 41, 38,
    37.3, 39.9, 45.2,
    41.1, 39.5, 38.9
  ),
  
  "Leonard_A11_SDW" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    18525.8, 19830.6, 16527.5,
    18761.5, 18989.9, 25174.8,
    15565.2, 16705.8, 17082.4,
    13678.7, 10463.2, 18203.3
  ),
  
  "Leonard_Rollete_IDC" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    2, 1, 1,
    2, 2.5, 1.5,
    1, 2, 1,
    1, 1, 1.5
  ),
  
  "Leonard_Rollete_SPAD" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    28, 45, 36.7,
    20.1, 34.5, 37,
    39.2, 31.6, 40,
    19, 43, 40.4
  ),
  
  "Leonard_Rollete_SDW" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    15183.2, 17307.2, 15798.6,
    4415.6, 12438.6, 13466.7,
    2638.4, 14515.4, 17798.7,
    11826.4, 19157.2, 10575.4
  ),
  
# Casselton
  
  "Casselton_ND18-17323_IDC" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    1, 1, 1,
    1, 1, 1,
    1, 1, 1,
    1, 1, 1
  ),
  
  "Casselton_ND18-17323_SPAD" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    37.6, 41.2, 38.1,
    39.5, 41.7, 40.6,
    NA, 36.4, 37.2,
    NA, NA, 33.7
  ),
  
  "Casselton_ND18-17323_SDW" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    3926.7, 6341.2, 2457.9,
    2483.6, 6710.8, 3778.2,
    NA, 2238.1, 10377.9,
    NA, NA, 1985.3
  ),
  
  "Casselton_ND17009GT_IDC" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    1, 1, 1,
    1, 1, 1,
    1, 1, 1,
    1, 1, 1
  ),
  
  "Casselton_ND17009GT_SPAD" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    34.9, 42.7, 47.8,
    38.3, 25.1, 43.8,
    31.8, 41.2, 35.1,
    47, 45.8, 36.6
  ),
  
  "Casselton_ND17009GT_SDW" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    2139.4, 1819.8, 7114.2,
    1849.3, 1507.3, 2369.4,
    6426.3, 7982.9, 1646.3,
    4784.3, 4782.3, 1996.8
  ),
  
  "Casselton_A11_IDC" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    1, 1, 1,
    1, 1, 1,
    1, 1, 1,
    1, 1, 1
  ),
  
  "Casselton_A11_SPAD" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    35.4, 38, 35.9,
    37.1, 38.7, 39.3,
    NA, 37.8, 35.8,
    NA, 39.6, NA
  ),
  
  "Casselton_A11_SDW" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    1612.5, 5193.4, 2827.6,
    2137.4, 5185.6, 4330.6,
    NA, 4230.4, 2832.7,
    NA, 6184.3, NA
  ),
  
  "Casselton_Rollete_IDC" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    1, 1, 1,
    1, 1, 1,
    1, 1, 1,
    1, 1, 1
  ),
  
  "Casselton_Rollete_SPAD" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    41.5, 40.9, 41.9,
    42.6, 38.7, 38.9,
    36.6, 42.1, 40,
    NA, 40.9, NA
  ),
  
  "Casselton_Rollete_SDW" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    4599.2, 4716.2, 3300.2,
    7719.8, 1885.6, 9482.3,
    5407.7, 1040.2, 1838.7,
    NA, 2283.6, NA
  )
)

# Shapiro-Wilk normality test
sw_all <- do.call(rbind, lapply(names(datasets), function(name) {
  
  dat_long <- datasets[[name]] %>%
    pivot_longer(cols = everything(), names_to = "Fertilizer", values_to = "Value") %>%
    drop_na(Value) %>%
    mutate(Fertilizer = factor(Fertilizer, levels = fert_order))
  
  do.call(rbind, lapply(fert_order, function(g) {
    
    vals <- dat_long$Value[dat_long$Fertilizer == g]
    
    if (length(vals) < 3) {
      return(data.frame(
        Dataset = name,
        Fertilizer = g,
        n = length(vals),
        W = NA,
        p = NA,
        Normal = "SKIP - too few"
      ))
    }
    
    if (length(unique(vals)) < 2) {
      return(data.frame(
        Dataset = name,
        Fertilizer = g,
        n = length(vals),
        W = NA,
        p = NA,
        Normal = "SKIP - all identical"
      ))
    }
    
    sw <- shapiro.test(vals)
    
    data.frame(
      Dataset = name,
      Fertilizer = g,
      n = length(vals),
      W = round(as.numeric(sw$statistic), 4),
      p = round(sw$p.value, 4),
      Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
    )
  }))
}))

print(sw_all, row.names = FALSE)

# Identify datasets that do not meet normality
failed <- sw_all %>%
  group_by(Dataset) %>%
  summarise(All_Normal = all(Normal == "YES"), .groups = "drop") %>%
  filter(!All_Normal)

cat("\nDatasets with at least one non-normal, identical, or too-small group:\n")
print(failed$Dataset)

# Statistical analysis and phenotype plots
run_phenotype_plot <- function(dat_wide, dataset_name) {
  
  dat_long <- dat_wide %>%
    pivot_longer(cols = everything(), names_to = "Fertilizer", values_to = "Value") %>%
    drop_na(Value) %>%
    mutate(Fertilizer = factor(Fertilizer, levels = fert_order))
  
  normality <- sw_all %>%
    filter(Dataset == dataset_name)
  
  all_normal <- all(normality$Normal == "YES")
  

# Skip inferential statistics when all values are identical
if (length(unique(dat_long$Value)) < 2) {
    
    stat.test <- data.frame()
    
# ANOVA and Tukey HSD when all groups meet normality
} else if (all_normal) {
    
    anova_result <- dat_long %>%
      anova_test(Value ~ Fertilizer)
    
    tukey_result <- dat_long %>%
      tukey_hsd(Value ~ Fertilizer)
    
    stat.test <- tukey_result %>%
      filter(group1 == "No fertilizer", group2 == "High fertilizer") %>%
      filter(p.adj <= 0.05) %>%
      mutate(plot_signif = ifelse(p.adj.signif == "****", "***", p.adj.signif))
    
# Kruskal-Wallis and pairwise Wilcoxon when normality is not met or cannot be assessed
} else {
    
    kw_result <- dat_long %>%
      kruskal_test(Value ~ Fertilizer)
    
    wilcox_result <- dat_long %>%
      wilcox_test(Value ~ Fertilizer, p.adjust.method = "BH", paired = FALSE) %>%
      add_significance("p.adj")
    
    stat.test <- wilcox_result %>%
      filter(group1 == "No fertilizer", group2 == "High fertilizer") %>%
      filter(p.adj <= 0.05) %>%
      mutate(plot_signif = ifelse(p.adj.signif == "****", "***", p.adj.signif))
    
  }
  
# Significance bracket position
y_max <- max(dat_long$Value, na.rm = TRUE)
  y_pad <- 0.25 * y_max
  
  if (nrow(stat.test) > 0) {
    stat.test$y.position <- y_max + 0.25 * y_pad
  }

# Phenotype label
trait <- sub(".*_", "", dataset_name)
  
  y_label <- switch(
    trait,
    IDC = "IDC Scores",
    SPAD = "SPAD Scores",
    SDW = "SDW"
  )
  
# Plot
p <- ggplot(dat_long, aes(x = Fertilizer, y = Value, fill = Fertilizer)) +
    geom_bar(stat = "summary", fun = "mean", width = 0.6, color = "black") +
    geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), width = 0.2) +
  geom_jitter(width = 0.05, alpha = 0.8, size = 15, color = "black") +
    {
      if (nrow(stat.test) > 0)
        stat_pvalue_manual(
          stat.test,
          label = "plot_signif",
          y.position = "y.position",
          tip.length = 0.02,
          size = 70,
          linewidth = 8,
          inherit.aes = FALSE
        )
    } +
    scale_x_discrete(limits = fert_order, drop = FALSE) +
    scale_fill_manual(values = fert_colors, limits = fert_order, drop = FALSE) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
    coord_cartesian(ylim = c(0, y_max + y_pad)) +
    labs(x = "", y = y_label) +
    theme_classic(base_size = 120) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 120, face = "bold", colour = "black"),
      axis.text.y = element_text(size = 100, face = "bold", colour = "black"),
      axis.title = element_text(size = 120, face = "bold", colour = "black"),
      legend.position = "none"
    )
  
# Save figure
trait_filename <- switch(
    trait,
    IDC = "IDCscores",
    SPAD = "SPADscores",
    SDW = "SDW"
  )
  
sample_name <- sub(paste0("_", trait, "$"), "", dataset_name)
  
ggsave(file.path(outdir, paste0("2023_", sample_name, "_", trait_filename, ".png")),
    p,
    width = 20,
    height = 25,
    dpi = 500,
    bg = "white"
  )
}

# Generate all phenotype figures
for (name in names(datasets)) {
  run_phenotype_plot(datasets[[name]], name)
}
