# Phenotypic responses across fertilizer treatments at Leonard and Casselton in 2022

# Load required packages
library(tidyverse)
library(rstatix)
library(ggpubr)

# Output directory
outdir <- "Leonard_Casselton_soybean_phenotype_2022"
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
    4.5, 3.5, 3,
    5, 1.5, 1,
    5, 3, 1.5,
    4.5, 1.5, 1.5
  ),
  
  "Leonard_ND18-17323_SPAD" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    8.5, 35.6, 25.9,
    6.9, 31.3, 38.5,
    11.1, 25.8, 30.9,
    10, 33.3, 37.5
  ),
  
  "Leonard_ND18-17323_SDW" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    764, 3077.8, 5177.5,
    362.2, 3307.8, 5387.6,
    771.4, 2132.5, 5032.4,
    602.1, 4255.7, 6408.7
  ),
  
  "Leonard_ND17009GT_IDC" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    5, 4, 3.5,
    5, 3.5, 3.5,
    5, 3.5, 3,
    5, 3.5, 2.5
  ),
  
  "Leonard_ND17009GT_SPAD" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    7.4, 31.9, 32.9,
    8.9, 37.9, 34.7,
    5.7, 34.2, 33.5,
    6.8, 24.1, 31.4
  ),
  
  "Leonard_ND17009GT_SDW" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    492.1, 2684.5, 3344.9,
    614.2, 4563.5, 5723.2,
    573.6, 2316.7, 4764.7,
    610.7, 2285.3, 3909.5
  ),
  
  "Leonard_A11_IDC" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    1.5, 1.5, 1,
    2, 1.5, 1,
    1.5, 1.5, 1,
    1, 1, 1
  ),
  
  "Leonard_A11_SPAD" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    30.9, 32.9, 35.6,
    22.2, 33, 36.6,
    33.3, 35.8, 37,
    33.8, 36.9, 38.4
  ),
  
  "Leonard_A11_SDW" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    3328.6, 4362.5, 6813.2,
    2292.7, 4431.4, 7533.2,
    4844.9, 4185.7, 3397.5,
    12354.8, 3601.6, 7202.5
  ),
  
  "Leonard_Rollete_IDC" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    4, 2, 1.5,
    3, 1.5, 1,
    4, 1.5, 1,
    4, 3, 1
  ),
  
  "Leonard_Rollete_SPAD" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    29.4, 27.7, 34.8,
    23.6, 28.4, 32.5,
    25.1, 38.1, 38.8,
    20.8, 30.9, 36.7
  ),
  
  "Leonard_Rollete_SDW" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    1772.8, 1904.5, 4427.5,
    1886.4, 2605.8, 8543.2,
    2063.7, 3635.4, 6378.3,
    1799.9, 2175.8, 4534.5
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
    41.7, 39.4, 35.5,
    33.2, 37.4, 44.2,
    40.5, 41.8, 48.3,
    40.3, 40, 44.7
  ),
  
  "Casselton_ND18-17323_SDW" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    14559.6, 6982.3, 9658.7,
    10857.8, 30566.1, 15992.7,
    11783.5, 20157.6, 16659.4,
    12477.3, 11820.5, 4328.7
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
    30.7, 35.9, 35.5,
    36.6, 43.5, 37.1,
    37.2, 33.8, 40.7,
    33.2, 37.4, 36.1
  ),
  
  "Casselton_ND17009GT_SDW" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    14850.3, 14789.8, 11153.6,
    6994.5, 2436.9, 7714.5,
    10318.6, 11081.2, 8663.2,
    12445.8, 7675.4, 6574.8
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
    40, 39.5, 38.3,
    41.9, 42.1, 44.1,
    38.8, 42.8, 33.5,
    38.5, 39.3, 41.6
  ),
  
  "Casselton_A11_SDW" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    18817.4, 28287.5, 7287.5,
    18944.3, 6231.7, 6426.4,
    13525.6, 19368.8, 8753.1,
    2345.8, 21785.3, 9085.8
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
    43.6, 39.1, 34.6,
    42.8, 43.9, 45.1,
    42.2, 37, 38.3,
    37.9, 42, 42.1
  ),
  
  "Casselton_Rollete_SDW" = tribble(
    ~`No fertilizer`, ~`Medium fertilizer`, ~`High fertilizer`,
    8955.9, 12355.8, 20683.5,
    17849, 30963.8, 24957.2,
    11061.3, 3565.8, 3979.6,
    5440, 13229.6, 9433.2
  )
)

# Shapiro-Wilk normality test
sw_all <- do.call(rbind, lapply(names(datasets), function(name) {
  
  dat_long <- datasets[[name]] %>%
    pivot_longer(cols = everything(), names_to = "Fertilizer", values_to = "Value") %>%
    mutate(Fertilizer = factor(Fertilizer, levels = fert_order))
  
  do.call(rbind, lapply(fert_order, function(g) {
    
    vals <- dat_long$Value[dat_long$Fertilizer == g]
    
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
      W = round(sw$statistic, 4),
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

cat("\nDatasets with at least one non-normal or identical group:\n")
print(failed$Dataset)

# Statistical analysis and phenotype plots
run_phenotype_plot <- function(dat_wide, dataset_name) {
  
  dat_long <- dat_wide %>%
    pivot_longer(cols = everything(), names_to = "Fertilizer", values_to = "Value") %>%
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
    
# Kruskal-Wallis and pairwise Wilcoxon when normality is not met
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
ggsave(file.path(outdir, paste0("2022_", sample_name, "_", trait_filename, ".png")),
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
