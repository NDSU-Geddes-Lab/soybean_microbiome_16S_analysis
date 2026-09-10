# IDC scores across soybean genotypes at Leonard under no fertilizer conditions 
# from 2021 to 2023

# Load required packages
library(tidyverse)
library(rstatix)
library(ggpubr)

# Output directory
outdir <- "Leonard_IDC_scores_across_genotypes_years"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Genotype order
geno_order <- c("A11_Early", "ND_Rolette", "ND17009GT", "ND18_17323")

# 2021 IDC scores

dat_wide_2021 <- tribble(
  ~A11_Early, ~ND_Rolette, ~ND17009GT, ~ND18_17323,
  1, 1.5, 3, 3,
  1, 1.5, 3, 2.5,
  1, 1.5, 3, 2.5,
  1, 1.5, 3, 2.5
)

# 2022 IDC scores

dat_wide_2022 <- tribble(
  ~A11_Early, ~ND_Rolette, ~ND17009GT, ~ND18_17323,
  1.5, 4, 5, 4.5,
  2, 3, 5, 5,
  1.5, 4, 5, 5,
  1, 4, 5, 4.5
)

# 2023 IDC scores

dat_wide_2023 <- tribble(
  ~A11_Early, ~ND_Rolette, ~ND17009GT, ~ND18_17323,
  1, 2, 4.5, 3,
  1, 2, 4.5, 3.5,
  1.5, 1, 4, 3,
  1.5, 1, 4.5, 1.5
)

# Statistical analysis and IDC plots
run_idc_plot <- function(dat_wide, year_label, fill_color, file_name) {
  
  dat_long <- dat_wide %>%
    pivot_longer(cols = everything(), names_to = "Genotype", values_to = "Score") %>%
    drop_na(Score) %>%
    mutate(Genotype = factor(Genotype, levels = geno_order))
  
# Shapiro-Wilk normality test
sw_results <- do.call(rbind, lapply(geno_order, function(g) {
    
    vals <- dat_long$Score[dat_long$Genotype == g]
    
    if (length(vals) < 3) {
      return(data.frame(
        Genotype = g,
        n = length(vals),
        W = NA,
        p = NA,
        Normal = "SKIP - too few"
      ))
    }
    
    if (length(unique(vals)) < 2) {
      return(data.frame(
        Genotype = g,
        n = length(vals),
        W = NA,
        p = NA,
        Normal = "SKIP - all identical"
      ))
    }
    
    sw <- shapiro.test(vals)
    
    data.frame(
      Genotype = g,
      n = length(vals),
      W = round(as.numeric(sw$statistic), 4),
      p = round(sw$p.value, 4),
      Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
    )
  }))
  
  cat("\n", year_label, "- Shapiro-Wilk normality test\n")
  print(sw_results, row.names = FALSE)
  
  all_normal <- all(sw_results$Normal == "YES")
  
# Statistical analysis
if (length(unique(dat_long$Score)) < 2) {
    
    cat("\n", year_label, "- No test: all values identical\n")
    
    stat.test <- data.frame()
    
# ANOVA and Tukey HSD when all groups meet normality
} else if (all_normal) {
    
    anova_result <- dat_long %>%
      anova_test(Score ~ Genotype)
    
    tukey_result <- dat_long %>%
      tukey_hsd(Score ~ Genotype)
    
    cat("\n", year_label, "- One-way ANOVA\n")
    print(anova_result)
    
    cat("\nTukey HSD\n")
    print(tukey_result)
    
    stat.test <- tukey_result %>%
      filter(
        (group1 == "A11_Early" & group2 == "ND17009GT") |
          (group1 == "ND17009GT" & group2 == "A11_Early")
      ) %>%
      filter(p.adj <= 0.05) %>%
      mutate(plot_signif = ifelse(p.adj.signif == "****", "***", p.adj.signif))
    
# Kruskal-Wallis and pairwise Wilcoxon when normality is not met or cannot be assessed
} else {
    
    kw_result <- dat_long %>%
      kruskal_test(Score ~ Genotype)
    
    comparisons <- combn(geno_order, 2, simplify = FALSE)
    
    wilcox_result <- do.call(rbind, lapply(comparisons, function(comp) {
      
      group1 <- comp[1]
      group2 <- comp[2]
      
      x <- dat_long$Score[dat_long$Genotype == group1]
      y <- dat_long$Score[dat_long$Genotype == group2]
      
      if (length(unique(c(x, y))) < 2) {
        statistic <- NA
        p <- 1
      } else {
        test <- wilcox.test(x, y, paired = FALSE, exact = FALSE)
        statistic <- as.numeric(test$statistic)
        p <- test$p.value
      }
      
      data.frame(
        .y. = "Score",
        group1 = group1,
        group2 = group2,
        n1 = length(x),
        n2 = length(y),
        statistic = statistic,
        p = p
      )
    }))
    
    wilcox_result$p.adj <- p.adjust(wilcox_result$p, method = "BH")
    
    wilcox_result$p.adj.signif <- case_when(
      wilcox_result$p.adj <= 0.001 ~ "***",
      wilcox_result$p.adj <= 0.01 ~ "**",
      wilcox_result$p.adj <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
    
    cat("\n", year_label, "- Kruskal-Wallis\n")
    print(kw_result)
    
    cat("\nPairwise Wilcoxon with BH correction\n")
    print(wilcox_result)
    
    stat.test <- wilcox_result %>%
      filter(
        (group1 == "A11_Early" & group2 == "ND17009GT") |
          (group1 == "ND17009GT" & group2 == "A11_Early")
      ) %>%
      filter(p.adj <= 0.05) %>%
      mutate(plot_signif = p.adj.signif)
    
  }
  
# Significance bracket positions
if (nrow(stat.test) > 0) {
    stat.test <- stat.test %>%
      arrange(p.adj)
    
    stat.test$y.position <- seq(5.3, 6.0, length.out = nrow(stat.test))
  }
  
# Genotype colors
geno_cols <- c(
    "A11_Early" = fill_color,
    "ND_Rolette" = fill_color,
    "ND17009GT" = fill_color,
    "ND18_17323" = fill_color
  )
  
# Plot
p <- ggplot(dat_long, aes(x = Genotype, y = Score, fill = Genotype)) +
    geom_bar(stat = "summary", fun = "mean", width = 0.6, color = "black") +
    geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), width = 0.2) +
  geom_jitter(width = 0.05, alpha = 0.8, size = 8, color = "black") +
    {
      if (nrow(stat.test) > 0) {
        stat_pvalue_manual(
          stat.test,
          label = "plot_signif",
          y.position = "y.position",
          tip.length = 0.02,
          size = 20,
          linewidth = 3.2,
          inherit.aes = FALSE
        )
      }
    } +
    scale_fill_manual(values = geno_cols) +
    scale_y_continuous(breaks = 1:5) +
    coord_cartesian(ylim = c(1, 6.2)) +
    labs(x = "Variety", y = "IDC Scores") +
    theme_classic(base_size = 60) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 60, face = "bold", colour = "black"),
      axis.text.y = element_text(size = 60, face = "bold", colour = "black"),
      axis.title = element_text(size = 60, face = "bold", colour = "black"),
      legend.position = "none"
    )
  
  
# Save figure
ggsave(
    file.path(outdir, file_name),
    p,
    width = 15,
    height = 13,
    dpi = 500,
    bg = "white"
  )
  
}

# 2021
p2021 <- run_idc_plot(
  dat_wide = dat_wide_2021,
  year_label = "2021",
  fill_color = "darkgreen",
  file_name = "Leonard_IDC_scores_across_genotypes_2021.png"
)

# 2022
p2022 <- run_idc_plot(
  dat_wide = dat_wide_2022,
  year_label = "2022",
  fill_color = "#6D712E",
  file_name = "Leonard_IDC_scores_across_genotypes_2022.png"
)

# 2023
p2023 <- run_idc_plot(
  dat_wide = dat_wide_2023,
  year_label = "2023",
  fill_color = "#7EBD01",
  file_name = "Leonard_IDC_scores_across_genotypes_2023.png"
)
