# Phenotypic responses of soybean genotypes under iron deficient and Soygreen treatments 
# in the greenhouse treatments

# Load required packages
library(tidyverse)
library(rstatix)
library(ggpubr)
library(cowplot)

# Output directory
out_dir <- "Greenhouse_Phenotype"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Colors
variety_colors <- c(
  "Sensitive" = "blue",
  "Resistant" = "red"
)

# Phenotypic data

datasets <- list(
  
# Nitrate/bicarbonate SPAD
  
  "NitrateBicarb_SPAD" = list(
    
    data = tribble(
      ~Treatment, ~Genotype, ~r1, ~r2, ~r3, ~r4, ~r5,
      "5Mm Nitrate 0mM Bicarbonate", "ND17009GT", 11.7, 14.3, 10.8, 13.4, 17.3,
      "5Mm Nitrate 0mM Bicarbonate", "A11", 25.8, 21.9, 25.2, 30.2, 28.3,
      "5Mm Nitrate 5mM Bicarbonate", "ND17009GT", 12, 13.2, 14.3, 15.2, 13.4,
      "5Mm Nitrate 5mM Bicarbonate", "A11", 24.7, 22.9, 21.8, 24.6, 26.4,
      "10 Mm Nitrate 0mM Bicarbonate", "ND17009GT", 16.2, 12, 11.5, 16.8, 13.1,
      "10 Mm Nitrate 0mM Bicarbonate", "A11", 23, 25, 22.8, 28.6, 22.6,
      "10 Mm Nitrate 5mM Bicarbonate", "ND17009GT", 8.6, 6.9, 7.7, 5.7, 8.3,
      "10 Mm Nitrate 5mM Bicarbonate", "A11", 25.1, 24.6, 24.3, 25.3, 26
    ),
    
    treatment_levels = c(
      "5Mm Nitrate 0mM Bicarbonate",
      "5Mm Nitrate 5mM Bicarbonate",
      "10 Mm Nitrate 0mM Bicarbonate",
      "10 Mm Nitrate 5mM Bicarbonate"
    ),
    
    treatment_labels = c(
      "5mM Nitrate\n0mM Bicarbonate",
      "5mM Nitrate\n5mM Bicarbonate",
      "10mM Nitrate\n0mM Bicarbonate",
      "10mM Nitrate\n5mM Bicarbonate"
    ),
    
    comparisons = "all",
    
    target_pairs = NULL,
    
    y_label = "SPAD Scores",
    
    step = 10,
    headroom = 30,
    tip_length = 0.01,
    signif_size = 70,
    bracket_width = 5,
    
    axis_line = 5,
    axis_title = 120,
    axis_x = 80,
    axis_y = 120,
    x_angle = 0,
    
    width = 40,
    height = 35,
    dpi = 300,
    
    file_name = "SPAD_NitrateBicarb_Grouped.png"
  ),
  
# Nitrate/bicarbonate shoot dry weight
  
  "NitrateBicarb_SDW" = list(
    
    data = tribble(
      ~Treatment, ~Genotype, ~r1, ~r2, ~r3, ~r4, ~r5,
      "5Mm Nitrate 0mM Bicarbonate", "A11", 543.1, 476.3, 446.1, 675.0, 717.5,
      "5Mm Nitrate 0mM Bicarbonate", "ND17009GT", 370.0, 269.7, 325.4, 438.2, 352.9,
      "5Mm Nitrate 5mM Bicarbonate", "A11", 777.3, 792.2, 545.2, 657.8, 558.4,
      "5Mm Nitrate 5mM Bicarbonate", "ND17009GT", 315.6, 332.4, 226.2, 157.4, 242.9,
      "10 Mm Nitrate 0mM Bicarbonate", "A11", 716.4, 684.4, 754.8, 672.1, 714.4,
      "10 Mm Nitrate 0mM Bicarbonate", "ND17009GT", 258.7, 234.7, 382.3, 435.8, 321.6,
      "10 Mm Nitrate 5mM Bicarbonate", "A11", 826.1, 676.1, 737.6, 611.3, 655.5,
      "10 Mm Nitrate 5mM Bicarbonate", "ND17009GT", 316.9, 238.4, 201.5, 146.9, 339.8
    ),
    
    treatment_levels = c(
      "5Mm Nitrate 0mM Bicarbonate",
      "5Mm Nitrate 5mM Bicarbonate",
      "10 Mm Nitrate 0mM Bicarbonate",
      "10 Mm Nitrate 5mM Bicarbonate"
    ),
    
    treatment_labels = c(
      "5mM Nitrate\n0mM Bicarbonate",
      "5mM Nitrate\n5mM Bicarbonate",
      "10mM Nitrate\n0mM Bicarbonate",
      "10mM Nitrate\n5mM Bicarbonate"
    ),
    
    comparisons = "selected",
    
    target_pairs = list(
      c("5mM Nitrate\n0mM Bicarbonate", "10mM Nitrate\n5mM Bicarbonate"),
      c("5mM Nitrate\n5mM Bicarbonate", "10mM Nitrate\n5mM Bicarbonate"),
      c("10mM Nitrate\n0mM Bicarbonate", "10mM Nitrate\n5mM Bicarbonate")
    ),
    
    y_label = "Shoot Dry Weight",
    
    step = 80,
    headroom = 200,
    tip_length = 0.02,
    signif_size = 70,
    bracket_width = 5,
    
    axis_line = 2,
    axis_title = 100,
    axis_x = 60,
    axis_y = 100,
    x_angle = 0,
    
    width = 40,
    height = 35,
    dpi = 600,
    
    file_name = "SDW_NitrateBicarb_Grouped.png"
  ),
  
# Soygreen shoot dry weight
  
  "Soygreen_SDW" = list(
    
    data = tribble(
      ~Treatment, ~Genotype, ~r1, ~r2, ~r3, ~r4, ~r5,
      "0ug soygreen", "A11", 501.0, 707.3, 535.3, 536.7, 683.4,
      "0ug soygreen", "ND17009GT", 178.5, 205.0, 330.2, 353.3, 240.0,
      "25ug soygreen", "A11", 660.0, 680.4, 728.5, 872.4, 983.2,
      "25ug soygreen", "ND17009GT", 718.6, 578.3, 585.0, 650.0, 825.3,
      "50ug soygreen", "A11", 1205.0, 1110.1, 1186.1, 992.5, 940.2,
      "50ug soygreen", "ND17009GT", 953.7, 980.1, 853.2, 881.5, 785.0
    ),
    
    treatment_levels = c(
      "0ug soygreen",
      "25ug soygreen",
      "50ug soygreen"
    ),
    
    treatment_labels = c(
      "No Fertilizer",
      "Medium Fertilizer",
      "High Fertilizer"
    ),
    
    comparisons = "selected",
    
    target_pairs = list(
      c("No Fertilizer", "Medium Fertilizer"),
      c("No Fertilizer", "High Fertilizer")
    ),
    
    y_label = "Shoot Dry Weight",
    
    step = 80,
    headroom = 140,
    tip_length = 0.02,
    signif_size = 50,
    bracket_width = 5,
    
    axis_line = 2,
    axis_title = 100,
    axis_x = 100,
    axis_y = 100,
    x_angle = 45,
    
    width = 28,
    height = 30,
    dpi = 600,
    
    file_name = "SDW_Soygreen_Grouped.png"
  ),
  
# Soygreen SPAD
  
  "Soygreen_SPAD" = list(
    
    data = tribble(
      ~Treatment, ~Genotype, ~r1, ~r2, ~r3, ~r4, ~r5,
      "0ug soygreen", "ND17009GT", 10.0, 14.0, 16.3, 10.2, 15.0,
      "0ug soygreen", "A11", 25.9, 29.2, 30.0, 26.0, 27.4,
      "25ug soygreen", "ND17009GT", 42.5, 33.0, 37.1, 38.8, 36.6,
      "25ug soygreen", "A11", 34.6, 32.7, 33.0, 43.1, 38.7,
      "50ug soygreen", "ND17009GT", 43.1, 38.7, 45.0, 40.6, 42.8,
      "50ug soygreen", "A11", 42.5, 41.9, 36.6, 40.3, 38.0
    ),
    
    treatment_levels = c(
      "0ug soygreen",
      "25ug soygreen",
      "50ug soygreen"
    ),
    
    treatment_labels = c(
      "No Fertilizer",
      "Medium Fertilizer",
      "High Fertilizer"
    ),
    
    comparisons = "selected",
    
    target_pairs = list(
      c("No Fertilizer", "Medium Fertilizer"),
      c("No Fertilizer", "High Fertilizer")
    ),
    
    y_label = "SPAD Scores",
    
    step = 10,
    headroom = 30,
    tip_length = 0.02,
    signif_size = 50,
    bracket_width = 5,
    
    axis_line = 2,
    axis_title = 100,
    axis_x = 100,
    axis_y = 100,
    x_angle = 45,
    
    width = 28,
    height = 30,
    dpi = 600,
    
    file_name = "SPAD_Soygreen_Grouped.png"
  )
)

# Shapiro-Wilk normality test
sw_all <- do.call(rbind, lapply(names(datasets), function(name) {
  
  info <- datasets[[name]]
  
  dat_long <- info$data %>%
    pivot_longer(starts_with("r"), names_to = "Replicate", values_to = "Value") %>%
    mutate(
      Variety = case_when(
        Genotype == "ND17009GT" ~ "Sensitive",
        Genotype == "A11" ~ "Resistant",
        TRUE ~ as.character(Genotype)
      ),
      Treatment = factor(
        Treatment,
        levels = info$treatment_levels,
        labels = info$treatment_labels
      )
    ) %>%
    filter(Variety == "Sensitive") %>%
    drop_na(Value)
  
  do.call(rbind, lapply(info$treatment_labels, function(g) {
    
    vals <- dat_long$Value[dat_long$Treatment == g]
    
    if (length(vals) < 3) {
      return(data.frame(
        Dataset = name,
        Treatment = g,
        n = length(vals),
        W = NA,
        p = NA,
        Normal = "SKIP - too few"
      ))
    }
    
    if (length(unique(vals)) < 2) {
      return(data.frame(
        Dataset = name,
        Treatment = g,
        n = length(vals),
        W = NA,
        p = NA,
        Normal = "SKIP - all identical"
      ))
    }
    
    sw <- shapiro.test(vals)
    
    data.frame(
      Dataset = name,
      Treatment = g,
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
run_phenotype_plot <- function(dataset_name) {
  
  info <- datasets[[dataset_name]]
  
  dat_long <- info$data %>%
    pivot_longer(starts_with("r"), names_to = "Replicate", values_to = "Value") %>%
    mutate(
      Variety = case_when(
        Genotype == "ND17009GT" ~ "Sensitive",
        Genotype == "A11" ~ "Resistant",
        TRUE ~ as.character(Genotype)
      ),
      Variety = factor(Variety, levels = c("Sensitive", "Resistant")),
      Treatment = factor(
        Treatment,
        levels = info$treatment_levels,
        labels = info$treatment_labels
      )
    ) %>%
    drop_na(Value)
  
# Summary statistics
sum_df <- dat_long %>%
    group_by(Treatment, Variety) %>%
    summarise(
      mean = mean(Value, na.rm = TRUE),
      se = sd(Value, na.rm = TRUE) / sqrt(sum(!is.na(Value))),
      .groups = "drop"
    )
  
# Sensitive genotype statistical analysis
one_df <- dat_long %>%
    filter(Variety == "Sensitive")
  
  normality <- sw_all %>%
    filter(Dataset == dataset_name)
  
  all_normal <- all(normality$Normal == "YES")
  
# Skip inferential statistics when all values are identical
if (length(unique(one_df$Value)) < 2) {
    
    cat("\n", dataset_name, "- No test: all values identical\n")
    
    posthoc_result <- data.frame()
    
# ANOVA and Tukey HSD when all groups meet normality
} else if (all_normal) {
    
    anova_result <- one_df %>%
      anova_test(Value ~ Treatment)
    
    tukey_result <- one_df %>%
      tukey_hsd(Value ~ Treatment)
    
    cat("\n", dataset_name, "- One-way ANOVA\n")
    print(anova_result)
    
    cat("\nTukey HSD\n")
    print(tukey_result)
    
    posthoc_result <- tukey_result %>%
      mutate(plot_signif = ifelse(p.adj.signif == "****", "***", p.adj.signif))
    
# Kruskal-Wallis and pairwise Wilcoxon when normality is not met
} else {
    
    kw_result <- one_df %>%
      kruskal_test(Value ~ Treatment)
    
    wilcox_result <- one_df %>%
      wilcox_test(Value ~ Treatment, p.adjust.method = "BH", paired = FALSE) %>%
      add_significance("p.adj")
    
    cat("\n", dataset_name, "- Kruskal-Wallis\n")
    print(kw_result)
    
    cat("\nPairwise Wilcoxon with BH correction\n")
    print(wilcox_result)
    
    posthoc_result <- wilcox_result %>%
      mutate(plot_signif = ifelse(p.adj.signif == "****", "***", p.adj.signif))
    
  }
  
# Select significant comparisons for the plot
if (nrow(posthoc_result) == 0) {
    
    stat.test <- data.frame()
    
  } else if (info$comparisons == "all") {
    
    stat.test <- posthoc_result %>%
      filter(p.adj <= 0.05)
    
  } else {
    
    stat.test <- posthoc_result %>%
      filter(p.adj <= 0.05) %>%
      filter(map2_lgl(group1, group2, function(g1, g2) {
        
        any(vapply(info$target_pairs, function(pr) {
          (g1 == pr[1] & g2 == pr[2]) |
            (g1 == pr[2] & g2 == pr[1])
        }, logical(1)))
        
      }))
    
  }
  
# Significance bracket positions
y_max <- max(sum_df$mean + sum_df$se, na.rm = TRUE)
  
  if (nrow(stat.test) > 0) {
    stat.test <- stat.test %>%
      arrange(p.adj)
    
    stat.test$y.position <- y_max + info$step * seq_len(nrow(stat.test))
    plot_ymax <- max(stat.test$y.position) + info$headroom
  } else {
    plot_ymax <- y_max * 1.15
  }
  
  # Plot
p <- ggplot(sum_df, aes(x = Treatment, y = mean, fill = Variety, color = Variety)) +
    geom_col(width = 0.62, alpha = 0.55, linewidth = 5, position = position_dodge(width = 0.75)) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 5, position = position_dodge(width = 0.75)) +
    geom_point(data = dat_long, aes(x = Treatment, y = Value, color = Variety), position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.75), size = 7, alpha = 0.95) +
    {
      if (nrow(stat.test) > 0) {
        stat_pvalue_manual(
          stat.test,
          label = "plot_signif",
          xmin = "group1",
          xmax = "group2",
          y.position = "y.position",
          tip.length = info$tip_length,
          size = info$signif_size,
          linewidth = info$bracket_width,
          inherit.aes = FALSE
        )
      }
    } +
    scale_fill_manual(values = variety_colors) +
    scale_color_manual(values = variety_colors) +
    labs(x = "", y = info$y_label, fill = NULL, color = NULL) +
    coord_cartesian(ylim = c(0, plot_ymax)) +
    theme_classic(base_size = 70) +
    theme(
      axis.line = element_line(linewidth = info$axis_line),
      text = element_text(face = "bold"),
      axis.title = element_text(size = info$axis_title),
      axis.text.x = element_text(size = info$axis_x, angle = info$x_angle, hjust = ifelse(info$x_angle == 45, 1, 0.5), colour = "black"),
      axis.text.y = element_text(size = info$axis_y, colour = "black"),
      legend.position = "none"
    )
  
  print(p)
  
  # Save figure
ggsave(
    file.path(out_dir, info$file_name),
    p,
    width = info$width,
    height = info$height,
    dpi = info$dpi,
    bg = "white"
  )
  
}

# Generate all phenotype figures
for (name in names(datasets)) {
  run_phenotype_plot(name)
}

# Genotype legend
legend_plot <- ggplot(
  data.frame(
    Variety = factor(
      c("Sensitive", "Resistant"),
      levels = c("Sensitive", "Resistant")
    )
  ),
  aes(x = Variety, y = 1, fill = Variety, color = Variety)
) +
  geom_col() +
  scale_fill_manual(values = variety_colors) +
  scale_color_manual(values = variety_colors) +
  labs(fill = NULL, color = NULL) +
  theme_classic(base_size = 70) +
  theme(
    text = element_text(face = "bold"),
    legend.text = element_text(size = 100, colour = "black"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.size = unit(3.5, "cm"),
    legend.key.width = unit(3.5, "cm"),
    legend.key.height = unit(3.5, "cm"),
    legend.spacing.y = unit(1.2, "cm")
  ) +
  guides(
    fill = guide_legend(ncol = 1, byrow = TRUE),
    color = guide_legend(ncol = 1, byrow = TRUE)
  )

leg <- cowplot::get_legend(legend_plot)

cowplot::save_plot(
  file.path(out_dir, "Soygreen_Legend.png"),
  leg,
  base_width = 20,
  base_height = 20,
  dpi = 600
)
