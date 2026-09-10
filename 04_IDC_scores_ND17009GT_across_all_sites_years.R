# IDC scores of  ND17009GT across the four sites over the years (2021-2023)

# Load required packages
library(tidyverse)
library(rstatix)
library(ggpubr)

# Output directory
out_dir <- "IDC_ND17009GT_across_sites_years"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Location order and colors
loc_order <- c("Prosper", "Casselton", "Colfax", "Leonard")

loc_colors <- c(
  "Prosper" = "orange",
  "Casselton" = "red",
  "Colfax" = "blue",
  "Leonard" = "darkgreen"
)

# IDC scores of ND17009GT
dat_wide <- tribble(
  ~Prosper, ~Casselton, ~Colfax, ~Leonard,
  1, 3, 1, 3.5,
  1, 3, 1.5, 4,
  1, 3, 2, 5,
  1, 3, 2.5, 4.5,
  1, 1, 5, 5,
  1, 1.5, 4.5, 5,
  1, 1.5, 4, 5,
  1, 1, 5, 5,
  NA, 1, NA, 4.5,
  NA, 1.5, NA, 4.5,
  NA, 1, NA, 4,
  NA, 1, NA, 4.5
)

dat_long <- dat_wide %>%
  pivot_longer(cols = everything(), names_to = "Location", values_to = "Score") %>%
  drop_na(Score) %>%
  mutate(Location = factor(Location, levels = loc_order))

# Shapiro-Wilk normality test
sw_results <- do.call(rbind, lapply(loc_order, function(g) {
  
  vals <- dat_long$Score[dat_long$Location == g]
  
  if (length(vals) < 3) {
    return(data.frame(
      Location = g,
      n = length(vals),
      W = NA,
      p = NA,
      Normal = "SKIP - too few"
    ))
  }
  
  if (length(unique(vals)) < 2) {
    return(data.frame(
      Location = g,
      n = length(vals),
      W = NA,
      p = NA,
      Normal = "SKIP - all identical"
    ))
  }
  
sw <- shapiro.test(vals)
  
  data.frame(
    Location = g,
    n = length(vals),
    W = round(as.numeric(sw$statistic), 4),
    p = round(sw$p.value, 4),
    Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
}))

print(sw_results, row.names = FALSE)

all_normal <- all(sw_results$Normal == "YES")

# Statistical analysis
if (length(unique(dat_long$Score)) < 2) {
  
  stat.test <- data.frame()
  
# ANOVA and Tukey HSD when all groups meet normality
} else if (all_normal) {
  
  anova_result <- dat_long %>%
    anova_test(Score ~ Location)
  
  tukey_result <- dat_long %>%
    tukey_hsd(Score ~ Location)
  
  print(anova_result)
  print(tukey_result)
  
  stat.test <- tukey_result %>%
    filter(
      (group1 == "Leonard" & group2 %in% c("Prosper", "Casselton", "Colfax")) |
        (group2 == "Leonard" & group1 %in% c("Prosper", "Casselton", "Colfax"))
    ) %>%
    filter(p.adj <= 0.05) %>%
    mutate(plot_signif = ifelse(p.adj.signif == "****", "***", p.adj.signif))
  
# Kruskal-Wallis and pairwise Wilcoxon when normality is not met or cannot be assessed
} else {
  
  kw_result <- dat_long %>%
    kruskal_test(Score ~ Location)
  
  wilcox_result <- dat_long %>%
    wilcox_test(Score ~ Location, p.adjust.method = "BH", paired = FALSE) %>%
    add_significance("p.adj")
  
  print(kw_result)
  print(wilcox_result)
  
  stat.test <- wilcox_result %>%
    filter(
      (group1 == "Leonard" & group2 %in% c("Prosper", "Casselton", "Colfax")) |
        (group2 == "Leonard" & group1 %in% c("Prosper", "Casselton", "Colfax"))
    ) %>%
    filter(p.adj <= 0.05) %>%
    mutate(plot_signif = ifelse(p.adj.signif == "****", "***", p.adj.signif))
}

# Significance bracket position
if (nrow(stat.test) > 0) {
  stat.test$y.position <- seq(from = 5.8, by = 0.5, length.out = nrow(stat.test))
}

# Figure
p <- ggplot(dat_long, aes(x = Location, y = Score, fill = Location)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.6, color = "black") +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), width = 0.2) +
  geom_jitter(width = 0.05, alpha = 0.8, size = 16) +
  {
    if (nrow(stat.test) > 0)
      stat_pvalue_manual(
        stat.test,
        label = "plot_signif",
        y.position = "y.position",
        size = 70,
        tip.length = 0.01,
        linewidth = 6,
        inherit.aes = FALSE
      )
  } +
  scale_x_discrete(limits = loc_order, drop = FALSE) +
  scale_fill_manual(values = loc_colors, limits = loc_order, drop = FALSE) +
  scale_y_continuous(breaks = 1:8) +
  coord_cartesian(ylim = c(1, 8)) +
  labs(x = "Location", y = "IDC Score") +
  theme_classic(base_size = 120) +
  theme(
    axis.line = element_line(linewidth = 2),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 120, face = "bold", colour = "black"),
    axis.text.y = element_text(size = 120, face = "bold", colour = "black"),
    axis.title = element_text(size = 120, face = "bold", colour = "black"),
    legend.position = "none"
  )

# Save figure
ggsave(
  file.path(out_dir, "ND17009GT_IDC_scores_across_sites_years.png"),
  p,
  width = 25,
  height = 30,
  dpi = 300,
  bg = "white"
)
