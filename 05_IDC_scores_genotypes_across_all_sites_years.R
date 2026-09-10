# IDC scores of soybean genotypes across four sites from 2021 to 2023 
# No fertilizer samples were used from year 2022 and 2023 to assess IDC across the sites

# Load required packages
library(tidyverse)
library(rstatix)
library(ggpubr)

# Output directory
out_dir <- "IDC_Scores_across_genotypes_sites_years"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Location colors
location_colors <- c(
  "Prosper" = "orange",
  "Casselton" = "red",
  "Colfax" = "blue",
  "Leonard" = "darkgreen"
)

# 2021 IDC data

dat_wide_2021 <- tribble(
  ~Location, ~A11, ~ND_Rolette, ~ND17009GT, ~ND16.7108,
  "Prosper", 1, 1, 1, 1,
  "Prosper", 1, 1, 1, 1,
  "Prosper", 1, 1, 1, 1,
  "Prosper", 1, 1, 1, 1,
  "Casselton", 1, 1.5, 3, 2,
  "Casselton", 1, 1.5, 3, 2,
  "Casselton", 2, 1, 3, 1.5,
  "Casselton", 2, 1, 3, 1.5,
  "Colfax", 1, 1, 1, 1,
  "Colfax", 1, 2.5, 1.5, 1.5,
  "Colfax", 1, 1.5, 2, 2.5,
  "Colfax", 1, 1, 2.5, 3,
  "Leonard", 1, 1.5, 3, 3,
  "Leonard", 1, 1.5, 3, 2.5,
  "Leonard", 1, 1.5, 3, 2.5,
  "Leonard", 1, 1.5, 3, 2.5
)

# 2022 IDC data

dat_wide_2022 <- tribble(
  ~Location, ~A11, ~ND_Rolette, ~ND17009GT, ~ND18_17323,
  "Casselton", 1, 1, 1, 1,
  "Casselton", 1, 1, 1, 1,
  "Casselton", 1, 1, 1, 1,
  "Casselton", 1, 1, 1, 1,
  "Leonard", 1.5, 4, 5, 4.5,
  "Leonard", 2, 3, 5, 5,
  "Leonard", 1.5, 4, 5, 5,
  "Leonard", 1, 4, 5, 4.5
)

# 2023 IDC data

dat_wide_2023 <- tribble(
  ~Location, ~A11, ~ND_Rolette, ~ND17009GT, ~ND18_17323,
  "Prosper", 1, 1, 1, 1,
  "Prosper", 1, 1, 1, 1,
  "Prosper", 1, 1, 1, 1,
  "Prosper", 1, 1, 1, 1,
  "Casselton", 1, 1, 1, 1,
  "Casselton", 1, 1, 1, 1,
  "Casselton", 1, 1, 1, 1,
  "Casselton", 1, 1, 1, 1,
  "Colfax", 1, 1, 5, 3,
  "Colfax", 1.5, 1, 4.5, 3.5,
  "Colfax", 1, 3, 4, 4,
  "Colfax", 2.5, 1.5, 5, 4,
  "Leonard", 1, 2, 5, 3,
  "Leonard", 1, 2, 4.5, 3.5,
  "Leonard", 1.5, 1, 4, 3,
  "Leonard", 1.5, 1, 5, 2.5
)


# Statistical analysis and IDC plots
run_idc_plot <- function(dat_wide, year_label, genotype_levels, location_levels, file_prefix) {
  
  dat_long <- dat_wide %>%
    pivot_longer(cols = -Location, names_to = "Genotype", values_to = "Value") %>%
    drop_na(Value) %>%
    mutate(
      Genotype = factor(Genotype, levels = genotype_levels),
      Location = factor(Location, levels = location_levels)
    )

  
# Shapiro-Wilk normality test
sw_all <- do.call(rbind, lapply(location_levels, function(loc) {
    
    do.call(rbind, lapply(genotype_levels, function(g) {
      
      vals <- dat_long$Value[
        dat_long$Location == loc &
          dat_long$Genotype == g
      ]
      
      if (length(vals) < 3) {
        return(data.frame(
          Location = loc,
          Genotype = g,
          n = length(vals),
          W = NA,
          p = NA,
          Normal = "SKIP - too few"
        ))
      }
      
      if (length(unique(vals)) < 2) {
        return(data.frame(
          Location = loc,
          Genotype = g,
          n = length(vals),
          W = NA,
          p = NA,
          Normal = "SKIP - all identical"
        ))
      }
      
      sw <- shapiro.test(vals)
      
      data.frame(
        Location = loc,
        Genotype = g,
        n = length(vals),
        W = round(as.numeric(sw$statistic), 4),
        p = round(sw$p.value, 4),
        Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
      )
    }))
  }))
  
  cat("\n", year_label, "- Shapiro-Wilk normality test\n")
  print(sw_all, row.names = FALSE)
  

  # Statistical tests
stat_all <- list()
  
  for (loc in location_levels) {
    
    dat_loc <- dat_long %>%
      filter(Location == loc)
    
    normality <- sw_all %>%
      filter(Location == loc)
    
    all_normal <- all(normality$Normal == "YES")
    
    # Skip inferential statistics when all values are identical
    
    if (length(unique(dat_loc$Value)) < 2) {
      
      cat("\n", year_label, "-", loc, "- No test: all values identical\n")
      
      stat.test <- data.frame()
      

# ANOVA and Tukey HSD when all groups meet normality
} else if (all_normal) {
      
      anova_result <- dat_loc %>%
        anova_test(Value ~ Genotype)
      
      tukey_result <- dat_loc %>%
        tukey_hsd(Value ~ Genotype)
      
      cat("\n", year_label, "-", loc, "- One-way ANOVA\n")
      print(anova_result)
      
      cat("\nTukey HSD\n")
      print(tukey_result)
      
      stat.test <- tukey_result %>%
        filter(
          (group1 == "A11" & group2 == "ND17009GT") |
            (group1 == "ND17009GT" & group2 == "A11")
        ) %>%
        filter(p.adj <= 0.05) %>%
        mutate(
          Location = loc,
          plot_signif = ifelse(p.adj.signif == "****", "***", p.adj.signif)
        )

            
# Kruskal-Wallis and pairwise Wilcoxon when normality is not met or cannot be assessed
} else {
      
      kw_result <- dat_loc %>%
        kruskal_test(Value ~ Genotype)
      
      comparisons <- combn(genotype_levels, 2, simplify = FALSE)
      
      wilcox_result <- do.call(rbind, lapply(comparisons, function(comp) {
        
        group1 <- comp[1]
        group2 <- comp[2]
        
        x <- dat_loc$Value[dat_loc$Genotype == group1]
        y <- dat_loc$Value[dat_loc$Genotype == group2]
        
        x <- x[!is.na(x)]
        y <- y[!is.na(y)]
        
        if (length(unique(c(x, y))) < 2) {
          
          statistic <- NA
          p <- 1
          
        } else {
          
          test <- wilcox.test(
            x,
            y,
            paired = FALSE,
            exact = FALSE
          )
          
          statistic <- as.numeric(test$statistic)
          p <- test$p.value
        }
        
        data.frame(
          .y. = "Value",
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
      
      cat("\n", year_label, "-", loc, "- Kruskal-Wallis\n")
      print(kw_result)
      
      cat("\nPairwise Wilcoxon with BH correction\n")
      print(wilcox_result)
      
      stat.test <- wilcox_result %>%
        filter(
          (group1 == "A11" & group2 == "ND17009GT") |
            (group1 == "ND17009GT" & group2 == "A11")
        ) %>%
        filter(p.adj <= 0.05) %>%
        mutate(
          Location = loc,
          plot_signif = ifelse(p.adj.signif == "****", "***", p.adj.signif)
        )
    }
    
    if (nrow(stat.test) > 0) stat_all[[loc]] <- stat.test
  }
  
stat_plot <- bind_rows(stat_all)

  
# Plot positions
x_levels <- unlist(lapply(location_levels, function(loc) paste(loc, genotype_levels, sep = "___")))
x_labels <- sub("^.*___", "", x_levels)
xpos <- setNames(seq_along(x_levels), x_levels)
  
  dat_long <- dat_long %>%
    mutate(XID = factor(paste(Location, Genotype, sep = "___"), levels = x_levels))
  
  
# Significance bracket positions
if (nrow(stat_plot) > 0) {
    
    bar_means <- dat_long %>%
      group_by(Location, Genotype) %>%
      summarise(m = mean(Value, na.rm = TRUE), .groups = "drop") %>%
      group_by(Location) %>%
      summarise(top = max(m, na.rm = TRUE), .groups = "drop")
    
    stat_plot <- stat_plot %>%
      mutate(
        xmin = xpos[paste(Location, "A11", sep = "___")],
        xmax = xpos[paste(Location, "ND17009GT", sep = "___")]
      ) %>%
      left_join(bar_means, by = "Location") %>%
      mutate(y.position = top + 0.6)
  }
  

# Plot
p <- ggplot(dat_long, aes(x = as.integer(XID), y = Value, fill = Location)) +
    geom_bar(stat = "summary", fun = "mean", width = 0.9, color = "black") +
    geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), width = 0.25) +
    geom_jitter(aes(x = as.integer(XID)), width = 0.15, height = 0, alpha = 1, size = 16, color = "black") +
    {
      if (nrow(stat_plot) > 0)
        stat_pvalue_manual(
          stat_plot,
          label = "plot_signif",
          xmin = "xmin",
          xmax = "xmax",
          y.position = "y.position",
          tip.length = 0.03,
          size = 70,
          linewidth = 6,
          inherit.aes = FALSE
        )
    } +
    scale_fill_manual(values = location_colors) +
    scale_x_continuous(
      breaks = seq_along(x_levels),
      labels = x_labels,
      expand = expansion(add = 0.6)
    ) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 6),
      limits = c(0, 7)
    ) +
    labs(
      x = "Soybean Genotype",
      y = "IDC Scores",
      title = year_label
    ) +
    theme_classic(base_size = 60) +
    theme(
      plot.title = element_text(size = 36, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 60, face = "bold", colour = "black"),
      axis.text.y = element_text(size = 60, face = "bold", colour = "black"),
      axis.title = element_text(size = 60, face = "bold", colour = "black"),
      legend.position = "none"
    )
  
# Save figure
ggsave(
    file.path(out_dir, paste0(file_prefix, ".png")),
    p,
    width = 24,
    height = 18,
    dpi = 500,
    bg = "white"
  )
  
}

# 2021
p2021 <- run_idc_plot(
  dat_wide = dat_wide_2021,
  year_label = "2021",
  genotype_levels = c("A11", "ND_Rolette", "ND17009GT", "ND16.7108"),
  location_levels = c("Prosper", "Casselton", "Colfax", "Leonard"),
  file_prefix = "IDCscores_across_genotypes_sites_2021"
)

# 2022
p2022 <- run_idc_plot(
  dat_wide = dat_wide_2022,
  year_label = "2022",
  genotype_levels = c("A11", "ND_Rolette", "ND17009GT", "ND18_17323"),
  location_levels = c("Casselton", "Leonard"),
  file_prefix = "IDCscores_across_genotypes_sites_no_fertilizer_2022"
)

# 2023
p2023 <- run_idc_plot(
  dat_wide = dat_wide_2023,
  year_label = "2023",
  genotype_levels = c("A11", "ND_Rolette", "ND17009GT", "ND18_17323"),
  location_levels = c("Prosper", "Casselton", "Colfax", "Leonard"),
  file_prefix = "IDCscores_across_genotypes_sites_no_fertilizer_2023"
)
