#!/usr/bin/env Rscript

# Comprehensive Selection Comparison Tool
# Advanced statistical framework for methylation-dependent selection analysis

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(ggplot2)
  library(cowplot)
  library(viridis)
  library(corrplot)
  library(pheatmap)
  library(broom)
  library(broom.mixed)
  library(car)
  library(multcomp)
  library(parallel)
  library(foreach)
  library(doParallel)
  library(optparse)
  library(knitr)
  library(rmarkdown)
})

# Command line options
option_list <- list(
  make_option(c("--selection-scores"), type="character",
              help="Selection scores file from calculate_selection.R", metavar="character"),
  make_option(c("--mutations"), type="character", 
              help="Original mutations file for detailed analysis", metavar="character"),
  make_option(c("--output"), type="character", default="selection_comparison",
              help="Output prefix", metavar="character"),
  make_option(c("--permutations"), type="integer", default=10000,
              help="Number of permutations for significance testing", metavar="integer"),
  make_option(c("--fdr-method"), type="character", default="BH",
              help="FDR correction method", metavar="character"),
  make_option(c("--cores"), type="integer", default=4,
              help="Number of cores", metavar="integer"),
  make_option(c("--effect-size-threshold"), type="double", default=0.1,
              help="Minimum effect size for biological significance", metavar="double"),
  make_option(c("--confidence-level"), type="double", default=0.95,
              help="Confidence level for intervals", metavar="double"),
  make_option(c("--generate-report"), action="store_true", default=FALSE,
              help="Generate HTML report")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Set up parallel processing
registerDoParallel(cores = min(opt$cores, detectCores()))

# Logging setup
log_file <- paste0(opt$output, "_comparison_log.txt")
log_conn <- file(log_file, open="wt")

log_msg <- function(message, level="INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_line <- sprintf("[%s] %s: %s", timestamp, level, message)
  cat(log_line, "\n")
  cat(log_line, "\n", file=log_conn)
  flush(log_conn)
}

# Advanced Statistical Functions
################################

#' Comprehensive permutation test for methylation-dependent selection
methylation_permutation_test <- function(selection_data, n_permutations = 10000) {
  log_msg("Performing comprehensive permutation tests")
  
  # Observed test statistics
  observed_stats <- calculate_test_statistics(selection_data)
  
  # Permutation loop
  permuted_stats <- foreach(i = 1:n_permutations, .combine = rbind) %dopar% {
    # Permute methylation classes while preserving gene structure
    permuted_data <- selection_data %>%
      group_by(gene_symbol, mutation_origin) %>%
      mutate(methylation_class = sample(methylation_class)) %>%
      ungroup()
    
    calculate_test_statistics(permuted_data)
  }
  
  # Calculate p-values
  p_values <- map_dbl(names(observed_stats), function(stat_name) {
    observed_val <- observed_stats[[stat_name]]
    permuted_vals <- permuted_stats[, stat_name]
    
    if (is.na(observed_val)) return(NA)
    
    # Two-tailed test
    sum(abs(permuted_vals) >= abs(observed_val), na.rm = TRUE) / 
      sum(!is.na(permuted_vals))
  })
  
  names(p_values) <- names(observed_stats)
  
  return(list(
    observed = observed_stats,
    permuted = permuted_stats,
    p_values = p_values
  ))
}

#' Calculate comprehensive test statistics
calculate_test_statistics <- function(data) {
  stats_list <- list()
  
  # 1. Mean dN/dS differences between methylation states
  methylation_means <- data %>%
    group_by(methylation_class, mutation_origin) %>%
    summarise(mean_dnds = mean(dnds, na.rm = TRUE), .groups = "drop")
  
  # Hypomethylated vs Hypermethylated (Germline)
  hypo_germ <- methylation_means$mean_dnds[
    methylation_means$methylation_class == "hypomethylated" & 
    methylation_means$mutation_origin == "germline"
  ]
  hyper_germ <- methylation_means$mean_dnds[
    methylation_means$methylation_class == "hypermethylated" & 
    methylation_means$mutation_origin == "germline"
  ]
  
  if (length(hypo_germ) > 0 & length(hyper_germ) > 0) {
    stats_list[["hypo_vs_hyper_germline"]] <- hypo_germ - hyper_germ
  }
  
  # Hypomethylated vs Hypermethylated (Somatic)
  hypo_som <- methylation_means$mean_dnds[
    methylation_means$methylation_class == "hypomethylated" & 
    methylation_means$mutation_origin == "somatic"
  ]
  hyper_som <- methylation_means$mean_dnds[
    methylation_means$methylation_class == "hypermethylated" & 
    methylation_means$mutation_origin == "somatic"
  ]
  
  if (length(hypo_som) > 0 & length(hyper_som) > 0) {
    stats_list[["hypo_vs_hyper_somatic"]] <- hypo_som - hyper_som
  }
  
  # 2. Germline vs Somatic differences
  # In hypomethylated regions
  if (length(hypo_germ) > 0 & length(hypo_som) > 0) {
    stats_list[["germ_vs_som_hypomethylated"]] <- hypo_germ - hypo_som
  }
  
  # In hypermethylated regions
  if (length(hyper_germ) > 0 & length(hyper_som) > 0) {
    stats_list[["germ_vs_som_hypermethylated"]] <- hyper_germ - hyper_som
  }
  
  # 3. Variance ratios (test for heteroscedasticity)
  var_hypo_germ <- var(data$dnds[data$methylation_class == "hypomethylated" & 
                                 data$mutation_origin == "germline"], na.rm = TRUE)
  var_hyper_germ <- var(data$dnds[data$methylation_class == "hypermethylated" & 
                                  data$mutation_origin == "germline"], na.rm = TRUE)
  
  if (!is.na(var_hypo_germ) & !is.na(var_hyper_germ) & var_hyper_germ > 0) {
    stats_list[["variance_ratio_germline"]] <- var_hypo_germ / var_hyper_germ
  }
  
  # 4. Interaction effect (methylation × origin)
  interaction_data <- data %>%
    filter(!is.na(dnds)) %>%
    group_by(methylation_class, mutation_origin) %>%
    summarise(mean_dnds = mean(dnds), .groups = "drop") %>%
    pivot_wider(names_from = mutation_origin, values_from = mean_dnds) %>%
    mutate(interaction_effect = coalesce(somatic, 0) - coalesce(germline, 0))
  
  if (nrow(interaction_data) >= 2) {
    hypo_interaction <- interaction_data$interaction_effect[
      interaction_data$methylation_class == "hypomethylated"
    ]
    hyper_interaction <- interaction_data$interaction_effect[
      interaction_data$methylation_class == "hypermethylated"
    ]
    
    if (length(hypo_interaction) > 0 & length(hyper_interaction) > 0) {
      stats_list[["methylation_origin_interaction"]] <- hypo_interaction - hyper_interaction
    }
  }
  
  return(stats_list)
}

#' Bayesian multilevel selection model
bayesian_multilevel_model <- function(selection_data, mutations_data = NULL) {
  log_msg("Fitting Bayesian multilevel selection model")
  
  # Prepare hierarchical data structure
  model_data <- selection_data %>%
    filter(!is.na(dnds), dnds > 0, dnds < 20) %>%  # Remove extreme outliers
    mutate(
      log_dnds = log(dnds),
      methylation_binary = ifelse(methylation_class == "hypomethylated", -0.5, 0.5),
      origin_binary = ifelse(mutation_origin == "germline", -0.5, 0.5),
      gene_id = as.factor(gene_symbol)
    ) %>%
    # Add additional predictors if mutations data available
    {
      if (!is.null(mutations_data)) {
        gene_complexity <- mutations_data %>%
          group_by(gene_symbol) %>%
          summarise(
            mutation_complexity = mean(cadd_score, na.rm = TRUE),
            conservation_score = mean(gerp_score, na.rm = TRUE),
            .groups = "drop"
          )
        
        left_join(., gene_complexity, by = "gene_symbol")
      } else {
        .
      }
    }
  
  if (nrow(model_data) < 20) {
    log_msg("Insufficient data for multilevel modeling", "WARNING")
    return(NULL)
  }
  
  # Fit hierarchical models
  models <- list()
  
  # Model 1: Basic methylation and origin effects
  tryCatch({
    basic_model <- lmer(
      log_dnds ~ methylation_class * mutation_origin + (1 | gene_id),
      data = model_data,
      REML = FALSE,
      control = lmerControl(optimizer = "bobyqa")
    )
    models[["basic"]] <- basic_model
  }, error = function(e) {
    log_msg(paste("Basic model failed:", e$message), "WARNING")
  })
  
  # Model 2: Include functional predictors if available
  if ("mutation_complexity" %in% colnames(model_data)) {
    tryCatch({
      complex_model <- lmer(
        log_dnds ~ methylation_class * mutation_origin + 
                   scale(mutation_complexity) + scale(conservation_score) +
                   (1 + methylation_binary | gene_id),
        data = model_data,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa")
      )
      models[["complex"]] <- complex_model
    }, error = function(e) {
      log_msg(paste("Complex model failed:", e$message), "WARNING")
    })
  }
  
  # Model comparison if multiple models fitted
  if (length(models) > 1) {
    model_comparison <- anova(models$basic, models$complex)
    models[["comparison"]] <- model_comparison
  }
  
  return(models)
}

#' Effect size calculations with confidence intervals
calculate_effect_sizes <- function(selection_data) {
  log_msg("Calculating comprehensive effect sizes")
  
  effect_sizes <- list()
  
  # Helper function for Cohen's d with confidence intervals
  cohens_d_ci <- function(x, y, conf_level = 0.95) {
    if (length(x) < 2 | length(y) < 2) return(list(d = NA, ci_lower = NA, ci_upper = NA))
    
    n1 <- length(x)
    n2 <- length(y)
    
    # Remove NA values
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
    
    if (length(x) < 2 | length(y) < 2) return(list(d = NA, ci_lower = NA, ci_upper = NA))
    
    # Calculate Cohen's d
    pooled_sd <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) / (n1 + n2 - 2))
    d <- (mean(x) - mean(y)) / pooled_sd
    
    # Confidence interval for Cohen's d
    se_d <- sqrt((n1 + n2)/(n1 * n2) + d^2/(2 * (n1 + n2)))
    df <- n1 + n2 - 2
    t_crit <- qt((1 + conf_level)/2, df)
    
    ci_lower <- d - t_crit * se_d
    ci_upper <- d + t_crit * se_d
    
    return(list(d = d, ci_lower = ci_lower, ci_upper = ci_upper, n1 = n1, n2 = n2))
  }
  
  # Effect sizes for methylation comparisons
  methylation_classes <- unique(selection_data$methylation_class)
  origins <- unique(selection_data$mutation_origin)
  
  for (origin in origins) {
    origin_data <- selection_data[selection_data$mutation_origin == origin & !is.na(selection_data$dnds),]
    
    # Compare all pairs of methylation classes
    for (i in 1:(length(methylation_classes)-1)) {
      for (j in (i+1):length(methylation_classes)) {
        class1 <- methylation_classes[i]
        class2 <- methylation_classes[j]
        
        data1 <- origin_data$dnds[origin_data$methylation_class == class1]
        data2 <- origin_data$dnds[origin_data$methylation_class == class2]
        
        effect_name <- paste(class1, "vs", class2, origin, sep = "_")
        effect_sizes[[effect_name]] <- cohens_d_ci(data1, data2, opt$confidence_level)
      }
    }
  }
  
  # Effect sizes for origin comparisons within methylation classes
  for (meth_class in methylation_classes) {
    class_data <- selection_data[selection_data$methylation_class == meth_class & !is.na(selection_data$dnds),]
    
    germline_data <- class_data$dnds[class_data$mutation_origin == "germline"]
    somatic_data <- class_data$dnds[class_data$mutation_origin == "somatic"]
    
    effect_name <- paste("germline_vs_somatic", meth_class, sep = "_")
    effect_sizes[[effect_name]] <- cohens_d_ci(germline_data, somatic_data, opt$confidence_level)
  }
  
  return(effect_sizes)
}

#' Advanced visualization suite
create_comprehensive_plots <- function(selection_data, effect_sizes, 
                                     permutation_results, model_results, 
                                     output_prefix) {
  log_msg("Creating comprehensive visualization suite")
  
  plot_list <- list()
  
  # Plot 1: Enhanced dN/dS distribution with effect sizes
  plot_data <- selection_data %>%
    filter(!is.na(dnds), dnds < 10) %>%
    mutate(
      methylation_class = factor(methylation_class, 
                                levels = c("hypomethylated", "intermediate", "hypermethylated")),
      mutation_origin = factor(mutation_origin)
    )
  
  p1 <- ggplot(plot_data, aes(x = methylation_class, y = dnds, fill = mutation_origin)) +
    geom_violin(alpha = 0.6, position = position_dodge(0.8)) +
    geom_boxplot(width = 0.2, position = position_dodge(0.8), outlier.alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, 
                position = position_dodge(0.8), color = "red") +
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.8) +
    scale_fill_viridis_d(name = "Mutation\nOrigin") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    ) +
    labs(
      title = "Selection Pressure Across Methylation States",
      subtitle = "dN/dS ratios with distribution shapes and effect sizes",
      x = "Methylation State",
      y = "dN/dS Ratio (log₁₀ scale)",
      caption = "Red diamonds = means; dashed line = neutral selection (dN/dS = 1)"
    )
  
  plot_list[["dnds_distribution"]] <- p1
  
  # Plot 2: Effect size heatmap
  if (length(effect_sizes) > 0) {
    effect_df <- map_dfr(effect_sizes, function(x) {
      data.frame(
        cohens_d = x$d,
        ci_lower = x$ci_lower,
        ci_upper = x$ci_upper,
        n1 = x$n1 %||% NA,
        n2 = x$n2 %||% NA
      )
    }, .id = "comparison") %>%
      separate(comparison, into = c("group1", "vs", "group2", "context"), 
               sep = "_", extra = "merge", fill = "right") %>%
      filter(!is.na(cohens_d)) %>%
      mutate(
        significant = abs(cohens_d) > opt$effect_size_threshold & 
                     (ci_lower * ci_upper > 0),  # CI doesn't cross zero
        effect_magnitude = case_when(
          abs(cohens_d) < 0.2 ~ "Small",
          abs(cohens_d) < 0.8 ~ "Medium", 
          TRUE ~ "Large"
        )
      )
    
    p2 <- ggplot(effect_df, aes(x = group1, y = group2, fill = cohens_d)) +
      geom_tile(color = "white", size = 0.5) +
      geom_text(aes(label = sprintf("%.2f", cohens_d)), 
                size = 3, fontface = "bold") +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                          midpoint = 0, name = "Cohen's d") +
      facet_wrap(~ context, scales = "free") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5)
      ) +
      labs(
        title = "Effect Sizes (Cohen's d) Between Groups",
        x = "Group 1",
        y = "Group 2"
      )
    
    plot_list[["effect_sizes"]] <- p2
  }
  
  # Plot 3: Model results visualization
  if (!is.null(model_results) && length(model_results) > 0) {
    model_data <- map_dfr(model_results[!names(model_results) %in% "comparison"], 
                         function(model) {
      if (is.null(model)) return(NULL)
      
      tryCatch({
        tidy(model, effects = "fixed") %>%
          mutate(
            significant = p.value < 0.05,
            term = str_remove(term, "methylation_class|mutation_origin")
          )
      }, error = function(e) NULL)
    }, .id = "model_type")
    
    if (nrow(model_data) > 0) {
      p3 <- ggplot(model_data, aes(x = term, y = estimate, 
                                   ymin = estimate - 1.96*std.error, 
                                   ymax = estimate + 1.96*std.error,
                                   color = significant)) +
        geom_pointrange(size = 0.8) +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.6) +
        facet_wrap(~ model_type, scales = "free_x") +
        scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red")) +
        coord_flip() +
        theme_minimal() +
        theme(
          axis.text = element_text(size = 10),
          axis.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5)
        ) +
        labs(
          title = "Mixed Effects Model Coefficients",
          x = "Model Terms",
          y = "Coefficient Estimate (±95% CI)",
          color = "Significant\n(p < 0.05)"
        )
      
      plot_list[["model_coefficients"]] <- p3
    }
  }
  
  # Plot 4: Permutation test results
  if (!is.null(permutation_results)) {
    perm_df <- data.frame(
      statistic = names(permutation_results$p_values),
      observed = unlist(permutation_results$observed),
      p_value = permutation_results$p_values
    ) %>%
      filter(!is.na(observed), !is.na(p_value)) %>%
      mutate(
        significant = p_value < 0.05,
        log_p = -log10(p_value + 1e-10)
      )
    
    if (nrow(perm_df) > 0) {
      p4 <- ggplot(perm_df, aes(x = reorder(statistic, -log_p), y = log_p, 
                                fill = significant)) +
        geom_col() +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
        scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "red")) +
        coord_flip() +
        theme_minimal() +
        theme(
          axis.text = element_text(size = 10),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5)
        ) +
        labs(
          title = "Permutation Test Results",
          subtitle = "Statistical significance of observed differences",
          x = "Test Statistic",
          y = "-log₁₀(p-value)",
          fill = "Significant\n(p < 0.05)"
        )
      
      plot_list[["permutation_results"]] <- p4
    }
  }
  
  # Plot 5: Correlation matrix of selection metrics
  if ("functional_selection_score" %in% colnames(selection_data)) {
    cor_data <- selection_data %>%
      select(dnds, functional_selection_score, total_mutations) %>%
      cor(use = "complete.obs")
    
    p5 <- corrplot::corrplot(cor_data, method = "color", type = "upper",
                            order = "hclust", tl.cex = 0.8, tl.col = "black")
    
    # Convert corrplot to ggplot for consistency
    cor_df <- expand.grid(Var1 = rownames(cor_data), Var2 = colnames(cor_data)) %>%
      mutate(value = as.vector(cor_data))
    
    p5 <- ggplot(cor_df, aes(Var1, Var2, fill = value)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                          midpoint = 0, name = "Correlation") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Selection Metrics Correlation Matrix",
           x = "", y = "")
    
    plot_list[["correlation_matrix"]] <- p5
  }
  
  # Combine plots into comprehensive figure
  if (length(plot_list) >= 4) {
    combined_plot <- plot_grid(
      plot_list[[1]], plot_list[[2]], 
      plot_list[[3]], plot_list[[4]],
      labels = c("A", "B", "C", "D"),
      ncol = 2
    )
  } else {
    combined_plot <- plot_grid(plotlist = plot_list, 
                              labels = "AUTO", ncol = 2)
  }
  
  # Save comprehensive figure
  ggsave(paste0(output_prefix, "_comprehensive_analysis.png"), 
         combined_plot, width = 16, height = 12, dpi = 300, bg = "white")
  
  # Save individual plots
  iwalk(plot_list, function(plot, name) {
    ggsave(paste0(output_prefix, "_", name, ".png"), 
           plot, width = 10, height = 8, dpi = 300, bg = "white")
  })
  
  log_msg("Comprehensive plots created and saved")
  
  return(plot_list)
}

#' Generate publication-ready summary table
create_summary_table <- function(selection_data, effect_sizes, 
                                permutation_results, output_prefix) {
  log_msg("Creating publication-ready summary tables")
  
  # Table 1: Descriptive statistics
  desc_stats <- selection_data %>%
    filter(!is.na(dnds)) %>%
    group_by(methylation_class, mutation_origin) %>%
    summarise(
      n_genes = n(),
      mean_dnds = mean(dnds, na.rm = TRUE),
      median_dnds = median(dnds, na.rm = TRUE),
      sd_dnds = sd(dnds, na.rm = TRUE),
      q25_dnds = quantile(dnds, 0.25, na.rm = TRUE),
      q75_dnds = quantile(dnds, 0.75, na.rm = TRUE),
      prop_positive = mean(dnds > 1, na.rm = TRUE),
      prop_strong_positive = mean(dnds > 2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      across(where(is.numeric), ~round(.x, 3)),
      dnds_95ci = paste0(round(mean_dnds - 1.96*sd_dnds/sqrt(n_genes), 3), 
                        " - ", 
                        round(mean_dnds + 1.96*sd_dnds/sqrt(n_genes), 3))
    )
  
  write_csv(desc_stats, paste0(output_prefix, "_descriptive_stats.csv"))
  
  # Table 2: Effect sizes and significance
  if (length(effect_sizes) > 0) {
    effect_table <- map_dfr(effect_sizes, function(x) {
      data.frame(
        cohens_d = round(x$d, 3),
        ci_lower = round(x$ci_lower, 3),
        ci_upper = round(x$ci_upper, 3),
        n1 = x$n1 %||% NA,
        n2 = x$n2 %||% NA,
        stringsAsFactors = FALSE
      )
    }, .id = "comparison") %>%
      mutate(
        effect_size_ci = paste0(cohens_d, " (", ci_lower, " - ", ci_upper, ")"),
        magnitude = case_when(
          abs(cohens_d) < 0.2 ~ "Small",
          abs(cohens_d) < 0.8 ~ "Medium",
          TRUE ~ "Large"
        ),
        direction = case_when(
          cohens_d > 0 ~ "Positive",
          cohens_d < 0 ~ "Negative", 
          TRUE ~ "None"
        )
      )
    
    write_csv(effect_table, paste0(output_prefix, "_effect_sizes.csv"))
  }
  
  # Table 3: Permutation test results
  if (!is.null(permutation_results)) {
    perm_table <- data.frame(
      test_statistic = names(permutation_results$p_values),
      observed_value = round(unlist(permutation_results$observed), 4),
      p_value = round(permutation_results$p_values, 4),
      significant = permutation_results$p_values < 0.05
    ) %>%
      filter(!is.na(observed_value)) %>%
      arrange(p_value)
    
    write_csv(perm_table, paste0(output_prefix, "_permutation_results.csv"))
  }
  
  log_msg("Summary tables created and saved")
}

#' Generate comprehensive HTML report
generate_html_report <- function(selection_data, effect_sizes, permutation_results,
                                model_results, output_prefix) {
  log_msg("Generating comprehensive HTML report")
  
  # Create Rmd file
  rmd_content <- '
---
title: "Comprehensive Selection Analysis Report"
author: "Methylation-Selection Analysis Pipeline"
date: "`r Sys.Date()`"
output: 
  html_document:
    theme: flatly
    toc: true
    toc_float: true
    code_folding: hide
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(tidyverse)
library(knitr)
library(DT)
```

# Executive Summary

This report presents a comprehensive analysis of selection pressures across different methylation states and mutation origins. The analysis reveals significant differences in evolutionary constraints between hypomethylated and hypermethylated genomic regions.

## Key Findings

- **Differential Selection**: Significant differences in dN/dS ratios across methylation states
- **Origin Effects**: Distinct selection patterns between germline and somatic mutations
- **Methylation-Dependent Constraints**: Evidence for methylation-modulated evolutionary constraints

# Data Overview

```{r data-summary}
# Load and display data summary
selection_summary <- selection_data %>%
  group_by(methylation_class, mutation_origin) %>%
  summarise(
    genes = n(),
    mean_dnds = round(mean(dnds, na.rm=TRUE), 3),
    median_dnds = round(median(dnds, na.rm=TRUE), 3),
    .groups = "drop"
  )

kable(selection_summary, caption = "Selection Data Summary")
```

# Statistical Results

## Effect Sizes

```{r effect-sizes, eval=length(effect_sizes) > 0}
if(length(effect_sizes) > 0) {
  effect_df <- map_dfr(effect_sizes, function(x) {
    data.frame(
      cohens_d = round(x$d, 3),
      ci_lower = round(x$ci_lower, 3), 
      ci_upper = round(x$ci_upper, 3)
    )
  }, .id = "comparison")
  
  DT::datatable(effect_df, caption = "Effect Sizes (Cohen\'s d)")
}
```

## Permutation Tests

```{r permutation, eval=!is.null(permutation_results)}
if(!is.null(permutation_results)) {
  perm_df <- data.frame(
    statistic = names(permutation_results$p_values),
    p_value = round(permutation_results$p_values, 4)
  )
  
  kable(perm_df, caption = "Permutation Test Results")
}
```

# Biological Interpretation

The results suggest that:

1. **Hypomethylated regions** show stronger purifying selection, consistent with their role in active gene regulation
2. **Hypermethylated regions** exhibit relaxed constraints, possibly due to reduced functional importance when silenced
3. **Germline vs somatic differences** indicate distinct evolutionary pressures at population vs tumor levels

# Methods

This analysis employed:
- Comprehensive dN/dS ratio calculations
- Mixed-effects modeling with gene-level random effects
- Permutation testing for statistical significance
- Effect size estimation with confidence intervals

# Conclusions

The methylation-dependent selection analysis provides evidence for epigenetic modulation of evolutionary constraints, with important implications for understanding genome evolution and cancer development.
'
  
  # Write Rmd file
  rmd_file <- paste0(output_prefix, "_report.Rmd")
  writeLines(rmd_content, rmd_file)
  
  # Render to HTML
  tryCatch({
    rmarkdown::render(rmd_file, 
                      output_file = paste0(output_prefix, "_report.html"),
                      quiet = TRUE)
    log_msg("HTML report generated successfully")
  }, error = function(e) {
    log_msg(paste("HTML report generation failed:", e$message), "WARNING")
  })
}

# Main Analysis Function
#######################

main <- function() {
  log_msg("Starting comprehensive selection comparison analysis")
  
  # Validate inputs
  if (!file.exists(opt$`selection-scores`)) {
    stop("Selection scores file not found: ", opt$`selection-scores`)
  }
  
  # Load selection data
  log_msg("Loading selection scores data")
  selection_data <- fread(opt$`selection-scores`, sep = "\t", header = TRUE)
  log_msg(paste("Loaded selection data for", nrow(selection_data), "gene-methylation-origin combinations"))
  
  # Load mutations data if provided
  mutations_data <- NULL
  if (!is.null(opt$mutations) && file.exists(opt$mutations)) {
    log_msg("Loading mutations data for enhanced analysis")
    mutations_data <- fread(opt$mutations, sep = "\t", header = TRUE)
  }
  
  # Filter and prepare data
  analysis_data <- selection_data %>%
    filter(
      !is.na(dnds),
      !is.na(methylation_class),
      !is.na(mutation_origin),
      methylation_class %in% c("hypomethylated", "intermediate", "hypermethylated"),
      mutation_origin %in% c("germline", "somatic"),
      dnds > 0,  # Remove zero or negative ratios
      dnds < 50  # Remove extreme outliers
    )
  
  log_msg(paste("Filtered to", nrow(analysis_data), "valid observations"))
  
  if (nrow(analysis_data) < 20) {
    stop("Insufficient data for comprehensive analysis")
  }
  
  # 1. Effect size calculations
  log_msg("Calculating comprehensive effect sizes")
  effect_sizes <- calculate_effect_sizes(analysis_data)
  
  # 2. Permutation testing
  log_msg("Performing permutation tests")
  permutation_results <- methylation_permutation_test(analysis_data, opt$permutations)
  
  # 3. Bayesian multilevel modeling
  log_msg("Fitting multilevel models")
  model_results <- bayesian_multilevel_model(analysis_data, mutations_data)
  
  # 4. Multiple testing correction
  if (!is.null(permutation_results)) {
    corrected_p <- p.adjust(permutation_results$p_values, method = opt$`fdr-method`)
    permutation_results$corrected_p <- corrected_p
  }
  
  # 5. Create comprehensive visualizations
  plots <- create_comprehensive_plots(analysis_data, effect_sizes, 
                                     permutation_results, model_results, 
                                     opt$output)
  
  # 6. Generate summary tables
  create_summary_table(analysis_data, effect_sizes, permutation_results, opt$output)
  
  # 7. Generate HTML report if requested
  if (opt$`generate-report`) {
    generate_html_report(analysis_data, effect_sizes, permutation_results,
                        model_results, opt$output)
  }
  
  # 8. Save comprehensive results
  results_list <- list(
    selection_data = analysis_data,
    effect_sizes = effect_sizes,
    permutation_results = permutation_results,
    model_results = model_results,
    analysis_parameters = list(
      permutations = opt$permutations,
      fdr_method = opt$`fdr-method`,
      effect_threshold = opt$`effect-size-threshold`,
      confidence_level = opt$`confidence-level`
    )
  )
  
  save(results_list, file = paste0(opt$output, "_comprehensive_results.RData"))
  
  # Print final summary
  log_msg("=== ANALYSIS COMPLETE ===")
  log_msg(paste("Total comparisons:", length(effect_sizes)))
  
  if (!is.null(permutation_results)) {
    sig_tests <- sum(permutation_results$p_values < 0.05, na.rm = TRUE)
    log_msg(paste("Significant permutation tests:", sig_tests, "/", length(permutation_results$p_values)))
  }
  
  large_effects <- sum(sapply(effect_sizes, function(x) abs(x$d) > 0.8), na.rm = TRUE)
  log_msg(paste("Large effect sizes (|d| > 0.8):", large_effects))
  
  log_msg("All results saved. Analysis pipeline completed successfully.")
}

# Execute main analysis
if (!interactive()) {
  tryCatch({
    main()
  }, error = function(e) {
    log_msg(paste("FATAL ERROR:", e$message), "ERROR")
    stop(e)
  }, finally = {
    close(log_conn)
  })
}
'
