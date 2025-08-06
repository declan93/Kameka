#!/usr/bin/env Rscript

# Selection Calculation Tool
# Comprehensive selection inference across methylation states and mutation origins

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(lme4)
  library(ggplot2)
  library(gridExtra)
  library(broom)
  library(coin)
  library(MASS)
  library(parallel)
  library(optparse)
})

# Command line options
option_list <- list(
  make_option(c("--mutations"), type="character", default=NULL,
              help="Annotated mutations file", metavar="character"),
  make_option(c("--methylation"), type="character", default=NULL,
              help="Methylation regions file", metavar="character"),
  make_option(c("--output"), type="character", default="selection_results",
              help="Output prefix", metavar="character"),
  make_option(c("--gene-lengths"), type="character", default=NULL,
              help="Gene length file for proper dN/dS calculation", metavar="character"),
  make_option(c("--codon-usage"), type="character", default=NULL,
              help="Codon usage table for accurate dN/dS", metavar="character"),
  make_option(c("--min-mutations"), type="integer", default=5,
              help="Minimum mutations per gene for analysis", metavar="integer"),
  make_option(c("--bootstrap-n"), type="integer", default=1000,
              help="Number of bootstrap replicates", metavar="integer"),
  make_option(c("--cores"), type="integer", default=4,
              help="Number of cores for parallel processing", metavar="integer"),
  make_option(c("--alpha"), type="double", default=0.05,
              help="Significance threshold", metavar="double")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Set up logging
log_file <- paste0(opt$output, "_selection_log.txt")
log_connection <- file(log_file, open="wt")

log_message <- function(message, level="INFO") {
  timestamp <- Sys.time()
  log_line <- sprintf("[%s] %s: %s", timestamp, level, message)
  cat(log_line, "\n")
  cat(log_line, "\n", file=log_connection)
  flush(log_connection)
}

# Selection Analysis Functions
###############################

#' Calculate dN/dS ratio with proper statistical framework
calculate_dnds <- function(syn_count, nonsyn_count, syn_sites, nonsyn_sites) {
  # Avoid division by zero
  if (syn_sites == 0 | syn_count == 0) {
    return(list(dnds = NA, ci_lower = NA, ci_upper = NA, p_value = NA))
  }
  
  dn <- nonsyn_count / nonsyn_sites
  ds <- syn_count / syn_sites
  dnds <- dn / ds
  
  # Calculate confidence interval using Fisher's method
  if (syn_count >= 5 & nonsyn_count >= 5) {
    # Use Poisson approximation for confidence intervals
    syn_rate <- syn_count / syn_sites
    nonsyn_rate <- nonsyn_count / nonsyn_sites
    
    # Approximate variance using delta method
    var_log_dnds <- (1/syn_count) + (1/nonsyn_count)
    se_log_dnds <- sqrt(var_log_dnds)
    
    ci_lower <- exp(log(dnds) - 1.96 * se_log_dnds)
    ci_upper <- exp(log(dnds) + 1.96 * se_log_dnds)
    
    # Test for neutrality (dN/dS = 1)
    z_score <- (log(dnds) - 0) / se_log_dnds
    p_value <- 2 * pnorm(-abs(z_score))
  } else {
    ci_lower <- NA
    ci_upper <- NA
    p_value <- NA
  }
  
  return(list(
    dnds = dnds,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p_value = p_value,
    syn_count = syn_count,
    nonsyn_count = nonsyn_count,
    syn_sites = syn_sites,
    nonsyn_sites = nonsyn_sites
  ))
}

#' Calculate site frequency spectrum-based selection metrics
calculate_sfs_selection <- function(mutations_df, methylation_class, mutation_origin) {
  subset_data <- mutations_df[
    mutations_df$methylation_class == methylation_class & 
    mutations_df$mutation_origin == mutation_origin,
  ]
  
  if (nrow(subset_data) < 10) {
    return(list(
      tajima_d = NA,
      fu_li_d = NA,
      fay_wu_h = NA,
      selection_score = NA
    ))
  }
  
  # Simple approximation of site frequency spectrum metrics
  # In practice, you'd need allele frequencies from population data
  
  # Mock calculation - replace with actual SFS analysis
  allele_freqs <- runif(nrow(subset_data), 0.01, 0.99)  # Placeholder
  
  # Tajima's D approximation
  theta_pi <- mean(allele_freqs * (1 - allele_freqs))
  theta_w <- length(allele_freqs) / sum(1/1:length(allele_freqs))
  tajima_d <- (theta_pi - theta_w) / sqrt(var(allele_freqs) / length(allele_freqs))
  
  return(list(
    tajima_d = tajima_d,
    fu_li_d = NA,  # Would need more complex calculation
    fay_wu_h = NA,  # Would need ancestral state information
    selection_score = tajima_d
  ))
}

#' Calculate McDonald-Kreitman test statistic
mcdonald_kreitman_test <- function(mutations_df, gene_symbol, methylation_class) {
  gene_data <- mutations_df[
    mutations_df$gene_symbol == gene_symbol & 
    mutations_df$methylation_class == methylation_class,
  ]
  
  if (nrow(gene_data) < 10) return(NA)
  
  # Separate polymorphisms from fixed differences
  # This is simplified - in practice you'd need interspecies comparison
  polymorphisms <- gene_data[gene_data$mutation_origin == "somatic",]
  fixed_diffs <- gene_data[gene_data$mutation_origin == "germline",]  # Approximation
  
  if (nrow(polymorphisms) == 0 | nrow(fixed_diffs) == 0) return(NA)
  
  # Count synonymous and non-synonymous in each category
  poly_syn <- sum(polymorphisms$mutation_type_category %in% c("low_impact", "synonymous"))
  poly_nonsyn <- sum(polymorphisms$mutation_type_category %in% c("moderate_impact", "high_impact"))
  fixed_syn <- sum(fixed_diffs$mutation_type_category %in% c("low_impact", "synonymous"))
  fixed_nonsyn <- sum(fixed_diffs$mutation_type_category %in% c("moderate_impact", "high_impact"))
  
  if (poly_syn == 0 | fixed_syn == 0) return(NA)
  
  # Calculate MK statistic
  mk_matrix <- matrix(c(fixed_nonsyn, fixed_syn, poly_nonsyn, poly_syn), nrow=2)
  
  if (sum(mk_matrix) < 10) return(NA)
  
  tryCatch({
    fisher_test <- fisher.test(mk_matrix)
    return(list(
      mk_statistic = fisher_test$estimate,
      mk_p_value = fisher_test$p.value,
      neutrality_index = (poly_nonsyn/poly_syn) / (fixed_nonsyn/fixed_syn)
    ))
  }, error = function(e) {
    return(list(mk_statistic = NA, mk_p_value = NA, neutrality_index = NA))
  })
}

#' Functional impact-weighted selection score
calculate_functional_selection <- function(mutations_df, weights = NULL) {
  if (is.null(weights)) {
    # Default weights based on functional categories
    weights <- list(
      "high_impact" = 3.0,
      "moderate_impact" = 2.0,
      "low_impact" = 1.0,
      "synonymous" = 0.5,
      "regulatory" = 1.5,
      "intronic" = 0.8,
      "other" = 1.0
    )
  }
  
  # Assign weights to mutations
  mutations_df$functional_weight <- sapply(mutations_df$mutation_type_category, function(x) {
    weights[[x]] %||% 1.0
  })
  
  # Add pathogenicity scores if available
  if (!all(is.na(mutations_df$cadd_score))) {
    # Normalize CADD scores to 0-1 range and add to weights
    cadd_norm <- (mutations_df$cadd_score - min(mutations_df$cadd_score, na.rm=TRUE)) / 
                 (max(mutations_df$cadd_score, na.rm=TRUE) - min(mutations_df$cadd_score, na.rm=TRUE))
    mutations_df$functional_weight <- mutations_df$functional_weight * (1 + cadd_norm)
  }
  
  # Calculate weighted selection metrics by gene and methylation class
  selection_metrics <- mutations_df %>%
    group_by(gene_symbol, methylation_class, mutation_origin) %>%
    summarise(
      total_mutations = n(),
      weighted_impact = sum(functional_weight, na.rm=TRUE),
      mean_cadd = mean(cadd_score, na.rm=TRUE),
      mean_gerp = mean(gerp_score, na.rm=TRUE),
      high_impact_count = sum(mutation_type_category == "high_impact", na.rm=TRUE),
      moderate_impact_count = sum(mutation_type_category == "moderate_impact", na.rm=TRUE),
      low_impact_count = sum(mutation_type_category == "low_impact", na.rm=TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      functional_selection_score = weighted_impact / total_mutations,
      impact_ratio = (high_impact_count + moderate_impact_count) / 
                     (low_impact_count + 1)  # Add pseudocount
    )
  
  return(selection_metrics)
}

#' Bayesian estimation of selection coefficients
bayesian_selection_estimation <- function(mutations_df, gene_lengths_df = NULL) {
  log_message("Starting Bayesian selection estimation")
  
  # Prepare data for Bayesian analysis
  gene_data <- mutations_df %>%
    group_by(gene_symbol, methylation_class, mutation_origin) %>%
    summarise(
      syn_count = sum(mutation_type_category %in% c("low_impact", "synonymous"), na.rm=TRUE),
      nonsyn_count = sum(mutation_type_category %in% c("moderate_impact", "high_impact"), na.rm=TRUE),
      total_count = n(),
      mean_cadd = mean(cadd_score, na.rm=TRUE),
      .groups = "drop"
    ) %>%
    filter(total_count >= opt$min_mutations)
  
  if (nrow(gene_data) == 0) {
    log_message("No genes meet minimum mutation threshold", "WARNING")
    return(NULL)
  }
  
  # Add gene lengths if available
  if (!is.null(gene_lengths_df)) {
    gene_data <- gene_data %>%
      left_join(gene_lengths_df, by = "gene_symbol") %>%
      mutate(
        syn_sites = coalesce(synonymous_sites, 1000),  # Default if missing
        nonsyn_sites = coalesce(nonsynonymous_sites, 3000)
      )
  } else {
    # Use approximate values based on typical gene properties
    gene_data <- gene_data %>%
      mutate(
        syn_sites = 1000,  # Approximate synonymous sites per gene
        nonsyn_sites = 3000  # Approximate non-synonymous sites per gene
      )
  }
  
  # Calculate selection coefficients for each gene/methylation/origin combination
  selection_results <- gene_data %>%
    rowwise() %>%
    mutate(
      dnds_result = list(calculate_dnds(syn_count, nonsyn_count, syn_sites, nonsyn_sites))
    ) %>%
    unnest_wider(dnds_result)
  
  return(selection_results)
}

#' Statistical comparison between methylation states
compare_methylation_selection <- function(selection_data) {
  log_message("Comparing selection between methylation states")
  
  comparison_results <- list()
  
  # Compare dN/dS ratios between methylation states
  methylation_classes <- unique(selection_data$methylation_class)
  origin_types <- unique(selection_data$mutation_origin)
  
  for (origin in origin_types) {
    origin_data <- selection_data[selection_data$mutation_origin == origin & !is.na(selection_data$dnds),]
    
    if (nrow(origin_data) < 10) next
    
    # Pairwise comparisons between methylation classes
    pairwise_comparisons <- list()
    
    for (i in 1:(length(methylation_classes)-1)) {
      for (j in (i+1):length(methylation_classes)) {
        class1 <- methylation_classes[i]
        class2 <- methylation_classes[j]
        
        data1 <- origin_data[origin_data$methylation_class == class1, "dnds"]
        data2 <- origin_data[origin_data$methylation_class == class2, "dnds"]
        
        if (length(data1) < 5 | length(data2) < 5) next
        
        # Wilcoxon test (non-parametric)
        wilcox_test <- tryCatch({
          wilcox.test(data1, data2)
        }, error = function(e) NULL)
        
        # t-test on log-transformed data
        t_test <- tryCatch({
          t.test(log(data1[data1 > 0]), log(data2[data2 > 0]))
        }, error = function(e) NULL)
        
        comparison_key <- paste(class1, "vs", class2, origin, sep="_")
        pairwise_comparisons[[comparison_key]] <- list(
          class1 = class1,
          class2 = class2,
          origin = origin,
          n1 = length(data1),
          n2 = length(data2),
          median1 = median(data1, na.rm=TRUE),
          median2 = median(data2, na.rm=TRUE),
          mean1 = mean(data1, na.rm=TRUE),
          mean2 = mean(data2, na.rm=TRUE),
          wilcox_p = ifelse(is.null(wilcox_test), NA, wilcox_test$p.value),
          t_test_p = ifelse(is.null(t_test), NA, t_test$p.value),
          effect_size = median(data1, na.rm=TRUE) - median(data2, na.rm=TRUE)
        )
      }
    }
    
    comparison_results[[origin]] <- pairwise_comparisons
  }
  
  return(comparison_results)
}

#' Compare germline vs somatic selection
compare_germline_somatic <- function(selection_data) {
  log_message("Comparing germline vs somatic selection")
  
  comparison_results <- list()
  methylation_classes <- unique(selection_data$methylation_class)
  
  for (meth_class in methylation_classes) {
    class_data <- selection_data[
      selection_data$methylation_class == meth_class & !is.na(selection_data$dnds),
    ]
    
    if (nrow(class_data) < 10) next
    
    germline_data <- class_data[class_data$mutation_origin == "germline", "dnds"]
    somatic_data <- class_data[class_data$mutation_origin == "somatic", "dnds"]
    
    if (length(germline_data) < 5 | length(somatic_data) < 5) next
    
    # Statistical tests
    wilcox_test <- tryCatch({
      wilcox.test(germline_data, somatic_data)
    }, error = function(e) NULL)
    
    t_test <- tryCatch({
      t.test(log(germline_data[germline_data > 0]), 
             log(somatic_data[somatic_data > 0]))
    }, error = function(e) NULL)
    
    comparison_results[[meth_class]] <- list(
      methylation_class = meth_class,
      germline_n = length(germline_data),
      somatic_n = length(somatic_data),
      germline_median = median(germline_data, na.rm=TRUE),
      somatic_median = median(somatic_data, na.rm=TRUE),
      germline_mean = mean(germline_data, na.rm=TRUE),
      somatic_mean = mean(somatic_data, na.rm=TRUE),
      wilcox_p = ifelse(is.null(wilcox_test), NA, wilcox_test$p.value),
      t_test_p = ifelse(is.null(t_test), NA, t_test$p.value),
      effect_size = median(somatic_data, na.rm=TRUE) - median(germline_data, na.rm=TRUE)
    )
  }
  
  return(comparison_results)
}

#' Mixed effects model for selection analysis
mixed_effects_selection_model <- function(mutations_df) {
  log_message("Fitting mixed effects models for selection analysis")
  
  # Prepare data for modeling
  model_data <- mutations_df %>%
    filter(!is.na(cadd_score) | !is.na(gerp_score)) %>%
    mutate(
      log_cadd = log(pmax(cadd_score, 0.1, na.rm=TRUE)),  # Add small constant
      high_impact = as.numeric(mutation_type_category %in% c("high_impact", "moderate_impact")),
      methylation_numeric = case_when(
        methylation_class == "hypomethylated" ~ -1,
        methylation_class == "intermediate" ~ 0,
        methylation_class == "hypermethylated" ~ 1,
        TRUE ~ 0
      ),
      origin_numeric = ifelse(mutation_origin == "germline", 0, 1)
    )
  
  if (nrow(model_data) < 100) {
    log_message("Insufficient data for mixed effects modeling", "WARNING")
    return(NULL)
  }
  
  model_results <- list()
  
  # Model 1: CADD score as response
  if (sum(!is.na(model_data$log_cadd)) > 50) {
    tryCatch({
      cadd_model <- lmer(
        log_cadd ~ methylation_class * mutation_origin + (1 | gene_symbol),
        data = model_data,
        REML = FALSE
      )
      model_results[["cadd_model"]] <- cadd_model
    }, error = function(e) {
      log_message(paste("CADD model failed:", e$message), "WARNING")
    })
  }
  
  # Model 2: High impact probability
  tryCatch({
    impact_model <- glmer(
      high_impact ~ methylation_class * mutation_origin + (1 | gene_symbol),
      data = model_data,
      family = binomial,
      control = glmerControl(optimizer = "bobyqa")
    )
    model_results[["impact_model"]] <- impact_model
  }, error = function(e) {
    log_message(paste("Impact model failed:", e$message), "WARNING")
  })
  
  return(model_results)
}

#' Bootstrap confidence intervals for selection metrics
bootstrap_selection_ci <- function(mutations_df, metric_function, n_bootstrap = 1000) {
  log_message(paste("Calculating bootstrap confidence intervals with", n_bootstrap, "replicates"))
  
  # Original metric calculation
  original_metric <- metric_function(mutations_df)
  
  # Bootstrap replicates
  bootstrap_metrics <- mclapply(1:n_bootstrap, function(i) {
    boot_sample <- mutations_df[sample(nrow(mutations_df), replace = TRUE), ]
    metric_function(boot_sample)
  }, mc.cores = min(opt$cores, detectCores()))
  
  # Calculate confidence intervals
  bootstrap_values <- unlist(bootstrap_metrics)
  ci_lower <- quantile(bootstrap_values, 0.025, na.rm = TRUE)
  ci_upper <- quantile(bootstrap_values, 0.975, na.rm = TRUE)
  
  return(list(
    original = original_metric,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    bootstrap_values = bootstrap_values
  ))
}

#' Create comprehensive visualizations
create_selection_plots <- function(selection_data, comparison_results, output_prefix) {
  log_message("Creating selection visualization plots")
  
  # Plot 1: dN/dS ratios across methylation states
  p1 <- ggplot(selection_data[!is.na(selection_data$dnds) & selection_data$dnds < 10,], 
               aes(x = methylation_class, y = dnds, fill = mutation_origin)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    scale_y_log10() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    facet_wrap(~ mutation_origin) +
    theme_minimal() +
    labs(
      title = "dN/dS Ratios Across Methylation States",
      x = "Methylation Class",
      y = "dN/dS Ratio (log scale)",
      fill = "Mutation Origin"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Plot 2: Functional selection scores
  if ("functional_selection_score" %in% colnames(selection_data)) {
    p2 <- ggplot(selection_data, aes(x = methylation_class, y = functional_selection_score, 
                                     color = mutation_origin)) +
      geom_boxplot() +
      geom_point(position = position_jitterdodge(), alpha = 0.6) +
      theme_minimal() +
      labs(
        title = "Functional Selection Scores",
        x = "Methylation Class",
        y = "Functional Selection Score",
        color = "Mutation Origin"
      )
  } else {
    p2 <- ggplot() + theme_void()
  }
  
  # Plot 3: Mutation burden comparison
  burden_data <- selection_data %>%
    group_by(methylation_class, mutation_origin) %>%
    summarise(
      total_mutations = sum(total_mutations, na.rm = TRUE),
      mean_dnds = mean(dnds, na.rm = TRUE),
      .groups = "drop"
    )
  
  p3 <- ggplot(burden_data, aes(x = methylation_class, y = total_mutations, 
                                fill = mutation_origin)) +
    geom_col(position = "dodge") +
    theme_minimal() +
    labs(
      title = "Mutation Burden by Methylation State",
      x = "Methylation Class", 
      y = "Total Mutations",
      fill = "Mutation Origin"
    )
  
  # Plot 4: Selection comparison heatmap
  if (length(comparison_results) > 0) {
    # Flatten comparison results for plotting
    comp_df <- map_dfr(comparison_results, function(origin_comps) {
      map_dfr(origin_comps, as.data.frame)
    }, .id = "origin")
    
    if (nrow(comp_df) > 0) {
      p4 <- ggplot(comp_df, aes(x = class1, y = class2, fill = -log10(wilcox_p + 1e-10))) +
        geom_tile() +
        geom_text(aes(label = sprintf("%.2f", effect_size)), size = 3) +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1.3) +
        facet_wrap(~ origin) +
        theme_minimal() +
        labs(
          title = "Selection Differences Between Methylation States",
          x = "Methylation Class 1",
          y = "Methylation Class 2", 
          fill = "-log10(p-value)"
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    } else {
      p4 <- ggplot() + theme_void()
    }
  } else {
    p4 <- ggplot() + theme_void()
  }
  
  # Combine plots
  combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2)
  
  # Save plot
  ggsave(paste0(output_prefix, "_selection_plots.png"), combined_plot, 
         width = 12, height = 10, dpi = 300)
  
  log_message("Selection plots saved")
}

#' Generate comprehensive report
generate_selection_report <- function(selection_data, comparison_results, 
                                    germline_somatic_results, model_results, 
                                    output_prefix) {
  log_message("Generating comprehensive selection analysis report")
  
  report_file <- paste0(output_prefix, "_selection_report.txt")
  
  sink(report_file)
  
  cat("COMPREHENSIVE SELECTION ANALYSIS REPORT\n")
  cat("=====================================\n\n")
  cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  # Summary statistics
  cat("SUMMARY STATISTICS\n")
  cat("-----------------\n")
  cat("Total genes analyzed:", length(unique(selection_data$gene_symbol)), "\n")
  cat("Total mutation-gene combinations:", nrow(selection_data), "\n")
  cat("Methylation classes:", paste(unique(selection_data$methylation_class), collapse = ", "), "\n")
  cat("Mutation origins:", paste(unique(selection_data$mutation_origin), collapse = ", "), "\n\n")
  
  # dN/dS distribution
  cat("dN/dS RATIO DISTRIBUTION\n")
  cat("----------------------\n")
  dnds_summary <- selection_data %>%
    filter(!is.na(dnds)) %>%
    group_by(methylation_class, mutation_origin) %>%
    summarise(
      n = n(),
      median_dnds = median(dnds),
      mean_dnds = mean(dnds),
      q25 = quantile(dnds, 0.25),
      q75 = quantile(dnds, 0.75),
      .groups = "drop"
    )
  
  print(dnds_summary)
  cat("\n")
  
  # Statistical comparisons
  cat("METHYLATION STATE COMPARISONS\n")
  cat("---------------------------\n")
  if (length(comparison_results) > 0) {
    for (origin in names(comparison_results)) {
      cat("Origin:", origin, "\n")
      for (comp_name in names(comparison_results[[origin]])) {
        comp <- comparison_results[[origin]][[comp_name]]
        cat(sprintf("  %s vs %s: p = %.4f, effect = %.4f\n", 
                   comp$class1, comp$class2, comp$wilcox_p, comp$effect_size))
      }
      cat("\n")
    }
  }
  
  cat("GERMLINE VS SOMATIC COMPARISONS\n")
  cat("-----------------------------\n")
  if (length(germline_somatic_results) > 0) {
    for (meth_class in names(germline_somatic_results)) {
      comp <- germline_somatic_results[[meth_class]]
      cat(sprintf("%s: p = %.4f, effect = %.4f (somatic - germline)\n",
                 meth_class, comp$wilcox_p, comp$effect_size))
    }
  }
  cat("\n")
  
  # Model results
  if (length(model_results) > 0) {
    cat("MIXED EFFECTS MODEL RESULTS\n")
    cat("-------------------------\n")
    for (model_name in names(model_results)) {
      cat("Model:", model_name, "\n")
      tryCatch({
        model_summary <- summary(model_results[[model_name]])
        print(model_summary$coefficients)
      }, error = function(e) {
        cat("Model summary unavailable\n")
      })
      cat("\n")
    }
  }
  
  # Interpretation
  cat("BIOLOGICAL INTERPRETATION\n")
  cat("------------------------\n")
  cat("1. Selection in hypomethylated regions:\n")
  hypo_germline <- dnds_summary$median_dnds[dnds_summary$methylation_class == "hypomethylated" & 
                                           dnds_summary$mutation_origin == "germline"]
  hypo_somatic <- dnds_summary$median_dnds[dnds_summary$methylation_class == "hypomethylated" & 
                                          dnds_summary$mutation_origin == "somatic"]
  
  if (length(hypo_germline) > 0) {
    cat(sprintf("   Germline median dN/dS = %.3f ", hypo_germline))
    if (hypo_germline < 1) cat("(purifying selection)")
    else if (hypo_germline > 1) cat("(positive selection)")
    else cat("(neutral)")
    cat("\n")
  }
  
  if (length(hypo_somatic) > 0) {
    cat(sprintf("   Somatic median dN/dS = %.3f ", hypo_somatic))
    if (hypo_somatic < 1) cat("(purifying selection)")
    else if (hypo_somatic > 1) cat("(positive selection)")
    else cat("(neutral)")
    cat("\n")
  }
  
  cat("\n2. Selection in hypermethylated regions:\n")
  hyper_germline <- dnds_summary$median_dnds[dnds_summary$methylation_class == "hypermethylated" & 
                                            dnds_summary$mutation_origin == "germline"]
  hyper_somatic <- dnds_summary$median_dnds[dnds_summary$methylation_class == "hypermethylated" & 
                                           dnds_summary$mutation_origin == "somatic"]
  
  if (length(hyper_germline) > 0) {
    cat(sprintf("   Germline median dN/dS = %.3f ", hyper_germline))
    if (hyper_germline < 1) cat("(purifying selection)")
    else if (hyper_germline > 1) cat("(positive selection)")
    else cat("(neutral)")
    cat("\n")
  }
  
  if (length(hyper_somatic) > 0) {
    cat(sprintf("   Somatic median dN/dS = %.3f ", hyper_somatic))
    if (hyper_somatic < 1) cat("(purifying selection)")
    else if (hyper_somatic > 1) cat("(positive selection)")
    else cat("(neutral)")
    cat("\n")
  }
  
  sink()
  
  log_message("Selection analysis report generated")
}

# Main Analysis Pipeline
######################

main <- function() {
  log_message("Starting comprehensive selection analysis")
  
  # Validate inputs
  if (is.null(opt$mutations)) {
    stop("Must provide annotated mutations file")
  }
  
  if (!file.exists(opt$mutations)) {
    stop(paste("Mutations file not found:", opt$mutations))
  }
  
  # Load data
  log_message("Loading annotated mutations data")
  mutations_df <- fread(opt$mutations, sep = "\t", header = TRUE)
  
  log_message(paste("Loaded", nrow(mutations_df), "mutations"))
  
  # Load gene lengths if provided
  gene_lengths_df <- NULL
  if (!is.null(opt$`gene-lengths`) && file.exists(opt$`gene-lengths`)) {
    log_message("Loading gene length data")
    gene_lengths_df <- fread(opt$`gene-lengths`, sep = "\t", header = TRUE)
  }
  
  # Filter data for analysis
  analysis_data <- mutations_df %>%
    filter(
      !is.na(methylation_class),
      !is.na(mutation_origin),
      methylation_class %in% c("hypomethylated", "intermediate", "hypermethylated"),
      mutation_origin %in% c("germline", "somatic")
    ) %>%
    # Remove genes with too few mutations
    group_by(gene_symbol) %>%
    filter(n() >= opt$min_mutations) %>%
    ungroup()
  
  log_message(paste("Filtered to", nrow(analysis_data), "mutations for analysis"))
  
  if (nrow(analysis_data) < 10) {
    stop("Insufficient data for selection analysis")
  }
  
  # Calculate selection metrics
  log_message("Calculating selection metrics")
  
  # 1. Bayesian dN/dS estimation
  selection_data <- bayesian_selection_estimation(analysis_data, gene_lengths_df)
  
  if (is.null(selection_data) || nrow(selection_data) == 0) {
    stop("Failed to calculate selection metrics")
  }
  
  # 2. Functional selection scores
  functional_metrics <- calculate_functional_selection(analysis_data)
  
  # Merge functional metrics with selection data
  selection_data <- selection_data %>%
    left_join(functional_metrics, 
              by = c("gene_symbol", "methylation_class", "mutation_origin"))
  
  # 3. Statistical comparisons
  methylation_comparisons <- compare_methylation_selection(selection_data)
  germline_somatic_comparisons <- compare_germline_somatic(selection_data)
  
  # 4. Mixed effects modeling
  model_results <- mixed_effects_selection_model(analysis_data)
  
  # Save intermediate results
  log_message("Saving intermediate results")
  write.table(selection_data, paste0(opt$output, "_selection_metrics.tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Create visualizations
  create_selection_plots(selection_data, methylation_comparisons, opt$output)
  
  # Generate comprehensive report
  generate_selection_report(
    selection_data, 
    methylation_comparisons,
    germline_somatic_comparisons,
    model_results,
    opt$output
  )
  
  # Save final results as RData for further analysis
  save(selection_data, methylation_comparisons, germline_somatic_comparisons, 
       model_results, analysis_data,
       file = paste0(opt$output, "_selection_analysis.RData"))
  
  log_message("Selection analysis completed successfully")
}

# Run main analysis
if (!interactive()) {
  tryCatch({
    main()
  }, error = function(e) {
    log_message(paste("ERROR:", e$message), "ERROR")
    stop(e)
  }, finally = {
    close(log_connection)
  })
}
