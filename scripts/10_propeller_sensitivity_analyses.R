# ============================================================
# 10_propeller_sensitivity_analyses.R
# Sensitivity analyses for immune-state differential abundance
# ============================================================

suppressPackageStartupMessages({
  library(speckle)
  library(limma)
  library(dplyr)
})

# Load primary workflow if required
if (
  !exists("run_matched_propeller") ||
  !exists("propeller_results") ||
  !exists("sample_info") ||
  !exists("meta_full")
) {
  source(
    "scripts/09_differential_abundance_propeller.R"
  )
}

dir.create(
  "results/tables/propeller",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/models",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Wrapper for the four primary contrasts
# ------------------------------------------------------------

run_four_contrasts <- function(
    cell_metadata,
    transform_method) {
  
  list(
    RRMS = run_matched_propeller(
      cell_metadata = cell_metadata,
      sample_metadata = sample_info,
      cohort_name = "RRMS_HI",
      reference_group = "HI2",
      test_group = "MS2",
      block_column = "batch_pair",
      contrast_label = "MS2_vs_HI2",
      transform_method = transform_method
    ),
    
    PPMS = run_matched_propeller(
      cell_metadata = cell_metadata,
      sample_metadata = sample_info,
      cohort_name = "PPMS_HI",
      reference_group = "HI3",
      test_group = "PPMS",
      block_column = "batch_pair",
      contrast_label = "PPMS_vs_HI3",
      transform_method = transform_method
    ),
    
    NAT_disease = run_matched_propeller(
      cell_metadata = cell_metadata,
      sample_metadata = sample_info,
      cohort_name = "NAT_HI",
      reference_group = "HI1",
      test_group = "MS1",
      block_column = "batch_pair",
      contrast_label = "MS1_vs_HI1",
      transform_method = transform_method
    ),
    
    NAT_treatment = run_matched_propeller(
      cell_metadata = cell_metadata,
      sample_metadata = sample_info,
      cohort_name = "NAT_HI",
      reference_group = "MS1",
      test_group = "MS1_nat",
      block_column = "donor",
      contrast_label = "MS1nat_vs_MS1",
      transform_method = transform_method
    )
  )
}

extract_results <- function(
    model_list,
    analysis_label) {
  
  bind_rows(
    lapply(
      model_list,
      function(x) x$results
    )
  ) %>%
    mutate(
      analysis_set = analysis_label
    )
}

# ------------------------------------------------------------
# 1. Full mapped dataset with arcsin transformation
# ------------------------------------------------------------

models_full_asin <- run_four_contrasts(
  cell_metadata = meta_full,
  transform_method = "asin"
)

results_full_asin <- extract_results(
  models_full_asin,
  "full_asin"
)

# ------------------------------------------------------------
# 2. Conservative mapped dataset
# Reference cells + transferred cells with score >= 0.50
# ------------------------------------------------------------

meta_confident <- meta_full %>%
  filter(
    state_mapping_origin == "reference" |
      (
        state_mapping_origin == "transferred" &
          state_prediction_score >= 0.50
      )
  )

models_confident_logit <- run_four_contrasts(
  cell_metadata = meta_confident,
  transform_method = "logit"
)

results_confident_logit <- extract_results(
  models_confident_logit,
  "confidence_050_logit"
)

# ------------------------------------------------------------
# 3. Original reference cells only
# ------------------------------------------------------------

meta_reference_only <- meta_full %>%
  filter(
    state_mapping_origin == "reference"
  )

models_reference_logit <- run_four_contrasts(
  cell_metadata = meta_reference_only,
  transform_method = "logit"
)

results_reference_logit <- extract_results(
  models_reference_logit,
  "reference_only_logit"
)

# ------------------------------------------------------------
# Combine primary and sensitivity analyses
# ------------------------------------------------------------

primary_results <- propeller_results %>%
  mutate(
    analysis_set = "full_logit"
  )

all_sensitivity_results <- bind_rows(
  primary_results,
  results_full_asin,
  results_confident_logit,
  results_reference_logit
)

write.csv(
  all_sensitivity_results,
  paste0(
    "results/tables/propeller/",
    "GSE144744_Propeller_all_sensitivity_analyses.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    full_asin = models_full_asin,
    confidence_050_logit =
      models_confident_logit,
    reference_only_logit =
      models_reference_logit
  ),
  paste0(
    "results/models/",
    "GSE144744_Propeller_sensitivity_models.rds"
  )
)

# ------------------------------------------------------------
# Robustness relative to the primary full-logit analysis
# ------------------------------------------------------------

primary_effects <- primary_results %>%
  select(
    contrast,
    immune_state_id,
    primary_difference =
      difference_percentage_points
  )

robustness_detail <- all_sensitivity_results %>%
  filter(
    analysis_set != "full_logit"
  ) %>%
  left_join(
    primary_effects,
    by = c(
      "contrast",
      "immune_state_id"
    )
  )

robustness_summary <- robustness_detail %>%
  group_by(
    analysis_set,
    contrast
  ) %>%
  summarise(
    pearson_effect_correlation = cor(
      primary_difference,
      difference_percentage_points,
      method = "pearson"
    ),
    
    spearman_effect_correlation = cor(
      primary_difference,
      difference_percentage_points,
      method = "spearman"
    ),
    
    direction_agreement_percent =
      100 * mean(
        sign(primary_difference) ==
          sign(difference_percentage_points)
      ),
    
    nominal_p_below_005 = sum(
      P.Value < 0.05
    ),
    
    FDR_below_005 = sum(
      FDR < 0.05
    ),
    
    .groups = "drop"
  )

write.csv(
  robustness_summary,
  paste0(
    "results/tables/propeller/",
    "GSE144744_Propeller_sensitivity_summary.csv"
  ),
  row.names = FALSE
)

# State 23 across all analyses
state23_sensitivity <- all_sensitivity_results %>%
  filter(
    immune_state_id == "23"
  ) %>%
  select(
    analysis_set,
    contrast,
    mean_reference_percent,
    mean_test_percent,
    difference_percentage_points,
    proportion_ratio,
    P.Value,
    FDR
  ) %>%
  arrange(
    contrast,
    analysis_set
  )

write.csv(
  state23_sensitivity,
  paste0(
    "results/tables/propeller/",
    "GSE144744_State23_sensitivity_results.csv"
  ),
  row.names = FALSE
)

print(
  tibble::as_tibble(robustness_summary),
  n = Inf,
  width = Inf
)

print(
  tibble::as_tibble(state23_sensitivity),
  n = Inf,
  width = Inf
)