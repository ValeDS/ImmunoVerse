# ============================================================
# 11_family_level_differential_abundance.R
# Matched differential-abundance analysis of immune-state families
# ============================================================

suppressPackageStartupMessages({
  library(speckle)
  library(limma)
  library(dplyr)
  library(tidyr)
})

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
# Load cell- and sample-level metadata
# ------------------------------------------------------------

meta_full <- readRDS(
  "data/processed/GSE144744_tcells_full_state_metadata.rds"
)

original_sample_metadata <- read.csv(
  "data/raw/GSE144744/GSE144744_metadata_per_sample.csv.gz",
  stringsAsFactors = FALSE
)

matching_metadata <- original_sample_metadata %>%
  select(
    sample,
    donor,
    group,
    cohort,
    match,
    batch_pair,
    age_sampling,
    sex
  ) %>%
  distinct()

stopifnot(
  !anyDuplicated(matching_metadata$sample)
)

sample_info <- meta_full %>%
  distinct(sample, donor, group, cohort) %>%
  left_join(
    matching_metadata,
    by = c(
      "sample",
      "donor",
      "group",
      "cohort"
    )
  )

stopifnot(
  !anyNA(sample_info$batch_pair)
)

# ------------------------------------------------------------
# A priori state-to-family mapping
# Defined from marker profiles and functional annotations,
# independently of the differential-abundance results.
# ------------------------------------------------------------

family_lookup <- tibble::tribble(
  ~immune_state_id, ~immune_family,
  "0",  "Cytotoxic",
  "1",  "Memory CD4-like",
  "2",  "Memory CD4-like",
  "3",  "Memory CD4-like",
  "4",  "Memory CD4-like",
  "5",  "Memory CD4-like",
  "6",  "Cytotoxic",
  "7",  "Cytotoxic",
  "8",  "Activated/stress-associated",
  "9",  "Innate-like",
  "10", "Activated/stress-associated",
  "11", "Donor-associated",
  "12", "Activated/stress-associated",
  "13", "Low-abundance candidate",
  "14", "Memory CD4-like",
  "15", "Regulatory-like",
  "16", "Cytotoxic",
  "17", "Low-abundance candidate",
  "18", "Low-abundance candidate",
  "19", "Low-abundance candidate",
  "20", "Memory CD4-like",
  "21", "Donor-associated",
  "22", "Activated/stress-associated",
  "23", "Memory CD4-like",
  "24", "Activated/stress-associated",
  "25", "Cytotoxic",
  "26", "Cytotoxic",
  "27", "Cytotoxic",
  "28", "Activated/stress-associated",
  "29", "Memory CD4-like",
  "30", "Activated/stress-associated",
  "31", "Low-abundance candidate",
  "32", "Low-abundance candidate",
  "33", "Activated/stress-associated",
  "34", "Cytotoxic",
  "35", "Low-abundance candidate"
)

expected_states <- as.character(0:35)

stopifnot(
  setequal(
    family_lookup$immune_state_id,
    expected_states
  ),
  !anyDuplicated(family_lookup$immune_state_id)
)

family_levels <- c(
  "Cytotoxic",
  "Memory CD4-like",
  "Activated/stress-associated",
  "Regulatory-like",
  "Innate-like",
  "Donor-associated",
  "Low-abundance candidate"
)

meta_full <- meta_full %>%
  mutate(
    immune_state_id = as.character(
      immune_state_id
    )
  ) %>%
  left_join(
    family_lookup,
    by = "immune_state_id"
  )

stopifnot(
  !anyNA(meta_full$immune_family)
)

meta_full$immune_family <- factor(
  meta_full$immune_family,
  levels = family_levels
)

write.csv(
  family_lookup,
  paste0(
    "results/tables/propeller/",
    "GSE144744_state_to_family_mapping.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Matched Propeller analysis at family level
# ------------------------------------------------------------

run_matched_family_propeller <- function(
    cell_metadata,
    sample_metadata,
    cohort_name,
    reference_group,
    test_group,
    block_column,
    contrast_label,
    transform_method = "logit") {
  
  analysis_samples <- sample_metadata %>%
    filter(
      cohort == cohort_name,
      group %in% c(
        reference_group,
        test_group
      )
    ) %>%
    mutate(
      .block = .data[[block_column]]
    ) %>%
    group_by(.block) %>%
    filter(
      n() == 2,
      n_distinct(group) == 2,
      all(
        c(reference_group, test_group) %in%
          group
      )
    ) %>%
    ungroup() %>%
    arrange(.block, group)
  
  if (nrow(analysis_samples) == 0) {
    stop(
      "No complete matched pairs found for ",
      contrast_label
    )
  }
  
  group_counts <- table(
    analysis_samples$group
  )
  
  if (
    group_counts[[reference_group]] !=
    group_counts[[test_group]]
  ) {
    stop(
      "Unequal groups after matched-pair selection: ",
      contrast_label
    )
  }
  
  cell_subset <- cell_metadata %>%
    semi_join(
      analysis_samples %>%
        select(sample),
      by = "sample"
    )
  
  cell_subset$immune_family <- factor(
    cell_subset$immune_family,
    levels = family_levels
  )
  
  propeller_input <- getTransformedProps(
    clusters = cell_subset$immune_family,
    sample = cell_subset$sample,
    transform = transform_method
  )
  
  sample_order <- colnames(
    propeller_input$TransformedProps
  )
  
  design_metadata <- analysis_samples[
    match(
      sample_order,
      analysis_samples$sample
    ),
    ,
    drop = FALSE
  ]
  
  stopifnot(
    !anyNA(design_metadata$sample),
    identical(
      sample_order,
      design_metadata$sample
    )
  )
  
  design_metadata$group_factor <- factor(
    design_metadata$group,
    levels = c(
      reference_group,
      test_group
    )
  )
  
  design_metadata$block_factor <- factor(
    design_metadata$.block
  )
  
  design <- model.matrix(
    ~ 0 + group_factor + block_factor,
    data = design_metadata
  )
  
  colnames(design) <- sub(
    "^group_factor",
    "group",
    colnames(design)
  )
  
  colnames(design) <- sub(
    "^block_factor",
    "block",
    colnames(design)
  )
  
  colnames(design) <- make.names(
    colnames(design),
    unique = TRUE
  )
  
  if (qr(design)$rank != ncol(design)) {
    stop(
      "Design matrix is not full rank: ",
      contrast_label
    )
  }
  
  contrast_expression <- paste0(
    "group",
    test_group,
    "-group",
    reference_group
  )
  
  contrast_matrix <- makeContrasts(
    contrasts = contrast_expression,
    levels = design
  )
  
  result <- propeller.ttest(
    prop.list = propeller_input,
    design = design,
    contrasts = contrast_matrix,
    robust = TRUE,
    trend = FALSE,
    sort = FALSE
  )
  
  result <- as.data.frame(result)
  result$immune_family <- rownames(result)
  
  proportion_matrix <-
    propeller_input$Proportions
  
  reference_columns <-
    design_metadata$group == reference_group
  
  test_columns <-
    design_metadata$group == test_group
  
  effect_sizes <- data.frame(
    immune_family = rownames(
      proportion_matrix
    ),
    
    mean_reference_percent =
      100 * rowMeans(
        proportion_matrix[
          ,
          reference_columns,
          drop = FALSE
        ]
      ),
    
    mean_test_percent =
      100 * rowMeans(
        proportion_matrix[
          ,
          test_columns,
          drop = FALSE
        ]
      ),
    
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      difference_percentage_points =
        mean_test_percent -
        mean_reference_percent,
      
      proportion_ratio = ifelse(
        mean_reference_percent > 0,
        mean_test_percent /
          mean_reference_percent,
        NA_real_
      )
    )
  
  result <- effect_sizes %>%
    left_join(
      result,
      by = "immune_family"
    ) %>%
    mutate(
      contrast = contrast_label,
      cohort = cohort_name,
      reference_group = reference_group,
      test_group = test_group,
      block_variable = block_column,
      n_pairs = nrow(
        analysis_samples
      ) / 2,
      transformation = transform_method
    ) %>%
    arrange(FDR, P.Value)
  
  list(
    results = result,
    proportions = propeller_input,
    design = design,
    design_metadata = design_metadata,
    contrast = contrast_matrix
  )
}

# ------------------------------------------------------------
# Wrapper for the four matched contrasts
# ------------------------------------------------------------

run_four_family_contrasts <- function(
    cell_metadata,
    transform_method) {
  
  list(
    RRMS = run_matched_family_propeller(
      cell_metadata = cell_metadata,
      sample_metadata = sample_info,
      cohort_name = "RRMS_HI",
      reference_group = "HI2",
      test_group = "MS2",
      block_column = "batch_pair",
      contrast_label = "MS2_vs_HI2",
      transform_method = transform_method
    ),
    
    PPMS = run_matched_family_propeller(
      cell_metadata = cell_metadata,
      sample_metadata = sample_info,
      cohort_name = "PPMS_HI",
      reference_group = "HI3",
      test_group = "PPMS",
      block_column = "batch_pair",
      contrast_label = "PPMS_vs_HI3",
      transform_method = transform_method
    ),
    
    NAT_disease = run_matched_family_propeller(
      cell_metadata = cell_metadata,
      sample_metadata = sample_info,
      cohort_name = "NAT_HI",
      reference_group = "HI1",
      test_group = "MS1",
      block_column = "batch_pair",
      contrast_label = "MS1_vs_HI1",
      transform_method = transform_method
    ),
    
    NAT_treatment = run_matched_family_propeller(
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

extract_family_results <- function(
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
# Primary and sensitivity analyses
# ------------------------------------------------------------

models_full_logit <- run_four_family_contrasts(
  cell_metadata = meta_full,
  transform_method = "logit"
)

models_full_asin <- run_four_family_contrasts(
  cell_metadata = meta_full,
  transform_method = "asin"
)

meta_confident <- meta_full %>%
  filter(
    state_mapping_origin == "reference" |
      (
        state_mapping_origin == "transferred" &
          state_prediction_score >= 0.50
      )
  )

models_confident_logit <- run_four_family_contrasts(
  cell_metadata = meta_confident,
  transform_method = "logit"
)

meta_reference_only <- meta_full %>%
  filter(
    state_mapping_origin == "reference"
  )

models_reference_logit <- run_four_family_contrasts(
  cell_metadata = meta_reference_only,
  transform_method = "logit"
)

results_full_logit <- extract_family_results(
  models_full_logit,
  "full_logit"
)

results_full_asin <- extract_family_results(
  models_full_asin,
  "full_asin"
)

results_confident_logit <- extract_family_results(
  models_confident_logit,
  "confidence_050_logit"
)

results_reference_logit <- extract_family_results(
  models_reference_logit,
  "reference_only_logit"
)

all_family_results <- bind_rows(
  results_full_logit,
  results_full_asin,
  results_confident_logit,
  results_reference_logit
)

# ------------------------------------------------------------
# Save models and complete results
# ------------------------------------------------------------

write.csv(
  all_family_results,
  paste0(
    "results/tables/propeller/",
    "GSE144744_family_level_all_analyses.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    full_logit = models_full_logit,
    full_asin = models_full_asin,
    confidence_050_logit =
      models_confident_logit,
    reference_only_logit =
      models_reference_logit
  ),
  paste0(
    "results/models/",
    "GSE144744_family_level_Propeller_models.rds"
  )
)

# ------------------------------------------------------------
# Summaries
# ------------------------------------------------------------

family_significance_summary <- all_family_results %>%
  group_by(
    analysis_set,
    contrast
  ) %>%
  summarise(
    tested_families = n(),
    nominal_p_below_005 = sum(
      P.Value < 0.05
    ),
    FDR_below_005 = sum(
      FDR < 0.05
    ),
    .groups = "drop"
  )

primary_family_results <- results_full_logit %>%
  arrange(
    contrast,
    FDR,
    P.Value
  )

primary_family_summary <- primary_family_results %>%
  select(
    contrast,
    immune_family,
    mean_reference_percent,
    mean_test_percent,
    difference_percentage_points,
    proportion_ratio,
    P.Value,
    FDR
  )

primary_effects <- results_full_logit %>%
  select(
    contrast,
    immune_family,
    primary_difference =
      difference_percentage_points
  )

family_robustness_summary <- all_family_results %>%
  filter(
    analysis_set != "full_logit"
  ) %>%
  left_join(
    primary_effects,
    by = c(
      "contrast",
      "immune_family"
    )
  ) %>%
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
    .groups = "drop"
  )

write.csv(
  family_significance_summary,
  paste0(
    "results/tables/propeller/",
    "GSE144744_family_level_significance_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  primary_family_summary,
  paste0(
    "results/tables/propeller/",
    "GSE144744_family_level_primary_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  family_robustness_summary,
  paste0(
    "results/tables/propeller/",
    "GSE144744_family_level_robustness_summary.csv"
  ),
  row.names = FALSE
)

print(
  tibble::as_tibble(
    family_significance_summary
  ),
  n = Inf,
  width = Inf
)

print(
  tibble::as_tibble(
    primary_family_summary
  ),
  n = Inf,
  width = Inf
)

print(
  tibble::as_tibble(
    family_robustness_summary
  ),
  n = Inf,
  width = Inf
)
