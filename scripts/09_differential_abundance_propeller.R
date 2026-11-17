# ============================================================
# 09_differential_abundance_propeller.R
# Matched differential-abundance analysis of immune states
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
# Load cell-level state metadata
# ------------------------------------------------------------

meta_full <- readRDS(
  "data/processed/GSE144744_tcells_full_state_metadata.rds"
)

meta_full$immune_state_id <- factor(
  meta_full$immune_state_id,
  levels = as.character(0:35)
)

# ------------------------------------------------------------
# Reconstruct sample-level matching information
# ------------------------------------------------------------

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
# Matched Propeller analysis
# ------------------------------------------------------------

run_matched_propeller <- function(
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
  
  cell_subset$immune_state_id <- factor(
    cell_subset$immune_state_id,
    levels = as.character(0:35)
  )
  
  propeller_input <- getTransformedProps(
    clusters = cell_subset$immune_state_id,
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
  result$immune_state_id <- rownames(result)
  
  proportion_matrix <-
    propeller_input$Proportions
  
  reference_columns <-
    design_metadata$group == reference_group
  
  test_columns <-
    design_metadata$group == test_group
  
  effect_sizes <- data.frame(
    immune_state_id = rownames(
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
      by = "immune_state_id"
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
# Primary analyses
# Positive effects indicate higher abundance in disease/treatment
# ------------------------------------------------------------

da_rrms <- run_matched_propeller(
  cell_metadata = meta_full,
  sample_metadata = sample_info,
  cohort_name = "RRMS_HI",
  reference_group = "HI2",
  test_group = "MS2",
  block_column = "batch_pair",
  contrast_label = "MS2_vs_HI2"
)

da_ppms <- run_matched_propeller(
  cell_metadata = meta_full,
  sample_metadata = sample_info,
  cohort_name = "PPMS_HI",
  reference_group = "HI3",
  test_group = "PPMS",
  block_column = "batch_pair",
  contrast_label = "PPMS_vs_HI3"
)

da_nat_disease <- run_matched_propeller(
  cell_metadata = meta_full,
  sample_metadata = sample_info,
  cohort_name = "NAT_HI",
  reference_group = "HI1",
  test_group = "MS1",
  block_column = "batch_pair",
  contrast_label = "MS1_vs_HI1"
)

da_nat_treatment <- run_matched_propeller(
  cell_metadata = meta_full,
  sample_metadata = sample_info,
  cohort_name = "NAT_HI",
  reference_group = "MS1",
  test_group = "MS1_nat",
  block_column = "donor",
  contrast_label = "MS1nat_vs_MS1"
)

# ------------------------------------------------------------
# Combine and save
# ------------------------------------------------------------

propeller_results <- bind_rows(
  da_rrms$results,
  da_ppms$results,
  da_nat_disease$results,
  da_nat_treatment$results
)

write.csv(
  propeller_results,
  "results/tables/propeller/GSE144744_Propeller_all_primary_contrasts.csv",
  row.names = FALSE
)

saveRDS(
  list(
    RRMS = da_rrms,
    PPMS = da_ppms,
    NAT_disease = da_nat_disease,
    NAT_treatment = da_nat_treatment
  ),
  "results/models/GSE144744_Propeller_primary_models.rds"
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

significance_summary <- propeller_results %>%
  group_by(contrast) %>%
  summarise(
    tested_states = n(),
    nominal_p_below_005 = sum(
      P.Value < 0.05
    ),
    FDR_below_005 = sum(
      FDR < 0.05
    ),
    .groups = "drop"
  )

top_results <- propeller_results %>%
  group_by(contrast) %>%
  slice_min(
    order_by = FDR,
    n = 5,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  select(
    contrast,
    immune_state_id,
    mean_reference_percent,
    mean_test_percent,
    difference_percentage_points,
    proportion_ratio,
    P.Value,
    FDR
  )

print(significance_summary)
print(top_results, n = 20)