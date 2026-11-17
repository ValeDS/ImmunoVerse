# ============================================================
# 12_immune_state_portability_audit.R
# Audit of recurrence, donor concentration, cohort distribution,
# and label-transfer quality for the 36 candidate immune states
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

dir.create(
  "results/tables/state_audit",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Load metadata and existing annotations
# ------------------------------------------------------------

meta_full <- readRDS(
  "data/processed/GSE144744_tcells_full_state_metadata.rds"
) %>%
  mutate(
    immune_state_id = as.character(
      immune_state_id
    )
  )

state_levels <- as.character(0:35)

meta_full$immune_state_id <- factor(
  meta_full$immune_state_id,
  levels = state_levels
)

sample_info <- meta_full %>%
  distinct(
    sample,
    donor,
    group,
    cohort
  )

stopifnot(
  !anyDuplicated(sample_info$sample)
)

atlas_annotations <- read.csv(
  paste0(
    "results/tables/",
    "ImmunoVerse_Immune_State_Atlas_Table1_filled.csv"
  ),
  stringsAsFactors = FALSE
) %>%
  transmute(
    immune_state_id = as.character(
      State_ID
    ),
    preliminary_family = Family,
    preliminary_state_name = State_Name,
    preliminary_confidence = Confidence,
    main_markers = Main_Markers
  )

mapping_qc <- read.csv(
  paste0(
    "results/tables/",
    "GSE144744_state_mapping_QC_by_state.csv"
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(
    immune_state_id = as.character(
      immune_state_id
    )
  )

# ------------------------------------------------------------
# Sample-by-state composition, including explicit zeroes
# ------------------------------------------------------------

sample_totals <- meta_full %>%
  count(
    sample,
    donor,
    group,
    cohort,
    name = "sample_total_tcells"
  )

sample_state <- meta_full %>%
  count(
    sample,
    donor,
    group,
    cohort,
    immune_state_id,
    name = "state_cells"
  ) %>%
  complete(
    nesting(
      sample,
      donor,
      group,
      cohort
    ),
    immune_state_id = factor(
      state_levels,
      levels = state_levels
    ),
    fill = list(
      state_cells = 0L
    )
  ) %>%
  left_join(
    sample_totals,
    by = c(
      "sample",
      "donor",
      "group",
      "cohort"
    )
  ) %>%
  mutate(
    sample_percent =
      100 * state_cells /
      sample_total_tcells,
    
    # Operational recurrence criterion:
    # at least 10 cells and at least 0.1% of T cells
    # within a biological sample.
    detected_010 =
      state_cells >= 10 &
      sample_percent >= 0.1
  )

write.csv(
  sample_state,
  paste0(
    "results/tables/state_audit/",
    "GSE144744_sample_by_state_portability_metrics.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Overall abundance and sample recurrence
# ------------------------------------------------------------

overall_summary <- meta_full %>%
  count(
    immune_state_id,
    name = "total_cells"
  ) %>%
  mutate(
    overall_percent =
      100 * total_cells /
      sum(total_cells)
  )

sample_recurrence <- sample_state %>%
  group_by(
    immune_state_id
  ) %>%
  summarise(
    total_samples = n(),
    samples_with_any_cells = sum(
      state_cells > 0
    ),
    samples_detected_010 = sum(
      detected_010
    ),
    sample_prevalence_any_percent =
      100 * mean(
        state_cells > 0
      ),
    sample_prevalence_010_percent =
      100 * mean(
        detected_010
      ),
    median_sample_percent =
      median(
        sample_percent
      ),
    mean_sample_percent =
      mean(
        sample_percent
      ),
    maximum_sample_percent =
      max(
        sample_percent
      ),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Donor recurrence and concentration
# ------------------------------------------------------------

donor_state <- meta_full %>%
  count(
    donor,
    immune_state_id,
    name = "state_cells"
  ) %>%
  complete(
    donor,
    immune_state_id = factor(
      state_levels,
      levels = state_levels
    ),
    fill = list(
      state_cells = 0L
    )
  )

donor_recurrence <- donor_state %>%
  group_by(
    immune_state_id
  ) %>%
  summarise(
    total_donors = n(),
    donors_with_any_cells = sum(
      state_cells > 0
    ),
    donor_prevalence_any_percent =
      100 * mean(
        state_cells > 0
      ),
    maximum_donor_contribution_percent =
      100 * max(state_cells) /
      sum(state_cells),
    .groups = "drop"
  )

top_donor <- donor_state %>%
  group_by(
    immune_state_id
  ) %>%
  slice_max(
    order_by = state_cells,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  transmute(
    immune_state_id,
    top_contributing_donor = donor,
    top_donor_state_cells = state_cells
  )

# ------------------------------------------------------------
# Cohort representation based on mean sample proportions
# ------------------------------------------------------------

cohort_state <- sample_state %>%
  group_by(
    immune_state_id,
    cohort
  ) %>%
  summarise(
    cohort_mean_sample_percent =
      mean(
        sample_percent
      ),
    cohort_median_sample_percent =
      median(
        sample_percent
      ),
    cohort_samples_detected_010 =
      sum(
        detected_010
      ),
    cohort_sample_prevalence_010_percent =
      100 * mean(
        detected_010
      ),
    .groups = "drop"
  )

write.csv(
  cohort_state,
  paste0(
    "results/tables/state_audit/",
    "GSE144744_state_by_cohort_portability_metrics.csv"
  ),
  row.names = FALSE
)

cohort_summary <- cohort_state %>%
  group_by(
    immune_state_id
  ) %>%
  summarise(
    cohorts_with_any_detected_samples = sum(
      cohort_samples_detected_010 > 0
    ),
    cohorts_with_prevalence_ge_20_percent = sum(
      cohort_sample_prevalence_010_percent >= 20
    ),
    maximum_cohort_mean_percent = max(
      cohort_mean_sample_percent
    ),
    minimum_cohort_mean_percent = min(
      cohort_mean_sample_percent
    ),
    max_to_min_cohort_ratio = ifelse(
      minimum_cohort_mean_percent > 0,
      maximum_cohort_mean_percent /
        minimum_cohort_mean_percent,
      NA_real_
    ),
    dominant_cohort_share_percent =
      100 * maximum_cohort_mean_percent /
      sum(cohort_mean_sample_percent),
    .groups = "drop"
  )

dominant_cohort <- cohort_state %>%
  group_by(
    immune_state_id
  ) %>%
  slice_max(
    order_by = cohort_mean_sample_percent,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  transmute(
    immune_state_id,
    dominant_cohort = cohort
  )

# ------------------------------------------------------------
# Integrated audit table
# ------------------------------------------------------------

state_audit <- tibble(
  immune_state_id = state_levels
) %>%
  left_join(
    atlas_annotations,
    by = "immune_state_id"
  ) %>%
  left_join(
    overall_summary %>%
      mutate(
        immune_state_id = as.character(
          immune_state_id
        )
      ),
    by = "immune_state_id"
  ) %>%
  left_join(
    sample_recurrence %>%
      mutate(
        immune_state_id = as.character(
          immune_state_id
        )
      ),
    by = "immune_state_id"
  ) %>%
  left_join(
    donor_recurrence %>%
      mutate(
        immune_state_id = as.character(
          immune_state_id
        )
      ),
    by = "immune_state_id"
  ) %>%
  left_join(
    top_donor %>%
      mutate(
        immune_state_id = as.character(
          immune_state_id
        )
      ),
    by = "immune_state_id"
  ) %>%
  left_join(
    cohort_summary %>%
      mutate(
        immune_state_id = as.character(
          immune_state_id
        )
      ),
    by = "immune_state_id"
  ) %>%
  left_join(
    dominant_cohort %>%
      mutate(
        immune_state_id = as.character(
          immune_state_id
        )
      ),
    by = "immune_state_id"
  ) %>%
  left_join(
    mapping_qc,
    by = "immune_state_id"
  ) %>%
  mutate(
    donor_concentration_flag =
      maximum_donor_contribution_percent >= 50,
    
    cohort_restriction_flag =
      cohorts_with_any_detected_samples == 1 |
      dominant_cohort_share_percent >= 80,
    
    low_recurrence_flag =
      sample_prevalence_010_percent < 20,
    
    low_transfer_confidence_flag =
      median_score < 0.75,
    
    provisional_portability_class =
      case_when(
        donor_concentration_flag ~
          "Donor-concentrated candidate",
        
        cohort_restriction_flag ~
          "Cohort-restricted candidate",
        
        low_recurrence_flag ~
          "Low-recurrence candidate",
        
        sample_prevalence_010_percent >= 50 &
          cohorts_with_prevalence_ge_20_percent == 3 ~
          "Broadly recurrent candidate",
        
        TRUE ~
          "Context-enriched recurrent candidate"
      )
  ) %>%
  arrange(
    factor(
      provisional_portability_class,
      levels = c(
        "Broadly recurrent candidate",
        "Context-enriched recurrent candidate",
        "Cohort-restricted candidate",
        "Donor-concentrated candidate",
        "Low-recurrence candidate"
      )
    ),
    desc(
      sample_prevalence_010_percent
    )
  )

write.csv(
  state_audit,
  paste0(
    "results/tables/state_audit/",
    "GSE144744_immune_state_portability_audit.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Compact summaries for review
# ------------------------------------------------------------

portability_class_summary <- state_audit %>%
  count(
    provisional_portability_class,
    name = "n_states"
  ) %>%
  arrange(
    desc(n_states)
  )

audit_review <- state_audit %>%
  select(
    immune_state_id,
    preliminary_state_name,
    total_cells,
    overall_percent,
    sample_prevalence_010_percent,
    donors_with_any_cells,
    maximum_donor_contribution_percent,
    cohorts_with_prevalence_ge_20_percent,
    dominant_cohort,
    dominant_cohort_share_percent,
    median_score,
    provisional_portability_class
  )

write.csv(
  portability_class_summary,
  paste0(
    "results/tables/state_audit/",
    "GSE144744_portability_class_summary.csv"
  ),
  row.names = FALSE
)

print(
  tibble::as_tibble(
    portability_class_summary
  ),
  n = Inf,
  width = Inf
)

print(
  tibble::as_tibble(
    audit_review
  ),
  n = Inf,
  width = Inf
)
