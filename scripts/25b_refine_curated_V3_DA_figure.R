# ============================================================
# 25b_refine_curated_V3_DA_figure.R
# Publication-oriented refinement of the curated V3
# differential-abundance figure. No model is refitted.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

INPUT_DIR <- paste0(
  "results/tables/",
  "curated_V3_differential_abundance"
)

FIGURE_DIR <- "results/figures"

ALL_RESULTS_FILE <- file.path(
  INPUT_DIR,
  "GSE144744_curated_V3_DA_all_analyses.csv"
)

PRIMARY_RESULTS_FILE <- file.path(
  INPUT_DIR,
  "GSE144744_curated_V3_DA_primary_results.csv"
)

required_files <- c(
  ALL_RESULTS_FILE,
  PRIMARY_RESULTS_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0L) {
  stop(
    "Required differential-abundance result file(s) not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

dir.create(
  FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

all_results <- read.csv(
  ALL_RESULTS_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

primary_results <- read.csv(
  PRIMARY_RESULTS_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "analysis_set",
  "contrast",
  "immune_state_id",
  "difference_percentage_points",
  "P.Value",
  "FDR",
  "n_pairs"
)

missing_columns <- setdiff(
  required_columns,
  colnames(all_results)
)

if (length(missing_columns) > 0L) {
  stop(
    "The all-analysis result table is missing column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

state_levels <- sprintf("IVS-%02d", 1:11)

contrast_levels <- c(
  "MS1_vs_HI1",
  "MS1nat_vs_MS1",
  "MS2_vs_HI2",
  "PPMS_vs_HI3"
)

contrast_labels <- c(
  MS1_vs_HI1 = paste0(
    "Untreated RRMS vs matched HI\n",
    "9 matched pairs"
  ),
  MS1nat_vs_MS1 = paste0(
    "Natalizumab vs pretreatment\n",
    "9 paired donors"
  ),
  MS2_vs_HI2 = paste0(
    "RRMS vs matched HI\n",
    "2 pairs; descriptive"
  ),
  PPMS_vs_HI3 = paste0(
    "PPMS vs matched HI\n",
    "10 matched pairs"
  )
)

analysis_levels <- c(
  "reference_only_logit",
  "confidence_050_logit",
  "full_asin"
)

analysis_labels <- c(
  reference_only_logit = "Reference cells only",
  confidence_050_logit =
    "Transferred score >= 0.50",
  full_asin = "Arcsin transformation"
)

# ------------------------------------------------------------
# Panel A: primary effects
# ------------------------------------------------------------

plot_primary <- primary_results %>%
  mutate(
    immune_state_id = factor(
      immune_state_id,
      levels = rev(state_levels)
    ),
    contrast = factor(
      contrast,
      levels = contrast_levels,
      labels = unname(
        contrast_labels[contrast_levels]
      )
    ),
    is_descriptive_rrms = (
      as.character(contrast) ==
        unname(
          contrast_labels["MS2_vs_HI2"]
        )
    ),
    significance_symbol = case_when(
      is_descriptive_rrms ~ "",
      FDR < 0.05 ~ "**",
      P.Value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    effect_label = paste0(
      sprintf(
        "%+.2f",
        difference_percentage_points
      ),
      significance_symbol
    )
  )

effect_limit <- max(
  abs(
    plot_primary$difference_percentage_points
  ),
  na.rm = TRUE
)

if (!is.finite(effect_limit) || effect_limit == 0) {
  effect_limit <- 1
}

panel_A <- ggplot(
  plot_primary,
  aes(
    x = contrast,
    y = immune_state_id,
    fill = difference_percentage_points
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = effect_label),
    size = 3.6,
    color = "#17202A"
  ) +
  scale_fill_gradient2(
    low = "#3C78B5",
    mid = "#F7F7F7",
    high = "#D64C4C",
    midpoint = 0,
    limits = c(-effect_limit, effect_limit),
    name = paste0(
      "Difference\n",
      "(percentage points)"
    )
  ) +
  labs(
    title = "Curated immune-state differential abundance",
    subtitle = paste0(
      "Positive values indicate enrichment in the disease or ",
      "treatment group. * nominal P < 0.05; ** FDR < 0.05. ",
      "Significance symbols are not shown for the two-pair ",
      "RRMS_HI comparison."
    ),
    x = NULL,
    y = "Curated immune state"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      hjust = 0.5,
      lineheight = 0.95
    ),
    axis.text.y = element_text(face = "bold"),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    plot.subtitle = element_text(
      color = "#5B5B5B",
      size = 10.2
    )
  )

# ------------------------------------------------------------
# Panel B: inferential-pattern sensitivity
# ------------------------------------------------------------

primary_for_sensitivity <- all_results %>%
  filter(analysis_set == "full_logit") %>%
  transmute(
    contrast,
    immune_state_id,
    primary_effect =
      difference_percentage_points,
    primary_P_Value = P.Value,
    primary_FDR = FDR,
    primary_signed_evidence =
      sign(difference_percentage_points) *
      -log10(pmax(P.Value, 1e-300)),
    primary_nominal = P.Value < 0.05
  )

sensitivity_state_comparison <- all_results %>%
  filter(analysis_set != "full_logit") %>%
  transmute(
    analysis_set,
    contrast,
    immune_state_id,
    sensitivity_effect =
      difference_percentage_points,
    sensitivity_P_Value = P.Value,
    sensitivity_FDR = FDR,
    sensitivity_signed_evidence =
      sign(difference_percentage_points) *
      -log10(pmax(P.Value, 1e-300)),
    sensitivity_nominal = P.Value < 0.05
  ) %>%
  left_join(
    primary_for_sensitivity,
    by = c(
      "contrast",
      "immune_state_id"
    )
  ) %>%
  mutate(
    direction_agreement = case_when(
      primary_effect == 0 &
        sensitivity_effect == 0 ~ TRUE,
      primary_effect == 0 |
        sensitivity_effect == 0 ~ FALSE,
      TRUE ~ sign(primary_effect) ==
        sign(sensitivity_effect)
    ),
    nominal_agreement = (
      primary_nominal ==
        sensitivity_nominal
    )
  )

safe_correlation <- function(x, y, method) {
  complete <- complete.cases(x, y)

  if (
    sum(complete) < 3L ||
    stats::sd(x[complete]) == 0 ||
    stats::sd(y[complete]) == 0
  ) {
    return(NA_real_)
  }

  stats::cor(
    x[complete],
    y[complete],
    method = method
  )
}

sensitivity_summary_refined <-
  sensitivity_state_comparison %>%
  group_by(analysis_set, contrast) %>%
  summarise(
    signed_evidence_spearman =
      safe_correlation(
        primary_signed_evidence,
        sensitivity_signed_evidence,
        method = "spearman"
      ),
    effect_pearson = safe_correlation(
      primary_effect,
      sensitivity_effect,
      method = "pearson"
    ),
    direction_agreement_percent =
      100 * mean(
        direction_agreement,
        na.rm = TRUE
      ),
    nominal_agreement_percent =
      100 * mean(
        nominal_agreement,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

write.csv(
  sensitivity_state_comparison,
  file.path(
    INPUT_DIR,
    paste0(
      "GSE144744_curated_V3_DA_",
      "refined_sensitivity_state_comparison.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  sensitivity_summary_refined,
  file.path(
    INPUT_DIR,
    paste0(
      "GSE144744_curated_V3_DA_",
      "refined_sensitivity_summary.csv"
    )
  ),
  row.names = FALSE
)

plot_sensitivity <- sensitivity_summary_refined %>%
  mutate(
    analysis_set = factor(
      analysis_set,
      levels = analysis_levels,
      labels = unname(
        analysis_labels[analysis_levels]
      )
    ),
    contrast = factor(
      contrast,
      levels = contrast_levels,
      labels = unname(
        contrast_labels[contrast_levels]
      )
    ),
    sensitivity_label = paste0(
      "\u03c1evidence = ",
      ifelse(
        is.na(signed_evidence_spearman),
        "NA",
        sprintf(
          "%.2f",
          signed_evidence_spearman
        )
      ),
      "\n",
      sprintf(
        "%.0f%% direction",
        direction_agreement_percent
      ),
      "; ",
      sprintf(
        "%.0f%% P-status",
        nominal_agreement_percent
      )
    )
  )

panel_B <- ggplot(
  plot_sensitivity,
  aes(
    x = contrast,
    y = analysis_set,
    fill = signed_evidence_spearman
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = sensitivity_label),
    size = 3.1,
    lineheight = 0.92,
    color = "#17202A"
  ) +
  scale_fill_gradient(
    low = "#F2F6FA",
    high = "#2E6FA8",
    limits = c(0, 1),
    na.value = "#E5E5E5",
    name = paste0(
      "Spearman \u03c1\n",
      "signed evidence"
    )
  ) +
  labs(
    title = "Sensitivity of inferential patterns",
    subtitle = paste0(
      "Agreement with the primary full-atlas logit analysis; ",
      "signed evidence combines effect direction and nominal P value"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      hjust = 0.5,
      lineheight = 0.95
    ),
    axis.text.y = element_text(face = "bold"),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    plot.subtitle = element_text(
      color = "#5B5B5B",
      size = 10.2
    )
  )

combined_figure <- (
  panel_A /
    panel_B
) +
  plot_layout(
    heights = c(1.85, 1)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 17
      )
    )
  )

figure_pdf <- file.path(
  FIGURE_DIR,
  paste0(
    "GSE144744_Tcells_V3_curated_state_",
    "differential_abundance_publication.pdf"
  )
)

figure_png <- file.path(
  FIGURE_DIR,
  paste0(
    "GSE144744_Tcells_V3_curated_state_",
    "differential_abundance_publication.png"
  )
)

ggsave(
  filename = figure_pdf,
  plot = combined_figure,
  width = 13.8,
  height = 11.2,
  units = "in",
  device = "pdf",
  useDingbats = FALSE
)

ggsave(
  filename = figure_png,
  plot = combined_figure,
  width = 13.8,
  height = 11.2,
  units = "in",
  dpi = 320,
  bg = "white"
)

message("\nRefined sensitivity summary:")
print(as.data.frame(sensitivity_summary_refined))

message("\nSaved refined PDF: ", figure_pdf)
message("Saved refined PNG: ", figure_png)
message(
  "Publication-oriented differential-abundance ",
  "figure completed."
)
