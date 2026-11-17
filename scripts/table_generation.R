# ------------------------------------------------------------
# Cell types vs group
# ------------------------------------------------------------

tab_group <- table(seu$basictype, seu$group)

write.csv(
  as.data.frame.matrix(tab_group),
  file = file.path(
    TAB_DIR,
    paste0(GEO_ACCESSION, "_basictype_vs_group.csv")
  )
)

tab_group

# ------------------------------------------------------------
# Cell types vs cohort
# ------------------------------------------------------------

tab_cohort <- table(seu$basictype, seu$cohort)

write.csv(
  as.data.frame.matrix(tab_cohort),
  file = file.path(
    TAB_DIR,
    paste0(GEO_ACCESSION, "_basictype_vs_cohort.csv")
  )
)

tab_cohort

# ------------------------------------------------------------
# Percentuali per group
# ------------------------------------------------------------

tab_group_pct <- prop.table(tab_group, margin = 2) * 100

write.csv(
  round(tab_group_pct, 2),
  file = file.path(
    TAB_DIR,
    paste0(GEO_ACCESSION, "_basictype_vs_group_percent.csv")
  )
)

tab_group_pct

# ------------------------------------------------------------
# Percentuali per cohort
# ------------------------------------------------------------

tab_cohort_pct <- prop.table(tab_cohort, margin = 2) * 100

write.csv(
  round(tab_cohort_pct, 2),
  file = file.path(
    TAB_DIR,
    paste0(GEO_ACCESSION, "_basictype_vs_cohort_percent.csv")
  )
)

tab_cohort_pct
