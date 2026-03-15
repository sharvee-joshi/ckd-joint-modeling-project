#02_missingness_summary.R
#Purpose: Summarize missingness for the three longitudinal biomarkers

# ------------------------------------------------------------------
# 0. Load required packages
# ------------------------------------------------------------------
source("00_libraries.R")

# ------------------------------------------------------------------
# 1. Load prepped data
# ------------------------------------------------------------------
load("ckd_prepped.RData")

# --------------------------------------------------
# 2. Create missingness summary table
# --------------------------------------------------
missing_summary <- data.frame(
  Variable = c("GFR", "Haematocrit", "Proteinuria"),
  Total_Obs = c(nrow(gfr), nrow(haem), nrow(prot)),
  Missing = c(
    sum(is.na(gfr$gfr_value)),
    sum(is.na(haem$haematocrit)),
    sum(is.na(prot$proteinuria))
  )
)

missing_summary$Percent_Missing <- round(
  100 * missing_summary$Missing / missing_summary$Total_Obs,
  1
)

# --------------------------------------------------
# 3. Print table
# --------------------------------------------------
missing_summary %>%
  kable(caption = "Missingness Summary for Longitudinal Biomarkers") %>%
  kable_styling(full_width = FALSE)

write.csv(missing_summary, "missing_summary.csv", row.names = FALSE)