#03_eda_plots.R
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
# 2. Kaplan-Meier curve by baseline proteinuria
# --------------------------------------------------
km_fit <- survfit(Surv(fuyears, failure) ~ proteinuria, data = surv_km)

p_km <- ggsurvplot(
  km_fit,
  data = surv_km,
  xlab = "Follow-Up Time (Years)",
  ylab = "Survival Probability",
  legend.title = "Proteinuria",
  legend.labs = c("Absent", "Present"),
  palette = c("#2C7BB6", "#D7191C"),
  legend = "right",
  ggtheme = theme_minimal()
)

# --------------------------------------------------
# 3. Smoothed trajectory line for GFR
# --------------------------------------------------
p_gfr <- ggplot(
  gfr %>% filter(!is.na(gfr_value)),
  aes(
    x = years,
    y = gfr_value,
    color = factor(
      failure,
      levels = c(0, 1),
      labels = c("No Graft Failure", "Graft Failure")
    )
  )
) +
  geom_smooth(se = FALSE, linewidth = 1.1) +
  scale_color_manual(values = c("#2C7BB6", "#D7191C")) +
  labs(
    x = "Follow-Up Time (Years)",
    y = "GFR",
    color = "Graft Failure Status"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

print(p_gfr)

# --------------------------------------------------
# 4. Smoothed trajectory line for haematocrit
# --------------------------------------------------
p_haem <- ggplot(
  haem %>% filter(!is.na(haematocrit)),
  aes(
    x = years,
    y = haematocrit,
    color = factor(
      failure,
      levels = c(0, 1),
      labels = c("No Graft Failure", "Graft Failure")
    )
  )
) +
  geom_smooth(se = FALSE, linewidth = 1.1) +
  scale_color_manual(values = c("#2C7BB6", "#D7191C")) +
  labs(
    x = "Follow-Up Time (Years)",
    y = "Haematocrit",
    color = "Graft Failure Status"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

print(p_haem)
