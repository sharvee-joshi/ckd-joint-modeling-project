library(dplyr)
library(survival)
library(GLMMadaptive)
library(JMbayes2)

# Load data
load("ckd.rdata")
prot <- ckd_data$prot
surv <- ckd_data$surv

# Basic cleaning
surv <- as.data.frame(lapply(surv, as.vector))
surv$gender <- factor(surv$gender)

prot <- prot %>%
  arrange(id, years) %>%
  mutate(
    proteinuria = ifelse(is.na(proteinuria), NA, as.integer(as.vector(proteinuria))),
    gender = factor(gender)
  )

# Baseline survival model
cox_base <- coxph(
  Surv(fuyears, failure) ~ age + weight + gender,
  data = surv,
  x = TRUE
)

# Working MAR model for missing proteinuria
imp_fit <- glm(
  proteinuria ~ years + age + weight + gender + failure + fuyears,
  family = binomial(),
  data = prot,
  subset = !is.na(proteinuria)
)

# Impute missing proteinuria under a delta-adjusted scenario
impute_proteinuria_delta <- function(prot_data, delta = 0, seed = 123) {
  set.seed(seed)
  
  dat <- prot_data
  miss_idx <- which(is.na(dat$proteinuria))
  
  if (length(miss_idx) == 0) return(dat)
  
  lp_mar <- predict(imp_fit, newdata = dat[miss_idx, ], type = "link")
  p_delta <- plogis(lp_mar + delta)
  
  dat$proteinuria[miss_idx] <- rbinom(length(miss_idx), 1, p_delta)
  dat
}

# Refit the joint model on the completed data
fit_jm_delta <- function(prot_completed) {
  glmm_prot_delta <- mixed_model(
    fixed = proteinuria ~ years + age + weight + gender,
    random = ~ years | id,
    family = binomial(),
    data = prot_completed
  )
  
  jm(cox_base, glmm_prot_delta, time_var = "years")
}

# Pull out the proteinuria association estimate
extract_assoc <- function(jm_fit, delta_label) {
  s <- summary(jm_fit)$Survival
  
  data.frame(
    delta = delta_label,
    estimate = s["value(proteinuria)", "Mean"],
    lower = s["value(proteinuria)", "2.5%"],
    upper = s["value(proteinuria)", "97.5%"],
    rhat = s["value(proteinuria)", "Rhat"]
  )
}

# Build a delta grid using +/- 1 SD on the logit scale
miss_idx <- which(is.na(prot$proteinuria))
lp_missing_mar <- predict(imp_fit, newdata = prot[miss_idx, ], type = "link")
delta_sd <- sd(lp_missing_mar, na.rm = TRUE)
delta_grid <- round(c(-delta_sd, 0, delta_sd), 3)

# Run the sensitivity analysis
assoc_results <- lapply(seq_along(delta_grid), function(k) {
  d <- delta_grid[k]
  cat("Running delta =", d, "\n")
  
  prot_delta <- impute_proteinuria_delta(prot, delta = d, seed = 100 + k)
  jm_fit_delta <- fit_jm_delta(prot_delta)
  extract_assoc(jm_fit_delta, d)
}) %>%
  bind_rows()

assoc_results
