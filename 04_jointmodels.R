#04_jointmodels.R
#Purpose: Joint models for GFR, haematocrit, and proteinuria
#compare selected models, extract association parameters
#generate diagnostic plots for longitudinal and survival submodels.

# ------------------------------------------------------------------
# 0. Load required packages
# ------------------------------------------------------------------
source("00_libraries.R")

# ------------------------------------------------------------------
# 1. Load prepped data
# ------------------------------------------------------------------
load("ckd_prepped.RData")

# ------------------------------------------------------------------
# 2. Final data formatting for modeling
# ------------------------------------------------------------------
#make sure survival data is stored as a standard data frame
surv <- as.data.frame(lapply(surv, as.vector))
surv$gender <- factor(surv$gender)

#sort data by subject and years
gfr  <- gfr  %>% arrange(id, years)
haem <- haem %>% arrange(id, years)
prot <- prot %>% arrange(id, years)

#make sure all variables are model ready
gfr$gfr_value <- as.numeric(as.vector(gfr$gfr))
haem$haematocrit <- as.numeric(as.vector(haem$haematocrit))
prot$proteinuria <- ifelse(
  is.na(prot$proteinuria),
  NA,
  as.integer(as.vector(prot$proteinuria))
)

#turning cateogirical into factors for easier time later
gfr$gender  <- factor(gfr$gender)
haem$gender <- factor(haem$gender)
prot$gender <- factor(prot$gender)

# ------------------------------------------------------------------
# 3. Shared survival submodel
# ------------------------------------------------------------------
#using this cox model as shared for the other longitudinal submodels
seed <- 293
set.seed(seed) #reporducibility purpose

cox_base <- coxph(
  Surv(fuyears, failure) ~ age + weight + gender,
  data = surv,
  x = TRUE
)

print(summary(cox_base)) #feel free to comment this out for less output

# ------------------------------------------------------------------
# 4. Joint model for GFR
# ------------------------------------------------------------------
#continuous outcome - LMM
set.seed(seed)

#added a control line
#feel free to mess around with this, but for the purpose of this project
#i kept the iterations low
#but make sure all iterations match for each model! 
lme_gfr <- lme(
  fixed = gfr_value ~ years + age + weight + gender,
  random = ~ years | id,
  data = gfr,
  na.action = na.exclude,
  control = lmeControl(
    opt = "optim",
    msMaxIter = 200,
    maxIter = 200,
    niterEM = 50,
    msVerbose = TRUE
  )
)

set.seed(seed)
jm_gfr <- jm(
  cox_base,
  lme_gfr,
  time_var = "years"
)

jm_gfr_sum <- summary(jm_gfr)

print(jm_gfr_sum)

# ------------------------------------------------------------------
# 5. Joint model for haematocrit
# ------------------------------------------------------------------
#continuous outcome - LMM
set.seed(seed)

#added a control line
#feel free to mess around with this, but for the purpose of this project
#i kept the iterations low
#but make sure all iterations match for each model! 
lme_haem <- lme(
  fixed = haematocrit ~ years + age + weight + gender,
  random = ~ years | id,
  data = haem,
  na.action = na.exclude,
  control = lmeControl(
    opt = "optim",
    msMaxIter = 200,
    maxIter = 200,
    niterEM = 50,
    msVerbose = TRUE
  )
)

set.seed(seed)
jm_haem <- jm(
  cox_base,
  lme_haem,
  time_var = "years"
)

jm_haem_sum <- summary(jm_haem)

print(jm_haem_sum)

# ------------------------------------------------------------------
# 6. Joint model for proteinuria
# ------------------------------------------------------------------
# Binary longitudinal outcome modeled with a generalized mixed model
set.seed(seed)

#had to get rid of control line
#not sure what i could do here to make sure iterations match the other two
glmm_prot <- mixed_model(
  fixed = proteinuria ~ years + age + weight + gender,
  random = ~ years | id,
  family = binomial(),
  data = prot,
  na.action = na.exclude
)

set.seed(seed)
jm_prot <- jm(
  cox_base,
  glmm_prot,
  time_var = "years"
)

jm_prot_sum <- summary(jm_prot)

print(jm_prot_sum)

# ------------------------------------------------------------------
# 7. Extract fit statistics from each joint model
# ------------------------------------------------------------------
#we will use this for model comparison
get_fitstats <- function(jm_fit) {
  s <- summary(jm_fit)
  
  data.frame(
    marginal_WAIC = s$fit_stats$marginal$WAIC,
    conditional_WAIC = s$fit_stats$conditional$WAIC,
    marginal_DIC = s$fit_stats$marginal$DIC,
    conditional_DIC = s$fit_stats$conditional$DIC,
    marginal_LPML = s$fit_stats$marginal$LPML,
    conditional_LPML = s$fit_stats$conditional$LPML
  )
}

stats_gfr  <- get_fitstats(jm_gfr)
stats_haem <- get_fitstats(jm_haem)
stats_prot <- get_fitstats(jm_prot)

# ------------------------------------------------------------------
# 8. WAIC comparison table for GFR and haematocrit only
# ------------------------------------------------------------------
#typically, we want a lower WAIC 
#not sure if I'm allowed to do this but its fine
model_comp <- data.frame(
  Biomarker = c("GFR", "Haematocrit", "Proteinuria"),
  WAIC = c(
    stats_gfr$marginal_WAIC,
    stats_haem$marginal_WAIC,
    stats_prot$marginal_WAIC
  ),
  DIC = c(
    stats_gfr$marginal_DIC,
    stats_haem$marginal_DIC,
    stats_prot$marginal_DIC
  ),
  LPML = c(
    stats_gfr$marginal_LPML,
    stats_haem$marginal_LPML,
    stats_prot$marginal_LPML
  )
) %>%
  arrange(WAIC)

print(model_comp)


# ------------------------------------------------------------------
# 9. Diagnostic plots for the GFR longitudinal model
# ------------------------------------------------------------------
#needed: resids and qq plots

gfr_diag <- data.frame(
  fitted = fitted(lme_gfr),
  resid = resid(lme_gfr, type = "pearson")
)

p_gfr_resid <- ggplot(gfr_diag, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Fitted Values",
    y = "Pearson Residuals",
    title = "GFR Model: Residuals vs Fitted"
  ) +
  theme_minimal()

p_gfr_qq <- ggplot(gfr_diag, aes(sample = resid)) +
  stat_qq() +
  stat_qq_line() +
  labs(
    title = "GFR Model: Normal Q-Q Plot"
  ) +
  theme_minimal()

print(p_gfr_resid)
print(p_gfr_qq)

# ------------------------------------------------------------------
# 10. Diagnostic plots for the haematocrit longitudinal model
# ------------------------------------------------------------------
#needed: resids and qq plots
haem_diag <- data.frame(
  fitted = fitted(lme_haem),
  resid = resid(lme_haem, type = "pearson")
)

p_haem_resid <- ggplot(haem_diag, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Fitted Values",
    y = "Pearson Residuals",
    title = "Haematocrit Mixed Model: Residuals vs Fitted"
  ) +
  theme_minimal()

p_haem_qq <- ggplot(haem_diag, aes(sample = resid)) +
  stat_qq() +
  stat_qq_line() +
  labs(
    title = "Haematocrit Mixed Model: Normal Q-Q Plot"
  ) +
  theme_minimal()

print(p_haem_resid)
print(p_haem_qq)

# ------------------------------------------------------------------
# 11. Diagnostic plots for the proteinuria longitudinal model
# ------------------------------------------------------------------
#no qqplot is needed as this is a glm and we're using binary data
prot_diag <- data.frame(
  fitted = fitted(glmm_prot, type = "mean_subject"),
  resid  = residuals(glmm_prot)
) %>%
  filter(complete.cases(.))

p_prot_resid <- ggplot(prot_diag, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Fitted Values",
    y = "Residuals",
    title = "Proteinuria Mixed Model: Residuals vs Fitted"
  ) +
  theme_minimal()

p_prot_hist <- ggplot(prot_diag, aes(x = resid)) +
  geom_histogram(bins = 30) +
  labs(
    x = "Residuals",
    y = "Count",
    title = "Proteinuria Mixed Model: Residual Histogram"
  ) +
  theme_minimal()

print(p_prot_resid)
# ------------------------------------------------------------------
# 12. Martingale residual plots for the shared Cox survival model
# ------------------------------------------------------------------
cox_diag <- data.frame(
  martingale = residuals(cox_base, type = "martingale"),
  lp = predict(cox_base, type = "lp"),
  age = surv$age,
  weight = surv$weight
)

p_cox_lp <- ggplot(cox_diag, aes(x = lp, y = martingale)) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE, method = "loess") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Linear Predictor",
    y = "Martingale Residuals",
    title = "Cox Model: Martingale Residuals vs Linear Predictor"
  ) +
  theme_minimal()

p_cox_age <- ggplot(cox_diag, aes(x = age, y = martingale)) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE, method = "loess") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Age",
    y = "Martingale Residuals",
    title = "Cox Model: Martingale Residuals vs Age"
  ) +
  theme_minimal()

p_cox_weight <- ggplot(cox_diag, aes(x = weight, y = martingale)) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE, method = "loess") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Weight",
    y = "Martingale Residuals",
    title = "Cox Model: Martingale Residuals vs Weight"
  ) +
  theme_minimal()

print(p_cox_lp)
print(p_cox_age)
print(p_cox_weight)
# ------------------------------------------------------------------
# 13. MCMC convergence diagnostics: trace plots
# ------------------------------------------------------------------
#we need to make sure we're also converging bc we've used Bayes techniques 
#alpha and betas presented, similar to notes given to us in class

#GFR joint model
traceplot(jm_gfr, parm = "alphas")
traceplot(jm_gfr, parm = "betas")

#Haematocrit joint model
traceplot(jm_haem, parm = "alphas")
traceplot(jm_haem, parm = "betas")

#Proteinuria joint model
traceplot(jm_prot, parm = "alphas")
traceplot(jm_prot, parm = "betas")
