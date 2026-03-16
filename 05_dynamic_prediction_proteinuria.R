############################################################
# Dynamic prediction for graft survival using proteinuria
############################################################

## Packages
library(dplyr)
library(survival)
library(splines)
library(GLMMadaptive)
library(JMbayes2)
library(ggplot2)

#############
## Load data
#############
load("ckd.rdata")

prot <- ckd_data$prot
surv <- ckd_data$surv

#####################
## Basic cleaning
####################

surv <- as.data.frame(lapply(surv, as.vector))
surv$gender <- factor(surv$gender)

prot <- prot %>%
  arrange(id, years) %>%
  mutate(proteinuria = ifelse(is.na(proteinuria), 
                              NA, as.integer(as.vector(proteinuria))))

##########################
## Fit or load final joint model
########################

use_saved_model <- FALSE

if (use_saved_model && file.exists("jm_prot.rds")) {
  jm_prot <- readRDS("jm_prot.rds")
} else {
  
  ## Survival submodel
  cox_base <- coxph(
    Surv(fuyears, failure) ~ age + weight + gender,
    data = surv,
    x = TRUE
  )
  
  ## Longitudinal submodel: proteinuria with linear time
  glmm_prot <- mixed_model(
    fixed = proteinuria ~ years + age + weight + gender,
    random = ~ years | id,
    family = binomial(),
    data = prot
  )
  
  ## Joint model
  jm_prot <- jm(
    cox_base,
    glmm_prot,
    time_var = "years"
  )
  saveRDS(jm_prot, "jm_prot.rds")
}
summary(jm_prot)

#######################
## Helper: create landmark histories
########################

make_landmark_data <- function(df, landmark_time) {
  df %>% filter(years <= landmark_time, !is.na(proteinuria)) %>%
     mutate(fuyears = landmark_time,failure = 0)
}

##################
## Patient 1: 5564
##################

prot_counts <- prot %>%
  filter(!is.na(proteinuria)) %>%
  count(id, sort = TRUE)

head(prot_counts, 10)

one_id <- as.character(prot_counts$id[1])

ND_prot <- prot %>%
  filter(as.character(id) == one_id) %>%
  arrange(years)

ND_surv <- surv %>%
  filter(as.character(id) == one_id)

#####################################
## Landmark histories for patient 5564
######################################

ND_2  <- make_landmark_data(ND_prot, 2)
ND_5  <- make_landmark_data(ND_prot, 5)
ND_8  <- make_landmark_data(ND_prot, 8)
ND_10 <- make_landmark_data(ND_prot, 10)

nrow(ND_2)
nrow(ND_5)
nrow(ND_8)
nrow(ND_10)

####################################
## Prediction times for patient 5564
#####################################

times_2  <- seq(2, 15, by = 0.5)
times_5  <- seq(5, 15, by = 0.5)
times_8  <- seq(8, 15, by = 0.5)
times_10 <- seq(10, 15, by = 0.5)

#######################################
## Dynamic predictions for patient 5564
#######################################

pred_2 <- predict(jm_prot, newdata = ND_2,
                  process = "event", times = times_2)

pred_5 <- predict(jm_prot,newdata = ND_5,
                  process = "event", times = times_5)

pred_8 <- predict(jm_prot, newdata = ND_8,
                  process = "event", times = times_8)

pred_10 <- predict(jm_prot, newdata = ND_10,
                   process = "event", times = times_10)

#####################################
## Convert predictions to plotting data: patient 5564
######################################

pred_df_2 <- data.frame(time = pred_2$times, surv_prob = 1 - pred_2$pred, 
                        low = 1 - pred_2$upp, upp = 1 - pred_2$low, 
                        landmark = "After year 2", landmark_time = 2)

pred_df_5 <- data.frame(time = pred_5$times, surv_prob = 1 - pred_5$pred, 
                        low = 1 - pred_5$upp, upp = 1 - pred_5$low, 
                        landmark = "After year 5", landmark_time = 5)

pred_df_8 <- data.frame(time = pred_8$times, surv_prob = 1 - pred_8$pred,
                        low = 1 - pred_8$upp, upp = 1 - pred_8$low,
                        landmark = "After year 8", landmark_time = 8)

pred_df_10 <- data.frame(time = pred_10$times, surv_prob = 1 - pred_10$pred,
                         low = 1 - pred_10$upp, upp = 1 - pred_10$low,
                         landmark = "After year 10", landmark_time = 10)

pred_all <- bind_rows(pred_df_2, pred_df_5, pred_df_8, pred_df_10)

pred_all$landmark <- factor(pred_all$landmark,
                            levels = c("After year 2", "After year 5", 
                                       "After year 8", "After year 10"))

landmark_lines <- data.frame(
  landmark = factor(
    c("After year 2", "After year 5", "After year 8", "After year 10"),
    levels = c("After year 2", "After year 5", "After year 8", "After year 10")
  ),
  landmark_time = c(2, 5, 8, 10)
)

#########################
##plot for patient 5564
##########################

p1 <- ggplot(pred_all, aes(x = time, y = surv_prob)) +
  geom_ribbon(aes(ymin = low, ymax = upp),
              fill = "#C96A6A", alpha = 0.18) +
  geom_line(color = "#8B1E1E", linewidth = 0.8) +
  geom_line(aes(y = low),linetype = "dashed", color = "#C96A6A",linewidth = 0.65) +
  geom_line(aes(y = upp), linetype = "dashed", color = "#C96A6A",linewidth = 0.65) +
  geom_rug(
    data = ND_prot %>% filter(!is.na(proteinuria)), aes(x = years), inherit.aes = FALSE,
    sides = "b", alpha = 0.15) +
  geom_vline(
    data = landmark_lines, aes(xintercept = landmark_time), linetype = "dashed",
    color = "grey35", linewidth = 0.8) +
  facet_wrap(~ landmark, ncol = 2) +
  labs(
    x = "Years since transplant",
    y = "Predicted graft survival probability"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 17, hjust = 0.5, face = "plain"),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.4),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 12)
  )

ggsave(
  filename = "patient5564_dynamic_prediction.png",
  plot = p1, width = 10, height = 7, dpi = 300)

########################
## Patient 2: 5500
#######################

second_id <- "5500"

ND2_prot <- prot %>%
  filter(as.character(id) == second_id) %>%
  arrange(years)

ND2_surv <- surv %>%
  filter(as.character(id) == second_id)

## Optional checks
ND2_prot
ND2_surv
range(ND2_prot$years, na.rm = TRUE)
table(ND2_prot$proteinuria, useNA = "ifany")

#####################################
## Landmark histories for patient 5500
#####################################

ND2_1 <- make_landmark_data(ND2_prot, 1)
ND2_2 <- make_landmark_data(ND2_prot, 2)
ND2_4 <- make_landmark_data(ND2_prot, 4)
ND2_5 <- make_landmark_data(ND2_prot, 5)

nrow(ND2_1)
nrow(ND2_2)
nrow(ND2_4)
nrow(ND2_5)

##################################
##times t for patient 5500
#################################

times2_1 <- seq(1, 5.9, by = 0.25)
times2_2 <- seq(2, 5.9, by = 0.25)
times2_4 <- seq(4, 5.9, by = 0.25)
times2_5 <- seq(5, 5.9, by = 0.25)

######################################
##Dynamic predictions for patient 5500
######################################

pred2_1 <- predict(jm_prot, newdata = ND2_1, 
                   process = "event", times = times2_1)

pred2_2 <- predict(jm_prot, newdata = ND2_2, 
                   process = "event", times = times2_2)

pred2_4 <- predict(jm_prot, newdata = ND2_4, 
                   process = "event", times = times2_4)

pred2_5 <- predict(jm_prot, newdata = ND2_5, 
                   process = "event", times = times2_5)

############################################################
## Convert predictions to plotting data: patient 5500
############################################################

pred2_df_1 <- data.frame(time = pred2_1$times, surv_prob = 1 - pred2_1$pred, 
                         low = 1 - pred2_1$upp, upp = 1 - pred2_1$low, 
                         landmark = "After year 1")

pred2_df_2 <- data.frame(time = pred2_2$times, surv_prob = 1 - pred2_2$pred, 
                         low = 1 - pred2_2$upp, upp = 1 - pred2_2$low, 
                         landmark = "After year 2")

pred2_df_4 <- data.frame(time = pred2_4$times, surv_prob = 1 - pred2_4$pred,
                         low = 1 - pred2_4$upp,upp = 1 - pred2_4$low,
                         landmark = "After year 4")

pred2_df_5 <- data.frame(time = pred2_5$times, surv_prob = 1 - pred2_5$pred, 
                         low = 1 - pred2_5$upp, upp = 1 - pred2_5$low, 
                         landmark = "After year 5")

pred_all_2 <- bind_rows(pred2_df_1, pred2_df_2, pred2_df_4, pred2_df_5)

pred_all_2$landmark <- factor(
  pred_all_2$landmark,
  levels = c("After year 1", "After year 2", "After year 4", "After year 5")
)

landmark_lines_2 <- data.frame(
  landmark = factor(
    c("After year 1", "After year 2", "After year 4", "After year 5"),
    levels = c("After year 1", "After year 2", "After year 4", "After year 5")
  ),
  landmark_time = c(1, 2, 4, 5)
)

##########################
##plot for patient 5500
#########################

p2 <- ggplot(pred_all_2, aes(x = time, y = surv_prob)) +
  geom_ribbon(
    aes(ymin = low, ymax = upp),
    fill = "#C96A6A", alpha = 0.18
  ) +
  geom_line(color = "#8B1E1E", linewidth = 0.8) +
  geom_line(
    aes(y = low),
    linetype = "dashed",
    color = "#C96A6A",
    linewidth = 0.65
  ) +
  geom_line(
    aes(y = upp),
    linetype = "dashed",
    color = "#C96A6A",
    linewidth = 0.65
  ) +
  geom_vline(
    data = landmark_lines_2,
    aes(xintercept = landmark_time),
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.8
  ) +
  
  facet_wrap(~ landmark, ncol = 2) +
  labs(
    x = "Years since transplant",
    y = "Predicted graft survival probability"
  ) +
  geom_rug(
    data = ND_prot %>% filter(!is.na(proteinuria)),
    aes(x = years),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.15
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 17, hjust = 0.5, face = "plain"),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.4),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 12)
  )

ggsave(
  filename = "patient5500_dynamic_prediction.png",
  plot = p2, width = 10, height = 7, dpi = 300)

