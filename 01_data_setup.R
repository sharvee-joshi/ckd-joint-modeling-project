#01_datasetup.R
#Purpose: Load and prepare CKD data for analyses

# ------------------------------------------------------------------
# 0. Load required packages
# ------------------------------------------------------------------
source("00_libraries.R")

# ------------------------------------------------------------------
# 1. Load raw data
# ------------------------------------------------------------------
load("ckd.rdata")

#Inspect available objects
ls()
names(ckd_data)

#Extract datasets from ckd_data
prot <- ckd_data$prot
haem <- ckd_data$haem
gfr  <- ckd_data$gfr
surv <- ckd_data$surv

#Inspect structure
str(prot)
str(haem)
str(gfr)
str(surv)

# ------------------------------------------------------------------
# 2. Clean baseline survival data
# ------------------------------------------------------------------
surv <- as.data.frame(lapply(surv, as.vector))
surv$gender <- factor(surv$gender)

# ------------------------------------------------------------------
# 3. Arrange longitudinal datasets
# ------------------------------------------------------------------
gfr  <- gfr  %>% arrange(id, years)
haem <- haem %>% arrange(id, years)
prot <- prot %>% arrange(id, years)

# ------------------------------------------------------------------
# 4. Standardize variable names
# ------------------------------------------------------------------
gfr  <- gfr  %>% rename(gfr_value = gfr)

haem <- haem %>% rename(haematocrit = haematocrit)

prot <- prot %>% rename(proteinuria = proteinuria)

# ------------------------------------------------------------------
# 5. Merge longitudinal outcomes into one long dataset
# ------------------------------------------------------------------
long_data <- full_join(
  gfr, haem,
  by = c("id", "years", "failure", "fuyears", "weight", "age", "gender")
)

long_data <- full_join(
  long_data, prot,
  by = c("id", "years", "failure", "fuyears", "weight", "age", "gender")
)

# ------------------------------------------------------------------
# 6. Quick checks
# ------------------------------------------------------------------
#Number of repeated measurements per subject for GFR
table(table(gfr$id))

#Check merged data structure
str(long_data)
summary(long_data)

# ------------------------------------------------------------------
# 8. Save cleaned objects for later scripts
# ------------------------------------------------------------------
save(
  gfr, haem, prot, surv, long_data,
  file = "data/ckd_prepped.RData"
)
