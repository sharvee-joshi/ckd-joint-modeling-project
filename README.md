# Chronic Kidney Disease: Joint and Dynamic Predictive Modeling on Risk of Graft Failure On Renal Transplant Patients
**Methods:** Survival Analysis, Joint Modeling, Dynamic Prediction  
**Language:** R

## Authors: Makena Grigsby & Sharvee Joshi

## Project Overview

This repository contains code for the analysis of a chronic kidney disease (CKD) dataset focusing on the progression toward kidney graft failure. The analysis explores longitudinal biomarkers collected during follow-up and their relationship with time-to-event outcomes.

Initial exploratory analyses include summaries of missingness patterns, Kaplan–Meier survival analysis, and visualization of biomarker trajectories over time. These exploratory steps help characterize the structure of the data and provide insight into how biomarkers evolve prior to graft failure.

Building on this exploratory work, the analysis implements joint models that link longitudinal biomarker measurements with the survival process for graft failure. This joint modeling framework allows the longitudinal trajectory of biomarkers to inform the risk of graft failure while properly accounting for repeated measurements within individuals.

Finally, dynamic prediction methods are applied using the fitted joint models to estimate patient-specific survival probabilities over time. These predictions update as new biomarker measurements become available, allowing individualized risk assessment throughout the follow-up period.

This repository contains reproducible code, the final R Markdown report, and all supporting functions.

## Data Description

The dataset used in this analysis contains longitudinal measurements of kidney function biomarkers collected from patients following kidney transplantation. Repeated measurements are available for each patient over follow-up time.

Key variables include:

- **id** – unique patient identifier  
- **years** – follow-up time (in years) since transplantation  
- **gfr** – estimated glomerular filtration rate (GFR), a measure of kidney function  
- **haematocrit** – percentage of red blood cells in blood, reflecting overall blood health  
- **proteinuria** – indicator of protein in urine, a clinical marker of kidney damage  
- **failure** – indicator variable for graft failure (1 = graft failure occurred, 0 = censored)  
- **fuyears** – total follow-up time for each patient  
- **age, weight, gender** – baseline patient characteristics used as covariates

The longitudinal biomarkers (GFR, haematocrit, and proteinuria) are observed repeatedly during follow-up, while graft failure is treated as the event outcome in the survival analysis.

## Files Included
```
00_libraries.R
01_data_setup.R
02_missingness_summary.R
03_eda_plots.R

#add more files here

ckd.rdata
```
---

## Required R Packages

The script `00_libraries.R` automatically installs and loads the required packages.  
These include:

- dplyr
- ggplot2
- survival
- survminer
- kableExtra
- nlme
- JMbayes2

No manual installation is required.

---

## How to Run the Code
To reproduce the analysis, please follow the steps outlined below, especially if this is your first time using GitHub.

### 1. Download the repository
Download the repository as a ZIP file from GitHub and unzip it on your computer. 
Make sure the folder contains the files presented above:

### 2. Open RStudio
Open **RStudio** and set the working directory to the folder cotaining your downloaded files. Be sure to unzip the file before completing this step. You can do this in R using:
```
setwd("path/to/the/unzipped/folder")
```

### 3. Run the scripts in order
Run the following commands in the R console:
```
source("01_data_setup.R")
source("02_missingness_summary.R")
source("03_eda_plots.R")

#add more files here
```

## Output

Running the scripts will generate:

- A missingness summary table
- Kaplan–Meier survival plot
- Smoothed biomarker trajectory plots
- Joint Modeling Results for all three biomarkers
- Dynamic Predictive Results and Grahpics


## References

Rizopoulos, D. (2012). *Joint Models for Longitudinal and Time-to-Event Data: With Applications in R.* Chapman & Hall/CRC.

Cox, D. R. (1972). Regression models and life-tables. *Journal of the Royal Statistical Society: Series B*, 34(2), 187–220.

Kaplan, E. L., & Meier, P. (1958). Nonparametric estimation from incomplete observations. *Journal of the American Statistical Association*, 53(282), 457–481.

Therneau, T. (2024). *A Package for Survival Analysis in R*. R package **survival**.

Kassambara, A., & Kosinski, M. (2024). *survminer: Drawing Survival Curves using 'ggplot2'*. R package.



## Team Contributions

| Team Member | Contributions |
|-------------|---------------|
| **Sharvee** | GitHub repository setup and organization, exploratory data analysis (EDA), missingness summaries, survival analysis visualizations, abstract writing, and drafting the joint modeling framework section of the report. |
| **Makena**  | Implementation and analysis of the dynamic prediction models, major contributions to report writing and proposal development. |
| **Sharvee & Makena** | Collaborative development of presentation slides and overall project interpretation. |

