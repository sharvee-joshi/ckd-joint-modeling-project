# Chronic Kidney Disease: Joint and Dynamic Predictive Modeling

## Authors: Makena Grigsby & Sharvee Joshi

## Project Overview

This repository contains code for the analysis of a chronic kidney disease (CKD) dataset focusing on the progression toward kidney graft failure. The analysis explores longitudinal biomarkers collected during follow-up and their relationship with time-to-event outcomes.

Initial exploratory analyses include summaries of missingness patterns, Kaplan–Meier survival analysis, and visualization of biomarker trajectories over time. These exploratory steps help characterize the structure of the data and provide insight into how biomarkers evolve prior to graft failure.

Building on this exploratory work, the analysis implements joint models that link longitudinal biomarker measurements with the survival process for graft failure. This joint modeling framework allows the longitudinal trajectory of biomarkers to inform the risk of graft failure while properly accounting for repeated measurements within individuals.

Finally, dynamic prediction methods are applied using the fitted joint models to estimate patient-specific survival probabilities over time. These predictions update as new biomarker measurements become available, allowing individualized risk assessment throughout the follow-up period.

This repository contains reproducible code, the final R Markdown report, and all supporting functions.

## Files Included
```
00_libraries.R
01_data_setup.R
02_missingness_summary.R
03_eda_plots.R

#add more files here

ckd.rdata
```


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

- 

## References



## Team Contributions 
| Team Member            | Contributions                                                                                                                                                                                                                                |
| ---------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Sharvee** | GitHub setup,  |
| **Makena**  | Stuff                                         |

