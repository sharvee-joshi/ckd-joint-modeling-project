#Load all required packages for the CKD project

packages <- c(
  "dplyr",
  "ggplot2",
  "survival",
  "survminer",
  "kableExtra",
  "nlme",
  "JMbayes2"
)

install_if_missing <- function(pkg){
  if(!pkg %in% installed.packages()){
    install.packages(pkg)
  }
}

lapply(packages, install_if_missing)
lapply(packages, library, character.only = TRUE)
