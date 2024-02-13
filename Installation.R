#renv::init()

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")
BiocManager::install("TCGAbiolinks")
BiocManager::install("edgeR")
BiocManager::install("biomaRt")
install.packages("DT")
install.packages("glmnet")
install.packages("pROC")
install.packages("caTools")
install.packages("survivalROC")
install.packages("risksetROC")
install.packages("lattice")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
if (!requireNamespace("biomaRt", quietly = TRUE))  install.packages("biomaRt")

renv::snapshot()

# Github configuration --> escrever no terminal
#git config --global user.email "you@example.com"
#git config --global user.name "Your Name"