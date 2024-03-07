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
install.packages("SIS")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
if (!requireNamespace("biomaRt", quietly = TRUE))  install.packages("biomaRt")

#renv::snapshot()

# Github configuration --> escrever no terminal
# em tools global git - criar ssh a colocar esse no github online
#git config --global user.email "beatnevelei@gmail.com"
#git config --global user.name "Beatriz Leitao"
#git remote add Monet https://github.com/BeatrizNL/Monet.git
# Depois de fazer commit das alterações que queres, voltar ao terminal
#git push Monet