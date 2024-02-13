if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")
BiocManager::install("TCGAbiolinks")
install.packages("DT")

install.packages("survminer") # survival

# Required libraries:
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)

library(survival)           # survival
library(survminer)          # survival

#----------------------------recolher os dados

#recolher os dados de GBM

query_GBM <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR XML"
)

GDCdownload(query_GBM)
follow_up_GBM <- GDCprepare_clinic(query_GBM,"follow_up")
#colnames(follow_up_GBM)
data_GBM <- select(follow_up_GBM, 'bcr_patient_barcode', 'lost_follow_up', 'vital_status', 'days_to_last_followup', 'days_to_death')
data_GBM['disease'] <- "TGCA_GBM"

#recolher os dados de LGG

query_LGG <- GDCquery(
  project = "TCGA-LGG",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR XML"
)
GDCdownload(query_LGG)
follow_up_LGG <- GDCprepare_clinic(query_LGG,"follow_up")
data_LGG <- select(follow_up_LGG,'bcr_patient_barcode', 'lost_follow_up', 'vital_status', 'days_to_last_followup', 'days_to_death')
data_LGG['disease'] <- "TGCA_LGG"

# Juntar os dados das duas doenças
data <- rbind(data_LGG, data_GBM)

#----------------------------limpar os dados

# colocar NA onde não há valor
data[data==""] <- NA

#table(data$vital_status) <- resumo de uma linha

# alterar o vital_status para 1 em caso de Dead e 0 se Alive
data <- data %>%
  mutate(vital_status = ifelse(vital_status == "Alive",0,1))
#data$vital_status <- as.numeric(data$vital_status)

# caso o vital status seja NA, days to death seja NA, e perca de follow up é NA, considerar o lost follow up YES
data$lost_follow_up[is.na(data$vital_status) & is.na(data$days_to_death) & is.na(data$lost_follow_up) ] <- "YES" 

# caso o vital status seja NA, days to death seja NA, e perca de follow up é YES, considerar o vital status Alive
data$vital_status[is.na(data$vital_status) & is.na(data$days_to_death) & data$lost_follow_up == 'YES' ] <- "0" 

# criar a coluna days, assumindo os dias até à morte em caso de Dead e os dias até ao último follow up em caso de Alive
data <- data %>% 
  mutate(days = case_when(data$vital_status == 1 ~ data$days_to_death, data$vital_status == 0 ~ data$days_to_last_followup))

# remover linhas que não tenham nem days_to_followup nem days_to_death
data <- data[!is.na(data$days), ]

# ordenar os dados de acordo com a coluna days
data <- data[order(data$days), ]

# apagar as linhas cujo valor da coluna days era negativo
data <- data[data$days >= 0, ]

# colocar o vital_status como númerico
data$vital_status <- as.numeric(data$vital_status)

# Keep only distinct rows based on the "bcr_patient_barcodet" column
data <- distinct(data, bcr_patient_barcode, .keep_all = TRUE)

#----------------------------Curvas de sobrevivência

# Kaplan-Meier sem separação por doença
fit_all <- survfit(Surv(data$days, data$vital_status) ~ 1)
#summary(fit)
ggsurvplot(fit_all, data = data)

# Kaplan-Meier com a separação por doença
fit_disease <- survfit(Surv(data$days, data$vital_status) ~ data$disease, data = data)
survdiff(Surv(data$days, data$vital_status) ~ data$disease, data = data)
ggsurvplot(fit_disease, data = data)

#----------------------------Comparação da ID dos pacientes

# deste código
#barcode_GBM <- select(data[data$disease == 'TGCA_GBM', ], 'bcr_patient_barcode')
#barcode_LGG <- select(data[data$disease == 'TGCA_LGG', ], 'bcr_patient_barcode')
#barcode<- select(data, 'bcr_patient_barcode')
barcode<- data[['bcr_patient_barcode']]

write.csv(barcode, file = "barcode_survival.csv")

# do código do Gabriel
library(readr)
Data_Gabriel <- read_delim("Glioma_merged_classification_mutations_survival.csv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

barcode_Gabriel<- Data_Gabriel[['Patient_ID']]


#diagrama de Venn

if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")

library("ggVennDiagram")
# Default plot
ggVennDiagram(x <- list(barcode, barcode_Gabriel),
              category.names = list("Beatriz", "Gabriel"))

