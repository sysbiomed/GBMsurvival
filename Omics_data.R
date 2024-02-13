library(tidyverse)
library(readxl)
library(writexl)

original_data <- read_excel("Glioma - data consolidation.xlsx")

# Remove duplicated rows
clean_data <- original_data %>%
  distinct()

# Extract unique sample IDs
sample_ids <- unique(unlist(clean_data))

# Create a data frame with Sample_ID (as character)
data <- data.frame(Sample_ID = as.character(sample_ids))

# Use pivot_longer to reshape the data from wide to long
long_data <- clean_data %>%
  pivot_longer(cols = everything(), names_to = "Characteristic", values_to = "Sample_ID") %>%
  mutate(Value = 1)

# Use distinct to ensure unique combinations
wide_data <- long_data %>%
  distinct(Sample_ID, Characteristic, .keep_all = TRUE) %>%
  pivot_wider(names_from = Characteristic, values_from = Value, values_fill = 0)

# Replace NA values in Sample_ID using coalesce
wide_data <- wide_data %>%
  mutate(Sample_ID = coalesce(Sample_ID, "0"))


# Write to Excel
write_xlsx(wide_data, "samples_omics.xlsx")

