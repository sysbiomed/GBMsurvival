## README for Glioma Biomarkers Analysis

### Project Overview

This project focuses on the development of **prognostic models for glioblastoma** by integrating **machine learning-based feature selection methods** with **regularised Cox regression**. The study leverages **mRNA expression and survival data** from **The Cancer Genome Atlas (TCGA)**, following the **WHO-2021 glioma classification guidelines**.

The analysis involves:\
- **Feature selection methods**, including **Boruta and Random Survival Forests (RSF)**, to enhance model interpretability.\
- **Regularised Cox regression models**, such as **LASSO, Elastic Net, and Ridge**, for survival prediction.\
- **Network-based regularisation techniques (HubCox and OrphanCox)** to incorporate biological network information into prognostic modelling.

By systematically evaluating these methodologies, this project aims to improve the identification of **reliable glioblastoma biomarkers** and contribute to precision oncology.

### Repository Structure

-   **R Scripts**:

    -   `Code_16.Rmd`: The main R Markdown file containing the entire analysis workflow.

    -   add network-regularisation script

-   **Data**:

    -   `GDCdata/`: Directory where TCGA data is downloaded and stored.

-   **Output**:

    -   `Results/`: Directory where results such as tables, plots, and models are stored.

### Dependencies

This project uses the `renv` package to manage dependencies. To ensure reproducibility, the exact package versions used in this analysis are recorded in the `renv.lock` file.

### Setup Instructions

1.  **Clone the Repository**:

    ``` bash
    git clone https://github.com/BeatrizNL/Monet.git
    cd Monet
    ```

2.  **Restore the R Environment**:

    -   Make sure R is installed on your system.

    -   Open R or RStudio in the project directory.

    -   Run the following command to restore the environment:

        ``` r
        renv::restore()
        ```

3.  **Run the Analysis**:

    -   Open the `Code_16.Rmd` file in RStudio.

    -   Before running the analysis, set the following variables:

        -   **Alpha for model fitting**:
            -   *Section*: 'Fit and Explore the Models (train ≠ test)' → Variable: `best_alpha`
            -   *Section*: 'Fit and Explore the Models (train = test)' → Variable: `best_alpha`
            -   *Options*: Any value from **0 to 1**
        -   **Select the model for detailed evaluation**:
            -   *Section*: 'Analyse the Best Model (train ≠ test)' → Variable: `best_model_index`
            -   *Section*: 'Analyse the Best Model (train = test)' → Variable: `best_model_index_all_equal`
            -   *Options*:
                -   `1` → Multivariate Cox regression\
                -   `2` → Multivariate Cox regression with features from the RSF model\
                -   `3` → Multivariate Cox regression with features from the Boruta model

    -   Knit the file to generate the HTML report:

        ``` r
        rmarkdown::render("Code_16.Rmd")
        ```

### Project Workflow

### **Project Workflow**

1.  **Setup Environment & Load Dependencies**
    -   Restore the R environment using `renv::restore()`.\
    -   Load necessary R packages for data handling, survival analysis, machine learning models, and enrichment analysis.
2.  **Data Download & Preprocessing**
    -   Retrieve gene expression and survival data from The Cancer Genome Atlas (TCGA) using `TCGAbiolinks`.\
    -   Prepare and clean survival data, ensuring correct formatting and missing data handling.\
    -   **Update classification labels** using `SIMPLIFIED_CLASSIFICATION_TCGA_2016_2021.csv`.\
    -   Integrate the updated classification into the survival dataset.
3.  **Feature Selection Methods**
    -   Apply **univariate Cox regression** to pre-select survival-associated genes.\
    -   Implement **machine learning-based selection methods**:
        -   **Boruta**.\
        -   **Random Survival Forest (RSF)**.
4.  **Prognostic Model Training & Evaluation**
    -   Train **regularised Cox regression models**:
        -   **LASSO**, **Elastic Net**, and **Ridge** using `glmnet`.\
        -   **HubCox and OrphanCox** for network-based regularisation.\
    -   Evaluate model performance across **10 replicates** using:
        -   **Concordance index (C-index)**\
        -   **Integrated Brier Score (IBS)**\
        -   **Feature count & convergence rates**
5.  **Risk Stratification & Kaplan-Meier Survival Analysis**
    -   Categorise patients into **high-risk and low-risk groups**.\
    -   Generate **Kaplan-Meier survival curves** to assess model stratification performance.\
    -   Evaluate statistical significance using **log-rank tests**.
6.  **Feature Consistency & Interpretability Analysis**
    -   Assess overlap in selected features across replicates using **UpSet plots**.\
    -   Compare gene selection stability between Boruta and RSF-based models.
7.  **Functional Enrichment Analysis**
    -   Perform **Gene Ontology (GO)** and **KEGG pathway enrichment** on selected biomarkers.\
    -   Identify biological processes and pathways associated with glioblastoma prognosis.

------------------------------------------------------------------------

### **Additional Notes**

-   The workflow allows flexible **disease selection** (e.g., glioblastoma vs astrocytoma).\
-   **Gene expression data can be filtered to include only protein-coding genes, lncRNA, or miRNA.**\
-   Outputs (plots, tables, models) are stored in the **Results/** directory.

### How to Cite

If you use this analysis in your research, please cite the following paper:

Leitão, B.N., Veríssimo, A., Carvalho, A.M., & Vinga, S. (2025).\
*"Enhancing Prognostic Signatures in Glioblastoma with Feature Selection and Regularised Cox Regression."*\
**Genes**. DOI: [Insert DOI here]

### Support

For any questions or issues, please open an issue in the repository or contact the authors.
