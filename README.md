## README for Glioma Biomarkers Analysis

### Project Overview

This project aims to identify and validate glioma biomarkers that can distinguish high-risk from low-risk patients. The analysis is based on gene expression data from The Cancer Genome Atlas (TCGA) and involves using multiple machine learning models, including Regularized Cox Regression, SIS, ISIS, Random Forest, Causal Forest, and Boruta.

### Repository Structure

-   **R Scripts**:

    -   `Code_16.Rmd`: The main R Markdown file containing the entire analysis workflow.

    -   `SIMPLIFIED_CLASSIFICATION_TCGA_2016_2021.csv`: A CSV file with simplified classification data used for analysis.

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

    -   Open the `Glioma_Biomarkers_Analysis.Rmd` file in RStudio.

    -   Knit the file to generate the HTML report:

        ``` r
        rmarkdown::render("Code_16.Rmd")
        ```

### Project Workflow

1.  **Data Download and Preparation**:

    -   TCGA data is downloaded using the `TCGAbiolinks` package.

    -   Survival and gene expression data are extracted and cleaned.

2.  **Normalization**:

    -   Gene expression data is normalized using the EdgeR package.

3.  **Model Fitting**:

    -   Multiple machine learning models are trained to identify biomarkers.

    -   The performance of each model is evaluated using the C-index.

4.  **Outlier Detection**:

    -   Outliers are identified using the rank product test, and models are re-evaluated after their removal.

5.  **Analysis of the Best Model**:

    -   The best model is selected based on the highest average C-index, and further analysis is performed to evaluate its performance.

### How to Cite

If you use this analysis in your research, please cite the following paper:

Leitão, B., & Vinga, S. (2024). "Glioma biomarkers to distinguish high-risk from low-risk patients."

### Support

For any questions or issues, please open an issue in the repository or contact the authors.
