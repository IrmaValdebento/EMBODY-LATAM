# Project EMBODY-LATAM (MOSAIC)

## Social Adversity, Epigenetic Aging, and Oral Cancer in Latin America

This repository contains the reproducible code for the **Two-Step Mendelian Randomization (MR)** analysis supporting the EMBODY-LATAM grant proposal.

### Project Overview
The study aims to establish a causal biological framework linking social adversity to oral carcinogenesis. Specifically, we investigate whether genetic liability to social disadvantage accelerates biological aging (PhenoAge/Hannum Age), thereby increasing the risk of Oral Squamous Cell Carcinoma (OSCC).

### Repository Contents
* **`01_analysis_pipeline.R`**: The main R script performing the two-step MR analysis:
    1.  **Step 1:** Education (Social Proxy) $\to$ Epigenetic Clocks.
    2.  **Step 2:** Epigenetic Clocks $\to$ Oral Cancer.
    3.  **Validation:** Education $\to$ Smoking Behavior.

### Requirements
The analysis runs in R and requires the following packages:
* `TwoSampleMR`
* `ieugwasr`
* `dplyr`
* `ggplot2`

### Note on Reproducibility & AI
* **Data Sources:** Public GWAS summary statistics from the OpenGWAS database (IEU).
* **Code Transparency:** The analysis pipeline was optimized with the assistance of LLM tools to ensure code cleanliness, reproducibility, and adherence to formatting standards. All statistical methodologies were selected by the principal investigator, and all outputs have been manually validated for accuracy.

### Contact
For questions regarding the methodology or preliminary data:
 Irma Valdebenito
