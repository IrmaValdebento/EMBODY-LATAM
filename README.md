# Project EMBODY-LATAM (MOSAIC)

## Social Adversity, Epigenetic Aging, and Oral Cancer in Latin America

This repository contains the reproducible code for the **Two-Step Mendelian Randomization (MR)** analysis supporting the EMBODY-LATAM grant proposal.

### Project Overview
The study aims to establish a causal biological framework linking social adversity to oral carcinogenesis. Specifically, we investigate whether genetic liability to social disadvantage accelerates biological aging (PhenoAge/Hannum Age), thereby increasing the risk of Oral Squamous Cell Carcinoma (OSCC).

### Repository Contents
This analysis is divided into two main pipelines:

* **`01_analysis_pipeline.R`**: The main R script performing the biological Two-Step MR analysis:
    1.  **Step 1:** Education (Social Proxy) $\to$ Epigenetic Clocks (PhenoAge/Hannum).
    2.  **Step 2:** Epigenetic Clocks $\to$ Oral Cancer Risk.

* **`02_behavioral_validation.R`**: The validation script for the behavioral pathway:
    * **Validation:** Education $\to$ Smoking Intensity (Cigarettes per Day).
    * *Purpose:* To confirm that the genetic instruments for education accurately replicate known behavioral patterns in this population context.

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
* **Principal Investigator:** Irma Valdebenito
