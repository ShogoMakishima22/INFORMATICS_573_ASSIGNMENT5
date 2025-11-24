# ADVANCED_R ASSIGNMENT_5
# NAME: VENKATESH JOSHI


# 🧬 Differential Gene Expression Analysis (Tumor vs Normal)

This repository contains a fully automated R-based pipeline for performing a complete differential gene expression workflow using **tumor vs. normal microarray samples**. The script dynamically adapts to your local directory structure and input files as long as the required files are present.

---

## 📁 Required Input Files

Ensure the following files are placed in the working directory (or any folder you specify):

* `Gene_Expression_Data.xlsx`
* `Gene_Information.csv`
* `Sample_Information.tsv`

The script automatically validates the presence of these files before execution.

---

## ⚙️ How to Use

### 1️⃣ Clone the Repository(Sample)

```bash
git clone https://github.com/your-username/differential-gene-expression.git
cd differential-gene-expression
```

### 2️⃣ Update Data Directory

In the script, update this line to your local directory path:

```r
data_dir <- "YOUR/DIRECTORY/PATH/"
```

Example:

```r
data_dir <- "D:/INFORMATICS_573/"
```

### 3️⃣ Run the Script

```r
source("DEG_pipeline.R")
```

---

## 🧪 Workflow Overview

The pipeline performs the following steps automatically:

### ✅ 1. Package Management

Automatically installs and loads:

* readxl
* readr
* dplyr
* tidyr
* stringr
* ggplot2
* pheatmap

---

### ✅ 2. Data Import & Validation

* Reads gene expression matrix
* Imports gene annotation data
* Loads GSM sample metadata
* Confirms file integrity and structure

---

### ✅ 3. Sample Metadata Cleaning

* Splits `tumor\tpatient: X` format
* Extracts patient ID
* Creates standardized labels: `tissue_patient`
* Renames matrix columns accordingly

---

### ✅ 4. Expression Profiling

Calculates:

* Mean Tumor Expression per gene
* Mean Normal Expression per gene

---

### ✅ 5. Fold Change Calculation

Formula used:

```r
log2((Tumor - Normal) / Normal)
```

Significance threshold:

```
|log2FC| > 5
```

---

### ✅ 6. DEG Annotation

* Merges gene info
* Classifies upregulation:

  * Tumor-up
  * Normal-up
* Filters chromosomes: 1–22, X, Y, MT

---

## 📊 Visualizations Generated

### 📍 Chromosome Distribution

* Total DEGs per chromosome
* Tumor-up vs Normal-up bar plots

### 📍 DEG Percentage Plot

Displays percentage of up/downregulated genes

### 🔥 Heatmaps

* Top 500 most variable genes
* Hierarchical clustering
* Teal-lime scientific palette
* Tumor vs Normal segregation

---

## 📈 Outputs

The script produces:

✔ DEG chromosome distribution plots
✔ Up/downregulated percentage plots
✔ Clustered heatmap of variable genes
✔ Enhanced clustermap
✔ Console summary statistics
✔ Preview of annotated DEGs

All plots are saved automatically in the `results/` directory.

---

## 🧠 Biological Insights

This workflow reveals:

* Clear tumor vs normal clustering
* Significant DEGs across multiple chromosomes
* Strong expression patterns
* Biologically meaningful fold changes

---

## 📂 Project Structure

```
📦 differential-gene-expression
 ┣ 📜 DEG_pipeline.R
 ┣ 📄 README.md
 ┣ 📂 data/
 ┃ ┣ Gene_Expression_Data.xlsx
 ┃ ┣ Gene_Information.csv
 ┃ ┣ Sample_Information.tsv
 ┗ 📂 results/
```

---

## 🔄 Auto-Adjusting Feature

✔ Automatically detects:

* Any valid working directory
* File presence and structure
* Sample naming formats

✔ Adapts dynamically to:

* Different dataset sizes
* Varying sample labels
* Alternate patient counts

---

## ✅ Example Console Output

```
✔ Packages Loaded
✔ Input Files Verified
✔ Sample Labels Processed
✔ Fold Change Computed
✔ 1,245 Differentially Expressed Genes Found
✔ Visualizations Generated
```
