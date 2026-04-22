# ImputeFISH

**ImputeFISH: A Spatially Aware and Scalable Framework for Gene Imputation in Imaging-Based Spatial Transcriptomics**

---

## 🔍 Overview

**ImputeFISH** (*Imputation for Fluorescence In-Situ Hybridization–based spatial transcriptomics*) is a spatially aware and scalable framework designed for **gene imputation** and **atlas-scale transcriptome enhancement** in imaging-based spatial transcriptomics (IST) datasets such as MERFISH, Xenium, and seqFISH.

ImputeFISH integrates spatial context and scRNA-seq references to recover unmeasured genes, enabling more complete transcriptome representations and improved biological interpretation.

---

## 🚀 Key Features

- 🧠 **Reference Selection:** Automatically identifies the most compatible scRNA-seq reference using a **gene co-expression network–based compatibility filter**.  
- 🌐 **Spatial Awareness:** Incorporates neighborhood-level spatial information inspired by **graph neural network message passing**.  
- 🔄 **Domain Adaptation:** Applies **SVD-based domain alignment** to correct batch effects between IST and scRNA-seq datasets.  
- ⚡ **Scalability:** Efficiently scales to **millions of cells**, suitable for atlas-scale datasets.  
- 🧬 **Accurate Gene Imputation:** Uses **non-negative ridge regression** to predict unmeasured gene expression.  

---

## 📦 Installation

You can install **ImputeFISH** using pip:

```bash
pip install imputefish
```
## Usage 
```bash
from ImputeFISH import GeneEnhancement
impfish_res = GeneEnhancement(sc_mat,sp_mat,g_list,sp_x,sp_y,
                                 referece_selection=True,
                                 batch_info=batch_info)
```
Argument Details 
sc_mat: scRNA-seq reference expression matrix (in pandas DataFrame format)
sp_mat: Spatial expression matrix (in pandas DataFrame format)
g_list: Target genes to impute (present in scRNA-seq, absent in spatial)
sp_x, sp_y: Spatial coordinates of cells/spots
reference_selection (optional): Enables co-expression-based selection of optimal scRNA-seq reference
batch_info (optional): Batch labels for scRNA-seq data (if multiple datasets are used)
