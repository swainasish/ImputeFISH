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
## Tutorials and Reproducibility section 
### Tutorial-1 (On MERFISH dataset) 
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1j2t6kzDEaMNihiTXc6lFZsK8MP4bGuq8?usp=sharing)
## Short-Usage 
```bash
from ImputeFISH import GeneEnhancement
impfish_res = GeneEnhancement(sc_mat,sp_mat,g_list,sp_x,sp_y,
                                 reference_selection=True,
                                 batch_info=batch_info)
```
Argument Details \
sc_mat: scRNA-seq reference expression matrix (in pandas DataFrame format) \
sp_mat: Spatial expression matrix (in pandas DataFrame format) \
g_list: Target genes to impute (present in scRNA-seq, absent in spatial) \
sp_x, sp_y: Spatial coordinates of cells/spots \
reference_selection (optional): Enables co-expression-based selection of optimal scRNA-seq reference \
batch_info (optional): Batch labels for scRNA-seq data (if multiple datasets are used) 

## Example
```bash
import numpy as np
import pandas as pd
from ImputeFISH import GeneEnhancement

# -----------------------------
# 1. scRNA-seq reference matrix
# -----------------------------
# Shape: cells × genes
sc_mat = pd.DataFrame(
    np.random.rand(200, 100),
    columns=[f"Gene{i}" for i in range(100)]
)

# ---------------------------------------
# 2. Spatial transcriptomics matrix
# ---------------------------------------
# Shape: spots/cells × measured genes
sp_mat = pd.DataFrame(
    np.random.rand(100, 50),
    columns=[f"Gene{i}" for i in range(50)]
)

# ---------------------------------------
# 3. Target genes to impute
# ---------------------------------------
# Must be present in sc_mat but not measured in spatial data
g_list = ["Gene88", "Gene89"]

# ---------------------------------------
# 4. Spatial coordinates
# ---------------------------------------
sp_x = np.random.rand(100)
sp_y = np.random.rand(100)

# ---------------------------------------
# 5. Optional: batch information (scRNA-seq)
# ---------------------------------------
batch_info = np.random.choice(["Batch1", "Batch2"], 200)

# ---------------------------------------
# 6. Run ImputeFISH
# ---------------------------------------
impfish_res = GeneEnhancement(
    sc_mat,
    sp_mat,
    g_list,
    sp_x,
    sp_y,
    reference_selection=True,
    batch_info=batch_info
)

# ---------------------------------------
# 7. Output
# ---------------------------------------
# impfish_res: imputed expression matrix for target genes
print(impfish_res.head())
```
