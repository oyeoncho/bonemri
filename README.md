# Prognostic Impact of Pelvic Bone MRI Features in Cervical Cancer Undergoing CCRT  
### Code Repository (Main + Supplementary Analyses)

This repository contains the complete implementation used in the manuscript:

**“Prognostic Impact of Pelvic Bone Magnetic Resonance Imaging Features in Patients with Cervical Cancer Undergoing Concurrent Chemoradiotherapy (CCRT)”**

It includes all analyses described in the Main Text and Supplementary Materials, including internal Monte Carlo validation, representative survival analyses, external validation (Supplementary), biological correlation analyses, visualization, decision curve analysis, and supplementary benchmarking experiments.

---

# 📁 Repository Structure

| Folder | Purpose | Corresponding Figures/Tables |
|--------|---------|-----------------------------|
| `preprocessing/`, `make_cache/` | MRI preprocessing, pelvic bone masking, slice caching | Supplementary Methods |
| `backbone_extract/` | BEiT-based backbone feature extraction | Supplementary Methods |
| `feature_selection/` | Stability grouping (n4–n7), recurrence statistics | Tables S3–S4 |
| `stable_feature_internal/` | Internal Monte Carlo survival modeling | Fig. S2–S3, Tables S5–S6 |
| `external/` | TCGA-CESC external validation | **Fig. S4, Table S9** |
| `km_cox/` | Representative model–based Kaplan–Meier & Cox analyses | **Fig. 1, Fig. 2, Table 1** |
| `visualization/` | Patch-level attention & spatial activation maps | **Fig. 3** |
| `biological_relevance/` | Correlation with exosomal RNAs & CBC indices | **Fig. S7–S8** |
| `dca/` | Decision curve analysis (5-year OS) | **Fig. S9** |
| `synthetic_benchmark/` | Mixture-SE benchmarking vs CoxPH/DeepSurv/Weibull | Supplementary Methods |
| `generalizable_analysis/` | Backbone robustness & pretraining comparisons | Supplementary Methods |
| `supplement_information/` | Scripts for Supplementary figures and tables | Fig. S1–S10, Tables S1–S9 |

---

# 🔧 Analysis Overview

## 1. MRI Preprocessing & Pelvic Bone Masking
Implements:

- N4 bias-field correction  
- Intensity standardization  
- 1-mm isotropic resampling  
- Pelvic cropping  
- Bone mask generation  
- Central / lower-pelvic slice caching  

Outputs bone-masked T1-weighted MRI slices for feature extraction.

---

## 2. BEiT Feature Extraction (Backbone)

- BEiT Base Patch16-224 vision transformer backbone  
- Pretrained on large natural-image datasets  
- Domain-adaptively tuned where specified  
- Frozen during survival modeling  
- Produces 768-dimensional slice-level embeddings  

Final stable imaging biomarkers: **Feature IDs 436 and 519 (n7 group)**.

---

## 3. Stability-Based Feature Selection (n4–n7)

- 30 independent 70:30 Monte Carlo splits  
- Recurrence-based grouping (n4–n7)  
- n7 features (436 & 519) selected as final stable biomarkers  
- No re-selection in validation or external cohorts  

Reproduces Tables S3–S4.

---

## 4. Internal Survival Modeling (Supplementary)

Implements the two-component Mixture Stretched-Exponential (Mixture-SE) model:

- Evaluated across 30 Monte Carlo runs  
- Time-dependent C-index and AUC (12–72 months)  

Reproduces:

- **Fig. S2–S3**
- **Tables S5–S6**

---

## 5. External Validation (Supplementary)

Location: `external/`

- Identical preprocessing and backbone configuration  
- n7 features applied without re-estimation  
- Internal models applied directly (no retraining)  
- 30 run-wise external C-index/AUC estimates  

Reproduces:

- **Fig. S4**
- **Table S9**

---

## 6. Representative Model: Kaplan–Meier & Cox Analyses (Main Text)

Location: `km_cox/`

Representative model selection:

- Backbone with highest mean external C-index  
- Internal run with validation C-index closest to median  

Using predefined n7 features:

- Compute 60-month predicted mortality risk  
- Kaplan–Meier survival curves (OS, PFS)  
- Log-rank testing  
- Multivariable Cox regression (age, stage, histology, imaging risk)  
- Stage × Imaging Risk 4-group stratification  

Reproduces:

- **Fig. 1 (Kaplan–Meier)**
- **Fig. 2 (Cox regression)**
- **Table 1**

---

## 7. Patch-Level Attention Visualization (Main Text)

Location: `visualization/`

- Attention maps over bone-masked slices  
- Long-term survivor vs early-death comparison  
- Hotspot frequency maps  

Reproduces:

- **Fig. 3**

---

## 8. Biological Correlation Analyses (Supplementary)

Location: `biological_relevance/`

- Correlation of n7 features with:
  - miR-574-3p  
  - LINC01003  
  - ACOT9  
  - CBC-derived indices  

Reproduces:

- **Fig. S7–S8**

---

## 9. Decision Curve Analysis (Supplementary)

Location: `dca/`

- 5-year overall survival net benefit  
- Image-only vs Image+Clinical vs Clinical-only  

Reproduces:

- **Fig. S9**

---

# 🔁 Reproducing Key Results

| Output | Folder |
|--------|--------|
| Fig. 1 (Kaplan–Meier) | `km_cox/` |
| Fig. 2 (Cox regression) | `km_cox/` |
| Fig. 3 (Visualization) | `visualization/` |
| Fig. S2–S3 (Internal MC) | `stable_feature_internal/` |
| Fig. S4 (External validation) | `external/` |
| Fig. S7–S8 (Biological correlation) | `biological_relevance/` |
| Fig. S9 (DCA) | `dca/` |
| Table 1 | `km_cox/` |
| Tables S1–S9 | `supplement_information/` |

---

# 💻 Software

Analyses were performed using:

- Python (PyTorch, lifelines, scikit-learn)  
- R (survival, survminer, ggplot2)

Detailed parameters and run configurations are described in the Supplementary Methods.

---

# 📄 Data Availability

Processed datasets supporting the manuscript findings will be made available upon journal request.  
Original imaging data require institutional review board approval.  
Code is provided to ensure full reproducibility of the analyses described in the manuscript.