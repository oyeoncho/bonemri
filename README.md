# Prognostic Impact of Pelvic Bone MRI Features in Cervical Cancer Undergoing CCRT  
### Code Repository (Main + Supplementary Analyses)

This repository contains the complete implementation used in the manuscript  

**“Prognostic Impact of Pelvic Bone MRI Features in Patients with Cervical Cancer Undergoing Concurrent Chemoradiotherapy (CCRT)”**  

including all analyses described in the **Main Text and Supplementary Materials**:

- Internal Monte-Carlo validation and stability-based feature selection (n4–n7, final n7 features 436 & 519)
- External validation in the **TCGA-CESC** cohort
- Correlation of n7 imaging features with exosomal RNAs (miR-574-3p, LINC01003, ACOT9) and CBC indices
- Patch-level attention visualization over pelvic bone–masked slices
- Decision curve analysis (DCA) for 5-year overall survival
:contentReference[oaicite:0]{index=0}

---

## 📁 Repository Structure

| Folder | Purpose | Corresponding Figures/Tables |
|--------|---------|-----------------------------|
| `stable_feature_internal/` | Internal Monte-Carlo (30×) validation, n4–n7 feature recurrence, OS/DM/DM+LP/LP model performance & progression patterns | **Fig. 1**, **Tables S3–S5** |
| `biological_relevance/` | Correlation of n7 imaging features with plasma exosomal RNAs (miR-574-3p, LINC01003, ACOT9) and hematologic indices (CBC0, CBCmin, CBC1, NLR/PLR/LMR) | **Fig. 3**, **Fig. S1** |
| `external/` | TCGA-CESC external validation workflow, mirroring internal MC structure | **Fig. 2**, **Table S6** |
| `dca/` | Decision curve analysis for 5-year OS (net benefit vs treat-all / treat-none) | **Fig. 5** |
| `visualization/` | Patch-level attention and spatial activation maps for n7 features over pelvic bone–masked slices | **Fig. 4**, **Fig. S2** |
| `supplement_information/` | Full Supplementary Methods implementation and scripts for tables/figures | **Supplementary Figs. S1–S2**, **Tables S1–S6** |
| `backbone_extract/` | BEiT backbone feature extraction & domain-adaptive pretraining utilities | Supplementary Methods (Backbone description) |
| `preprocessing/`, `make_cache/` | MRI preprocessing, bone mask application, slice caching (central / lower pelvis) | Supplementary Methods (Image preprocessing) |
| `feature_selection/` | n4–n7 stability grouping, recurrence statistics, reproducibility across datasets/regions | **Tables S3–S4** |
| `synthetic_benchmark/` | Mixture stretched-exponential (Mixture-SE) benchmarking against CoxPH / DeepSurv / Weibull on synthetic data | Supplementary Methods (Model comparison) |
| `generalizable_analysis/` | Comparison of VT backbones / pretraining variants (BEiT & others) and robustness checks | Supplementary Methods (Backbone comparison) |

---

# 🔧 Analysis Overview

## 1. Preprocessing & Pelvic Bone Masking

**Location:**  
`preprocessing/`, `make_cache/`

Implements the MRI pipeline described in the Methods and Supplementary Materials:

- N4 bias-field correction  
- Intensity standardization (e.g., robust Z-score)  
- Resampling to 1-mm isotropic voxels  
- Cropping/padding to a standardized pelvic field of view  
- Slice caching (central or lower-pelvic slices per patient)  
- Pelvic bone mask generation (manual labels → pseudo-labels → DL propagation)

These scripts create the bone-masked T1W-MRI inputs used for BEiT feature extraction and subsequent survival modeling.

---

## 2. BEiT Feature Extraction (Backbone Model)

**Location:**  
`backbone_extract/`

Implements the **BEiT-based vision transformer** used in the manuscript:

- BEiT Base Patch16-224 backbone  
- Initialized from publicly available weights, with domain-adaptive tuning on MRI when specified  
- Frozen backbone during survival modeling; selected blocks can be fine-tuned in pretraining scripts  
- Produces 768-dim slice embeddings from **bone-masked pretreatment T1W-MRI**  
- Supports multiple preprocessing variants (native 224×224, resized / alternative pretraining schemes, lower-pelvic slice set)

These slice-level features are the basis for stability-based feature selection and the final n7 feature set (IDs **436** and **519**).

---

## 3. Stability-Based Feature Selection (n4–n7)

**Location:**  
`feature_selection/`, `stable_feature_internal/`

Implements the **Monte-Carlo (MC) stability analysis**:

1. 30 independent 70:30 train/validation splits per dataset (MC cross-validation).  
2. Within each split, survival-associated imaging features are screened.  
3. Recurrence counts across 30 runs define the stability groups: **n4**, **n5**, **n6**, **n7**.  
4. Across six datasets, features with recurrence 4–7 times form n4–n7 groups.  
5. Reproducibility is checked in the lower-pelvic slice analysis;  
   - n7 and n6 features are fully replicated,  
   - n4/n5 show partial replication.  
6. The **n7 features (IDs 436 and 519)** constitute the **final, most stable imaging biomarker set**, used in all downstream modeling (no re-selection in internal validation or external cohort).

Reproduces:

- Stability summaries and cross-dataset reproducibility  
- **Fig. 1 (part of panel A–C via performance by feature set)**  
- **Tables S3–S4**

---

## 4. Internal Survival Modeling (Mixture Stretched-Exponential)

**Location:**  
`stable_feature_internal/`

Implements the **two-component Mixture Stretched-Exponential (Mixture-SE) survival model**:

- Survival function modeled as a weighted mixture of exponential components  
- Captures heterogeneous risk and time-varying hazards  
- For each of the 30 MC runs:  
  - Slice-level parameters are predicted and aggregated to patient-level features.  
  - Mixture-SE model is trained using selected features (n7, or alternative nX groups in supplementary analyses).  
  - Time-dependent **C-index** and **AUC** are computed at 12–72 months.  
  - Mean values across time are used as run-level performance metrics.  
- Comparisons:  
  - **Clinical-only**  
  - **Image-only (n7 features)**  
  - **Image + Clinical**

Reproduces:

- Internal validation performance for **OS, DM, DM+LP, LP**  
- Progression pattern distribution by survival status (LP, DM, DM+LP)  
- **Fig. 1**  
- **Table S5**

---

## 5. External Validation (TCGA-CESC)

**Location:**  
`external/`

Implements the **TCGA-CESC external validation** as described in the manuscript:

- Includes only CCRT-treated cases with adequate follow-up.  
- Pretreatment T1W-MRI is processed with the **identical preprocessing pipeline** and frozen BEiT backbone used for the discovery cohort.  
- For each of the 30 internal runs:  
  - External slice-level features are extracted using the corresponding run-specific aggregation.  
  - The **predefined n7 features** (436, 519) are selected **without re-estimation**.  
  - Mixture-SE models trained on internal runs are applied **directly** to external data (no retraining).  
- Produces 30 run-wise external **C-index** and **AUC** values, mirroring the internal MC design.

Reproduces:

- External C-index distribution across 30 runs (Image-only, n7).  
- Comparison of **Image-only**, **Clinical-only**, and **Image + Clinical** models for the representative best run.  
- **Fig. 2**  
- **Table S6**

---

## 6. Biological Relevance: Exosomal RNAs & Hematologic Indices

**Location:**  
`biological_relevance/`

Assesses the **biological interpretation of n7 features** using the 41-patient subset with:

- Plasma exosomal RNA log₂ fold change (baseline → 2 weeks post-RT) for:  
  - **miR-574-3p**, **LINC01003**, **ACOT9**  
- Serial CBC metrics (CBC0, minimum during treatment, CBC1, and derived ratios):  
  - ANC, ALC, PLT, Hb, Mo  
  - NLR, PLR, LMR (0 and 1)

Pipeline:

1. Compute Pearson correlations between n7 features and biomarker variables.  
2. Select top variables per run (e.g., |r| ≈ 0.4–0.6).  
3. Aggregate across 30 runs → 150 feature–biomarker associations.  
4. Summarize recurrence frequencies to identify dominant axes.

Key outputs:

- Dominant **miR-574-3p–LINC01003–ACOT9** exosomal RNA axis associated with n7 features.  
- Lower and less consistent recurrence for CBC-derived indices (e.g., Hb, NLR).  

Reproduces:

- **Fig. 3**  
- **Fig. S1**

---

## 7. Decision Curve Analysis (DCA)

**Location:**  
`dca/`

Implements the DCA described in the manuscript:

- Uses n7 feature set from **Run 4** and a representative survival model from **Run 6** (external C-index closest to the median).  
- Computes predicted 5-year OS risk (T = 60 months) per patient.  
- Net Benefit is evaluated across threshold probabilities (approximately 0–0.8):  
  - **Image-only** model  
  - **Image + Clinical** model  
  - **Clinical-only** model  
  - “Treat all” and “Treat none” strategies  
- Provides separate analyses for:  
  - Internal validation cohort (n = 149)  
  - External TCGA-CESC cohort (n = 38)

Reproduces:

- Net benefit curves for 5-year OS  
- Demonstrates higher net benefit for image-based models in threshold ranges ~0.2–0.6  
- **Fig. 5**

---

## 8. Patch-Level Attention Visualization

**Location:**  
`visualization/`

Visualizes **patch-level attention and spatial activation differences** for n7 features:

- Generates bone-masked T1W-MRI slices for:  
  - Long-term survivors (>60 months)  
  - Early-death cases (<12 months)  
- Overlays BEiT-derived attention maps and activation patterns, showing:  
  - High attention localized to **posterior iliac bone marrow** and **upper sacral marrow**, regions rich in active hematopoietic tissue.  
  - Distinct spatial activation patterns between survival-extreme groups.  

Reproduces:

- **Fig. 4** (attention maps and activation differences)  
- **Fig. S2** (additional visualization details, if applicable)

---

## 9. Synthetic Benchmark & Generalizability Analyses

**Location:**  
- `synthetic_benchmark/`  
- `generalizable_analysis/`

These folders contain supplementary experiments:

- **Synthetic benchmark:**  
  - Compares the Mixture-SE model against CoxPH, DeepSurv, and Weibull-AFT models on synthetic data with nonlinear and time-varying hazards.  

- **Generalizable analysis:**  
  - Explores alternative backbones and pretraining variants to assess robustness of the n7 feature concept and overall modeling framework.

They support the **methodological context** described in the Supplementary Methods.

---

## 10. Reproducing the Manuscript Results

To reproduce the main figures and tables:

- **Fig. 1 & Tables S3–S5** → `stable_feature_internal/`, `feature_selection/`  
- **Fig. 2 & Table S6** → `external/`  
- **Fig. 3 & Fig. S1** → `biological_relevance/`  
- **Fig. 4 & Fig. S2** → `visualization/`  
- **Fig. 5** → `dca/`  
- **Supplementary tables/figures** → `supplement_information/`

Each subfolder contains R/Python scripts with comments indicating:

- Required input files (discovery cohort, TCGA-CESC, biomarker subset)  
- Output CSV/figure paths  
- Key hyperparameters (e.g., number of MC runs, thresholds, bootstrap iterations)

Please refer to the **Supplementary Methods** for detailed parameter choices and additional notes on preprocessing, model training, and evaluation.
