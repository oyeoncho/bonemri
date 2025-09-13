# Pelvic Bone MRI Survival – Code Repository

> End-to-end code and notebooks for **deep learning–based radiomics of pelvic bone T1‑weighted MRI** and survival prediction with a **Mixture Stretched‑Exponential (Mixture‑SE)** model.

This README summarizes how to set up the environment, prepare data, reproduce the main experiments (Model A / Model B), select robust features across BEiT variants, run external validation, and generate visualizations.

> 📄 Related manuscript: “Deep Learning–Based Radiomics of Pelvic Bone T1‑weighted MRI for Cervical Cancer Survival Prediction Using a Mixture Stretched‑Exponential Model.”  (Please cite when using this code.)

---

## 1) Repository Map (top-level)

```
code_paper/
├─ environment.json
├─ environment.txt
├─ nvidia-smi.txt
├─ pip-freeze.txt
├─ external/                    # External validation & cohort harmonization (R, ipynb)
│  ├─ ex_analysis.R
│  ├─ ex_integrate.R
│  ├─ ex_test.ipynb
│  ├─ preprocessing.ipynb
│  └─ ...
├─ feature_selection/           # Feature recurrence & refinement across runs/variants (R)
│  ├─ selection1.R
│  └─ selection2.R
├─ generalizable_analysis/      # Cross-variant common-feature analysis (R)
│  └─ overlap_features_analysis.R
├─ each_group/                  # Per‑group utilities (R)
│  └─ ...
├─ interpretation/              # Imaging attention & molecular links (R)
│  └─ interpretation.R
├─ IXI_pretrain/                # Domain-adaptive pretraining utilities / notes (ipynb)
│  └─ _pretrained_models.ipynb
├─ make_cache/                  # Caching helpers for preprocessing/features (py)
│  └─ make_cache.py (if present)
├─ beit0.py                     # Model A config: native 224×224, domain-adapted BEiT
├─ beit.py                      # Model A config: 196→224 padded, domain-adapted BEiT
├─ beit_resize.py               # Model A config: 196→224 resized, domain-adapted BEiT
├─ beit0_o/                     # ImageNet-only pretrain (no domain adaptation), native
├─ beit_o/                      # ImageNet-only, padded
├─ beit_resize_o/               # ImageNet-only, resized
├─ beit1/                       # Inferior 20-slice variant (domain-adapted)
├─ model_A_backbone_extract/    # Slice encoder & feature extraction
│  └─ beit/
│     ├─ beit_extract.py        # Extract 768-dim slice features; average→patient vector
│     └─ beit.py                # Model A core (blocks 4–11 fine-tuned)
├─ model_B/                     # Survival heads
│  ├─ model_B.ipynb             # Mixture‑SE training/eval (image / clinical / combined)
│  └─ cox_model_B_analysis.R    # CoxPH benchmarking
├─ preprocessing/               # MRI & mask preprocessing (ipynb)
│  └─ preprocessing.ipynb
├─ synthetic_benchmark/         # Simulation study
│  ├─ synthetic_benchmark.ipynb
│  └─ synthetic_analysis.R
├─ visualization/               # Attention maps & activation differences (py)
│  ├─ activation_diff.py
│  ├─ patch_level_attention.py
│  ├─ simlav.py / vis.py (if present)
│  └─ ...
└─ README.md
```

> The exact filenames may differ slightly by commit; use the closest variant if you don’t see an exact match.

---

## 2) Environment

- **Python** ≥ 3.10, **CUDA** (optional, recommended)
- Install from snapshot:
  ```bash
  # suggested: create a clean env first, then:
  pip install -r environment.txt
  # or replicate fully:
  pip install -r pip-freeze.txt
  ```
- For GPU diagnostics see `nvidia-smi.txt`.
- R scripts require **R ≥ 4.2** with the tidyverse ecosystem.

---

## 3) Data Layout (expected)

- **Inputs**
  - Pretreatment **T1‑weighted pelvic MRI** per patient.
  - Corresponding **pelvic bone masks** (NIfTI). (Masks can be generated via the paper’s pseudo‑labeling pipeline; see notes in `preprocessing/`.)
- **Preprocessing** (performed in `preprocessing/preprocessing.ipynb`):
  - N4 bias correction, Z‑score normalization.
  - Resample to **1.0×1.0×1.0 mm**.
  - Crop/pad to **[128, 224, 224]** (or **[128, 196, 196]** depending on variant).
  - Extract **20 slices** (central or inferior set by variant) **within the bone mask**.

> External cohorts with heterogeneous FOVs should be resampled to a fixed **physical FOV 224×224×128 mm³** to ensure complete pelvic coverage (see `external/` notes).

---

## 4) BEiT Variants (Model A)

| Variant dir/file | Resize strategy | Pretraining |
|---|---|---|
| `beit0.py`       | **Native 224×224** (no resize) | **Domain‑adaptive** (ImageNet→IXI T1) |
| `beit.py`        | **196→224 padded** (no stretch) | **Domain‑adaptive** |
| `beit_resize.py` | **196→224 resized** (stretched) | **Domain‑adaptive** |
| `beit0_o/`       | Native 224×224 | ImageNet‑only |
| `beit_o/`        | 196→224 padded | ImageNet‑only |
| `beit_resize_o/` | 196→224 resized | ImageNet‑only |
| `beit1/`         | **Inferior 20 slices** (native 224) | Domain‑adaptive |

Each Model A run encodes 20 masked slices with BEiT → **768‑dim** slice features → **patient‑level average** vector. During fine‑tuning, only **BEiT blocks 4–11** are updated to preserve pretraining.

**Run Model A (example)**

```bash
# Example: domain‑adapted, native 224
python model_A_backbone_extract/beit/beit_extract.py   --variant beit0   --images /path/to/mri/   --masks  /path/to/masks/   --outdir ./features/beit0/run_XX/
```

- Repeat **30×** (different seeds) per variant. Each run saves:
  - patient‑level features (CSV or NPZ),
  - split indices (train/val),
  - logs/plots.

---

## 5) Feature Refinement (recurrence across runs)

1. For each run, compute survival correlation and keep **top‑100** backbone features.
2. Stack across 30 runs (**3,000** entries) and form **nX** groups = features recurring ≥ **X** times (X=4..9).
3. Use **`feature_selection/selection1.R`** and **`selection2.R`** to materialize refined sets.
4. **`generalizable_analysis/overlap_features_analysis.R`** derives **common features across BEiT variants** (nested: **n4 ⊃ n5 ⊃ n6 ⊃ n7**).

---

## 6) Survival Modeling (Model B) & Baselines

- **Model B (Mixture‑SE)** — Train on refined features (image‑only, clinical‑only, or combined):
  - Notebook: `model_B/model_B.ipynb`
  - Repeats **30×** (70/30 MC‑CV). Primary metric: **mean C‑index** over {12,24,36,48,60,72} months.
- **CoxPH baseline** — `model_B/cox_model_B_analysis.R`
- Outputs:
  - Per‑run metrics (AUC, C‑index), distribution plots, and tables.

**CLI‑style pseudocode**

```bash
# Mixture‑SE (image‑only)
jupyter nbconvert --to notebook --execute model_B/model_B.ipynb   --TagRemovePreprocessor.remove_input_tags='["long-run"]'   --output outputs/model_B_mixtureSE_image_only.ipynb
```

---

## 7) External Validation

- Data harmonization & follow‑up definitions: `external/` (R & ipynb).
- Ensure resampling to **224×224×128 mm³** physical FOV before slice extraction.
- Reuse internal **run‑specific indices** (n4–n7) to extract matching features from external patients.
- Evaluate with the same MC‑CV protocol (30×) and report **C‑index/AUC**.

---

## 8) Visualization & Interpretation

- **Patch‑level attention**: `visualization/patch_level_attention.py`
- **Group activation differences** (e.g., >60‑month survivors vs <12‑month deaths): `visualization/activation_diff.py`
- **Molecular associations** (exosomal RNA / CBC): `interpretation/interpretation.R`

---

## 9) Synthetic Benchmark (optional but recommended)

- Reproduce Fig. 2‑style comparisons (Mixture‑SE vs CoxPH/WeibullAFT/DeepSurv):
  - `synthetic_benchmark/synthetic_benchmark.ipynb`
  - `synthetic_benchmark/synthetic_analysis.R`

---

## 10) End‑to‑End Quickstart

```bash
# 0) Create env
pip install -r environment.txt

# 1) Preprocess MRI + masks
jupyter lab preprocessing/preprocessing.ipynb

# 2) Train Model A (30 runs) for chosen BEiT variants
python model_A_backbone_extract/beit/beit_extract.py --variant beit0 ...

# 3) Select recurrent features (n4–n7) per variant
Rscript feature_selection/selection1.R
Rscript feature_selection/selection2.R

# 4) Derive cross‑variant common features
Rscript generalizable_analysis/overlap_features_analysis.R

# 5) Train Model B (Mixture‑SE) on refined/common features
jupyter lab model_B/model_B.ipynb

# 6) Benchmark with CoxPH
Rscript model_B/cox_model_B_analysis.R

# 7) External validation
jupyter lab external/ex_analysis.R  # or run via RStudio
```

---

## 11) Key Results (for orientation)

- **Internal validation (image‑only)**: mean **C‑index ≈ 0.829**, **AUC ≈ 0.852** across runs; adding clinical variables did **not** significantly improve overall survival prediction.
- **External validation (TCGA‑CESC, CCRT‑only)**: **C‑index ≈ 0.703**, **AUC ≈ 0.732** (image‑only, n4–n7); **n7** (2 features) reached **C‑index ≈ 0.744**, **AUC ≈ 0.803**, comparable to image+clinical.

See manuscript for full statistics and confidence intervals.

---

## 12) Tips & Troubleshooting

- **Slice coverage**: If the left/right pelvis is cut off, increase FOV or adjust centering before cropping to [128,224,224].
- **Reproducibility**: Fix seeds per run; log train/val splits.
- **Overfitting**: Prefer **n6/n7** in very small external cohorts to keep events‑per‑variable reasonable.
- **Domain adaptation**: Models with **IXI T1** domain‑adaptive pretraining tend to yield more reproducible features and cleaner attention maps.

---

## 13) Citation

If you use this code, please cite the accompanying manuscript.

Cho O, El Naqa I. *Deep Learning–Based Radiomics of Pelvic Bone T1‑weighted MRI for Cervical Cancer Survival Prediction Using a Mixture Stretched‑Exponential Model.*

---

## 14) License

Unless specified elsewhere, this code is released for **research use**. Please contact the authors for other uses.

---

## 15) Acknowledgments

Thanks to collaborators and the institutions supporting this work.
