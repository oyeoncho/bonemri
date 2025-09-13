import os
import torch
import torch.nn.functional as F
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

# ⚙️ 경로 설정
FEATURE_A_CSV = "./features/beit0_features/beit_base_patch16_224_backbone_features_run0.csv"
FEATURE_B_CSV = "./features/beit0_o_features/beit_base_patch16_224_backbone_features_run0.csv"
SAVE_DIR = "./feature_similarity"
os.makedirs(SAVE_DIR, exist_ok=True)

# 📥 Feature 불러오기 (CSV → DataFrame → Tensor)
df_a = pd.read_csv(FEATURE_A_CSV, index_col=0)
df_b = pd.read_csv(FEATURE_B_CSV, index_col=0)
feat_a = torch.tensor(df_a.values, dtype=torch.float32)
feat_b = torch.tensor(df_b.values, dtype=torch.float32)

# ✅ 유효성 검사
assert feat_a.shape == feat_b.shape, f"❌ Shape mismatch: {feat_a.shape} vs {feat_b.shape}"

# ✅ Cosine similarity 계산
def compute_feature_cosine(feat_a, feat_b):
    feat_a = F.normalize(feat_a, dim=0)
    feat_b = F.normalize(feat_b, dim=0)
    return (feat_a * feat_b).mean(dim=0).cpu().numpy()

# ✅ Pearson correlation 계산
def compute_feature_pearson(feat_a, feat_b):
    if feat_a.shape[0] < 2:
        print("⚠️ N=1이므로 Pearson correlation 생략 (NaN으로 채움)")
        return np.full(feat_a.shape[1], np.nan)
    corr = []
    for i in range(feat_a.shape[1]):
        a = feat_a[:, i].cpu()
        b = feat_b[:, i].cpu()
        r, _ = pearsonr(a, b)
        corr.append(r)
    return np.array(corr)

# 🔍 유사도 계산
cos_sim = compute_feature_cosine(feat_a, feat_b)         # [768]
pearson_corr = compute_feature_pearson(feat_a, feat_b)   # [768]

# ✅ 저장: 수치 결과 (npy)
np.save(f"{SAVE_DIR}/cosine_similarity.npy", cos_sim)
np.save(f"{SAVE_DIR}/pearson_correlation.npy", pearson_corr)

# ✅ 저장: 수치 결과 (CSV)
df_similarity = pd.DataFrame({
    "Feature_Index": np.arange(len(cos_sim)),
    "Cosine_Similarity": cos_sim,
    "Pearson_Correlation": pearson_corr
})
df_similarity.to_csv(f"{SAVE_DIR}/feature_similarity_values.csv", index=False)

# 📊 히스토그램 저장
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
sns.histplot(cos_sim, bins=40, kde=True)
plt.title("Cosine Similarity per Feature Index")
plt.xlabel("Cosine Similarity")

plt.subplot(1, 2, 2)
if np.isnan(pearson_corr).all():
    plt.text(0.5, 0.5, "No Pearson Correlation (N=1)", ha="center", va="center")
else:
    sns.histplot(pearson_corr[~np.isnan(pearson_corr)], bins=40, kde=True)
    plt.title("Pearson Correlation per Feature Index")
    plt.xlabel("Pearson Correlation")
plt.tight_layout()
plt.savefig(f"{SAVE_DIR}/feature_similarity_hist.png")
plt.close()

# 🔝 분산 기준 상위 feature index overlap
top_k = 20
top_a = np.argsort(np.var(feat_a.numpy(), axis=0))[-top_k:]
top_b = np.argsort(np.var(feat_b.numpy(), axis=0))[-top_k:]
overlap = len(set(top_a).intersection(set(top_b)))

# ✅ 요약 텍스트 저장
with open(f"{SAVE_DIR}/summary.txt", "w") as f:
    f.write(f"Top-{top_k} Feature Index Overlap:\n")
    f.write(f"Model A Top-{top_k}: {sorted(top_a.tolist())}\n")
    f.write(f"Model B Top-{top_k}: {sorted(top_b.tolist())}\n")
    f.write(f"Overlap Count: {overlap} / {top_k}\n")

# ✅ Top-K feature index CSV 저장
pd.DataFrame({
    "Top_A": pd.Series(sorted(top_a.tolist())),
    "Top_B": pd.Series(sorted(top_b.tolist()))
}).to_csv(f"{SAVE_DIR}/topk_indices.csv", index=False)

# 🔍 Cosine vs Pearson 산점도 저장
plt.figure(figsize=(6, 6))
if np.isnan(pearson_corr).all():
    plt.text(0.5, 0.5, "No Pearson Correlation", ha="center", va="center")
else:
    plt.scatter(cos_sim, pearson_corr, alpha=0.6)
    plt.xlabel("Cosine Similarity")
    plt.ylabel("Pearson Correlation")
    plt.title("Feature-wise Similarity Between Model A and B")
    plt.grid(True)
plt.savefig(f"{SAVE_DIR}/feature_scatter.png")
plt.close()
