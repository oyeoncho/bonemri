import os
import torch
import numpy as np
import matplotlib.pyplot as plt
import torch.nn.functional as F
import timm
from torch import nn
import seaborn as sns
import gc

# ----------------------------
# 설정
# -----------------------------
CACHE_DIR = "./slice_cache1_224"  # 마스킹된 슬라이스 텐서들이 저장된 디렉토리
MODEL_PATH = "./models/beit0/best_model_beit_run0.pt"  # BEiT 모델 체크포인트
SAVE_DIR = "./models/beit0/visualization/group_spatial_diff_all_beit0_n7"  # 결과 저장 폴더
os.makedirs(SAVE_DIR, exist_ok=True)

PATCH_SIZE = 16
IMAGE_RESOLUTION = (224, 224)
PATCH_GRID = IMAGE_RESOLUTION[0] // PATCH_SIZE
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

important_features = [436, 519]

# 비교 그룹 정의
group_a = ["VP005_pre_0000", "VP007_pre_0000", "VP013_pre_0000", "VP032_pre_0000", "VP035_pre_0000",
           "VP046_pre_0000", "VP047_pre_0000", "VP052_pre_0000", "VP074_pre_0000", "VP079_pre_0000"]
group_b = ["VP025_pre_0000", "VP034_pre_0000", "VP067_pre_0000", "VP078_pre_0000", "VP147_pre_0000",
           "VP174_pre_0000", "VP197_pre_0000", "VP201_pre_0000", "VP203_pre_0000", "VP228_pre_0000"]

# -----------------------------
# 모델 정의
# -----------------------------
class BEiTFeatureExtractor(nn.Module):
    def __init__(self):
        super().__init__()
        self.model = timm.create_model("beit_base_patch16_224", pretrained=False, num_classes=0)
    def forward(self, x):
        return self.model.forward_features(x)

model = BEiTFeatureExtractor().to(DEVICE)
weights = torch.load(MODEL_PATH, map_location=DEVICE)
model.load_state_dict({k.replace("backbone.", "model."): v 
                       for k, v in weights["state_dict"].items() 
                       if k.startswith("backbone.")}, strict=False)
model.eval()

# -----------------------------
# 환자별 패치 맵 평균 계산
# -----------------------------
def compute_patch_map_mean(pid, feat_idx, max_slices=5):
    path = os.path.join(CACHE_DIR, f"{pid}.pt")
    if not os.path.exists(path):
        print(f"❌ 캐시 없음: {pid}")
        return None
    slices = torch.load(path)
    slices = slices[:max_slices]
    maps = []
    for img_tensor in slices:
        input_tensor = img_tensor.unsqueeze(0).to(DEVICE)
        with torch.no_grad():
            features = model(input_tensor)
        patch_feats = features[:, 1:, :].cpu().numpy()  # [1, 196, 768]
        patch_map = patch_feats[0, :, feat_idx].reshape(PATCH_GRID, PATCH_GRID)
        patch_map = (patch_map - patch_map.min()) / (np.ptp(patch_map) + 1e-6)
        maps.append(patch_map)
        del input_tensor, features
        gc.collect()
        torch.cuda.empty_cache()
    return np.mean(maps, axis=0)

# -----------------------------
# 메인 루프
# -----------------------------
overall_diff_maps = []

for feat_idx in important_features:
    print(f"▶ Feature {feat_idx} 처리 중...")

    group_a_maps = [compute_patch_map_mean(pid, feat_idx) for pid in group_a]
    group_b_maps = [compute_patch_map_mean(pid, feat_idx) for pid in group_b]

    group_a_maps = [m for m in group_a_maps if m is not None]
    group_b_maps = [m for m in group_b_maps if m is not None]

    if len(group_a_maps) == 0 or len(group_b_maps) == 0:
        print(f"⚠️ 그룹 데이터 부족 → feature {feat_idx} 스킵")
        continue

    group_a_mean = np.mean(group_a_maps, axis=0)
    group_b_mean = np.mean(group_b_maps, axis=0)
    diff_map = group_a_mean - group_b_mean
    overall_diff_maps.append(diff_map)

    # 시각화 저장
    plt.figure(figsize=(6, 5))
    sns.heatmap(diff_map, cmap="bwr", center=0, cbar=True, square=True)
    plt.title(f"Spatial Activation Difference (Feature {feat_idx})")
    plt.tight_layout()
    save_path = os.path.join(SAVE_DIR, f"spatial_diff_feat{feat_idx}.png")
    plt.savefig(save_path)
    plt.close()
    print(f"✅ 저장 완료: {save_path}")

# -----------------------------
# 전체 평균 차이 맵 저장
# -----------------------------
if overall_diff_maps:
    mean_diff_map = np.mean(overall_diff_maps, axis=0)
    plt.figure(figsize=(6, 5))
    sns.heatmap(mean_diff_map, cmap="bwr", center=0, cbar=True, square=True)
    plt.title("Average Spatial Activation Difference (All Features)")
    plt.tight_layout()
    plt.savefig(os.path.join(SAVE_DIR, "spatial_diff_all_features.png"))
    plt.close()
    print("✅ 전체 평균 차이 맵 저장 완료")
else:
    print("❌ 유효한 feature 차이 맵이 없음")
