import os
import torch
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------
# 설정
# -----------------------------
CACHE_DIR0 = "./slice_cache0_224"
CACHE_DIR1 = "./slice_cache1_224"
SAVE_DIR = "./visualization"

SAVE_DIR0 = os.path.join(SAVE_DIR, "a0")
SAVE_DIR1 = os.path.join(SAVE_DIR, "a1")
os.makedirs(SAVE_DIR0, exist_ok=True)
os.makedirs(SAVE_DIR1, exist_ok=True)

patient_id = "VP005_pre_0000"
num_slices = 20

# ----------------------------
# 시각화 함수
# -----------------------------
def visualize_and_save(cache_dir0, cache_dir1, patient_id, num_slices, save_dir0, save_dir1, comparison_dir):
    path0 = os.path.join(cache_dir0, f"{patient_id}.pt")
    path1 = os.path.join(cache_dir1, f"{patient_id}.pt")

    tensor0 = torch.load(path0)[:num_slices]  # (N, 3, 224, 224)
    tensor1 = torch.load(path1)[:num_slices]

    fig, axs = plt.subplots(num_slices, 2, figsize=(6, num_slices * 1.5))

    for i in range(num_slices):
        img0 = tensor0[i][0].numpy()
        img1 = tensor1[i][0].numpy()

        # 정규화
        img0 = (img0 - img0.min()) / (img0.max() - img0.min() + 1e-6)
        img1 = (img1 - img1.min()) / (img1.max() - img1.min() + 1e-6)

        # PNG로 각각 저장
        plt.imsave(os.path.join(save_dir0, f"{patient_id}_slice{i:02d}.png"), img0, cmap='gray', vmin=0, vmax=1)
        plt.imsave(os.path.join(save_dir1, f"{patient_id}_slice{i:02d}.png"), img1, cmap='gray', vmin=0, vmax=1)

        # 비교 이미지 시각화용
        axs[i, 0].imshow(img0, cmap='gray', vmin=0, vmax=1)
        axs[i, 0].axis('off')
        axs[i, 0].set_title(f"Cache0 Slice {i:02d}", fontsize=8)

        axs[i, 1].imshow(img1, cmap='gray', vmin=0, vmax=1)
        axs[i, 1].axis('off')
        axs[i, 1].set_title(f"Cache1 Slice {i:02d}", fontsize=8)

    # 전체 비교 이미지 저장
    plt.tight_layout()
    comparison_path = os.path.join(comparison_dir, f"{patient_id}_comparison.png")
    plt.savefig(comparison_path, dpi=300)
    plt.close()
    print(f"✅ 저장 완료: 슬라이스 PNGs → {save_dir0}, {save_dir1}")
    print(f"✅ 저장 완료: 비교 이미지 → {comparison_path}")

# ----------------------------
# 실행
# -----------------------------
visualize_and_save(CACHE_DIR0, CACHE_DIR1, patient_id, num_slices, SAVE_DIR0, SAVE_DIR1, SAVE_DIR)
