import os
import torch
import nibabel as nib
import numpy as np
from PIL import Image
from torchvision import transforms
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
import random

# -----------------------------
# 설정
# -----------------------------
IMAGE_DIR = "./image_data/output/images"
MASK_DIR = "./image_data/output/masks"
SAVE_DIR = "./slice_cache_224_resize"
MAX_SLICES = 20
MODEL_RESOLUTION = (224, 224)
VALID_AREA_THRESHOLD = 50

os.makedirs(SAVE_DIR, exist_ok=True)

# -----------------------------
# 함수 정의
# -----------------------------
def slice_to_tensor(slice_2d):
    normed = (slice_2d - np.min(slice_2d)) / (np.ptp(slice_2d) + 1e-6)
    img = Image.fromarray(np.uint8(normed * 255)).convert("RGB")
    resized = img.resize(MODEL_RESOLUTION, Image.BILINEAR)
    return transforms.Compose([
        transforms.ToTensor(),
        transforms.Normalize([0.485]*3, [0.229]*3)
    ])(resized)

def extract_masked_slices(vol, mask, valid_area_threshold=VALID_AREA_THRESHOLD):
    mask = (mask > 0).astype(np.float32)
    slices = []
    for i in range(vol.shape[2]):
        m = mask[:, :, i]
        if np.sum(m) >= valid_area_threshold:
            masked = vol[:, :, i] * m
            slices.append(masked)

    num_valid = len(slices)
    if num_valid == 0:
        return [np.zeros_like(vol[:, :, 0])] * MAX_SLICES, num_valid

    if num_valid >= MAX_SLICES:
        center = num_valid // 2
        half = MAX_SLICES // 2
        slices = slices[center - half:center - half + MAX_SLICES]
    else:
        pad = MAX_SLICES - num_valid
        slices += [np.zeros_like(vol[:, :, 0])] * pad

    return slices, num_valid

def process_patient(pid):
    img_path = os.path.join(IMAGE_DIR, f"{pid}.nii.gz")
    mask_path = os.path.join(MASK_DIR, f"{pid}.nii.gz")

    if not os.path.exists(img_path) or not os.path.exists(mask_path):
        return None

    vol = nib.load(img_path).get_fdata().astype(np.float32)
    mask = nib.load(mask_path).get_fdata().astype(np.float32)

    slices, valid_count = extract_masked_slices(vol, mask)
    tensors = torch.stack([slice_to_tensor(s) for s in slices])  # (20, 3, 242, 242)
    torch.save(tensors, os.path.join(SAVE_DIR, f"{pid}.pt"))
    return pid, valid_count

# -----------------------------
# 캐싱 실행
# -----------------------------
all_pids = [f.replace(".nii.gz", "") for f in os.listdir(IMAGE_DIR) if f.endswith(".nii.gz")]
log_path = os.path.join(SAVE_DIR, "slice_cache_log.csv")

records = []
for pid in tqdm(all_pids, desc="🔄 Caching slices (242x242)"):
    result = process_patient(pid)
    if result:
        records.append(result)

df = pd.DataFrame(records, columns=["PatientID", "ValidSlices"])
df.to_csv(log_path, index=False)
print(f"✅ 총 {len(df)}명 캐싱 완료 (242x242), 로그 저장됨: {log_path}")

# -----------------------------
# 캐싱된 .pt 파일 중 3명 시각화
# -----------------------------
NUM_PATIENTS_TO_VISUALIZE = 3
SLICE_INDICES = [4, 9, 14]  # 시각화할 슬라이스 번호

cached_files = [f for f in os.listdir(SAVE_DIR) if f.endswith(".pt")]
sampled_files = random.sample(cached_files, min(NUM_PATIENTS_TO_VISUALIZE, len(cached_files)))

for file in sampled_files:
    tensor = torch.load(os.path.join(SAVE_DIR, file))  # (20, 3, 242, 242)
    pid = file.replace(".pt", "")

    fig, axes = plt.subplots(1, len(SLICE_INDICES), figsize=(12, 4))
    fig.suptitle(f"📊 Patient ID: {pid}", fontsize=14)

    for i, idx in enumerate(SLICE_INDICES):
        slice_tensor = tensor[idx]
        img = slice_tensor * 0.229 + 0.485  # Unnormalize
        img = img.permute(1, 2, 0).numpy()
        img = np.clip(img, 0, 1)

        axes[i].imshow(img)
        axes[i].set_title(f"Slice {idx+1}")
        axes[i].axis("off")

    plt.tight_layout()
    plt.show()
