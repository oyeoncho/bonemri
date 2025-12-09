import os
import torch
import nibabel as nib
import numpy as np
from PIL import Image
from torchvision import transforms
from tqdm import tqdm
import pandas as pd

# -----------------------------
# 설정
# -----------------------------
IMAGE_DIR = "./image_data/output_196/images"
MASK_DIR = "./image_data/output_196/masks"
SAVE_DIR = "./slice_cache_224"
MAX_SLICES = 20
MODEL_RESOLUTION = (224, 224)  # ✅ 해상도 변경
VALID_AREA_THRESHOLD = 50

os.makedirs(SAVE_DIR, exist_ok=True)

# -----------------------------
# 함수 정의
# -----------------------------
def slice_to_tensor(slice_2d):
    normed = (slice_2d - np.min(slice_2d)) / (np.ptp(slice_2d) + 1e-6)
    img = Image.fromarray(np.uint8(normed * 255)).convert("RGB")
    padded = transforms.Pad((0, 0, max(0, MODEL_RESOLUTION[1] - img.width), max(0, MODEL_RESOLUTION[0] - img.height)))(img)
    return transforms.Compose([
        transforms.Resize(MODEL_RESOLUTION),
        transforms.ToTensor(),
        transforms.Normalize([0.485]*3, [0.229]*3)
    ])(padded)

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

    tensors = torch.stack([slice_to_tensor(s) for s in slices])  # (20, 3, 192, 192)
    torch.save(tensors, os.path.join(SAVE_DIR, f"{pid}.pt"))
    return pid, valid_count

# -----------------------------
# 실행
# -----------------------------
all_pids = [f.replace(".nii.gz", "") for f in os.listdir(IMAGE_DIR) if f.endswith(".nii.gz")]
log_path = os.path.join(SAVE_DIR, "slice_cache_log.csv")

records = []
for pid in tqdm(all_pids, desc="🔄 Caching slices (192x192)"):
    result = process_patient(pid)
    if result:
        records.append(result)

# -----------------------------
# 결과 저장
# -----------------------------
df = pd.DataFrame(records, columns=["PatientID", "ValidSlices"])
df.to_csv(log_path, index=False)
print(f"✅ 총 {len(df)}명 캐싱 완료 (192x192), 로그 저장됨: {log_path}")
