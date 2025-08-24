import os
import torch
import nibabel as nib
import numpy as np
import pandas as pd
from tqdm import tqdm
from PIL import Image
from torchvision import transforms

# -----------------------------
# 설정
# -----------------------------
IMAGE_DIR = "./image_data/output_224/images"
MASK_DIR = "./image_data/output_224/masks"
SAVE_DIR = "./slice_cache1_224"  # ✅ 폴더명 변경
MAX_SLICES = 20
VALID_AREA_THRESHOLD = 50
MODEL_RESOLUTION = None  # (H, W) 유지. 필요시 (224, 224) 등 설정 가능

os.makedirs(SAVE_DIR, exist_ok=True)

# -----------------------------
# 함수 정의
# -----------------------------
def slice_to_tensor(slice_2d):
    # 1. Min-max 정규화 (0–255)
    normed = (slice_2d - np.min(slice_2d)) / (np.ptp(slice_2d) + 1e-6)
    img = Image.fromarray(np.uint8(normed * 255)).convert("RGB")  # 3채널 복제

    # 2. Transform 정의
    transform = transforms.Compose([
        transforms.ToTensor(),  # (H, W, 3) → (3, H, W), [0,1]
        transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225])  # ImageNet
    ])

    if MODEL_RESOLUTION:
        img = transforms.Resize(MODEL_RESOLUTION)(img)

    return transform(img)

def extract_bottom_masked_slices(vol, mask, valid_area_threshold=VALID_AREA_THRESHOLD):
    mask = (mask > 0).astype(np.float32)
    valid_slices = []
    for i in range(vol.shape[2]):
        m = mask[:, :, i]
        if np.sum(m) >= valid_area_threshold:
            masked = vol[:, :, i] * m
            valid_slices.append(masked)

    num_valid = len(valid_slices)
    if num_valid == 0:
        return [np.zeros_like(vol[:, :, 0])] * MAX_SLICES, num_valid

    if num_valid >= MAX_SLICES:
        selected = valid_slices[-MAX_SLICES:]  # ✅ 뒤쪽(아래쪽)부터 선택
    else:
        selected = [np.zeros_like(vol[:, :, 0])] * (MAX_SLICES - num_valid) + valid_slices

    return selected, num_valid

def process_patient(pid):
    img_path = os.path.join(IMAGE_DIR, f"{pid}.nii.gz")
    mask_path = os.path.join(MASK_DIR, f"{pid}.nii.gz")

    if not os.path.exists(img_path) or not os.path.exists(mask_path):
        return None

    vol = nib.load(img_path).get_fdata().astype(np.float32)
    mask = nib.load(mask_path).get_fdata().astype(np.float32)

    slices, valid_count = extract_bottom_masked_slices(vol, mask)

    tensors = torch.stack([slice_to_tensor(s) for s in slices])  # (20, 3, H, W)
    torch.save(tensors, os.path.join(SAVE_DIR, f"{pid}.pt"))

    return pid, valid_count

# -----------------------------
# 실행
# -----------------------------
all_pids = [
    f.replace(".nii.gz", "")
    for f in os.listdir(IMAGE_DIR)
    if f.endswith(".nii.gz") and not f.startswith("._")
]
log_path = os.path.join(SAVE_DIR, "bottom_slice_log.csv")

records = []
for pid in tqdm(all_pids, desc="🔄 Caching 20 bottom slices (RGB + norm)"):
    result = process_patient(pid)
    if result:
        records.append(result)

# -----------------------------
# 결과 저장
# -----------------------------
df = pd.DataFrame(records, columns=["PatientID", "ValidSlices"])
df.to_csv(log_path, index=False)
print(f"✅ 총 {len(df)}명 캐싱 완료 (3채널, ImageNet 정규화), 로그 저장됨: {log_path}")
