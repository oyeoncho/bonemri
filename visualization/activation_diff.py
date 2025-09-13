import os
import torch
import numpy as np
import matplotlib.pyplot as plt
import torch.nn.functional as F
import timm
from torch import nn
from tqdm import tqdm
from torchvision.transforms import ToPILImage
import imageio
import gc

# -----------------------------
# 설정
# -----------------------------
CACHE_DIR = "./slice_cache1_224"
MODEL_PATH = "./models/beit0/best_model_beit_run0.pt"
SAVE_DIR = "./models/beit0/visualization/patient_patch_vis_beit0_n5"
os.makedirs(SAVE_DIR, exist_ok=True)

PATCH_SIZE = 16
IMAGE_RESOLUTION = (224, 224)
PATCH_GRID = IMAGE_RESOLUTION[0] // PATCH_SIZE
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
#important_features = [436, 519]
#important_features = [213, 215, 266, 436, 519]
important_features = [2, 107, 109, 137, 173, 209, 213, 215, 223, 266, 327, 374, 391, 436, 499, 519, 554, 715]
#important_features = [2, 10, 15, 25, 55, 107, 109, 123, 137, 143, 152, 163, 167, 169, 173, 181, 194, 209, 210, 213, 215, 220, 223, 240, 255, 266, 289, 309, 327, 328, 352, 361, 368, 374, 391, 402, 407, 420, 436, 439, 451, 498, 499, 507, 514, 519, 554, 560, 565, 576, 577, 583, 605, 617, 633, 653, 656, 657, 666, 710, 715, 747]

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
model.load_state_dict({k.replace("backbone.", "model."): v for k, v in weights["state_dict"].items() if k.startswith("backbone.")}, strict=False)
model.eval()

# -----------------------------
# 유틸: 16의 배수로 padding
# -----------------------------
def pad_to_multiple_of_16(image):
    h, w, c = image.shape
    new_h = (h + 15) // 16 * 16
    new_w = (w + 15) // 16 * 16
    pad_h = new_h - h
    pad_w = new_w - w
    return np.pad(image, ((0, pad_h), (0, pad_w), (0, 0)), mode='constant')

# -----------------------------
# 시각화 함수 (회전 없음, 기본 여백 유지)
# -----------------------------
def visualize_patient(pid, max_slices=20):
    path = os.path.join(CACHE_DIR, f"{pid}.pt")
    if not os.path.exists(path):
        print(f"\u274c 캐시 없음: {pid}")
        return

    slices = torch.load(path)[:max_slices]
    frames = []

    for i, img_tensor in enumerate(tqdm(slices, desc=f"Processing {pid}")):
        input_tensor = img_tensor.unsqueeze(0).to(DEVICE)
        with torch.no_grad():
            features = model(input_tensor)

        patch_feats = features[:, 1:, :].cpu().numpy()
        patch_maps = [patch_feats[0, :, idx].reshape(PATCH_GRID, PATCH_GRID) for idx in important_features]
        avg_map = np.mean(patch_maps, axis=0)
        avg_map = (avg_map - avg_map.min()) / (np.ptp(avg_map) + 1e-6)

        avg_map_up = torch.tensor(avg_map).unsqueeze(0).unsqueeze(0)
        cam_up = F.interpolate(avg_map_up, size=IMAGE_RESOLUTION, mode='bilinear', align_corners=False)[0, 0].numpy()
        cam_up = (cam_up - cam_up.min()) / (np.ptp(cam_up) + 1e-6)

        img_np = img_tensor.numpy()
        if img_np.shape[0] == 3:
            img_np = img_np.transpose(1, 2, 0)
        elif img_np.shape[0] == 1:
            img_np = img_np.squeeze(0)
        img_np = (img_np - img_np.min()) / (np.ptp(img_np) + 1e-6)

        fig, ax = plt.subplots(figsize=(4, 4))
        ax.imshow(img_np, cmap='gray')
        ax.imshow(cam_up, cmap='jet', alpha=0.6)
        ax.axis('off')

        fig.canvas.draw()
        img = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
        img = img.reshape(fig.canvas.get_width_height()[::-1] + (3,))
        img_padded = pad_to_multiple_of_16(img)
        frames.append(img_padded)

        plt.close(fig)
        gc.collect()

    mp4_path = os.path.join(SAVE_DIR, f"{pid}_heatmap.mp4")
    imageio.mimsave(mp4_path, frames, fps=3)
    print(f"\u2705 저장 완료: {mp4_path}")

# -----------------------------
# 실행
# -----------------------------
visualize_patient("VP005_pre_0000")
visualize_patient("VP025_pre_0000")
