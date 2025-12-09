import os
import torch
import torch.nn as nn
import pandas as pd
import numpy as np
import timm
from tqdm import tqdm

# 설정
CACHE_DIR = "./slice_cache0_224"
MODEL_DIR = "./models/beit0"
CSV_PATH = "./cx.csv"
FEATURE_DIR = "./features/beit0_features"
NUM_RUNS = 30
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# 디렉토리 준비
os.makedirs(FEATURE_DIR, exist_ok=True)

# 모델 정의 (BEiT 백본 사용)
class BEiTBackbone(nn.Module):
    def __init__(self):
        super().__init__()
        self.backbone = timm.create_model("beit_base_patch16_224", pretrained=False, num_classes=0)

    def forward(self, x):
        return self.backbone(x)  # (B, D)

# DataFrame 로딩
df = pd.read_csv(CSV_PATH)
df = df[df['PatientID'].str.contains('pre')]
df = df[df['PatientID'].apply(lambda x: os.path.exists(os.path.join(CACHE_DIR, f"{x}.pt")))]
df = df.reset_index(drop=True)

print(f"\n👉 [DEBUG] 최종 df PatientID 수: {len(df)}")
if len(df) > 0:
    print(df["PatientID"].head())
else:
    print("⚠️ df 가 비어있음! → CSV 또는 캐시 확인 필요")

# 전체 Run 반복
for run in range(NUM_RUNS):
    print(f"\n🚀 Run {run}/{NUM_RUNS} → Feature 추출 시작...")

    model_path = os.path.join(MODEL_DIR, f"best_model_beit_run{run}.pt")
    if not os.path.exists(model_path):
        print(f"❌ Run {run} 모델 없음 → skip: {model_path}")
        continue

    # 모델 로딩
    model = BEiTBackbone().to(DEVICE)
    state_dict = torch.load(model_path, map_location=DEVICE, weights_only=False)["state_dict"]
    model.backbone.load_state_dict(
        {k.replace("backbone.", ""): v for k, v in state_dict.items() if k.startswith("backbone.")}, strict=False
    )
    model.eval()

    # 피처 추출
    features = []
    with torch.no_grad():
        for _, row in tqdm(df.iterrows(), total=len(df)):
            pid = row["PatientID"]
            cache_path = os.path.join(CACHE_DIR, f"{pid}.pt")
            if not os.path.exists(cache_path):
                print(f"⚠️ {pid} 캐시 없음 → skip")
                continue
            imgs = torch.load(cache_path).to(DEVICE)  # (20, 3, 224, 224)
            if imgs.shape[1] == 1:
                imgs = imgs.repeat(1, 3, 1, 1)
            feats = model(imgs)  # (20, D)
            mean_feat = feats.mean(dim=0).cpu().numpy()
            features.append([pid] + mean_feat.tolist())

    # 저장
    if features:
        feature_dim = len(features[0]) - 1
        columns = ["PatientID"] + [f"feat_{i}" for i in range(feature_dim)]
        df_feat = pd.DataFrame(features, columns=columns)
        save_path = os.path.join(FEATURE_DIR, f"beit_base_patch16_224_backbone_features_run{run}.csv")
        df_feat.to_csv(save_path, index=False)
        print(f"✅ Run {run} Feature 저장 완료 → {save_path}")
    else:
        print(f"⚠️ Run {run} → 추출된 feature 없음 → 저장 skip")

# 병합
print("\n📢 전체 Run Feature 병합 시작...")
dfs = []
for run in range(NUM_RUNS):
    file_path = os.path.join(FEATURE_DIR, f"beit_base_patch16_224_backbone_features_run{run}.csv")
    if os.path.exists(file_path):
        df_run = pd.read_csv(file_path)
        df_run["Run"] = run
        dfs.append(df_run)
    else:
        print(f"⚠️ Run {run} feature 파일 없음 → skip")

if dfs:
    df_merged = pd.concat(dfs, ignore_index=True)
    merged_path = os.path.join(FEATURE_DIR, "beit_base_patch16_224_backbone_features_merged.csv")
    df_merged.to_csv(merged_path, index=False)
    print(f"🎉 전체 Feature 병합 저장 완료 → {merged_path}")
    print("\n📊 Run별 Patient 수:")
    print(df_merged['Run'].value_counts().sort_index())
else:
    print("⚠️ 병합할 Feature 파일이 없음!")
