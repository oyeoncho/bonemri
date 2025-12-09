# beit_base_patch16_224

import os
import gc
import random
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torchvision import transforms
from PIL import Image
import nibabel as nib
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from lifelines.utils import concordance_index
from torch.utils.data import Dataset, DataLoader
import matplotlib.pyplot as plt
import timm

# -----------------------------
# 설정
# -----------------------------
IMAGE_DIR = "./image_data/output_224/images"
MASK_DIR = "./image_data/output_224/masks"
CSV_PATH = "./cx.csv"
MODEL_PATH = "./pretrained_models/beit/best_model_overall.pt"
SAVE_DIR = "./models/beit0"
CACHE_DIR = "./slice_cache0_224"

NUM_RUNS = 30
BATCH_SIZE = 8
NUM_WORKERS = 16
TIME_POINTS = [12, 24, 36, 48, 60, 72]
MODEL_RESOLUTION = (224, 224)
MAX_SLICES = 20
MAX_EPOCHS = 500
USE_AUTOMIXED_PRECISION = True

os.makedirs(SAVE_DIR, exist_ok=True)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
torch.backends.cuda.matmul.allow_tf32 = True
torch.backends.cudnn.allow_tf32 = True

def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


class SurvivalSliceDataset(Dataset):
    def __init__(self, df):
        self.df = df

    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):
        row = self.df.iloc[idx]
        pid, time, event = row['PatientID'], row['time'], row['event']

        cache_path = os.path.join(CACHE_DIR, f"{pid}.pt")
        if not os.path.exists(cache_path):
            raise FileNotFoundError(f"❌ 캐시 파일이 존재하지 않습니다: {cache_path}")

        imgs = torch.load(cache_path)  # (N, C, H, W)
        return imgs, torch.tensor(time, dtype=torch.float32), torch.tensor(event, dtype=torch.float32)

class BEiTMixtureSurvival(nn.Module):
    def __init__(self, num_components=2):
        super().__init__()
        self.backbone = timm.create_model("beit_base_patch16_224", pretrained=False, num_classes=0)
        self.head = nn.Sequential(
            nn.Linear(self.backbone.num_features, 256), nn.ReLU(),
            nn.Linear(256, 128), nn.ReLU()
        )
        self.pi_layer = nn.Linear(128, num_components)
        self.lam_layer = nn.Linear(128, num_components)
        self.alpha_layer = nn.Linear(128, num_components)

    def forward(self, x):
        B, N, C, H, W = x.shape
        x = x.view(-1, C, H, W)
        features = self.backbone(x)
        h = self.head(features)
        pi = F.softmax(self.pi_layer(h), dim=1)
        lam = F.softplus(self.lam_layer(h)) + 1e-3
        alpha = F.softplus(self.alpha_layer(h)) + 1e-3
        pi = pi.view(B, N, -1).mean(1)
        lam = lam.view(B, N, -1).mean(1)
        alpha = alpha.view(B, N, -1).mean(1)
        return pi, lam, alpha


def mixture_stretched_nll(t, e, pi, lam, a, eps=1e-8):
    t = t.view(-1, 1)
    t_a = torch.pow(t + eps, a)
    S_k = torch.exp(-lam * t_a)
    f_k = lam * a * torch.pow(t + eps, a - 1) * S_k
    f = torch.sum(pi * f_k, dim=1) + eps
    S = torch.sum(pi * S_k, dim=1) + eps
    loglik = e * torch.log(f) + (1 - e) * torch.log(S)
    return -loglik.mean()

def predict_survival(model, x, times):
    model.eval()
    pi, lam, a = x[:, 0, :], x[:, 1, :], x[:, 2, :]
    surv = []
    for t in times:
        t_tensor = torch.tensor([t], dtype=torch.float32, device=pi.device)
        t_a = torch.pow(t_tensor + 1e-8, a)
        S_k = torch.exp(-lam * t_a)
        S = torch.sum(pi * S_k, dim=1)
        surv.append(S.cpu().numpy())
    return np.vstack(surv)

def safe_concordance_index(times, risks, events):
    times, risks, events = np.asarray(times), np.asarray(risks), np.asarray(events)
    mask = ~(np.isnan(times) | np.isnan(risks) | np.isnan(events))
    if np.sum(mask) < 2 or np.std(risks[mask]) < 1e-6:
        return np.nan
    return concordance_index(times[mask], risks[mask], events[mask])

def save_raw_metrics(cidx_df, auc_df, save_dir):
    cidx_path = os.path.join(save_dir, "raw_cindex_per_time.csv")
    auc_path = os.path.join(save_dir, "raw_auc_per_time.csv")
    cidx_df.to_csv(cidx_path, index=False)
    auc_df.to_csv(auc_path, index=False)
    print(f"💾 Raw C-index saved to {cidx_path}")
    print(f"💾 Raw AUC saved to {auc_path}")

def plot_mean_metrics(cidx_df, auc_df, time_points, save_dir):
    plt.figure(figsize=(10, 5))
    for t in time_points:
        plt.plot(cidx_df['run'], cidx_df[t], label=f"C-index@{t}", alpha=0.6)
    plt.plot(cidx_df['run'], cidx_df['mean_cidx'], label="Mean C-index", linewidth=2, color='black')
    plt.title("Run-wise C-index over Time Points")
    plt.xlabel("Run")
    plt.ylabel("C-index")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, "cindex_runs.png"))
    plt.close()

    plt.figure(figsize=(10, 5))
    for t in time_points:
        plt.plot(auc_df['run'], auc_df[t], label=f"AUC@{t}", alpha=0.6)
    plt.plot(auc_df['run'], auc_df['mean_auc'], label="Mean AUC", linewidth=2, color='black')
    plt.title("Run-wise AUC over Time Points")
    plt.xlabel("Run")
    plt.ylabel("AUC")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, "auc_runs.png"))
    plt.close()
    print("📊 평균 C-index와 AUC 시각화 저장 완료!")

def append_raw_metrics(raw_list, time_points, save_dir, filename):
    df = pd.DataFrame(raw_list)
    if not df.empty:
        if 'mean_cidx' not in df.columns and 'mean_auc' not in df.columns:
            if 'mean_cidx' in filename:
                df['mean_cidx'] = df[[t for t in time_points]].mean(axis=1)
            elif 'mean_auc' in filename:
                df['mean_auc'] = df[[t for t in time_points]].mean(axis=1)
        df.to_csv(os.path.join(save_dir, filename), index=False)
        print(f"💾 Partial raw metrics saved to {filename} (Run {len(raw_list)})")
        
def evaluate(model, loader):
    pi_all, lam_all, alpha_all, t_all, e_all = [], [], [], [], []
    with torch.no_grad():
        for imgs, t, e in loader:
            imgs = imgs[:, :MAX_SLICES].to(device)
            if imgs.shape[2] == 1:
                imgs = imgs.repeat(1, 1, 3, 1, 1)
            B, N, C, H, W = imgs.shape
            x = imgs.view(-1, C, H, W)
            features = model.backbone(x)
            h = model.head(features)
            pi = F.softmax(model.pi_layer(h), dim=1).view(B, N, -1).mean(1)
            lam = F.softplus(model.lam_layer(h)).view(B, N, -1).mean(1) + 1e-3
            alpha = F.softplus(model.alpha_layer(h)).view(B, N, -1).mean(1) + 1e-3
            pi_all.append(pi)
            lam_all.append(lam)
            alpha_all.append(alpha)
            t_all.append(t)
            e_all.append(e)
    params = torch.stack([torch.cat(pi_all), torch.cat(lam_all), torch.cat(alpha_all)], dim=1)
    t_arr = torch.cat(t_all).cpu().numpy()
    e_arr = torch.cat(e_all).cpu().numpy()
    surv = predict_survival(model, params, TIME_POINTS)
    cidx = [safe_concordance_index(t_arr, 1 - surv[i], e_arr) for i in range(len(TIME_POINTS))]
    aucs = {}
    for i, t in enumerate(TIME_POINTS):
        true = ((e_arr == 1) & (t_arr <= t)).astype(int)
        pred = 1 - surv[i, :]
        try:
            aucs[t] = roc_auc_score(true, pred)
        except:
            aucs[t] = np.nan
    return cidx, aucs

# ===============================
# 학습 루프
# ===============================
df_all = pd.read_csv(CSV_PATH)
df_all = df_all.dropna(subset=["fu_date", "survival"])
df_all = df_all[df_all['PatientID'].str.contains("pre")]
df_all = df_all[df_all['PatientID'].apply(
    lambda x: os.path.exists(os.path.join(IMAGE_DIR, f"{x}.nii.gz")) and os.path.exists(os.path.join(MASK_DIR, f"{x}.nii.gz"))
)]
df_all = df_all.rename(columns={"fu_date": "time", "survival": "event"})

raw_cidx_list = []
raw_auc_list = []
raw_train_cidx_list = []
raw_train_auc_list = []
best_val_score = -np.inf
best_overall_model_state = None

for run in range(NUM_RUNS):
    print(f"\n🚀 Run {run+1}/{NUM_RUNS}")
    set_seed(run)
    train_df, val_df = train_test_split(df_all, test_size=0.3, random_state=run)

    train_loader = DataLoader(SurvivalSliceDataset(train_df), batch_size=BATCH_SIZE, shuffle=True, num_workers=NUM_WORKERS)
    val_loader = DataLoader(SurvivalSliceDataset(val_df), batch_size=BATCH_SIZE, shuffle=False, num_workers=NUM_WORKERS)

    model = BEiTMixtureSurvival().to(device)
    weights = torch.load(MODEL_PATH, map_location=device)
    if "state_dict" in weights:
        model.load_state_dict(weights["state_dict"], strict=False)
    else:
        model.load_state_dict(weights, strict=False)

    for name, param in model.backbone.named_parameters():
        if any(f"blocks.{i}" in name for i in range(4, 12)):
            param.requires_grad = True
        else:
            param.requires_grad = False

    for p in list(model.head.parameters()) + list(model.pi_layer.parameters()) + \
             list(model.lam_layer.parameters()) + list(model.alpha_layer.parameters()):
        p.requires_grad = True

    optimizer = torch.optim.Adam(filter(lambda p: p.requires_grad, model.parameters()), lr=1e-4)
    scaler = torch.cuda.amp.GradScaler()
    best_loss, patience, counter = float('inf'), 10, 0

    for epoch in range(MAX_EPOCHS):
        model.train()
        epoch_loss = 0
        for imgs, t, e in train_loader:
            imgs = imgs.to(device)
            if imgs.shape[2] == 1: 
                imgs = imgs.repeat(1, 1, 3, 1, 1)
            with torch.cuda.amp.autocast(enabled=USE_AUTOMIXED_PRECISION):
                pi, lam, alpha = model(imgs)
                loss = mixture_stretched_nll(t.to(device), e.to(device), pi, lam, alpha)
            optimizer.zero_grad()
            scaler.scale(loss).backward()
            scaler.step(optimizer)
            scaler.update()
            epoch_loss += loss.item()
        print(f"Epoch {epoch+1}: Loss={epoch_loss:.4f}")
        if epoch_loss < best_loss - 1e-4:
            best_loss = epoch_loss
            best_model_state = model.state_dict()
            counter = 0
            print("🔁 개선됨. Best loss 갱신.")
            
            torch.save({
                "state_dict": best_model_state,
                "run": run,
                "val_score": best_val_score
            }, os.path.join(SAVE_DIR, f"best_model_beit_run{run}.pt"))

        else:
            counter += 1
            if counter >= patience:
                print("Early stopping.")
                break

        torch.cuda.empty_cache()
        torch.cuda.ipc_collect()

    model.load_state_dict(best_model_state)
    model.eval()

    

    cidx_val, auc_val = evaluate(model, val_loader)
    raw_cidx_list.append({"run": run, "mean_cidx": np.nanmean(cidx_val), **{t: v for t, v in zip(TIME_POINTS, cidx_val)}})
    raw_auc_list.append({"run": run, "mean_auc": np.nanmean(list(auc_val.values())), **auc_val})

    cidx_train, auc_train = evaluate(model, train_loader)
    raw_train_cidx_list.append({"run": run, "mean_cidx": np.nanmean(cidx_train), **{t: v for t, v in zip(TIME_POINTS, cidx_train)}})
    raw_train_auc_list.append({"run": run, "mean_auc": np.nanmean(list(auc_train.values())), **auc_train})

 # Run loop 내부 끝부분 evaluate 이후:

    append_raw_metrics(raw_cidx_list, TIME_POINTS, SAVE_DIR, "raw_cindex_per_time_val_partial.csv")
    append_raw_metrics(raw_auc_list, TIME_POINTS, SAVE_DIR, "raw_auc_per_time_val_partial.csv")
    append_raw_metrics(raw_train_cidx_list, TIME_POINTS, SAVE_DIR, "raw_cindex_per_time_train_partial.csv")
    append_raw_metrics(raw_train_auc_list, TIME_POINTS, SAVE_DIR, "raw_auc_per_time_train_partial.csv")


    if np.nanmean(cidx_val) > best_val_score:
        best_val_score = np.nanmean(cidx_val)
        best_overall_model_state = best_model_state

    del train_loader, val_loader
    import gc
    gc.collect()


# 전체 베스트 모델 저장
if best_overall_model_state is not None:
    torch.save(best_overall_model_state, os.path.join(SAVE_DIR, "best_model_beit_overall.pt"))
    print("✅ 전체 최고 모델 저장 완료")

# 결과 저장 및 시각화
cidx_df = pd.DataFrame(raw_cidx_list)
cidx_df['mean_cidx'] = cidx_df[[t for t in TIME_POINTS]].mean(axis=1)
auc_df = pd.DataFrame(raw_auc_list)
auc_df['mean_auc'] = auc_df[[t for t in TIME_POINTS]].mean(axis=1)
train_cidx_df = pd.DataFrame(raw_train_cidx_list)
train_cidx_df['mean_cidx'] = train_cidx_df[[t for t in TIME_POINTS]].mean(axis=1)
train_auc_df = pd.DataFrame(raw_train_auc_list)
train_auc_df['mean_auc'] = train_auc_df[[t for t in TIME_POINTS]].mean(axis=1)

save_raw_metrics(train_cidx_df, train_auc_df, SAVE_DIR)
save_raw_metrics(cidx_df, auc_df, SAVE_DIR)
plot_mean_metrics(train_cidx_df, train_auc_df, TIME_POINTS, SAVE_DIR)
plot_mean_metrics(cidx_df, auc_df, TIME_POINTS, SAVE_DIR)

print("🎉 전체 학습, 검증, 저장, AUC/C-index 시각화 완료!")