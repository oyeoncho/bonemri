#!/bin/bash
#SBATCH --job-name=beit_original_ex
#SBATCH --partition=red
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# logs 폴더 없으면 생성
mkdir -p logs

# Conda 환경 활성화
source ~/.bashrc
conda activate my_env

# 디버그 정보 출력
echo "======================================"
echo "Running on host: $(hostname)"
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"
echo "Start time: $(date)"
echo "======================================"
nvidia-smi
echo "======================================"

# run.py 실행
python beit_original_ex.py

echo "======================================"
echo "End time: $(date)"
echo "Run completed!"
echo "======================================"

