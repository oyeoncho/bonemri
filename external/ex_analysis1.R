######## analysis (Image only 전용 버전)

# ---- 설정: file01..file30가 들어있는 상위 폴더 경로로 바꾸세요 ----
base_dir <- "image_data/external/TCGA_CCTH/fu_new/n4_30_30/test1_1"

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
})

# -------------------------
# 1) 파일 경로
# -------------------------
ids <- sprintf("%02d", 1:30)

files_cindex <- file.path(base_dir, paste0("file", ids),
                          paste0("external_cindex_ALL_runs_from_file", ids, ".csv"))
files_auc    <- file.path(base_dir, paste0("file", ids),
                          paste0("external_auc_ALL_runs_from_file",    ids, ".csv"))

exist_cidx <- file.exists(files_cindex)
exist_auc  <- file.exists(files_auc)

if (!any(exist_cidx)) stop("C-index CSV가 없습니다. base_dir를 확인하세요: ", base_dir)
if (!any(exist_auc))  message("⚠️ AUC CSV가 없습니다. (AUC는 NA로 처리)")

ids_cidx <- ids[exist_cidx]
ids_auc  <- ids[exist_auc]

# -------------------------
# 2) 로더
# -------------------------
read_cidx_one <- function(path, id2) {
  read_csv(path, show_col_types = FALSE) %>%
    rename(
      runfile     = `RunFile`,
      feature_set = `Feature Set`,
      time_raw    = `Time (Months)`,
      cindex      = `C-index (External)`
    ) %>%
    mutate(
      time    = as.character(time_raw),
      file_id = paste0("file", id2)
    ) %>%
    select(file_id, runfile, feature_set, time, cindex)
}

read_auc_one <- function(path, id2) {
  read_csv(path, show_col_types = FALSE) %>%
    rename(
      runfile     = `RunFile`,
      feature_set = `Feature Set`,
      time_raw    = `Time (Months)`,
      auc_raw     = `AUC (External)`
    ) %>%
    mutate(
      time    = as.character(time_raw),
      auc     = suppressWarnings(as.numeric(auc_raw)),
      file_id = paste0("file", id2)
    ) %>%
    select(file_id, runfile, feature_set, time, auc)
}

# -------------------------
# 3) 읽기 + 조인
# -------------------------
df_cidx <- map2_dfr(files_cindex[exist_cidx], ids_cidx, read_cidx_one)

if (any(exist_auc)) {
  df_auc <- map2_dfr(files_auc[exist_auc], ids_auc, read_auc_one)
  all_df <- df_cidx %>%
    left_join(df_auc, by = c("file_id", "runfile", "feature_set", "time"))
  message("✅ Joined C-index + AUC by keys: file_id, runfile, feature_set, time")
} else {
  all_df <- df_cidx %>% mutate(auc = NA_real_)
}

# -------------------------
# 4) Overall만 추출
# -------------------------
overall_df_all <- all_df %>%
  filter(time == "Overall") %>%
  select(file_id, runfile, feature_set, time, cindex, auc)

# -------------------------
# 5) Image only만 유지 + 전체 저장
# -------------------------
img_only_df <- overall_df_all %>%
  filter(feature_set == "Image only")

if (nrow(img_only_df) == 0) {
  stop("`Image only` 레코드가 없습니다. 원본 CSV의 Feature Set 값을 확인하세요.")
}

summary_all_path <- file.path(base_dir, "external_overall_summary_image_only.csv")
write_csv(img_only_df %>% arrange(desc(cindex)), summary_all_path)
cat("✅ Overall 요약(Image only):", summary_all_path, "\n")

# -------------------------
# 6) file_id별 요약 (Image only): mean ± 95% CI (C-index, AUC)
# -------------------------
metrics_image_only <- img_only_df %>%
  group_by(file_id) %>%
  summarise(
    # C-index
    n_cidx         = sum(!is.na(cindex)),
    mean_cindex    = mean(cindex, na.rm = TRUE),
    sd_cindex      = sd(cindex, na.rm = TRUE),
    se_cindex      = sd_cindex / sqrt(pmax(n_cidx, 1)),
    tcrit_cindex   = qt(0.975, df = pmax(n_cidx - 1, 1)),
    lower95_cindex = if_else(n_cidx > 1, mean_cindex - tcrit_cindex * se_cindex, NA_real_),
    upper95_cindex = if_else(n_cidx > 1, mean_cindex + tcrit_cindex * se_cindex, NA_real_),
    
    # AUC
    n_auc         = sum(!is.na(auc)),
    mean_auc      = mean(auc, na.rm = TRUE),
    sd_auc        = sd(auc, na.rm = TRUE),
    se_auc        = sd_auc / sqrt(pmax(n_auc, 1)),
    tcrit_auc     = qt(0.975, df = pmax(n_auc - 1, 1)),
    lower95_auc   = if_else(n_auc > 1, mean_auc - tcrit_auc * se_auc, NA_real_),
    upper95_auc   = if_else(n_auc > 1, mean_auc + tcrit_auc * se_auc, NA_real_),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_cindex))

write_csv(metrics_image_only,
          file.path(base_dir, "external_overall_metrics_image_only.csv"))
cat("✅ saved CSV:\n  - external_overall_metrics_image_only.csv\n", sep = "")

# -------------------------
# 7) 그림: Image only (그대로 유지)
# -------------------------
img_df <- img_only_df %>%
  mutate(
    run_num = as.integer(gsub("^run", "", runfile)),
    file_id = factor(file_id),
    runfile = factor(runfile)
  )

file_order_img <- img_df %>%
  group_by(file_id) %>%
  summarise(mean_ci = mean(cindex, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_ci)) %>% pull(file_id)

run_order_img <- img_df %>%
  arrange(run_num) %>% pull(runfile) %>% unique()

img_df <- img_df %>%
  mutate(
    file_id    = factor(file_id, levels = file_order_img),
    runfile    = factor(runfile, levels = run_order_img),
    file_label = paste0("Run", as.integer(gsub("file", "", file_id))),
    file_label = factor(file_label, levels = paste0("Run", as.integer(gsub("file", "", levels(file_id)))))
  )

p_file <- ggplot(img_df, aes(x = file_label, y = cindex)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.4, size = 1.7) +
  geom_hline(yintercept = 0.6, linetype = 2, color = "red", lwd = 2) +
  labs(
    title = "External validation on TCGA-CESC + CCTH (N = 54)\nC-index distribution across 30 runs (n4; Image only)",
    x = "Run ID", y = "External C-index"
  ) +
  scale_y_continuous(breaks = seq(0.3, 0.8, 0.1), limits = c(0.3, 0.8)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15),
    axis.text.y  = element_text(size = 15),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    plot.title   = element_text(size = 20, face = "bold", hjust = 0.5)
  )


# ---- Top 5 Run 기준 mean CV 계산 ----
top5_runs <- img_df %>%
  group_by(file_label) %>%
  summarise(mean_ci = mean(cindex, na.rm = TRUE),
            sd_ci   = sd(cindex, na.rm = TRUE),
            cv_ci   = if_else(is.finite(mean_ci) & mean_ci > 0, sd_ci / mean_ci, NA_real_),
            .groups = "drop") %>%
  arrange(desc(mean_ci)) %>%
  slice_head(n = 5)

top5_cv <- mean(top5_runs$cv_ci, na.rm = TRUE)

# 라벨 위치: x 중앙, y 맨 위쪽
label_x <- length(levels(img_df$file_label)) / 2
label_y <- 0.8

p_file <- p_file +
  annotate("text",
           x = label_x, y = label_y,
           label = paste0("Top 5 Runs Mean Coefficient of Variation = ",
                          scales::percent(top5_cv, accuracy = 0.1)),
           color = "blue", fontface = "bold", size = 7)

print(p_file)
