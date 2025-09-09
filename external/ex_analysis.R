######## analysis — 전체 스크립트 (정렬·저장 순서 수정판)

# ---- 설정: file01..file30가 들어있는 상위 폴더 경로로 바꾸세요 ----
base_dir <- "image_data/external/TCGA_CCRT/fu_new/n7_30_30/test1_1"

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(tidyr)
  library(scales)   # ← percent() 사용
})

# 옵션: pairwise p-values 계산 여부
compute_pairwise <- TRUE

# -------------------------
# 0) 유틸
# -------------------------
slug_fs <- function(x) {
  x %>% trimws() %>% tolower() %>%
    gsub("[^a-z0-9]+", "_", .) %>% gsub("^_|_$", "", .)
}

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
      file_id = paste0("file", id2),
      feature_set = trimws(feature_set)
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
      file_id = paste0("file", id2),
      feature_set = trimws(feature_set)
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
# 4) Overall만 추출 (모든 Feature Set 포함)
# -------------------------
overall_df_all <- all_df %>%
  filter(time == "Overall") %>%
  mutate(feature_set = trimws(feature_set)) %>%
  select(file_id, runfile, feature_set, time, cindex, auc)

# 전체 저장 (전체 Feature Set 합본)
summary_all_path <- file.path(base_dir, "external_overall_summary_allfeatures.csv")
write_csv(overall_df_all %>% arrange(feature_set, desc(cindex)), summary_all_path)
cat("✅ Overall 요약(전체 Feature Set):", summary_all_path, "\n")

# -------------------------
# 5) file_id × Feature Set: mean ± 95% CI (C-index & AUC)
# -------------------------
sets_use <- c("Clinical only", "Image + Clinical", "Image only")
fs_map <- c("Clinical only" = "clinical",
            "Image + Clinical" = "img_clin",
            "Image only" = "image")

metrics_by_feature <- overall_df_all %>%
  filter(feature_set %in% sets_use) %>%
  group_by(file_id, feature_set) %>%
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
  arrange(file_id, factor(feature_set, levels = sets_use))

# -------------------------
# 6) (선택) pairwise p-values (paired t-test; runfile 매칭)
# -------------------------
if (isTRUE(compute_pairwise)) {
  pair_list  <- list(c("clinical","image"),
                     c("img_clin","image"),
                     c("img_clin","clinical"))
  pair_names <- c("Clinical vs Image",
                  "Image+Clinical vs Image",
                  "Image+Clinical vs Clinical")
  
  pairwise_fun <- function(dat, value_col, metric_label) {
    w <- dat %>%
      select(runfile, fs_short, value = {{ value_col }}) %>%
      pivot_wider(names_from = fs_short, values_from = value)
    
    purrr::imap_dfr(pair_list, function(p, i) {
      dfp <- w %>% select(all_of(p))
      ok  <- stats::complete.cases(dfp)
      x   <- dfp[[p[1]]][ok]; y <- dfp[[p[2]]][ok]
      n   <- length(x)
      if (n >= 2) {
        tt <- t.test(x, y, paired = TRUE)
        tibble::tibble(
          pair_label = pair_names[i],
          metric     = metric_label,
          x_set = p[1], y_set = p[2],
          n_pairs = n,
          mean_x  = mean(x), mean_y = mean(y),
          mean_diff = mean(x - y),
          t_stat  = unname(tt$statistic),
          p_value = tt$p.value
        )
      } else {
        tibble::tibble(
          pair_label = pair_names[i],
          metric     = metric_label,
          x_set = p[1], y_set = p[2],
          n_pairs = n,
          mean_x = NA_real_, mean_y = NA_real_, mean_diff = NA_real_,
          t_stat = NA_real_, p_value = NA_real_
        )
      }
    })
  }
  
  pairwise_pvals <- overall_df_all %>%
    filter(feature_set %in% sets_use) %>%
    mutate(fs_short = recode(feature_set, !!!fs_map)) %>%
    group_by(file_id) %>%
    group_modify(~{
      dat <- .x
      bind_rows(
        pairwise_fun(dat, cindex, "cindex"),
        pairwise_fun(dat, auc,    "auc")
      )
    }) %>%
    ungroup() %>%
    group_by(file_id, metric) %>%
    mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
    ungroup()
}

# -------------------------
# 7) 저장 — 계산 후 저장(순서 중요!)
# -------------------------
# 통합 저장
write_csv(metrics_by_feature,
          file.path(base_dir, "external_overall_metrics_by_feature.csv"))
cat("✅ saved: external_overall_metrics_by_feature.csv\n")

if (isTRUE(compute_pairwise)) {
  write_csv(pairwise_pvals,
            file.path(base_dir, "external_overall_pairwise_pvalues.csv"))
  cat("✅ saved: external_overall_pairwise_pvalues.csv\n")
}

# Feature Set별 분할 저장 (metrics_by_feature)
unique_metrics_fs <- sort(unique(metrics_by_feature$feature_set))
for (fs in unique_metrics_fs) {
  fs_slug <- slug_fs(fs)
  outp <- file.path(base_dir, paste0("external_overall_metrics_by_feature_", fs_slug, ".csv"))
  metrics_by_feature %>%
    filter(feature_set == fs) %>%
    arrange(desc(mean_cindex), file_id) %>%
    write_csv(outp)
  cat("  └─ saved (metrics_by_feature split):", outp, "\n")
}

# Summary도 Feature Set별로 분할 저장
unique_fs <- sort(unique(overall_df_all$feature_set))
for (fs in unique_fs) {
  fs_slug <- slug_fs(fs)
  p <- file.path(base_dir, paste0("external_overall_summary_", fs_slug, ".csv"))
  overall_df_all %>%
    filter(feature_set == fs) %>%
    arrange(desc(cindex)) %>%
    write_csv(p)
  cat("  └─ saved (summary split):", p, "\n")
}

# -------------------------
# 8) 그림: Image only만 유지 (평균 C-index 내림차순 정렬)
# -------------------------
set.seed(1234)  # 지터 재현성 고정

img_df <- overall_df_all %>% filter(feature_set == "Image only")

if (nrow(img_df) > 0) {
  img_df <- img_df %>%
    mutate(
      run_num  = suppressWarnings(as.integer(gsub("^run", "", runfile))),
      file_id  = factor(file_id),
      runfile  = factor(runfile)
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
  
  # 기본 박스플롯
  p_file <- ggplot(img_df, aes(x = file_label, y = cindex)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.4, size = 1.7) +
    geom_hline(yintercept = 0.7, linetype = 2, color = "red", lwd = 2) +
    labs(
      title = "External validation on TCGA-CESC (N=38)\nC-index distribution across 30 runs (n7; Image only)",
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
  
} else {
  message("ℹ️ 'Image only' 레코드가 없어 그림을 건너뜁니다.")
}


