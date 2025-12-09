# ================================
# 공통 로딩 파트
#  - all_results : 모든 Feature set (Clinical / Image / Image+Clinical)
#  - all_results_img : Image only만 필터링
# ================================

library(dplyr)
library(tidyr)
library(purrr)
library(readr)

# ---------------------------
# 1) 설정
# ---------------------------

# 그룹: n4_30_30 ~ n7_30_30
groups <- paste0("n", 7, "_30_30")

# base 경로 (필요하면 여기만 바꾸면 됨)
base_dir <- "image_data/survival_model/mixture_non_fix/non_nest/beit0/results/generalization/test2/dl0"

# ---------------------------
# 2) 한 그룹(grp) 로딩 함수
# ---------------------------

load_one_group <- function(grp) {
  message("=== 그룹 로딩: ", grp, " ===")
  
  # run 01~30 반복
  map_dfr(1:30, function(i) {
    auc_path <- file.path(
      base_dir, grp,
      sprintf("raw_auc_per_time_run%02d.csv", i)
    )
    cindex_path <- file.path(
      base_dir, grp,
      sprintf("raw_cindex_per_time_run%02d.csv", i)
    )
    
    # 파일 존재 체크 (없으면 그 run 스킵)
    if (!file.exists(auc_path) || !file.exists(cindex_path)) {
      message("  - run ", sprintf("%02d", i), " 스킵 (파일 없음)")
      return(NULL)
    }
    
    # AUC 읽기
    a <- read.csv(auc_path) %>%
      as_tibble() %>%
      rename(
        Feature   = Feature.Set,
        Time      = Time..Months.,
        AUC_train = AUC..Train.,
        AUC_val   = AUC..Val.
      )
    
    # C-index 읽기
    b <- read.csv(cindex_path) %>%
      as_tibble() %>%
      rename(
        Feature       = Feature.Set,
        Time          = Time..Months.,
        C_index_train = C.index..Train.,
        C_index_val   = C.index..Val.
      )
    
    # "Overall" 제외한 time-point만
    a_time <- a %>% filter(Time != "Overall")
    b_time <- b %>% filter(Time != "Overall")
    
    if (nrow(a_time) == 0 || nrow(b_time) == 0) {
      message("  - run ", sprintf("%02d", i), " 스킵 (time-point 행 없음)")
      return(NULL)
    }
    
    # 안전하게 Feature + Time 기준으로 조인
    ab <- inner_join(
      b_time %>% select(Feature, Time, C_index_train, C_index_val),
      a_time %>% select(Feature, Time, AUC_train, AUC_val),
      by = c("Feature", "Time")
    )
    
    if (nrow(ab) == 0) {
      message("  - run ", sprintf("%02d", i), " 스킵 (조인 결과 0행)")
      return(NULL)
    }
    
    # Feature별로 time-point 평균 내기
    combined <- ab %>%
      group_by(Feature) %>%
      summarise(
        mean_AUC_train    = mean(AUC_train,    na.rm = TRUE),
        mean_AUC_val      = mean(AUC_val,      na.rm = TRUE),
        mean_cindex_train = mean(C_index_train,na.rm = TRUE),
        mean_cindex_val   = mean(C_index_val,  na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        diff_AUC    = mean_AUC_train    - mean_AUC_val,
        diff_cindex = mean_cindex_train - mean_cindex_val,
        run   = i,
        group = grp
      )
    
    combined
  })
}

# ---------------------------
# 3) 모든 그룹(n4~n7) 로딩 후 하나로 병합
# ---------------------------

all_results_list <- lapply(groups, load_one_group)
names(all_results_list) <- groups

# 전체 Feature 세트 포함 (Clinical / Image / Image+Clinical)
all_results <- bind_rows(all_results_list) %>%
  mutate(
    group = gsub("_30_30", "", group),  # n4_30_30 → n4
    Feature = as.character(Feature)
  )

# 확인
dplyr::count(all_results, group, Feature)


# ================================
# Distant Metastasis: n7 (Feature set 비교)
#  - 그림: C-index / AUC by Feature set (n7)
#  - 표  : Feature별 mean ± SD + bootstrap 95% CI + p-value
# ================================

library(dplyr)
library(ggplot2)
library(ggpubr)
library(boot)

## 1) Feature factor 순서 고정
feature_levels <- c("Clinical only","Image only","Image + Clinical")
all_results$Feature <- factor(all_results$Feature, levels = feature_levels)

## 2) n7만 필터 (DM endpoint의 n7)
all_results_n7 <- all_results %>%
  filter(group == "n7",
         Feature %in% feature_levels)

# sanity check
dplyr::count(all_results_n7, Feature)

# ----------------------------------------
# (A) 그림 1: C-index by Feature set (n7, DM)
# ----------------------------------------
p_cindex_n7_dm <- ggplot(all_results_n7,
                         aes(x = Feature, y = mean_cindex_val, fill = Feature)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  scale_fill_brewer(palette = "Set2", name = "Feature Set") +
  labs(
    title    = "Validation C-index by Feature Set (DM)",
    y        = "Mean C-index (validation)",
    x        = ""
  ) +
  coord_cartesian(ylim = c(0.4, 1)) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "top",        # 🔹 범례 위쪽
    legend.box      = "horizontal",
    axis.text.x     = element_text(angle = 30, hjust = 1, size = 20),
    axis.text.y     = element_text(size = 25),
    axis.title.y    = element_text(size = 25, face = "bold"),
    strip.text      = element_text(size = 25, face = "bold"),
    plot.title      = element_text(face = "bold", size = 35),
    plot.subtitle   = element_text(size = 29, margin = margin(b = 10), face = "bold"),
    legend.title    = element_text(face = "bold", size = 25),
    legend.text     = element_text(size = 25)
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      c("Image only",    "Clinical only"),
      c("Image only",    "Image + Clinical"),
      c("Clinical only", "Image + Clinical")
    ),
    label = "p.signif",
    size  = 6,
    step.increase = 0.05
  )

print(p_cindex_n7_dm)

# ----------------------------------------
# (B) 그림 2: AUC by Feature set (n7, DM)
# ----------------------------------------
p_auc_n7_dm <- ggplot(all_results_n7,
                      aes(x = Feature, y = mean_AUC_val, fill = Feature)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  scale_fill_brewer(palette = "Set2", name = "Feature Set") +
  labs(
    title    = "Validation AUC by Feature Set (DM)",
    y        = "Mean AUC (validation)",
    x        = ""
  ) +
  coord_cartesian(ylim = c(0.4, 1)) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "top",        # 🔹 범례 위쪽
    legend.box      = "horizontal",
    axis.text.x     = element_text(angle = 30, hjust = 1, size = 20),
    axis.text.y     = element_text(size = 25),
    axis.title.y    = element_text(size = 25, face = "bold"),
    strip.text      = element_text(size = 25, face = "bold"),
    plot.title      = element_text(face = "bold", size = 35),
    plot.subtitle   = element_text(size = 29, margin = margin(b = 10), face = "bold"),
    legend.title    = element_text(face = "bold", size = 25),
    legend.text     = element_text(size = 25)
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      c("Image only",    "Clinical only"),
      c("Image only",    "Image + Clinical"),
      c("Clinical only", "Image + Clinical")
    ),
    label = "p.signif",
    size  = 6,
    step.increase = 0.05
  )

print(p_auc_n7_dm)

# 필요하면 저장
# ggsave("dm_n7_cindex_by_feature.png", p_cindex_n7_dm, width = 8, height = 6, dpi = 300)
# ggsave("dm_n7_auc_by_feature.png",    p_auc_n7_dm,    width = 8, height = 6, dpi = 300)


# ----------------------------------------
# (C) 표: Feature별 mean ± SD + bootstrap 95% CI
# ----------------------------------------

# 1) 요약 (mean / SD)
summary_n7_dm <- all_results_n7 %>%
  group_by(Feature) %>%
  summarise(
    mean_cindex = mean(mean_cindex_val, na.rm = TRUE),
    sd_cindex   = sd(mean_cindex_val,   na.rm = TRUE),
    mean_auc    = mean(mean_AUC_val,    na.rm = TRUE),
    sd_auc      = sd(mean_AUC_val,      na.rm = TRUE),
    .groups = "drop"
  )

# 2) 부트스트랩용 함수
boot_mean <- function(data, indices) {
  d <- data[indices]
  mean(d, na.rm = TRUE)
}

# 3) Feature별 bootstrap CI (C-index / AUC)
boot_ci_n7_dm <- all_results_n7 %>%
  group_by(Feature) %>%
  summarise(
    # C-index bootstrap
    boot_cindex = list(boot(mean_cindex_val, statistic = boot_mean, R = 5000)),
    cindex_ci_lower = {
      bc <- tryCatch(
        boot.ci(boot_cindex[[1]], type = "perc"),
        error = function(e) NULL
      )
      if (is.null(bc) || is.null(bc$percent)) NA_real_ else bc$percent[4]
    },
    cindex_ci_upper = {
      bc <- tryCatch(
        boot.ci(boot_cindex[[1]], type = "perc"),
        error = function(e) NULL
      )
      if (is.null(bc) || is.null(bc$percent)) NA_real_ else bc$percent[5]
    },
    # AUC bootstrap
    boot_auc = list(boot(mean_AUC_val, statistic = boot_mean, R = 5000)),
    auc_ci_lower = {
      bc <- tryCatch(
        boot.ci(boot_auc[[1]], type = "perc"),
        error = function(e) NULL
      )
      if (is.null(bc) || is.null(bc$percent)) NA_real_ else bc$percent[4]
    },
    auc_ci_upper = {
      bc <- tryCatch(
        boot.ci(boot_auc[[1]], type = "perc"),
        error = function(e) NULL
      )
      if (is.null(bc) || is.null(bc$percent)) NA_real_ else bc$percent[5]
    },
    .groups = "drop"
  ) %>%
  select(
    Feature,
    cindex_ci_lower, cindex_ci_upper,
    auc_ci_lower,    auc_ci_upper
  )

# 4) 최종 요약 테이블 (논문용 포맷)
summary_with_ci_n7_dm <- summary_n7_dm %>%
  left_join(boot_ci_n7_dm, by = "Feature") %>%
  mutate(
    Cindex_mean_sd = sprintf("%.3f ± %.3f", mean_cindex, sd_cindex),
    AUC_mean_sd    = sprintf("%.3f ± %.3f", mean_auc,    sd_auc),
    Cindex_CI      = ifelse(
      is.na(cindex_ci_lower),
      "NA",
      sprintf("[%.3f, %.3f]", cindex_ci_lower, cindex_ci_upper)
    ),
    AUC_CI         = ifelse(
      is.na(auc_ci_lower),
      "NA",
      sprintf("[%.3f, %.3f]", auc_ci_lower, auc_ci_upper)
    )
  ) %>%
  select(
    Feature,
    mean_cindex, sd_cindex, Cindex_CI, Cindex_mean_sd,
    mean_auc,    sd_auc,    AUC_CI,    AUC_mean_sd
  )

print(summary_with_ci_n7_dm)

# 필요하면 CSV 저장
write.csv(summary_with_ci_n7_dm,
           "dm_n7_feature_summary_with_bootstrap_CI.csv",
           row.names = FALSE)


# ----------------------------------------
# (D) n7 Feature-set 간 p-value (Wilcoxon, Holm 보정)
# ----------------------------------------

# C-index에 대한 쌍별 Wilcoxon
pw_cindex_n7_dm <- pairwise.wilcox.test(
  all_results_n7$mean_cindex_val,
  all_results_n7$Feature,
  p.adjust.method = "holm"
)

# AUC에 대한 쌍별 Wilcoxon
pw_auc_n7_dm <- pairwise.wilcox.test(
  all_results_n7$mean_AUC_val,
  all_results_n7$Feature,
  p.adjust.method = "holm"
)

# pairwise.wilcox.test 결과에서 관심 쌍만 추출
p_cindex_clin_img    <- pw_cindex_n7_dm$p.value["Image only",       "Clinical only"]
p_cindex_img_imgClin <- pw_cindex_n7_dm$p.value["Image + Clinical", "Image only"]

p_auc_clin_img       <- pw_auc_n7_dm$p.value["Image only",          "Clinical only"]
p_auc_img_imgClin    <- pw_auc_n7_dm$p.value["Image + Clinical",    "Image only"]

p_table_n7_dm <- tibble::tibble(
  Comparison = c("Clinical only vs Image only",
                 "Image only vs Image + Clinical"),
  p_Cindex   = c(p_cindex_clin_img,   p_cindex_img_imgClin),
  p_AUC      = c(p_auc_clin_img,      p_auc_img_imgClin)
) %>%
  mutate(
    p_Cindex_fmt = sprintf("%.3f", p_Cindex),
    p_AUC_fmt    = sprintf("%.3f", p_AUC)
  )

print(p_table_n7_dm)

# 필요하면 CSV로 저장
 write.csv(p_table_n7_dm,
           "dm_n7_pairwise_pvalues.csv",
           row.names = FALSE)

