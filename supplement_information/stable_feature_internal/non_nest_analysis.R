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
groups <- paste0("n", 4:7, "_30_30")

# base 경로 (필요하면 여기만 바꾸면 됨)
base_dir <- "image_data/survival_model/mixture_non_fix/non_nest/beit0/results/generalization/test0/dl0"

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

# ---------------------------
# 4) Image-only만 별도 객체로
# ---------------------------

all_results_img <- all_results %>%
  filter(Feature == "Image only")

# 확인
dplyr::count(all_results_img, group)

# ---------------------------
# 2) run 단위 요약 (각 group당 30개 값)
# ---------------------------

run_summary_img <- all_results_img %>%
  group_by(group, run) %>%
  summarise(
    mean_cindex_val_run = mean(mean_cindex_val, na.rm = TRUE),
    mean_AUC_val_run    = mean(mean_AUC_val,   na.rm = TRUE),
    .groups = "drop"
  )

# ---------------------------
# 3) 그룹별 요약 (mean ± SD, 95% CI, CV)
# ---------------------------

# 부트스트랩용 함수
boot_mean <- function(data, indices) {
  d <- data[indices]
  mean(d, na.rm = TRUE)
}

# C-index/AUC에 대해 group별 mean, SD, 95% CI, CV 계산
summary_run_table_img <- run_summary_img %>%
  group_by(group) %>%
  summarise(
    mean_cindex_val_mean = mean(mean_cindex_val_run, na.rm = TRUE),
    mean_cindex_val_sd   = sd(mean_cindex_val_run,   na.rm = TRUE),
    mean_AUC_val_mean    = mean(mean_AUC_val_run,    na.rm = TRUE),
    mean_AUC_val_sd      = sd(mean_AUC_val_run,      na.rm = TRUE),
    .groups = "drop"
  )

bootstrap_run_img <- run_summary_img %>%
  group_by(group) %>%
  summarise(
    # C-index bootstrap CI
    boot_cindex = list(boot(mean_cindex_val_run, statistic = boot_mean, R = 5000)),
    cindex_ci_lower = boot.ci(boot_cindex[[1]], type = "perc")$percent[4],
    cindex_ci_upper = boot.ci(boot_cindex[[1]], type = "perc")$percent[5],
    # AUC bootstrap CI
    boot_auc = list(boot(mean_AUC_val_run, statistic = boot_mean, R = 5000)),
    auc_ci_lower = boot.ci(boot_auc[[1]], type = "perc")$percent[4],
    auc_ci_upper = boot.ci(boot_auc[[1]], type = "perc")$percent[5],
    .groups = "drop"
  ) %>%
  select(group, cindex_ci_lower, cindex_ci_upper, auc_ci_lower, auc_ci_upper)

# CV 추가
summary_run_with_boot_img <- summary_run_table_img %>%
  left_join(bootstrap_run_img, by = "group") %>%
  mutate(
    cv_cindex = mean_cindex_val_sd / mean_cindex_val_mean,
    cv_auc    = mean_AUC_val_sd    / mean_AUC_val_mean
  ) %>%
  mutate(
    Cindex_mean_sd = sprintf("%.3f ± %.3f", mean_cindex_val_mean, mean_cindex_val_sd),
    AUC_mean_sd    = sprintf("%.3f ± %.3f", mean_AUC_val_mean,   mean_AUC_val_sd),
    Cindex_CI      = sprintf("[%.3f, %.3f]", cindex_ci_lower, cindex_ci_upper),
    AUC_CI         = sprintf("[%.3f, %.3f]", auc_ci_lower,    auc_ci_upper),
    CV_Cindex      = sprintf("%.3f", cv_cindex),
    CV_AUC         = sprintf("%.3f", cv_auc)
  ) %>%
  select(
    group,
    Cindex_mean_sd, Cindex_CI, CV_Cindex,
    AUC_mean_sd,    AUC_CI,    CV_AUC
  ) %>%
  arrange(group)

# 🔹 표 1개: summary_run_with_boot_img 객체로 사용
print(summary_run_with_boot_img)

# 원하면 CSV로 저장 (필요 없으면 주석 처리)
write.csv(summary_run_with_boot_img,
           "supp_table_imageonly_n4_n7_summary.csv",
           row.names = FALSE)

# ---------------------------
# 4) 그룹간 차이 & 안정성 검정 (p-value는 본문에서 언급용)
# ---------------------------

# C-index 차이 (Kruskal–Wallis)
kw_cindex_img <- kruskal.test(mean_cindex_val_run ~ group, data = run_summary_img)

# 안정성 (분산 차이, Levene)
levene_cindex_img <- car::leveneTest(mean_cindex_val_run ~ group, data = run_summary_img)

cat("\n[Image only] Kruskal–Wallis (C-index): p-value =",
    signif(kw_cindex_img$p.value, 3), "\n")
cat("[Image only] Levene's test (C-index): p-value =",
    signif(as.numeric(levene_cindex_img[1, "Pr(>F)"]), 3), "\n")

# ---------------------------
# 5) 그림 1개: run-level C-index boxplot (Image only, n4–n7)
# ---------------------------
# p-value 미리 뽑아두기
lev_p  <- levene_cindex_img[1, "Pr(>F)"]
kw_p   <- kw_cindex_img$p.value

p_box_cindex_img <- ggplot(run_summary_img, aes(x = group, y = mean_cindex_val_run)) +
  geom_boxplot(outlier.shape = NA, fill = "#4C72B0", alpha = 0.6) +
  # ✅ 점(지터) 크기/투명도 조정
  geom_jitter(width = 0.1, alpha = 0.7, size = 3) +
  labs(
    title = "Run-level Validation C-index (Image-only model)",
    subtitle = sprintf("Groups n4–n7; Kruskal–Wallis p = %.3f, Levene p = %.3f",
                       kw_p, lev_p),
    x     = "Feature Group",
    y     = "Mean C-index (Validation, per run)"
  ) +
  coord_cartesian(ylim = c(0.7, 0.95)) +
  theme_bw(base_size = 18) +   # ✅ 전체 기본 글자 크기 키우기
  theme(
    plot.title    = element_text(size = 22, face = "bold"),
    plot.subtitle = element_text(size = 18),
    axis.title.x  = element_text(size = 18),
    axis.title.y  = element_text(size = 18),
    axis.text.x   = element_text(size = 16),
    axis.text.y   = element_text(size = 16)
  )

print(p_box_cindex_img)



# 원하면 PNG로 저장 (필요 없으면 주석 처리)
# ggsave("supp_figure_imageonly_cindex_boxplot.png",
#        p_box_cindex_img, width = 6, height = 4, dpi = 300)

# 콘솔에서 summary_run_with_boot_img 확인하면
# 각 그룹별:
#  - run-level mean ± SD
#  - 부트스트랩 95% CI
#  - → "성능 비슷 / 안정성 비슷"을 정량적으로 보여주는 근거로 사용 가능
library(dplyr)
library(ggplot2)
library(ggpubr)

# 1) Feature factor 순서 고정
feature_levels <- c("Clinical only","Image only","Image + Clinical")
all_results$Feature <- factor(all_results$Feature, levels = feature_levels)

# 2) n7만 필터
all_results_n7 <- all_results %>%
  filter(group == "n7")

# ----------------------------------------
# (1) C-index: Feature set 비교 플롯
# ----------------------------------------
p_cindex_n7 <- ggplot(all_results_n7, aes(x = Feature, y = mean_cindex_val, fill = Feature)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  scale_fill_brewer(palette = "Set2", name = "Feature Set") +
  labs(
    title    = "Validation C-index by Feature Set (OS)",
    y        = "Mean C-index (validation)",
    x        = ""
  ) +
  coord_cartesian(ylim = c(0.4, 1)) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "top",          # 🔹 범례를 위로
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

print(p_cindex_n7)

# ----------------------------------------
# (2) AUC: Feature set 비교 플롯
# ----------------------------------------
p_auc_n7 <- ggplot(all_results_n7, aes(x = Feature, y = mean_AUC_val, fill = Feature)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  scale_fill_brewer(palette = "Set2", name = "Feature Set") +
  labs(
    title    = "Validation AUC by Feature Set (OS)",
    y        = "Mean AUC (validation)",
    x        = ""
  ) +
  coord_cartesian(ylim = c(0.4, 1)) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "top",          # 🔹 범례를 위로
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

print(p_auc_n7)

# 필요하면 각각 저장
# ggsave("supp_figure_n7_cindex_by_feature.png", p_cindex_n7,
#        width = 8, height = 6, dpi = 300)
# ggsave("supp_figure_n7_auc_by_feature.png",    p_auc_n7,
#        width = 8, height = 6, dpi = 300)


library(dplyr)
library(boot)

## 0) 전제: all_results 이미 존재
## Feature factor 순서 고정
feature_levels <- c("Clinical only","Image only","Image + Clinical")
all_results$Feature <- factor(all_results$Feature, levels = feature_levels)

## 1) n7만 필터
all_results_n7 <- all_results %>%
  filter(group == "n7",
         Feature %in% feature_levels)

## 2) Feature별 평균 / SD 계산  (C-index, AUC)
summary_n7 <- all_results_n7 %>%
  group_by(Feature) %>%
  summarise(
    mean_cindex = mean(mean_cindex_val, na.rm = TRUE),
    sd_cindex   = sd(mean_cindex_val,   na.rm = TRUE),
    mean_auc    = mean(mean_AUC_val,    na.rm = TRUE),
    sd_auc      = sd(mean_AUC_val,      na.rm = TRUE),
    .groups = "drop"
  )

## 3) 부트스트랩 CI 계산 함수
boot_mean <- function(data, indices) {
  d <- data[indices]
  mean(d, na.rm = TRUE)
}

## 4) 부트스트랩 CI (에러/degenerate case 방어 처리)
boot_ci_n7 <- all_results_n7 %>%
  group_by(Feature) %>%
  summarise(
    # C-index bootstrap object
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
    
    # AUC bootstrap object
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

## 5) 최종 요약 테이블 (mean ± SD + CI 포맷)
summary_with_ci_n7 <- summary_n7 %>%
  left_join(boot_ci_n7, by = "Feature") %>%
  mutate(
    Cindex_mean_sd = sprintf("%.3f ± %.3f", mean_cindex, sd_cindex),
    AUC_mean_sd    = sprintf("%.3f ± %.3f", mean_auc,    sd_auc),
    Cindex_CI      = ifelse(
      is.na(cindex_ci_lower),
      # degenerate case: CI 계산 실패하면 "NA" 혹은 mean으로 고정
      sprintf("NA"),
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

print(summary_with_ci_n7)
write.csv(summary_with_ci_n7, "internal_table.csv")

## 5) n7에서 Feature-set 간 p-value (Wilcoxon) -------------------

# C-index에 대한 쌍별 Wilcoxon
pw_cindex_n7 <- pairwise.wilcox.test(
  all_results_n7$mean_cindex_val,
  all_results_n7$Feature,
  p.adjust.method = "holm"
)

# AUC에 대한 쌍별 Wilcoxon
pw_auc_n7 <- pairwise.wilcox.test(
  all_results_n7$mean_AUC_val,
  all_results_n7$Feature,
  p.adjust.method = "holm"
)

# pairwise.wilcox.test 결과에서 원하는 쌍만 뽑기
# levels = c("Clinical only","Image only","Image + Clinical") 이라고 가정
p_cindex_clin_img     <- pw_cindex_n7$p.value["Image only",        "Clinical only"]
p_cindex_img_imgClin  <- pw_cindex_n7$p.value["Image + Clinical",  "Image only"]

p_auc_clin_img        <- pw_auc_n7$p.value["Image only",           "Clinical only"]
p_auc_img_imgClin     <- pw_auc_n7$p.value["Image + Clinical",     "Image only"]

p_table_n7 <- tibble::tibble(
  Comparison = c("Clinical only vs Image only",
                 "Image only vs Image + Clinical"),
  p_Cindex   = c(p_cindex_clin_img,    p_cindex_img_imgClin),
  p_AUC      = c(p_auc_clin_img,       p_auc_img_imgClin)
) %>%
  mutate(
    p_Cindex_fmt = sprintf("%.3f", p_Cindex),
    p_AUC_fmt    = sprintf("%.3f", p_AUC)
  )

print(p_table_n7)

write.csv(p_table_n7, "internal_table_pval.csv")
