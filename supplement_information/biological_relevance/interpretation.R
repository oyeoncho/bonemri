library(readxl)
library(dplyr)
library(xlsx)
library(tidyverse)
library(corrplot)
library(Hmisc)
library(Hmisc)  # rcorr 사용

cx_test <- readxl::read_excel(path= "interpretation.xlsx", sheet="interpretation", col_names = TRUE)

## ─────────────────────────────────────────────
## 2. BEiT feature 로드 + n7 feature만 사용 (feat_436, feat_519)
## ─────────────────────────────────────────────

# 전체 이미지 feature (Run 컬럼 포함)
beit0i <- read.csv(
  "image_data/image_features/tuning/final_features/slice_20_dl0/beit0_backbone_features_merged.csv"
) %>%
  as_tibble()

# n7 only features
sel_features <- c("feat_436", "feat_519")

# run index: 0 ~ 29 (총 30 run이라고 가정)
N_RUN <- 30

# 임상/치료/혈액 변수 중에서 상관분석에서 제외할 것들 (필요시 수정 가능)
exclude_vars <- c(
  "Age", "pathology", "stage0", "recur_date", "recur", "recur1", "fu_date",
  "survival", "RTF", "duration", "EQD2",
  "SCC0", "SCC2", "Cyfra0", "Cyfra2"
)

## ─────────────────────────────────────────────
## 3. Run별 상관분석: feat_436, feat_519 vs RNA/임상 (p < 0.05만 저장)
## ─────────────────────────────────────────────

all_cor_results <- list()

for (k in 0:(N_RUN - 1)) {
  cat("Processing Run:", k, "\n")
  
  beit0_k <- beit0i %>% filter(Run == k)
  
  merged <- cx_test %>%
    inner_join(beit0_k %>% select(PatientID, all_of(sel_features)), by = "PatientID")
  
  # 상관분석 대상: numeric 임상/바이오마커 + feat_436, feat_519
  cor_input <- merged %>%
    select(-any_of(c("ID", "ID2", "RTID", exclude_vars))) %>%
    select(where(is.numeric), all_of(sel_features)) %>%
    drop_na()
  
  if (nrow(cor_input) < 5) {
    # 표본이 너무 적으면 스킵
    all_cor_results[[k + 1]] <- tibble(
      Run        = integer(0),
      feature    = character(0),
      variable   = character(0),
      correlation = numeric(0),
      p_value    = numeric(0)
    )
    next
  }
  
  res <- Hmisc::rcorr(as.matrix(cor_input), type = "pearson")
  r_mat <- res$r
  p_mat <- res$P
  
  cor_df_list <- list()
  
  for (feat in sel_features) {
    clinical_vars <- setdiff(colnames(cor_input), sel_features)
    
    r_vals <- r_mat[feat, clinical_vars]
    p_vals <- p_mat[feat, clinical_vars]
    
    df <- tibble(
      Run        = k,
      feature    = feat,
      variable   = clinical_vars,
      correlation = as.numeric(r_vals),
      p_value    = as.numeric(p_vals)
    ) %>%
      filter(!is.na(correlation), p_value < 0.05)
    
    cor_df_list[[feat]] <- df
  }
  
  all_cor_results[[k + 1]] <- bind_rows(cor_df_list)
}

final_correlation_all_runs <- bind_rows(all_cor_results)

## ─────────────────────────────────────────────
## 4. n7 전용 heatmap: Run × feature × variable
## ─────────────────────────────────────────────

# Run을 1~30으로 (보기 좋게)
K1 <- 5  # 각 Run당 상관계수 상위 K개만 표시 (원하면 변경/삭제)
top_k_plot_data <- final_correlation_all_runs %>%
  group_by(Run) %>%
  slice_max(order_by = abs(correlation), n = K1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Run = Run + 1)  # 0-based → 1-based

# miR-574-3p-LINC01003-ACOT9 표기 줄바꿈
top_k_plot_data <- top_k_plot_data %>%
  mutate(
    variable = if_else(
      variable == "miR-574-3p-LINC01003-ACOT9",
      "miR-574-3p\n-LINC01003\n-ACOT9",
      variable
    )
  )

# Heatmap: 상관계수에 따라 색을 주고, 값도 텍스트로 표시
p_heat_n7 <- ggplot(top_k_plot_data,
                    aes(x = variable, y = feature)) +
  geom_tile(aes(fill = correlation), color = "white") +
  geom_text(aes(label = round(correlation, 2)),
            size = 3.5, fontface = "bold", color = "black") +
  scale_fill_gradient2(
    name = "r",
    low = "blue", mid = "white", high = "red",
    midpoint = 0
  ) +
  facet_wrap(~ Run, scales = "free") +
  theme_minimal(base_size = 16) +
  labs(
    title = paste0("Top ", K1, " significant correlations per run"),
    subtitle = "n7 features: feat_436 & feat_519",
    x = "Plasma exosomal RNAs / hemogram",
    y = "n7 features"
  ) +
  theme(
    legend.position   = "top",
    axis.text.x       = element_text(angle = 45, hjust = 1),
    strip.text        = element_text(face = "bold")
  )

print(p_heat_n7)

## ─────────────────────────────────────────────
## 5. 변수별 등장 빈도 (n7에서 어떤 biomarker가 자주 뜨는지)
## ─────────────────────────────────────────────

TOTAL_TOP_ENTRIES <- length(unique(top_k_plot_data$Run)) * K1

top_variable_freq_n7 <- top_k_plot_data %>%
  count(variable, name = "count") %>%
  mutate(
    percentage = 100 * count / TOTAL_TOP_ENTRIES
  ) %>%
  arrange(desc(percentage))

# bar plot: 어떤 변수들이 얼마나 자주 등장했는지
p_var_freq_n7 <- top_variable_freq_n7 %>%
  ggplot(aes(x = reorder(variable, percentage), y = percentage)) +
  geom_col(fill = "gray60") +
  coord_flip() +
  labs(
    title = paste0("Top variables associated with n7 (IDs 436 and 519)\nTop ", K1, " entries per run"),
    x = "Variables (exosomal RNA / hemogram)",
    y = "Relative frequency (%)"
  ) +
  scale_y_continuous(
    limits = c(0, 40),
    breaks = seq(0, 50, 10),
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_text(angle = 20, hjust = 1, face = "bold"),  # ★ 수정
    axis.text.x = element_text(face = "bold"),
    plot.title  = element_text(face = "bold", hjust = 0.5)
  )


print(p_var_freq_n7)

## ─────────────────────────────────────────────
## 6. miR-574-3p–LINC01003–ACOT9 (axis) 에 대한 feat_436 vs feat_519 비교
## ─────────────────────────────────────────────
target_var <- "miR-574-3p\n-LINC01003\n-ACOT9"

dat_target_n7 <- top_k_plot_data %>%
  filter(variable == target_var)

summary_n7_target <- dat_target_n7 %>%
  group_by(feature) %>%
  summarise(
    freq_runs = n(),                     # 유의하게 잡힌 run 수
    pct_runs  = freq_runs / N_RUN,       # 전체 30 run 대비 비율
    mean_abs_r = mean(abs(correlation), na.rm = TRUE),   # 절대값 r 평균
    sd_abs_r   = sd(abs(correlation), na.rm = TRUE),     # 절대값 r 표준편차
    median_abs_r = median(abs(correlation), na.rm = TRUE),  # 절대값 r 중앙값
    .groups = "drop"
  ) %>%
  arrange(desc(freq_runs), desc(mean_abs_r))

print(summary_n7_target)

# bar plot: feat_436 vs feat_519의 run 등장 비율
p_feat_bar_n7 <- summary_n7_target %>%
  ggplot(aes(x = feature, y = pct_runs, fill = feature)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = percent(pct_runs, accuracy = 1)),
            vjust = -0.3, size = 6, fontface = "bold") +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title    = "Association with miR-574-3p–LINC01003–ACOT9 (log2FC)",
    subtitle = "Relative frequency across 30 runs (n7 features)",
    x = "n7 feature", y = "% of runs"
  ) +
  theme_bw(base_size = 18) +
  theme(
    axis.text  = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(face = "bold")
  )

print(p_feat_bar_n7)
