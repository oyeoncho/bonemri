library(readxl)
library(dplyr)
library(xlsx)
library(tidyverse)
library(corrplot)
library(Hmisc)
library(Hmisc)  # rcorr 사용

cx_test <- readxl::read_excel(path= "interpretation.xlsx", sheet="interpretation", col_names = TRUE)

N_RUN <- 29

#전체 데이터 불러오기
beit0i <- read.csv("image_data/image_features/tuning/final_features/slice_20_dl0/beit0_backbone_features_merged.csv") %>%
  as_tibble()

exclude_vars <- c("Age", "pathology", "stage0", "recur_date", "recur", "recur1", "fu_date",
                  "survival", "RTF", "duration", "EQD2", 
                  "SCC0", "SCC2",  "Cyfra0", "Cyfra2")


sel_features =c("feat_213", "feat_266", "feat_499", "feat_436", "feat_2", "feat_327", "feat_391", "feat_519", "feat_173", "feat_715", 
                "feat_107", "feat_80", "feat_137", "feat_209", "feat_215", "feat_374", "feat_55", "feat_223", "feat_554", "feat_577", 
                "feat_109", "feat_583", "feat_657")


#리스트
all_cor_results <- list()

for (k in 0:N_RUN) {
  cat("Processing Run:", k, "\n")
  
  beit0 <- beit0i %>% filter(Run == k)
  
  merged <- cx_test %>%
    inner_join(beit0 %>% select(PatientID, all_of(sel_features)), by = "PatientID")
  
  # 상관분석 대상: clinical + sel_features (NA 제거)
  cor_input <- merged %>%
    select(-any_of(c("ID", "ID2", "RTID", exclude_vars))) %>%
    select(where(is.numeric), all_of(sel_features)) %>%
    drop_na()
  
  result <- Hmisc::rcorr(as.matrix(cor_input), type = "pearson")
  r_mat <- result$r
  p_mat <- result$P
  
  cor_df_list <- list()
  
  for (feat in sel_features) {
    # 대상 임상 변수만 필터
    target_vars <- setdiff(colnames(cor_input), sel_features)
    
    corr_vals <- r_mat[feat, target_vars]
    p_vals <- p_mat[feat, target_vars]
    
    df <- tibble(
      Run = k,
      feature = feat,
      variable = target_vars,
      correlation = as.numeric(corr_vals),
      p_value = as.numeric(p_vals)
    ) %>%
      filter(!is.na(correlation)) %>%
      filter(p_value < 0.05)
    
    cor_df_list[[feat]] <- df
  }
  
  cor_df_run <- bind_rows(cor_df_list)
  all_cor_results[[k + 1]] <- cor_df_run
}

# 전체 결합
final_correlation_all_runs <- bind_rows(all_cor_results)




# ─────────────────────────────────────────────
# ✅ 설정
K1 <- 10                # Run당 선택할 상위 변수 개수
N_RUN <- 30            # 총 Run 수 (0~29일 경우 N_RUN = 30)
TOTAL_TOP_ENTRIES <- N_RUN * K1  # 전체 entry 수

# ─────────────────────────────────────────────
# ✅ 1. Run별 상위 K개의 feature-variable 상관관계 추출
top_k_plot_data <- final_correlation_all_runs %>%
  group_by(Run) %>%
  slice_max(order_by = abs(correlation), n = K1) %>%
  ungroup() %>% mutate(Run = Run + 1)

# 1. 강조 대상 feature 지정
highlight_feats0 <- c("feat_436", "feat_519")
highlight_feats1 = c("feat_213", "feat_266", "feat_215")

# 긴 변수명을 줄바꿈으로 대체
# 박스색 지정: group0 = 검정배경, group1 = 빨강배경, 나머지 = 회색배경
top_k_plot_data <- top_k_plot_data %>%
  mutate(variable = if_else(variable == "miR-574-3p-LINC01003-ACOT9", 
                            "miR-574-3p\n-LINC01003\n-ACOT9", variable)) %>%
  mutate(
    fill_group = case_when(
      feature %in% highlight_feats0 ~ "Black: feat_436 / feat_519 (n7 only)",
      feature %in% highlight_feats1 ~ "Red: feat_213 / feat_266 / feat_215 (n6 only)",
      TRUE ~ "n5 only"
    ),
    font_face = if_else(feature %in% c(highlight_feats0, highlight_feats1), "bold", "plain")
  )

# ─────────────────────────────────────────────
# ✅ 2. Heatmap 시각화 (Run별 상위 K개)
ggplot(top_k_plot_data, aes(x = variable, y = feature)) +
  geom_tile(aes(fill = fill_group), color = "white") +
  geom_text(aes(label = round(correlation, 2), fontface = font_face), 
            color = "blue", size = 3.5) +
  
  scale_fill_manual(
    name = "Highlighted Features",
    values = c(
      "Black: feat_436 / feat_519 (n7 only)" = "black",
      "Red: feat_213 / feat_266 / feat_215 (n6 only)" = "red",
      "n5 only" = "grey90"
    ),
    breaks = c(
      "Black: feat_436 / feat_519 (n7 only)", 
      "Red: feat_213 / feat_266 / feat_215 (n6 only)", 
      "n5 only"
    )
  ) +
  
  facet_wrap(~Run, scales = "free") +
  theme_minimal(base_size = 20) +
  labs(
    title = paste("Top", K1, "Correlations per Run"),
    x = "Plasma exosomal RNAs and hemogram", y = "Common features (n5, ≥5 reps) from 6 variants of 494 images"
  ) +
  theme(
    legend.position = "top",  # ✅ 범례 위쪽으로 이동
    legend.title = element_text(size = 20, face="bold"),
    legend.text = element_text(size = 20, face="bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face ="bold"),
    strip.text = element_text(size = 20),
    plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.title.x = element_text(size = 30, face = "bold")
  )

# ─────────────────────────────────────────────
# ✅ 3. 변수별 등장 횟수 및 비율 계산
top_variable_freq <- top_k_plot_data %>%
  count(variable, name = "count") %>%
  mutate(percentage = 100 * count / TOTAL_TOP_ENTRIES) %>%
  arrange(desc(percentage))

# ─────────────────────────────────────────────
# ✅ 4. 바 차트: Top 변수 등장 비율 시각화
top_variable_freq %>%
  ggplot(aes(x = reorder(variable, percentage), y = percentage)) +
  geom_col(fill = "gray") +
  coord_flip() +
  labs(
    title = paste("Top feature related variables in Top", K1, "per Run \n(", TOTAL_TOP_ENTRIES, " entries)"),
    x = "Plasma exosomal RNAs \nand hemogram", y = "Relative Frequency (%)"
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 20),     # ✅ 20% 단위 표시
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 12)+
  theme(
    legend.position = "top",  # ✅ 범례 위쪽으로 이동
    axis.text = element_text(angle = 45, hjust = 1, size = 25, face="bold"),
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 35, face = "bold"),
    axis.title.x = element_text(size = 35, face = "bold")
  )

library(dplyr)
library(ggplot2)
library(scales)

# 그룹 정의
n5_only <- c("feat_499","feat_2","feat_327","feat_391", "feat_173","feat_715","feat_107","feat_80","feat_137","feat_209",
             "feat_374","feat_55","feat_223","feat_554","feat_577", "feat_109","feat_583","feat_657") #18

n6_only <- c("feat_215","feat_213","feat_266") #3
n7_only <- c("feat_436","feat_519")# 2



# 분석할 variable 원래 표기 맞춰 설정
target_var <- "miR-574-3p\n-LINC01003\n-ACOT9" 

# 데이터 필터링
dat_target <- top_k_plot_data %>%
  filter(variable == target_var,
         feature %in% c(n5_only, n6_only, n7_only)) %>%
  mutate(group = case_when(
    feature %in% n7_only ~ "n7 only",
    feature %in% n6_only ~ "n6 only",
    feature %in% n5_only ~ "n5 only"
  ))

n_runs <- length(unique(top_k_plot_data$Run))

# 그룹별 feature 빈도 및 평균 r
summary_table <- dat_target %>%
  group_by(group, feature) %>%
  summarise(
    freq_runs = n(),  # 해당 variable과 함께 등장한 횟수
    pct_runs  = freq_runs / n_runs,
    mean_r    = mean(correlation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(group, desc(freq_runs))

print(summary_table)

p_bar <- summary_table %>%
  group_by(group) %>%
  mutate(feature = fct_reorder(feature, pct_runs, .desc = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = feature, y = pct_runs, fill = group)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = percent(pct_runs, accuracy = 1)),
            hjust = 0.7, size = 10, fontface = "bold") +  # 글자 크기 & 굵기
  scale_y_continuous(labels = function(x) round(x * 100, 1), 
                     expand = expansion(mult = c(0, 0.05))) +
  facet_wrap(~ group, scales = "free_x") +
  labs(title = "Association with miR-574-3p-LINC01003-ACOT9 (log2FC)",
       subtitle = "% of runs (relative frequency across 30 runs)",
       x = "Feature", y = "% of runs") +
  coord_flip() +
  theme_bw(base_size = 20) +  # 기본 글자 크기 증가
  theme(
    axis.text  = element_text(face = "bold", size =25),     # 축 숫자 굵게
    strip.text = element_text(face = "bold", size =30),     # facet 제목 굵게
    axis.title = element_text(face = "bold", size =35),     # 축 라벨 굵게
    plot.title = element_text(face = "bold", size = 35), # 메인 타이틀 굵게/크게
    plot.subtitle = element_text(size = 30, face = "bold") # 부제목 굵게/크게
  )

print(p_bar)


# 그룹별 합산 비율 계산
group_summary <- dat_target %>%
  group_by(group) %>%
  summarise(
    freq_runs = n_distinct(Run),   # 해당 그룹에서 등장한 run 수
    pct_runs = freq_runs / 30,     # 전체 30 run 대비 비율
    .groups = "drop"
  )

print(group_summary)

# 시각화
p_group <- ggplot(group_summary, aes(x = group, y = pct_runs, fill = group)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = percent(pct_runs, accuracy = 1)),
            vjust = -0.3, size = 10, fontface = "bold") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Association with miR-574-3p-LINC01003-ACOT9 (log2FC)\n in group level",
    subtitle = "% of runs (relative frequency across 30 runs)",
    x = "Group", y = "% of runs"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text  = element_text(face = "bold", size=30),
    axis.title = element_text(face = "bold", size=30),
    plot.title = element_text(face = "bold", size = 35),
    plot.subtitle = element_text(size = 30, face = "bold")
  )

print(p_group)





