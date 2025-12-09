library(readxl)
library(dplyr)
library(tidyverse)
library(moonBook)
library(ggpubr)
# 그룹 정의
groups <- paste0("n", 4:9, "_30_30")

# 결과 저장 리스트
all_results_all <- list()
all_results_long_all <- list()

for (grp in groups) {
  group_results <- map_dfr(1:30, function(i) {
    auc_path <- sprintf("image_data/survival_model/mixture_non_fix/tune/beit_original/results/generalization/test/dl0/%s/raw_auc_per_time_run%02d.csv", grp, i)
    cindex_path <- sprintf("image_data/survival_model/mixture_non_fix/tune/beit_original/results/generalization/test/dl0/%s/raw_cindex_per_time_run%02d.csv", grp, i)
    
    a <- read.csv(auc_path) %>%
      as_tibble() %>%
      rename(Feature = Feature.Set, Time = Time..Months., AUC_train = AUC..Train., AUC_val = AUC..Val.)
    
    b <- read.csv(cindex_path) %>%
      as_tibble() %>%
      rename(Feature = Feature.Set, Time = Time..Months., C_index_train = C.index..Train., C_index_val = C.index..Val.)
    
    # "Overall" 제외한 시점만 필터링
    a_time <- a %>% filter(Time != "Overall")
    b_time <- b %>% filter(Time != "Overall")
    
    combined <- bind_cols(b_time %>% select(C_index_train, C_index_val), a_time) %>%
      group_by(Feature) %>%
      summarise(
        mean_AUC_train = mean(AUC_train),
        mean_AUC_val = mean(AUC_val),
        mean_cindex_train = mean(C_index_train),
        mean_cindex_val = mean(C_index_val),
        .groups = "drop"
      ) %>%
      mutate(
        diff_AUC = mean_AUC_train - mean_AUC_val,
        diff_cindex = mean_cindex_train - mean_cindex_val,
        run = i,
        group = grp
      )
  })
  
  group_results_long <- group_results %>%
    select(group, run, Feature, diff_cindex, mean_cindex_val) %>%
    pivot_longer(cols = c(diff_cindex, mean_cindex_val), names_to = "Metric", values_to = "Value")
  
  all_results_all[[grp]] <- group_results
  all_results_long_all[[grp]] <- group_results_long
}

# 전체 결합
all_results <- bind_rows(all_results_all)
all_results_long <- bind_rows(all_results_long_all)
all_results$group <- gsub("_30_30", "", all_results$group)
all_results_long$group <- gsub("_30_30", "", all_results_long$group)

feature_levels <- c("Image only", "Clinical only", "Image + Clinical")
all_results$Feature <- factor(all_results$Feature, levels = feature_levels)

ggplot(all_results, aes(x = Feature, y = mean_cindex_val, fill = Feature)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  facet_wrap(~ group) +
  scale_fill_brewer(palette = "Set2", name = "Feature Set") +  # 범례 제목 추가
  labs(
    title = "Validation C-index by Feature Set and Group \n(Beit_resize_o)",
    subtitle = expression("nX group (≥ X repetitions): n4, n5, n6, n7, n8, n9 from 494 images"),
    y = "Mean C-index (val)",
    x = ""
  ) +
  coord_cartesian(ylim = c(0.4, 1)) +
  theme_minimal(base_size = 16) +  # 전체 기본 글꼴 크기 키움
  theme(
    plot.title = element_text(size = 30, face = "bold"),
    plot.subtitle = element_text(size = 25, margin = margin(b = 10)),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 20, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 20),
    strip.text = element_text(size = 25, face = "bold"),  # facet 제목
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20)
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      c("Image only", "Clinical only"),
      c("Image only", "Image + Clinical"),
      c("Clinical only", "Image + Clinical")
    ),
    label = "p.signif", size=7,
    y.position = c(0.92, 0.95, 0.98),
    step.increase = 0.05
  )

# ✅ Threshold 별 feature별 비율
thresholds <- seq(0.7, 0.8, 0.05)

df_plot <- map_dfr(thresholds, function(thr) {
  all_results %>%
    filter(mean_cindex_val > thr) %>%
    group_by(group, Feature) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(threshold = thr, n_r = n / 30)
})

df_plot$threshold <- as.numeric(df_plot$threshold)

ggplot(df_plot, aes(x = group, y = n_r, fill = Feature)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  facet_wrap(~ threshold, labeller = label_bquote(plain(">") * .(threshold))) +
  labs(
    title = "Proportion of Runs Above Threshold by Feature Set \n(Beit_resize_o)",
    subtitle = "Validation C-index of selected features",
    y = "Proportion of Runs", x = expression(bold("nX group (≥ X repetitions): n4, n5, n6, n7, n8, n9"))
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 16) +  # 전체 기본 글꼴 크기 키움
  theme(
    legend.position="top",
    plot.title = element_text(size = 30, face = "bold"),
    plot.subtitle = element_text(size = 25, margin = margin(b = 10)),
    axis.title.x = element_text(size = 30, face = "bold"),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.text.x = element_text(size = 25, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 25),
    strip.text = element_text(size = 25, face = "bold"),  # facet 제목
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20)
  ) +
  ylim(0, 1)

###
library(moonBook)
out=mytable(Feature ~ mean_cindex_train+mean_cindex_val + diff_cindex + mean_AUC_train+ mean_AUC_val + diff_AUC, data = all_results,digits=3, method=3)
#mycsv(out, file="time_wise_mean.csv")

#########

target_feature <- "Image only"

filtered_plot_data <- all_results_long %>%
  filter(Feature == target_feature) %>%
  filter(Metric %in% c("mean_cindex_val", "diff_cindex"))


ggplot(filtered_plot_data, aes(x = group, y = Value, fill = Metric)) + 
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0.8, linetype = "solid", color = "red", linewidth = 1) +  
  facet_wrap(~ run, ncol = 5) +  
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = paste("Run-wise Metrics for", target_feature, "(Beit_resize_o)"),
    subtitle = "Validation C-index of selected features",
    y = "Value", 
    x = expression(bold("nX group (≥ X repetitions): n4, n5, n6, n7, n8, n9"))
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.title.x = element_text(size = 30, face = "bold"),
    axis.title.y = element_text(size = 30, face = "bold"),
    plot.title = element_text(size = 30, face = "bold"),
    plot.subtitle = element_text(size = 25),
    strip.text = element_text(size = 25, face = "bold"),
    legend.title = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 20)
  )

#####
#########
# 📦 필요한 패키지 로드
library(readxl)
library(dplyr)
library(tidyverse)
library(moonBook)
library(ggpubr)

# 📁 그룹 정의
groups <- paste0("n", 4:9, "_30_30")

# 📌 결과 저장 리스트
all_results_all <- list()

# 🔁 그룹 루프
for (grp in groups) {
  cat("Processing group:", grp, "\n")
  
  group_results <- map_dfr(1:30, function(i) {
    # 📁 파일 경로 지정
    auc_path <- sprintf("image_data/survival_model/mixture_non_fix/tune/beit_original/results/generalization/test/dl0/%s/raw_auc_per_time_run%02d.csv", grp, i)
    cindex_path <- sprintf("image_data/survival_model/mixture_non_fix/tune/beit_original/results/generalization/test/dl0/%s/raw_cindex_per_time_run%02d.csv", grp, i)
    
    # 🔄 데이터 로딩 및 컬럼명 정리
    a <- read.csv(auc_path) %>%
      as_tibble() %>%
      rename(
        Feature = Feature.Set,
        Time = Time..Months.,
        AUC_train = AUC..Train.,
        AUC_val = AUC..Val.
      )
    
    b <- read.csv(cindex_path) %>%
      as_tibble() %>%
      rename(
        Feature = Feature.Set,
        Time = Time..Months.,
        C_index_train = C.index..Train.,
        C_index_val = C.index..Val.
      )
    
    
    
    # 🧩 시간 + 피처별 평균 계산
    combined <- bind_cols(
      b %>% select(Time, C_index_train, C_index_val),
      a %>% select(Feature, Time, AUC_train, AUC_val)
    ) %>%
      # Time 중복 열 제거
      select(Time = Time...1, Feature, C_index_train, C_index_val, AUC_train, AUC_val) %>%
      group_by(Time, Feature) %>%
      summarise(
        mean_AUC_train = mean(AUC_train, na.rm = TRUE),
        mean_AUC_val   = mean(AUC_val, na.rm = TRUE),
        mean_cindex_train = mean(C_index_train, na.rm = TRUE),
        mean_cindex_val   = mean(C_index_val, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        diff_AUC = mean_AUC_train - mean_AUC_val,
        diff_cindex = mean_cindex_train - mean_cindex_val,
        run = i,
        group = grp
      )
    
    return(combined)
  })
  
  # 💾 그룹별 결과 저장
  all_results_all[[grp]] <- group_results
}

# 📦 최종 결과 통합
final_results <- bind_rows(all_results_all)
image  <- final_results %>% filter(Feature=="Image only")
out =mytable(Time ~ mean_cindex_train+mean_cindex_val + diff_cindex + mean_AUC_train+ mean_AUC_val + diff_AUC, data = image, digits=3, method=3)
mycsv(out, file="image.csv")
i_c  <- final_results %>% filter(Feature=="Image + Clinical")
out=mytable(Time ~ mean_cindex_train+mean_cindex_val + diff_cindex + mean_AUC_train+ mean_AUC_val + diff_AUC, data = i_c, digits=3, method=3)
mycsv(out, file="i_c.csv")
