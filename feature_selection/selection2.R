library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)

# ------------------------------
# 설정
# ------------------------------
n_min_values <- 4:7
run_file <- "dh11_run02.csv"

# 그룹 정의 (필요시 Group3 추가 가능)
group1 <- c("beit0", "beit_resize", "beit")
group2 <- c("beit_o", "beit0_o", "beit_original")
group3 <- c()  # 비워두면 무시됨

groups <- list(Group1 = group1, Group2 = group2)
if (length(group3) > 0) {
  groups$Group3 <- group3
}

# ------------------------------
# 그룹별 공통 feat 계산
# ------------------------------
group_common_feats_by_n <- list()

for (n_min in n_min_values) {
  group_common_feats_by_n[[paste0("n", n_min)]] <- map(groups, function(dir_group) {
    feat_cols_list <- map(dir_group, function(base_dir) {
      load_dir <- sprintf("./dataset/%s/test/dl0/n%d_30_30", base_dir, n_min)
      file_path <- file.path(load_dir, run_file)
      df <- read.csv(file_path)
      grep("feat", names(df), value = TRUE)
    })
    Reduce(intersect, feat_cols_list)
  })
}

# ------------------------------
# 비교 및 시각화 함수
# ------------------------------
compare_and_plot <- function(n_key, top_k = 30) {
  cat("\n\n===== 📊 분석 for ", n_key, " =====\n")
  
  feats_list <- group_common_feats_by_n[[n_key]]
  group_names <- names(feats_list)
  common_all <- Reduce(intersect, feats_list)
  
  # 고유 feature 계산
  diffs <- list()
  for (i in seq_along(group_names)) {
    current_group <- group_names[i]
    others <- setdiff(group_names, current_group)
    diffs[[paste0("Only_", current_group)]] <- setdiff(
      feats_list[[current_group]],
      unique(unlist(feats_list[others]))
    )
  }
  
  # 공통 및 고유 feature 출력
  cat("✅ 공통 feature (모든 그룹):\n")
  print(common_all)
  
  cat("❌ 그룹별 고유 feature:\n")
  print(diffs)
  
  # 표 생성
  all_feats <- sort(unique(unlist(feats_list)))
  feat_df <- tibble(feat = all_feats)
  for (g in group_names) {
    feat_df[[g]] <- all_feats %in% feats_list[[g]]
  }
  
  feat_df <- feat_df %>%
    mutate(count = rowSums(select(., all_of(group_names)))) %>%
    arrange(desc(count))
  
  # ✅/❌ 변환
  feat_display <- feat_df %>%
    mutate(across(all_of(group_names), ~ ifelse(.x, "✅", "❌")))
  
  cat("📋 ✅/❌ 표시 표:\n")
  print(feat_display)
  
  # 시각화
  p <- ggplot(feat_df[1:min(top_k, nrow(feat_df)), ], aes(x = reorder(feat, count), y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = paste0("Top ", top_k, " Feats in ", n_key), x = "Feature", y = "Included Group Count") +
    theme_minimal(base_size = 14)
  
  print(p)
}

# ------------------------------
# 실행
# ------------------------------
for (n_key in names(group_common_feats_by_n)) {
  compare_and_plot(n_key)
}

for (n_key in names(group_common_feats_by_n)) {
  feats_list <- group_common_feats_by_n[[n_key]]
  common_all <- Reduce(intersect, feats_list)
  
  # feat_숫자 패턴만 추출
  numeric_feats <- grep("^feat_\\d+$", common_all, value = TRUE)
  
  # 숫자 부분만 추출
  feat_nums <- as.integer(gsub("feat_", "", numeric_feats))
  feat_nums <- feat_nums[!is.na(feat_nums)]
  
  cat(paste0("# ", n_key, " overlapping features across 6 groups: [", paste(sort(feat_nums), collapse = ", "), "]\n"))
}

for (n_key in names(group_common_feats_by_n)) {
  feats_list <- group_common_feats_by_n[[n_key]]
  common_all <- Reduce(intersect, feats_list)
  
  quoted_feats <- paste0('"', common_all, '"')
  cat("\n🔍 ", n_key, "overlapping features across 6 groups:\n")
  cat("[", paste(quoted_feats, collapse = ", "), "]\n", sep = "")
}

