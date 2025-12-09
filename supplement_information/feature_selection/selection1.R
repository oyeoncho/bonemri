library(readxl)
library(dplyr)
library(xlsx)
library(tidyverse)
library(corrplot)
library(Hmisc)
library(moonBook)

cx <- read.csv(file = "cx.csv") %>% as_tibble() %>% rename("ID1"="PatientID")



N_RUN <- 29

#전체 데이터 불러오기
beit1i <- read.csv("./features/slice_20_dl0/beit1_backbone_features_merged.csv") %>%
  as_tibble()


# 1️⃣ 전체 수치형 피처만 추출해서 공통 스케일 기준 계산
numeric_cols <- beit1i %>% select(where(is.numeric)) %>% select(-Run)  # Run 제외
global_mean <- apply(numeric_cols, 2, mean, na.rm = TRUE)
global_sd   <- apply(numeric_cols, 2, sd, na.rm = TRUE)

# 2️⃣ 스케일 함수 정의
scale_global <- function(df) {
  df_scaled <- as_tibble(sweep(df, 2, global_mean, "-"))
  df_scaled <- as_tibble(sweep(df_scaled, 2, global_sd, "/"))
  return(df_scaled)
}
library(tidyverse)
library(Hmisc)  # rcorr 사용
library(ggplot2)

# 1️⃣ 전체 수치형 피처만 추출해서 공통 스케일 기준 계산
numeric_cols <- beit1i %>% select(where(is.numeric)) %>% select(-Run)  # Run 제외
global_mean <- apply(numeric_cols, 2, mean, na.rm = TRUE)
global_sd   <- apply(numeric_cols, 2, sd, na.rm = TRUE)

# 2️⃣ 스케일 함수 정의
scale_global <- function(df) {
  df_scaled <- as_tibble(sweep(df, 2, global_mean, "-"))
  df_scaled <- as_tibble(sweep(df_scaled, 2, global_sd, "/"))
  return(df_scaled)
}

# 3️⃣ 반복문 시작
b_list <- list()

for (k in 0:N_RUN) {
  
  # Run별 데이터 추출
  beit1 <- beit1i %>% filter(Run == k)
  
  # 수치형 피처만 추출하고 공통 스케일 적용
  numeric_feats <- beit1 %>% select(where(is.numeric)) %>% select(-Run)
  numeric_feats_scaled <- scale_global(numeric_feats)
  
  # PatientID 유지
  numeric_feats_scaled <- numeric_feats_scaled %>%
    mutate(PatientID = beit1$PatientID)
  
  # pre 추출
  pre <- numeric_feats_scaled %>%
    filter(str_detect(PatientID, "pre")) %>%
    rename("ID1" = "PatientID") %>%
    mutate(ID = ifelse(substr(ID1, 1, 1) == "P", substr(ID1, 1, 4), substr(ID1, 1, 5))) %>%
    select(ID, ID1, where(is.numeric))
  
  # 병합
  dh1 <- cx %>%
    inner_join(pre, by = "ID1") 
  
  # 상관 분석용 준비
  all_feats <- numeric_feats_scaled %>% select(-PatientID)
  a <- dh1 %>% select(survival, colnames(all_feats))
  
  # NA/표준편차 0 제거
  numeric_data_clean <- a[, sapply(a, function(x) {
    sd_val <- tryCatch(sd(x, na.rm = TRUE), error = function(e) NA)
    !is.na(sd_val) && sd_val != 0
  })]
  
  print(paste0("Run = ", k, " / 남은 열 수: ", ncol(numeric_data_clean)))
  
  # 상관 분석
  result <- Hmisc::rcorr(as.matrix(numeric_data_clean), type = "pearson")
  cor_matrix <- result$r
  cor_matrix[result$P >= 0.05] <- 0
  
  corm <- cor_matrix %>%
    as_tibble() %>%
    mutate(Name = colnames(cor_matrix)) %>%
    select(Name, survival) %>%
    filter(!Name %in% c("Age", "pathology", "stage", "fu_date", "recur_date", "survival", "recur", "css"))
  
  b <- corm %>%
    arrange(desc(abs(survival))) %>% slice_head(n = 100) %>% mutate(order = 1:100)
  
  b_list[[k + 1]] <- b %>% mutate(Run = k)
}

# 전체 통합
b_all <- bind_rows(b_list)

library(dplyr)
library(ggplot2)

b_all %>%
  group_by(Run) %>% 
  filter(abs(survival) > 0) %>%
  summarise(survival_mean = mean(abs(survival)), .groups = "drop") %>%
  mutate(Run = Run + 1) %>%  # 🔹 Run을 n+1로 변환
  arrange(desc(survival_mean)) %>%
  ggplot(aes(x = reorder(Run, -survival_mean), y = survival_mean)) +
  geom_col() +
  geom_hline(yintercept = 0.7, color = "red", size = 2, linetype=1) +
  geom_hline(yintercept = 0.5, color = "blue", size = 2, linetype=2) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  labs(x = "Run index", y = "Mean |R| of top 100 features") +
  ggtitle("Beit1 - No resize, Domain-Adaptive, Lower") +
  theme_minimal(base_size = 30) +
  theme(
    axis.text.x = element_text(size = 25, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 25),
    axis.title.y = element_text(size = 35, face="bold"),
    axis.title.x = element_text(size = 35, face="bold"),
    plot.title = element_text(size = 35, face = "bold", hjust = 0.5)
  )


top_runs <- b_all %>%
  group_by(Run) %>%
  summarise(survival_mean = mean(abs(survival)), .groups = "drop") %>%
  arrange(desc(survival_mean)) 

b_top <- b_all %>% filter(abs(survival)>0) %>%
  semi_join(top_runs, by = "Run")  


library(ggplot2)

b_top %>%
  count(Name, sort = TRUE) %>%
  top_n(20, n) %>%
  ggplot(aes(x = reorder(Name, n), y = n)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Feature Name", y = "Frequency in Top 20 Runs")

name_count <- b_top %>%
  count(Name, sort = TRUE)

#final <- name_count %>% filter(n >= 8) 
#write.csv(final, "./image_data/dataset/beit1/test/n8_f_final.csv")


# 4~9까지 필final# 4~9까지 필터링 조건별로 행 수 요약
result <- sapply(4:9, function(k) {
  name_count %>%
    filter(n >= k) %>%
    nrow()
})

# 결과를 데이터프레임으로 보기 좋게 정리
table <- data.frame(
  n_cutoff = 4:9,
  row_count = result
)

write.csv(table, "./dataset/beit1/test/count_num.csv")


library(purrr)

for (n_min in c(4:9)) {
  # 1. sel 정의 (n회 이상 등장한 feature만 추출)
  
  sel <- b_top %>%
    count(Name, sort = TRUE) %>%
    filter(n >= n_min)
  
  # 2. top feature 이름 추출
  top_names <- sel %>% pull(Name)
  
  # 3. dh11 리스트 초기화
  dh11_list <- list()
  
  for (k in 0:N_RUN) {
    # beit1에서 해당 Run 데이터만 추출
    beit1 <- beit1i %>% filter(Run == k)
    
    beit1 <- beit1 %>%
      select(where(is.numeric)) %>%
      as_tibble() %>%
      mutate(PatientID = beit1$PatientID)
    
    # pre 데이터 추출 및 ID 생성
    pre <- beit1 %>%
      filter(str_detect(PatientID, "pre")) %>%
      rename("ID1" = "PatientID") %>%
      mutate(ID = ifelse(substr(ID1,1,1) == "P", substr(ID1, 1, 4), substr(ID1, 1, 5))) %>%
      select(ID, ID1, where(is.numeric))
    
    # clinical 데이터와 병합
    dh1 <- cx %>%
      inner_join(pre, by = "ID1") %>%
      select(-ID.x, -ID.y, -ID2)
    
    # 유효한 변수만 선택
    valid_names <- intersect(top_names, colnames(dh1))
    
    # 필요한 열만 선택하여 저장
    dh11_k <- dh1 %>% select(ID1:EQD2, all_of(valid_names))
    dh11_list[[k + 1]] <- dh11_k
  }
  
  # 4. 결과 저장 경로 생성
  save_dir <- sprintf("./image_data/dataset/beit1/test/dl0/n%d_30_30", n_min)
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 5. CSV 파일 저장
  walk2(
    dh11_list,
    seq_along(dh11_list),
    ~ write.csv(.x, file = file.path(save_dir, sprintf("dh11_run%02d.csv", .y)), row.names = FALSE)
  )
}

