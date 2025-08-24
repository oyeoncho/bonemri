library(readxl)
library(dplyr)
library(xlsx)
library(tidyverse)
library(corrplot)
library(Hmisc)
library(moonBook)
library(Hmisc)  # rcorr 사용

cx <- read.csv(file = "cx.csv") %>% as_tibble() %>% rename("ID1"="PatientID")



N_RUN <- 29

#전체 데이터 불러오기
beit0i <- read.csv("image_data/image_features/tuning/final_features/slice_20_dl0/beit0_backbone_features_merged.csv") %>%
  as_tibble()

cx_pre %>% filter(fraction_size < 4 | is.na(fraction_size)==T) %>% select(ICR, fraction_size)

# 1️⃣ 전체 수치형 피처만 추출해서 공통 스케일 기준 계산
numeric_cols <- beit0i %>% select(where(is.numeric)) %>% select(-Run)  # Run 제외
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
numeric_cols <- beit0i %>% select(where(is.numeric)) %>% select(-Run)  # Run 제외
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
  beit0 <- beit0i %>% filter(Run == k)
  
  # 수치형 피처만 추출하고 공통 스케일 적용
  numeric_feats <- beit0 %>% select(where(is.numeric)) %>% select(-Run)
  numeric_feats_scaled <- scale_global(numeric_feats)
  
  # PatientID 유지
  numeric_feats_scaled <- numeric_feats_scaled %>%
    mutate(PatientID = beit0$PatientID)
  
  # pre 추출
  pre <- numeric_feats_scaled %>%
    filter(str_detect(PatientID, "pre")) %>%
    rename("ID1" = "PatientID") %>%
    mutate(ID = ifelse(substr(ID1, 1, 1) == "P", substr(ID1, 1, 4), substr(ID1, 1, 5))) %>%
    select(ID, ID1, where(is.numeric))
  
  # 병합
  dh1 <- cx %>%
    inner_join(pre, by = "ID1") %>%
    select(-ID.x, -ID.y, -ID2)
  
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
library(purrr)
library(tidyr)
library(ggplot2)

# 1) 합의집합 정의
consensus_sets <- list(
  n4 = c("feat_213","feat_194","feat_163","feat_266","feat_407","feat_468","feat_499","feat_169","feat_436","feat_560",
         "feat_2","feat_327","feat_391","feat_519","feat_173","feat_181","feat_389","feat_715","feat_107","feat_203",
         "feat_361","feat_439","feat_451","feat_565","feat_747","feat_80","feat_10","feat_123","feat_137","feat_15",
         "feat_209","feat_215","feat_289","feat_368","feat_374","feat_55","feat_576","feat_578","feat_121","feat_125",
         "feat_143","feat_223","feat_240","feat_25","feat_309","feat_498","feat_514","feat_554","feat_577","feat_617",
         "feat_653","feat_710","feat_109","feat_210","feat_220","feat_352","feat_420","feat_507","feat_583","feat_605",
         "feat_657","feat_666","feat_152","feat_167","feat_255","feat_328","feat_378","feat_402","feat_633","feat_656"),
  n5 = c("feat_213","feat_266","feat_499","feat_436","feat_2","feat_327","feat_391","feat_519","feat_173","feat_715",
         "feat_107","feat_80","feat_137","feat_209","feat_215","feat_374","feat_55","feat_223","feat_554","feat_577",
         "feat_109","feat_583","feat_657"),
  n6 = c("feat_213","feat_266","feat_436","feat_519","feat_215"),
  n7 = c("feat_436","feat_519")
)

# 2) Jaccard 함수
jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  u <- union(a, b); i <- intersect(a, b)
  if (length(u) == 0) return(NA_real_)
  length(i) / length(u)
}

# 3) run별 Top100 리스트 만들기
runs_df <- b_all %>%
  group_by(Run) %>%
  summarise(feats = list(Name), .groups = "drop")

# 4) 합의집합별 Jaccard 계산
jac_by_run <- runs_df %>%
  mutate(
    J_n4 = map_dbl(feats, ~ jaccard(.x, consensus_sets$n4)),
    J_n5 = map_dbl(feats, ~ jaccard(.x, consensus_sets$n5)),
    J_n6 = map_dbl(feats, ~ jaccard(.x, consensus_sets$n6)),
    J_n7 = map_dbl(feats, ~ jaccard(.x, consensus_sets$n7))
  )

# 5) 대표 run 후보 정렬 (보수적 합의집합 n5 우선, 동률 시 n4→n6→n7)
selection_table <- jac_by_run %>%
  arrange(desc(J_n5), desc(J_n4), desc(J_n6), desc(J_n7)) %>%
  select(Run, starts_with("J_"))

selection_table %>% arrange(desc(J_n4+J_n5+J_n6+J_n7))

