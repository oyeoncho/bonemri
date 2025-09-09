library(dplyr)
library(readr)
library(forcats)
library(ggplot2)
# 0) 임상 데이터 1번만 로드
a <- read.csv("image_data/external/clinical1.csv", stringsAsFactors = FALSE)

# 공통 피처 목록
feat_cols <- c("feat_213", "feat_194", "feat_163", "feat_266", "feat_407", "feat_468", "feat_499", "feat_169", "feat_436", "feat_560", 
               "feat_2", "feat_327", "feat_391", "feat_519", "feat_173", "feat_181", "feat_389", "feat_715", "feat_107", "feat_203", 
               "feat_361", "feat_439", "feat_451", "feat_565", "feat_747", "feat_80", "feat_10", "feat_123", "feat_137", "feat_15", 
               "feat_209", "feat_215", "feat_289", "feat_368", "feat_374", "feat_55", "feat_576", "feat_578", "feat_121", "feat_125", 
               "feat_143", "feat_223", "feat_240", "feat_25", "feat_309", "feat_498", "feat_514", "feat_554", "feat_577", "feat_617", 
               "feat_653", "feat_710", "feat_109", "feat_210", "feat_220", "feat_352", "feat_420", "feat_507", "feat_583", "feat_605", 
               "feat_657", "feat_666", "feat_152", "feat_167", "feat_255", "feat_328", "feat_378", "feat_402", "feat_633", "feat_656")

# 1) 한 개 인덱스를 처리하는 함수
make_external <- function(k) {
  f_in  <- sprintf("image_data/external/beit_base_patch16_224_features_%d.csv", k)
  f_out <- sprintf("image_data/external/external%d.csv", k)
  
  if (!file.exists(f_in)) {
    message("⚠️ 파일 없음: ", f_in, " → 스킵")
    return(invisible(NULL))
  }
  
  b <- read.csv(f_in, stringsAsFactors = FALSE) %>%
    mutate(PatientID = sub("_0000$", "", PatientID)) %>%
    distinct(PatientID, .keep_all = TRUE)
  
  df <- a %>% inner_join(b, by = "PatientID")
  
  out <- df %>% select(any_of(c("PatientID", "survival", "fu_date", "Age",	"pathology",	"stage0", feat_cols)))
  
  write_csv(out, f_out)
  message("✅ 저장: ", f_out, " (n=", nrow(out), ")")
}

# 2) 3 ~ 30 반복 실행
for (k in 1:30) make_external(k)

