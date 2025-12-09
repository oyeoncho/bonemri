library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# ---- 경로 설정 ----
tsv_path <- "image_data/external/clinical.tsv"   # GDC clinical TSV
csv_path <- "image_data/external/clinical_origin.csv"   # 기존 clinical.csv

# ---- ID 표준화 ----
normalize_id <- function(x) {
  x %>% trimws() %>% toupper() %>% 
    str_replace("_0000$", "") %>% str_replace_all("_", "-")
}

# ---- 1) 기존 clinical.csv 로드 + ID 표준화 ----
a <- read.csv(csv_path, stringsAsFactors = FALSE) %>%
  mutate(PatientID = normalize_id(PatientID))

# ---- 2) clinical.tsv 로드 ----
clin_tsv <- read_tsv(tsv_path, show_col_types = FALSE, progress = FALSE)

# ---- 3) 필요한 컬럼만 선별 ----
req_cols <- c(
  "cases.submitter_id",
  "cases.disease_type",
  "diagnoses.figo_stage",
  "treatments.days_to_treatment_end",
  "diagnoses.figo_staging_edition_year",
  "diagnoses.diagnosis_is_primary_disease",
  "diagnoses.synchronous_malignancy",
  "treatments.chemo_concurrent_to_radiation",
  "treatments.number_of_fractions",
  "treatments.reason_treatment_not_given",
  "treatments.therapeutic_agents",
  "treatments.timepoint_category",
  "treatments.treatment_anatomic_sites",
  "treatments.treatment_dose",
  "treatments.treatment_intent_type",
  "treatments.treatment_type",
  "cases.lost_to_followup",
  "demographic.vital_status",
  "diagnoses.classification_of_tumor",
  "diagnoses.days_to_last_follow_up"
)

present_cols <- intersect(req_cols, names(clin_tsv))
missing_cols <- setdiff(req_cols, present_cols)
if (length(missing_cols) > 0) {
  message("⚠️ clinical.tsv에 없는 컬럼: ", paste(missing_cols, collapse = ", "))
}

clin_sel <- clin_tsv %>%
  select(any_of(present_cols)) %>%
  rename(submitter_id = `cases.submitter_id`) %>%
  mutate(PatientID = normalize_id(submitter_id)) %>%
  select(-submitter_id)

# ---- 4) 환자 단위로 '; ' 병합 (원자료 유지) ----
collapse_semicolon <- function(x) {
  u <- unique(na.omit(as.character(x)))
  if (length(u) == 0) NA_character_ else paste(u, collapse = "; ")
}
clin_agg <- clin_sel %>%
  group_by(PatientID) %>%
  summarise(across(everything(), collapse_semicolon), .groups = "drop")

# ---- 5) 조인 → TCGA-VS/CCTH만 선별 (원본 문자열은 그대로 유지) ----
joined <- a %>% left_join(clin_agg, by = "PatientID")

joined_sel <- joined %>% 
  filter(!is.na(PatientID)) %>%
  filter(str_detect(PatientID, "^(TCGA-VS|CCTH)"))

# ─────────────────────────────────────────────────────────────────────
# ⬇ 문자열 정리(클린 컬럼) 없이, 원본에서 바로 파생 플래그 계산
# ─────────────────────────────────────────────────────────────────────

# 6) treatment_type 기반 플래그: brachytherapy / hysterectomy
tt_long <- joined_sel %>%
  select(PatientID, tt = `treatments.treatment_type`) %>%
  mutate(tt = coalesce(tt, "")) %>%
  separate_rows(tt, sep = ";") %>%
  mutate(tt = str_trim(tt)) %>%
  filter(tt != "", !tt %in% c("'--","--","Not Reported","Unknown","UNKNOWN")) %>%
  mutate(tt_l = str_to_lower(tt))

tt_flags <- tt_long %>%
  mutate(
    # brachytherapy 키워드: brachy / HDR / high-dose(rate)
    is_brachy = str_detect(tt_l, "brachy") |
      str_detect(tt_l, "(^|\\s)hdr(\\s|$)") |
      str_detect(tt_l, "high[- ]dose( rate)?"),
    
    # hysterectomy 키워드: hysterect*, 구체 표현, 흔한 약어
    is_hyst = str_detect(tt_l, "hysterect") |                         # hysterectomy, hysterectomies...
      str_detect(tt_l, "radical\\s+hyst") |                   # radical hysterectomy
      str_detect(tt_l, "total\\s+hyst") |                     # total hysterectomy
      str_detect(tt_l, "simple\\s+hyst|subtotal\\s+hyst") |   # simple/subtotal
      str_detect(tt_l, "\\b(tah|tlh|lavh|vh|lrh|rrh|tlrh|tvrh)\\b") # 약어들
  ) %>%
  group_by(PatientID) %>%
  summarise(
    brachytherapy = as.integer(any(is_brachy, na.rm = TRUE)),
    hysterectomy  = as.integer(any(is_hyst,   na.rm = TRUE)),
    .groups = "drop"
  )


# 7) intent_adjuvant 플래그: treatment_intent_type에 Adjuvant 포함 시 1
intent_long <- joined_sel %>%
  select(PatientID, intent = `treatments.treatment_intent_type`) %>%
  mutate(intent = coalesce(intent, "")) %>%
  separate_rows(intent, sep = ";") %>%
  mutate(intent = str_trim(intent)) %>%
  filter(intent != "", !intent %in% c("'--","--","Not Reported","Unknown","UNKNOWN")) %>%
  mutate(intent_l = str_to_lower(intent))

intent_flag <- intent_long %>%
  group_by(PatientID) %>%
  summarise(intent_adjuvant = as.integer(any(str_detect(intent_l, "\\badjuvant\\b"), na.rm = TRUE)),
            .groups = "drop")

# 8) chemo_concurrent 플래그: Yes/True/Y/1 → 1, No/False/N/0 → 0, 그 외 NA
cc_long <- joined_sel %>%
  select(PatientID, cc = `treatments.chemo_concurrent_to_radiation`) %>%
  mutate(cc = coalesce(cc, "")) %>%
  separate_rows(cc, sep = ";") %>%
  mutate(cc = str_trim(cc),
         cc_l = str_to_lower(cc)) %>%
  filter(cc != "", !cc %in% c("'--","--","Not Reported","Unknown","UNKNOWN"))

cc_flag <- cc_long %>%
  group_by(PatientID) %>%
  summarise(
    chemo_concurrent = case_when(
      any(cc_l %in% c("yes","y","true","t","1")) ~ 1L,
      any(cc_l %in% c("no","n","false","f","0")) ~ 0L,
      TRUE ~ NA_integer_
    ),
    .groups = "drop"
  )

# 9) 플래그들 결합 (원본 문자열 컬럼은 그대로 유지)
joined_flags <- joined_sel %>%
  left_join(tt_flags,   by = "PatientID") %>%   # ← brachy + hyst 둘 다 포함
  left_join(intent_flag, by = "PatientID") %>%
  left_join(cc_flag,     by = "PatientID") %>%
  mutate(
    brachytherapy   = ifelse(is.na(brachytherapy), 0L, brachytherapy),
    hysterectomy    = ifelse(is.na(hysterectomy),  0L, hysterectomy),  # ← 추가
    intent_adjuvant = ifelse(is.na(intent_adjuvant), 0L, intent_adjuvant)
    # chemo_concurrent 는 기록 없으면 NA 유지
  )

joined_flags <-joined_flags %>% 
  filter(signal =="SE" | signal == "FSE")

# ─────────────────────────────────────────────────────────────
# NEW) days_to_treatment_end: 42일 초과값 중 '최소값' 선택 → 결측은 중앙값으로 대체
#      fu_date를 치료종료일 기준으로 업데이트 (진단기준에서 d/30 '빼기')
# ─────────────────────────────────────────────────────────────
library(purrr)


# 세미콜론 문자열 → 숫자 벡터
extract_nums <- function(x) {
  toks <- unlist(str_split(coalesce(as.character(x), ""), ";"))
  toks <- str_trim(toks)
  toks <- toks[!(toks %in% c("'--","--","Not Reported","Unknown","UNKNOWN","", NA))]
  as.numeric(toks[grepl("^-?\\d+(\\.\\d+)?$", toks)])
}

LOW  <- 42
HIGH <- 120  # (42, 120) 사이 값만 사용

# 환자별: (42일 초과 & 120일 미만) 값 중 '최소값'
dte_tbl <- joined_flags %>%
  transmute(
    PatientID,
    dte_vals = map(`treatments.days_to_treatment_end`, extract_nums),
    dte_days_raw = map_dbl(dte_vals, ~{
      v <- .x[.x > LOW & .x < HIGH]
      if (length(v) > 0) min(v) else NA_real_
    })
  )

# 결측을 전체 중앙값으로 대체
median_dte <- median(dte_tbl$dte_days_raw, na.rm = TRUE)
# (옵션) 모두 NA인 경우 대비 기본값 설정
if (is.na(median_dte)) median_dte <- 90

dte_tbl <- dte_tbl %>%
  mutate(
    dte_days_imputed   = ifelse(is.na(dte_days_raw), median_dte, dte_days_raw),
    dte_months_imputed = dte_days_imputed / 30
  )

# fu_date 업데이트: 진단기준 fu에서 (치료종료일/월) 감산 → post-treatment 추적
joined_flags <- joined_flags %>%
  mutate(fu_date = as.numeric(fu_date)) %>%
  left_join(dte_tbl %>% select(PatientID, dte_days_imputed, dte_months_imputed),
            by = "PatientID") %>%
  mutate(
    fu_date_orig = fu_date,
    fu_date      = ifelse(!is.na(dte_months_imputed), fu_date - dte_months_imputed, fu_date),
    fu_date      = pmax(fu_date, 0)
  )


# (선택) QC 저장
readr::write_csv(
  dte_tbl %>% mutate(median_used = median_dte),
  file.path(dirname(csv_path), "dte_selection_qc_after_join_42_120.csv")
)

# (필요 시) 최종 저장
out_path <- "image_data/external/clinical_joined_from_tsv_flags.csv"
readr::write_csv(joined_flags, out_path)
cat("✅ fu_date updated AFTER JOIN & saved:", out_path, "\n")
