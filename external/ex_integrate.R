# ========= n4~n7 통합: file_id × Feature Set별 mean ±95%CI + pairwise p =========
root_base  <- "image_data/external/TCGA_CCRT/fu_new"
groups     <- paste0("n", 4:7, "_30_30")
subfolder  <- "test1_1"
sets_use   <- c("Clinical only", "Image + Clinical", "Image only")

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(purrr); library(stringr); library(tidyr)
})

# ---------- 유틸 ----------
ci95 <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n <= 1) return(list(n=n, mean=mean(x), sd=sd(x), se=NA_real_, lcl=NA_real_, ucl=NA_real_))
  m  <- mean(x); s <- sd(x); se <- s/sqrt(n); tcrit <- qt(0.975, df=n-1)
  list(n=n, mean=m, sd=s, se=se, lcl=m - tcrit*se, ucl=m + tcrit*se)
}
fs_map <- c("Clinical only" = "clinical",
            "Image + Clinical" = "img_clin",
            "Image only" = "image")

# ---------- 한 그룹 읽기 ----------
load_one_group <- function(grp){
  base_dir <- file.path(root_base, grp, subfolder)
  ids <- sprintf("%02d", 1:30)
  files_cindex <- file.path(base_dir, paste0("file", ids),
                            paste0("external_cindex_ALL_runs_from_file", ids, ".csv"))
  files_auc    <- file.path(base_dir, paste0("file", ids),
                            paste0("external_auc_ALL_runs_from_file",    ids, ".csv"))
  
  exist_cidx <- file.exists(files_cindex)
  exist_auc  <- file.exists(files_auc)
  if (!any(exist_cidx)) stop("C-index CSV가 없습니다. base_dir를 확인하세요: ", base_dir)
  
  read_cidx_one <- function(path, id2) {
    read_csv(path, show_col_types = FALSE) %>%
      rename(runfile=`RunFile`, feature_set=`Feature Set`, time_raw=`Time (Months)`,
             cindex=`C-index (External)`) %>%
      mutate(time=as.character(time_raw), group=grp,
             file_id=paste0("file", id2), feature_set=trimws(feature_set)) %>%
      select(group, file_id, runfile, feature_set, time, cindex)
  }
  read_auc_one <- function(path, id2) {
    read_csv(path, show_col_types = FALSE) %>%
      rename(runfile=`RunFile`, feature_set=`Feature Set`, time_raw=`Time (Months)`,
             auc_raw=`AUC (External)`) %>%
      mutate(time=as.character(time_raw), group=grp,
             auc=suppressWarnings(as.numeric(auc_raw)),
             file_id=paste0("file", id2), feature_set=trimws(feature_set)) %>%
      select(group, file_id, runfile, feature_set, time, auc)
  }
  
  df_cidx <- map2_dfr(files_cindex[exist_cidx], ids[exist_cidx], read_cidx_one)
  if (any(exist_auc)) {
    df_auc <- map2_dfr(files_auc[exist_auc], ids[exist_auc], read_auc_one)
    all_df <- df_cidx %>% left_join(df_auc, by=c("group","file_id","runfile","feature_set","time"))
  } else {
    all_df <- df_cidx %>% mutate(auc = NA_real_)
  }
  
  # Overall만
  all_df %>%
    filter(time == "Overall", feature_set %in% sets_use) %>%
    select(group, file_id, runfile, feature_set, cindex, auc)
}

# ---------- 1) n4~n7 모두 로드 & 통합 ----------
overall_all <- map_dfr(groups, load_one_group)

# ---------- 2) file_id × Feature Set별 요약(통합) ----------
metrics_by_file_all <- overall_all %>%
  group_by(file_id, feature_set) %>%
  summarise(
    # C-index
    {
      ci <- ci95(cindex)
      n_cidx         <- ci$n; mean_cindex <- ci$mean; sd_cindex <- ci$sd
      se_cindex      <- ci$se; lower95_cindex <- ci$lcl; upper95_cindex <- ci$ucl
      # AUC
      cia <- ci95(auc)
      n_auc         <- cia$n; mean_auc <- cia$mean; sd_auc <- cia$sd
      se_auc        <- cia$se; lower95_auc <- cia$lcl; upper95_auc <- cia$ucl
      tibble(n_cidx, mean_cindex, sd_cindex, se_cindex, lower95_cindex, upper95_cindex,
             n_auc, mean_auc, sd_auc, se_auc, lower95_auc, upper95_auc)
    },
    .groups = "drop"
  ) %>%
  arrange(file_id, factor(feature_set, levels = sets_use))

# ---------- 3) file_id 내 Feature Set간 pairwise p-values (통합 데이터로, run 매칭) ----------
pair_list  <- list(c("clinical","image"),
                   c("img_clin","image"),
                   c("img_clin","clinical"))
pair_names <- c("Clinical vs Image",
                "Image+Clinical vs Image",
                "Image+Clinical vs Clinical")

pairwise_fun <- function(dat, value_col, metric_label) {
  # 같은 file_id 안에서 group+runfile을 매칭 키로, 세트별 wide 구성
  w <- dat %>%
    mutate(fs_short = recode(feature_set, !!!fs_map)) %>%
    select(group, runfile, fs_short, value = {{ value_col }}) %>%
    pivot_wider(names_from = fs_short, values_from = value)
  
  purrr::imap_dfr(pair_list, function(p, i) {
    # 필요한 세트가 모두 있어야 함
    if (!all(p %in% names(w))) {
      tibble(pair_label = pair_names[i], metric = metric_label,
             x_set = p[1], y_set = p[2], n_pairs = 0,
             mean_x = NA_real_, mean_y = NA_real_, mean_diff = NA_real_,
             t_stat = NA_real_, p_value = NA_real_)
    } else {
      dfp <- w %>% select(all_of(p))
      ok  <- stats::complete.cases(dfp)
      x   <- dfp[[p[1]]][ok]; y <- dfp[[p[2]]][ok]
      n   <- length(x)
      if (n >= 2 && sd(x - y, na.rm = TRUE) > 0) {
        tt <- t.test(x, y, paired = TRUE)
        tibble(pair_label = pair_names[i], metric = metric_label,
               x_set = p[1], y_set = p[2], n_pairs = n,
               mean_x = mean(x), mean_y = mean(y),
               mean_diff = mean(x - y), t_stat = unname(tt$statistic),
               p_value = unname(tt$p.value))
      } else {
        tibble(pair_label = pair_names[i], metric = metric_label,
               x_set = p[1], y_set = p[2], n_pairs = n,
               mean_x = if (n>0) mean(x) else NA_real_,
               mean_y = if (n>0) mean(y) else NA_real_,
               mean_diff = if (n>0) mean(x - y) else NA_real_,
               t_stat = NA_real_, p_value = NA_real_)
      }
    }
  })
}

pairwise_pvalues_by_file_all <- overall_all %>%
  group_by(file_id) %>%
  group_modify(~{
    dat <- .x
    bind_rows(
      pairwise_fun(dat, cindex, "cindex"),
      pairwise_fun(dat, auc,    "auc")
    ) %>%
      mutate(file_id = unique(dat$file_id))
  }) %>%
  ungroup() %>%
  group_by(file_id, metric) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  select(file_id, metric, pair_label, x_set, y_set,
         n_pairs, mean_x, mean_y, mean_diff, t_stat, p_value, p_adj_BH)

# ---------- 4) 저장 ----------
out_metrics <- file.path(root_base, "metrics_by_file_all.csv")
out_pvals   <- file.path(root_base, "pairwise_pvalues_by_file_all.csv")
write_csv(metrics_by_file_all, out_metrics)
write_csv(pairwise_pvalues_by_file_all, out_pvals)
cat("✅ saved:\n  - ", out_metrics, "\n  - ", out_pvals, "\n", sep = "")


# ========= 통합 리포트 생성 =========
library(dplyr)
library(tidyr)
library(readr)

root_base  <- "image_data/external/TCGA_CCRT/fu_new"

# 기존에 저장된 CSV 불러오기
metrics_by_file_all <- read_csv(file.path(root_base, "metrics_by_file_all.csv"),
                                show_col_types = FALSE)
pairwise_pvalues_by_file_all <- read_csv(file.path(root_base, "pairwise_pvalues_by_file_all.csv"),
                                         show_col_types = FALSE)

# ---- 1) p-value 테이블을 wide로 변환 ----
pvals_wide <- pairwise_pvalues_by_file_all %>%
  mutate(metric_pair = paste(metric, pair_label, sep = "_")) %>%
  select(file_id, metric_pair, p_value, p_adj_BH) %>%
  pivot_wider(
    names_from = metric_pair,
    values_from = c(p_value, p_adj_BH),
    names_sep = "."
  )

# ---- 2) metrics와 p-value 병합 ----
report_all <- metrics_by_file_all %>%
  left_join(pvals_wide, by = "file_id") %>%
  arrange(file_id, feature_set)

# ---- 3) 저장 ----
out_report <- file.path(root_base, "report_metrics_pvalues_by_file.csv")
write_csv(report_all, out_report)

cat("✅ saved:", out_report, "\n")
