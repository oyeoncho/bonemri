#####
## Paired t-tests on outer-run (1..30) using inner-run (0..29) summaries
## - Mean test: compare vectors of inner-mean across 30 outer runs (paired by outer run id)
## - SD   test: compare vectors of inner-SD   across 30 outer runs (paired by outer run id)
## - Plots with mean±SD over outer runs + significance stars per time/group
#####

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)

# -------------------------------
# 0) Configs
# -------------------------------
groups_dirs <- paste0("n", 4:7, "_30_30")           # folder names
base_dir    <- "image_data/survival_model/mixture_non_fix/tune/beit0/results/generalization"
model_dirs  <- c(beit0 = "test0", beit0_cox = "test_cox")
dl_dir      <- "dl0"

keep_times     <- c(12, 24, 36, 48, 60, 72)
target_groups  <- paste0("n", 4:7)                  # n4, n5, n6, n7 (표시용)
target_feature <- "Image only"

# -------------------------------
# 1) Safe loader
# -------------------------------
safe_read_csv <- function(path) {
  if (!file.exists(path)) {
    message("❌ not found: ", path)
    return(NULL)
  }
  suppressWarnings(read_csv(path, show_col_types = FALSE))
}

# -------------------------------
# 2) Load all (model → group → outer_run)
#    Keep inner Run (0..29) to summarize inside each outer_run
# -------------------------------
all_outer <- list()

for (model in names(model_dirs)) {
  exp_dir <- model_dirs[[model]]
  for (grp in groups_dirs) {
    # 1..30 : outer run id (파일 이름의 2자리 번호)
    this_list <- map_dfr(1:30, function(outer_id) {
      auc_path    <- file.path(base_dir, exp_dir, dl_dir, grp, sprintf("raw_auc_per_time_run%02d.csv", outer_id))
      cindex_path <- file.path(base_dir, exp_dir, dl_dir, grp, sprintf("raw_cindex_per_time_run%02d.csv", outer_id))
      
      a <- safe_read_csv(auc_path)
      b <- safe_read_csv(cindex_path)
      if (is.null(a) || is.null(b)) return(tibble())
      
      # 표준화
      a <- a %>%
        rename(
          Feature   = `Feature Set`,
          Time      = `Time (Months)`,
          AUC_train = `AUC (Train)`,
          AUC_val   = `AUC (Val)`
        )
      b <- b %>%
        rename(
          Feature        = `Feature Set`,
          Time           = `Time (Months)`,
          C_index_train  = `C-index (Train)`,
          C_index_val    = `C-index (Val)`
        )
      
      # Scope가 있다면 Time-wise만 남기고, 없다면 무시
      if ("Scope" %in% names(a)) a <- a %>% filter(Scope == "Time-wise")
      if ("Scope" %in% names(b)) b <- b %>% filter(Scope == "Time-wise")
      
      # 내부런 Run 열이 반드시 있어야 함
      if (!("Run" %in% names(a)) || !("Run" %in% names(b))) {
        message("⚠️ Missing 'Run' column in: ", dirname(auc_path))
      }
      
      inner_join(
        a %>% select(Feature, Time, Run, AUC_train, AUC_val),
        b %>% select(Feature, Time, Run, C_index_train, C_index_val),
        by = c("Feature", "Time", "Run")
      ) %>%
        mutate(
          outer_run = outer_id,
          group_dir = grp,
          group     = gsub("_30_30$", "", grp),  # n4..n7
          model     = model
        )
    })
    
    if (nrow(this_list) > 0) {
      all_outer[[paste(model, grp, sep = "::")]] <- this_list
    }
  }
}

raw_all <- bind_rows(all_outer)

# -------------------------------
# 3) Filter: Image only, keep times, groups
# -------------------------------
raw_use <- raw_all %>%
  filter(
    Feature == target_feature,
    suppressWarnings(as.numeric(Time)) %in% keep_times
  ) %>%
  mutate(
    Time  = as.numeric(Time),
    group = factor(group, levels = target_groups),
    model = recode(model, "beit0" = "Beit0", "beit0_cox" = "Beit0_cox")
  )

# -------------------------------
# 4) Inner summaries per Outer Run (핵심)
#    For each group×Time×model×outer_run:
#      - inner-mean of AUC_val / C_index_val over inner Run(0..29)
#      - inner-SD   of AUC_val / C_index_val over inner Run(0..29)
# -------------------------------
outer_summ <- raw_use %>%
  group_by(group, Time, model, outer_run) %>%
  summarise(
    mean_auc_val    = mean(AUC_val,     na.rm = TRUE),
    sd_auc_val      = sd(AUC_val,       na.rm = TRUE),
    mean_cindex_val = mean(C_index_val, na.rm = TRUE),
    sd_cindex_val   = sd(C_index_val,   na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------------------------------
# 5) Paired tests across outer_run (30 pairs)
#    Compare Beit0 vs Beit0_cox (paired on outer_run)
#    - mean difference test
#    - sd difference test
# -------------------------------------------------
# helper: paired p-value for one variable name
paired_p <- function(df, value_col) {
  df %>%
    select(group, Time, model, outer_run, !!sym(value_col)) %>%
    pivot_wider(names_from = model, values_from = !!sym(value_col)) %>%
    filter(is.finite(Beit0), is.finite(Beit0_cox)) %>%
    group_by(group, Time) %>%
    summarise(
      n_pairs = n(),
      p_val   = ifelse(n_pairs > 1,
                       t.test(Beit0_cox, Beit0, paired = TRUE)$p.value,
                       NA_real_),
      .groups = "drop"
    )
}

p_to_star <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )
}

# mean 유의성
p_mean_auc  <- paired_p(outer_summ, "mean_auc_val")    %>% mutate(metric = "AUC (Val)",      kind = "Mean")
p_mean_cidx <- paired_p(outer_summ, "mean_cindex_val") %>% mutate(metric = "C-index (Val)",  kind = "Mean")
# sd 유의성
p_sd_auc    <- paired_p(outer_summ, "sd_auc_val")      %>% mutate(metric = "AUC (Val)",      kind = "SD")
p_sd_cidx   <- paired_p(outer_summ, "sd_cindex_val")   %>% mutate(metric = "C-index (Val)",  kind = "SD")

sig_tbl <- bind_rows(p_mean_auc, p_mean_cidx, p_sd_auc, p_sd_cidx) %>%
  mutate(stars = p_to_star(p_val))

# -------------------------------------------------
# 6) Plot data (line = outer_run 평균, ribbon/errorbar = outer_run SD)
#    (a) Mean curves
#    (b) SD curves
#    + stars at Beit0_cox only
# -------------------------------------------------
# (a) mean curves: 먼저 outer_run 평균/표준편차(= across outer runs) 계산
mean_plot_df <- outer_summ %>%
  select(group, Time, model, outer_run, mean_auc_val, mean_cindex_val) %>%
  pivot_longer(cols = c(mean_auc_val, mean_cindex_val),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         "mean_auc_val"    = "AUC (Val)",
                         "mean_cindex_val" = "C-index (Val)")) %>%
  group_by(group, Time, model, metric) %>%
  summarise(
    mean_over_outer = mean(value, na.rm = TRUE),
    sd_over_outer   = sd(value,   na.rm = TRUE),
    .groups = "drop"
  )

# (b) sd curves: 내부 SD의 outer_run 평균/표준편차
sd_plot_df <- outer_summ %>%
  select(group, Time, model, outer_run, sd_auc_val, sd_cindex_val) %>%
  pivot_longer(cols = c(sd_auc_val, sd_cindex_val),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         "sd_auc_val"     = "AUC (Val)",
                         "sd_cindex_val"  = "C-index (Val)")) %>%
  group_by(group, Time, model, metric) %>%
  summarise(
    mean_over_outer = mean(value, na.rm = TRUE),
    sd_over_outer   = sd(value,   na.rm = TRUE),
    .groups = "drop"
  )

# 별표 위치(y) 계산: 각 패널에서 y의 최대치보다 조금 위
ypos_mean <- mean_plot_df %>%
  group_by(group, Time, metric) %>%
  summarise(ypos = max(mean_over_outer + sd_over_outer, na.rm = TRUE), .groups = "drop") %>%
  mutate(ypos = ypos * 1.02)
ypos_sd <- sd_plot_df %>%
  group_by(group, Time, metric) %>%
  summarise(ypos = max(mean_over_outer + sd_over_outer, na.rm = TRUE), .groups = "drop") %>%
  mutate(ypos = ypos * 1.02)

sig_mean <- sig_tbl %>% filter(kind == "Mean") %>%
  left_join(ypos_mean, by = c("group", "Time", "metric"))
sig_sd   <- sig_tbl %>% filter(kind == "SD") %>%
  left_join(ypos_sd,   by = c("group", "Time", "metric"))


# 1) 데이터 결합
mean_plot_df2 <- mean_plot_df %>% mutate(type = "Mean", sig = FALSE)
sd_plot_df2   <- sd_plot_df   %>% mutate(type = "SD",   sig = FALSE)

plot_df <- bind_rows(mean_plot_df2, sd_plot_df2)

# 2) 유의성 데이터 결합
sig_mean2 <- sig_mean %>% mutate(type = "Mean")
sig_sd2   <- sig_sd   %>% mutate(type = "SD")
sig_df    <- bind_rows(sig_mean2, sig_sd2)

# 3) 통합 플롯
ggplot(plot_df,
       aes(x = factor(Time, levels = keep_times),
           y = mean_over_outer,
           color = group,
           linetype = model, shape = model,
           group = interaction(group, model))) +
  geom_line(size = 1.3) +
  geom_point(position = pos, size = 3.5) +
  geom_errorbar(aes(ymin = mean_over_outer - sd_over_outer,
                    ymax = mean_over_outer + sd_over_outer),
                position = pos, width = 0.18, linewidth = 0.9) +
  geom_text(
    data = sig_df %>% mutate(Time = factor(Time, levels = keep_times)),
    aes(x = Time, y = ypos, label = stars, color = group),
    position = pos, vjust = 1,
    fontface = "bold", size =7, show.legend = FALSE
  ) +
  facet_grid(type ~ metric + group, scales = "free_y") +
  scale_linetype_manual(values = linetypes) +
  scale_shape_manual(values = shapes) +
  labs(title = "Mean & Standard deviation (SD) of metrics over 30 runs — with paired t-test",
       x = "Time (Months)", y = "Value (mean ± SD)") +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(face = "bold", size = 20),
    strip.text = element_text(face = "bold", size = 20),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 20),
    axis.text.y = element_text(face = "bold", size = 20),
    axis.title.x = element_text(face = "bold", size = 30),
    axis.title.y = element_text(face = "bold", size = 30),
    plot.title = element_text(face = "bold", size = 35, hjust = 0.5)
  )



