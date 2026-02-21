# install.packages(c("survival","survminer","ggplot2","dplyr","readr","stringr","broom","forestplot"))
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

OUT_DIR <- "./image_data/survival_model/mixture_non_fix/non_nest/beit0/results/reviewer_km_pack/dl0/n7_30_30/file04_run06"

df_stage <- read_csv(file.path(OUT_DIR, "DATA_INTERNAL_ALL_stagebin_long.csv"), show_col_types = FALSE)

df_stage <- df_stage %>%
  mutate(
    time  = as.numeric(time),
    event = as.numeric(event),
    group = factor(group, levels = c("IB–IIIC1 (stage0 1–2)", "IIIC2–IVB (stage0 3–4)"))
  ) %>%
  filter(!is.na(time), !is.na(event), !is.na(group))

fit_stage <- survfit(Surv(time, event) ~ group, data = df_stage)

p_stage <- ggsurvplot(
  fit_stage, data = df_stage,
  conf.int = TRUE,          # CI on/off
  risk.table = TRUE,        # risk table on/off
  pval = TRUE,              # log-rank p 표시
  pval.coord = c(10, 0.15), # p-value 위치 (x,y) 필요시 조정
  break.time.by = 12,       # x축 눈금 간격
  xlim = c(0, 84),         # x축 범위 (원하면 수정)
  ylim = c(0, 1),           # y축 범위
  legend.title = "",
  legend.labs = c("Stage IB–IIIC1", "Stage IIIC2–IVB"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
  surv.median.line = NULL   # 중앙생존선 (원치 않으면 NULL)
)

# ggplot 스타일 수정(폰트/라인/여백 등)
p_stage$plot <- p_stage$plot +
  theme_classic(base_size = 13) +
  theme(
    legend.position = c(0.75, 1),
    legend.background = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(face="bold")
  ) +
  ggtitle("FIGO 2018 Stage (IB–IIIC1 vs IIIC2–IVB)")

# 저장
p_stage

#################
df_risk <- read_csv(file.path(OUT_DIR, "DATA_INTERNAL_ALL_imagingRisk_T60_withGroup.csv"), show_col_types = FALSE)

df_risk <- df_risk %>%
  mutate(
    time  = as.numeric(time),
    event = as.numeric(event),
    risk  = as.numeric(risk),
    risk_group = factor(risk_group, levels=c("Low","High"))
  ) %>%
  filter(!is.na(time), !is.na(event), !is.na(risk_group))

cutoff_val <- unique(df_risk$cutoff)[1]

fit_risk <- survfit(Surv(time, event) ~ risk_group, data=df_risk)

p_risk <- ggsurvplot(
  fit_risk, data=df_risk,
  conf.int = TRUE,
  risk.table = TRUE,
  pval = TRUE,
  pval.coord = c(10, 0.15),
  break.time.by = 12,
  xlim = c(0, 84),
  ylim = c(0, 1),
  legend.title = "",
  legend.labs = c("Low imaging risk", "High imaging risk"),
  xlab = "Time (months)",
  ylab = "Overall survival probability"
)

p_risk$plot <- p_risk$plot +
  theme_classic(base_size = 13) +
  ggtitle("Imaging risk (High vs Low)") +
  labs(subtitle = paste0("Cutoff (internal median) = ", signif(cutoff_val, 4))) +
  theme(
    legend.position = c(0.75, 0.3),
    plot.title = element_text(face="bold")
  )

p_risk

###
library(readr)
library(dplyr)
library(stringr)
library(survival)
library(survminer)
library(broom)
library(ggplot2)
library(ggpubr)
library(ggrepel)

# ----------------------------
# 0) Load
# ----------------------------
df4 <- read_csv(
  file.path(OUT_DIR, "DATA_INTERNAL_ALL_stagebin_x_imagingRisk_T60_long.csv"),
  show_col_types = FALSE
)

# ----------------------------
# 1) Group label cleanup
# ----------------------------
df4$group <- gsub("IB–IIIC1 \\(stage0 1–2\\) \\| LowRisk",   "IB-IIIC1 + LowRisk",   df4$group)
df4$group <- gsub("IB–IIIC1 \\(stage0 1–2\\) \\| HighRisk",  "IB-IIIC1 + HighRisk",  df4$group)
df4$group <- gsub("IIIC2–IVB \\(stage0 3–4\\) \\| LowRisk",  "IIIC2-IVB + LowRisk",  df4$group)
df4$group <- gsub("IIIC2–IVB \\(stage0 3–4\\) \\| HighRisk", "IIIC2-IVB + HighRisk", df4$group)

df4$group <- str_replace_all(df4$group, "[\u2012\u2013\u2014\u2212]", "-")
df4$group <- str_replace_all(df4$group, "[\u00A0\u2007\u202F]", " ")
df4$group <- str_squish(df4$group)

# ----------------------------
# 2) Factor levels
# ----------------------------
levels4 <- c(
  "IB-IIIC1 + LowRisk",
  "IB-IIIC1 + HighRisk",
  "IIIC2-IVB + LowRisk",
  "IIIC2-IVB + HighRisk"
)

df4 <- df4 %>%
  mutate(
    time  = as.numeric(time),
    event = as.numeric(event),
    group = factor(group, levels = levels4)
  ) %>%
  filter(!is.na(time), !is.na(event), !is.na(group))

print(table(df4$group))

# ----------------------------
# 3) KM fit
# ----------------------------
fit4 <- survfit(Surv(time, event) ~ group, data = df4)
summary(fit4, time=60)
# ----------------------------
# 4) Palette (levels4 순서와 매칭)
# ----------------------------
pal4 <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
names(pal4) <- paste0("group=", levels4)

# ----------------------------
# 5) Main KM (legend 제거)
# ----------------------------
p_main <- ggsurvplot(
  fit4, data = df4,
  conf.int = FALSE,
  risk.table = FALSE,
  pval = FALSE,
  break.time.by = 12,
  xlim = c(0, 84),
  ylim = c(0, 1),
  legend = "none",
  palette = pal4,
  xlab = "Time (months)",
  ylab = "Overall survival probability"
)

# ----------------------------
# 6) Curve labels (겹침 해결: ggrepel)
# ----------------------------
km_df <- broom::tidy(fit4)
t_label <- 60

label_df <- km_df %>%
  group_by(strata) %>%
  filter(time <= t_label) %>%
  slice_max(time, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    label = gsub("^group=", "", strata),
    x = t_label + 2,
    y = estimate
  )

p_main$plot <- p_main$plot +
  ggtitle("Stage × Imaging risk (4 groups)") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face="bold"),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(0, 84), breaks = seq(0, 84, 12)) +
  guides(colour = "none", fill = "none") +
  ggrepel::geom_label_repel(
    data = label_df,
    aes(x = x, y = y, label = label, colour = strata),
    inherit.aes = FALSE,
    direction = "y",        # ✅ y방향으로만 밀어서 top 라벨 겹침 해결
    nudge_x = 0,            # x는 고정 (원하면 +0.5 정도 가능)
    box.padding = 0.25,
    point.padding = 0.1,
    min.segment.length = 0,
    segment.size = 0.25,
    size = 4,
    fontface = "bold",
    fill = "white",
    alpha = 0.9,
    show.legend = FALSE     # ✅ 범례 다시 안 생기게
  )

# ----------------------------
# 7) Risk table 직접 생성 (정확 매칭 유지)
# ----------------------------
times_tbl <- seq(0, 84, 12)
sfit <- summary(fit4, times = times_tbl)

risk_df <- data.frame(
  strata = sfit$strata,
  time   = sfit$time,
  n_risk = sfit$n.risk
)

risk_df$strata <- factor(risk_df$strata, levels = paste0("group=", levels4))

p_tbl <- ggplot(risk_df, aes(x = time, y = strata, label = n_risk, colour = strata)) +
  geom_text(size = 5, show.legend = FALSE) +  # ✅ 범례 제거
  scale_colour_manual(values = pal4, guide = "none") +
  scale_x_continuous(limits = c(0, 84), breaks = times_tbl) +
  scale_y_discrete(labels = function(x) gsub("^group=", "", x)) +
  labs(
    title = "Number at risk",
    x = "Time (months)",
    y = "Strata"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face="bold", hjust = 0.5, size = 22),
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14)
  )

# ----------------------------
# 8) Combine
# ----------------------------
ggarrange(
  p_main$plot,
  p_tbl,
  ncol = 1,
  heights = c(3.2, 1.3),
  align = "v"
)

library(dplyr)
library(survival)
library(survminer)
library(tidyr)

# ---------------------------------
# 1) Pairwise log-rank test
# ---------------------------------
pw <- pairwise_survdiff(
  Surv(time, event) ~ group,
  data = df4,
  p.adjust.method = "none"   # 먼저 raw
)

# ---------------------------------
# 2) Long format 변환
# ---------------------------------
pw_long <- as.data.frame(as.table(pw$p.value)) %>%
  rename(group1 = Var1,
         group2 = Var2,
         p_raw  = Freq) %>%
  filter(!is.na(p_raw)) %>%
  mutate(
    group1 = gsub("^group=", "", group1),
    group2 = gsub("^group=", "", group2)
  )

# ---------------------------------
# 3) Multiple comparison correction
# ---------------------------------
pw_long <- pw_long %>%
  mutate(
    p_bonf = p.adjust(p_raw, method = "bonferroni"),
    p_fdr  = p.adjust(p_raw, method = "fdr")
  )

# ---------------------------------
# 4) 보기 좋게 formatting
# ---------------------------------
format_p <- function(p){
  ifelse(p < 0.001, "<0.001",
         formatC(p, format="f", digits=3))
}

pw_table <- pw_long %>%
  mutate(
    p_raw_txt  = format_p(p_raw),
    p_bonf_txt = format_p(p_bonf),
    p_fdr_txt  = format_p(p_fdr)
  ) %>%
  arrange(p_raw) %>%   # ✅ 먼저 정렬
  select(group1, group2, p_raw_txt, p_bonf_txt, p_fdr_txt)

pw_table
##

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

# ============================================================
# CONFIG
# ============================================================
OUT_DIR <- "./image_data/survival_model/mixture_non_fix/non_nest/beit0/results/reviewer_km_pack/dl0/n7_30_30/file04_run06"

uv_path <- file.path(OUT_DIR, "Cox_INTERNAL_ALL_T60_stageBIN_imgrisk_per0p1_UV.csv")
mv_path <- file.path(OUT_DIR, "Cox_INTERNAL_ALL_T60_stageBIN_imgrisk_per0p1_MV.csv")

uv <- read_csv(uv_path, show_col_types = FALSE)
mv <- read_csv(mv_path, show_col_types = FALSE)

# ============================================================
# 1) Loader -> standardize columns: term, HR, lo, hi, p, model
# ============================================================
prep_cox_for_forest <- function(df, model_label = "MV") {
  d <- df
  
  # 케이스 1) UV에서 우리가 만든 컬럼(variable, HR, HR_95low, HR_95high, p)
  if (all(c("variable","HR","HR_95low","HR_95high","p") %in% names(d))) {
    out <- d %>%
      transmute(
        term = as.character(variable),
        HR   = as.numeric(HR),
        lo   = as.numeric(HR_95low),
        hi   = as.numeric(HR_95high),
        p    = as.numeric(p)
      )
    out$model <- model_label
    return(out)
  }
  
  # 케이스 2) lifelines summary 그대로 저장된 MV (열 이름이 길 수 있음)
  term_col <- if ("variable" %in% names(d)) "variable" else if ("covariate" %in% names(d)) "covariate" else names(d)[1]
  
  p_col  <- intersect(names(d), c("p", "p_value", "p-value"))[1]
  hr_col <- intersect(names(d), c("exp(coef)", "HR", "exp_coef"))[1]
  lo_col <- intersect(names(d), c("exp(coef) lower 95%", "HR_95low", "lower 0.95", "lower 95%"))[1]
  hi_col <- intersect(names(d), c("exp(coef) upper 95%", "HR_95high", "upper 0.95", "upper 95%"))[1]
  
  out <- d %>%
    transmute(
      term = as.character(.data[[term_col]]),
      HR   = as.numeric(.data[[hr_col]]),
      lo   = as.numeric(.data[[lo_col]]),
      hi   = as.numeric(.data[[hi_col]]),
      p    = if (!is.na(p_col) && !is.null(p_col)) as.numeric(.data[[p_col]]) else NA_real_
    )
  out$model <- model_label
  out
}

# ============================================================
# 2) AGE scaling: per 10-year increase
#    HR_10 = HR_1^10; CI also ^10
# ============================================================
scale_age_to_10yr <- function(df_forest, age_term = "Age", years = 10) {
  df_forest %>%
    mutate(
      HR = ifelse(term == age_term, HR^years, HR),
      lo = ifelse(term == age_term, lo^years, lo),
      hi = ifelse(term == age_term, hi^years, hi)
    )
}

# ============================================================
# 3) Pretty labels
# ============================================================
pretty_term <- function(term) {
  term %>%
    str_replace("^Age$", "Age (per 10-year increase)") %>%
    str_replace("^stage_bin$", "Stage (IIIC2–IVB)") %>%
    str_replace("^imaging_risk$", "Imaging risk (per 10% increase)") %>%
    str_replace("^pathology_1$", "Pathology (non-SqCC)") %>%
    str_replace("_", " ")
}

# ============================================================
# 4) Single forest plot
# ============================================================
plot_forest <- function(df_forest,
                        title = "Cox proportional hazards model",
                        out_png = NULL) {
  
  d <- df_forest %>%
    filter(is.finite(HR), is.finite(lo), is.finite(hi)) %>%
    mutate(
      term_pretty = pretty_term(as.character(term)),
      p_txt = case_when(
        is.na(p) ~ "",
        p < 0.001 ~ "p < 0.001",
        TRUE ~ paste0("p = ", formatC(p, format="f", digits=3))
      ),
      hr_txt = paste0(
        formatC(HR, format="f", digits=2),
        " (", formatC(lo, format="f", digits=2), "–", formatC(hi, format="f", digits=2), ")"
      ),
      label_right = ifelse(p_txt == "", hr_txt, paste0(hr_txt, ", ", p_txt))
    ) %>%
    arrange(desc(HR)) %>%
    mutate(term_pretty = factor(term_pretty, levels = rev(unique(term_pretty))))
  
  xmin <- min(d$lo, na.rm = TRUE)
  xmax <- max(d$hi, na.rm = TRUE)
  
  p <- ggplot(d, aes(x = HR, y = term_pretty)) +
    geom_vline(xintercept = 1, linetype = 2) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15) +
    geom_point(size = 2.2) +
    geom_text(aes(x = xmax * 1.15, label = label_right),
              hjust = 0, size = 3.3) +
    scale_x_log10(limits = c(xmin * 0.8, xmax * 2.2)) +
    labs(title = title, x = "Hazard ratio (log scale)", y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 11),
      plot.margin = margin(5.5, 90, 5.5, 5.5)
    )
  
  if (!is.null(out_png)) {
    ggsave(out_png, p, width = 8.5, height = 4.8, dpi = 300)
  }
  p
}

# ============================================================
# 5) UV / MV build + Age 10y scaling
# ============================================================
uv_f <- prep_cox_for_forest(uv, "UV") %>% scale_age_to_10yr(age_term = "Age", years = 10)
mv_f <- prep_cox_for_forest(mv, "MV") %>% scale_age_to_10yr(age_term = "Age", years = 10)

# ============================================================
# 6) Plot: UV, MV separately
# ============================================================
p_uv <- plot_forest(
  uv_f,
  title = "Univariable Cox (Internal) — Age per 10-year increase",
  out_png = file.path(OUT_DIR, "FOREST_UV_stageBIN_Age10.png")
)

p_mv <- plot_forest(
  mv_f,
  title = "Multivariable Cox (Internal) — Age per 10-year increase",
  out_png = file.path(OUT_DIR, "FOREST_MV_stageBIN_Age10.png")
)


# ============================================================
# 7) Combined UV vs MV forest in one panel
# ============================================================
cox_all <- bind_rows(
  uv_f %>% mutate(model = "Univariate"),
  mv_f %>% mutate(model = "Multivariate")
) %>%
  filter(is.finite(HR), is.finite(lo), is.finite(hi)) %>%
  mutate(term_pretty = pretty_term(as.character(term)))

p_both <- ggplot(cox_all, aes(x = HR, y = term_pretty)) +
  geom_vline(xintercept = 1, linetype = 4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, color = model),
                 height = 0.15,
                 position = position_dodge(width = 0.55)) +
  geom_point(aes(color = model),
             size = 6,
             position = position_dodge(width = 0.55)) +
  scale_x_log10() +
  labs(title = "Cox model (Overall survival)",
       x = "Hazard ratio (log scale)", y = NULL) +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(face="bold")) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5.5, 20, 5.5, 5.5))


p_both

#######

# ============================================================
# Risk-stratified Table 1
# ============================================================

library(readr)
library(dplyr)
library(moonBook)
library(tidyverse)
OUT_DIR <- "./image_data/survival_model/mixture_non_fix/non_nest/beit0/results/reviewer_km_pack/dl0/n7_30_30/file04_run06"

# ----------------------------
# 1. Load data
# ----------------------------
df <- read_csv(
  file.path(OUT_DIR, "table.csv"),
  show_col_types = FALSE
)

summary(df)
out=mytable(imaging_risk_group ~ Age+pathology+stage0+recur1+survival+fu_date+brachy+duration+EQD2+feat_436+feat_519, data=df, digits=2)

mycsv(out, file="./image_data/survival_model/mixture_non_fix/non_nest/beit0/results/reviewer_km_pack/dl0/n7_30_30/file04_run06/table2.csv")


