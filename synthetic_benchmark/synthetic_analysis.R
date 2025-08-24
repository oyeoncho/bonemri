library(readxl)
library(dplyr)
library(xlsx)
library(tidyverse)
library(corrplot)
library(Hmisc)
library(moonBook)
library(RColorBrewer)
library(viridis)
library(ggrepel)

# ===== 1. 데이터 불러오기 & DeepHit 제외 =====
all <- read.csv("raw_data/survival_model_comparison_50runs_RAW1.csv") %>%
  as_tibble() %>%
  select(where(~ !any(is.na(.)))) %>%
  filter(Model != "DeepHit")  # 🔹 DeepHit 제거


# ===== 공통 색상 팔레트 (Okabe-Ito 4색) =====
model_colors <- c(
  "CoxPH"      = "#E69F00", # 주황
  "WeibullAFT" = "#56B4E9", # 하늘
  "DeepSurv"   = "#009E73", # 녹색
  "Mixture"    = "#D55E00"  # 진한 주황-빨강
)

# ===== 2. IBS 막대그래프 =====
ggplot(all, aes(x = Model, y = IBS, fill = Model)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(0.8), width = 0.6) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, position = position_dodge(0.8)) +
  stat_summary(fun = mean, geom = "text", aes(label = round(after_stat(y), 2)),
               position = position_dodge(0.8), vjust = 1.5, color = "white", fontface = "bold", size = 9) +
  facet_wrap(~ Scenario) +
  scale_fill_manual(values = model_colors) +
  theme_bw(base_size = 25) +
  labs(title = "Integrated Brier Score (IBS) by Model and Scenario",
       x = " ", y = "IBS") +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 40, face = "bold"),
    axis.text.x = element_text(size = 40, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 40, face = "bold"),
    xis.title.x  = element_text(size = 40, face = "bold"),      # x축 이름
    axis.title.y  = element_text(size = 40, face = "bold"),      # y축 이름
    plot.title    = element_text(size = 40, face = "bold", hjust = 0.5)
  )

# ===== 3. AUC 시계열 =====
auc_long <- all %>%
  pivot_longer(cols = starts_with("AUC"), names_to = "Time", values_to = "AUC") %>%
  mutate(Time = str_replace(Time, "AUC_t_", "t = 0."))

auc_heatmap <- auc_long %>%
  mutate(Time = as.numeric(gsub("AUC_t\\.", "", Time))) %>%
  group_by(Model, Scenario, Time) %>%
  summarise(AUC = mean(AUC), .groups = "drop")

end_labels_auc <- auc_heatmap %>%
  group_by(Scenario, Model) %>%
  filter(Time == max(Time)) %>%
  ungroup()

y_min <- 0.6; y_max <- 0.9

ggplot(auc_heatmap, aes(x = Time, y = AUC, color = Model, group = Model)) +
  geom_line(linewidth = 2) +
  geom_point(size = 2) +
  geom_text_repel(data = end_labels_auc, aes(label = Model),
                  size = 9, fontface = "bold", nudge_x = 0.06, hjust = 0,
                  direction = "y", segment.color = NA, box.padding = 0.2,
                  point.padding = 0.2, show.legend = FALSE) +
  facet_wrap(~ Scenario, nrow = 1, scales = "fixed") +
  scale_color_manual(values = model_colors) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.18)),
                     breaks = seq(0.1, 1.0, by = 0.1)) +
  scale_y_continuous(limits = c(y_min, y_max), breaks = seq(y_min, y_max, by = 0.05)) +
  theme_bw(base_size = 25) +
  labs(title = "Time-dependent Area Under Curve (AUC) Trends by Model",
       x = "Time", y = "AUC") +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 40, face = "bold"),
    axis.text.x = element_text(size = 33, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 35, face = "bold"),
    axis.title.x  = element_text(size = 40, face = "bold"),      # x축 이름
    axis.title.y  = element_text(size = 40, face = "bold"),      # y축 이름
    plot.title    = element_text(size = 35, face = "bold", hjust = 0.5),
        plot.margin = margin(10, 50, 10, 10))

# ===== 4. C-index 시계열 =====
cindex_long <- all %>%
  pivot_longer(cols = starts_with("Cindex"), names_to = "Time", values_to = "Cindex") %>%
  mutate(Time = as.numeric(str_replace(Time, "Cindex_t\\.\\.", "")))

cindex_heatmap <- cindex_long %>%
  group_by(Model, Scenario, Time) %>%
  summarise(Cindex = mean(Cindex, na.rm = TRUE), .groups = "drop")

end_labels_cindex <- cindex_heatmap %>%
  group_by(Scenario, Model) %>%
  filter(Time == max(Time, na.rm = TRUE)) %>%
  ungroup()

y_min <- 0.6; y_max <- 0.90

ggplot(cindex_heatmap, aes(x = Time, y = Cindex, color = Model, group = Model)) +
  geom_line(linewidth = 2) +
  geom_point(size = 2) +
  geom_text_repel(data = end_labels_cindex, aes(label = Model),
                  size = 9, fontface = "bold", nudge_x = 0.06, hjust = 0,
                  direction = "y", segment.color = NA,
                  box.padding = 0.2, point.padding = 0.2, show.legend = FALSE) +
  facet_wrap(~ Scenario, nrow = 1, scales = "fixed") +
  scale_color_manual(values = model_colors) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.18)),
                     breaks = seq(0.1, 1.0, by = 0.1)) +
  scale_y_continuous(limits = c(y_min, y_max), breaks = seq(y_min, y_max, by = 0.05)) +
  theme_bw(base_size = 25) +
  labs(title = "Time-dependent Concordance Index (C-index) Trends by Model",
       x = "Time", y = "C-index") +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 40, face = "bold"),
    axis.text.x = element_text(size = 33, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 35, face = "bold"),
    axis.title.x  = element_text(size = 40, face = "bold"),      # x축 이름
    axis.title.y  = element_text(size = 40, face = "bold"),      # y축 이름
    plot.title    = element_text(size = 34, face = "bold", hjust = 0.5),
    plot.margin = margin(10, 50, 10, 10))

library(tidyverse)

raw <- read_csv("raw_data/survival_model_comparison_50runs_RAW1.csv")  # 파일명 맞춰 변경

auc_cols  <- names(raw)[startsWith(names(raw), "AUC_t=")]
cidx_cols <- names(raw)[startsWith(names(raw), "Cindex_t>=")]

df <- raw %>%
  mutate(
    AUC_global = rowMeans(select(., all_of(auc_cols)), na.rm = TRUE),
    CIDX_global = rowMeans(select(., all_of(cidx_cols)), na.rm = TRUE)
  )

summary_tbl <- df %>%
  group_by(Scenario, Model) %>%
  summarise(
    runs = n(),
    IBS_mean  = mean(IBS, na.rm = TRUE),
    IBS_sd    = sd(IBS,   na.rm = TRUE),
    AUCg_mean = mean(AUC_global,  na.rm = TRUE),
    AUCg_sd   = sd(AUC_global,    na.rm = TRUE),
    CIDXg_mean = mean(CIDX_global, na.rm = TRUE),
    CIDXg_sd   = sd(CIDX_global,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Scenario, match(Model, c("CoxPH","Mixture","DeepSurv","WeibullAFT")))
print(summary_tbl)

# df: 위에서 만든 raw + AUC_global + CIDX_global 이 있다고 가정
# DeepHit 제외하려면 filter(Model != "DeepHit") 추가

# df: 위에서 만든 raw + AUC_global + CIDX_global 이 있다고 가정
# DeepHit 제외하려면 filter(Model != "DeepHit") 추가

pair_tbl <- df %>%
  filter(Model != "DeepHit") %>%
  select(Scenario, Model, Run, IBS, AUC_global, CIDX_global) %>%
  group_by(Scenario) %>%
  group_modify(~{
    mix  <- filter(.x, Model == "Mixture") %>%
      select(Run, IBS_mix = IBS, AUCg_mix = AUC_global, CIDXg_mix = CIDX_global)
    comp <- filter(.x, Model != "Mixture") %>%
      select(Model, Run, IBS_cmp = IBS, AUCg_cmp = AUC_global, CIDXg_cmp = CIDX_global)
    
    inner_join(mix, comp, by = "Run") %>%
      group_by(Model) %>%
      summarise(
        n = n(),
        
        # --- 쌍대비 차이 (Mixture - 비교모델)
        dIBS  = IBS_mix  - IBS_cmp,
        dAUC  = AUCg_mix - AUCg_cmp,
        dCIDX = CIDXg_mix - CIDXg_cmp,
        
        # 요약 통계 + 95% CI(평균의 정규근사)
        dIBS_mean  = mean(dIBS),  dIBS_sd  = sd(dIBS),
        dAUC_mean  = mean(dAUC),  dAUC_sd  = sd(dAUC),
        dCIDX_mean = mean(dCIDX), dCIDX_sd = sd(dCIDX),
        dIBS_ci_l  = dIBS_mean  - 1.96 * dIBS_sd  / sqrt(n),
        dIBS_ci_u  = dIBS_mean  + 1.96 * dIBS_sd  / sqrt(n),
        dAUC_ci_l  = dAUC_mean  - 1.96 * dAUC_sd  / sqrt(n),
        dAUC_ci_u  = dAUC_mean  + 1.96 * dAUC_sd  / sqrt(n),
        dCIDX_ci_l = dCIDX_mean - 1.96 * dCIDX_sd / sqrt(n),
        dCIDX_ci_u = dCIDX_mean + 1.96 * dCIDX_sd / sqrt(n),
        
        # Wilcoxon signed-rank (차이의 중앙값=0?)
        p_IBS  = {x <- na.omit(dIBS);  if (length(x)<2 || all(x==0)) NA_real_
        else wilcox.test(x, mu=0, exact=FALSE, correct=FALSE)$p.value},
        p_AUC  = {x <- na.omit(dAUC);  if (length(x)<2 || all(x==0)) NA_real_
        else wilcox.test(x, mu=0, exact=FALSE, correct=FALSE)$p.value},
        p_CIDX = {x <- na.omit(dCIDX); if (length(x)<2 || all(x==0)) NA_real_
        else wilcox.test(x, mu=0, exact=FALSE, correct=FALSE)$p.value},
        
        # --- win rate (Mixture가 이긴 비율)
        win_IBS  = mean(dIBS  < 0, na.rm = TRUE),  # IBS는 작을수록 좋음
        win_AUC  = mean(dAUC  > 0, na.rm = TRUE),  # AUC는 클수록 좋음
        win_CIDX = mean(dCIDX > 0, na.rm = TRUE),  # C-index도 클수록 좋음
        
        # 동시 승률(IBS↓ & AUC↑), 세 지표 모두 승 (참고용)
        win_both_IBS_AUC = mean((dIBS<0) & (dAUC>0), na.rm = TRUE),
        win_all3         = mean((dIBS<0) & (dAUC>0) & (dCIDX>0), na.rm = TRUE),
        
        # 효과크기 대용: probability of superiority (=win rate), Cliff's delta
        prob_sup_AUC  = win_AUC,
        prob_sup_CIDX = win_CIDX,
        cliffs_AUC    = 2*win_AUC  - 1,
        cliffs_CIDX   = 2*win_CIDX - 1,
        .groups = "drop"
      )
  }) %>%
  group_by(Scenario) %>%                       # 시나리오 내 비교들에 Holm 보정
  mutate(
    p_IBS_adj  = p.adjust(p_IBS,  method = "holm"),
    p_AUC_adj  = p.adjust(p_AUC,  method = "holm"),
    p_CIDX_adj = p.adjust(p_CIDX, method = "holm")
  ) %>%
  ungroup() %>%
  arrange(Scenario, factor(Model, levels = c("CoxPH","WeibullAFT","DeepSurv")))


readr::write_csv(pair_tbl, "pairwise_mixture_vs_others_with_winrates.csv")

summary(pair_tbl)
out = mytable(Model~., data=pair_tbl, digits=3, method=3)
mycsv(out, file="suppl_table.csv")
