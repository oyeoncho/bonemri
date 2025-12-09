library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(stringr)

## -----------------------------
## 1) 설정
## -----------------------------
T0 <- 60
base_dir <- "image_data/survival_model/mixture_non_fix/non_nest/beit0/results/generalization/test1_1/dl0/n7_30_30/file04"

fname_internal_image   <- file.path(base_dir, "dca_internal_file04_run06_Image_only_T60.csv")
fname_internal_clin    <- file.path(base_dir, "dca_internal_file04_run06_Clinical_only_T60.csv")
fname_internal_both    <- file.path(base_dir, "dca_internal_file04_run06_Image_+_Clinical_T60.csv")

fname_external_image   <- file.path(base_dir, "dca_external_file04_run06_Image_only_T60.csv")
fname_external_clin    <- file.path(base_dir, "dca_external_file04_run06_Clinical_only_T60.csv")
fname_external_both    <- file.path(base_dir, "dca_external_file04_run06_Image_+_Clinical_T60.csv")

thresholds <- seq(0.01, 0.80, by = 0.01)

## -----------------------------
## 2) DCA 계산 함수
## -----------------------------
dca_from_df <- function(df, model, T0 = 60, threshs = seq(0.01, 0.80, by = 0.01)) {
  df <- df %>%
    mutate(event_T0 = ifelse(time <= T0 & event == 1, 1, 0))
  n <- nrow(df)
  prev <- mean(df$event_T0 == 1)
  
  out_list <- lapply(threshs, function(pt) {
    pred_pos <- df$risk >= pt
    TP <- sum(df$event_T0 == 1 & pred_pos)
    FP <- sum(df$event_T0 == 0 & pred_pos)
    
    nb_model <- (TP / n) - (FP / n) * (pt / (1 - pt))
    nb_all   <- prev - (1 - prev) * (pt / (1 - pt))
    nb_none  <- 0
    
    tibble(
      threshold = pt,
      model = model,
      NB_model = nb_model,
      NB_all   = nb_all,
      NB_none  = nb_none
    )
  })
  
  bind_rows(out_list)
}

dca_from_file <- function(path, model, T0 = 60, threshs = seq(0.01, 0.80, by = 0.01)) {
  df <- read_csv(path, show_col_types = FALSE) %>%
    select(time, event, risk)
  dca_from_df(df, model = model, T0 = T0, threshs = threshs)
}

## -----------------------------
## 3) Internal DCA
## -----------------------------
dca_int_img  <- dca_from_file(fname_internal_image, "Image only",    T0, thresholds)
dca_int_cli  <- dca_from_file(fname_internal_clin,  "Clinical only", T0, thresholds)
dca_int_both <- dca_from_file(fname_internal_both,  "Image + Clinical", T0, thresholds)

dca_internal_all <- bind_rows(dca_int_img, dca_int_cli, dca_int_both)

# treat-all / none를 long 포맷으로 변환해서 같이 플롯
nb_ref_internal_long <- dca_int_img %>%
  select(threshold, NB_all, NB_none) %>%
  distinct() %>%
  pivot_longer(cols = c(NB_all, NB_none),
               names_to = "model",
               values_to = "NB") %>%
  mutate(model = recode(model,
                        NB_all  = "Treat all",
                        NB_none = "Treat none"))

plot_internal_models <- dca_internal_all %>%
  select(threshold, model, NB_model) %>%
  rename(NB = NB_model)

# 모델 + 참조선 모두 합치기
plot_internal_all <- bind_rows(plot_internal_models, nb_ref_internal_long)

# 색/선형 정의 (원하는 색으로 조정 가능)
model_levels <- c("Image only", "Clinical only", "Image + Clinical",
                  "Treat all", "Treat none")

plot_internal_all$model <- factor(plot_internal_all$model,
                                  levels = model_levels)

cols <- c("Image only"       = "blue",
          "Clinical only"    = "red",
          "Image + Clinical" = "green4",
          "Treat all"        = "black",
          "Treat none"       = "black")

lts  <- c("Image only"       = "solid",
          "Clinical only"    = "solid",
          "Image + Clinical" = "solid",
          "Treat all"        = "dashed",
          "Treat none"       = "dotted")

## -----------------------------
## Internal plot
## -----------------------------
p_int <- ggplot(plot_internal_all,
                aes(x = threshold, y = NB,
                    color = model, linetype = model)) +
  geom_line(size = 1) +
  scale_color_manual(values = cols, name = "Model / Strategy") +
  scale_linetype_manual(values = lts, name = "Model / Strategy") +
  # 🔹 x, y 축 간격 0.1로 설정
  scale_x_continuous(breaks = seq(0, 0.8, by = 0.1), limits = c(0, 0.8)) +
  scale_y_continuous(breaks = seq(-0.2, 0.15, by = 0.05),
                     limits = c(-0.2, 0.15)) +
  labs(
    title = sprintf("Decision Curve Analysis (Internal, T = %d months)", T0),
    x     = "Threshold probability",
    y     = "Net benefit"
  ) +
  # 🔹 글자 크게
  theme_bw(base_size = 18) +
  theme(
    legend.position   = "bottom",
    axis.title.x      = element_text(size = 18),
    axis.title.y      = element_text(size = 18),
    axis.text.x       = element_text(size = 16),
    axis.text.y       = element_text(size = 16),
    legend.title      = element_blank(),
    legend.text       = element_text(size = 14),
    plot.title        = element_text(size = 20, face = "bold", hjust = 0.5)
  )

print(p_int)


## -----------------------------
## 4) External DCA (동일 방식)
## -----------------------------
dca_ext_img  <- dca_from_file(fname_external_image, "Image only",    T0, thresholds)
dca_ext_cli  <- dca_from_file(fname_external_clin,  "Clinical only", T0, thresholds)
dca_ext_both <- dca_from_file(fname_external_both,  "Image + Clinical", T0, thresholds)

dca_external_all <- bind_rows(dca_ext_img, dca_ext_cli, dca_ext_both)

nb_ref_external_long <- dca_ext_img %>%
  select(threshold, NB_all, NB_none) %>%
  distinct() %>%
  pivot_longer(cols = c(NB_all, NB_none),
               names_to = "model",
               values_to = "NB") %>%
  mutate(model = recode(model,
                        NB_all  = "Treat all",
                        NB_none = "Treat none"))

plot_external_models <- dca_external_all %>%
  select(threshold, model, NB_model) %>%
  rename(NB = NB_model)

plot_external_all <- bind_rows(plot_external_models, nb_ref_external_long)
plot_external_all$model <- factor(plot_external_all$model,
                                  levels = model_levels)

## -----------------------------
## External plot
## -----------------------------
p_ext <- ggplot(plot_external_all,
                aes(x = threshold, y = NB,
                    color = model, linetype = model)) +
  geom_line(size = 1) +
  scale_color_manual(values = cols, name = "Model / Strategy") +
  scale_linetype_manual(values = lts, name = "Model / Strategy") +
  # 🔹 x, y 축 간격 0.1로 설정
  scale_x_continuous(breaks = seq(0, 0.8, by = 0.1), limits = c(0, 0.8)) +
  scale_y_continuous(breaks = seq(-0.2, 0.4, by = 0.05),
                     limits = c(-0.2, 0.4)) +
  labs(
    title = sprintf("Decision Curve Analysis (External, T = %d months)", T0),
    x     = "Threshold probability",
    y     = "Net benefit"
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.position   = "bottom",
    axis.title.x      = element_text(size = 18),
    axis.title.y      = element_text(size = 18),
    axis.text.x       = element_text(size = 16),
    axis.text.y       = element_text(size = 16),
    legend.title      = element_blank(),
    legend.text       = element_text(size = 14),
    plot.title        = element_text(size = 20, face = "bold", hjust = 0.5)
  )

print(p_ext)
