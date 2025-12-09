
library(tidyverse)
library(Hmisc)  # rcorr 사용
library(dplyr)
library(ggplot2)
library(scales)

cx <- read.csv("image_data/cx.csv")


library(dplyr)
library(ggplot2)
library(scales)

# 1) 데이터 준비
m <- cx_pre %>%
  select(survival, recur1) %>%
  filter(recur1 > 0) %>%
  mutate(
    survival_f = factor(survival, levels = c(0,1), labels = c("Alive","Dead")),
    recur_f = factor(recur1, levels = c(1,2,3),
                     labels = c("LP","DM","LP+DM"))
  )

plot_df <- m %>%
  count(survival_f, recur_f) %>%
  group_by(survival_f) %>%
  mutate(
    prop  = n / sum(n),
    label = paste0(n, " (", percent(prop, accuracy = 0.1), ")")
  )

# 2) Alive vs Dead 에서 LP/DM/LP+DM 분포 비교 (카이제곱)
tab <- table(m$survival_f, m$recur_f)
chisq_res <- chisq.test(tab)
p_val <- chisq_res$p.value

# p-value → 별표 변환 함수
p_to_stars <- function(p){
  if (p < 0.0001) "****"
  else if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "ns"
}
p_star <- p_to_stars(p_val)

p_star    # 확인용

# 3) Figure 2D + 별표
ggplot(plot_df, aes(x = survival_f, y = prop, fill = recur_f)) +
  geom_col(width = 0.8, color = "white") +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 7,
            color = "black") +
  # y축 살짝 위로 늘려서 별표 공간 확보
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.1),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("LP"="#F8766D",
                               "DM"="#00BA38",
                               "LP+DM"="#619CFF"),
                    name = "Progression") +
  # 🔹 두 막대 위를 잇는 선 + 중앙에 별표
  geom_text(aes(x = 1.5, y = 1.07, label = p_star),
            inherit.aes = FALSE, size = 8) +
  labs(
    title = "Progression pattern by survival status after progression",
    x = "",
    y = "Percent within group"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title  = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title  = element_text(size = 18, face = "bold"),
    axis.text   = element_text(size = 18, face ="bold"),
    legend.position = "top",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold")
  )
