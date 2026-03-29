# ============================================================
# 공간적 심슨의 역설 시뮬레이션 (R)
# 
# Step 1: Fotheringham(2023) 방식으로 역설 재현
# Step 2: Z(교란변수) DGP로 재구성
#
# 필요 패키지: GWmodel, ggplot2, dplyr, patchwork
# install.packages(c("GWmodel", "ggplot2", "dplyr", "patchwork"))
# ============================================================

library(GWmodel)
library(ggplot2)
library(dplyr)
library(patchwork)

set.seed(42)

# ─────────────────────────────────────────────────────────────
# STEP 1: Fotheringham(2023) 방식 역설 재현
# 10개 집단, 각 100개 포인트
# 집단 내: x와 y 음의 관계 (beta 음수)
# 전체: x와 y 양의 관계 (집단 간 x 분포 분리)
# ─────────────────────────────────────────────────────────────

n_groups    <- 10
n_per_group <- 100
n_total     <- n_groups * n_per_group

# 집단별 x 평균: 15~35 (집단 간 분리)
x_means   <- seq(15, 35, length.out = n_groups)
# 집단별 beta: -15~-55 (모두 음수)
beta_true_group <- seq(-15, -55, length.out = n_groups)

# 데이터 생성
dat_list <- lapply(1:n_groups, function(g) {
  x_g   <- rnorm(n_per_group, mean = x_means[g],   sd = 1.0)
  eps_g <- rnorm(n_per_group, mean = 0,             sd = 1.0)
  y_g   <- beta_true_group[g] * x_g + eps_g

  # 집단별 y 수준을 분리하기 위해 intercept 조정
  # (집단 간 x가 클수록 y도 커야 역설 발생)
  # intercept_g를 집단 번호에 비례하게 설정
  intercept_g <- g * 600
  y_g <- y_g + intercept_g

  # 공간 좌표: 집단별로 y축 방향으로 분리
  cx <- runif(n_per_group, 0, 10)
  cy <- runif(n_per_group, (g-1)*10, g*10)

  data.frame(
    x = x_g, y = y_g,
    group = g,
    cx = cx, cy = cy
  )
})

df_foto <- do.call(rbind, dat_list)

# OLS (전역 모형)
ols_foto <- lm(y ~ x, data = df_foto)
ols_beta_foto <- coef(ols_foto)["x"]

# 집단별 OLS
group_betas <- df_foto %>%
  group_by(group) %>%
  summarise(beta = coef(lm(y ~ x))[2], .groups = "drop")

cat("========================================================\n")
cat("STEP 1: Fotheringham 방식 역설 재현\n")
cat("========================================================\n")
cat(sprintf("OLS beta_X (전역) = %.4f  (%s)\n",
            ols_beta_foto,
            ifelse(ols_beta_foto > 0, "양(+) -> 역설 조건", "음(-) -> 역설 미발생")))
cat(sprintf("집단별 beta 범위: %.1f ~ %.1f (모두 음(-))\n",
            min(group_betas$beta), max(group_betas$beta)))
cat(sprintf("=> 심슨의 역설 재현: %s\n",
            ifelse(ols_beta_foto > 0 && all(group_betas$beta < 0), "성공", "실패")))

# ─────────────────────────────────────────────────────────────
# STEP 2: Z(교란변수) DGP로 재구성
#
# DGP:
#   Z_i ~ 집단별 군집 (공간적으로 클러스터)
#   X_i = alpha * Z_i + eta_i
#   Y_i = beta_true * X_i + gamma * Z_i + epsilon_i
#   beta_true = 0 (X->Y 진짜 효과 없음)
# ─────────────────────────────────────────────────────────────

cat("\n========================================================\n")
cat("STEP 2: Z(교란변수) DGP로 재구성\n")
cat("========================================================\n")

beta_true <- 0.0   # X->Y 진짜 효과 없음
gamma_z   <- 10.0  # Z->Y 효과
alpha_z   <- 7.0   # Z->X 효과

# Z를 집단별로 다른 수준으로 설정 (공간적 군집)
z_group_means <- seq(-3, 3, length.out = n_groups)

dat_list_z <- lapply(1:n_groups, function(g) {
  z_g   <- rnorm(n_per_group, mean = z_group_means[g], sd = 0.3)
  eta_g <- rnorm(n_per_group, sd = 0.5)
  eps_g <- rnorm(n_per_group, sd = 0.5)

  x_g <- alpha_z * z_g + eta_g
  y_g <- beta_true * x_g + gamma_z * z_g + eps_g

  cx <- runif(n_per_group, 0, 10)
  cy <- runif(n_per_group, (g-1)*10, g*10)

  data.frame(
    x = x_g, y = y_g, z = z_g,
    group = g,
    cx = cx, cy = cy
  )
})

df_z <- do.call(rbind, dat_list_z)

# OLS - Z 미포함 (현실적 상황: Z 관측 불가)
ols_no_z  <- lm(y ~ x, data = df_z)
beta_no_z <- coef(ols_no_z)["x"]

# OLS - Z 포함 (참조: Z를 통제하면 어떻게 되는가)
ols_with_z  <- lm(y ~ x + z, data = df_z)
beta_with_z <- coef(ols_with_z)["x"]

cat(sprintf("beta_true = %.1f (X->Y 진짜 효과 없음)\n", beta_true))
cat(sprintf("OLS beta_X (Z 미포함) = %.4f  (%s)\n",
            beta_no_z,
            ifelse(beta_no_z > 0, "spurious 양(+)", "spurious 음(-)")))
cat(sprintf("OLS beta_X (Z 포함)   = %.4f  (참조: beta_true=0에 가까워야 함)\n",
            beta_with_z))

# ─────────────────────────────────────────────────────────────
# GWR 적합 (GWmodel 패키지)
# ─────────────────────────────────────────────────────────────

# SpatialPointsDataFrame 생성
coords_z <- as.matrix(df_z[, c("cx", "cy")])
sp_z <- SpatialPointsDataFrame(
  coords = coords_z,
  data   = df_z[, c("x", "y", "z", "group")]
)

# 최적 bandwidth 탐색 (고정 커널)
cat("\nGWR 최적 bandwidth 탐색 중...\n")
bw_opt <- bw.gwr(
  y ~ x,
  data       = sp_z,
  approach   = "AICc",
  kernel     = "gaussian",
  adaptive   = FALSE
)
cat(sprintf("최적 bandwidth (AICc): %.4f\n", bw_opt))

# GWR 적합
cat("GWR 적합 중...\n")
gwr_result <- gwr.basic(
  y ~ x,
  data      = sp_z,
  bw        = bw_opt,
  kernel    = "gaussian",
  adaptive  = FALSE
)

# GWR 로컬 계수 추출
gwr_betas <- as.data.frame(gwr_result$SDF)
df_z$gwr_beta_x <- gwr_betas$x

# 결과 요약
cat("\n========================================================\n")
cat("GWR 결과 요약\n")
cat("========================================================\n")
cat(sprintf("OLS beta_X (Z 미포함) = %.4f  (spurious)\n", beta_no_z))
cat(sprintf("GWR beta_X 전체 평균  = %.4f\n", mean(df_z$gwr_beta_x)))
cat(sprintf("GWR beta_X 음(-) 비율 = %.2f%%\n",
            mean(df_z$gwr_beta_x < 0) * 100))
cat(sprintf("OLS beta_X (Z 포함)   = %.4f  (참조)\n", beta_with_z))

# 집단별 GWR 평균
group_gwr <- df_z %>%
  group_by(group) %>%
  summarise(gwr_mean = mean(gwr_beta_x), .groups = "drop")

cat("\n집단별 GWR beta_X 평균:\n")
print(group_gwr)

# 심슨의 역설 판정
sp_ols_gwr <- (beta_no_z > 0) & (mean(df_z$gwr_beta_x < 0) > 0.5)
cat(sprintf("\n[OLS vs GWR 심슨의 역설]: %s\n",
            ifelse(sp_ols_gwr, "발생", "미발생")))
cat(sprintf("  OLS(no Z)=%.3f vs GWR 음(-) %.1f%%\n",
            beta_no_z, mean(df_z$gwr_beta_x < 0)*100))

# ─────────────────────────────────────────────────────────────
# 시각화
# ─────────────────────────────────────────────────────────────

# 1. Step1: X-Y 산점도 (Fotheringham)
p1 <- ggplot(df_foto, aes(x = x, y = y, color = factor(group))) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  geom_abline(intercept = coef(ols_foto)[1],
              slope     = coef(ols_foto)[2],
              color = "red", linewidth = 1.5, linetype = "dashed") +
  annotate("text", x = max(df_foto$x) * 0.7, y = max(df_foto$y) * 0.1,
           label = sprintf("Global OLS beta=%.2f", ols_beta_foto),
           color = "red", size = 3.5) +
  scale_color_viridis_d() +
  labs(title = "Step1: Fotheringham 방식\n(집단별 음(-), 전체 양(+))",
       x = "X", y = "Y", color = "Group") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

# 2. Step1: 집단별 beta 막대그래프
p2 <- ggplot(group_betas, aes(x = factor(group), y = beta)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_hline(yintercept = ols_beta_foto,
             color = "red", linewidth = 1.5, linetype = "dashed") +
  annotate("text", x = 5, y = ols_beta_foto * 0.8,
           label = sprintf("Global OLS=%.1f", ols_beta_foto),
           color = "red", size = 3.5) +
  labs(title = "Step1: 집단별 OLS beta_X\nvs Global OLS",
       x = "Group", y = "beta_X") +
  theme_minimal(base_size = 10)

# 3. Step2: X-Y 산점도 (Z로 색상)
p3 <- ggplot(df_z, aes(x = x, y = y, color = z)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_abline(intercept = coef(ols_no_z)[1],
              slope     = coef(ols_no_z)[2],
              color = "red", linewidth = 1.5) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, name = "Z") +
  annotate("text",
           x = min(df_z$x) + diff(range(df_z$x)) * 0.05,
           y = max(df_z$y) * 0.95,
           label = sprintf("OLS(no Z) beta=%.3f\nbeta_true=0",
                           beta_no_z),
           color = "red", size = 3.5, hjust = 0) +
  labs(title = "Step2: Z-DGP X-Y scatter\n(color=Z, beta_true=0)",
       x = "X", y = "Y") +
  theme_minimal(base_size = 10)

# 4. Step2: Z 공간 분포
p4 <- ggplot(df_z, aes(x = cx, y = cy, color = z)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, name = "Z") +
  labs(title = "Step2: Z 공간 분포\n(집단별 군집)",
       x = "x coord", y = "y coord") +
  theme_minimal(base_size = 10)

# 5. GWR 로컬 계수 분포
p5 <- ggplot(df_z, aes(x = gwr_beta_x)) +
  geom_histogram(bins = 30, fill = "steelblue",
                 color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0,
             color = "black", linewidth = 1.5,
             linetype = "dashed") +
  geom_vline(xintercept = beta_no_z,
             color = "red", linewidth = 1.5,
             linetype = "solid") +
  geom_vline(xintercept = mean(df_z$gwr_beta_x),
             color = "blue", linewidth = 1.5) +
  geom_vline(xintercept = beta_with_z,
             color = "darkgreen", linewidth = 1.5,
             linetype = "dotted") +
  annotate("text", x = beta_no_z,
           y = Inf, vjust = 1.5,
           label = sprintf("OLS(no Z)\n%.3f", beta_no_z),
           color = "red", size = 3) +
  annotate("text", x = mean(df_z$gwr_beta_x),
           y = Inf, vjust = 3,
           label = sprintf("GWR mean\n%.3f", mean(df_z$gwr_beta_x)),
           color = "blue", size = 3) +
  annotate("text", x = beta_with_z,
           y = Inf, vjust = 5,
           label = sprintf("OLS(with Z)\n%.3f", beta_with_z),
           color = "darkgreen", size = 3) +
  labs(title = "Step2: GWR local beta_X 분포",
       x = "beta_X", y = "count") +
  theme_minimal(base_size = 10)

# 6. GWR 로컬 계수 공간 분포
gwr_lim <- max(abs(df_z$gwr_beta_x))
p6 <- ggplot(df_z, aes(x = cx, y = cy, color = gwr_beta_x)) +
  geom_point(alpha = 0.7, size = 1.2) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0,
                        limits = c(-gwr_lim, gwr_lim),
                        name = "GWR\nbeta_X") +
  labs(title = "Step2: GWR local beta_X\n공간 분포",
       x = "x coord", y = "y coord") +
  theme_minimal(base_size = 10)

# 7. 요약 비교
summary_df <- data.frame(
  model = factor(
    c("OLS\n(no Z)", "GWR\nmean", "OLS\n(with Z)", "beta_true"),
    levels = c("OLS\n(no Z)", "GWR\nmean", "OLS\n(with Z)", "beta_true")
  ),
  beta  = c(beta_no_z, mean(df_z$gwr_beta_x), beta_with_z, beta_true),
  color = c("red", "steelblue", "darkgreen", "black")
)

p7 <- ggplot(summary_df, aes(x = model, y = beta, fill = color)) +
  geom_col(alpha = 0.8) +
  geom_hline(yintercept = 0,
             linetype = "dashed", color = "black") +
  scale_fill_identity() +
  geom_text(aes(label = sprintf("%.3f", beta)),
            vjust = ifelse(summary_df$beta >= 0, -0.3, 1.3),
            size = 3.5) +
  labs(title = "Summary: beta_X 추정치 비교",
       x = "", y = "beta_X") +
  theme_minimal(base_size = 10)

# 전체 레이아웃
layout <- "
AABB
CCDD
EEFG
"

final_plot <- p1 + p2 + p3 + p4 + p5 + p6 + p7 +
  plot_layout(design = layout) +
  plot_annotation(
    title = "Spatial Simpson's Paradox via Confounder Z",
    subtitle = sprintf(
      "DGP: Y = %.0f*X + %.0f*Z + e,  X = %.0f*Z + eta  |  beta_true = %.0f",
      beta_true, gamma_z, alpha_z, beta_true
    ),
    theme = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  )

ggsave("simpson_simulation_R.png",
       final_plot,
       width = 14, height = 12,
       dpi = 150)

cat("\n그림 저장 완료: simpson_simulation_R.png\n")
cat("\n========================================================\n")
cat("실험 파라미터 요약\n")
cat("========================================================\n")
cat(sprintf("  beta_true = %.1f  (X->Y 진짜 효과)\n", beta_true))
cat(sprintf("  gamma     = %.1f  (Z->Y 효과)\n", gamma_z))
cat(sprintf("  alpha     = %.1f  (Z->X 효과)\n", alpha_z))
cat(sprintf("  GWR bw    = %.4f  (최적 bandwidth)\n", bw_opt))
cat(sprintf("\n  OLS(no Z)   beta = %.4f  <- spurious\n", beta_no_z))
cat(sprintf("  GWR mean    beta = %.4f\n", mean(df_z$gwr_beta_x)))
cat(sprintf("  OLS(with Z) beta = %.4f  <- Z 통제 시\n", beta_with_z))
cat(sprintf("  beta_true        = %.4f\n", beta_true))
cat("\n  핵심 확인사항:\n")
cat(sprintf("  1. Z 미포함 OLS가 spurious 양(+): %s\n",
            ifelse(beta_no_z > 0, "YES", "NO")))
cat(sprintf("  2. Z 포함 OLS가 beta_true에 근접: %s (|diff|=%.4f)\n",
            ifelse(abs(beta_with_z) < 0.1, "YES", "NO"),
            abs(beta_with_z - beta_true)))
cat(sprintf("  3. GWR이 역설 완화(음(-) 다수): %s\n",
            ifelse(mean(df_z$gwr_beta_x < 0) > 0.5, "YES", "NO")))

