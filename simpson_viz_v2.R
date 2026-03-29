# ============================================================
# Fotheringham Figure 5 스타일 시각화 v2
#
# 왼쪽: Y vs X 산점도 (집단별 색상 + 회귀선)
#        집단별: 음의 기울기 (beta_true=-1)
#        전체 OLS: 양의 기울기 (spurious)
#        => 심슨의 역설을 직접 보여줌
#
# 오른쪽: beta_hat * x vs x
#        회색: beta_true * x (진짜 관계)
#        빨강: OLS 추정선 (spurious)
#        파랑: GWR 로컬 추정 점 (beta_gwr_i * x_i)
# ============================================================

library(GWmodel)
library(ggplot2)
library(dplyr)
library(patchwork)

set.seed(42)

# ── 파라미터 설정 ─────────────────────────────────────────────
n_groups    <- 10
n_per_group <- 100
n_total     <- n_groups * n_per_group

beta_true <- -1.0   # X->Y 진짜 효과 (음수)
gamma_z   <- 10.0   # Z->Y 효과
alpha_z   <- 7.0    # Z->X 효과

# 집단별 Z 수준 (집단 간 분리)
z_group_means <- seq(-3, 3, length.out = n_groups)

# ── DGP 생성 ──────────────────────────────────────────────────
dat_list <- lapply(1:n_groups, function(g) {
  z_g   <- rnorm(n_per_group, mean = z_group_means[g], sd = 0.3)
  eta_g <- rnorm(n_per_group, sd = 0.5)
  eps_g <- rnorm(n_per_group, sd = 0.5)

  x_g <- alpha_z * z_g + eta_g
  y_g <- beta_true * x_g + gamma_z * z_g + eps_g

  # 공간 좌표: 집단별 y축 방향 분리 (간격 50)
  cx <- runif(n_per_group, 0, 10)
  cy <- runif(n_per_group, (g-1)*50, g*50)

  data.frame(
    x = x_g, y = y_g, z = z_g,
    group = factor(g),
    cx = cx, cy = cy
  )
})

df <- do.call(rbind, dat_list)

# ── 모형 적합 ─────────────────────────────────────────────────
# OLS (Z 미포함)
ols_fit  <- lm(y ~ x, data = df)
beta_ols <- coef(ols_fit)["x"]
int_ols  <- coef(ols_fit)["(Intercept)"]

# OLS (Z 포함, 참조)
ols_full      <- lm(y ~ x + z, data = df)
beta_ols_full <- coef(ols_full)["x"]

# GWR
coords_mat <- as.matrix(df[, c("cx", "cy")])
sp_df <- SpatialPointsDataFrame(
  coords = coords_mat,
  data   = df[, c("x", "y", "z", "group")]
)

cat("GWR 최적 bandwidth 탐색 중...\n")
bw_opt <- bw.gwr(
  y ~ x,
  data     = sp_df,
  approach = "AICc",
  kernel   = "gaussian",
  adaptive = FALSE
)
cat(sprintf("최적 bandwidth: %.4f\n", bw_opt))

cat("GWR 적합 중...\n")
gwr_fit  <- gwr.basic(
  y ~ x,
  data     = sp_df,
  bw       = bw_opt,
  kernel   = "gaussian",
  adaptive = FALSE
)

gwr_betas   <- as.data.frame(gwr_fit$SDF)
df$gwr_beta <- gwr_betas$x

# ── 결과 출력 ─────────────────────────────────────────────────
bias_ols <- abs(beta_ols - beta_true)
bias_gwr <- abs(mean(df$gwr_beta) - beta_true)

cat("\n========================================================\n")
cat("결과 요약\n")
cat("========================================================\n")
cat(sprintf("beta_true         = %.4f\n", beta_true))
cat(sprintf("OLS(no Z)  beta   = %.4f  (부호 역전: %s)\n",
            beta_ols,
            ifelse(sign(beta_ols) != sign(beta_true), "YES", "NO")))
cat(sprintf("GWR mean   beta   = %.4f\n", mean(df$gwr_beta)))
cat(sprintf("OLS(with Z) beta  = %.4f  (참조)\n", beta_ols_full))
cat(sprintf("\n편향:\n"))
cat(sprintf("  OLS(no Z): %.4f\n", bias_ols))
cat(sprintf("  GWR mean:  %.4f\n", bias_gwr))
cat(sprintf("  편향 감소율: %.1f%%\n",
            (bias_ols - bias_gwr) / bias_ols * 100))

# ── 시각화 ────────────────────────────────────────────────────

group_colors <- scales::hue_pal()(n_groups)

# 전체 OLS 추세선의 y 범위 계산
x_range <- range(df$x)
y_ols_line <- data.frame(
  x = x_range,
  y = int_ols + beta_ols * x_range
)

# ── 왼쪽: Y vs X 산점도 ───────────────────────────────────────
# 집단별 색상 + 집단별 회귀선(음) + 전체 OLS선(양)
# => 심슨의 역설 직접 시각화

p_left <- ggplot(df, aes(x = x, y = y)) +
  # 집단별 점
  geom_point(aes(color = group), alpha = 0.35, size = 1.0) +
  # 집단별 회귀선 (모두 음의 기울기)
  geom_smooth(aes(color = group, group = group),
              method = "lm", se = FALSE,
              linewidth = 0.9, alpha = 0.8) +
  # 전체 OLS (spurious, 빨강 점선)
  geom_abline(intercept = int_ols, slope = beta_ols,
              color = "red", linewidth = 1.5,
              linetype = "dashed") +
  annotate("text",
           x = max(df$x) * 0.55,
           y = min(df$y) + diff(range(df$y)) * 0.08,
           label = sprintf("Global OLS beta = %.3f\n(spurious: sign reversed)",
                           beta_ols),
           color = "red", size = 3.2, hjust = 0) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Simulated conditioned relationship\nbetween x and y",
    x = "Predictor variable (x)",
    y = "y"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# ── 오른쪽: beta_hat * x vs x ─────────────────────────────────
# 회색: beta_true * x (진짜 관계)
# 파랑: GWR 로컬 추정 (beta_gwr_i * x_i)
# 빨강: OLS 추정선 (beta_ols * x)

df$true_betax <- beta_true * df$x
df$gwr_betax  <- df$gwr_beta * df$x

# GWR 점의 추세선 (집단별로 묶어서)
p_right <- ggplot(df, aes(x = x)) +
  # 진짜 관계 (회색)
  geom_point(aes(y = true_betax),
             color = "gray70", alpha = 0.4, size = 1.0) +
  # GWR 로컬 추정 (파랑)
  geom_point(aes(y = gwr_betax),
             color = "#2166AC", alpha = 0.55, size = 1.2) +
  # OLS 추정 (빨강 선)
  geom_abline(intercept = 0, slope = beta_ols,
              color = "#D6604D", linewidth = 2.0) +
  # 수동 범례
  annotate("point", x = max(df$x) * 0.55,
           y = max(df$true_betax) * 0.92,
           color = "#D6604D", size = 3) +
  annotate("text",
           x = max(df$x) * 0.58,
           y = max(df$true_betax) * 0.92,
           label = sprintf("OLS  beta = %.3f (spurious)", beta_ols),
           color = "#D6604D", size = 3.2, hjust = 0) +
  annotate("point", x = max(df$x) * 0.55,
           y = max(df$true_betax) * 0.80,
           color = "#2166AC", size = 3) +
  annotate("text",
           x = max(df$x) * 0.58,
           y = max(df$true_betax) * 0.80,
           label = sprintf("GWR  mean = %.3f", mean(df$gwr_beta)),
           color = "#2166AC", size = 3.2, hjust = 0) +
  annotate("point", x = max(df$x) * 0.55,
           y = max(df$true_betax) * 0.68,
           color = "gray50", size = 3) +
  annotate("text",
           x = max(df$x) * 0.58,
           y = max(df$true_betax) * 0.68,
           label = sprintf("Simulated  beta_true = %.1f", beta_true),
           color = "gray50", size = 3.2, hjust = 0) +
  labs(
    title = "Modeled conditioned relationship between x and y\nusing OLS (Red) and GWR (Blue)",
    x = "Predictor variable (x)",
    y = expression(hat(beta) %.% x)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# ── 결합 ──────────────────────────────────────────────────────
final <- p_left + p_right +
  plot_annotation(
    title = "Spatial Simpson's Paradox: Simulated vs Modeled",
    subtitle = sprintf(
      "DGP: Y = %.0f·X + %.0f·Z + ε,  X = %.0f·Z + η  |  beta_true = %.1f  |  GWR bw = %.2f\nBias: OLS=%.3f, GWR=%.3f  |  Bias Reduction=%.1f%%",
      beta_true, gamma_z, alpha_z, beta_true, bw_opt,
      bias_ols, bias_gwr,
      (bias_ols - bias_gwr) / bias_ols * 100
    ),
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "gray40")
    )
  )

ggsave("simpson_figure5_v2.png",
       final,
       width = 13, height = 6,
       dpi = 150)

cat("\n그림 저장 완료: simpson_figure5_v2.png\n")
