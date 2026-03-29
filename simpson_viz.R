# ============================================================
# Fotheringham Figure 5 스타일 비교 시각화
#
# 왼쪽: Simulated conditioned relationship (집단별 색상)
# 오른쪽: Modeled relationship - OLS(Red) vs GWR(Blue)
#
# 본 연구 DGP:
#   Z_i ~ 집단별 군집 (공간적 교란변수)
#   X_i = alpha * Z_i + eta_i
#   Y_i = beta_true * X_i + gamma * Z_i + epsilon_i
#   beta_true = -1 (X->Y 진짜 효과: 음수)
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

# ── DGP 생성 ─────────────────────────────────────────────────
dat_list <- lapply(1:n_groups, function(g) {
  # Z: 집단별로 다른 수준 (공간적 군집)
  z_g   <- rnorm(n_per_group, mean = z_group_means[g], sd = 0.3)
  eta_g <- rnorm(n_per_group, sd = 0.5)
  eps_g <- rnorm(n_per_group, sd = 0.5)

  # X = alpha*Z + eta
  x_g <- alpha_z * z_g + eta_g

  # Y = beta_true*X + gamma*Z + epsilon
  y_g <- beta_true * x_g + gamma_z * z_g + eps_g

  # 공간 좌표: 집단별로 y축 방향 분리 (간격 50)
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
ols_fit   <- lm(y ~ x, data = df)
beta_ols  <- coef(ols_fit)["x"]

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
gwr_fit   <- gwr.basic(
  y ~ x,
  data     = sp_df,
  bw       = bw_opt,
  kernel   = "gaussian",
  adaptive = FALSE
)

gwr_betas    <- as.data.frame(gwr_fit$SDF)
df$gwr_beta  <- gwr_betas$x
df$gwr_yhat  <- gwr_betas$yhat

# ── 결과 출력 ─────────────────────────────────────────────────
cat("\n========================================================\n")
cat("결과 요약\n")
cat("========================================================\n")
cat(sprintf("beta_true        = %.4f\n", beta_true))
cat(sprintf("OLS(no Z)  beta  = %.4f  (spurious: %s)\n",
            beta_ols,
            ifelse(sign(beta_ols) != sign(beta_true), "YES", "NO")))
cat(sprintf("GWR mean   beta  = %.4f\n", mean(df$gwr_beta)))
cat(sprintf("OLS(with Z) beta = %.4f  (참조)\n", beta_ols_full))
cat(sprintf("\n편향 (|beta - beta_true|):\n"))
cat(sprintf("  OLS(no Z): %.4f\n", abs(beta_ols - beta_true)))
cat(sprintf("  GWR mean:  %.4f\n", abs(mean(df$gwr_beta) - beta_true)))
cat(sprintf("  편향 감소율: %.1f%%\n",
            (abs(beta_ols - beta_true) - abs(mean(df$gwr_beta) - beta_true)) /
              abs(beta_ols - beta_true) * 100))

# ── 시각화: Fotheringham Figure 5 스타일 ─────────────────────

# beta*x 계산 (conditioned relationship)
# Fotheringham Figure 5의 y축은 beta*x
df$true_beta_x <- beta_true * df$x          # True: beta_true * x
df$ols_beta_x  <- beta_ols  * df$x          # OLS:  beta_ols * x
df$gwr_beta_x_val <- df$gwr_beta * df$x     # GWR:  local_beta * x

# 집단 색상 (Fotheringham과 유사하게 무지개)
group_colors <- scales::hue_pal()(n_groups)

# ── 왼쪽 패널: Simulated (집단별 색상) ───────────────────────
# Fotheringham: 집단별 음의 기울기, 전체 양의 기울기
# 본 연구: 집단별 양의 기울기(spurious), 진짜는 음의 기울기

# 왼쪽에 표시할 것:
# - 집단별 색상 점
# - 각 집단의 진짜 조건부 관계선 (beta_true * x, 집단별)
# - 전체 OLS 선 (회색, spurious)

p_left <- ggplot(df, aes(x = x, y = true_beta_x)) +
  # 집단별 점
  geom_point(aes(color = group), alpha = 0.4, size = 1.2) +
  # 집단별 진짜 관계선 (beta_true * x이므로 모두 동일 기울기)
  geom_smooth(aes(color = group, group = group),
              method = "lm", se = FALSE, linewidth = 0.8) +
  # 전체 OLS (spurious, 회색)
  geom_abline(intercept = 0, slope = beta_true,
              color = "gray40", linewidth = 1.2,
              linetype = "solid") +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Simulated conditioned relationship\nbetween x and y",
    x = "Predictor variable (x)",
    y = expression(beta[true] %.% x)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# ── 오른쪽 패널: Modeled - OLS(Red) vs GWR(Blue) ─────────────
# Fotheringham과 동일한 구성:
# - 회색: simulated (진짜 beta*x)
# - 파란색: GWR 로컬 추정치 (beta_gwr_i * x_i)
# - 빨간색: OLS 추정치 (beta_ols * x)

p_right <- ggplot(df, aes(x = x)) +
  # Simulated (회색 배경)
  geom_point(aes(y = true_beta_x),
             color = "gray70", alpha = 0.4, size = 1.0) +
  # GWR 로컬 추정 (파란색)
  geom_point(aes(y = gwr_beta_x_val),
             color = "#2166AC", alpha = 0.5, size = 1.2) +
  # OLS 추정 (빨간색 선)
  geom_abline(intercept = 0, slope = beta_ols,
              color = "#D6604D", linewidth = 2.0) +
  # 범례용 dummy
  geom_point(data = data.frame(x = NA, y = NA, type = "MGWR"),
             aes(x = x, y = y, shape = type),
             color = "#2166AC", size = 3) +
  annotate("text",
           x = max(df$x) * 0.6,
           y = max(df$true_beta_x) * 0.85,
           label = sprintf("OLS beta=%.3f\n(spurious)", beta_ols),
           color = "#D6604D", size = 3.5, hjust = 0) +
  annotate("text",
           x = max(df$x) * 0.6,
           y = max(df$true_beta_x) * 0.65,
           label = sprintf("GWR mean=%.3f\nbeta_true=%.1f",
                           mean(df$gwr_beta), beta_true),
           color = "#2166AC", size = 3.5, hjust = 0) +
  labs(
    title = sprintf(
      "Modeled conditioned relationship between x and y\nusing OLS (Red) and GWR (Blue)"
    ),
    x = "Predictor variable (x)",
    y = expression(hat(beta) %.% x)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# ── 최종 결합 ─────────────────────────────────────────────────
final <- p_left + p_right +
  plot_annotation(
    title = "Spatial Simpson's Paradox: Simulated vs Modeled",
    subtitle = sprintf(
      "DGP: Y = %.0f·X + %.0f·Z + ε,  X = %.0f·Z + η  |  beta_true = %.1f  |  GWR bw = %.2f",
      beta_true, gamma_z, alpha_z, beta_true, bw_opt
    ),
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "gray40")
    )
  )

ggsave("simpson_figure5_style.png",
       final,
       width = 12, height = 5.5,
       dpi = 150)

cat("\n그림 저장 완료: simpson_figure5_style.png\n")
