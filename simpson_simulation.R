# ============================================================
# 공간적 심슨의 역설 시뮬레이션 (R)
#
# Z-DGP 시나리오: GRF 기반 파라미터 시나리오 (OLS / GWR / SGWR 비교)
#
# DGP (연구계획서 7절):
#   1단계: Y_i = beta_true * X_i + gamma * Z_i + epsilon_i
#   2단계: Z_i ~ GRF(h)   [rho_s = h: 공간적 군집 정도]
#   3단계: X_i = alpha * Z_i + eta_i   [alpha: X-Z 상관 강도]
#   4단계: W_i = delta * Z_i + nu_i    [delta: SGWR 대리변수 품질]
#
# 평가 지표: Bias, Bias Reduction, 부호 역전 비율
# ============================================================

# install.packages(c("GWmodel", "ggplot2", "dplyr", "patchwork", "MASS"))

library(GWmodel)
library(ggplot2)
library(dplyr)
library(patchwork)
library(MASS)

set.seed(42)

# ─────────────────────────────────────────────────────────────
# 유틸리티 함수 정의
# ─────────────────────────────────────────────────────────────

# 가우시안 랜덤 필드 생성 (squared exponential kernel)
# Omega(h) = exp(-0.5 * d_ij^2 / h^2)  [연구계획서 각주 3]
gen_grf <- function(coords, h) {
  n     <- nrow(coords)
  D     <- as.matrix(dist(coords))
  Sigma <- exp(-0.5 * D^2 / h^2) + diag(1e-6, n)
  L     <- chol(Sigma)
  z     <- as.vector(crossprod(L, rnorm(n)))
  as.vector(scale(z))
}

# SGWR: 지리적 + 속성 유사성 결합 가중 회귀
# w_comb = w_geo + impact_ratio * w_attr  (Lessani & Li 2024)
run_sgwr <- function(y, X_mat, coords, w_attr_var,
                     bw_geo, bw_attr, impact_ratio = 1.0) {
  n      <- nrow(X_mat)
  D_geo  <- as.matrix(dist(coords))
  D_attr <- as.matrix(dist(matrix(w_attr_var, ncol = 1)))
  W_geo  <- exp(-0.5 * D_geo^2  / bw_geo^2)
  W_attr <- exp(-0.5 * D_attr^2 / bw_attr^2)
  W_comb <- W_geo + impact_ratio * W_attr

  betas <- matrix(NA_real_, nrow = n, ncol = ncol(X_mat))
  colnames(betas) <- colnames(X_mat)

  for (i in seq_len(n)) {
    sq_w <- sqrt(W_comb[i, ])
    Xw   <- X_mat * sq_w
    yw   <- y     * sq_w
    tryCatch(
      betas[i, ] <- solve(crossprod(Xw), crossprod(Xw, yw)),
      error = function(e) {}
    )
  }
  betas
}

# 평가 지표 계산 (연구계획서 7절)
# Bias: mean |beta_hat_i - beta_true|
# Bias Reduction: (Bias_OLS - Bias_local) / Bias_OLS
# 부호 역전 비율: sign(beta_local) != sign(beta_OLS) 인 비율
compute_metrics <- function(beta_ols_scalar, beta_local_vec, beta_true_val) {
  bias_ols   <- abs(beta_ols_scalar - beta_true_val)
  bias_local <- mean(abs(beta_local_vec - beta_true_val), na.rm = TRUE)
  bias_red   <- ifelse(bias_ols > 1e-10,
                       (bias_ols - bias_local) / bias_ols, NA_real_)
  sign_rev   <- mean(sign(beta_local_vec) != sign(beta_ols_scalar),
                     na.rm = TRUE)
  list(bias_ols       = bias_ols,
       bias_local     = bias_local,
       bias_reduction = bias_red,
       sign_reversal  = sign_rev)
}

# ─────────────────────────────────────────────────────────────
# Z-DGP 파라미터 시나리오 시뮬레이션
# ─────────────────────────────────────────────────────────────
cat("========================================================\n")
cat("Z-DGP 파라미터 시나리오 시뮬레이션\n")
cat("========================================================\n")

# 30 × 30 격자 공간 좌표 (연구계획서 각주 3)
grid_size   <- 30
coords_grid <- as.matrix(expand.grid(x = 1:grid_size, y = 1:grid_size))
n_total     <- nrow(coords_grid)

# 공통 DGP 파라미터
# beta_true = 0: toy example과 동일하게 X→Y 진짜 효과 없음
# spurious 상관이 순수하게 Z 혼재에서만 발생
beta_true <- 0.0   # X→Y 진짜 효과 없음
gamma_z   <- 5.0   # Z→Y 효과
noise_sd  <- 0.3

# GWR/SGWR 공간 bandwidth 고정 (이론적 근거)
#   h=6(Cond1): bw < h → GWR 이웃 내 Z 분산 작음 → Z 암묵 통제 가능
#   h=0.5(Cond2/3): bw >> h → GWR 이웃 내 Z 분산 큼 → Z 통제 불가
# * AICc 최적화는 예측 정확도 기준이므로 beta_X 편향 감소와 무관
#   → 이론 검증 목적의 시뮬레이션에서는 고정 bandwidth 사용
bw_fixed <- 2.0

# 시나리오 정의 (연구계획서 7절 표)
# rho_s 높음(h=6), alpha 높음       → GWR 역설 완화  (bw < h)
# rho_s 낮음(h=0.5), alpha↑, delta↓ → GWR, SGWR 모두 완화 실패 (bw >> h)
# rho_s 낮음(h=0.5), alpha↑, delta↑ → SGWR만 역설 완화 (W가 Z의 대리변수)
# alpha = 0                          → X-Z 독립 → spurious correlation 없음
scenarios <- list(
  list(h = 6.0, alpha = 5.0, delta = 0.0,
       label = "Cond1: rho_s↑ alpha↑\n(GWR 완화 예상)"),
  list(h = 0.5, alpha = 5.0, delta = 0.5,
       label = "Cond2: rho_s↓ alpha↑ delta↓\n(완화 실패 예상)"),
  list(h = 0.5, alpha = 5.0, delta = 5.0,
       label = "Cond3: rho_s↓ alpha↑ delta↑\n(SGWR 완화 예상)"),
  list(h = 2.0, alpha = 0.0, delta = 0.0,
       label = "Cond4: alpha=0\n(역설 없음 예상)")
)

results_list <- vector("list", length(scenarios))

for (s_idx in seq_along(scenarios)) {
  sc <- scenarios[[s_idx]]
  set.seed(100 + s_idx)

  cat(sprintf("\n[시나리오 %d] h=%.1f, alpha=%.1f, delta=%.1f\n",
              s_idx, sc$h, sc$alpha, sc$delta))

  # 2단계: Z ~ GRF(h)
  z <- gen_grf(coords_grid, h = sc$h)

  # 3단계: X = alpha * Z + eta
  x <- sc$alpha * z + rnorm(n_total, sd = noise_sd)

  # 1단계: Y = beta_true * X + gamma * Z + epsilon
  y <- beta_true * x + gamma_z * z + rnorm(n_total, sd = noise_sd)

  # 4단계: W = delta * Z + nu  (SGWR용 속성 대리변수)
  w <- sc$delta * z + rnorm(n_total, sd = noise_sd)

  df_sc <- data.frame(x = x, y = y, z = z, w = w,
                      cx = coords_grid[, 1], cy = coords_grid[, 2])

  # OLS (Z 미포함)
  ols_fit  <- lm(y ~ x, data = df_sc)
  beta_ols <- coef(ols_fit)["x"]

  # OLS (Z 포함, 참조)
  beta_ols_z <- coef(lm(y ~ x + z, data = df_sc))["x"]

  # GWR (고정 bandwidth 사용)
  sp_sc <- SpatialPointsDataFrame(
    coords = coords_grid,
    data   = df_sc[, c("x", "y", "z", "w")]
  )
  cat(sprintf("  GWR 적합 중... (bw=%.1f 고정)\n", bw_fixed))
  gwr_fit  <- gwr.basic(y ~ x, data = sp_sc,
                        bw = bw_fixed, kernel = "gaussian",
                        adaptive = FALSE)
  beta_gwr <- as.data.frame(gwr_fit$SDF)$x

  # SGWR (지리 bandwidth = bw_fixed, 속성 bandwidth = sd(w)/2로 타이트하게)
  # sd(w)/2: W 유사성 기반 그룹이 더 좁게 형성 → Z 통제 효과 강화
  bw_attr       <- sd(df_sc$w) / 2
  X_mat         <- model.matrix(~ x, data = df_sc)
  cat("  SGWR 적합 중...\n")
  beta_sgwr_mat <- run_sgwr(
    y          = df_sc$y,
    X_mat      = X_mat,
    coords     = coords_grid,
    w_attr_var = df_sc$w,
    bw_geo     = bw_fixed,
    bw_attr    = bw_attr,
    impact_ratio = 2.0
  )
  beta_sgwr <- beta_sgwr_mat[, "x"]

  # 평가 지표
  m_gwr  <- compute_metrics(beta_ols, beta_gwr,  beta_true)
  m_sgwr <- compute_metrics(beta_ols, beta_sgwr, beta_true)

  # 결과 저장
  results_list[[s_idx]] <- list(
    sc = sc, df = df_sc,
    beta_ols = beta_ols, beta_ols_z = beta_ols_z,
    beta_gwr = beta_gwr, beta_sgwr = beta_sgwr,
    bw_gwr = bw_fixed, m_gwr = m_gwr, m_sgwr = m_sgwr
  )

  # 콘솔 출력
  cat(sprintf("  beta_true        = %.1f\n", beta_true))
  cat(sprintf("  OLS(no Z)        = %+.4f\n", beta_ols))
  cat(sprintf("  OLS(with Z)      = %+.4f  (참조)\n", beta_ols_z))
  cat(sprintf("  GWR  mean beta   = %+.4f  | BiasReduction=%+.1f%%  | 부호역전=%.1f%%\n",
              mean(beta_gwr,  na.rm = TRUE),
              m_gwr$bias_reduction  * 100, m_gwr$sign_reversal  * 100))
  cat(sprintf("  SGWR mean beta   = %+.4f  | BiasReduction=%+.1f%%  | 부호역전=%.1f%%\n",
              mean(beta_sgwr, na.rm = TRUE),
              m_sgwr$bias_reduction * 100, m_sgwr$sign_reversal * 100))
}

# ─────────────────────────────────────────────────────────────
# 시각화
# ─────────────────────────────────────────────────────────────

# ── 시나리오별 Z 공간 분포 & beta 분포 ──────────────────────

sc_z_plots <- lapply(seq_along(results_list), function(s) {
  res      <- results_list[[s]]
  sc_short <- sprintf("h=%.1f α=%.1f δ=%.1f", res$sc$h, res$sc$alpha, res$sc$delta)
  ggplot(res$df, aes(x = cx, y = cy, color = z)) +
    geom_point(size = 1.0, alpha = 0.8) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red",
                          midpoint = 0, guide = "none") +
    labs(title = sprintf("Cond%d Z 분포\n%s", s, sc_short),
         x = "", y = "") +
    theme_minimal(base_size = 8) +
    theme(axis.text = element_blank(), axis.ticks = element_blank())
})

sc_beta_plots <- lapply(seq_along(results_list), function(s) {
  res     <- results_list[[s]]
  n_pts   <- length(res$beta_gwr)
  beta_df <- data.frame(
    beta  = c(res$beta_gwr, res$beta_sgwr),
    model = rep(c("GWR", "SGWR"), each = n_pts)
  )
  ggplot(beta_df, aes(x = beta, fill = model)) +
    geom_histogram(bins = 25, alpha = 0.6, position = "identity") +
    geom_vline(xintercept = beta_true,
               linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_vline(xintercept = res$beta_ols,
               color = "red", linewidth = 1.0) +
    scale_fill_manual(values = c("GWR" = "steelblue", "SGWR" = "darkorange")) +
    labs(title = sprintf("Cond%d Local beta_X", s),
         subtitle = sprintf("OLS=%.2f | beta_true=%.1f",
                            res$beta_ols, beta_true),
         x = "beta_X", y = "count", fill = NULL) +
    theme_minimal(base_size = 8) +
    theme(legend.position = "bottom")
})

# ── 요약 지표 비교 플롯 ──────────────────────────────────────

summary_metrics <- do.call(rbind, lapply(seq_along(results_list), function(s) {
  res <- results_list[[s]]
  data.frame(
    scenario       = factor(paste0("Cond", s)),
    model          = c("GWR", "SGWR"),
    bias_reduction = c(res$m_gwr$bias_reduction,
                       res$m_sgwr$bias_reduction) * 100,
    sign_reversal  = c(res$m_gwr$sign_reversal,
                       res$m_sgwr$sign_reversal) * 100
  )
}))

p_bias <- ggplot(summary_metrics,
                 aes(x = scenario, y = bias_reduction, fill = model)) +
  geom_col(position = "dodge", alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("GWR" = "steelblue", "SGWR" = "darkorange")) +
  labs(title = "Bias Reduction (OLS 대비, %)",
       subtitle = "양수 = OLS보다 beta_true에 가까움",
       x = "", y = "Bias 감소율 (%)", fill = NULL) +
  theme_minimal(base_size = 10)

p_sign <- ggplot(summary_metrics,
                 aes(x = scenario, y = sign_reversal, fill = model)) +
  geom_col(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = c("GWR" = "steelblue", "SGWR" = "darkorange")) +
  labs(title = "부호 역전 비율 (%)",
       subtitle = "Local beta 부호 ≠ OLS beta 부호",
       x = "", y = "부호 역전 (%)", fill = NULL) +
  theme_minimal(base_size = 10)

# ── 최종 레이아웃 조합 ───────────────────────────────────────

# 그림 1: Z 공간 분포 (4 시나리오)
fig1 <- (sc_z_plots[[1]] | sc_z_plots[[2]] | sc_z_plots[[3]] | sc_z_plots[[4]]) +
  plot_annotation(
    title    = "시나리오별 Z 공간 분포 (GRF)",
    subtitle = sprintf("DGP: Y = %.0f*X + %.0f*Z + e,  X = alpha*Z + eta", beta_true, gamma_z),
    theme    = theme(plot.title    = element_text(size = 12, face = "bold"),
                     plot.subtitle = element_text(size = 10))
  )

ggsave("simpson_step2_Z.png", fig1, width = 12, height = 3.5, dpi = 150)
cat("그림 저장: simpson_step2_Z.png\n")

# 그림 2: Local beta 분포 (4 시나리오)
fig2 <- (sc_beta_plots[[1]] | sc_beta_plots[[2]] |
         sc_beta_plots[[3]] | sc_beta_plots[[4]]) +
  plot_annotation(
    title = "시나리오별 GWR/SGWR Local beta_X 분포",
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  )

ggsave("simpson_step2_beta.png", fig2, width = 12, height = 4, dpi = 150)
cat("그림 저장: simpson_step2_beta.png\n")

# 그림 3: 요약 지표
fig3 <- (p_bias | p_sign) +
  plot_annotation(
    title    = "모형별 평가 지표 비교 (4 시나리오)",
    subtitle = "GWR(지리적 근접성) vs SGWR(지리적 + 속성 유사성)",
    theme    = theme(plot.title    = element_text(size = 12, face = "bold"),
                     plot.subtitle = element_text(size = 10))
  )

ggsave("simpson_step2_metrics.png", fig3, width = 10, height = 5, dpi = 150)
cat("그림 저장: simpson_step2_metrics.png\n")

# ─────────────────────────────────────────────────────────────
# 최종 요약 출력
# ─────────────────────────────────────────────────────────────
cat("\n========================================================\n")
cat("최종 요약: 시나리오별 결과\n")
cat("========================================================\n")
cat(sprintf("%-8s %-8s %10s %12s %12s\n",
            "시나리오", "모형", "beta_mean", "BiasRed(%)", "부호역전(%)"))
cat(strrep("-", 58), "\n")
for (s in seq_along(results_list)) {
  res <- results_list[[s]]
  cat(sprintf("Cond%-4d OLS      %+10.4f %12s %12s\n",
              s, res$beta_ols, "-", "-"))
  cat(sprintf("         GWR      %+10.4f %12.1f %12.1f\n",
              mean(res$beta_gwr,  na.rm = TRUE),
              res$m_gwr$bias_reduction  * 100,
              res$m_gwr$sign_reversal   * 100))
  cat(sprintf("         SGWR     %+10.4f %12.1f %12.1f\n",
              mean(res$beta_sgwr, na.rm = TRUE),
              res$m_sgwr$bias_reduction * 100,
              res$m_sgwr$sign_reversal  * 100))
  cat(sprintf("         beta_true=%+.1f\n", beta_true))
  cat(strrep("-", 58), "\n")
}
