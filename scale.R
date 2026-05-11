# =============================================================================
# Scale Alignment between GWR Bandwidth and Spatial Confounder Scale
# -----------------------------------------------------------------------------
# 질문: 공간 교란 Z의 스케일 ℓ_Z 와 GWR 대역폭 h가 일치할 때,
#       GWR이 OLS 대비 mean bias를 줄일 수 있는가?
#
# 이론적 예측 (Paciorek): h ≈ sqrt(ℓ_X · ℓ_Z) 근방에서 bias 최소
# =============================================================================

library(GWmodel)
library(sp)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(viridis)
library(patchwork)

set.seed(20260508)

# ---- 1. 설정 -----------------------------------------------------------------

CFG <- list(
  n_side       = 25,
  domain_size  = 10,
  scale_X      = 0.5,
  alpha_XZ     = 0.6,
  gamma_true   = 1.5,
  beta_true    = 1.0,
  sigma_eps    = 0.3,
  
  # 실험 그리드
  scale_Z_grid = c(0.3, 0.7, 1.5, 3.0, 6.0),
  h_grid       = c(20, 50, 120, 300, 600),
  n_reps       = 10
)

# ---- 2. Helper functions ----------------------------------------------------

sim_grf <- function(coords, range_param, sigma2 = 1) {
  d <- as.matrix(dist(coords))
  Sigma <- sigma2 * exp(-d / range_param) + 1e-6 * diag(nrow(d))
  L <- chol(Sigma)
  as.numeric(t(L) %*% rnorm(nrow(d)))
}

generate_dgp_B <- function(coords, scale_Z, cfg) {
  Z      <- sim_grf(coords, scale_Z)
  X_fine <- sim_grf(coords, cfg$scale_X)
  X      <- cfg$alpha_XZ * Z + X_fine
  y      <- cfg$beta_true * X + cfg$gamma_true * Z +
            rnorm(nrow(coords), sd = cfg$sigma_eps)
  data.frame(coords, X = X, Z = Z, y = y)
}

#' (scale_Z, h) 한 셀에서 한 번의 시뮬레이션
single_run <- function(coords, scale_Z, h, cfg) {
  df <- generate_dgp_B(coords, scale_Z, cfg)
  sp_df <- SpatialPointsDataFrame(coords = df[, c("u", "v")], data = df)
  
  # GWR 적합
  fit <- tryCatch(
    gwr.basic(y ~ X, data = sp_df, bw = h,
              kernel = "bisquare", adaptive = TRUE),
    error = function(e) NULL
  )
  
  # OLS 참조
  ols_bias <- coef(lm(y ~ X, data = df))["X"] - cfg$beta_true
  
  if (is.null(fit)) {
    return(data.frame(scale_Z = scale_Z, h = h,
                      mean_beta = NA, sd_beta = NA,
                      gwr_bias = NA, ols_bias = ols_bias))
  }
  
  beta_hat <- fit$SDF$X
  data.frame(
    scale_Z   = scale_Z,
    h         = h,
    mean_beta = mean(beta_hat),
    sd_beta   = sd(beta_hat),
    gwr_bias  = mean(beta_hat) - cfg$beta_true,
    ols_bias  = ols_bias
  )
}

# ---- 3. 메인 시뮬레이션 ------------------------------------------------------

cat("[1/3] Setting up grid...\n")
s_seq  <- seq(0, CFG$domain_size, length.out = CFG$n_side)
coords <- expand.grid(u = s_seq, v = s_seq)

# 모든 (scale_Z, h, rep) 조합
grid_design <- expand.grid(
  scale_Z = CFG$scale_Z_grid,
  h       = CFG$h_grid,
  rep     = seq_len(CFG$n_reps)
)
n_total <- nrow(grid_design)
cat(sprintf("Total simulations: %d\n", n_total))

cat("[2/3] Running simulations (this may take 30-60 min)...\n")
pb <- txtProgressBar(min = 0, max = n_total, style = 3)
results <- map_dfr(seq_len(n_total), function(i) {
  out <- single_run(coords,
                    grid_design$scale_Z[i],
                    grid_design$h[i],
                    CFG)
  out$rep <- grid_design$rep[i]
  setTxtProgressBar(pb, i)
  out
})
close(pb)

# ---- 4. 결과 집계 ------------------------------------------------------------

cat("[3/3] Aggregating results...\n")

# 각 (scale_Z, h) 셀에서 reps에 걸친 중앙값
agg <- results %>%
  group_by(scale_Z, h) %>%
  summarize(
    median_gwr_bias = median(gwr_bias, na.rm = TRUE),
    mad_gwr_bias   = mad(gwr_bias, na.rm = TRUE),
    median_sd_beta = median(sd_beta, na.rm = TRUE),
    median_ols_bias = median(ols_bias, na.rm = TRUE),
    bias_reduction = median(ols_bias - gwr_bias, na.rm = TRUE),
    n              = sum(!is.na(gwr_bias)),
    .groups = "drop"
  ) %>%
  mutate(
    scale_Z_label = factor(scale_Z),
    h_label       = factor(h),
    # Paciorek 예측: h* ≈ sqrt(ℓ_X · ℓ_Z), 단위는 격자 NN 수가 아니라 거리이므로 환산 필요
    # 격자점 간 거리: domain_size/n_side, NN 수와 거리의 근사 관계
    # 단순화: log scale로 비교
    log_h = log(h),
    log_scaleZ = log(scale_Z)
  )

print(agg)

# ---- 5. 시각화 --------------------------------------------------------------

# Panel 1: GWR bias의 heatmap
p_heat <- agg %>%
  ggplot(aes(x = h_label, y = scale_Z_label, fill = abs(median_gwr_bias))) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%+.2f", median_gwr_bias)),
            color = "white", size = 3.5) +
  scale_fill_viridis_c(name = "|bias|", option = "magma", direction = -1) +
  labs(title = "GWR mean bias across (ℓ_Z, h) grid",
       subtitle = sprintf("ℓ_X = %.1f fixed; numbers = signed median bias",
                          CFG$scale_X),
       x = "Bandwidth h (adaptive NN)",
       y = expression("Z scale " ~ ell[Z])) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())

# Panel 2: 각 ℓ_Z 별 bias의 h-궤적 + OLS 참조선
p_lines <- agg %>%
  ggplot(aes(x = h, color = scale_Z_label, group = scale_Z_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_line(aes(y = median_gwr_bias), linewidth = 0.8) +
  geom_point(aes(y = median_gwr_bias), size = 2.5) +
  geom_line(aes(y = median_ols_bias), linewidth = 0.5,
            linetype = "dotted") +
  scale_x_log10() +
  scale_color_viridis_d(name = expression(ell[Z]), option = "viridis") +
  labs(title = "GWR bias vs bandwidth, by Z scale",
       subtitle = "Solid = GWR bias; dotted = OLS bias (reference)",
       x = "Bandwidth h (log scale)",
       y = "Median bias") +
  theme_minimal(base_size = 11)

# Panel 3: bias 감소 (OLS - GWR) heatmap — 양수면 GWR이 OLS보다 좋음
p_reduction <- agg %>%
  ggplot(aes(x = h_label, y = scale_Z_label, fill = bias_reduction)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%+.2f", bias_reduction)),
            color = "black", size = 3.5) +
  scale_fill_gradient2(name = "OLS - GWR bias",
                       low = "#d73027", mid = "white", high = "#1a9850",
                       midpoint = 0) +
  labs(title = "Bias reduction: OLS bias minus GWR bias",
       subtitle = "Green = GWR helps; red = GWR hurts (Hodges-Reich region)",
       x = "Bandwidth h", y = expression(ell[Z])) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())

# Panel 4: 거짓 변동(sd_beta) heatmap
p_sd <- agg %>%
  ggplot(aes(x = h_label, y = scale_Z_label, fill = median_sd_beta)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", median_sd_beta)),
            color = "white", size = 3.5) +
  scale_fill_viridis_c(name = "sd(β̂)", option = "plasma") +
  labs(title = "False spatial variation: sd of β̂(s)",
       subtitle = "진실 sd = 0 (β is constant)",
       x = "Bandwidth h", y = expression(ell[Z])) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())

# 결합
final_plot <- (p_heat | p_reduction) / (p_lines | p_sd) +
  plot_annotation(
    title = "Scale alignment between GWR bandwidth and Z's spatial scale",
    subtitle = sprintf("Does h ≈ ℓ_Z reduce bias under DGP B? (γ=%.1f, α_XZ=%.1f, ℓ_X=%.1f)",
                       CFG$gamma_true, CFG$alpha_XZ, CFG$scale_X),
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

print(final_plot)
ggsave("scale_alignment.png", final_plot, width = 13, height = 10, dpi = 150)

# ---- 6. 핵심 진단 통계 -------------------------------------------------------

cat("\n=== 핵심 발견 ===\n\n")

# 각 scale_Z 별 최적 h (bias 최소)
optimal_h <- agg %>%
  group_by(scale_Z) %>%
  slice_min(abs(median_gwr_bias), n = 1) %>%
  select(scale_Z, h, median_gwr_bias, median_ols_bias) %>%
  mutate(predicted_h_paciorek = sqrt(CFG$scale_X * scale_Z))

cat("1. 각 scale_Z 별 bias-최소 h:\n")
print(optimal_h)

# Paciorek 예측과의 비교
cat("\n2. Paciorek 예측 (h* ≈ sqrt(ℓ_X · ℓ_Z)) vs 관측 최적 h:\n")
cat("   상관:",
    cor(log(optimal_h$h), log(optimal_h$predicted_h_paciorek)), "\n")

# Hodges-Reich 영역 식별 (GWR이 OLS보다 나쁜 셀)
hr_cells <- agg %>% filter(bias_reduction < 0)
cat(sprintf("\n3. Hodges-Reich 영역 (GWR이 OLS보다 나쁨): %d / %d 셀\n",
            nrow(hr_cells), nrow(agg)))
if (nrow(hr_cells) > 0) {
  cat("   해당 셀들:\n")
  print(hr_cells %>% select(scale_Z, h, median_gwr_bias, median_ols_bias))
}

cat("\n=== Done ===\n")