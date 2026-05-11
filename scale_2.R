# =============================================================================
# Standard GWR vs Mixed GWR — DGP B에서 OVB 보정 비교
# -----------------------------------------------------------------------------
# Standard GWR: y ~ β(s)·x + ε        (절편, 기울기 모두 변동)
# Mixed GWR:    y ~ α(s) + β·x + ε    (절편만 변동, 기울기 전역 상수)
#
# Jiho 직관: 공간 인접성 = 국지 confounder 통제 ≈ panel fixed effects
# 이 직관의 *깨끗한 구현*이 Mixed GWR (β를 전역 상수로 강제)
# =============================================================================

library(GWmodel)
library(sp)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(patchwork)
library(viridis)

set.seed(20260508)

# ---- 1. 설정 -----------------------------------------------------------------

CFG <- list(
  n_side       = 15,          # 25→15: 625→225 포인트 (GWR은 O(n²), 핵심 절감)
  domain_size  = 10,
  scale_X      = 0.5,
  scale_Z      = 3.0,
  alpha_XZ     = 0.6,
  gamma_true   = 1.5,
  beta_true    = 1.0,
  sigma_eps    = 0.3,

  bw_grid      = c(20, 50, 120, 300),   # 5→4개
  n_reps       = 4            # 8→4회
)

# ---- 2. Helpers --------------------------------------------------------------

sim_grf <- function(coords, range_param, sigma2 = 1) {
  d <- as.matrix(dist(coords))
  Sigma <- sigma2 * exp(-d / range_param) + 1e-6 * diag(nrow(d))
  L <- chol(Sigma)
  as.numeric(t(L) %*% rnorm(nrow(d)))
}

generate_dgp_B <- function(coords, cfg) {
  Z      <- sim_grf(coords, cfg$scale_Z)
  X_fine <- sim_grf(coords, cfg$scale_X)
  X      <- cfg$alpha_XZ * Z + X_fine
  y      <- cfg$beta_true * X + cfg$gamma_true * Z +
            rnorm(nrow(coords), sd = cfg$sigma_eps)
  data.frame(coords, X = X, Z = Z, y = y)
}

# ---- 3. 한 시뮬레이션에서 Standard + Mixed 적합 ------------------------------

fit_both <- function(sp_df, df, bw, cfg) {
  # Standard GWR: 절편 + 기울기 모두 국지
  fit_std <- tryCatch(
    gwr.basic(y ~ X, data = sp_df, bw = bw,
              kernel = "bisquare", adaptive = TRUE),
    error = function(e) NULL
  )
  
  # Mixed GWR: X 계수를 전역 상수로 고정 (절편만 국지)
  fit_mix <- tryCatch(
    gwr.mixed(y ~ X, data = sp_df,
              fixed.vars = c("X"),
              bw = bw,
              kernel = "bisquare", adaptive = TRUE),
    error = function(e) NULL
  )
  
  # 추출
  beta_std_surface <- if (!is.null(fit_std)) fit_std$SDF$X else NA
  
  # gwr.mixed 출력 구조:
  #   fixed(전역) 계수 → fit_mix$fixed  (named vector)
  #   varying(국지) 계수 → fit_mix$SDF  (SpatialPointsDataFrame)
  beta_mixed_global <- if (!is.null(fit_mix)) {
    fx <- fit_mix$fixed
    if (!is.null(fx) && "X" %in% names(fx)) as.numeric(fx["X"]) else NA
  } else NA

  # Mixed GWR의 국지 절편이 Z를 따라가는지 진단
  intercept_mixed <- if (!is.null(fit_mix)) {
    int_col <- grep("Intercept", names(fit_mix$SDF), value = TRUE, ignore.case = TRUE)
    if (length(int_col) > 0) fit_mix$SDF[[int_col[1]]] else NA
  } else NA
  
  cor_int_Z <- if (length(intercept_mixed) > 1 && !any(is.na(intercept_mixed))) {
    cor(intercept_mixed, df$Z)
  } else NA
  
  data.frame(
    bw = bw,
    std_mean       = if (length(beta_std_surface) > 1) mean(beta_std_surface) else NA,
    std_sd         = if (length(beta_std_surface) > 1) sd(beta_std_surface) else NA,
    std_rmse       = if (length(beta_std_surface) > 1)
                       sqrt(mean((beta_std_surface - cfg$beta_true)^2)) else NA,
    mixed_beta     = beta_mixed_global,
    mixed_bias     = beta_mixed_global - cfg$beta_true,
    cor_int_Z      = cor_int_Z
  )
}

# ---- 4. Monte Carlo 반복 -----------------------------------------------------

cat("[1/3] Generating grid and running simulations...\n")
s_seq  <- seq(0, CFG$domain_size, length.out = CFG$n_side)
coords <- expand.grid(u = s_seq, v = s_seq)

run_one_rep <- function(rep_id, cfg) {
  cat(sprintf("  rep %d/%d\n", rep_id, cfg$n_reps))
  df <- generate_dgp_B(coords, cfg)
  sp_df <- SpatialPointsDataFrame(coords = df[, c("u", "v")], data = df)
  
  # OLS 참조
  ols_beta <- coef(lm(y ~ X, data = df))["X"]
  
  # 모든 bw에서 두 모형 적합
  res <- map_dfr(cfg$bw_grid, function(bw) fit_both(sp_df, df, bw, cfg))
  res$ols_beta <- ols_beta
  res$rep      <- rep_id
  res
}

all_results <- map_dfr(seq_len(CFG$n_reps), ~ run_one_rep(.x, CFG))

# ---- 5. 집계 -----------------------------------------------------------------

cat("[2/3] Aggregating...\n")

agg <- all_results %>%
  group_by(bw) %>%
  summarize(
    std_mean_med   = median(std_mean, na.rm = TRUE),
    std_sd_med     = median(std_sd, na.rm = TRUE),
    std_rmse_med   = median(std_rmse, na.rm = TRUE),
    mixed_bias_med = median(mixed_bias, na.rm = TRUE),
    mixed_beta_med = median(mixed_beta, na.rm = TRUE),
    cor_int_Z_med  = median(cor_int_Z, na.rm = TRUE),
    ols_beta_med   = median(ols_beta, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    std_bias_med   = std_mean_med - CFG$beta_true,
    ols_bias       = ols_beta_med - CFG$beta_true
  )

cat("\n=== 집계 결과 ===\n")
print(agg %>% select(bw, std_bias_med, std_sd_med, mixed_bias_med, cor_int_Z_med, ols_bias))

# ---- 6. 시각화 --------------------------------------------------------------

cat("[3/3] Plotting...\n")

theme_clean <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

# Panel 1: bias 비교 — Standard mean vs Mixed
p_bias <- agg %>%
  select(bw, `Standard GWR (mean β̂(s))` = std_bias_med,
              `Mixed GWR (global β̂)` = mixed_bias_med,
              OLS = ols_bias) %>%
  pivot_longer(-bw, names_to = "model", values_to = "bias") %>%
  ggplot(aes(x = bw, y = bias, color = model, group = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  scale_x_log10() +
  scale_color_manual(values = c("Standard GWR (mean β̂(s))" = "#d95f02",
                                "Mixed GWR (global β̂)" = "#1b9e77",
                                "OLS" = "gray50")) +
  labs(title = "Bias of estimated β across models",
       subtitle = sprintf("진실 β = %.1f, ℓ_Z = %.1f, ℓ_X = %.1f",
                          CFG$beta_true, CFG$scale_Z, CFG$scale_X),
       x = "Bandwidth (log scale)", y = "Bias (β̂ − β_true)",
       color = NULL) +
  theme_clean

# Panel 2: Standard GWR의 sd(β̂(s)) — 거짓 변동
p_std_sd <- agg %>%
  ggplot(aes(x = bw, y = std_sd_med)) +
  geom_line(color = "#d95f02", linewidth = 0.9) +
  geom_point(color = "#d95f02", size = 2.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_x_log10() +
  labs(title = "Standard GWR의 거짓 공간 변동",
       subtitle = "진실 sd = 0 (β는 상수). Mixed GWR은 정의상 sd = 0",
       x = "Bandwidth (log scale)",
       y = expression(sd[s](hat(beta)(s)))) +
  theme_clean

# Panel 3: Mixed GWR 절편이 Z를 따라가는가 (메커니즘 진단)
p_int_Z <- agg %>%
  ggplot(aes(x = bw, y = cor_int_Z_med)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 1, linetype = "dotted", color = "gray60") +
  geom_line(color = "#1b9e77", linewidth = 0.9) +
  geom_point(color = "#1b9e77", size = 2.8) +
  scale_x_log10() +
  ylim(c(-0.2, 1.05)) +
  labs(title = "Mixed GWR의 국지 절편 α(s)가 Z(s)를 따라가는가",
       subtitle = "1.0에 가까우면 → 직관 검증 (절편이 Z를 흡수)",
       x = "Bandwidth (log scale)",
       y = expression(cor(hat(alpha)(s), Z(s)))) +
  theme_clean

# Panel 4: 한 줄 요약 표
summary_df <- agg %>%
  select(bw,
         `Std. mean β̂` = std_mean_med,
         `Std. sd β̂`   = std_sd_med,
         `Mixed β̂`     = mixed_beta_med,
         `cor(α̂, Z)`   = cor_int_Z_med) %>%
  mutate(across(-bw, ~ round(.x, 3)))

cat("\n=== 한 줄 요약 ===\n")
print(summary_df)

# 결합
final_plot <- (p_bias / p_std_sd / p_int_Z) +
  plot_annotation(
    title = "Standard GWR vs Mixed GWR — DGP B에서 OVB 보정 비교",
    subtitle = sprintf("Mixed GWR이 Jiho 직관의 깨끗한 구현 (β 전역 상수, 절편 국지) — γ=%.1f, α_XZ=%.1f",
                       CFG$gamma_true, CFG$alpha_XZ),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

print(final_plot)
ggsave("std_vs_mixed.png", final_plot, width = 10, height = 12, dpi = 150)

cat("\n=== Done ===\n")