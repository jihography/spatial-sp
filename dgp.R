# =============================================================================
# Multi-Bandwidth GWR Diagnostic for Spatial Confounding vs Nonstationarity
# -----------------------------------------------------------------------------
# Purpose:
#   DGP A (true nonstationarity)와 DGP B (전역 상수 + 공간 교란요인)에서
#   각각 multi-bandwidth GWR을 적합하고, stability profile signature가
#   두 DGP를 식별 가능하게 분리하는지 검증.
#
# Author: 곽지호
# Date:   2026-05-08
# =============================================================================

# ---- 0. Libraries -----------------------------------------------------------

# install.packages(c("GWmodel", "sp", "fields", "ggplot2", "dplyr", 
#                    "tidyr", "purrr", "patchwork", "viridis"))

library(GWmodel)
library(sp)
library(fields)       # for image.plot, GRF utilities
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(patchwork)
library(viridis)

set.seed(20260508)

# ---- 1. Configuration -------------------------------------------------------

CONFIG <- list(
  # 공간 도메인
  n_side       = 30,           # 30x30 = 900 격자점
  domain_size  = 10,           # [0, 10] x [0, 10]
  
  # GRF 스케일 (range 모수 = exp covariance)
  scale_Z      = 3.0,          # 교란요인 스케일 (큼)
  scale_X_fine = 0.6,          # X의 미세 변동 스케일 (작음)
  scale_beta   = 2.0,          # DGP A에서 β(s)의 스케일
  
  # X = alpha_XZ * Z + fine_scale_noise  (X와 Z의 상관 제어)
  alpha_XZ     = 0.6,          # X-Z 결합 강도
  
  # 진짜 모수
  beta_true    = 1.0,          # DGP A: β(s) 평균; DGP B: β
  gamma_true   = 1.5,          # DGP B: Z 계수
  sigma_eps    = 0.3,          # 잔차 표준편차
  
  # GWR 대역폭 격자 (adaptive: nearest neighbors 수)
  bw_grid      = c(30, 60, 120, 240, 480, 900),
  
  # Monte Carlo 반복 (계산비용 고려: 처음엔 1, 검증 시 50+)
  n_replicates = 1
)

# ---- 2. Helper: Gaussian Random Field 시뮬레이션 -----------------------------

#' 지정된 좌표와 range 모수로 GRF 표본 생성 (지수 공분산)
#' 
#' @param coords matrix of (x, y) coordinates
#' @param range_param exp covariance의 range
#' @param sigma2 분산
#' @param nugget 수치 안정성을 위한 작은 nugget
sim_grf <- function(coords, range_param, sigma2 = 1, nugget = 1e-6) {
  d <- as.matrix(dist(coords))
  Sigma <- sigma2 * exp(-d / range_param) + nugget * diag(nrow(d))
  L <- chol(Sigma)
  as.numeric(t(L) %*% rnorm(nrow(d)))
}

# ---- 3. 공간 설계 및 데이터 생성 함수 ----------------------------------------

build_spatial_grid <- function(cfg) {
  s <- seq(0, cfg$domain_size, length.out = cfg$n_side)
  expand.grid(u = s, v = s)
}

#' 한 번의 시뮬레이션: DGP A와 DGP B를 동일한 X 분포 위에서 생성
#' 
#' 같은 X와 Z를 공유하는 두 결과변수 y_A, y_B를 만들어
#' GWR 적합 결과가 DGP에 의존함을 직접 비교 가능하게 함.
generate_data <- function(coords, cfg) {
  # 공간 잠재 필드들
  Z         <- sim_grf(coords, cfg$scale_Z)        # 큰 스케일 교란요인
  X_fine    <- sim_grf(coords, cfg$scale_X_fine)   # 작은 스케일 X 성분
  X         <- cfg$alpha_XZ * Z + X_fine           # X = αZ + fine
  beta_surf <- cfg$beta_true + sim_grf(coords, cfg$scale_beta) * 0.5  # DGP A
  
  # 결과변수
  eps_A <- rnorm(nrow(coords), sd = cfg$sigma_eps)
  eps_B <- rnorm(nrow(coords), sd = cfg$sigma_eps)
  
  # DGP A: 진짜 비정상성 (Z 효과 없음)
  y_A <- beta_surf * X + eps_A
  
  # DGP B: 전역 상수 + 공간 교란
  y_B <- cfg$beta_true * X + cfg$gamma_true * Z + eps_B
  
  data.frame(
    coords,
    X = X, Z = Z,
    beta_surf = beta_surf,
    y_A = y_A, y_B = y_B
  )
}

# ---- 4. GWR 적합 및 profile 추출 --------------------------------------------

#' 한 대역폭에서 GWR 적합 후 계수 표면 통계 추출
fit_gwr_profile <- function(formula, sp_df, bw, dgp_label, bw_label) {
  fit <- tryCatch(
    gwr.basic(formula, data = sp_df, bw = bw,
              kernel = "bisquare", adaptive = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(data.frame(
      dgp = dgp_label, bw = bw,
      mean_beta = NA, sd_beta = NA,
      median_beta = NA, q025 = NA, q975 = NA,
      moran_beta = NA, n_local = NA
    ))
  }
  
  # 계수 표면 추출
  beta_var <- as.character(formula[[3]])  # X (단일 공변량 가정)
  beta_hat <- fit$SDF[[beta_var]]
  
  data.frame(
    dgp = dgp_label, bw = bw,
    mean_beta   = mean(beta_hat, na.rm = TRUE),
    sd_beta     = sd(beta_hat, na.rm = TRUE),
    median_beta = median(beta_hat, na.rm = TRUE),
    q025        = quantile(beta_hat, 0.025, na.rm = TRUE),
    q975        = quantile(beta_hat, 0.975, na.rm = TRUE),
    moran_beta  = compute_moran_simple(beta_hat, sp_df),
    n_local     = length(beta_hat)
  )
}

#' 단순한 Moran's I (k-nearest neighbor 가중치)
compute_moran_simple <- function(z, sp_df, k = 8) {
  if (any(is.na(z))) return(NA)
  coords <- coordinates(sp_df)
  d <- as.matrix(dist(coords))
  W <- matrix(0, nrow(d), ncol(d))
  for (i in 1:nrow(d)) {
    nn <- order(d[i, ])[2:(k + 1)]  # 자기 자신 제외
    W[i, nn] <- 1 / k
  }
  z_c <- z - mean(z)
  num <- sum(W * outer(z_c, z_c))
  den <- sum(z_c^2)
  (length(z) / sum(W)) * (num / den)
}

#' 한 시뮬레이션 데이터셋에서 모든 (DGP, 대역폭) 조합 적합
run_full_profile <- function(df, cfg) {
  sp_df <- SpatialPointsDataFrame(
    coords = df[, c("u", "v")], data = df
  )
  
  # DGP A profile
  profile_A <- map_dfr(cfg$bw_grid, function(bw) {
    fit_gwr_profile(y_A ~ X, sp_df, bw, "A", as.character(bw))
  })
  
  # DGP B profile
  profile_B <- map_dfr(cfg$bw_grid, function(bw) {
    fit_gwr_profile(y_B ~ X, sp_df, bw, "B", as.character(bw))
  })
  
  bind_rows(profile_A, profile_B)
}

# ---- 5. 단일 시뮬레이션 실행 -------------------------------------------------

cat("[1/4] Building spatial grid...\n")
coords <- build_spatial_grid(CONFIG)

cat("[2/4] Generating data (DGP A & DGP B)...\n")
df <- generate_data(coords, CONFIG)

# OLS 참조
ols_A <- lm(y_A ~ X, data = df)
ols_B <- lm(y_B ~ X, data = df)

cat("OLS estimates:\n")
cat(sprintf("  DGP A: beta_hat = %.3f (true mean = %.3f)\n",
            coef(ols_A)["X"], mean(df$beta_surf)))
cat(sprintf("  DGP B: beta_hat = %.3f (true beta = %.3f, expected bias > 0)\n",
            coef(ols_B)["X"], CONFIG$beta_true))

cat("[3/4] Fitting GWR across bandwidth grid...\n")
profile <- run_full_profile(df, CONFIG)
print(profile)

# ---- 6. 시각화 --------------------------------------------------------------

cat("[4/4] Generating diagnostic plots...\n")

theme_clean <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

# Panel 1: mean(β̂) profile
p_mean <- profile %>%
  ggplot(aes(x = bw, y = mean_beta, color = dgp, group = dgp)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = CONFIG$beta_true,
             linetype = "dashed", color = "black", alpha = 0.6) +
  geom_hline(yintercept = coef(ols_A)["X"],
             linetype = "dotted", color = "#1b9e77") +
  geom_hline(yintercept = coef(ols_B)["X"],
             linetype = "dotted", color = "#d95f02") +
  scale_x_log10() +
  scale_color_manual(values = c("A" = "#1b9e77", "B" = "#d95f02"),
                     labels = c("A: nonstationarity", "B: confounding")) +
  labs(title = "Mean of GWR coefficient across bandwidth",
       subtitle = "Dashed = true β; Dotted = OLS estimate per DGP",
       x = "Bandwidth (k-nearest neighbors, log scale)",
       y = expression(bar(hat(beta))(h)),
       color = "DGP") +
  theme_clean

# Panel 2: sd(β̂) profile - 가장 진단적인 그림
p_sd <- profile %>%
  ggplot(aes(x = bw, y = sd_beta, color = dgp, group = dgp)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_log10() +
  scale_color_manual(values = c("A" = "#1b9e77", "B" = "#d95f02"),
                     labels = c("A: nonstationarity", "B: confounding")) +
  labs(title = "Spatial sd of GWR coefficient — diagnostic signature",
       subtitle = "DGP A: high at small h (true variation revealed)\nDGP B: low at small h (true β is constant)",
       x = "Bandwidth (k-nearest neighbors, log scale)",
       y = expression(sd[s](hat(beta)(s, h))),
       color = "DGP") +
  theme_clean

# Panel 3: 95% 분포 폭 (mean ± 2sd 대신 quantile)
p_range <- profile %>%
  ggplot(aes(x = bw, color = dgp, fill = dgp, group = dgp)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.2, color = NA) +
  geom_line(aes(y = median_beta), linewidth = 0.8) +
  geom_point(aes(y = median_beta), size = 2.5) +
  geom_hline(yintercept = CONFIG$beta_true, linetype = "dashed", alpha = 0.6) +
  scale_x_log10() +
  scale_color_manual(values = c("A" = "#1b9e77", "B" = "#d95f02")) +
  scale_fill_manual(values = c("A" = "#1b9e77", "B" = "#d95f02")) +
  labs(title = "Median and 95% spatial range of β̂(s; h)",
       x = "Bandwidth (log scale)",
       y = expression(hat(beta)(s, h)),
       color = "DGP", fill = "DGP") +
  theme_clean

# Panel 4: 입력 데이터 시각화
p_X <- ggplot(df, aes(u, v, fill = X)) +
  geom_raster() + scale_fill_viridis_c() +
  coord_equal() + labs(title = "X(s)") + theme_clean +
  theme(legend.position = "right")

p_Z <- ggplot(df, aes(u, v, fill = Z)) +
  geom_raster() + scale_fill_viridis_c() +
  coord_equal() + labs(title = "Z(s) — confounder") + theme_clean +
  theme(legend.position = "right")

p_beta_surf <- ggplot(df, aes(u, v, fill = beta_surf)) +
  geom_raster() + scale_fill_viridis_c(option = "plasma") +
  coord_equal() + labs(title = expression(beta(s) ~ "(DGP A truth)")) +
  theme_clean + theme(legend.position = "right")

# 결합 출력
diagnostic_plot <- (p_X | p_Z | p_beta_surf) /
                    (p_mean | p_sd) /
                    p_range +
  plot_annotation(
    title = "Multi-Bandwidth GWR Diagnostic: DGP A vs DGP B",
    subtitle = sprintf("scale(Z) = %.1f, scale(X_fine) = %.1f, alpha_XZ = %.1f, gamma = %.1f",
                       CONFIG$scale_Z, CONFIG$scale_X_fine,
                       CONFIG$alpha_XZ, CONFIG$gamma_true),
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

print(diagnostic_plot)
# ggsave("multi_bw_diagnostic.png", diagnostic_plot, 
#        width = 12, height = 14, dpi = 150)

# ---- 7. 다중 반복 시뮬레이션 (선택) ------------------------------------------

# Monte Carlo replicate를 통해 profile signature의 안정성 검증
# 처음에는 cfg$n_replicates = 1로 단일 결과 점검,
# 이후 50+회로 늘려 신뢰구간 확보

run_monte_carlo <- function(cfg) {
  coords <- build_spatial_grid(cfg)
  results <- map_dfr(seq_len(cfg$n_replicates), function(rep_id) {
    cat(sprintf("  Replicate %d / %d\n", rep_id, cfg$n_replicates))
    df_rep <- generate_data(coords, cfg)
    profile_rep <- run_full_profile(df_rep, cfg)
    profile_rep$rep <- rep_id
    profile_rep
  })
  results
}

# 다음과 같이 실행:
# CONFIG$n_replicates <- 50
# mc_results <- run_monte_carlo(CONFIG)
# 
# mc_results %>%
#   group_by(dgp, bw) %>%
#   summarize(mean_sd_beta = mean(sd_beta, na.rm = TRUE),
#             se_sd_beta   = sd(sd_beta, na.rm = TRUE) / sqrt(n()),
#             .groups = "drop") %>%
#   ggplot(aes(bw, mean_sd_beta, color = dgp)) +
#     geom_ribbon(aes(ymin = mean_sd_beta - 2*se_sd_beta,
#                     ymax = mean_sd_beta + 2*se_sd_beta,
#                     fill = dgp), alpha = 0.2, color = NA) +
#     geom_line() + geom_point() + scale_x_log10()

# ---- 8. 진단 결정 규칙 (탐색) ------------------------------------------------

#' DGP B (confounding) 진단 신호 점수
#' 
#' 작은 h에서 sd가 낮고 (안정), h가 커지며 sd가 증가하면 신호 증가.
#' 이는 stability profile의 핵심 부등식 점수화.
diagnostic_score <- function(profile_one_dgp) {
  prof <- profile_one_dgp %>% arrange(bw)
  sd_small <- prof$sd_beta[1]
  sd_large <- prof$sd_beta[nrow(prof)]
  # 음수면 DGP A 시그니처, 양수면 DGP B 시그니처
  log(sd_large / sd_small)
}

cat("\n--- Diagnostic scores (log sd_large / sd_small) ---\n")
profile %>%
  group_split(dgp) %>%
  walk(function(p) {
    s <- diagnostic_score(p)
    cat(sprintf("  DGP %s: score = %+.3f  (positive → confounding signature)\n",
                unique(p$dgp), s))
  })

# 해석:
#   score > 0: small h에서 sd 작고 large h에서 sd 큼 → DGP B (confounding)
#   score < 0: small h에서 sd 크고 large h에서 sd 작음 → DGP A (nonstationarity)

cat("\n=== Done ===\n")