# =============================================================================
# DGP B: Z의 공간 구조가 GWR 편향에 미치는 체계적 영향
# -----------------------------------------------------------------------------
# 실험 설계:
#   Exp 1. scale_Z 변화  (range 효과 — 비단조성 검증)
#   Exp 2. gamma_true 변화 (효과 크기 — 선형 증폭 검증)
#   Exp 3. alpha_XZ 변화  (X–Z 결합 강도)
#
# 진단 지표 (각 조건마다):
#   sd(beta_hat)        : GWR 계수의 거짓 공간 변동
#   cor(beta_hat, Z)    : 계수가 confounder를 얼마나 반영하는지
#   RMSE vs beta_true   : 절대적 추정 오차
#   |OLS_bias|          : OLS 편향 (비교 기준)
# =============================================================================

library(GWmodel)
library(sp)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

set.seed(20260508)

# ---- 공통 설정 ---------------------------------------------------------------

N_SIDE      <- 20          # 격자 크기 (20×20=400, 속도 위해 30에서 축소)
DOMAIN_SIZE <- 10
BETA_TRUE   <- 1.0
SIGMA_EPS   <- 0.3

# 기준 파라미터 (각 실험에서 나머지는 고정)
BASE_SCALE_Z  <- 2.0
BASE_GAMMA    <- 3.0
BASE_ALPHA_XZ <- 0.7
BASE_SCALE_X_FINE <- 0.5

# ---- 헬퍼 함수 ---------------------------------------------------------------

sim_grf <- function(coords, range_param, sigma2 = 1) {
  d     <- as.matrix(dist(coords))
  Sigma <- sigma2 * exp(-d / range_param) + 1e-6 * diag(nrow(d))
  L     <- chol(Sigma)
  as.numeric(t(L) %*% rnorm(nrow(d)))
}

run_one <- function(scale_Z   = BASE_SCALE_Z,
                    gamma     = BASE_GAMMA,
                    alpha_XZ  = BASE_ALPHA_XZ,
                    scale_X_fine = BASE_SCALE_X_FINE) {

  s      <- seq(0, DOMAIN_SIZE, length.out = N_SIDE)
  coords <- expand.grid(u = s, v = s)

  Z      <- sim_grf(coords, scale_Z)
  X_fine <- sim_grf(coords, scale_X_fine)
  X      <- alpha_XZ * Z + X_fine
  y      <- BETA_TRUE * X + gamma * Z + rnorm(nrow(coords), sd = SIGMA_EPS)

  sp_df  <- SpatialPointsDataFrame(coords = coords,
                                   data   = data.frame(X = X, y = y, Z = Z))

  # 대역폭 선택
  bw <- bw.gwr(y ~ X, data = sp_df, kernel = "bisquare",
               adaptive = TRUE, approach = "AICc")

  fit       <- gwr.basic(y ~ X, data = sp_df, bw = bw,
                         kernel = "bisquare", adaptive = TRUE)
  beta_hat  <- fit$SDF$X
  ols_beta  <- coef(lm(y ~ X, data = as.data.frame(sp_df@data)))["X"]

  list(
    sd_beta    = sd(beta_hat),
    cor_Z      = cor(beta_hat, Z),
    rmse       = sqrt(mean((beta_hat - BETA_TRUE)^2)),
    ols_bias   = abs(ols_beta - BETA_TRUE),
    opt_bw     = bw
  )
}

# ---- Exp 1: scale_Z 변화 (range 효과) ----------------------------------------

cat("=== Exp 1: Varying scale_Z ===\n")

scale_Z_grid <- c(0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0)
res_scaleZ   <- lapply(scale_Z_grid, function(sz) {
  cat(sprintf("  scale_Z = %.1f ... ", sz))
  r <- run_one(scale_Z = sz)
  cat(sprintf("sd=%.3f  cor=%.3f  bw=%d\n", r$sd_beta, r$cor_Z, r$opt_bw))
  data.frame(scale_Z   = sz,
             sd_beta   = r$sd_beta,
             cor_Z     = r$cor_Z,
             rmse      = r$rmse,
             ols_bias  = r$ols_bias,
             opt_bw    = r$opt_bw)
})
df_scaleZ <- bind_rows(res_scaleZ)

# ---- Exp 2: gamma_true 변화 --------------------------------------------------

cat("\n=== Exp 2: Varying gamma ===\n")

gamma_grid <- c(0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0)
res_gamma  <- lapply(gamma_grid, function(g) {
  cat(sprintf("  gamma = %.1f ... ", g))
  r <- run_one(gamma = g)
  cat(sprintf("sd=%.3f  cor=%.3f\n", r$sd_beta, r$cor_Z))
  data.frame(gamma    = g,
             sd_beta  = r$sd_beta,
             cor_Z    = r$cor_Z,
             rmse     = r$rmse,
             ols_bias = r$ols_bias)
})
df_gamma <- bind_rows(res_gamma)

# ---- Exp 3: alpha_XZ 변화 ----------------------------------------------------

cat("\n=== Exp 3: Varying alpha_XZ ===\n")

alpha_grid <- c(0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
res_alpha  <- lapply(alpha_grid, function(a) {
  cat(sprintf("  alpha_XZ = %.2f ... ", a))
  r <- run_one(alpha_XZ = a)
  cat(sprintf("sd=%.3f  cor=%.3f\n", r$sd_beta, r$cor_Z))
  data.frame(alpha_XZ = a,
             sd_beta  = r$sd_beta,
             cor_Z    = r$cor_Z,
             rmse     = r$rmse,
             ols_bias = r$ols_bias)
})
df_alpha <- bind_rows(res_alpha)

# ---- 요약 출력 ----------------------------------------------------------------

cat("\n--- Exp 1 summary (scale_Z sensitivity) ---\n")
print(df_scaleZ)
cat("\n--- Exp 2 summary (gamma sensitivity) ---\n")
print(df_gamma)
cat("\n--- Exp 3 summary (alpha_XZ sensitivity) ---\n")
print(df_alpha)

# ---- 시각화 ------------------------------------------------------------------

theme_pub <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

# (1) scale_Z — 비단조성 및 bandwidth와의 관계
p1a <- ggplot(df_scaleZ, aes(x = scale_Z)) +
  geom_line(aes(y = sd_beta, color = "sd(β̂)"), linewidth = 0.9) +
  geom_point(aes(y = sd_beta, color = "sd(β̂)"), size = 2.5) +
  geom_line(aes(y = rmse, color = "RMSE"), linewidth = 0.9, linetype = "dashed") +
  geom_point(aes(y = rmse, color = "RMSE"), size = 2.5) +
  scale_color_manual(values = c("sd(β̂)" = "#d95f02", "RMSE" = "#7570b3")) +
  labs(title = "Exp 1: scale_Z → GWR 편향 지표",
       x = "Z의 공간 범위 (scale_Z)", y = "값",
       color = NULL) +
  theme_pub

p1b <- ggplot(df_scaleZ, aes(x = scale_Z, y = cor_Z)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(color = "#1b9e77", linewidth = 0.9) +
  geom_point(color = "#1b9e77", size = 2.5) +
  labs(title = "Exp 1: cor(β̂, Z) — 계수가 confounder를 따라가는 정도",
       x = "Z의 공간 범위 (scale_Z)", y = "cor(β̂(s), Z(s))") +
  theme_pub

# bandwidth도 함께 — scale_Z와 opt_bw의 관계
p1c <- ggplot(df_scaleZ, aes(x = scale_Z, y = opt_bw)) +
  geom_line(color = "#666666", linewidth = 0.9) +
  geom_point(color = "#666666", size = 2.5) +
  labs(title = "Exp 1: AICc 최적 대역폭 vs scale_Z",
       subtitle = "scale_Z ≈ bandwidth일 때 편향 극대화 예상",
       x = "Z의 공간 범위 (scale_Z)", y = "최적 adaptive BW (# 이웃)") +
  theme_pub

# (2) gamma — 선형 증폭
p2 <- ggplot(df_gamma, aes(x = gamma)) +
  geom_line(aes(y = sd_beta, color = "GWR sd(β̂)"), linewidth = 0.9) +
  geom_point(aes(y = sd_beta, color = "GWR sd(β̂)"), size = 2.5) +
  geom_line(aes(y = ols_bias, color = "|OLS bias|"), linewidth = 0.9, linetype = "dashed") +
  geom_point(aes(y = ols_bias, color = "|OLS bias|"), size = 2.5) +
  scale_color_manual(values = c("GWR sd(β̂)" = "#d95f02", "|OLS bias|" = "#377eb8")) +
  labs(title = "Exp 2: γ 효과 크기 → 편향 증폭",
       subtitle = "선형 증폭이면 β̂ 분산이 γ에 비례",
       x = "γ (Z의 효과 크기)", y = "값",
       color = NULL) +
  theme_pub

# (3) alpha_XZ — 단조 증가
p3 <- ggplot(df_alpha, aes(x = alpha_XZ)) +
  geom_line(aes(y = sd_beta, color = "GWR sd(β̂)"), linewidth = 0.9) +
  geom_point(aes(y = sd_beta, color = "GWR sd(β̂)"), size = 2.5) +
  geom_line(aes(y = cor_Z, color = "cor(β̂, Z)"), linewidth = 0.9, linetype = "dashed") +
  geom_point(aes(y = cor_Z, color = "cor(β̂, Z)"), size = 2.5) +
  scale_color_manual(values = c("GWR sd(β̂)" = "#d95f02", "cor(β̂, Z)" = "#1b9e77")) +
  labs(title = "Exp 3: X–Z 결합 강도 α → 편향",
       subtitle = "α → 1 : X와 Z 거의 공선, 분리 불가",
       x = "α (X–Z 결합 강도)", y = "값",
       color = NULL) +
  theme_pub

# 결합
final_plot <- (p1a | p1b) / (p1c | p2) / (p3 + plot_spacer()) +
  plot_annotation(
    title = "DGP B: Z의 공간 구조 → GWR 편향 민감도 분석",
    subtitle = sprintf(
      "격자 %d×%d, β_true=%.1f, 기준 파라미터: scale_Z=%.1f, γ=%.1f, α=%.1f",
      N_SIDE, N_SIDE, BETA_TRUE, BASE_SCALE_Z, BASE_GAMMA, BASE_ALPHA_XZ
    ),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

print(final_plot)
ggsave("dgp_B_sensitivity.png", final_plot, width = 12, height = 12, dpi = 150)
cat("\n저장 완료: dgp_B_sensitivity.png\n")

# ---- 핵심 수치 요약 -----------------------------------------------------------

cat("\n============================================================\n")
cat("핵심 결과 요약\n")
cat("------------------------------------------------------------\n")

cat("\n[Exp 1] scale_Z 효과:\n")
cat(sprintf("  sd(β̂) 최대: scale_Z=%.1f 에서 %.3f\n",
            df_scaleZ$scale_Z[which.max(df_scaleZ$sd_beta)],
            max(df_scaleZ$sd_beta)))
cat(sprintf("  cor(β̂,Z) 최대: scale_Z=%.1f 에서 %.3f\n",
            df_scaleZ$scale_Z[which.max(df_scaleZ$cor_Z)],
            max(df_scaleZ$cor_Z)))
cat(sprintf("  → 비단조성 여부: sd_beta 최솟값=%.3f (scale_Z=%.1f), 최댓값=%.3f (scale_Z=%.1f)\n",
            min(df_scaleZ$sd_beta), df_scaleZ$scale_Z[which.min(df_scaleZ$sd_beta)],
            max(df_scaleZ$sd_beta), df_scaleZ$scale_Z[which.max(df_scaleZ$sd_beta)]))

cat("\n[Exp 2] gamma 효과:\n")
lm_gamma <- lm(sd_beta ~ gamma, data = df_gamma)
cat(sprintf("  sd(β̂) ~ γ 선형 회귀: 기울기=%.4f (R²=%.3f)\n",
            coef(lm_gamma)["gamma"],
            summary(lm_gamma)$r.squared))

cat("\n[Exp 3] alpha_XZ 효과:\n")
cat(sprintf("  α=0.1 → sd(β̂)=%.3f,  α=0.95 → sd(β̂)=%.3f\n",
            df_alpha$sd_beta[df_alpha$alpha_XZ == 0.1],
            df_alpha$sd_beta[df_alpha$alpha_XZ == 0.95]))
cat("============================================================\n")