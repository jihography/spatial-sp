# =============================================================================
# DGP B에서 GWR의 misspecification 시각적 시연
# -----------------------------------------------------------------------------
# 진실: y = beta * x + gamma * Z(s) + eps,  beta = 1.0 (상수)
# 적합: y ~ x via GWR  →  결과적으로 beta_hat(s)가 공간상 변동
#       이 변동은 진실과 무관, Z(s)의 그림자
# =============================================================================

library(GWmodel)
library(sp)
library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)

set.seed(20260508)

# ---- 1. 설정 -----------------------------------------------------------------

n_side       <- 30
domain_size  <- 10
scale_Z      <- 3.0      # 큰 스케일 confounder
scale_X_fine <- 0.6      # 작은 스케일 X 변동
alpha_XZ     <- 0.85      # X-Z 결합 강도
beta_true    <- 2.0      # 진짜 효과 (상수)
gamma_true   <- 5.0      # Z의 효과
sigma_eps    <- 0.2

# ---- 2. 공간 필드 생성 -------------------------------------------------------

sim_grf <- function(coords, range_param, sigma2 = 1) {
  d <- as.matrix(dist(coords))
  Sigma <- sigma2 * exp(-d / range_param) + 1e-6 * diag(nrow(d))
  L <- chol(Sigma)
  as.numeric(t(L) %*% rnorm(nrow(d)))
}

s <- seq(0, domain_size, length.out = n_side)
coords <- expand.grid(u = s, v = s)
n <- nrow(coords)

Z      <- sim_grf(coords, scale_Z)
X_fine <- sim_grf(coords, scale_X_fine)
X      <- alpha_XZ * Z + X_fine

# DGP B: 상수 beta + 공간 confounder
y <- beta_true * X + gamma_true * Z + rnorm(n, sd = sigma_eps)

df <- data.frame(coords, X = X, Z = Z, y = y)
sp_df <- SpatialPointsDataFrame(coords = coords, data = df)

# ---- 3. GWR 적합 (데이터 기반 최적 대역폭) -----------------------------------

cat("Selecting optimal bandwidth (AICc)...\n")
opt_bw <- bw.gwr(y ~ X, data = sp_df, kernel = "bisquare",
                 adaptive = TRUE, approach = "AICc")
cat(sprintf("Optimal adaptive bandwidth: %d nearest neighbors\n", opt_bw))

fit <- gwr.basic(y ~ X, data = sp_df, bw = opt_bw,
                 kernel = "bisquare", adaptive = TRUE)

beta_hat <- fit$SDF$X
df$beta_hat <- beta_hat

# ---- 4. 진단 통계량 ----------------------------------------------------------

ols <- lm(y ~ X, data = df)

cat("\n--- Estimation summary ---\n")
cat(sprintf("True beta:                   %.3f (constant)\n", beta_true))
cat(sprintf("OLS beta_hat:                %.3f (biased due to omitted Z)\n",
            coef(ols)["X"]))
cat(sprintf("GWR mean(beta_hat(s)):       %.3f\n", mean(beta_hat)))
cat(sprintf("GWR sd(beta_hat(s)):         %.3f  (truth: 0)\n", sd(beta_hat)))
cat(sprintf("cor(beta_hat(s), Z(s)):      %+.3f  (truth: 0; if non-zero, beta_hat tracks confounder)\n",
            cor(beta_hat, Z)))
cat(sprintf("RMSE of beta_hat vs truth:   %.3f\n",
            sqrt(mean((beta_hat - beta_true)^2))))

# ---- 5. 시각화 --------------------------------------------------------------

theme_map <- theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        legend.position = "right")

# (a) Z(s) 지도
p_Z <- ggplot(df, aes(u, v, fill = Z)) +
  geom_raster() + coord_equal() +
  scale_fill_viridis_c(name = "Z(s)") +
  labs(title = "Z(s) — 미관측 confounder (진실)") +
  theme_map

# (b) GWR이 추정한 beta_hat(s) 지도
p_beta_hat <- ggplot(df, aes(u, v, fill = beta_hat)) +
  geom_raster() + coord_equal() +
  scale_fill_viridis_c(name = expression(hat(beta)(s)), option = "plasma") +
  labs(title = expression(hat(beta)(s) ~ "— GWR 추정 (진실: 상수 1.0)")) +
  theme_map

# (c) beta_hat(s) 분포 vs 진실
p_dist <- ggplot(df, aes(x = beta_hat)) +
  geom_histogram(bins = 30, fill = "#d95f02", alpha = 0.7) +
  geom_vline(xintercept = beta_true, linetype = "dashed",
             color = "black", linewidth = 0.8) +
  geom_vline(xintercept = mean(beta_hat), linetype = "dotted",
             color = "#d95f02", linewidth = 0.8) +
  annotate("text", x = beta_true, y = Inf, label = " 진실 β = 1.0",
           hjust = 0, vjust = 1.5, size = 3.5) +
  annotate("text", x = mean(beta_hat), y = Inf,
           label = sprintf(" mean(β̂) = %.2f", mean(beta_hat)),
           hjust = 0, vjust = 3.5, size = 3.5, color = "#d95f02") +
  labs(title = expression(hat(beta)(s) ~ "분포 — 진실은 점선"),
       x = expression(hat(beta)(s)), y = "frequency") +
  theme_minimal(base_size = 11)

# (d) beta_hat(s) vs Z(s) 산점도 — "추정된 효과 변동이 사실 confounder를 따라감"
p_scatter <- ggplot(df, aes(x = Z, y = beta_hat)) +
  geom_point(alpha = 0.5, color = "#d95f02") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
  geom_hline(yintercept = beta_true, linetype = "dashed", color = "gray40") +
  labs(title = expression("산점도: " ~ hat(beta)(s) ~ "vs" ~ Z(s)),
       subtitle = sprintf("cor = %+.3f  (진실대로라면 0이어야 함)",
                          cor(beta_hat, Z)),
       x = "Z(s)", y = expression(hat(beta)(s))) +
  theme_minimal(base_size = 11)

# 결합
final_plot <- (p_Z | p_beta_hat) / (p_dist | p_scatter) +
  plot_annotation(
    title = "DGP B에서 GWR의 misspecification",
    subtitle = sprintf("진실: β = %.1f (상수). GWR: 변동하는 β̂(s) — 그러나 그 변동은 Z의 그림자",
                       beta_true),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

print(final_plot)
ggsave("dgp_B_misspecification.png", final_plot, width = 11, height = 9, dpi = 150)
