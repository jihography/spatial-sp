# ══════════════════════════════════════════════════════════════════
#  실증분석: 출산율 GWR — 누락변수(x7)의 영향
#
#  Y  : 합계출산율
#  X  : x1 — 아동 1000명당 어린이집 수  (관심 변수)
#  Z  : x7 — 여성 경제활동참가율 (%)    (누락변수)
#
#  분석 흐름:
#    1. 공간 가중행렬 W 구성
#    2. Moran's I — x1, x7의 공간 스케일 비교
#    3. Local-L   — x1-x7 국지 이변량 자기상관
#    4. Short GWR : y ~ x1
#    5. Long GWR  : y ~ x1 + x7
#    6. β̂_x1 지도 비교 (short vs long) + Local-L 지도
#
#  패키지: sp, GWmodel, spdep, ggplot2, dplyr, patchwork, viridis
# ══════════════════════════════════════════════════════════════════

library(sp); library(GWmodel); library(spdep)
library(ggplot2); library(dplyr); library(patchwork); library(viridis)

# ──────────────────────────────────────────────────────────────────
# 0. 데이터 로드
# ──────────────────────────────────────────────────────────────────
df <- read.csv("sgwr_data.csv")   # Longitude, Latitude, y, x1~x7
n  <- nrow(df)
cat(sprintf("n = %d 시군구\n", n))

coords <- as.matrix(df[, c("Longitude", "Latitude")])

# sp 객체 (GWmodel용)
dat_sp <- df
coordinates(dat_sp) <- ~Longitude + Latitude

# ──────────────────────────────────────────────────────────────────
# 1. 공간 가중행렬 W (k=8 nearest neighbor, row-standardized)
# ──────────────────────────────────────────────────────────────────
knn8  <- knearneigh(coords, k = 8)
nb8   <- knn2nb(knn8)
W_lst <- nb2listw(nb8, style = "W")   # listw for spdep

# ──────────────────────────────────────────────────────────────────
# 2. Moran's I — 각 변수의 공간 자기상관
#    값이 클수록 공간 클러스터링이 강함 → 큰 스케일로 변동
# ──────────────────────────────────────────────────────────────────
vars_of_interest <- c("y", "x1", "x7")
moran_res <- sapply(vars_of_interest, function(v) {
  mt <- moran.test(df[[v]], W_lst)
  c(I = round(mt$estimate["Moran I statistic"], 4),
    p = round(mt$p.value, 4))
})
cat("\n── Moran's I ──────────────────────────────\n")
print(t(moran_res))

# ──────────────────────────────────────────────────────────────────
# 3. GWR Bandwidth 탐색
#    Short BW (y ~ x1) / Long BW (y ~ x1 + x7) 각각 탐색
# ──────────────────────────────────────────────────────────────────
cat("\nShort GWR BW 탐색 중...\n")
bw_short <- bw.gwr(y ~ x1, data = dat_sp,
                   kernel = "bisquare", adaptive = TRUE, approach = "AICc")
cat(sprintf("Short BW = %d\n", bw_short))

cat("Long GWR BW 탐색 중...\n")
bw_long  <- bw.gwr(y ~ x1 + x7, data = dat_sp,
                   kernel = "bisquare", adaptive = TRUE, approach = "AICc")
cat(sprintf("Long  BW = %d\n", bw_long))

# ──────────────────────────────────────────────────────────────────
# 4. Short GWR: y ~ x1  (x7 누락)
# ──────────────────────────────────────────────────────────────────
gwr_short <- gwr.basic(y ~ x1, data = dat_sp,
                       bw = bw_short, kernel = "bisquare", adaptive = TRUE)
cat("\n── Short GWR (y ~ x1) ─────────────────────\n")
print(gwr_short)

# ──────────────────────────────────────────────────────────────────
# 5. Long GWR: y ~ x1 + x7  (x7 포함)
# ──────────────────────────────────────────────────────────────────
gwr_long  <- gwr.basic(y ~ x1 + x7, data = dat_sp,
                       bw = bw_long,  kernel = "bisquare", adaptive = TRUE)
cat("\n── Long GWR (y ~ x1 + x7) ─────────────────\n")
print(gwr_long)

# ──────────────────────────────────────────────────────────────────
# 6. Local-L: x1-x7 국지 이변량 자기상관
#    GWR short BW 재활용
# ──────────────────────────────────────────────────────────────────
x1_std <- scale(df$x1)[, 1]
x7_std <- scale(df$x7)[, 1]

local_L <- vapply(seq_len(n), function(i) {
  d <- sqrt(rowSums(sweep(coords, 2, coords[i, ], "-")^2))
  h <- sort(d)[bw_short]
  w <- (1 - (d / h)^2)^2 * (d < h)
  w[i] <- 0
  if (sum(w) == 0) return(NA_real_)

  xw  <- sum(w * x1_std) / sum(w)
  zw  <- sum(w * x7_std) / sum(w)
  x_c <- x1_std - xw
  z_c <- x7_std - zw

  num   <- sum(w * x_c * z_c)
  denom <- sqrt(sum(w * x_c^2) * sum(w * z_c^2))
  if (denom > 0) num / denom else NA_real_
}, numeric(1))

cat(sprintf("\nLocal-L(x1, x7): mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
            mean(local_L, na.rm=TRUE), sd(local_L, na.rm=TRUE),
            min(local_L, na.rm=TRUE), max(local_L, na.rm=TRUE)))

# ──────────────────────────────────────────────────────────────────
# 7. 결과 데이터프레임 통합
# ──────────────────────────────────────────────────────────────────
res <- data.frame(
  lon         = df$Longitude,
  lat         = df$Latitude,
  y           = df$y,
  x1          = df$x1,
  x7          = df$x7,
  beta_x1_short = as.numeric(gwr_short$SDF$x1),
  beta_x1_long  = as.numeric(gwr_long$SDF$x1),
  beta_x7_long  = as.numeric(gwr_long$SDF$x7),
  local_L       = local_L
)
write.csv(res, "empirical_results.csv", row.names = FALSE)

# ──────────────────────────────────────────────────────────────────
# 8. 시각화
# ──────────────────────────────────────────────────────────────────

# 공통 테마
map_theme <- theme_minimal(base_size = 10) +
  theme(
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    panel.grid   = element_blank(),
    plot.title   = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle= element_text(size = 8.5, color = "grey40", hjust = 0.5),
    legend.position = "right"
  )

# 공통 스케일 (short/long 비교를 위해 동일 범위)
beta_lim <- max(abs(c(res$beta_x1_short, res$beta_x1_long))) * 1.05

# ── 8-1. Short GWR β̂_x1
p_short <- ggplot(res, aes(x = lon, y = lat, color = beta_x1_short)) +
  geom_point(size = 2.2, alpha = 0.9) +
  scale_color_gradient2(
    low = "#2166AC", mid = "white", high = "#D6604D",
    midpoint = 0, limits = c(-beta_lim, beta_lim),
    name = expression(hat(beta)[x1])
  ) +
  labs(
    title    = expression(bold("Short GWR")~~hat(beta)[x1]~~(y~"~"~x1)),
    subtitle = "x7 (여성경활률) 누락"
  ) +
  coord_fixed(ratio = 1.3) + map_theme

# ── 8-2. Long GWR β̂_x1
p_long <- ggplot(res, aes(x = lon, y = lat, color = beta_x1_long)) +
  geom_point(size = 2.2, alpha = 0.9) +
  scale_color_gradient2(
    low = "#2166AC", mid = "white", high = "#D6604D",
    midpoint = 0, limits = c(-beta_lim, beta_lim),
    name = expression(hat(beta)[x1])
  ) +
  labs(
    title    = expression(bold("Long GWR")~~hat(beta)[x1]~~(y~"~"~x1+x7)),
    subtitle = "x7 (여성경활률) 포함"
  ) +
  coord_fixed(ratio = 1.3) + map_theme

# ── 8-3. 차이 지도: Δβ̂_x1 = short - long
res$delta_beta <- res$beta_x1_short - res$beta_x1_long
delta_lim <- max(abs(res$delta_beta)) * 1.05

p_delta <- ggplot(res, aes(x = lon, y = lat, color = delta_beta)) +
  geom_point(size = 2.2, alpha = 0.9) +
  scale_color_gradient2(
    low = "#2166AC", mid = "white", high = "#D6604D",
    midpoint = 0, limits = c(-delta_lim, delta_lim),
    name = expression(Delta~hat(beta)[x1])
  ) +
  labs(
    title    = expression(bold("OVB 지도")~~Delta~hat(beta)[x1]),
    subtitle = "Short − Long  (누락변수로 인한 편향)"
  ) +
  coord_fixed(ratio = 1.3) + map_theme

# ── 8-4. Local-L 지도
ll_lim <- max(abs(res$local_L), na.rm = TRUE) * 1.05

p_ll <- ggplot(res, aes(x = lon, y = lat, color = local_L)) +
  geom_point(size = 2.2, alpha = 0.9) +
  scale_color_gradient2(
    low = "#2166AC", mid = "white", high = "#D6604D",
    midpoint = 0, limits = c(-ll_lim, ll_lim),
    name = "Local-L"
  ) +
  labs(
    title    = "Local-L (x1, x7)",
    subtitle = "커널 기반 국지 이변량 공간 자기상관"
  ) +
  coord_fixed(ratio = 1.3) + map_theme

# ── 8-5. Long GWR β̂_x7 지도
p_x7 <- ggplot(res, aes(x = lon, y = lat, color = beta_x7_long)) +
  geom_point(size = 2.2, alpha = 0.9) +
  scale_color_gradient2(
    low = "#2166AC", mid = "white", high = "#D6604D",
    midpoint = 0,
    name = expression(hat(beta)[x7])
  ) +
  labs(
    title    = expression(bold("Long GWR")~~hat(beta)[x7]),
    subtitle = "여성경활률의 국지 효과"
  ) +
  coord_fixed(ratio = 1.3) + map_theme

# ── 통합 저장
p_main <- (p_short | p_long | p_delta) /
          (p_ll    | p_x7  | plot_spacer()) +
  plot_annotation(
    title    = "출산율 GWR 실증분석  |  누락변수(여성경활률)의 효과",
    subtitle = expression(
      "Y = 합계출산율,  X = 어린이집 수 (x1),  Z = 여성경활률 (x7)"
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9.5, color = "grey40")
    )
  )

ggsave("fig_empirical_main.png", p_main, width = 16, height = 11, dpi = 150)
cat("저장: fig_empirical_main.png\n")

# ── Short vs Long 계수 산점도 (Δβ 크기 확인)
p_scatter <- ggplot(res, aes(x = beta_x1_long, y = beta_x1_short)) +
  geom_point(aes(color = local_L), size = 2, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey50") +
  scale_color_gradient2(
    low = "#2166AC", mid = "white", high = "#D6604D",
    midpoint = 0, name = "Local-L"
  ) +
  labs(
    title    = expression(hat(beta)[x1]~~": Short vs Long GWR"),
    subtitle = "대각선(점선) = bias 없음  /  이탈 정도 = OVB 크기",
    x = expression("Long GWR  "~hat(beta)[x1]~"  (x7 포함)"),
    y = expression("Short GWR  "~hat(beta)[x1]~"  (x7 누락)")
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave("fig_empirical_scatter.png", p_scatter,
       width = 7, height = 6, dpi = 150)
cat("저장: fig_empirical_scatter.png\n")

# empirical_analysis.R 실행 후
cat("Short GWR β̂_x1 std:", sd(gwr_short$SDF$x1), "\n")
cat("Long  GWR β̂_x1 std:", sd(gwr_long$SDF$x1),  "\n")
cat("감소율:", 1 - sd(gwr_long$SDF$x1)/sd(gwr_short$SDF$x1), "\n")

library(spdep)
lee_result <- lee(df$x1, df$x7, W_lst, n, zero.policy = TRUE)
cat("전역 Lee's L (x1, x7):", lee_result$L, "\n")