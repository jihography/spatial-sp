# ══════════════════════════════════════════════════════════════════
#  Bivariate Coregionalization GP — 파라미터별 표면 시각화
#
#  같은 기본 노이즈(seed 고정)를 써서 파라미터 변화의 효과만 보여줌
#
#  Panel A: l_Z 변화 → Z, X 표면 (l_X=0.6, rho=0.5 고정)
#  Panel B: l_X 변화 → Z, X 표면 (l_Z=0.6, rho=0.5 고정)
#  Panel C: rho 변화 → Z, X 표면 (l_Z=0.6, l_X=0.2 고정)
# ══════════════════════════════════════════════════════════════════

library(ggplot2); library(dplyr); library(patchwork)

set.seed(42)

# ── 격자 설정 ─────────────────────────────────────────────────
N      <- 20
g1d    <- seq(0, 1, length.out = N)
grid   <- expand.grid(x = g1d, y = g1d)
coords <- as.matrix(grid)
n      <- nrow(coords)

# ── Matern(nu=1.5) 공분산 행렬 ───────────────────────────────
matern15 <- function(coords, l) {
  D <- as.matrix(dist(coords))
  r <- sqrt(3) * D / l
  K <- (1 + r) * exp(-r)
  diag(K) <- diag(K) + 1e-6
  K
}

# ── 기본 노이즈 고정 (한 번만 생성) ──────────────────────────
xi_Z  <- rnorm(n)
xi_X0 <- rnorm(n)

# ── Z, X 생성 함수 ───────────────────────────────────────────
gen_ZX <- function(l_z, l_x, rho) {
  Kz <- matern15(coords, l_z)
  Kx <- matern15(coords, l_x)
  Lz <- t(chol(Kz))   # lower Cholesky
  Lx <- t(chol(Kx))

  Z_raw <- as.vector(Lz %*% xi_Z)

  # X|Z: rho * Lx * Lz^{-1} * Z + sqrt(1-rho^2) * Lx * xi_X0
  Lz_inv_Z <- forwardsolve(Lz, Z_raw)
  X_raw <- rho * as.vector(Lx %*% Lz_inv_Z) +
           sqrt(max(1 - rho^2, 0)) * as.vector(Lx %*% xi_X0)

  # 표준화
  Z <- (Z_raw - mean(Z_raw)) / sd(Z_raw)
  X <- (X_raw - mean(X_raw)) / sd(X_raw)

  data.frame(x = coords[,1], y = coords[,2], Z = Z, X = X)
}

# ── 표면 단일 플롯 함수 ──────────────────────────────────────
surf_plot <- function(df, var, title, show_legend = FALSE) {
  cmap <- if (var == "Z") "BrBG" else "RdYlBu"
  ggplot(df, aes(x = x, y = y, fill = .data[[var]])) +
    geom_raster(interpolate = FALSE) +
    scale_fill_distiller(palette = cmap, limits = c(-2.5, 2.5),
                         direction = 1, name = var) +
    labs(title = title) +
    coord_equal() +
    theme_void(base_size = 9) +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 9,
                                     face = "bold"),
      legend.position = if (show_legend) "right" else "none"
    )
}

# ── Panel 빌더 ────────────────────────────────────────────────
build_panel <- function(param_vals, param_name,
                        fixed_lz, fixed_lx, fixed_rho,
                        panel_label) {
  plots <- list()
  for (pv in param_vals) {
    lz  <- if (param_name == "l_Z")  pv else fixed_lz
    lx  <- if (param_name == "l_X")  pv else fixed_lx
    rho <- if (param_name == "rho")  pv else fixed_rho

    df <- gen_ZX(lz, lx, rho)

    lab <- switch(param_name,
                  "l_Z"  = sprintf("l[Z]==%s", pv),
                  "l_X"  = sprintf("l[X]==%s", pv),
                  "rho"  = sprintf("rho==%s",  pv))

    is_last <- (pv == tail(param_vals, 1))
    p_Z <- surf_plot(df, "Z",
                     parse(text = paste0("Z~~(", lab, ")")))
    p_X <- surf_plot(df, "X",
                     parse(text = paste0("X~~(", lab, ")")),
                     show_legend = is_last)
    plots <- c(plots, list(p_Z, p_X))
  }

  # 3열 × 2행 배치 (열=파라미터, 행=Z/X)
  wrap_plots(plots, nrow = 2, ncol = length(param_vals),
             byrow = FALSE) +
    plot_annotation(
      title = panel_label,
      theme = theme(plot.title = element_text(face = "bold", size = 11))
    )
}

# ── 세 패널 생성 ──────────────────────────────────────────────
panelA <- build_panel(
  param_vals  = c(0.2, 0.6, 1.5),
  param_name  = "l_Z",
  fixed_lz    = 0.6,   # 미사용
  fixed_lx    = 0.6,
  fixed_rho   = 0.5,
  panel_label = "(A)  l_Z 변화  |  l_X=0.6, rho=0.5 고정"
)

panelB <- build_panel(
  param_vals  = c(0.2, 0.6, 1.5),
  param_name  = "l_X",
  fixed_lz    = 0.6,
  fixed_lx    = 0.6,   # 미사용
  fixed_rho   = 0.5,
  panel_label = "(B)  l_X 변화  |  l_Z=0.6, rho=0.5 고정"
)

panelC <- build_panel(
  param_vals  = c(0.0, 0.5, 0.9),
  param_name  = "rho",
  fixed_lz    = 0.6,
  fixed_lx    = 0.2,
  fixed_rho   = 0.5,   # 미사용
  panel_label = "(C)  rho 변화  |  l_Z=0.6, l_X=0.2 고정"
)

# ── 통합 저장 ─────────────────────────────────────────────────
p_all <- panelA / panelB / panelC +
  plot_annotation(
    title    = "Bivariate Coregionalization GP — 파라미터별 Z, X 표면",
    subtitle = "동일 seed(xi_Z, xi_X0) 사용: 파라미터 변화의 순수 효과",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "grey40")
    )
  )

ggsave("fig_gp_surfaces.png", p_all,
       width = 13, height = 15, dpi = 150)
cat("저장: fig_gp_surfaces.png\n")