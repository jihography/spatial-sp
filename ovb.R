# ══════════════════════════════════════════════════════════════════
#  GWR OVB Simulation  (R version, bivariate coregionalization GP)
#
#  DGP:
#    Z        ~ GP(0, K_Z)                    K_Z = Matern(l_Z)
#    X | Z    ~ GP(rho * L_X %*% solve(L_Z, Z),  (1-rho^2) * K_X)
#    y = beta*X + gamma*Z + eps,  beta=1 (constant), gamma=1
#
#  팩토리얼:
#    l_Z in {0.2, 0.6, 1.5}
#    l_X in {0.2, 0.6, 1.5}
#    rho in {0.0, 0.5, 0.9}
#
#  BW: 매 반복마다 adaptive AICc로 재탐색
#  DGP 검증: X, Z 각각의 Moran's I (rook adjacency, 행 표준화)
# ══════════════════════════════════════════════════════════════════

library(sp)
library(GWmodel)
library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
library(viridis)

dir.create("output_ovb", showWarnings = FALSE)
set.seed(2025)

# ──────────────────────────────────────────────────────────────────
# 0. 설정
# ──────────────────────────────────────────────────────────────────
N_GRID     <- 15
N_SIM      <- 30
BETA_TRUE  <- 1.0
GAMMA_TRUE <- 1.0
SIGMA_EPS  <- 0.3

L_Z_VALUES <- c(0.05, 0.2, 0.5)
L_X_VALUES <- c(0.05, 0.2, 0.5)
RHO_VALUES <- c(0.0, 0.4, 0.8)

# ──────────────────────────────────────────────────────────────────
# 1. 격자 좌표
# ──────────────────────────────────────────────────────────────────
grid_df <- expand.grid(x = seq(0, 1, length.out = N_GRID),
                       y = seq(0, 1, length.out = N_GRID))
coords  <- as.matrix(grid_df)
n       <- nrow(coords)

# ──────────────────────────────────────────────────────────────────
# 2. Moran's I 사전 설정 (rook adjacency, 행 표준화)
#    격자 간격 = 1/(N_GRID-1), 최소 거리 셀들이 이웃
# ──────────────────────────────────────────────────────────────────
D_mat  <- as.matrix(dist(coords))
min_d  <- min(D_mat[D_mat > 0])
W_rook <- (abs(D_mat - min_d) < 1e-10) * 1.0
diag(W_rook) <- 0
W_rook <- W_rook / rowSums(W_rook)   # 행 표준화
S0     <- sum(W_rook)                 # = n (행 표준화 후)

moran_i <- function(x) {
  xz <- x - mean(x)
  as.numeric(n / S0 * (t(xz) %*% W_rook %*% xz) / (t(xz) %*% xz))
}

# ──────────────────────────────────────────────────────────────────
# 3. Matern(nu=1.5) 공분산 행렬
# ──────────────────────────────────────────────────────────────────
matern15_cov <- function(coords, len_scale) {
  D <- as.matrix(dist(coords))
  r <- sqrt(3) * D / len_scale
  (1 + r) * exp(-r)
}

# ──────────────────────────────────────────────────────────────────
# 4. Bivariate coregionalization GP 샘플링
# ──────────────────────────────────────────────────────────────────
precompute_chol <- function(len_scale, coords) {
  K <- matern15_cov(coords, len_scale)
  diag(K) <- diag(K) + 1e-6
  t(chol(K))
}

sample_bivariate_gp <- function(L_Z, L_X, rho, seed) {
  set.seed(seed)
  xi_Z <- rnorm(n);  xi_X <- rnorm(n)
  Z_raw    <- as.vector(L_Z %*% xi_Z)
  Lz_inv_Z <- forwardsolve(L_Z, Z_raw)
  X_raw    <- rho * as.vector(L_X %*% Lz_inv_Z) +
              sqrt(max(1 - rho^2, 0)) * as.vector(L_X %*% xi_X)
  list(Z = (Z_raw - mean(Z_raw)) / sd(Z_raw),
       X = (X_raw - mean(X_raw)) / sd(X_raw))
}

# ──────────────────────────────────────────────────────────────────
# 5. Cholesky 캐싱
# ──────────────────────────────────────────────────────────────────
chol_cache <- list()
get_chol <- function(len_scale) {
  key <- as.character(len_scale)
  if (is.null(chol_cache[[key]]))
    chol_cache[[key]] <<- precompute_chol(len_scale, coords)
  chol_cache[[key]]
}

# ──────────────────────────────────────────────────────────────────
# 6. 단일 시뮬레이션
# ──────────────────────────────────────────────────────────────────
run_one <- function(l_z, l_x, rho, seed_base,
                    return_surface = FALSE) {
  L_Z <- get_chol(l_z);  L_X <- get_chol(l_x)
  gp  <- sample_bivariate_gp(L_Z, L_X, rho, seed = seed_base)
  Z   <- gp$Z;  X <- gp$X

  y <- BETA_TRUE * X + GAMMA_TRUE * Z + rnorm(n, 0, SIGMA_EPS)

  # ── Moran's I (DGP 검증) ────────────────────────────────────────
  moran_X <- moran_i(X)
  moran_Z <- moran_i(Z)

  # ── Short GWR (Z 누락) ──────────────────────────────────────────
  dat <- data.frame(y = y, X = X, cx = coords[, 1], cy = coords[, 2])
  coordinates(dat) <- ~cx + cy

  tryCatch({
    bw       <- bw.gwr(y ~ X, data = dat,
                       kernel = "bisquare", adaptive = TRUE,
                       approach = "AICc")
    gwr_res  <- gwr.basic(y ~ X, data = dat,
                          bw = bw, kernel = "bisquare", adaptive = TRUE)
    beta_hat <- gwr_res$SDF$X

    # ── Long GWR (Z 포함) ───────────────────────────────────────
    dat_long <- data.frame(y = y, X = X, Z = Z,
                           cx = coords[, 1], cy = coords[, 2])
    coordinates(dat_long) <- ~cx + cy
    bw_long      <- bw.gwr(y ~ X + Z, data = dat_long,
                           kernel = "bisquare", adaptive = TRUE,
                           approach = "AICc")
    gwr_long_res <- gwr.basic(y ~ X + Z, data = dat_long,
                              bw = bw_long, kernel = "bisquare",
                              adaptive = TRUE)
    beta_long    <- gwr_long_res$SDF$X

    list(
      gwr_std   = sd(beta_hat),
      gwr_bias  = mean(beta_hat)  - BETA_TRUE,
      lgwr_std  = sd(beta_long),
      lgwr_bias = mean(beta_long) - BETA_TRUE,
      moran_X   = moran_X,
      moran_Z   = moran_Z,
      bw        = bw,
      beta_hat  = if (return_surface) beta_hat  else NULL,
      beta_long = if (return_surface) beta_long else NULL,
      Z_surf    = if (return_surface) Z          else NULL,
      X_surf    = if (return_surface) X          else NULL
    )
  }, error = function(e) NULL)
}

# ──────────────────────────────────────────────────────────────────
# 7. 팩토리얼 실험
# ──────────────────────────────────────────────────────────────────
param_grid   <- expand.grid(l_z = L_Z_VALUES, l_x = L_X_VALUES,
                            rho = RHO_VALUES)
summary_recs <- vector("list", nrow(param_grid))
surface_recs <- vector("list", nrow(param_grid))

for (i in seq_len(nrow(param_grid))) {
  l_z <- param_grid$l_z[i]
  l_x <- param_grid$l_x[i]
  rho <- param_grid$rho[i]
  cat(sprintf("[%2d/%d] l_Z=%.1f  l_X=%.1f  rho=%.1f\n",
              i, nrow(param_grid), l_z, l_x, rho))

  sims        <- vector("list", N_SIM)
  rep_surface <- NULL;  rep_Z <- NULL;  rep_X <- NULL

  for (k in seq_len(N_SIM)) {
    r <- run_one(l_z, l_x, rho,
                 seed_base      = k * 100,
                 return_surface = (k == 1))
    if (!is.null(r)) {
      sims[[k]] <- r
      if (k == 1) {
        rep_surface <- r$beta_hat
        rep_Z       <- r$Z_surf
        rep_X       <- r$X_surf
        cat(sprintf("  k=1: BW=%d (%.0f%%)  Moran_X=%.3f  Moran_Z=%.3f\n",
                    r$bw, 100*r$bw/n, r$moran_X, r$moran_Z))
      }
    }
  }

  sims <- Filter(Negate(is.null), sims)

  if (length(sims) > 0) {
    summary_recs[[i]] <- data.frame(
      l_z         = l_z, l_x = l_x, rho = rho,
      gwr_std     = mean(sapply(sims, `[[`, "gwr_std")),
      gwr_bias    = mean(sapply(sims, `[[`, "gwr_bias")),
      lgwr_std    = mean(sapply(sims, `[[`, "lgwr_std")),
      lgwr_bias   = mean(sapply(sims, `[[`, "lgwr_bias")),
      moran_X     = mean(sapply(sims, `[[`, "moran_X")),
      moran_Z     = mean(sapply(sims, `[[`, "moran_Z")),
      bw_mean     = mean(sapply(sims, `[[`, "bw")),
      bw_sd       = sd(sapply(sims,   `[[`, "bw")),
      n_ok        = length(sims)
    )
    if (!is.null(rep_surface)) {
      surface_recs[[i]] <- data.frame(
        x         = coords[, 1],
        y         = coords[, 2],
        beta_hat  = rep_surface,
        beta_long = sims[[1]]$beta_long,
        Z_surf    = rep_Z,
        X_surf    = rep_X,
        moran_X   = sims[[1]]$moran_X,   # 대표 k=1 값
        moran_Z   = sims[[1]]$moran_Z,
        l_z = l_z, l_x = l_x, rho = rho
      )
    }
  }
}

df_summary  <- do.call(rbind, Filter(Negate(is.null), summary_recs))
df_surfaces <- do.call(rbind, Filter(Negate(is.null), surface_recs)) %>%
  left_join(df_summary %>%
              select(l_z, l_x, rho, gwr_std, gwr_bias,
                     lgwr_std, lgwr_bias, moran_X, moran_Z, bw_mean),
            by = c("l_z", "l_x", "rho"),
            suffix = c("_rep", "_mean"))

write.csv(df_summary,  "output_ovb/sim_results.csv",  row.names = FALSE)
write.csv(df_surfaces, "output_ovb/sim_surfaces.csv", row.names = FALSE)
cat("\n결과 저장 완료\n")
print(df_summary %>% select(l_z, l_x, rho,
                             gwr_std, gwr_bias,
                             lgwr_std, lgwr_bias,
                             moran_X, moran_Z))

# ──────────────────────────────────────────────────────────────────
# 공통 설정
# ──────────────────────────────────────────────────────────────────
lz_labeller  <- as_labeller(
  setNames(paste0("l[Z]==", L_Z_VALUES), as.character(L_Z_VALUES)),
  label_parsed)
rho_labeller <- as_labeller(
  setNames(paste0("rho==", RHO_VALUES), as.character(RHO_VALUES)),
  label_parsed)

make_heatmap <- function(df, metric, title_expr,
                         diverging = TRUE, midpt = 0) {
  sub <- df %>%
    mutate(
      l_z_f = factor(l_z, levels = L_Z_VALUES,
                     labels = paste0("l[Z]==", L_Z_VALUES)),
      l_x_f = factor(l_x, levels = L_X_VALUES,
                     labels = paste0("l[X]==", L_X_VALUES)),
      rho_f = factor(rho, levels = RHO_VALUES,
                     labels = paste0("rho==", RHO_VALUES))
    )
  ggplot(sub, aes(x = rho_f, y = l_z_f, fill = .data[[metric]])) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = round(.data[[metric]], 3)), size = 3) +
    facet_wrap(~l_x_f, nrow = 1, labeller = label_parsed) +
    { if (diverging)
        scale_fill_gradient2(low = "#2166AC", mid = "white",
                             high = "#D6604D", midpoint = midpt,
                             name = "value")
      else
        scale_fill_viridis(option = "plasma", direction = -1,
                           name = "value") } +
    scale_x_discrete(labels = \(x) parse(text = x)) +
    scale_y_discrete(labels = \(x) parse(text = x)) +
    labs(title = title_expr,
         x = expression(paste("X-Z correlation (", rho, ")")),
         y = expression(paste("Z spatial scale (", l[Z], ")"))) +
    theme_minimal(base_size = 10) +
    theme(plot.title      = element_text(face = "bold", size = 10),
          strip.text      = element_text(face = "bold", size = 10),
          legend.position = "right")
}

# ──────────────────────────────────────────────────────────────────
# Fig 1: Short GWR β̂_X 표면 매트릭스
# ──────────────────────────────────────────────────────────────────
surf_lim <- max(abs(df_surfaces$beta_hat - 1)) * 1.05

plot_beta <- vector("list", length(L_X_VALUES))
for (idx in seq_along(L_X_VALUES)) {
  lx  <- L_X_VALUES[idx]
  sub <- df_surfaces %>%
    filter(l_x == lx) %>%
    mutate(l_z_f = factor(l_z, levels = L_Z_VALUES),
           rho_f = factor(rho, levels = RHO_VALUES))
  ann <- sub %>%
    group_by(l_z_f, rho_f) %>%
    summarise(gwr_std  = first(gwr_std),
              gwr_bias = first(gwr_bias),
              bw_mean  = first(bw_mean), .groups = "drop") %>%
    mutate(label = paste0("std=",  round(gwr_std,  3),
                          "\nbias=", round(gwr_bias, 3),
                          "\nbw=",   round(bw_mean)))
  plot_beta[[idx]] <- ggplot(sub, aes(x = x, y = y, fill = beta_hat)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D6604D",
                         midpoint = 1, limits = c(1-surf_lim, 1+surf_lim),
                         name = expression(hat(beta)[X](u[i]))) +
    geom_text(data = ann, aes(x = Inf, y = -Inf, label = label),
              inherit.aes = FALSE, hjust = 1.05, vjust = -0.2,
              size = 2.0, color = "grey25") +
    facet_grid(l_z_f ~ rho_f,
               labeller = labeller(l_z_f = lz_labeller,
                                   rho_f = rho_labeller)) +
    coord_equal() +
    labs(title = bquote(l[X] == .(lx)), x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(plot.title      = element_text(face = "bold", size = 11,
                                         hjust = 0.5),
          strip.text      = element_text(size = 9, face = "bold"),
          legend.position = if (idx == length(L_X_VALUES)) "right"
                            else "none",
          axis.text       = element_blank(),
          axis.ticks      = element_blank(),
          panel.grid      = element_blank(),
          panel.border    = element_rect(color = "grey80", fill = NA,
                                         linewidth = 0.5))
}
p_beta <- wrap_plots(plot_beta, nrow = 1, guides = "collect") +
  plot_annotation(
    title    = expression(bold("Short GWR")~~hat(beta)[X](u[i])~~
               "Surface  |  "~beta[true]==1~"(constant)"),
    subtitle = "셀 하단: std / bias / 평균BW",
    theme    = theme(plot.title    = element_text(face = "bold", size = 12),
                     plot.subtitle = element_text(size = 9, color = "grey40")))
ggsave("output_ovb/fig_beta_surfaces_all.png", p_beta,
       width = 22, height = 7, dpi = 150)
cat("저장: fig_beta_surfaces_all.png\n")

# ──────────────────────────────────────────────────────────────────
# Fig 2: Long GWR β̂_X 표면 매트릭스
# ──────────────────────────────────────────────────────────────────
long_lim <- max(abs(df_surfaces$beta_long - BETA_TRUE)) * 1.05

plot_long <- vector("list", length(L_X_VALUES))
for (idx in seq_along(L_X_VALUES)) {
  lx  <- L_X_VALUES[idx]
  sub <- df_surfaces %>%
    filter(l_x == lx) %>%
    mutate(l_z_f = factor(l_z, levels = L_Z_VALUES),
           rho_f = factor(rho, levels = RHO_VALUES))
  ann <- sub %>%
    group_by(l_z_f, rho_f) %>%
    summarise(lgwr_std  = first(lgwr_std),
              lgwr_bias = first(lgwr_bias), .groups = "drop") %>%
    mutate(label = paste0("std=",  round(lgwr_std,  3),
                          "\nbias=", round(lgwr_bias, 3)))
  plot_long[[idx]] <- ggplot(sub, aes(x = x, y = y, fill = beta_long)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D6604D",
                         midpoint = BETA_TRUE,
                         limits = c(BETA_TRUE-long_lim,
                                    BETA_TRUE+long_lim),
                         oob = scales::squish,
                         name = expression(hat(beta)[X]^{long}(u[i]))) +
    geom_text(data = ann, aes(x = Inf, y = -Inf, label = label),
              inherit.aes = FALSE, hjust = 1.05, vjust = -0.2,
              size = 2.0, color = "grey25") +
    facet_grid(l_z_f ~ rho_f,
               labeller = labeller(l_z_f = lz_labeller,
                                   rho_f = rho_labeller)) +
    coord_equal() +
    labs(title = bquote(l[X] == .(lx)), x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(plot.title      = element_text(face = "bold", size = 11,
                                         hjust = 0.5),
          strip.text      = element_text(size = 9, face = "bold"),
          legend.position = if (idx == length(L_X_VALUES)) "right"
                            else "none",
          axis.text       = element_blank(),
          axis.ticks      = element_blank(),
          panel.grid      = element_blank(),
          panel.border    = element_rect(color = "grey80", fill = NA,
                                         linewidth = 0.5))
}
p_long <- wrap_plots(plot_long, nrow = 1, guides = "collect") +
  plot_annotation(
    title    = expression(bold("Long GWR")~~hat(beta)[X]^{long}(u[i])~~
               "Surface  (y ~ X + Z)  |  "~beta[true]==1),
    subtitle = "Z 포함 후 β̂_X 이질성 감소 확인 / 셀 하단: std / bias",
    theme    = theme(plot.title    = element_text(face = "bold", size = 12),
                     plot.subtitle = element_text(size = 9, color = "grey40")))
ggsave("output_ovb/fig_long_gwr_surfaces.png", p_long,
       width = 22, height = 7, dpi = 150)
cat("저장: fig_long_gwr_surfaces.png\n")

# ──────────────────────────────────────────────────────────────────
# Fig 3: Heatmap — Short vs Long GWR
# ──────────────────────────────────────────────────────────────────
p1 <- make_heatmap(df_summary, "gwr_std",
  expression(bold("Short GWR")~~hat(beta)[X]~~
             "Spatial Std  (spurious heterogeneity)"),
  diverging = FALSE)
p2 <- make_heatmap(df_summary, "gwr_bias",
  expression(bold("Short GWR")~~hat(beta)[X]~~
             "Mean Bias  ("*beta[true]==1*")"),
  diverging = TRUE, midpt = 0)
p3 <- make_heatmap(df_summary, "lgwr_std",
  expression(bold("Long GWR")~~hat(beta)[X]~~
             "Spatial Std  (Z 포함 후 이질성)"),
  diverging = FALSE)
p4 <- make_heatmap(df_summary, "lgwr_bias",
  expression(bold("Long GWR")~~hat(beta)[X]~~
             "Mean Bias  ("*beta[true]==1*")"),
  diverging = TRUE, midpt = 0)

p_heatmap <- p1 / p2 / p3 / p4 +
  plot_annotation(
    title = "Short vs Long GWR  |  beta_true=1, gamma=1  |  BW 반복 재탐색",
    theme = theme(plot.title = element_text(face = "bold", size = 13)))
ggsave("output_ovb/fig_heatmap_4row.png", p_heatmap,
       width = 13, height = 16, dpi = 150)
cat("저장: fig_heatmap_4row.png\n")

# ──────────────────────────────────────────────────────────────────
# Fig 4: Z, X 표면 매트릭스 — 27개 조합 + Moran's I 캡션
# ──────────────────────────────────────────────────────────────────
xz_lim <- max(abs(c(df_surfaces$Z_surf, df_surfaces$X_surf)),
              na.rm = TRUE) * 1.05

make_xz_panel <- function(surf_var, cmap, leg_name,
                           moran_col, moran_label) {
  plot_list <- vector("list", length(L_X_VALUES))
  for (idx in seq_along(L_X_VALUES)) {
    lx  <- L_X_VALUES[idx]
    sub <- df_surfaces %>%
      filter(l_x == lx) %>%
      mutate(l_z_f = factor(l_z, levels = L_Z_VALUES),
             rho_f = factor(rho, levels = RHO_VALUES))
    # Moran's I 캡션: 각 셀마다 대표 k=1 값
    ann <- sub %>%
      group_by(l_z_f, rho_f) %>%
      summarise(moran_val = first(.data[[moran_col]]), .groups = "drop") %>%
      mutate(label = paste0(moran_label, "=", round(moran_val, 3)))
    p <- ggplot(sub, aes(x = x, y = y, fill = .data[[surf_var]])) +
      geom_raster(interpolate = TRUE) +
      scale_fill_gradient2(low = "#2166AC", mid = "white",
                           high = "#D6604D", midpoint = 0,
                           limits = c(-xz_lim, xz_lim),
                           oob = scales::squish,
                           name = leg_name) +
      geom_text(data = ann, aes(x = Inf, y = -Inf, label = label),
                inherit.aes = FALSE, hjust = 1.05, vjust = -0.3,
                size = 2.2, color = "grey20", fontface = "bold") +
      facet_grid(l_z_f ~ rho_f,
                 labeller = labeller(l_z_f = lz_labeller,
                                     rho_f = rho_labeller)) +
      coord_equal() +
      labs(title = bquote(l[X] == .(lx)), x = NULL, y = NULL) +
      theme_minimal(base_size = 9) +
      theme(plot.title      = element_text(face = "bold", size = 11,
                                           hjust = 0.5),
            strip.text      = element_text(size = 9, face = "bold"),
            legend.position = if (idx == length(L_X_VALUES)) "right"
                              else "none",
            axis.text       = element_blank(),
            axis.ticks      = element_blank(),
            panel.grid      = element_blank(),
            panel.border    = element_rect(color = "grey80", fill = NA,
                                           linewidth = 0.5))
    plot_list[[idx]] <- p
  }
  plot_list
}

pl_Z <- make_xz_panel("Z_surf", "BrBG", "Z",
                       "moran_Z_rep", "I(Z)")
pl_X <- make_xz_panel("X_surf", "RdBu_r", "X",
                       "moran_X_rep", "I(X)")

p_Z_mat <- wrap_plots(pl_Z, nrow = 1, guides = "collect") +
  plot_annotation(
    title = "Z 표면  |  Z ~ GP(K_Z)  |  셀 우하단: Moran's I(Z)",
    theme = theme(plot.title = element_text(face = "bold", size = 12)))
p_X_mat <- wrap_plots(pl_X, nrow = 1, guides = "collect") +
  plot_annotation(
    title = "X 표면  |  X|Z ~ GP(...)  |  셀 우하단: Moran's I(X)",
    theme = theme(plot.title = element_text(face = "bold", size = 12)))

p_xz_all <- p_Z_mat / p_X_mat +
  plot_annotation(
    title    = "Z, X 표면 매트릭스  |  27개 조합  (대표 시뮬레이션, k=1)",
    subtitle = paste0("DGP 검증: l_Z↑ → I(Z)↑,  l_X↑ → I(X)↑,",
                      "  rho: Z·X 공간 구조 공유 강도"),
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "grey40")))
ggsave("output_ovb/fig_XZ_surfaces.png", p_xz_all,
       width = 22, height = 14, dpi = 150)
cat("저장: fig_XZ_surfaces.png\n")

# ──────────────────────────────────────────────────────────────────
# Fig 5: Moran's I Heatmap (DGP 검증 요약)
# ──────────────────────────────────────────────────────────────────
p_mI_X <- make_heatmap(df_summary, "moran_X",
  expression("Moran's"~~italic(I)(X)~~"(Monte Carlo 평균)"),
  diverging = TRUE, midpt = 0)
p_mI_Z <- make_heatmap(df_summary, "moran_Z",
  expression("Moran's"~~italic(I)(Z)~~"(Monte Carlo 평균)"),
  diverging = TRUE, midpt = 0)

p_moran <- p_mI_Z / p_mI_X +
  plot_annotation(
    title    = "DGP 검증: Moran's I  |  l_Z/l_X에 따른 공간 스케일 확인",
    subtitle = "I(Z)는 l_Z에만 의존, I(X)는 l_X에 의존  (rho 무관)",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "grey40")))
ggsave("output_ovb/fig_moran_heatmap.png", p_moran,
       width = 13, height = 8, dpi = 150)
cat("저장: fig_moran_heatmap.png\n")