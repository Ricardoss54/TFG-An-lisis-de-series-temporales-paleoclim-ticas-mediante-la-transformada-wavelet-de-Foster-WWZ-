#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Weighted Wavelet Z-Transform (WWZ) -- funciones (v3)
# Foster, G. (1996). Wavelets for period analysis of unevenly
# sampled time series. The Astronomical Journal, 112(4), 1709-1729.
#
# (C) 2026. Ricardo Sanchez Santos & Josue M. Polanco-Martinez
# GNU General Public License v3 -- https://www.gnu.org/licenses/
#
# Este fichero contiene SOLO funciones (sin experimentos), para poder
# cargarlas con source() desde los experimentos y desde los tests.
#
# Cambios respecto a WWZfinal_v2.R:
#  1. Las celdas no calculables devuelven NA (antes 0). Un 0 significa
#     "calculado, sin senal"; un NA significa "no se pudo calcular".
#  2. neff < neff_min -> NA (antes se forzaba a 3 + 1e-6).
#  3. n_sim = 200 por defecto, argumento `seed` y simulaciones generadas
#     en el proceso principal: resultado reproducible e independiente
#     del numero de nucleos.
#  4. wwz_core() es ahora una funcion de primer nivel (testeable) y las
#     tolerancias absolutas (1e-12) pasan a ser relativas a la escala de x.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggplot2)
library(gridExtra)
library(viridis)
library(parallel)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HELPER: guardar figuras en figures/
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
save_plot <- function(filename, plot, width = 10, height = 6, dpi = 300) {
  fig_dir <- "figures"
  dir.create(fig_dir, showWarnings = FALSE)
  path <- file.path(fig_dir, filename)
  ggsave(path, plot = plot, width = width, height = height, dpi = dpi)
  message(paste("Guardado en:", path))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NUCLEO: estadistico Z de la WWZ en una celda (tau0, w)
#
# Valor devuelto:
#   numero >= 0 : estadistico Z calculado (0 = sin varianza explicada)
#   NA          : celda no calculable (sin datos cerca de tau0,
#                 neff insuficiente, sistema mal condicionado o
#                 ventana con x constante)
#
# Argumentos:
#   c_decay   constante de decaimiento de la ventana gaussiana
#   sw_min    suma minima de pesos. 1e-12 reproduce la v2. Con 1 se
#             exige al menos "un dato completo" de peso en la ventana;
#             evita celdas dentro de huecos grandes cuyo ajuste solo
#             usa los datos de los extremos del hueco. Conviene explorarlo.
#   neff_min  neff minimo. La v2 equivale a 3 + 1e-6. Con 3 parametros
#             ajustados, neff - 3 son los grados de libertad del
#             denominador de la F(2, neff - 3).
#   tol_const ventana "constante": sd ponderada <= tol_const * max|x|
#   tol_fit   ajuste "perfecto": Vr se acota inferiormente por
#             tol_fit * Vx, de modo que Z es enorme pero finito (la v2
#             devolvia 0 justo donde el ajuste era mas perfecto).
#   full      TRUE devuelve c(Z =, neff =) (para tests de calibracion)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wwz_core <- function(t, x, w, tau0,
                     c_decay   = 0.0125,
                     neff_min  = 4,
                     sw_min    = 1e-12,
                     tol_const = 1e-10,
                     tol_fit   = 1e-8,
                     full      = FALSE) {
  out <- function(z, neff = NA_real_) if (full) c(Z = z, neff = neff) else z

  wt <- exp(-c_decay * w^2 * (t - tau0)^2)
  Sw <- sum(wt)
  if (!is.finite(Sw) || Sw < sw_min) return(out(NA_real_))

  neff <- Sw^2 / sum(wt^2)
  if (!is.finite(neff) || neff < neff_min) return(out(NA_real_, neff))

  X   <- cbind(1, cos(w * t), sin(w * t))
  XtW <- t(X * wt)
  A   <- XtW %*% X
  if (rcond(A) < 1e-12) return(out(NA_real_, neff))

  beta <- solve(A, XtW %*% x)
  yhat <- drop(X %*% beta)

  xbar <- sum(wt * x) / Sw
  Vx   <- sum(wt * (x - xbar)^2)    # varianza ponderada total
  Vr   <- sum(wt * (x - yhat)^2)    # varianza residual
  if (sqrt(Vx / Sw) <= tol_const * max(abs(x))) return(out(NA_real_, neff))

  # Vx - Vr = varianza explicada por el modelo (>= 0 salvo redondeo)
  z <- ((neff - 3) / 2) * max(Vx - Vr, 0) / max(Vr, tol_fit * Vx)
  out(z, neff)
}

# Matriz WWZ [nf x nt] de una serie (x, time) sobre la malla (freq, tau)
wwz_matrix <- function(x, time, freq, tau,
                       c_decay = 0.0125, neff_min = 4, sw_min = 1e-12) {
  M <- matrix(NA_real_, length(freq), length(tau))
  for (i in seq_along(freq)) {
    w <- 2 * pi * freq[i]
    for (j in seq_along(tau))
      M[i, j] <- wwz_core(time, x, w, tau[j],
                          c_decay = c_decay, neff_min = neff_min,
                          sw_min = sw_min)
  }
  M
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCION PRINCIPAL
# Los argumentos nuevos (seed, neff_min, sw_min) van al final para no
# romper llamadas posicionales antiguas.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wwz_robust <- function(time, x,
                       freq = seq(1e-04, 0.001, length.out = 200),
                       tau_n = 120,
                       c = 0.0125,
                       sig_level = 0.95,
                       n_sim = 200,
                       use_mc = TRUE,
                       return_plot = TRUE,
                       show_plot = TRUE,
                       seed = NULL,
                       neff_min = 4,
                       sw_min = 1e-12) {

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 0. CLEAN DATA
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ok <- is.finite(time) & is.finite(x)
  time <- time[ok]; x <- x[ok]

  ord <- order(time)
  time <- time[ord]; x <- x[ord]

  n   <- length(x)
  tau <- seq(min(time), max(time), length.out = tau_n)
  nf  <- length(freq)
  nt  <- length(tau)

  if (use_mc && n_sim < 100)
    warning("n_sim = ", n_sim, ": con pocas simulaciones el cuantil ",
            sig_level, " por celda es muy ruidoso (recomendado >= 200).")

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 1. MATRIZ WWZ (NA = celda no calculable)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  WWZ <- wwz_matrix(x, time, freq, tau,
                    c_decay = c, neff_min = neff_min, sw_min = sw_min)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2. AR(1) -- blindado contra series deterministas (sigma2 ~ 0)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  fit <- tryCatch(
    arima(x, order = c(1, 0, 0), method = "CSS-ML"),
    error = function(e) NULL
  )
  if (is.null(fit) || !is.finite(fit$sigma2) || fit$sigma2 < 1e-10) {
    phi   <- 0
    sigma <- max(sd(x) * 0.01, 1e-6)
  } else {
    phi   <- as.numeric(fit$coef["ar1"])
    sigma <- sqrt(fit$sigma2)
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 3. SIGNIFICANCIA (MC AR1)
  # Las series simuladas se generan AQUI, en el proceso principal, con
  # la semilla dada: el resultado no depende del numero de nucleos.
  # Los workers solo calculan la matriz WWZ de cada serie.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  WWZ_sig <- matrix(NA_real_, nf, nt)

  if (use_mc) {
    if (!is.null(seed)) set.seed(seed)
    sims <- lapply(seq_len(n_sim), function(k)
      as.numeric(arima.sim(n = n, list(ar = phi), sd = sigma)))

    n_cores <- max(1L, min(detectCores() - 1L, n_sim))
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl), add = TRUE)   # se cierra aunque haya error
    clusterExport(cl, c("wwz_core", "wwz_matrix"),
                  envir = environment(wwz_core))

    sim_list <- parLapply(cl, sims, wwz_matrix,
                          time = time, freq = freq, tau = tau,
                          c_decay = c, neff_min = neff_min, sw_min = sw_min)

    # Cuantil por celda; NA si menos de la mitad de las sims es valida
    for (i in seq_len(nf)) {
      for (j in seq_len(nt)) {
        v <- vapply(sim_list, `[`, 0, i, j)
        v <- v[is.finite(v)]
        if (length(v) >= 0.5 * n_sim)
          WWZ_sig[i, j] <- quantile(v, sig_level, names = FALSE)
      }
    }
  } else {
    # Umbral empirico simple (como en v2). Ojo: el TFG (ec. 3.1) usa
    # media + k*sd con k = 1; aqui es solo la media. A decidir.
    WWZ_sig[] <- mean(WWZ, na.rm = TRUE)
  }

  # Una celda es significativa solo si ambas cantidades existen
  signif_mask <- !is.na(WWZ) & !is.na(WWZ_sig) & (WWZ > WWZ_sig)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 4. COI
  # OJO: aqui edge = periodo / sqrt(2); el apendice del TFG tiene
  # periodo * sqrt(2). Hay que decidir cual es y, idealmente, derivarlo
  # de c (el e-folding de la ventana es 1 / (w * sqrt(c))).
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  period <- 1 / freq
  coi <- matrix(0, nf, nt)
  for (i in 1:nf) {
    edge <- 1 / (2 * pi * freq[i] * sqrt(c))
    for (j in 1:nt) {
      dist_edge <- min(tau[j] - min(time), max(time) - tau[j])
      coi[i, j] <- min(dist_edge / edge, 1)
    }
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 5. DATAFRAME
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  df <- data.frame(
    tau  = rep(tau, each = nf),
    freq = rep(freq, times = nt),
    WWZ  = as.vector(WWZ),
    sig  = as.vector(signif_mask),
    coi  = as.vector(coi)
  )
  df$WWZ_plot <- df$WWZ
  df$WWZ_plot[!df$sig] <- df$WWZ_plot[!df$sig] * 0.5   # NA sigue siendo NA
  df$coi[is.na(df$coi)] <- 0

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 6. PLOT (celdas NA en gris claro)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  plt_wwz <- ggplot(df, aes(x = tau, y = freq)) +
    geom_tile(aes(fill = WWZ_plot, alpha = coi)) +
    scale_fill_viridis(option = "C", na.value = "grey85") +
    scale_alpha(range = c(0.2, 1), guide = "none") +
    labs(title = "WWZ robust (AR1 + COI + significance)",
         x = "Tiempo", y = "Frecuencia") +
    theme_minimal(base_size = 12)

  if (show_plot) print(plt_wwz)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 7. RETURN
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  return(list(
    WWZ             = WWZ,
    WWZ_significant = signif_mask,
    COI             = coi,
    AR1_phi         = phi,
    AR1_sigma       = sigma,
    plot            = plt_wwz,
    data            = df,
    params = list(c = c, n = n, n_sim = if (use_mc) n_sim else 0L,
                  seed = seed, neff_min = neff_min, sw_min = sw_min,
                  n_cells = length(WWZ), n_cells_NA = sum(is.na(WWZ)))
  ))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT AUXILIAR (celdas NA en gris claro)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
make_wwz_plot <- function(df, titulo, zlim,
                          contour_col = "white",
                          contour_lwd = 0.7) {
  df$sig_num <- as.numeric(df$sig)
  ggplot(df, aes(x = tau, y = freq)) +
    geom_tile(aes(fill = WWZ_plot, alpha = coi)) +
    scale_fill_viridis(option = "C", limits = zlim, name = "WWZ",
                       na.value = "grey85") +
    scale_alpha(range = c(0.2, 1), guide = "none") +
    geom_contour(aes(z = sig_num), breaks = 0.5,
                 colour = contour_col, linewidth = contour_lwd) +
    labs(title = titulo, x = "Tiempo", y = "Frecuencia") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AUXILIAR: eliminar un % de datos al azar
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
create_missing <- function(x, t, porcentaje, seed = 123) {
  set.seed(seed)
  N         <- length(x)
  n_missing <- round(porcentaje / 100 * N)
  idx_missing <- sort(sample(N, n_missing))
  idx_present <- setdiff(seq_len(N), idx_missing)
  x_na              <- x
  x_na[idx_missing] <- NA
  cat(sprintf("  [%3d%%]  %d datos eliminados | %d conservados\n",
              porcentaje, n_missing, length(idx_present)))
  list(
    porcentaje  = porcentaje,
    n_missing   = n_missing,
    idx_missing = idx_missing,
    idx_present = idx_present,
    t_obs       = t[idx_present],
    x_obs       = x[idx_present],
    x_na        = x_na
  )
}
