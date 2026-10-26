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
#
# Cambios v4 (esta version):
#  5. wwz_matrix estandariza x una vez (Z es invariante a a + s*x).
#  6. Monte Carlo por lotes (WWZ_batch.R): pesos y A^-1 una vez por celda,
#     K simulaciones por producto de matrices; reparto entre nucleos POR
#     FRECUENCIAS; resultado identico a la v3 con la misma semilla.
#  7. noise = "ou": ruido rojo en tiempo continuo ajustado sobre los
#     residuos tras prewhitening (wwz_prewhiten, wwz_fit_ou, wwz_sim_ou).
#     noise = "ar1" (defecto) reproduce la v3.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggplot2)
library(gridExtra)
library(viridis)
library(parallel)
source("WWZ_batch.R")   # wwz_matrix_batch(): MC por lotes

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
#             (x se espera ya estandarizada, como hace wwz_matrix;
#             si se llama a wwz_core directamente con x cruda, un offset
#             grande frente a su variacion pierde precision numerica)
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
#
# x se estandariza UNA vez aqui: x' = (x - media) / sd. Z es exactamente
# invariante a x -> a + s*x (el modelo incluye termino constante, y Vx y
# Vr escalan por s^2), pero en doble precision un offset grande frente a
# la variacion de x degrada el ajuste. Asi la serie real y las simuladas
# (wwz_matrix_batch estandariza igual) se calculan de la misma manera.
# Serie constante (sd = 0 o no finita) -> toda la matriz NA.
wwz_matrix <- function(x, time, freq, tau,
                       c_decay = 0.0125, neff_min = 4, sw_min = 1e-12) {
  M <- matrix(NA_real_, length(freq), length(tau))
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(M)
  x <- (x - mean(x)) / s
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
# RUIDO ROJO EN TIEMPO CONTINUO (OU) SOBRE RESIDUOS
#
# El AR(1) por indice ajustado a la serie ENTERA tiene dos problemas:
#  (a) la senal entra en el ajuste: con un seno limpio y muestreo denso
#      phi ~ 1 y el nulo se traga la senal;
#  (b) en series irregulares phi mezcla pasos de distinto tamano.
# Aqui: (1) se blanquean (prewhitening) las componentes sinusoidales
# significativas, (2) se ajusta sobre los residuos un proceso OU
#   x_i = a_i x_{i-1} + sqrt(1 - a_i^2) e_i,  a_i = exp(-dt_i / tau)
# por maxima verosimilitud, y (3) se simula con los tiempos reales.
# Z es invariante a la escala de x, asi que solo importa tau.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SS explicada por una sinusoide (+ const.) a la frecuencia f, para cada
# columna de R (n x K, ya centrada)
.ls_explained <- function(time, R, f) {
  n  <- nrow(R)
  Fm <- qr.fitted(qr(cbind(1, cos(2*pi*f*time), sin(2*pi*f*time))), R)
  colSums(Fm^2) - colSums(Fm)^2 / n
}

# Prewhitening SECUENCIAL contra el nulo de ruido rojo: en cada paso se
# estima tau sobre los residuos actuales, se simulan K series OU con los
# tiempos reales y se compara la mayor fraccion de varianza explicada por
# una sinusoide con la distribucion del maximo bajo ese nulo (corrige de
# paso las multiples frecuencias). Si no es significativa, se para.
# (Un test F de ruido blanco quitaba componentes espurias a todo ruido rojo.)
wwz_prewhiten <- function(time, x, max_comp = 5, alpha = 0.01,
                          K = 199, oversample = 2) {
  n <- length(x); Tspan <- diff(range(time))
  fmax <- 0.5 / median(diff(time))
  fg <- seq(1 / (oversample * Tspan), fmax,
            length.out = max(20, ceiling(oversample * Tspan * fmax)))
  freqs <- numeric(0); pvals <- numeric(0)
  res <- x - mean(x); tau <- wwz_fit_ou(time, res)
  design <- function(fr) cbind(1, do.call(cbind, lapply(fr, function(f)
    cbind(cos(2*pi*f*time), sin(2*pi*f*time)))))
  for (k in seq_len(max_comp)) {
    if (n - 3 - 2 * length(freqs) <= 5) break
    obs <- vapply(fg, function(f) .ls_explained(time, matrix(res, n, 1), f), 0)
    ib <- which.max(obs); step <- fg[2] - fg[1]
    f_ref <- optimize(function(f) .ls_explained(time, matrix(res, n, 1), f),
                      c(max(fg[ib] - step, fg[1] / 2), fg[ib] + step),
                      maximum = TRUE)$maximum
    frac_obs <- max(obs, .ls_explained(time, matrix(res, n, 1), f_ref)) / sum(res^2)
    # tau con el candidato YA retirado: si no, la propia senal (p. ej. un
    # seno limpio) da tau enorme y un nulo tipo paseo aleatorio que se la traga
    res_c <- x - qr.fitted(qr(design(c(freqs, f_ref))), x)
    tau_c <- wwz_fit_ou(time, res_c)
    S <- wwz_sim_ou(time, tau_c, K); S <- sweep(S, 2, colMeans(S))
    tot <- colSums(S^2)
    mx  <- rep(0, K)
    for (f in fg) mx <- pmax(mx, .ls_explained(time, S, f))
    p <- (1 + sum(mx / tot >= frac_obs)) / (K + 1)
    if (p > alpha) break
    freqs <- c(freqs, f_ref); pvals <- c(pvals, p)
    res <- res_c; tau <- tau_c
  }
  list(resid = res, freq = freqs, p = pvals, tau = tau)
}

# tau (maxima verosimilitud) de un OU con muestreo irregular
wwz_fit_ou <- function(time, r) {
  r  <- (r - mean(r)) / sd(r)
  dt <- pmax(diff(time), 1e-8 * diff(range(time)))
  nll <- function(lt) {
    a  <- exp(-dt / exp(lt))
    v  <- 1 - a^2
    -sum(dnorm(r[-1], a * r[-length(r)], sqrt(v), log = TRUE))
  }
  lo <- log(min(dt) / 10); hi <- log(10 * diff(range(time)))
  o  <- optimize(nll, c(lo, hi))
  exp(o$minimum)
}

# K series OU (matriz n x K) con los tiempos reales
wwz_sim_ou <- function(time, tau, K) {
  n  <- length(time)
  a  <- exp(-pmax(diff(time), 0) / tau)
  sd_i <- sqrt(1 - a^2)
  Z <- matrix(rnorm(n * K), n, K)
  X <- Z
  for (i in 2:n) X[i, ] <- a[i - 1] * X[i - 1, ] + sd_i[i - 1] * Z[i, ]
  X
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
                       sw_min = 1e-12,
                       n_cores = NULL,
                       noise = c("ar1", "ou"),
                       prewhiten_max = 5) {
  noise <- match.arg(noise)

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
  pw <- NULL; tau_hat <- NA_real_
  if (use_mc && noise == "ou") {
    if (!is.null(seed)) set.seed(seed)
    pw <- wwz_prewhiten(time, x, max_comp = prewhiten_max)
    tau_hat <- pw$tau
    phi <- NA_real_; sigma <- NA_real_
  } else {
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
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 3. SIGNIFICANCIA (MC AR1)
  # Las series simuladas se generan AQUI, en el proceso principal, con
  # la semilla dada: el resultado no depende del numero de nucleos.
  # Los workers solo calculan la matriz WWZ de cada serie.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  WWZ_sig <- matrix(NA_real_, nf, nt)

  if (use_mc) {
    if (!is.null(seed) && noise == "ar1") set.seed(seed)
    # Matriz n x n_sim (mismo orden de RNG que la v3: misma semilla,
    # mismas series)
    Xsim <- if (noise == "ou") wwz_sim_ou(time, tau_hat, n_sim) else
      vapply(seq_len(n_sim), function(k)
        as.numeric(arima.sim(n = n, list(ar = phi), sd = sigma)), numeric(n))

    # Reparto POR FRECUENCIAS (intercaladas, para equilibrar la carga):
    # los pesos y A^-1 de cada celda solo dependen de los tiempos, asi que
    # cada celda se prepara una sola vez y se aplica a las n_sim series.
    # Cada worker devuelve ya los cuantiles [nf_bloque x nt].
    min_valid <- 0.5 * n_sim
    block_fun <- function(idx, Xsim, time, freq, tau, c_decay, neff_min,
                          sw_min, sig_level, min_valid) {
      B  <- wwz_matrix_batch(Xsim, time, freq[idx], tau, c_decay = c_decay,
                             neff_min = neff_min, sw_min = sw_min)
      B2 <- matrix(B, length(idx) * length(tau), dim(B)[3])
      q  <- apply(B2, 1, function(v) {
        v <- v[is.finite(v)]
        if (length(v) >= min_valid) quantile(v, sig_level, names = FALSE)
        else NA_real_
      })
      matrix(q, length(idx), length(tau))
    }

    if (is.null(n_cores)) n_cores <- detectCores() - 1L
    n_cores <- max(1L, min(as.integer(n_cores), nf))
    blocks  <- lapply(seq_len(n_cores), function(b) seq(b, nf, by = n_cores))
    args <- list(Xsim = Xsim, time = time, freq = freq, tau = tau,
                 c_decay = c, neff_min = neff_min, sw_min = sw_min,
                 sig_level = sig_level, min_valid = min_valid)

    if (n_cores == 1L) {
      res <- lapply(blocks, function(idx) do.call(block_fun, c(list(idx), args)))
    } else {
      cl <- makeCluster(n_cores)
      on.exit(stopCluster(cl), add = TRUE)   # se cierra aunque haya error
      clusterExport(cl, c("wwz_cells_batch", "wwz_matrix_batch"),
                    envir = environment(wwz_matrix_batch))
      res <- do.call(parLapply, c(list(cl, blocks, block_fun), args))
    }
    for (b in seq_len(n_cores)) WWZ_sig[blocks[[b]], ] <- res[[b]]
  } else {
    # Umbral empirico simple (como en v2). Ojo: el TFG (ec. 3.1) usa
    # media + k*sd con k = 1; aqui es solo la media. A decidir.
    WWZ_sig[] <- mean(WWZ, na.rm = TRUE)
  }

  # Una celda es significativa solo si ambas cantidades existen
  signif_mask <- !is.na(WWZ) & !is.na(WWZ_sig) & (WWZ > WWZ_sig)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 4. COI
  # Distancia de e-folding de la ventana gaussiana
  #   exp(-c w^2 dt^2) = e^-1  =>  dt = 1 / (w sqrt(c)) = 1 / (2 pi f sqrt(c)).
  # coi = min(distancia al borde / dt, 1): 1 = lejos del borde.
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
    noise           = list(model = noise, tau = tau_hat,
                           removed_freq = if (is.null(pw)) NULL else pw$freq),
    plot            = plt_wwz,
    data            = df,
    params = list(c = c, n = n, n_sim = if (use_mc) n_sim else 0L,
                  seed = seed, n_cores_used = if (use_mc) n_cores else 0L, neff_min = neff_min, sw_min = sw_min,
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
