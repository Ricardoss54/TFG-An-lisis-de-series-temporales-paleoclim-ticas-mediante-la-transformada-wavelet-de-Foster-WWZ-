
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Weighted Wavelet Z-Transform (WWZ)
# Implementation in R
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description:
#   This script implements the Weighted Wavelet Z-Transform (WWZ),
#   a time-frequency analysis method specifically designed for
#   unevenly spaced time series. The method was originally proposed
#   by Foster (1996) and is particularly useful in paleoclimatology,
#   astronomy, and any field dealing with irregularly sampled data.
#
# Key reference:
#   Foster, G. (1996). Wavelets for period analysis of unevenly
#   sampled time series. The Astronomical Journal, 112(4), 1709-1729.
#   https://doi.org/10.1086/118137

# Other important references: (@Ricardo: estudialas)
#   Witt, A. & Schumann, A. Y. (2005). Holocene climate variability
#     on millennial scales recorded in Greenland ice cores.
#     Nonlinear Processes in Geophysics, 12, 345-352.
# 
#   Grinsted, A., Moore, J. C. & Jevrejeva, S. (2004). Application
#     of the cross wavelet transform and wavelet coherence to
#     geophysical time series. Nonlinear Processes in Geophysics,
#     11, 561-566. https://doi.org/10.5194/npg-11-561-2004
# 
#   Schulz, M., & Mudelsee, M. (2002). REDFIT: estimating red-noise 
#     spectra directly from unevenly spaced paleoclimatic time series. 
#     Computers & Geosciences, 28(3), 421-426. Este paper es un must! 
# 
#   Torrence, C., & Compo, G. P. (1998). A practical guide to wavelet 
#     analysis. Bulletin of the American Meteorological society, 79(1), 
# 61-78. Este paper es un must! 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (C) 2026. Ricardo Sánchez Santos & Josue M. Polanco-Martinez
# ricardosanchezsantos@usal.es | 
# josue.polanco@usal.es | josue.m.polanco@gmail.com
# Dept. of Statistics, Fac. of Sciences, Univ. of Salamanca, Spain
# March 2026, Salamanca, Spain
# GNU General Public License v3 — https://www.gnu.org/licenses/
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCIÓN ENVOLTORIO: PREPARADOR INTELIGENTE PARA WWZ
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggplot2)
library(gridExtra)
library(viridis)
library(parallel)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PARÁMETROS GLOBALES — comunes a todas las pruebas
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
N  <- 1000
p1 <- 23   # periodo 1 (único periodo en pruebas de un seno)
p2 <- 41   # periodo 2 (sólo en pruebas de dos senos)
if (p1 > p2) { temp <- p1; p1 <- p2; p2 <- temp }

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

wwz_robust <- function(time, x,
                       freq = seq(1e-04, 0.001, length.out = 200),
                       #freq = seq(0.01, 0.5, length.out = 80),
                       tau_n = 120,
                       c = 0.0125,
                       sig_level = 0.95,
                       n_sim = 50,
                       use_mc = TRUE,
                       return_plot = TRUE,
                       show_plot = TRUE){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 0. CLEAN DATA
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ok <- is.finite(time) & is.finite(x)
  time <- time[ok]; x <- x[ok]
  
  ord <- order(time)
  time <- time[ord]; x <- x[ord]
  
  n <- length(x)
  tau <- seq(min(time), max(time), length.out = tau_n)
  
  nf <- length(freq)
  nt <- length(tau)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 1. WWZ CORE
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  wwz_core <- function(t, x, w, tau0){
    
    wt <- exp(-c * w^2 * (t - tau0)^2)
    
    Sw <- sum(wt)
    if(Sw < 1e-12) return(0)
    
    X <- cbind(1, cos(w*t), sin(w*t))
    
    XtW <- t(X * wt)
    A <- XtW %*% X
    
    if(rcond(A) < 1e-12) return(0)
    
    beta <- solve(A, XtW %*% x)
    yhat <- X %*% beta
    
    xbar <- sum(wt * x) / Sw
    Vx <- sum(wt * (x - xbar)^2)
    Vy <- sum(wt * (x - yhat)^2)
    
    if(Vy < 1e-12) return(0)
    
    neff <- max(sum(wt)^2 / sum(wt^2), 3 + 1e-6)
    if((Vx - Vy) < 1e-12) return(0)
    ((neff - 3) / 2) * (Vx - Vy) / Vy
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2. WWZ MATRIX (SAFE)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  WWZ <- matrix(0, nf, nt)
  
  for(i in 1:nf){
    w <- 2*pi*freq[i]
    for(j in 1:nt){
      WWZ[i, j] <- wwz_core(time, x, w, tau[j])
    }
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 3. AR(1) — blindado contra series deterministas (sigma2 ~ 0)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  fit <- tryCatch(
    arima(x, order = c(1, 0, 0), method = "CSS-ML"),
    error = function(e) NULL
  )
  if (is.null(fit) || !is.finite(fit$sigma2) || fit$sigma2 < 1e-10) {
    # Serie casi determinista: referencia de ruido blanco mínimo
    phi   <- 0
    sigma <- max(sd(x) * 0.01, 1e-6)
  } else {
    phi   <- as.numeric(fit$coef["ar1"])
    sigma <- sqrt(fit$sigma2)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 4. SIGNIFICANCE (MC AR1) — paralelizado
  #   Cada worker calcula UNA simulación completa [nf × nt].
  #   Se evita guardar el array gordo [nf × nt × n_sim]:
  #   el quantile se calcula celda a celda al final.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  WWZ_sig <- matrix(0, nf, nt)
  
  if (use_mc) {
    n_cores <- max(1L, detectCores() - 1L)   # deja 1 core libre al SO
    cl      <- makeCluster(n_cores)
    
    
    clusterExport(cl,
                  varlist = c("wwz_core", "c", "freq", "tau", "time",
                              "nf", "nt", "n", "phi", "sigma"),
                  envir   = environment())
    
    # Cada elemento de sim_list es una matriz [nf × nt] de una simulación AR1
    sim_list <- parLapply(cl, seq_len(n_sim), function(k) {
      sim <- as.numeric(arima.sim(n = n, list(ar = phi), sd = sigma))
      mat <- matrix(0, nf, nt)
      for (i in seq_len(nf)) {
        w <- 2 * pi * freq[i]
        for (j in seq_len(nt))
          mat[i, j] <- wwz_core(time, sim, w, tau[j])
      }
      mat
    })
    
    stopCluster(cl)
    
    # Quantile 95% celda a celda (sin el array [nf × nt × n_sim] en memoria)
    for (i in seq_len(nf))
      for (j in seq_len(nt))
        WWZ_sig[i, j] <- quantile(vapply(sim_list, `[`, 0, i, j),
                                  sig_level, na.rm = TRUE)
    
  } else {
    WWZ_sig <- matrix(mean(WWZ), nf, nt)
  }
  
  signif_mask <- WWZ > WWZ_sig
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 5. COI (SMOOTH)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  period <- 1 / freq
  coi <- matrix(0, nf, nt)
  
  for(i in 1:nf){
    edge <- period[i] / sqrt(2)
    
    for(j in 1:nt){
      dist_edge <- min(tau[j] - min(time), max(time) - tau[j])
      coi[i, j] <- min(dist_edge / edge, 1)
    }
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 6. DATAFRAME SAFE (NO expand.grid BUGS)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  df <- data.frame(
    tau = rep(tau, each = nf),
    freq = rep(freq, times = nt),
    WWZ = as.vector(WWZ),
    sig = as.vector(signif_mask),
    coi = as.vector(coi)
  )
  
  df$WWZ_plot <- df$WWZ
  df$WWZ_plot[!df$sig] <- df$WWZ_plot[!df$sig] * 0.5
  
  # safety
  df$coi[is.na(df$coi)] <- 0
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 7. PLOT (ROBUST + ALWAYS WORKS)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  plt_wwz <- ggplot(df, aes(x = tau, y = freq)) +
    geom_tile(aes(fill = WWZ_plot, alpha = coi)) +
    scale_fill_viridis(option = "C") +
    scale_alpha(range = c(0.2, 1), guide = "none") +
    labs(
      title = "WWZ robust (AR1 + COI + significance)",
      x = "Tiempo",
      y = "Frecuencia"
    ) +
    theme_minimal(base_size = 12)
  
  if(show_plot) print(plt_wwz)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 8. RETURN
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  return(list(
    WWZ = WWZ,
    WWZ_significant = signif_mask,
    COI = coi,
    AR1_phi = phi,
    AR1_sigma = sigma,
    plot = plt_wwz,
    data = df
  ))
}

# FUNCIÓN AUXILIAR DE PLOT WWZ — definida una sola vez, usada en todas las pruebas
# contour_col / contour_lwd: controlan el aspecto de la línea de significancia al 95%
make_wwz_plot <- function(df, titulo, zlim,
                          contour_col = "white",
                          contour_lwd = 0.7) {
  # sig_num convierte el logical sig (TRUE/FALSE) a numérico para geom_contour.
  # breaks = 0.5 traza exactamente la frontera entre zona significativa (1) y no (0).
  df$sig_num <- as.numeric(df$sig)
  ggplot(df, aes(x = tau, y = freq)) +
    geom_tile(aes(fill = WWZ_plot, alpha = coi)) +
    scale_fill_viridis(option = "C", limits = zlim, name = "WWZ") +
    scale_alpha(range = c(0.2, 1), guide = "none") +
    # Línea de contorno del nivel de significancia al 95% (Test de Monte Carlo)
    geom_contour(aes(z = sig_num),
                 breaks    = 0.5,
                 colour    = contour_col,
                 linewidth = contour_lwd) +
    labs(title = titulo, x = "Tiempo", y = "Frecuencia") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))
}

# FUNCIÓN AUXILIAR create_missing — definida una sola vez, usada en todas las pruebas
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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. PRUEBA SINTÉTICA: UN SENO, VARIOS % DE DATOS ELIMINADOS ####
# Genera sin(2*pi*t/p), elimina distintos porcentajes y compara
# visualmente las series y los WWZ resultantes.
# (C) 2026. Ricardo Sánchez Santos & Josue M. Polanco-Martinez
# GNU General Public License v3 — https://www.gnu.org/licenses/
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. PARÁMETROS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# N, p1, p2 definidos globalmente al inicio del script
pct <- c(10, 30, 50, 90)

t <- seq(0, 200, length.out = N)
x <- sin(2 * pi * t / p1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. GENERAR SERIES CON DATOS FALTANTES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("\n--- Generando series con datos faltantes (Prueba 1) ---\n")
resultados_list <- lapply(pct, function(p_val) {
  create_missing(x, t, porcentaje = p_val)
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3A. PLOTS SERIES TEMPORALES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_completa <- data.frame(t = t, x = x)

p_ts_completa <- ggplot(df_completa, aes(t, x)) +
  geom_line(color = "red", linewidth = 0.7) +
  geom_point(size = 0.8, color = "red", alpha = 0.6) +
  labs(title    = "Serie completa",
       subtitle = paste0("N = ", N, " puntos | p = ", p1),
       x = "Tiempo", y = "Amplitud") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey40"))

plots_ts <- lapply(seq_along(resultados_list), function(j) {
  res_j <- resultados_list[[j]]
  df_j  <- data.frame(t = t, x = res_j$x_na)
  ggplot(df_j, aes(t, x)) +
    geom_line(color = "red", linewidth = 0.7) +
    geom_point(size = 0.8, color = "red", alpha = 0.6, na.rm = TRUE) +
    labs(title    = paste0("Faltan el ", pct[j], "% de los datos"),
         subtitle = paste0("N observados = ", length(res_j$t_obs),
                           " | N eliminados = ", res_j$n_missing),
         x = "Tiempo", y = "Amplitud") +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, color = "grey40"))
})

cat("\n--- Comparativa series temporales (Prueba 1) ---\n")
print(p_ts_completa)
save_plot("p2_completa.png", p_ts_completa)
p_ts_grid <- arrangeGrob(grobs = plots_ts, ncol = 2)
grid.arrange(grobs = plots_ts, ncol = 2)
save_plot("p2_incompletas.png", p_ts_grid)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3B. WWZ
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

freq_wwz <- seq(1/p2 * 0.5, 1/p1 * 2, length.out = 100)

cat("\nEjecutando WWZ serie completa (Prueba 1)...\n")
res_completa <- wwz_robust(t, x,
                           freq   = freq_wwz,
                           use_mc =FALSE,
                           show_plot = FALSE)

cat("Ejecutando WWZ series incompletas (Prueba 1)...\n")
res_wwz_list <- lapply(resultados_list, function(lst_r) {
  wwz_robust(lst_r$t_obs, lst_r$x_obs,
             freq   = freq_wwz,
             use_mc = TRUE,
             n_sim=1000,
             show_plot = FALSE)
})


zlim <- range(c(res_completa$data$WWZ_plot,
                unlist(lapply(res_wwz_list, function(lst_r) lst_r$data$WWZ_plot))),
              na.rm = TRUE)

p_wwz1 <- make_wwz_plot(res_completa$data, "WWZ — Serie completa", zlim)

plots_wwz <- lapply(seq_along(res_wwz_list), function(j) {
  make_wwz_plot(res_wwz_list[[j]]$data,
                paste0("WWZ — Serie con ", pct[j], "% datos eliminados"),
                zlim)
})

print(p_wwz1)
save_plot("p2_wwz_completaMC.png", p_wwz1)
p_wwz_grid <- arrangeGrob(grobs = plots_wwz, ncol = 2)
grid.arrange(grobs = plots_wwz, ncol = 2)
save_plot("p2_wwz_incompletasMC.png", p_wwz_grid)

#Aparte de las conclusiones que sacaremos a la hora de ir quitando datos
#Comentar la diferencia entre poner un número de datos u otro. Cuando pongo
# 500 puntos, la función va bien con un periodo de 23000 años. Pero con 1000
#datos la función completa deja de funcionar. Esto se debe a que al poner un 
#número de datos tan elevado,no hay casi distancia entre cada punto. Al calcular
#el AR1, el AR1 calcula un sigma muy alto ya que los puntos están juntos pero lo
#suficientemente separados. Sin embargo con 1000 pasos, están muy juntos cada
#punto es prácticamente idéntico al anterior por lo que la sigma es muy cercana a 1
# y el AR1 identifica la función como una no estacionaria. Sin embargo, cuando le quitamos
#el 80% de los datos, los datos ya están lo suficientemente separados.


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. PRUEBA SINTÉTICA: DOS SENOS, UN % DE DATOS ELIMINADOS ####
# Genera sin(2*pi*t/p1) + sin(2*pi*t/p2), elimina un porcentaje y compara.
# (C) 2026. Ricardo Sánchez Santos & Josue M. Polanco-Martinez
# GNU General Public License v3 — https://www.gnu.org/licenses/
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. PARÁMETROS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# N, p1, p2 definidos globalmente al inicio del script
pct <- c(10, 30, 50, 90)

t <- seq(0, 200, length.out = N)
x <- sin(2 * pi * t / p1) + sin(2 * pi * t / p2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. GENERAR SERIES CON DATOS FALTANTES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("\n--- Generando series con datos faltantes (Prueba 2) ---\n")
resultados_list <- lapply(pct, function(p_val) {
  create_missing(x, t, porcentaje = p_val)
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3A. PLOTS SERIES TEMPORALES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_completa <- data.frame(t = t, x = x)

p_ts_completa <- ggplot(df_completa, aes(t, x)) +
  geom_line(color = "red", linewidth = 0.7) +
  geom_point(size = 0.8, color = "red", alpha = 0.6) +
  labs(title    = "Serie completa",
       subtitle = paste0("N = ", N, " puntos | p1 = ", p1, " | p2 = ", p2),
       x = "Tiempo", y = "Amplitud") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey40"))

plots_ts <- lapply(seq_along(resultados_list), function(j) {
  res_j <- resultados_list[[j]]
  df_j  <- data.frame(t = t, x = res_j$x_na)
  ggplot(df_j, aes(t, x)) +
    geom_line(color = "red", linewidth = 0.7) +
    geom_point(size = 0.8, color = "red", alpha = 0.6, na.rm = TRUE) +
    labs(title    = paste0("Faltan el ", pct[j], "% de los datos"),
         subtitle = paste0("N observados = ", length(res_j$t_obs),
                           " | N eliminados = ", res_j$n_missing),
         x = "Tiempo", y = "Amplitud") +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, color = "grey40"))
})

cat("\n--- Comparativa series temporales (Prueba 1) ---\n")
print(p_ts_completa)
save_plot("p2_completa.png", p_ts_completa)
p_ts_grid <- arrangeGrob(grobs = plots_ts, ncol = 2)
grid.arrange(grobs = plots_ts, ncol = 2)
save_plot("p2_incompletas.png", p_ts_grid)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3B. WWZ
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

freq_wwz <- seq(1/p2 * 0.5, 1/p1 * 2, length.out = 100)

cat("\nEjecutando WWZ serie completa (Prueba 1)...\n")
res_completa <- wwz_robust(t, x,
                           freq   = freq_wwz,
                           use_mc =FALSE,
                           show_plot = FALSE)

cat("Ejecutando WWZ series incompletas (Prueba 1)...\n")
res_wwz_list <- lapply(resultados_list, function(lst_r) {
  wwz_robust(lst_r$t_obs, lst_r$x_obs,
             freq   = freq_wwz,
             use_mc = FALSE,
             show_plot = FALSE)
})

zlim <- range(c(res_completa$data$WWZ_plot,
                unlist(lapply(res_wwz_list, function(lst_r) lst_r$data$WWZ_plot))),
              na.rm = TRUE)

p_wwz1 <- make_wwz_plot(res_completa$data, "WWZ — Serie completa", zlim)

plots_wwz <- lapply(seq_along(res_wwz_list), function(j) {
  make_wwz_plot(res_wwz_list[[j]]$data,
                paste0("WWZ — Serie con ", pct[j], "% datos eliminados"),
                zlim)
})

print(p_wwz1)
save_plot("p2_wwz_completaTS.png", p_wwz1)
p_wwz_grid <- arrangeGrob(grobs = plots_wwz, ncol = 2)
grid.arrange(grobs = plots_wwz, ncol = 2)
save_plot("p2_wwz_incompletasTS.png", p_wwz_grid)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. PRUEBA SINTÉTICA: UN SENO + RUIDO ROJO, UN % DE DATOS ELIMINADOS ####
# Genera sin(2*pi*t/p) + AR1, elimina un porcentaje y compara.
# (C) 2026. Ricardo Sánchez Santos & Josue M. Polanco-Martinez
# GNU General Public License v3 — https://www.gnu.org/licenses/
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. PARÁMETROS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# N, p1 definidos globalmente al inicio del script
pct <- c(10, 30, 50, 90)

t        <- seq(0, 200, length.out = N)
sd_senal <- sd(sin(2 * pi * t / p1))
ruido    <- as.numeric(arima.sim(model = list(ar = 0.7), n = N, sd = 0.5 * sd_senal))
x        <- sin(2 * pi * t / p1) + ruido

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. GENERAR SERIES CON DATOS FALTANTES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("\n--- Generando series con datos faltantes (Prueba 3) ---\n")
resultados_list <- lapply(pct, function(p_val) {
  create_missing(x, t, porcentaje = p_val)
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3A. PLOTS SERIES TEMPORALES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_completa <- data.frame(t = t, x = x)

p_ts_completa <- ggplot(df_completa, aes(t, x)) +
  geom_line(color = "red", linewidth = 0.7) +
  geom_point(size = 0.8, color = "red", alpha = 0.6) +
  labs(title    = "Serie completa",
       subtitle = paste0("N = ", N, " puntos | p = ", p1),
       x = "Tiempo", y = "Amplitud") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey40"))

plots_ts <- lapply(seq_along(resultados_list), function(j) {
  res_j <- resultados_list[[j]]
  df_j  <- data.frame(t = t, x = res_j$x_na)
  ggplot(df_j, aes(t, x)) +
    geom_line(color = "red", linewidth = 0.7) +
    geom_point(size = 0.8, color = "red", alpha = 0.6, na.rm = TRUE) +
    labs(title    = paste0("Faltan el ", pct[j], "% de los datos"),
         subtitle = paste0("N observados = ", length(res_j$t_obs),
                           " | N eliminados = ", res_j$n_missing),
         x = "Tiempo", y = "Amplitud") +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, color = "grey40"))
})

cat("\n--- Comparativa series temporales (Prueba 1) ---\n")
print(p_ts_completa)
save_plot("p3_completa.png", p_ts_completa)
p_ts_grid <- arrangeGrob(grobs = plots_ts, ncol = 2)
grid.arrange(grobs = plots_ts, ncol = 2)
save_plot("p3_incompletas.png", p_ts_grid)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3B. WWZ
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

freq_wwz <- seq(1/p1 * 0.5, 1/p1 * 2, length.out = 100)

cat("\nEjecutando WWZ serie completa (Prueba 3)...\n")
res_completa <- wwz_robust(t, x,
                           freq   = freq_wwz,
                           use_mc =FALSE,
                           show_plot = FALSE)

cat("Ejecutando WWZ series incompletas (Prueba 3)...\n")
res_wwz_list <- lapply(resultados_list, function(lst_r) {
  wwz_robust(lst_r$t_obs, lst_r$x_obs,
             freq   = freq_wwz,
             use_mc = FALSE,
             show_plot = FALSE)
})

zlim <- range(c(res_completa$data$WWZ_plot,
                unlist(lapply(res_wwz_list, function(lst_r) lst_r$data$WWZ_plot))),
              na.rm = TRUE)

p_wwz1 <- make_wwz_plot(res_completa$data, "WWZ — Serie completa", zlim)

plots_wwz <- lapply(seq_along(res_wwz_list), function(j) {
  make_wwz_plot(res_wwz_list[[j]]$data,
                paste0("WWZ — Serie con ", pct[j], "% datos eliminados"),
                zlim)
})

print(p_wwz1)
save_plot("p3_wwz_completaTS.png", p_wwz1)
p_wwz_grid <- arrangeGrob(grobs = plots_wwz, ncol = 2)
grid.arrange(grobs = plots_wwz, ncol = 2)
save_plot("p3_wwz_incompletasTS.png", p_wwz_grid)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. PRUEBA SINTÉTICA: DOS SENOS + RUIDO ROJO, UN % DE DATOS ELIMINADOS ####
# Genera sin(2*pi*t/p1) + sin(2*pi*t/p2) + AR1, elimina un porcentaje y compara.
# (C) 2026. Ricardo Sánchez Santos & Josue M. Polanco-Martinez
# GNU General Public License v3 — https://www.gnu.org/licenses/
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. PARÁMETROS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# N, p1, p2 definidos globalmente al inicio del script
pct <- c(10, 30, 50, 90)

t        <- seq(0, 200, length.out = N)
sd_senal <- sd(sin(2 * pi * t / p1) + sin(2 * pi * t / p2))
ruido    <- as.numeric(arima.sim(model = list(ar = 0.7), n = N, sd = 0.5 * sd_senal))
x        <- sin(2 * pi * t / p1) + sin(2 * pi * t / p2) + ruido


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. GENERAR SERIES CON DATOS FALTANTES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("\n--- Generando series con datos faltantes (Prueba 4) ---\n")
resultados_list <- lapply(pct, function(p_val) {
  create_missing(x, t, porcentaje = p_val)
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3A. PLOTS SERIES TEMPORALES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_completa <- data.frame(t = t, x = x)

p_ts_completa <- ggplot(df_completa, aes(t, x)) +
  geom_line(color = "red", linewidth = 0.7) +
  geom_point(size = 0.8, color = "red", alpha = 0.6) +
  labs(title    = "Serie completa",
       subtitle = paste0("N = ", N, " puntos | p1 = ", p1, " | p2 = ", p2),
       x = "Tiempo", y = "Amplitud") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey40"))

plots_ts <- lapply(seq_along(resultados_list), function(j) {
  res_j <- resultados_list[[j]]
  df_j  <- data.frame(t = t, x = res_j$x_na)
  ggplot(df_j, aes(t, x)) +
    geom_line(color = "red", linewidth = 0.7) +
    geom_point(size = 0.8, color = "red", alpha = 0.6, na.rm = TRUE) +
    labs(title    = paste0("Faltan el ", pct[j], "% de los datos"),
         subtitle = paste0("N observados = ", length(res_j$t_obs),
                           " | N eliminados = ", res_j$n_missing),
         x = "Tiempo", y = "Amplitud") +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, color = "grey40"))
})

cat("\n--- Comparativa series temporales (Prueba 4) ---\n")
print(p_ts_completa)
save_plot("p4_completa.png", p_ts_completa)
p_ts_grid <- arrangeGrob(grobs = plots_ts, ncol = 2)
grid.arrange(grobs = plots_ts, ncol = 2)
save_plot("p4_incompletas.png", p_ts_grid)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3B. WWZ
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

freq_wwz <- seq(1/p2 * 0.5, 1/p1 * 2, length.out = 100)

cat("\nEjecutando WWZ serie completa (Prueba 1)...\n")
res_completa <- wwz_robust(t, x,
                           freq   = freq_wwz,
                           use_mc =FALSE,,
                           show_plot = FALSE)

cat("Ejecutando WWZ series incompletas (Prueba 4)...\n")
res_wwz_list <- lapply(resultados_list, function(lst_r) {
  wwz_robust(lst_r$t_obs, lst_r$x_obs,
             freq   = freq_wwz,
             use_mc = FALSE,
             show_plot = FALSE)
})



zlim <- range(c(res_completa$data$WWZ_plot,
                unlist(lapply(res_wwz_list, function(lst_r) lst_r$data$WWZ_plot))),
              na.rm = TRUE)

p_wwz1 <- make_wwz_plot(res_completa$data, "WWZ — Serie completa", zlim)

plots_wwz <- lapply(seq_along(res_wwz_list), function(j) {
  make_wwz_plot(res_wwz_list[[j]]$data,
                paste0("WWZ — Serie con ", pct[j], "% datos eliminados"),
                zlim)
})

print(p_wwz1)
save_plot("p4_wwz_completaTS.png", p_wwz1)
p_wwz_grid <- arrangeGrob(grobs = plots_wwz, ncol = 2)
grid.arrange(grobs = plots_wwz, ncol = 2)
save_plot("p4_wwz_incompletasTS.png", p_wwz_grid)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. PRUEBA SINTÉTICA: UN SENO + RUIDO BLANCO, UN % DE DATOS ELIMINADOS ####
# Genera sin(2*pi*t/p) + Normal, elimina un porcentaje y compara.
# (C) 2026. Ricardo Sánchez Santos & Josue M. Polanco-Martinez
# GNU General Public License v3 — https://www.gnu.org/licenses/
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. PARÁMETROS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# N, p1 definidos globalmente al inicio del script
pct <- c(10, 30, 50, 90)

t        <- seq(0, 200, length.out = N)
sd_senal <- sd(sin(2 * pi * t / p1))
ruido    <- rnorm(N, mean = 0, sd = 0.5 * sd_senal)
x        <- sin(2 * pi * t / p1) + ruido


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. GENERAR SERIES CON DATOS FALTANTES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("\n--- Generando series con datos faltantes (Prueba 5) ---\n")
resultados_list <- lapply(pct, function(p_val) {
  create_missing(x, t, porcentaje = p_val)
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3A. PLOTS SERIES TEMPORALES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_completa <- data.frame(t = t, x = x)

p_ts_completa <- ggplot(df_completa, aes(t, x)) +
  geom_line(color = "red", linewidth = 0.7) +
  geom_point(size = 0.8, color = "red", alpha = 0.6) +
  labs(title    = "Serie completa",
       subtitle = paste0("N = ", N, " puntos | p = ", p1),
       x = "Tiempo", y = "Amplitud") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey40"))

plots_ts <- lapply(seq_along(resultados_list), function(j) {
  res_j <- resultados_list[[j]]
  df_j  <- data.frame(t = t, x = res_j$x_na)
  ggplot(df_j, aes(t, x)) +
    geom_line(color = "red", linewidth = 0.7) +
    geom_point(size = 0.8, color = "red", alpha = 0.6, na.rm = TRUE) +
    labs(title    = paste0("Faltan el ", pct[j], "% de los datos"),
         subtitle = paste0("N observados = ", length(res_j$t_obs),
                           " | N eliminados = ", res_j$n_missing),
         x = "Tiempo", y = "Amplitud") +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, color = "grey40"))
})

cat("\n--- Comparativa series temporales (Prueba 5) ---\n")
print(p_ts_completa)
save_plot("p5_completa.png", p_ts_completa)
p_ts_grid <- arrangeGrob(grobs = plots_ts, ncol = 2)
grid.arrange(grobs = plots_ts, ncol = 2)
save_plot("p5_incompletas.png", p_ts_grid)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3B. WWZ
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

freq_wwz <- seq(1/p1 * 0.5, 1/p1 * 2, length.out = 100)

cat("\nEjecutando WWZ serie completa (Prueba 5)...\n")
res_completa <- wwz_robust(t, x,
                           freq   = freq_wwz,
                           use_mc =FALSE,
                           show_plot = FALSE)

cat("Ejecutando WWZ series incompletas (Prueba 5)...\n")
res_wwz_list <- lapply(resultados_list, function(lst_r) {
  wwz_robust(lst_r$t_obs, lst_r$x_obs,
             freq   = freq_wwz,
             use_mc = FALSE,
             show_plot = FALSE)
})


zlim <- range(c(res_completa$data$WWZ_plot,
                unlist(lapply(res_wwz_list, function(lst_r) lst_r$data$WWZ_plot))),
              na.rm = TRUE)

p_wwz1 <- make_wwz_plot(res_completa$data, "WWZ — Serie completa", zlim)

plots_wwz <- lapply(seq_along(res_wwz_list), function(j) {
  make_wwz_plot(res_wwz_list[[j]]$data,
                paste0("WWZ — Serie con ", pct[j], "% datos eliminados"),
                zlim)
})

print(p_wwz1)
save_plot("p5_wwz_completaTS.png", p_wwz1)
p_wwz_grid <- arrangeGrob(grobs = plots_wwz, ncol = 2)
grid.arrange(grobs = plots_wwz, ncol = 2)
save_plot("p5_wwz_incompletasTS.png", p_wwz_grid)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6. PRUEBA SINTÉTICA: TRES SENOS + UN % DE DATOS ELIMINADOS ####
# Genera  + Normal, elimina un porcentaje y compara.
# (C) 2026. Ricardo Sánchez Santos & Josue M. Polanco-Martinez
# GNU General Public License v3 — https://www.gnu.org/licenses/
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. PARÁMETROS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# N, p1 definidos globalmente al inicio del script
pct <- c(10, 30, 50, 90)

p3<-5
p4<-25
p5<-100
t        <- seq(0, 200, length.out = N)
x        <- sin(2 * pi * t / p3) + sin(2 * pi * t / p4)+sin(2 * pi * t / p5)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. GENERAR SERIES CON DATOS FALTANTES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("\n--- Generando series con datos faltantes (Prueba 5) ---\n")
resultados_list <- lapply(pct, function(p_val) {
  create_missing(x, t, porcentaje = p_val)
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3A. PLOTS SERIES TEMPORALES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_completa <- data.frame(t = t, x = x)

p_ts_completa <- ggplot(df_completa, aes(t, x)) +
  geom_line(color = "red", linewidth = 0.7) +
  geom_point(size = 0.8, color = "red", alpha = 0.6) +
  labs(title    = "Serie completa",
       subtitle = paste0("N = ", N, " puntos | p = ", p3,",",p4,",",p5),
       x = "Tiempo", y = "Amplitud") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey40"))

plots_ts <- lapply(seq_along(resultados_list), function(j) {
  res_j <- resultados_list[[j]]
  df_j  <- data.frame(t = t, x = res_j$x_na)
  ggplot(df_j, aes(t, x)) +
    geom_line(color = "red", linewidth = 0.7) +
    geom_point(size = 0.8, color = "red", alpha = 0.6, na.rm = TRUE) +
    labs(title    = paste0("Faltan el ", pct[j], "% de los datos"),
         subtitle = paste0("N observados = ", length(res_j$t_obs),
                           " | N eliminados = ", res_j$n_missing),
         x = "Tiempo", y = "Amplitud") +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, color = "grey40"))
})

cat("\n--- Comparativa series temporales (Prueba 5) ---\n")
print(p_ts_completa)
save_plot("p6_completa.png", p_ts_completa)
p_ts_grid <- arrangeGrob(grobs = plots_ts, ncol = 2)
grid.arrange(grobs = plots_ts, ncol = 2)
save_plot("p6_incompletas.png", p_ts_grid)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3B. WWZ
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

freq_wwz <- seq(1/p5 * 0.5, 1/p3 * 2, length.out = 100)

cat("\nEjecutando WWZ serie completa (Prueba 6)...\n")
res_completa <- wwz_robust(t, x,
                           freq   = freq_wwz,
                           use_mc =FALSE,
                           show_plot = FALSE)

cat("Ejecutando WWZ series incompletas (Prueba 6)...\n")
res_wwz_list <- lapply(resultados_list, function(lst_r) {
  wwz_robust(lst_r$t_obs, lst_r$x_obs,
             freq   = freq_wwz,
             use_mc = FALSE,
             show_plot = FALSE)
})


zlim <- range(c(res_completa$data$WWZ_plot,
                unlist(lapply(res_wwz_list, function(lst_r) lst_r$data$WWZ_plot))),
              na.rm = TRUE)

p_wwz1 <- make_wwz_plot(res_completa$data, "WWZ — Serie completa", zlim)

plots_wwz <- lapply(seq_along(res_wwz_list), function(j) {
  make_wwz_plot(res_wwz_list[[j]]$data,
                paste0("WWZ — Serie con ", pct[j], "% datos eliminados"),
                zlim)
})

print(p_wwz1)
save_plot("p6_wwz_completaTS.png", p_wwz1)
p_wwz_grid <- arrangeGrob(grobs = plots_wwz, ncol = 2)
grid.arrange(grobs = plots_wwz, ncol = 2)
save_plot("p6_wwz_incompletasTS.png", p_wwz_grid)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7. COMPARATIVA: RUIDO ROJO vs RUIDO BLANCO (serie completa, sin datos faltantes) ####
# sin(2*pi*t/p) + AR1 vs sin(2*pi*t/p) + Normal
# (C) 2026. Ricardo Sánchez Santos & Josue M. Polanco-Martinez
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. PARÁMETROS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# N, p1 definidos globalmente al inicio del script
t <- seq(0, 200, length.out = N)

set.seed(123)
sd_senal <- sd(sin(2 * pi * t / p1))

senal_base   <- sin(2 * pi * t / p1)
ruido_rojo   <- as.numeric(arima.sim(model = list(ar = 0.7), n = N, sd = 0.5 * sd_senal))
ruido_blanco <- rnorm(N, mean = 0, sd = 0.5 * sd_senal)

x_rojo   <- senal_base + ruido_rojo
x_blanco <- senal_base + ruido_blanco

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. COMPARATIVA SERIES TEMPORALES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_rojo   <- data.frame(t = t, x = x_rojo)
df_blanco <- data.frame(t = t, x = x_blanco)

p_ts_rojo <- ggplot(df_rojo, aes(t, x)) +
  geom_line(color = "firebrick", linewidth = 0.7) +
  geom_point(size = 0.6, color = "firebrick", alpha = 0.5) +
  labs(title    = "Serie completa — Ruido rojo (AR1, phi = 0.7)",
       subtitle = paste0("p = ", p1, " | N = ", N),
       x = "Tiempo", y = "Amplitud") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey40"))

p_ts_blanco <- ggplot(df_blanco, aes(t, x)) +
  geom_line(color = "steelblue", linewidth = 0.7) +
  geom_point(size = 0.6, color = "steelblue", alpha = 0.5) +
  labs(title    = "Serie completa — Ruido blanco (Normal)",
       subtitle = paste0("p = ", p1, " | N = ", N),
       x = "Tiempo", y = "Amplitud") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey40"))

cat("\n--- Comparativa series temporales: ruido rojo vs ruido blanco ---\n")
plot1<-grid.arrange(p_ts_rojo, p_ts_blanco, ncol = 1,
                    top = grid::textGrob(
                      paste0("Comparativa ruido rojo vs ruido blanco (p = ", p1, ")"),
                      gp = grid::gpar(fontsize = 13, fontface = "bold")))

save_plot("p7plot.png", plot1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. WWZ DE CADA UNA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

freq_wwz <- seq(1/p1 * 0.5, 1/p1 * 2, length.out = 100)

cat("\nEjecutando WWZ — Ruido rojo (Prueba 7)...\n")
res_rojo <- wwz_robust(t, x_rojo,
                       freq   = freq_wwz,
                       use_mc = TRUE,
                       n_sim  = 1000,
                       show_plot = FALSE)

cat("Ejecutando WWZ — Ruido blanco (Prueba 7)...\n")
res_blanco <- wwz_robust(t, x_blanco,
                         freq   = freq_wwz,
                         use_mc = TRUE,
                         n_sim  = 1000,
                         show_plot = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. PLOTS WWZ CON ESCALA COMÚN
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

zlim_comp <- range(c(res_rojo$data$WWZ_plot,
                     res_blanco$data$WWZ_plot), na.rm = TRUE)

p_wwz_rojo   <- make_wwz_plot(res_rojo$data,
                              paste0("WWZ — Ruido rojo (AR1, phi=0.7) | p = ", p1),
                              zlim_comp)
p_wwz_blanco <- make_wwz_plot(res_blanco$data,
                              paste0("WWZ — Ruido blanco (Normal) | p = ", p1),
                              zlim_comp)

pwwz<-grid.arrange(p_wwz_rojo, p_wwz_blanco, ncol = 1,
                   top = grid::textGrob(
                     "Comparativa WWZ: ruido rojo vs ruido blanco (serie completa)",
                     gp = grid::gpar(fontsize = 13, fontface = "bold")))
save_plot("p7_wwz_MC.png", pwwz)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 8. DEMOSTRACIÓN DE FALSOS POSITIVOS CON RUIDO ROJO ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

N_p8     <- 500                              # N local: no sobreescribe el global N=1000
t_p8     <- seq(0, 100, length.out = N_p8)
phi_vals <- c(0.3, 0.7, 0.95)
freq_wwz <- seq(1/60, 1/10, length.out = 100)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. BUCLE: calcular todo y guardar
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

resultados <- list()

for(i in seq_along(phi_vals)){
  
  phi <- phi_vals[i]
  set.seed(100 + i)
  x_ruido <- as.numeric(arima.sim(model = list(ar = phi), n = N_p8, sd = 1))
  phi_estimado <- round(arima(x_ruido, c(1,0,0))$coef["ar1"], 3)
  
  cat(sprintf("\nEjecutando WWZ — phi = %.2f (Prueba 8)...\n", phi))
  res_i <- wwz_robust(t_p8, x_ruido,
                      freq      = freq_wwz,
                      use_mc    = TRUE,
                      n_sim     = 50,
                      show_plot = FALSE)
  
  n_sig   <- sum(res_i$WWZ_significant, na.rm = TRUE)
  pct_sig <- round(100 * n_sig / length(res_i$WWZ_significant), 1)
  
  resultados[[i]] <- list(
    phi          = phi,
    phi_estimado = phi_estimado,
    x_ruido      = x_ruido,
    res          = res_i,
    n_sig        = n_sig,
    pct_sig      = pct_sig
  )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. PLOTS desde resultados guardados
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plots_ts  <- list()
plots_wwz <- list()

for(i in seq_along(phi_vals)){
  r_i <- resultados[[i]]
  
  df_ts <- data.frame(t = t_p8, x = r_i$x_ruido)
  
  plots_ts[[i]] <- ggplot(df_ts, aes(t, x)) +
    geom_line(color = "firebrick", linewidth = 0.6) +
    labs(title    = paste0("Ruido AR1 puro — phi real = ", r_i$phi,
                           " | phi estimado = ", r_i$phi_estimado),
         subtitle = "Sin ninguna señal senoidal",
         x = "Tiempo", y = "Amplitud") +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, color = "grey40"))
  
  plots_wwz[[i]] <- make_wwz_plot(
    r_i$res$data,
    paste0("WWZ — phi real = ", r_i$phi,
           " | Celdas significativas: ", r_i$n_sig, " (", r_i$pct_sig, "%)"),
    zlim = range(r_i$res$data$WWZ_plot, na.rm = TRUE)
  )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. PANELES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

grid.arrange(grobs = plots_ts, ncol = 1,
             top = grid::textGrob(
               "Ruido rojo puro (sin señal): phi = 0.3, 0.7, 0.95",
               gp = grid::gpar(fontsize = 13, fontface = "bold")))

grid.arrange(grobs = plots_wwz, ncol = 1,
             top = grid::textGrob(
               "WWZ sobre ruido rojo puro — falsos positivos con MC",
               gp = grid::gpar(fontsize = 13, fontface = "bold")))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. RESUMEN NUMÉRICO
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("\n--- Resumen: falsos positivos por phi ---\n")
cat(sprintf("%-12s %-15s %-15s %-20s\n",
            "phi real", "phi estimado", "Celdas sig.", "% significativas"))
cat(strrep("-", 65), "\n")

for(i in seq_along(phi_vals)){
  r_i <- resultados[[i]]
  cat(sprintf("%-12.2f %-15.3f %-15d %-20s\n",
              r_i$phi, r_i$phi_estimado, r_i$n_sig, paste0(r_i$pct_sig, "%")))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 9. FALSOS POSITIVOS SOBRE SEÑAL SINUSOIDAL ####
# sin(2*pi*t/p) + AR1 con phi creciente
# Demostración de frecuencias espurias
# (C) 2026. Ricardo Sánchez Santos & Josue M. Polanco-Martinez
# GNU General Public License v3 — https://www.gnu.org/licenses/
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. PARÁMETROS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# p1 definido globalmente; N local para esta prueba
N_p9     <- 500                              # N local: no sobreescribe el global N=1000
t_p9     <- seq(0, 100, length.out = N_p9)
phi_vals <- c(0.3, 0.7, 0.95)
freq_wwz <- seq(1/60, 1/10, length.out = 100)

senal_base <- sin(2 * pi * t_p9 / p1)
sd_senal   <- sd(senal_base)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. BUCLE: calcular todo y guardar
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

resultados <- list()

for(i in seq_along(phi_vals)){
  
  phi <- phi_vals[i]
  set.seed(100 + i)
  ruido <- as.numeric(arima.sim(model = list(ar = phi), n = N_p9,
                                sd = 2.0 * sd_senal))
  x <- senal_base + ruido
  phi_estimado <- round(arima(x, c(1,0,0))$coef["ar1"], 3)
  
  cat(sprintf("\nEjecutando WWZ — phi = %.2f (Prueba 9)...\n", phi))
  res_i <- wwz_robust(t_p9, x,
                      freq      = freq_wwz,
                      use_mc    = TRUE,
                      n_sim     = 50,
                      show_plot = FALSE)
  
  n_sig    <- sum(res_i$WWZ_significant, na.rm = TRUE)
  pct_sig  <- round(100 * n_sig / length(res_i$WWZ_significant), 1)
  f_real   <- 1 / p1
  df_sig   <- res_i$data[!is.na(res_i$data$sig) & res_i$data$sig, ]
  n_falsas <- sum(abs(df_sig$freq - f_real) > (f_real * 0.2), na.rm = TRUE)
  
  resultados[[i]] <- list(
    phi          = phi,
    phi_estimado = phi_estimado,
    x            = x,
    ruido        = ruido,
    res          = res_i,
    n_sig        = n_sig,
    pct_sig      = pct_sig,
    n_falsas     = n_falsas
  )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. PLOTS desde resultados guardados
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plots_ts  <- list()
plots_wwz <- list()

for(i in seq_along(phi_vals)){
  r_i <- resultados[[i]]
  
  df_ts <- data.frame(
    t    = rep(t_p9, 2),
    x    = c(senal_base, r_i$x),
    tipo = rep(c("Señal pura", "Señal + ruido AR1"), each = N_p9)
  )
  
  plots_ts[[i]] <- ggplot(df_ts, aes(t, x, color = tipo)) +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = c("Señal pura"       = "steelblue",
                                  "Señal + ruido AR1" = "firebrick")) +
    labs(title    = paste0("phi = ", r_i$phi,
                           " | phi estimado = ", r_i$phi_estimado),
         subtitle = paste0("Ruido = 2x amplitud señal | p = ", p1),
         x = "Tiempo", y = "Amplitud",
         color = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title      = element_text(face = "bold", size = 11),
          plot.subtitle   = element_text(size = 9, color = "grey40"),
          legend.position = "bottom")
  
  plots_wwz[[i]] <- make_wwz_plot(
    r_i$res$data,
    paste0("WWZ — phi = ", r_i$phi,
           " | Sig. totales: ", r_i$n_sig,
           " | Frecuencias falsas: ", r_i$n_falsas),
    zlim = range(r_i$res$data$WWZ_plot, na.rm = TRUE)
  )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. PANELES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

grid.arrange(grobs = plots_ts, ncol = 1,
             top = grid::textGrob(
               paste0("Señal sin(2\u03c0 t/", p1, ") + ruido AR1 creciente"),
               gp = grid::gpar(fontsize = 13, fontface = "bold")))

grid.arrange(grobs = plots_wwz, ncol = 1,
             top = grid::textGrob(
               "WWZ — aparición de frecuencias espurias con ruido rojo fuerte",
               gp = grid::gpar(fontsize = 13, fontface = "bold")))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. RESUMEN NUMÉRICO
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("\n--- Resumen: señal real p =", p1, "| ruido = 2x señal ---\n")
cat(sprintf("%-12s %-15s %-15s %-20s %-15s\n",
            "phi real", "phi estimado", "Celdas sig.", "%  sig.", "Frec. falsas"))
cat(strrep("-", 80), "\n")

for(i in seq_along(phi_vals)){
  r_i <- resultados[[i]]
  cat(sprintf("%-12.2f %-15.3f %-15d %-20s %-15d\n",
              r_i$phi, r_i$phi_estimado, r_i$n_sig,
              paste0(r_i$pct_sig, "%"), r_i$n_falsas))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 10. ROBUSTEZ FRENTE A AUTOCORRELACIÓN PERSISTENTE (RUIDO ROJO) ####
# Evalúa la capacidad de la WWZ para detectar periodicidades reales
# bajo ruido AR(1) con distintos niveles de persistencia (phi).
# A medida que phi aumenta, el ruido concentra más potencia en bajas
# frecuencias, dificultando la discriminación señal/ruido.
# phi ∈ {0, 0.2, 0.5, 0.7, 0.95}
# (C) 2026. Ricardo Sánchez Santos & Josue M. Polanco-Martinez
# GNU General Public License v3 — https://www.gnu.org/licenses/
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. PARÁMETROS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# N, p1 definidos globalmente al inicio del script
N_p10    <- 1000
t_p10    <- seq(0, 200, length.out = N_p10)
phi_vals <- c(0, 0.2, 0.5, 0.7, 0.95)
freq_wwz <- seq(1/p1 * 0.5, 1/p1 * 2, length.out = 100)

senal_base <- sin(2 * pi * t_p10 / p1)
sd_senal   <- sd(senal_base)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. BUCLE: calcular todo y guardar
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("\n--- Generando series y WWZ (Prueba 10) ---\n")
resultados_p10 <- list()

for (i in seq_along(phi_vals)) {
  
  phi <- phi_vals[i]
  set.seed(200 + i)
  
  # AR(1) con phi=0 equivale a ruido blanco puro
  if (phi == 0) {
    ruido <- rnorm(N_p10, mean = 0, sd = sd_senal)
  } else {
    ruido <- as.numeric(arima.sim(model = list(ar = phi), n = N_p10,
                                  sd = sd_senal))
  }
  
  x <- senal_base + ruido
  phi_estimado <- round(arima(x, order = c(1, 0, 0), method = "CSS-ML")$coef["ar1"], 3)
  
  cat(sprintf("  Ejecutando WWZ — phi = %.2f (Prueba 10)...\n", phi))
  
  res_i <- wwz_robust(t_p10, x,
                      freq      = freq_wwz,
                      use_mc    = TRUE,
                      n_sim     = 1000,
                      show_plot = FALSE)
  
  n_sig    <- sum(res_i$WWZ_significant, na.rm = TRUE)
  pct_sig  <- round(100 * n_sig / length(res_i$WWZ_significant), 1)
  f_real   <- 1 / p1
  
  # Celdas significativas alejadas >20% de la frecuencia real → falsos positivos
  df_sig   <- res_i$data[!is.na(res_i$data$sig) & res_i$data$sig, ]
  n_falsas <- sum(abs(df_sig$freq - f_real) > (f_real * 0.2), na.rm = TRUE)
  
  resultados_p10[[i]] <- list(
    phi          = phi,
    phi_estimado = phi_estimado,
    x            = x,
    ruido        = ruido,
    res          = res_i,
    n_sig        = n_sig,
    pct_sig      = pct_sig,
    n_falsas     = n_falsas
  )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3A. PLOTS SERIES TEMPORALES — uno por uno
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("\n--- Series temporales (Prueba 10) ---\n")

# phi = 0 (ruido blanco)
r <- resultados_p10[[1]]
df_ts <- data.frame(t = rep(t_p10, 2), x = c(senal_base, r$x),
                    tipo = rep(c("Señal pura", "Señal + ruido AR1"), each = N_p10))
p_ts_phi0 <- ggplot(df_ts, aes(t, x, color = tipo)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("Señal pura" = "steelblue", "Señal + ruido AR1" = "firebrick")) +
  labs(title = paste0("phi = 0 (ruido blanco) | phi estimado = ", r$phi_estimado),
       subtitle = paste0("Ruido = 1x amplitud señal | p = ", p1),
       x = "Tiempo", y = "Amplitud", color = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        legend.position = "bottom")
print(p_ts_phi0)
save_plot("p10_serie_phi0.png", p_ts_phi0, width = 10, height = 5)

# phi = 0.2
r <- resultados_p10[[2]]
df_ts <- data.frame(t = rep(t_p10, 2), x = c(senal_base, r$x),
                    tipo = rep(c("Señal pura", "Señal + ruido AR1"), each = N_p10))
p_ts_phi02 <- ggplot(df_ts, aes(t, x, color = tipo)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("Señal pura" = "steelblue", "Señal + ruido AR1" = "firebrick")) +
  labs(title = paste0("phi = 0.2 | phi estimado = ", r$phi_estimado),
       subtitle = paste0("Ruido = 1x amplitud señal | p = ", p1),
       x = "Tiempo", y = "Amplitud", color = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        legend.position = "bottom")
print(p_ts_phi02)
save_plot("p10_serie_phi02.png", p_ts_phi02, width = 10, height = 5)

# phi = 0.5
r <- resultados_p10[[3]]
df_ts <- data.frame(t = rep(t_p10, 2), x = c(senal_base, r$x),
                    tipo = rep(c("Señal pura", "Señal + ruido AR1"), each = N_p10))
p_ts_phi05 <- ggplot(df_ts, aes(t, x, color = tipo)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("Señal pura" = "steelblue", "Señal + ruido AR1" = "firebrick")) +
  labs(title = paste0("phi = 0.5 | phi estimado = ", r$phi_estimado),
       subtitle = paste0("Ruido = 1x amplitud señal | p = ", p1),
       x = "Tiempo", y = "Amplitud", color = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        legend.position = "bottom")
print(p_ts_phi05)
save_plot("p10_serie_phi05.png", p_ts_phi05, width = 10, height = 5)

# phi = 0.7
r <- resultados_p10[[4]]
df_ts <- data.frame(t = rep(t_p10, 2), x = c(senal_base, r$x),
                    tipo = rep(c("Señal pura", "Señal + ruido AR1"), each = N_p10))
p_ts_phi07 <- ggplot(df_ts, aes(t, x, color = tipo)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("Señal pura" = "steelblue", "Señal + ruido AR1" = "firebrick")) +
  labs(title = paste0("phi = 0.7 | phi estimado = ", r$phi_estimado),
       subtitle = paste0("Ruido = 1x amplitud señal | p = ", p1),
       x = "Tiempo", y = "Amplitud", color = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        legend.position = "bottom")
print(p_ts_phi07)
save_plot("p10_serie_phi07.png", p_ts_phi07, width = 10, height = 5)

# phi = 0.95
r <- resultados_p10[[5]]
df_ts <- data.frame(t = rep(t_p10, 2), x = c(senal_base, r$x),
                    tipo = rep(c("Señal pura", "Señal + ruido AR1"), each = N_p10))
p_ts_phi095 <- ggplot(df_ts, aes(t, x, color = tipo)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("Señal pura" = "steelblue", "Señal + ruido AR1" = "firebrick")) +
  labs(title = paste0("phi = 0.95 | phi estimado = ", r$phi_estimado),
       subtitle = paste0("Ruido = 1x amplitud señal | p = ", p1),
       x = "Tiempo", y = "Amplitud", color = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        legend.position = "bottom")
print(p_ts_phi095)
save_plot("p10_serie_phi095.png", p_ts_phi095, width = 10, height = 5)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3B. PLOTS WWZ — uno por uno
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

zlim_p10 <- range(
  unlist(lapply(resultados_p10, function(r) r$res$data$WWZ_plot)),
  na.rm = TRUE
)

cat("\n--- WWZ (Prueba 10) ---\n")

# phi = 0 (ruido blanco)
r <- resultados_p10[[1]]
p_wwz_phi0 <- make_wwz_plot(r$res$data,
                            paste0("WWZ — phi = 0 (ruido blanco) | Sig.: ", r$n_sig, " (", r$pct_sig, "%) | Falsas: ", r$n_falsas),
                            zlim_p10)
print(p_wwz_phi0)
save_plot("p10_wwz_phi0.png", p_wwz_phi0, width = 10, height = 6)

# phi = 0.2
r <- resultados_p10[[2]]
p_wwz_phi02 <- make_wwz_plot(r$res$data,
                             paste0("WWZ — phi = 0.2 | Sig.: ", r$n_sig, " (", r$pct_sig, "%) | Falsas: ", r$n_falsas),
                             zlim_p10)
print(p_wwz_phi02)
save_plot("p10_wwz_phi02.png", p_wwz_phi02, width = 10, height = 6)

# phi = 0.5
r <- resultados_p10[[3]]
p_wwz_phi05 <- make_wwz_plot(r$res$data,
                             paste0("WWZ — phi = 0.5 | Sig.: ", r$n_sig, " (", r$pct_sig, "%) | Falsas: ", r$n_falsas),
                             zlim_p10)
print(p_wwz_phi05)
save_plot("p10_wwz_phi05.png", p_wwz_phi05, width = 10, height = 6)

# phi = 0.7
r <- resultados_p10[[4]]
p_wwz_phi07 <- make_wwz_plot(r$res$data,
                             paste0("WWZ — phi = 0.7 | Sig.: ", r$n_sig, " (", r$pct_sig, "%) | Falsas: ", r$n_falsas),
                             zlim_p10)
print(p_wwz_phi07)
save_plot("p10_wwz_phi07.png", p_wwz_phi07, width = 10, height = 6)

# phi = 0.95
r <- resultados_p10[[5]]
p_wwz_phi095 <- make_wwz_plot(r$res$data,
                              paste0("WWZ — phi = 0.95 | Sig.: ", r$n_sig, " (", r$pct_sig, "%) | Falsas: ", r$n_falsas),
                              zlim_p10)
print(p_wwz_phi095)
save_plot("p10_wwz_phi095.png", p_wwz_phi095, width = 10, height = 6)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 11. PRUEBA SINTÉTICA: TRES SENO + RUIDO ROJO, UN % DE DATOS ELIMINADOS ####
# Genera sin(2*pi*t/p) + AR1, elimina un porcentaje y compara.
# (C) 2026. Ricardo Sánchez Santos & Josue M. Polanco-Martinez
# GNU General Public License v3 — https://www.gnu.org/licenses/
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. PARÁMETROS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# N, p1 definidos globalmente al inicio del script
pct <- c(10, 30, 50, 90)

p3<-5
p4<-25
p5<-100
t        <- seq(0, 200, length.out = N)
sd_senal <- sd(sin(2 * pi * t / p3) + sin(2 * pi * t / p4)+sin(2 * pi * t / p5))
ruido    <- as.numeric(arima.sim(model = list(ar = 0.7), n = N, sd = 0.5 * sd_senal))
x        <- sin(2 * pi * t / p3) + sin(2 * pi * t / p4)+sin(2 * pi * t / p5)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. GENERAR SERIES CON DATOS FALTANTES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("\n--- Generando series con datos faltantes (Prueba 3) ---\n")
resultados_list <- lapply(pct, function(p_val) {
  create_missing(x, t, porcentaje = p_val)
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3A. PLOTS SERIES TEMPORALES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_completa <- data.frame(t = t, x = x)

p_ts_completa <- ggplot(df_completa, aes(t, x)) +
  geom_line(color = "red", linewidth = 0.7) +
  geom_point(size = 0.8, color = "red", alpha = 0.6) +
  labs(title    = "Serie completa",
       subtitle = paste0("N = ", N, " puntos | p = ", p1),
       x = "Tiempo", y = "Amplitud") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey40"))

plots_ts <- lapply(seq_along(resultados_list), function(j) {
  res_j <- resultados_list[[j]]
  df_j  <- data.frame(t = t, x = res_j$x_na)
  ggplot(df_j, aes(t, x)) +
    geom_line(color = "red", linewidth = 0.7) +
    geom_point(size = 0.8, color = "red", alpha = 0.6, na.rm = TRUE) +
    labs(title    = paste0("Faltan el ", pct[j], "% de los datos"),
         subtitle = paste0("N observados = ", length(res_j$t_obs),
                           " | N eliminados = ", res_j$n_missing),
         x = "Tiempo", y = "Amplitud") +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, color = "grey40"))
})

cat("\n--- Comparativa series temporales (Prueba 1) ---\n")
print(p_ts_completa)
save_plot("p11_completa.png", p_ts_completa)
p_ts_grid <- arrangeGrob(grobs = plots_ts, ncol = 2)
grid.arrange(grobs = plots_ts, ncol = 2)
save_plot("p11_incompletas.png", p_ts_grid)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3B. WWZ
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

freq_wwz <- seq(1/p5 * 0.5, 1/p3 * 2, length.out = 100)

cat("\nEjecutando WWZ serie completa (Prueba 3)...\n")
res_completa <- wwz_robust(t, x,
                           freq   = freq_wwz,
                           use_mc =TRUE,
                           n_sim = 1000,
                           show_plot = FALSE)

cat("Ejecutando WWZ series incompletas (Prueba 3)...\n")
res_wwz_list <- lapply(resultados_list, function(lst_r) {
  wwz_robust(lst_r$t_obs, lst_r$x_obs,
             freq   = freq_wwz,
             use_mc = TRUE,
             n_sim=1000,
             show_plot = FALSE)
})

zlim <- range(c(res_completa$data$WWZ_plot,
                unlist(lapply(res_wwz_list, function(lst_r) lst_r$data$WWZ_plot))),
              na.rm = TRUE)

p_wwz1 <- make_wwz_plot(res_completa$data, "WWZ — Serie completa", zlim)

plots_wwz <- lapply(seq_along(res_wwz_list), function(j) {
  make_wwz_plot(res_wwz_list[[j]]$data,
                paste0("WWZ — Serie con ", pct[j], "% datos eliminados"),
                zlim)
})

print(p_wwz1)
save_plot("p11_wwz_completaMC.png", p_wwz1)
p_wwz_grid <- arrangeGrob(grobs = plots_wwz, ncol = 2)
grid.arrange(grobs = plots_wwz, ncol = 2)
save_plot("p11_wwz_incompletasMC.png", p_wwz_grid)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 12. PRUEBA SINTÉTICA: TRES SENO + RUIDO BLANCO, UN % DE DATOS ELIMINADOS ####
# Genera sin(2*pi*t/p) + Normal, elimina un porcentaje y compara.
# (C) 2026. Ricardo Sánchez Santos & Josue M. Polanco-Martinez
# GNU General Public License v3 — https://www.gnu.org/licenses/
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. PARÁMETROS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# N, p1 definidos globalmente al inicio del script
pct <- c(10, 30, 50, 90)

p3<-5
p4<-25
p5<-100
t        <- seq(0, 200, length.out = N)
sd_senal <- sd(sin(2 * pi * t / p3) + sin(2 * pi * t / p4)+sin(2 * pi * t / p5))
ruido    <- rnorm(N, mean = 0, sd = 0.5 * sd_senal)
x        <- sin(2 * pi * t / p3) + sin(2 * pi * t / p4)+sin(2 * pi * t / p5) + ruido


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. GENERAR SERIES CON DATOS FALTANTES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("\n--- Generando series con datos faltantes (Prueba 5) ---\n")
resultados_list <- lapply(pct, function(p_val) {
  create_missing(x, t, porcentaje = p_val)
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3A. PLOTS SERIES TEMPORALES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_completa <- data.frame(t = t, x = x)

p_ts_completa <- ggplot(df_completa, aes(t, x)) +
  geom_line(color = "red", linewidth = 0.7) +
  geom_point(size = 0.8, color = "red", alpha = 0.6) +
  labs(title    = "Serie completa",
       subtitle = paste0("N = ", N, " puntos | p = ", p1),
       x = "Tiempo", y = "Amplitud") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey40"))

plots_ts <- lapply(seq_along(resultados_list), function(j) {
  res_j <- resultados_list[[j]]
  df_j  <- data.frame(t = t, x = res_j$x_na)
  ggplot(df_j, aes(t, x)) +
    geom_line(color = "red", linewidth = 0.7) +
    geom_point(size = 0.8, color = "red", alpha = 0.6, na.rm = TRUE) +
    labs(title    = paste0("Faltan el ", pct[j], "% de los datos"),
         subtitle = paste0("N observados = ", length(res_j$t_obs),
                           " | N eliminados = ", res_j$n_missing),
         x = "Tiempo", y = "Amplitud") +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, color = "grey40"))
})

cat("\n--- Comparativa series temporales (Prueba 5) ---\n")
print(p_ts_completa)
save_plot("p12_completa.png", p_ts_completa)
p_ts_grid <- arrangeGrob(grobs = plots_ts, ncol = 2)
grid.arrange(grobs = plots_ts, ncol = 2)
save_plot("p12_incompletas.png", p_ts_grid)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3B. WWZ
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

freq_wwz <- seq(1/p5 * 0.5, 1/p3 * 2, length.out = 100)

cat("\nEjecutando WWZ serie completa (Prueba 5)...\n")
res_completa <- wwz_robust(t, x,
                           freq   = freq_wwz,
                           use_mc =TRUE,
                           n_sim=1000,
                           show_plot = FALSE)

cat("Ejecutando WWZ series incompletas (Prueba 5)...\n")
res_wwz_list <- lapply(resultados_list, function(lst_r) {
  wwz_robust(lst_r$t_obs, lst_r$x_obs,
             freq   = freq_wwz,
             use_mc = TRUE,
             n_sim=1000,
             show_plot = FALSE)
})


zlim <- range(c(res_completa$data$WWZ_plot,
                unlist(lapply(res_wwz_list, function(lst_r) lst_r$data$WWZ_plot))),
              na.rm = TRUE)

p_wwz1 <- make_wwz_plot(res_completa$data, "WWZ — Serie completa", zlim)

plots_wwz <- lapply(seq_along(res_wwz_list), function(j) {
  make_wwz_plot(res_wwz_list[[j]]$data,
                paste0("WWZ — Serie con ", pct[j], "% datos eliminados"),
                zlim)
})

print(p_wwz1)
save_plot("p12_wwz_completaMC.png", p_wwz1)
p_wwz_grid <- arrangeGrob(grobs = plots_wwz, ncol = 2)
grid.arrange(grobs = plots_wwz, ncol = 2)
save_plot("p12_wwz_incompletasMC.png", p_wwz_grid)
