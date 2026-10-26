#####################################################################
# Comentarios: al final vamos a usar este código, como veo que el 
# el tiempo lo tenemos encima me he puesto serio hoy por la noche 
# y he generado el codigo R para la WWZ, usando lo que algna vez 
# desarrollé (aunque poco), aunque básicamente utlicé el código que 
# hiciste en R translated from Python, se lo pase a Claude interactuamos 
# varias veces, no me convencia pero se agoto el tiempo gratis y me 
# pasé a Chatgpt con dudas y la verdad que me ha sorprendido, le he 
# pasado el código que había conseguido funcional con Claude basado en 
# tuyo, pero lo ha re-escrito la tía, y ademas me dio opciones para 
# mejorarlo, y aprovechando le dije que implemente el COI y un test 
# de significancia estadísita de dos modos diferentes: 
# mediante un modelo teórico (test paramétrico) de ruido rojo y 
# mediante simualciones de Monte Carlo. Lee este paper para esto último: https://www.sciencedirect.com/science/article/pii/S0098300401000449 es la base de todos estos rollos del analisis espectral via la transtornada de wavelet,
# todo inicio con la de Fourier! En mi paper review explico mucho 
# de eso: https://rclimatol.eu/wp-content/uploads/2022/04/Articulo22e.pdf 

#####################################################################
# Weighted Wavelet Z-Transform (WWZ)
# Implementation in R
#####################################################################
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

#####################################################################

#####################################################################
# (C) 2026. Ricardo Sánchez Santos & Josue M. Polanco-Martinez
# ricardosanchezsantos@usal.es | 
# josue.polanco@usal.es | josue.m.polanco@gmail.com
# Dept. of Statistics, Fac. of Sciences, Univ. of Salamanca, Spain
# March 2026, Salamanca, Spain
# GNU General Public License v3 — https://www.gnu.org/licenses/
#####################################################################


wwz_robustn <- function(time, x,
                       freq = seq(1e-04, 0.001, length.out = 200),
                       #freq = seq(0.01, 0.5, length.out = 80),
                       tau_n = 120,
                       c = 0.0125,
                       sig_level = 0.95,
                       n_sim = 50,
                       use_mc = TRUE,
                       return_plot = TRUE,
                       show_plot = TRUE){

  # =========================
  # 0. CLEAN DATA
  # =========================
  ok <- is.finite(time) & is.finite(x)
  time <- time[ok]; x <- x[ok]

  ord <- order(time)
  time <- time[ord]; x <- x[ord]

  n <- length(x)
  tau <- seq(min(time), max(time), length.out = tau_n)

  nf <- length(freq)
  nt <- length(tau)

  # =========================
  # 1. WWZ CORE
  # =========================
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

    ((n - 3)/2) * (Vx - Vy) / Vy
  }

  # =========================
  # 2. WWZ MATRIX (SAFE)
  # =========================
  WWZ <- matrix(0, nf, nt)

  for(i in 1:nf){
    w <- 2*pi*freq[i]
    for(j in 1:nt){
      WWZ[i, j] <- wwz_core(time, x, w, tau[j])
    }
  }

  # =========================
  # 3. AR(1)
  # =========================
  fit <- arima(x, c(1,0,0))
  phi <- fit$coef["ar1"]
  sigma <- sqrt(fit$sigma2)

  # =========================
  # 4. SIGNIFICANCE (MC AR1)
  # =========================
  WWZ_sig <- matrix(0, nf, nt)

  if(use_mc){

    WWZ_sim <- array(0, c(nf, nt, n_sim))

    for(k in 1:n_sim){

      sim <- arima.sim(n = n, list(ar = phi), sd = sigma)

      for(i in 1:nf){
        w <- 2*pi*freq[i]
        for(j in 1:nt){
          WWZ_sim[i, j, k] <- wwz_core(time, sim, w, tau[j])
        }
      }
    }

    for(i in 1:nf){
      for(j in 1:nt){
        WWZ_sig[i, j] <- quantile(WWZ_sim[i, j, ], sig_level, na.rm = TRUE)
      }
    }

  } else {
    WWZ_sig <- matrix(mean(WWZ), nf, nt)
  }

  signif_mask <- WWZ > WWZ_sig

  # =========================
  # 5. COI (SMOOTH)
  # =========================
  period <- 1 / freq
  coi <- matrix(0, nf, nt)

  for(i in 1:nf){
    edge <- period[i] / sqrt(2)

    for(j in 1:nt){
      dist_edge <- min(tau[j] - min(time), max(time) - tau[j])
      coi[i, j] <- min(dist_edge / edge, 1)
    }
  }

  # =========================
  # 6. DATAFRAME SAFE (NO expand.grid BUGS)
  # =========================
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

  # =========================
  # 7. PLOT (ROBUST + ALWAYS WORKS)
  # =========================
  library(ggplot2)
  library(viridis)

  p <- ggplot(df, aes(x = tau, y = freq)) +

    geom_tile(aes(fill = WWZ_plot, alpha = coi)) +

    scale_fill_viridis(option = "C") +
    scale_alpha(range = c(0.2, 1), guide = "none") +

    labs(
      title = "WWZ robust (AR1 + COI + significance)",
      x = "Time",
      y = "Frequency"
    ) +

    theme_minimal(base_size = 12)

  if(show_plot) print(p)

  # =========================
  # 8. RETURN
  # =========================
  return(list(
    WWZ = WWZ,
    WWZ_significant = signif_mask,
    COI = coi,
    AR1_phi = phi,
    AR1_sigma = sigma,
    plot = p,
    data = df
  ))
}

wwz_robustneff <- function(time, x,
                       freq = seq(1e-04, 0.001, length.out = 200),
                       #freq = seq(0.01, 0.5, length.out = 80),
                       tau_n = 120,
                       c = 0.0125,
                       sig_level = 0.95,
                       n_sim = 50,
                       use_mc = TRUE,
                       return_plot = TRUE,
                       show_plot = TRUE){
  
  # =========================
  # 0. CLEAN DATA
  # =========================
  ok <- is.finite(time) & is.finite(x)
  time <- time[ok]; x <- x[ok]
  
  ord <- order(time)
  time <- time[ord]; x <- x[ord]
  
  n <- length(x)
  tau <- seq(min(time), max(time), length.out = tau_n)
  
  nf <- length(freq)
  nt <- length(tau)
  
  # =========================
  # 1. WWZ CORE
  # =========================
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
    if ((Vx-Vy)<1e-12)return(0)
    
    neff <- max(sum(wt)^2 / sum(wt^2), 3 + 1e-6)
    ((neff - 3) / 2) * (Vx - Vy) / Vy
  }
  
  # =========================
  # 2. WWZ MATRIX (SAFE)
  # =========================
  WWZ <- matrix(0, nf, nt)
  
  for(i in 1:nf){
    w <- 2*pi*freq[i]
    for(j in 1:nt){
      WWZ[i, j] <- wwz_core(time, x, w, tau[j])
    }
  }
  
  # =========================
  # 3. AR(1)
  # =========================
  fit <- arima(x, c(1,0,0))
  phi <- fit$coef["ar1"]
  sigma <- sqrt(fit$sigma2)
  
  # =========================
  # 4. SIGNIFICANCE (MC AR1)
  # =========================
  WWZ_sig <- matrix(0, nf, nt)
  
  if(use_mc){
    
    WWZ_sim <- array(0, c(nf, nt, n_sim))
    
    for(k in 1:n_sim){
      
      sim <- arima.sim(n = n, list(ar = phi), sd = sigma)
      
      for(i in 1:nf){
        w <- 2*pi*freq[i]
        for(j in 1:nt){
          WWZ_sim[i, j, k] <- wwz_core(time, sim, w, tau[j])
        }
      }
    }
    
    for(i in 1:nf){
      for(j in 1:nt){
        WWZ_sig[i, j] <- quantile(WWZ_sim[i, j, ], sig_level, na.rm = TRUE)
      }
    }
    
  } else {
    WWZ_sig <- matrix(mean(WWZ), nf, nt)
  }
  
  signif_mask <- WWZ > WWZ_sig
  
  # =========================
  # 5. COI (SMOOTH)
  # =========================
  period <- 1 / freq
  coi <- matrix(0, nf, nt)
  
  for(i in 1:nf){
    edge <- period[i] / sqrt(2)
    
    for(j in 1:nt){
      dist_edge <- min(tau[j] - min(time), max(time) - tau[j])
      coi[i, j] <- min(dist_edge / edge, 1)
    }
  }
  
  # =========================
  # 6. DATAFRAME SAFE (NO expand.grid BUGS)
  # =========================
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
  
  # =========================
  # 7. PLOT (ROBUST + ALWAYS WORKS)
  # =========================
  library(ggplot2)
  library(viridis)
  
  p <- ggplot(df, aes(x = tau, y = freq)) +
    
    geom_tile(aes(fill = WWZ_plot, alpha = coi)) +
    
    scale_fill_viridis(option = "C") +
    scale_alpha(range = c(0.2, 1), guide = "none") +
    
    labs(
      title = "WWZ robust (AR1 + COI + significance)",
      x = "Time",
      y = "Frequency"
    ) +
    
    theme_minimal(base_size = 12)
  
  if(show_plot) print(p)
  
  # =========================
  # 8. RETURN
  # =========================
  return(list(
    WWZ = WWZ,
    WWZ_significant = signif_mask,
    COI = coi,
    AR1_phi = phi,
    AR1_sigma = sigma,
    plot = p,
    data = df
  ))
}


library(gridExtra)
set.seed(123)
time <- seq(0, 100, by = 0.2)
x    <- sin(2*pi*time/10) + rnorm(length(time), 0, 0.3)

freq_uso <- seq(0.01, 0.5, length.out = 80)

res_n    <- wwz_robustn(time, x, freq = freq_uso, use_mc = TRUE, n_sim = 50, show_plot = FALSE)
res_neff <- wwz_robustneff(time, x, freq = freq_uso, use_mc = TRUE, n_sim = 50, show_plot = FALSE)

grid.arrange(res_n$plot, res_neff$plot, ncol = 1,
             top = grid::textGrob("n global (arriba) vs neff local (abajo)",
                                  gp = grid::gpar(fontsize = 13, fontface = "bold")))

