#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WWZ -- Monte Carlo por lotes (K series a la vez)
#
# Para una celda (tau0, w) los pesos y la matriz A = X'WX solo dependen
# de los TIEMPOS, no de x. Se calculan una vez y se aplican a las K
# series simuladas mediante productos de matrices.
#
# Equivalencia con wwz_core (normal equations, beta = A^-1 b):
#   b  = X'W x                 (3 x K)
#   Q  = sum(w * x^2)          (K)
#   Vx = Q - b[1]^2 / Sw       (varianza ponderada total)
#   Vr = Q - b' A^-1 b         (varianza residual)
#   Vy = Vx - Vr = b' A^-1 b - b[1]^2 / Sw
#   Z  = (neff - 3)/2 * Vy / Vr
#
# Estabilidad: la formula Q - b'A^-1 b sufre cancelacion si x tiene un
# offset grande, asi que cada columna se centra y escala ANTES (igual que
# wwz_matrix). Z es invariante a x -> a + s*x.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Una celda para K series ya estandarizadas. Devuelve vector de longitud K.
wwz_cells_batch <- function(wt, X, Xs, Xs2, thr,
                            neff_min = 4, sw_min = 1e-12,
                            tol_fit = 1e-8) {
  K  <- ncol(Xs)
  na <- rep(NA_real_, K)

  Sw <- sum(wt)
  if (!is.finite(Sw) || Sw < sw_min) return(na)
  neff <- Sw^2 / sum(wt^2)
  if (!is.finite(neff) || neff < neff_min) return(na)

  XW <- X * wt                       # n x 3
  A  <- crossprod(XW, X)             # 3 x 3 (solo depende de los tiempos)
  if (rcond(A) < 1e-12) return(na)

  B    <- crossprod(XW, Xs)          # 3 x K
  Q    <- drop(crossprod(wt, Xs2))   # K
  quad <- colSums(B * solve(A, B))   # b' A^-1 b
  b1s  <- B[1, ]^2 / Sw

  Vx <- Q - b1s
  Vy <- pmax(quad - b1s, 0)
  Vr <- pmax(Q - quad, tol_fit * Vx)
  z  <- ((neff - 3) / 2) * Vy / Vr
  z[sqrt(pmax(Vx, 0) / Sw) <= thr] <- NA_real_   # ventana constante
  z
}

# Array [nf x nt x K] con la WWZ de cada columna de Xmat (n x K)
wwz_matrix_batch <- function(Xmat, time, freq, tau,
                             c_decay = 0.0125, neff_min = 4,
                             sw_min = 1e-12, tol_const = 1e-10,
                             tol_fit = 1e-8) {
  Xmat <- as.matrix(Xmat)
  K <- ncol(Xmat)
  mu <- colMeans(Xmat)
  s  <- apply(Xmat, 2, sd)
  s[!is.finite(s) | s == 0] <- NA_real_
  Xs  <- sweep(sweep(Xmat, 2, mu), 2, s, "/")
  Xs[is.na(Xs)] <- 0
  Xs2 <- Xs^2
  # criterio de x constante de wwz_core, que ahora recibe x estandarizada:
  # sd ponderada <= tol_const * max|x_est|
  thr <- tol_const * apply(abs(Xs), 2, max)
  thr[is.na(s)] <- Inf                               # sd = 0 -> todo NA

  out <- array(NA_real_, c(length(freq), length(tau), K))
  for (i in seq_along(freq)) {
    w <- 2 * pi * freq[i]
    X <- cbind(1, cos(w * time), sin(w * time))
    for (j in seq_along(tau)) {
      wt <- exp(-c_decay * w^2 * (time - tau[j])^2)
      out[i, j, ] <- wwz_cells_batch(wt, X, Xs, Xs2, thr,
                                     neff_min, sw_min, tol_fit)
    }
  }
  out
}
