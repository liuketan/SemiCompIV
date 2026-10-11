## ===========================================================================
##  Data generation, estimator wrapper and population-limit solver for the
##  counting-process IV estimator of manuscript v2.
##
##  Designs (all share b01 = 0.20, b02 = 0.03, b12 = 0.12):
##    A  entry separability (manuscript Prop. 2 example): alpha01 tuned so that
##       alpha01 + b01 E(X|G,U) factors as a(U) c(G).  U confounds everything.
##    B  U-free, onset-time-free post-onset baseline hazard.  Entry
##       factorization FAILS; balance should still hold.
##    C  design A with 0.15 U added to alpha01: factorization broken, baseline
##       IV validity and additive homogeneity retained.
##    D  as B, but the post-onset baseline hazard is onset-time DEPENDENT
##       (still U-free).  Tests whether "U-free" alone suffices for balance.
##    E  as A but with the exposure located near 1.5 rather than 0.5, so the
##       weight exp{X A(T1)} is much more sensitive to the estimated
##       first-event curve.  Entry separability still holds.
## ===========================================================================

## design-D tuning constants: alpha01 = 0.20 + D_aU * U (entry factorization
## fails) and alpha12 = D_a0 + D_bT * T1 (U-free but onset-time dependent)
D_aU <- 0.45; D_a0 <- 0.08; D_bT <- 0.70

alpha12_of <- function(design, U, T1) {
  switch(design,
         A = 0.15 + 0.20 * U,
         B = rep(0.25, length(U)),
         C = 0.15 + 0.20 * U,
         D = D_a0 + D_bT * T1,
         E = 0.15 + 0.20 * U,
         stop("unknown design"))
}
alpha01_of <- function(design, mu, U) {
  switch(design,
         A = 0.2 * (1.2 - mu),
         B = 0.20 + 0.15 * U,
         C = 0.2 * (1.2 - mu) + 0.15 * U,
         D = 0.20 + D_aU * U,
         E = 0.50 - 0.20 * mu,
         stop("unknown design"))
}

cp_dgp <- function(n, design = "A", gammaG = 0.40, eps_a = 0.50,
                   cens_rate = 0.15, cens_G = 0.0,
                   tau_e = 2.0, tau_g = 1.5,
                   b01 = 0.20, b02 = 0.03, b12 = 0.12) {
  tau_max <- tau_e + tau_g
  G  <- rbinom(n, 1, 0.5)
  U  <- rbinom(n, 1, 0.5)
  mu <- if (design == "E") 1.0 + 0.2 * U else 0.2 + 0.2 * U
  X  <- mu + gammaG * G + runif(n, -eps_a, eps_a)

  a01 <- alpha01_of(design, mu, U)
  a02 <- 0.10 + 0.08 * U
  l01 <- a01 + b01 * X
  l02 <- a02 + b02 * X
  stopifnot(all(l01 > 0), all(l02 > 0))

  ltot  <- l01 + l02
  Tfirst <- rexp(n, ltot)
  onset  <- rbinom(n, 1, l01 / ltot) == 1

  a12 <- alpha12_of(design, U, Tfirst)
  l12 <- a12 + b12 * X
  stopifnot(all(l12[onset] > 0))
  W  <- rexp(n, l12)

  T1 <- ifelse(onset, Tfirst, Inf)
  T2 <- ifelse(onset, Tfirst + W, Tfirst)

  hC <- cens_rate * exp(cens_G * G)
  Ct <- rexp(n, hC)

  dfe <- ifelse(Ct < Tfirst, 0L, ifelse(onset, 1L, 2L))
  Tfe <- pmin(Tfirst, Ct)
  isEnt <- as.integer(dfe == 1L & Tfe <= tau_e)

  gapCap <- pmin(Ct - T1, tau_max - T1)          # censoring + administrative
  Wobs <- ifelse(isEnt == 1L, pmin(W, gapCap), 0)
  dgap <- as.integer(isEnt == 1L & W <= gapCap)

  # censoring-model likelihood data
  Vc  <- pmin(Ct, T2, tau_max)
  dC  <- as.integer(Ct < pmin(T2, tau_max))

  data.frame(X = X, G = G, U = U, Tfe = Tfe, dfe = dfe,
             Wobs = Wobs, dgap = dgap, isEnt = isEnt,
             hC_true = hC, Vc = Vc, dC = dC)
}

## --------------------------------------------------------------- nuisances
cp_nuisance <- function(d, cens_G_model = FALSE, known_K = FALSE) {
  n <- nrow(d)
  mhat <- mean(d$G)
  xi_m <- d$G - mhat                                   # P(L = 1) = 1
  if (known_K) {
    return(list(mvec = rep(mhat, n), hC = d$hC_true,
                Zc = matrix(0, n, 1), XiMat = cbind(xi_m), nstrat = 1L,
                use_K = TRUE))
  }
  Zc <- if (cens_G_model) cbind(1, d$G) else matrix(1, n, 1)
  ## exponential censoring MLE: rate is constant within each column pattern
  if (cens_G_model) {
    r0 <- sum(d$dC[d$G == 0]) / sum(d$Vc[d$G == 0])
    r1 <- sum(d$dC[d$G == 1]) / sum(d$Vc[d$G == 1])
    th <- c(log(r0), log(r1) - log(r0))
  } else {
    th <- log(sum(d$dC) / sum(d$Vc))
  }
  hC <- as.vector(exp(Zc %*% th))
  S  <- crossprod(Zc, Zc * (hC * d$Vc)) / n            # E[Z Z' h V]
  res <- d$dC - hC * d$Vc
  xi_th <- t(solve(S, t(Zc * res)))
  list(mvec = rep(mhat, n), hC = hC, Zc = Zc,
       XiMat = cbind(xi_m, xi_th), nstrat = 1L, use_K = TRUE)
}

cp_fit <- function(d, tau_e = 2.0, tau_g = 1.5,
                   eval_t = c(0.5, 1.0, 1.5), eval_s = c(0.5, 1.0, 1.5),
                   do_if = TRUE, use_aL = TRUE, use_H = TRUE,
                   cens_G_model = FALSE, known_K = FALSE) {
  nu <- cp_nuisance(d, cens_G_model, known_K)
  cp_fit_cpp(d$X, d$G, integer(nrow(d)), nu$nstrat,
             d$Tfe, as.integer(d$dfe), d$Wobs, as.integer(d$dgap),
             as.integer(d$isEnt),
             nu$mvec, nu$hC, nu$Zc, nu$XiMat,
             tau_e, tau_g, do_if, nu$use_K,
             eval_t, eval_s, use_aL, use_H)
}

cp_se <- function(fit, n) {
  list(se01 = sqrt(colMeans(fit$phi01at^2) / n),
       se02 = sqrt(colMeans(fit$phi02at^2) / n),
       se12 = sqrt(colMeans(fit$phi12at^2) / n))
}

## ------------------------------------------------- population-limit solver
## Solves dBbar/ds = b12 + rho(s; Bbar)/D12(s; Bbar) (manuscript eq. bias).
## Constant-in-s transition hazards, so the calendar integral is done by
## quadrature over onset time and the exposure error by quadrature over eps.
cp_poplim <- function(design, gammaG = 0.40, eps_a = 0.50, tau_e = 2.0,
                      tau_g = 1.5, b01 = 0.20, b02 = 0.03, b12 = 0.12,
                      n_eps = 121, n_t = 400, n_s = 601) {
  ge <- expand.grid(G = 0:1, U = 0:1)
  pe <- rep(0.25, nrow(ge))
  ee <- (seq_len(n_eps) - 0.5) / n_eps * (2 * eps_a) - eps_a
  we <- rep(1 / n_eps, n_eps)
  tt <- (seq_len(n_t) - 0.5) / n_t * tau_e
  wt <- rep(tau_e / n_t, n_t)

  cell <- expand.grid(k = seq_len(nrow(ge)), j = seq_len(n_eps))
  Gc <- ge$G[cell$k]; Uc <- ge$U[cell$k]
  pc <- pe[cell$k] * we[cell$j]
  muc <- if (design == "E") 1.0 + 0.2 * Uc else 0.2 + 0.2 * Uc
  Xc  <- muc + gammaG * Gc + ee[cell$j]
  a01 <- alpha01_of(design, muc, Uc)
  a02 <- 0.10 + 0.08 * Uc
  entry <- outer(a01 + a02, tt, function(r, t) exp(-r * t)) *
           (a01 + b01 * Xc)                       # cells x t
  entry <- sweep(entry, 2, wt, "*")
  A12 <- if (design == "D") outer(rep(1, length(Xc)), D_a0 + D_bT * tt)
         else matrix(alpha12_of(design, Uc, rep(0, length(Uc))),
                     nrow = length(Xc), ncol = n_t)

  step <- function(s, B) {
    delta <- B - b12 * s
    gapf  <- exp(-A12 * s)
    base  <- entry * gapf
    mass  <- rowSums(base) * exp(Xc * delta)
    amass <- rowSums(base * A12) * exp(Xc * delta)
    pm <- pc * mass
    g  <- sum(pm * Gc) / sum(pm)
    rho <- sum(pc * amass * (Gc - g))
    D12 <- sum(pm * (Gc - g) * Xc)
    c(rho = rho, D12 = D12, g = g, drift = rho / D12)
  }

  sg <- seq(0, tau_g, length.out = n_s)
  h  <- sg[2] - sg[1]
  B  <- numeric(n_s); rho <- numeric(n_s); D12 <- numeric(n_s); gb <- numeric(n_s)
  q <- step(0, 0); rho[1] <- q["rho"]; D12[1] <- q["D12"]; gb[1] <- q["g"]
  for (i in 2:n_s) {                                # RK4
    s0 <- sg[i-1]; y0 <- B[i-1]
    k1 <- b12 + step(s0,         y0)["drift"]
    k2 <- b12 + step(s0 + h/2, y0 + h/2*k1)["drift"]
    k3 <- b12 + step(s0 + h/2, y0 + h/2*k2)["drift"]
    k4 <- b12 + step(s0 + h,   y0 + h*k3)["drift"]
    B[i] <- y0 + h/6*(k1 + 2*k2 + 2*k3 + k4)
    q <- step(sg[i], B[i]); rho[i] <- q["rho"]; D12[i] <- q["D12"]; gb[i] <- q["g"]
  }
  data.frame(s = sg, Bbar = B, rho = rho, D12 = D12, gbar = gb,
             drift = rho / D12, Btrue = b12 * sg)
}
