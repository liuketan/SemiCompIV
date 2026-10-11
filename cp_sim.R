## ===========================================================================
##  ADEMP simulation study for the counting-process IV estimator.
##  Aims      identification under / without weighted entry balance; validity
##            of the Appendix-A influence functions; behaviour under weak
##            instruments, heavier censoring and modelled censoring.
##  Estimands B01(t), B02(t) at t = 0.5, 1.0, 1.5 and B12(s) at s = 0.5, 1.0, 1.5.
##  Methods   proposed recursion with influence-function SE; naive Aalen.
##  Perform.  bias, empirical SE, model SE, coverage of nominal 95% pointwise
##            intervals, all with Monte Carlo standard errors.
##  n_sim     2000, matching Martinussen and Vansteelandt (2020, Biostatistics).
## ===========================================================================
suppressMessages({library(parallel)})
Rcpp::sourceCpp("cp_core.cpp"); source("cp_lib.R")

EV     <- c(0.5, 1.0, 1.5)
TAU_E  <- 2.0; TAU_G <- 1.5
B01T <- 0.20 * EV; B02T <- 0.03 * EV; B12T <- 0.12 * EV
NSIM   <- as.integer(Sys.getenv("NSIM", "2000"))
NCORE  <- as.integer(Sys.getenv("NCORE", "7"))

run_one <- function(seed, n, design, gammaG, eps_a, cens_rate, cens_G,
                    known_K, cens_G_model, ablation, nboot) {
  set.seed(seed)
  d <- cp_dgp(n, design, gammaG = gammaG, eps_a = eps_a,
              cens_rate = cens_rate, cens_G = cens_G,
              tau_e = TAU_E, tau_g = TAU_G)
  f <- cp_fit(d, TAU_E, TAU_G, EV, EV, do_if = TRUE,
              cens_G_model = cens_G_model, known_K = known_K)
  se <- cp_se(f, n)
  nv <- aalen_naive_cpp(d$X, d$G, d$Tfe, as.integer(d$dfe), d$Wobs,
                        as.integer(d$dgap), as.integer(d$isEnt),
                        TAU_E, TAU_G, EV, EV)
  out <- c(B01 = f$B01at, B02 = f$B02at, B12 = f$B12at,
           se01 = se$se01, se02 = se$se02, se12 = se$se12,
           nv01 = nv$B01, nv02 = nv$B02, nv12 = nv$B12,
           n_ent = f$n_ent, n_gap = f$n_gap,
           minD0 = min(abs(f$D0)), minD12 = min(abs(f$D12)),
           cens = mean(d$dfe == 0))
  if (ablation) {
    fa <- cp_fit(d, TAU_E, TAU_G, EV, EV, use_aL = FALSE,
                 cens_G_model = cens_G_model, known_K = known_K)
    fh <- cp_fit(d, TAU_E, TAU_G, EV, EV, use_H = FALSE,
                 cens_G_model = cens_G_model, known_K = known_K)
    out <- c(out, se12_noaL = cp_se(fa, n)$se12, se12_noH = cp_se(fh, n)$se12)
  }
  if (nboot > 0) {
    bs <- matrix(NA_real_, nboot, 3)
    for (b in seq_len(nboot)) {
      idx <- sample.int(n, n, replace = TRUE)
      fb <- try(cp_fit(d[idx, , drop = FALSE], TAU_E, TAU_G, EV, EV,
                       do_if = FALSE, cens_G_model = cens_G_model,
                       known_K = known_K), silent = TRUE)
      if (!inherits(fb, "try-error")) bs[b, ] <- fb$B12at
    }
    out <- c(out, bse12 = apply(bs, 2, sd, na.rm = TRUE))
  }
  out
}

run_scenario <- function(label, n, design = "A", gammaG = 0.40, eps_a = 0.50,
                         cens_rate = 0.15, cens_G = 0.0, known_K = FALSE,
                         cens_G_model = FALSE, ablation = FALSE, nboot = 0,
                         nsim = NSIM, seed0 = 20260909) {
  t0 <- Sys.time()
  res <- mclapply(seq_len(nsim), function(k)
           run_one(seed0 + k, n, design, gammaG, eps_a, cens_rate, cens_G,
                   known_K, cens_G_model, ablation, nboot),
           mc.cores = NCORE)
  ok <- vapply(res, function(z) is.numeric(z) && all(is.finite(z[1:9])), TRUE)
  M <- do.call(rbind, res[ok])
  cat(sprintf("[%s] n=%d design=%s  %d/%d usable  %.1f min\n", label, n, design,
              sum(ok), nsim, as.numeric(Sys.time() - t0, units = "mins")))
  list(label = label, n = n, design = design, gammaG = gammaG,
       cens_rate = cens_rate, cens_G = cens_G, known_K = known_K,
       cens_G_model = cens_G_model, nsim = nsim, n_ok = sum(ok), M = M)
}

## ------------------------------------------------------------- ADEMP summary
summ <- function(est, se, truth) {
  k <- sum(is.finite(est)); est <- est[is.finite(est)]
  b  <- mean(est) - truth
  sd_ <- sd(est)
  ase <- mean(se, na.rm = TRUE)
  cp  <- mean(est - 1.96*se <= truth & est + 1.96*se >= truth, na.rm = TRUE)
  c(truth = truth, bias = b, mcse_bias = sd_/sqrt(k),
    empSE = sd_, mcse_empSE = sd_/sqrt(2*(k-1)),
    modSE = ase,
    mcse_modSE = sd(se, na.rm=TRUE)/sqrt(k),
    relerr = 100*(ase/sd_ - 1),
    cover = cp, mcse_cover = sqrt(cp*(1-cp)/k),
    rmse = sqrt(b^2 + sd_^2), nrep = k)
}

summarise_scenario <- function(S, poplim = NULL) {
  M <- S$M; rows <- NULL
  for (j in 1:3) {
    for (cv in c("B01","B02","B12")) {
      tr <- switch(cv, B01 = B01T[j], B02 = B02T[j], B12 = B12T[j])
      r <- summ(M[, paste0(cv, j)], M[, paste0("se", substr(cv,2,3), j)], tr)
      rows <- rbind(rows, data.frame(scenario = S$label, n = S$n,
                 design = S$design, curve = cv, time = EV[j],
                 estimator = "IV", t(r)))
      nvcol <- paste0("nv", substr(cv,2,3), j)
      e <- M[, nvcol]
      rows <- rbind(rows, data.frame(scenario = S$label, n = S$n,
                 design = S$design, curve = cv, time = EV[j],
                 estimator = "Naive-Aalen",
                 truth = tr, bias = mean(e) - tr, mcse_bias = sd(e)/sqrt(length(e)),
                 empSE = sd(e), mcse_empSE = sd(e)/sqrt(2*(length(e)-1)),
                 modSE = NA, mcse_modSE = NA, relerr = NA, cover = NA,
                 mcse_cover = NA, rmse = sqrt((mean(e)-tr)^2 + sd(e)^2),
                 nrep = length(e)))
      if (cv == "B12" && !is.null(poplim)) {
        pb <- poplim[j]
        r2 <- summ(M[, paste0("B12", j)], M[, paste0("se12", j)], pb)
        rows <- rbind(rows, data.frame(scenario = S$label, n = S$n,
                   design = S$design, curve = "B12", time = EV[j],
                   estimator = "IV vs pop-limit", t(r2)))
      }
    }
  }
  rows
}

show <- function(df, cols = c("curve","time","estimator","truth","bias","mcse_bias",
                              "empSE","modSE","relerr","cover")) {
  d <- df[, cols, drop = FALSE]
  num <- vapply(d, is.numeric, TRUE)
  d[num] <- lapply(d[num], function(z) round(z, 4))
  print(d, row.names = FALSE)
}
