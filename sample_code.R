## Minimal reader/reviewer example: one reproducible data set and one fit.
## The full Monte Carlo driver is run_manifest.sh; this file is intentionally
## short enough to serve as a supplementary-code example.

suppressPackageStartupMessages(Rcpp::sourceCpp("cp_core.cpp"))
source("cp_lib.R")

set.seed(20260909)
rho <- 0.50
eps_a <- 0.50
sigma_R2 <- 0.2^2 / 4 + eps_a^2 / 3
gammaG <- 2 * rho * sqrt(sigma_R2 / (1 - rho^2))

dat <- cp_dgp(n = 2000, design = "A", gammaG = gammaG,
              eps_a = eps_a, cens_rate = 0.15,
              tau_e = 2.0, tau_g = 1.5)

fit <- cp_fit(dat, tau_e = 2.0, tau_g = 1.5,
              eval_t = c(0.5, 1.0, 1.5),
              eval_s = c(0.5, 1.0, 1.5),
              do_if = TRUE)

se <- cp_se(fit, nrow(dat))
result <- data.frame(
  clock = rep(c("calendar", "calendar", "gap"), each = 3),
  transition = rep(c("01", "02", "12"), each = 3),
  time = rep(c(0.5, 1.0, 1.5), 3),
  estimate = c(fit$B01at, fit$B02at, fit$B12at),
  standard_error = c(se$se01, se$se02, se$se12)
)
print(result, row.names = FALSE)

## A compact diagnostic that should be reported with an application:
cat("Observed entrants:", fit$n_ent, "\n")
cat("Observed post-onset deaths:", fit$n_gap, "\n")
cat("Smallest |D0|:", min(abs(fit$D0)), "\n")
cat("Smallest |D12|:", min(abs(fit$D12)), "\n")

## For an empirical 95% pointwise interval at a selected grid point:
result$lower <- result$estimate - 1.96 * result$standard_error
result$upper <- result$estimate + 1.96 * result$standard_error
print(result, row.names = FALSE)
