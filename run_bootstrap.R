source("cp_sim.R")
gamma_for_rho <- function(rho, eps_a=.50) {
  sigma_R2 <- .2^2/4 + eps_a^2/3
  2*rho*sqrt(sigma_R2/(1-rho^2))
}
set.seed(20260909)
nsim <- as.integer(Sys.getenv("BOOT_NSIM", "100"))
nboot <- as.integer(Sys.getenv("BOOT_NBOOT", "99"))
NCORE <- 1L
z <- run_scenario("S7_boot_rho50_n1000", n=1000, design="A",
                  gammaG=gamma_for_rho(.50), nboot=nboot, nsim=nsim,
                  seed0=20262909)
saveRDS(z, "scenario_S7_boot_rho50_n1000.rds")
cat("BOOTSTRAP DONE: datasets=",nsim,"resamples=",nboot,"\n")
