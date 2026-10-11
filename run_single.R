source("cp_sim.R")
gamma_for_rho <- function(rho, eps_a = 0.50) {
  sigma_R2 <- 0.2^2 / 4 + eps_a^2 / 3
  2 * rho * sqrt(sigma_R2 / (1 - rho^2))
}
specs <- list(
  S1a_rho30_n1000=list(n=1000,design="A",gammaG=gamma_for_rho(.30)),
  S1b_rho30_n2000=list(n=2000,design="A",gammaG=gamma_for_rho(.30)),
  S1c_rho30_n4000=list(n=4000,design="A",gammaG=gamma_for_rho(.30)),
  S1d_rho50_n1000=list(n=1000,design="A",gammaG=gamma_for_rho(.50)),
  S1e_rho50_n2000=list(n=2000,design="A",gammaG=gamma_for_rho(.50)),
  S1f_rho50_n4000=list(n=4000,design="A",gammaG=gamma_for_rho(.50)),
  S2_ufree_n2000=list(n=2000,design="B",gammaG=gamma_for_rho(.50)),
  S3a_viol_n2000=list(n=2000,design="C",gammaG=gamma_for_rho(.50)),
  S3b_viol_n4000=list(n=4000,design="C",gammaG=gamma_for_rho(.50)),
  S4a_onset_n2000=list(n=2000,design="D",gammaG=gamma_for_rho(.50)),
  S4b_onset_n4000=list(n=4000,design="D",gammaG=gamma_for_rho(.50)),
  S5a_weak_rho20=list(n=2000,design="A",gammaG=gamma_for_rho(.20)),
  S5b_weak_rho10=list(n=4000,design="A",gammaG=gamma_for_rho(.10)),
  S6a_cens40=list(n=2000,design="A",gammaG=gamma_for_rho(.50),cens_rate=.40),
  S6b_censG=list(n=2000,design="A",gammaG=gamma_for_rho(.50),cens_G=.60,cens_G_model=TRUE),
  S6c_knownK=list(n=2000,design="A",gammaG=gamma_for_rho(.50),known_K=TRUE)
)
label <- Sys.getenv("SCENARIO", "S1e_rho50_n2000")
nsim <- as.integer(Sys.getenv("NSIM", "2000"))
ncore <- as.integer(Sys.getenv("NCORE", "2"))
if (!label %in% names(specs)) stop("Unknown SCENARIO")
assign("NCORE", ncore, envir=.GlobalEnv)
s <- specs[[label]]
z <- do.call(run_scenario, c(list(label=label, nsim=nsim, seed0=20260909), s))
saveRDS(z, paste0("scenario_",label,".rds"))
saveRDS(z$M, paste0("raw_",label,".rds"))
cat("DONE", label, "usable", z$n_ok, "of", nsim, "\n")
