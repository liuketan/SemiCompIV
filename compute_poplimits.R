## Probability limits of the post-onset recursion (manuscript equation for the
## limit without entry balance), Designs A-D at Corr(G, X) = 0.5.
source("cp_lib.R")
gamma_for_rho <- function(rho, eps_a=.50) {
  sigma_R2 <- .2^2/4 + eps_a^2/3
  2*rho*sqrt(sigma_R2/(1-rho^2))
}
g <- gamma_for_rho(.50)
ans <- lapply(c("A","B","C","D"), function(d) {
  z <- cp_poplim(d, gammaG=g)
  data.frame(design=d, rho=.50,
             s=c(.5,1,1.5), Btrue=.12*c(.5,1,1.5),
             Bbar=sapply(c(.5,1,1.5), function(s) z$Bbar[which.min(abs(z$s-s))]))
})
ans <- do.call(rbind, ans); ans$bias_limit <- ans$Bbar-ans$Btrue
write.csv(ans, "population_limits_rho50.csv", row.names=FALSE)
saveRDS(ans, "population_limits_rho50.rds")
print(ans, row.names=FALSE)
