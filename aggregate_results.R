source("cp_sim.R")
files <- list.files(pattern="^scenario_.*\\.rds$")
SC <- lapply(files, readRDS)
names(SC) <- vapply(SC, `[[`, character(1), "label")
if (!length(SC)) stop("No scenario files found")
truth_at <- function(curve, j) switch(curve, B01=.20*c(.5,1,1.5)[j], B02=.03*c(.5,1,1.5)[j], B12=.12*c(.5,1,1.5)[j])
summ <- function(est, se, truth, total=length(est)) {
  ok <- is.finite(est) & is.finite(se) & se > 0
  est <- est[ok]; se <- se[ok]; k <- length(est)
  if (k < 2) return(data.frame(truth=truth,bias=NA,ESE=NA,ASE=NA,relerr=NA,CP=NA,MCSE_bias=NA,MCSE_CP=NA,valid_pct=100*k/total))
  cp <- mean(abs(est-truth) <= 1.96*se)
  data.frame(truth=truth,bias=mean(est)-truth,ESE=sd(est),ASE=mean(se),
             relerr=100*(mean(se)/sd(est)-1),CP=cp,
             MCSE_bias=sd(est)/sqrt(k),MCSE_CP=sqrt(cp*(1-cp)/k),
             valid_pct=100*k/total)
}
rows <- list()
for (S in SC) for (cv in c("B01","B02","B12")) for (j in 1:3) {
  M <- S$M; est <- M[,paste0(cv,j)]; se <- M[,paste0("se",substr(cv,2,3),j)]
  z <- summ(est,se,truth_at(cv,j), total=nrow(M)); z$scenario <- S$label; z$n <- S$n; z$design <- S$design
  z$curve <- cv; z$time <- c(.5,1,1.5)[j]; z$method <- "IV"
  rows[[length(rows)+1]] <- z
  nv <- M[,paste0("nv",substr(cv,2,3),j)]; ok <- is.finite(nv); e <- nv[ok]
  rows[[length(rows)+1]] <- data.frame(truth=truth_at(cv,j),bias=mean(e)-truth_at(cv,j),ESE=sd(e),ASE=NA,relerr=NA,CP=NA,MCSE_bias=sd(e)/sqrt(length(e)),MCSE_CP=NA,valid_pct=100*length(e)/nrow(M),scenario=S$label,n=S$n,design=S$design,curve=cv,time=c(.5,1,1.5)[j],method="Naive-Aalen")
}
out <- do.call(rbind,rows); rownames(out) <- NULL
write.csv(out,"performance_summary.csv",row.names=FALSE)
saveRDS(out,"performance_summary.rds")
writeLines(capture.output(print(out,row.names=FALSE)),"performance_summary_print.txt")
cat("Wrote",nrow(out),"rows from",length(SC),"scenario files.\n")
