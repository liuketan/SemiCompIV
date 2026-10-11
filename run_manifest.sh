#!/bin/sh
## Runs the 16 main simulation scenarios, one output file per scenario.
## NSIM and NCORE may be set in the environment; completed scenarios are skipped.
set -eu
NSIM="${NSIM:-2000}"
NCORE="${NCORE:-2}"
for scenario in \
  S1a_rho30_n1000 S1b_rho30_n2000 S1c_rho30_n4000 \
  S1d_rho50_n1000 S1e_rho50_n2000 S1f_rho50_n4000 S2_ufree_n2000 \
  S3a_viol_n2000 S3b_viol_n4000 S4a_onset_n2000 S4b_onset_n4000 \
  S5a_weak_rho20 S5b_weak_rho10 S6a_cens40 S6b_censG S6c_knownK
do
  if [ -f "scenario_${scenario}.rds" ]; then
    echo "SKIP ${scenario}"
    continue
  fi
  echo "START ${scenario}"
  SCENARIO="$scenario" NSIM="$NSIM" NCORE="$NCORE" Rscript run_single.R > "${scenario}.log" 2>&1
  tail -2 "${scenario}.log"
done
