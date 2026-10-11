// ============================================================================
//  Counting-process IV estimator for semi-competing risks.
//  Implements the three recursions of manuscript Section 3.3 (eq. est0, est12)
//  and the joint influence functions of Appendix A (eq. if0, if12).
//
//  Data layout (one row per subject):
//    Tfe    first-event / censoring time on the calendar clock
//    dfe    0 = censored first, 1 = onset observed, 2 = direct death observed
//    Wobs   observed gap time = min(T2-T1, C-T1); 0 for non-entrants
//    dgap   1 = post-onset death observed, else 0
//    isEnt  1 if dfe==1 and Tfe <= tau_e
//
//  Censoring is exponential with subject-specific rate hC[i], so
//    K_i(v) = exp(-hC[i]*v),   -d/dtheta log K_i(v) = v * hC[i] * Zc(i,.).
//
//  XiMat is n x q, the influence function of eta-hat: its first `nstrat`
//  columns are the stratum instrument means m_ell, the remaining p columns
//  the censoring parameters theta.
// ============================================================================
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
List cp_fit_cpp(NumericVector X, NumericVector G, IntegerVector Lidx, int nstrat,
                NumericVector Tfe, IntegerVector dfe,
                NumericVector Wobs, IntegerVector dgap, IntegerVector isEnt,
                NumericVector mvec, NumericVector hC, NumericMatrix Zc,
                NumericMatrix XiMat,
                double tau_e, double tau_g, bool do_if, bool use_K,
                NumericVector eval_t, NumericVector eval_s,
                bool use_aL, bool use_H)
{
  const int n = X.size();
  const int q = XiMat.ncol();
  const int p = Zc.ncol();
  const double dn = (double) n;

  std::vector<double> r(n);
  for (int i = 0; i < n; ++i) r[i] = G[i] - mvec[i];

  // ---- first-event jump times, increasing ---------------------------------
  std::vector<int> ordfe;
  for (int i = 0; i < n; ++i) if (dfe[i] > 0 && Tfe[i] <= tau_e) ordfe.push_back(i);
  std::sort(ordfe.begin(), ordfe.end(),
            [&](int a, int b){ return Tfe[a] < Tfe[b]; });
  const int m = (int) ordfe.size();

  // ---- entrants ------------------------------------------------------------
  std::vector<int> ent;
  for (int i = 0; i < n; ++i) if (isEnt[i] == 1) ent.push_back(i);
  const int ne = (int) ent.size();
  std::vector<int> entpos(n, ne);
  for (int j = 0; j < ne; ++j) entpos[ ent[j] ] = j;

  NumericVector fe_time(m), dB01(m), dB02(m), D0v(m), B01(m), B02(m);
  std::vector<double> ph1(n, 0.0), ph2(n, 0.0), w0(n, 0.0), Yfe(n, 0.0);
  std::vector<double> PhiEnt;                       // n x ne, column major
  if (do_if) PhiEnt.assign((size_t) n * (size_t) ne, 0.0);

  double A = 0.0;                                   // A(t-)
  std::vector<double> dG1(q, 0.0), dG2(q, 0.0);

  const int net = eval_t.size();
  NumericVector B01at(net), B02at(net);
  NumericMatrix phi01at(do_if ? n : 1, do_if ? net : 1),
                phi02at(do_if ? n : 1, do_if ? net : 1);
  int esnap = 0;
  double Bc01 = 0.0, Bc02 = 0.0;

  for (int k = 0; k < m; ++k) {
    const int ik = ordfe[k];
    const double t = Tfe[ik];
    fe_time[k] = t;
    while (esnap < net && eval_t[esnap] < t) {
      B01at[esnap] = Bc01; B02at[esnap] = Bc02;
      if (do_if) for (int i = 0; i < n; ++i) {
        phi01at(i, esnap) = ph1[i]; phi02at(i, esnap) = ph2[i];
      }
      ++esnap;
    }

    double D0 = 0.0;
    for (int i = 0; i < n; ++i) {
      Yfe[i] = (Tfe[i] >= t) ? 1.0 : 0.0;
      if (Yfe[i] > 0.0) {
        double lw = X[i] * A;
        if (use_K) lw += hC[i] * t;                 // 1 / K_i(t)
        w0[i] = std::exp(lw);
      } else w0[i] = 0.0;
      D0 += r[i] * w0[i] * Yfe[i] * X[i];
    }
    D0 /= dn;
    D0v[k] = D0;

    const double num = r[ik] * w0[ik] / dn;
    double d01 = 0.0, d02 = 0.0;
    if (dfe[ik] == 1) d01 = num / D0; else d02 = num / D0;
    dB01[k] = d01; dB02[k] = d02;

    if (do_if) {
      // store the first-event influence at this entrant's onset time (left
      // limit): the jump at t has not been applied yet
      if (dfe[ik] == 1 && isEnt[ik] == 1) {
        const size_t off = (size_t) entpos[ik] * (size_t) n;
        for (int i = 0; i < n; ++i) PhiEnt[off + i] = ph1[i] + ph2[i];
      }
      double k01 = 0.0, k02 = 0.0;
      std::fill(dG1.begin(), dG1.end(), 0.0);
      std::fill(dG2.begin(), dG2.end(), 0.0);
      for (int i = 0; i < n; ++i) {
        if (w0[i] == 0.0) continue;
        const double e1 = ((i == ik && dfe[ik] == 1) ? 1.0 : 0.0) - Yfe[i]*X[i]*d01;
        const double e2 = ((i == ik && dfe[ik] == 2) ? 1.0 : 0.0) - Yfe[i]*X[i]*d02;
        k01 += r[i]*w0[i]*X[i]*e1;
        k02 += r[i]*w0[i]*X[i]*e2;
        dG1[ Lidx[i] ] += (-w0[i]) * e1;            // d/dm_ell of (r w0)
        dG2[ Lidx[i] ] += (-w0[i]) * e2;
        if (use_K) {
          const double c = r[i]*w0[i]*t*hC[i];      // d/dtheta of (r w0)
          for (int b = 0; b < p; ++b) {
            dG1[nstrat + b] += c * Zc(i,b) * e1;
            dG2[nstrat + b] += c * Zc(i,b) * e2;
          }
        }
      }
      k01 /= (dn*D0); k02 /= (dn*D0);
      for (int b = 0; b < q; ++b) { dG1[b] /= (dn*D0); dG2[b] /= (dn*D0); }

      for (int i = 0; i < n; ++i) {
        const double e1 = ((i == ik && dfe[ik] == 1) ? 1.0 : 0.0) - Yfe[i]*X[i]*d01;
        const double e2 = ((i == ik && dfe[ik] == 2) ? 1.0 : 0.0) - Yfe[i]*X[i]*d02;
        const double dz1 = (w0[i] == 0.0) ? 0.0 : r[i]*w0[i]/D0*e1;
        const double dz2 = (w0[i] == 0.0) ? 0.0 : r[i]*w0[i]/D0*e2;
        double g1 = 0.0, g2 = 0.0;
        for (int b = 0; b < q; ++b) { g1 += dG1[b]*XiMat(i,b); g2 += dG2[b]*XiMat(i,b); }
        const double prop = ph1[i] + ph2[i];        // phi(u-)
        ph1[i] += dz1 + g1 + k01*prop;
        ph2[i] += dz2 + g2 + k02*prop;
      }
    }

    A += d01 + d02;
    Bc01 += d01; Bc02 += d02;
    B01[k] = Bc01; B02[k] = Bc02;
  }
  while (esnap < net) {
    B01at[esnap] = Bc01; B02at[esnap] = Bc02;
    if (do_if) for (int i = 0; i < n; ++i) {
      phi01at(i, esnap) = ph1[i]; phi02at(i, esnap) = ph2[i];
    }
    ++esnap;
  }

  // A(T1-) for each entrant
  NumericVector AatT1(ne);
  for (int j = 0; j < ne; ++j) {
    const double t1 = Tfe[ ent[j] ];
    double a = 0.0;
    for (int k = 0; k < m; ++k) { if (fe_time[k] >= t1) break; a += dB01[k] + dB02[k]; }
    AatT1[j] = a;
  }

  // ---- gap-time recursion --------------------------------------------------
  std::vector<int> ordg;
  for (int j = 0; j < ne; ++j) {
    const int i = ent[j];
    if (dgap[i] == 1 && Wobs[i] <= tau_g) ordg.push_back(j);
  }
  std::sort(ordg.begin(), ordg.end(),
            [&](int a, int b){ return Wobs[ ent[a] ] < Wobs[ ent[b] ]; });
  const int mg = (int) ordg.size();

  NumericVector gap_time(mg), dB12(mg), B12(mg), D12v(mg), gbar_v(mg),
                Qmass(mg), dAbase_v(mg);
  std::vector<double> ph12(n, 0.0), wg(ne, 0.0), Y12(ne, 0.0);
  std::vector<double> gs(nstrat, 0.0), qs(nstrat, 0.0), numA(nstrat, 0.0), dAb(nstrat, 0.0);
  std::vector<double> vvec(ne, 0.0), dzv(ne, 0.0), dG12(q, 0.0), hterm(n, 0.0);

  double Bc = 0.0;                                  // B12(s-)
  // The propagation factor (1 + k12) is a scalar, so the recursion
  //   phi_l = (1+k12_l) phi_{l-1} + u_l
  // has closed form phi_l = Prod_l * sum_{l'<=l} u_{l'}/Prod_{l'}.  Keeping the
  // deflated accumulators lets the n x ne product be formed only at the
  // snapshot times instead of at every gap event.
  double Prod = 1.0; double Prod_min = 1.0, Prod_max = 1.0;
  std::vector<double> Zacc(n, 0.0), Vacc(ne, 0.0);
  const int nes = eval_s.size();
  NumericVector B12at(nes);
  NumericMatrix phi12at(do_if ? n : 1, do_if ? nes : 1);
  int ssnap = 0;

  for (int l = 0; l < mg; ++l) {
    const int jl = ordg[l];
    const double s = Wobs[ ent[jl] ];
    gap_time[l] = s;
    while (ssnap < nes && eval_s[ssnap] < s) {
      B12at[ssnap] = Bc;
      if (do_if) {
        for (int i = 0; i < n; ++i) phi12at(i, ssnap) = Prod * Zacc[i];
        for (int j = 0; j < ne; ++j) {
          if (Vacc[j] == 0.0) continue;
          const double v = Prod * Vacc[j];
          const size_t off = (size_t) j * (size_t) n;
          for (int i = 0; i < n; ++i) phi12at(i, ssnap) += v * PhiEnt[off + i];
        }
      }
      ++ssnap;
    }

    std::fill(gs.begin(), gs.end(), 0.0);
    std::fill(qs.begin(), qs.end(), 0.0);
    for (int j = 0; j < ne; ++j) {
      const int i = ent[j];
      Y12[j] = (Wobs[i] >= s) ? 1.0 : 0.0;
      if (Y12[j] > 0.0) {
        double lw = X[i] * (AatT1[j] + Bc);
        if (use_K) lw += hC[i] * (Tfe[i] + s);
        wg[j] = std::exp(lw);
        qs[ Lidx[i] ] += wg[j];
        gs[ Lidx[i] ] += wg[j]*G[i];
      } else wg[j] = 0.0;
    }
    for (int e = 0; e < nstrat; ++e) if (qs[e] > 0.0) gs[e] /= qs[e];

    double D12 = 0.0;
    for (int j = 0; j < ne; ++j) {
      if (Y12[j] == 0.0) continue;
      const int i = ent[j];
      D12 += wg[j] * (G[i] - gs[ Lidx[i] ]) * X[i];
    }
    D12 /= dn;
    D12v[l] = D12;
    gbar_v[l] = gs[ Lidx[ ent[jl] ] ];
    Qmass[l]  = qs[ Lidx[ ent[jl] ] ] / dn;

    const double d12 = wg[jl] * (G[ ent[jl] ] - gs[ Lidx[ ent[jl] ] ]) / (dn * D12);
    dB12[l] = d12;

    if (do_if) {
      // feasible weighted baseline post-onset increment (eq. amean)
      std::fill(numA.begin(), numA.end(), 0.0);
      for (int j = 0; j < ne; ++j) {
        if (Y12[j] == 0.0) continue;
        const int i = ent[j];
        numA[ Lidx[i] ] += wg[j] * ( ((j == jl) ? 1.0 : 0.0) - X[i]*d12 );
      }
      for (int e = 0; e < nstrat; ++e)
        dAb[e] = (use_aL && qs[e] > 0.0) ? numA[e]/qs[e] : 0.0;
      dAbase_v[l] = dAb[ Lidx[ ent[jl] ] ];

      double k12 = 0.0;
      std::fill(dG12.begin(), dG12.end(), 0.0);
      std::fill(vvec.begin(), vvec.end(), 0.0);
      std::fill(dzv.begin(),  dzv.end(),  0.0);
      for (int j = 0; j < ne; ++j) {
        if (Y12[j] == 0.0) continue;
        const int i = ent[j];
        const double dR  = ((j == jl) ? 1.0 : 0.0) - X[i]*d12 - dAb[ Lidx[i] ];
        const double cen = wg[j] * (G[i] - gs[ Lidx[i] ]);
        k12    += cen * X[i] * dR;
        vvec[j] = cen * X[i] * dR / (dn * D12);
        dzv[j]  = cen * dR / D12;
        if (use_K) {
          const double c = cen * (Tfe[i] + s) * hC[i] * dR;
          for (int b = 0; b < p; ++b) dG12[nstrat + b] += c * Zc(i,b);
        }
      }
      k12 /= (dn*D12);
      for (int b = 0; b < q; ++b) dG12[b] /= (dn*D12);

      Prod *= (1.0 + k12);
      if (Prod < Prod_min) Prod_min = Prod;
      if (Prod > Prod_max) Prod_max = Prod;
      const double inv = 1.0 / Prod;
      if (use_H) for (int j = 0; j < ne; ++j)
        if (vvec[j] != 0.0) Vacc[j] += vvec[j] * inv;
      for (int i = 0; i < n; ++i) {
        double gterm = 0.0;
        for (int b = 0; b < q; ++b) gterm += dG12[b]*XiMat(i,b);
        Zacc[i] += gterm * inv;
      }
      for (int j = 0; j < ne; ++j)
        if (dzv[j] != 0.0) Zacc[ ent[j] ] += dzv[j] * inv;
    }

    Bc += d12;
    B12[l] = Bc;
  }
  while (ssnap < nes) {
    B12at[ssnap] = Bc;
    if (do_if) {
      for (int i = 0; i < n; ++i) phi12at(i, ssnap) = Prod * Zacc[i];
      for (int j = 0; j < ne; ++j) {
        if (Vacc[j] == 0.0) continue;
        const double v = Prod * Vacc[j];
        const size_t off = (size_t) j * (size_t) n;
        for (int i = 0; i < n; ++i) phi12at(i, ssnap) += v * PhiEnt[off + i];
      }
    }
    ++ssnap;
  }
  if (do_if) for (int i = 0; i < n; ++i) ph12[i] = Prod * Zacc[i];

  IntegerVector entv(ne);
  for (int j = 0; j < ne; ++j) entv[j] = ent[j] + 1;

  List out = List::create(
    _["fe_time"] = fe_time, _["dB01"] = dB01, _["dB02"] = dB02,
    _["B01"] = B01, _["B02"] = B02, _["D0"] = D0v,
    _["gap_time"] = gap_time, _["dB12"] = dB12, _["B12"] = B12,
    _["D12"] = D12v, _["gbar"] = gbar_v, _["Qmass"] = Qmass,
    _["dAbase"] = dAbase_v, _["AatT1"] = AatT1, _["entrant"] = entv,
    _["n_fe"] = m, _["n_ent"] = ne, _["n_gap"] = mg,
    _["B01at"] = B01at, _["B02at"] = B02at, _["B12at"] = B12at,
    _["prod_range"] = NumericVector::create(Prod_min, Prod_max));
  if (do_if) {
    out["phi01"] = NumericVector(ph1.begin(), ph1.end());
    out["phi02"] = NumericVector(ph2.begin(), ph2.end());
    out["phi12"] = NumericVector(ph12.begin(), ph12.end());
    out["phi01at"] = phi01at; out["phi02at"] = phi02at; out["phi12at"] = phi12at;
  }
  return out;
}

// ---------------------------------------------------------------------------
//  Naive Aalen additive-hazards comparator: cause-specific regression on
//  (1, X, G) for the two first-event transitions, and on (1, X, G) among
//  observed entrants for the post-onset transition.  Returns the cumulative
//  X coefficient.  The 3x3 risk-set cross-product is updated incrementally.
// ---------------------------------------------------------------------------
static void solve3(const double* S, const double* b, double* out, bool& ok) {
  double a[9]; for (int i = 0; i < 9; ++i) a[i] = S[i];
  const double det =
      a[0]*(a[4]*a[8]-a[5]*a[7]) - a[3]*(a[1]*a[8]-a[2]*a[7])
    + a[6]*(a[1]*a[5]-a[2]*a[4]);
  if (!R_finite(det) || std::fabs(det) < 1e-12) { ok = false; return; }
  ok = true;
  const double i0 = ( a[4]*a[8]-a[5]*a[7])/det;
  const double i1 = (-(a[1]*a[8]-a[2]*a[7]))/det;
  const double i2 = ( a[1]*a[5]-a[2]*a[4])/det;
  const double i3 = (-(a[3]*a[8]-a[5]*a[6]))/det;
  const double i4 = ( a[0]*a[8]-a[2]*a[6])/det;
  const double i5 = (-(a[0]*a[5]-a[2]*a[3]))/det;
  const double i6 = ( a[3]*a[7]-a[4]*a[6])/det;
  const double i7 = (-(a[0]*a[7]-a[1]*a[6]))/det;
  const double i8 = ( a[0]*a[4]-a[1]*a[3])/det;
  out[0] = i0*b[0]+i3*b[1]+i6*b[2];
  out[1] = i1*b[0]+i4*b[1]+i7*b[2];
  out[2] = i2*b[0]+i5*b[1]+i8*b[2];
}

// [[Rcpp::export]]
List aalen_naive_cpp(NumericVector X, NumericVector G,
                     NumericVector Tfe, IntegerVector dfe,
                     NumericVector Wobs, IntegerVector dgap, IntegerVector isEnt,
                     double tau_e, double tau_g,
                     NumericVector eval_t, NumericVector eval_s)
{
  const int n = X.size();
  NumericVector out01(eval_t.size()), out02(eval_t.size()), out12(eval_s.size());
  // ---- first events --------------------------------------------------------
  {
    std::vector<int> ord(n); for (int i = 0; i < n; ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](int a,int b){ return Tfe[a] < Tfe[b]; });
    double S[9] = {0,0,0,0,0,0,0,0,0};
    for (int i = 0; i < n; ++i) {
      const double z[3] = {1.0, X[i], G[i]};
      for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) S[a*3+b] += z[a]*z[b];
    }
    int ptr = 0; double B01 = 0.0, B02 = 0.0; int e = 0;
    std::vector<int> ev;
    for (int i = 0; i < n; ++i) if (dfe[i] > 0 && Tfe[i] <= tau_e) ev.push_back(i);
    std::sort(ev.begin(), ev.end(), [&](int a,int b){ return Tfe[a] < Tfe[b]; });
    for (size_t k = 0; k < ev.size(); ++k) {
      const int ik = ev[k]; const double t = Tfe[ik];
      while (ptr < n && Tfe[ ord[ptr] ] < t) {          // remove those already gone
        const int j = ord[ptr];
        const double z[3] = {1.0, X[j], G[j]};
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) S[a*3+b] -= z[a]*z[b];
        ++ptr;
      }
      const double b[3] = {1.0, X[ik], G[ik]};
      double coef[3]; bool ok = false; solve3(S, b, coef, ok);
      if (ok) { if (dfe[ik] == 1) B01 += coef[1]; else B02 += coef[1]; }
      while (e < eval_t.size() && (k+1 == ev.size() || Tfe[ ev[k+1] ] > eval_t[e])) {
        if (t <= eval_t[e]) { out01[e] = B01; out02[e] = B02; ++e; } else break;
      }
    }
    while (e < eval_t.size()) { out01[e] = B01; out02[e] = B02; ++e; }
  }
  // ---- post-onset ----------------------------------------------------------
  {
    std::vector<int> ent;
    for (int i = 0; i < n; ++i) if (isEnt[i] == 1) ent.push_back(i);
    const int ne = (int) ent.size();
    if (ne > 3) {
      std::vector<int> ord(ne); for (int j = 0; j < ne; ++j) ord[j] = j;
      std::sort(ord.begin(), ord.end(),
                [&](int a,int b){ return Wobs[ ent[a] ] < Wobs[ ent[b] ]; });
      double S[9] = {0,0,0,0,0,0,0,0,0};
      for (int j = 0; j < ne; ++j) {
        const int i = ent[j]; const double z[3] = {1.0, X[i], G[i]};
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) S[a*3+b] += z[a]*z[b];
      }
      std::vector<int> ev;
      for (int j = 0; j < ne; ++j)
        if (dgap[ ent[j] ] == 1 && Wobs[ ent[j] ] <= tau_g) ev.push_back(j);
      std::sort(ev.begin(), ev.end(),
                [&](int a,int b){ return Wobs[ ent[a] ] < Wobs[ ent[b] ]; });
      int ptr = 0; double B12 = 0.0; int e = 0;
      for (size_t k = 0; k < ev.size(); ++k) {
        const int jl = ev[k]; const double s = Wobs[ ent[jl] ];
        while (ptr < ne && Wobs[ ent[ ord[ptr] ] ] < s) {
          const int i = ent[ ord[ptr] ]; const double z[3] = {1.0, X[i], G[i]};
          for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) S[a*3+b] -= z[a]*z[b];
          ++ptr;
        }
        const int il = ent[jl];
        const double b[3] = {1.0, X[il], G[il]};
        double coef[3]; bool ok = false; solve3(S, b, coef, ok);
        if (ok) B12 += coef[1];
        while (e < eval_s.size() && (k+1 == ev.size() || Wobs[ ent[ ev[k+1] ] ] > eval_s[e])) {
          if (s <= eval_s[e]) { out12[e] = B12; ++e; } else break;
        }
      }
      while (e < eval_s.size()) { out12[e] = B12; ++e; }
    }
  }
  return List::create(_["B01"] = out01, _["B02"] = out02, _["B12"] = out12);
}
