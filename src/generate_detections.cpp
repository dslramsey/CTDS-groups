// generate_detections.cpp
// Fast C++ implementation of generate_detections() via Rcpp.
//
// r == rho for uniform sector points (polar-coordinate identity)
//
// gr() is called ONCE per group size on a batched vector of distances,
//
// Detection modes :
//   any=FALSE (default): group detected if closest individual detected through gr();
//                        return closest or all distances depending on `closest`.
//   any=TRUE           : each individual detected independently via gr(r_k);
//                        return closest detected distance or all
//                        distances depending on `closest`.

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List generate_detections_cpp(
    IntegerVector entries_per_group,
    Function      gr,
    double        w,
    bool          closest = true,
    bool          any     = false
) {
  int G = entries_per_group.size();
  List out(G);

  for (int g = 0; g < G; g++) {
    int n = g + 1;   // group size (1-indexed)
    int E = entries_per_group[g];

    if (E <= 0) {
      out[g] = R_NilValue;
      continue;
    }

    std::vector<int>    entry_vec;
    std::vector<double> dist_vec;

    // ==================================================================
    // MODE A: any=FALSE  — detection decided by closest individual only
    // ==================================================================
    if (!any) {

      NumericVector r_min(E);
      NumericMatrix all_rho;   // only allocated when !closest && n > 1

      if (!closest && n > 1) {
        all_rho = NumericMatrix(E, n);
        for (int j = 0; j < E; j++) {
          NumericVector rho = sqrt(runif(n)) * w;
          double mn = rho[0];
          for (int k = 1; k < n; k++) if (rho[k] < mn) mn = rho[k];
          r_min[j] = mn;
          for (int k = 0; k < n; k++) all_rho(j, k) = rho[k];
        }
      } else {
        for (int j = 0; j < E; j++) {
          NumericVector rho = sqrt(runif(n)) * w;
          double mn = rho[0];
          for (int k = 1; k < n; k++) if (rho[k] < mn) mn = rho[k];
          r_min[j] = mn;
        }
      }

      // One vectorised gr() call for all E minimum distances
      NumericVector p_vec = as<NumericVector>(gr(r_min));
      NumericVector u     = runif(E);

      for (int j = 0; j < E; j++) {
        if (u[j] >= p_vec[j]) continue;   // not detected

        if (closest || n == 1) {
          entry_vec.push_back(j + 1);
          dist_vec.push_back(r_min[j]);
        } else {
          // return all n distances sorted ascending
          std::vector<double> r_entry(n);
          for (int k = 0; k < n; k++) r_entry[k] = all_rho(j, k);
          std::sort(r_entry.begin(), r_entry.end());
          for (int k = 0; k < n; k++) {
            entry_vec.push_back(j + 1);
            dist_vec.push_back(r_entry[k]);
          }
        }
      }

    // ==================================================================
    // MODE B: any=TRUE  — each individual detected independently at gr(r_k)
    // ==================================================================
    } else {

      int total = E * n;

      // Generate all rho values in one flat vector (E*n values)
      NumericVector rho_all = sqrt(runif(total)) * w;

      // One vectorised gr() call across all individuals
      NumericVector p_all = as<NumericVector>(gr(rho_all));
      NumericVector u_all = runif(total);

      for (int j = 0; j < E; j++) {
        int base = j * n;

        // Collect detected distances for this entry
        std::vector<double> det_r;
        for (int k = 0; k < n; k++) {
          if (u_all[base + k] < p_all[base + k]) {
            det_r.push_back(rho_all[base + k]);
          }
        }
        if (det_r.empty()) continue;

        if (closest) {
          // return only the closest detected individual
          entry_vec.push_back(j + 1);
          dist_vec.push_back(*std::min_element(det_r.begin(), det_r.end()));
        } else {
          // return all n distances sorted — encounter is triggered, all visible
          std::vector<double> r_all(n);
          for (int k = 0; k < n; k++) r_all[k] = rho_all[base + k];
          std::sort(r_all.begin(), r_all.end());
          for (int k = 0; k < n; k++) {
            entry_vec.push_back(j + 1);
            dist_vec.push_back(r_all[k]);
          }
        }
      }
    }

    if (entry_vec.empty()) {
      out[g] = R_NilValue;
    } else {
      out[g] = DataFrame::create(
        Named("entry")    = IntegerVector(entry_vec.begin(), entry_vec.end()),
        Named("group")    = g + 1,
        Named("distance") = NumericVector(dist_vec.begin(), dist_vec.end())
      );
    }
  }

  return out;
}

