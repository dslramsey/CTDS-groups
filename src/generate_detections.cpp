// generate_detections.cpp
// Fast C++ implementation of generate_detections() via Rcpp.
//
// r == rho for uniform points on circle (polar-coordinate identity)
//
// gr() is called once per group size on a batched vector of distances,
//


#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List generate_detections_cpp(
    IntegerVector entries_per_group,
    Function      gr,
    double        w,
    bool          closest = true
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

