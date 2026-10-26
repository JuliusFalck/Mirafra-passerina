#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double sliding_window_cor_cpp_penalized(const NumericMatrix& A,
                                        const NumericMatrix& B,
                                        const NumericVector& window_sizes,
                                        int time_step = 12,
                                        int freq_max_shift = 0,
                                        int freq_step = 1,
                                        double lambda_t = 0.0,
                                        double lambda_f = 0.0) {
  int lenA = A.ncol(), lenB = B.ncol();
  
  const NumericMatrix *big=&A, *small=&B;
  int c_max=lenA, c_min=lenB;
  if (lenB > lenA) { big=&B; small=&A; c_max=lenB; c_min=lenA; }
  
  int nrow = big->nrow();
  if (small->nrow() != nrow) stop("A and B must have same number of rows.");
  
  const double* bigPtr = big->begin();
  const double* smallPtr = small->begin();
  
  double best_score = -1e9;
  bool found = false;
  
  for (int wi=0; wi<window_sizes.size(); ++wi) {
    double w = window_sizes[wi];
    if (w<=0.0 || w>1.0) continue;
    
    int win_len = (int)std::floor(w * c_min);
    if (win_len < 2) continue;
    
    int max_i = c_max - win_len;
    int max_j = c_min - win_len;
    if (max_i < 0 || max_j < 0) continue;
    
    for (int s = -freq_max_shift; s <= freq_max_shift; s += freq_step) {
      int abs_s = std::abs(s);
      int nrow_eff = nrow - abs_s;
      if (nrow_eff <= 0) continue;
      
      int r_big0   = (s > 0) ? s : 0;
      int r_small0 = (s < 0) ? -s : 0;
      
      // normalized freq shift in [0,1]
      double df = (nrow > 1) ? (double)abs_s / (double)(nrow - 1) : 0.0;
      
      for (int i=0; i<=max_i; i+=time_step) {
        for (int j=0; j<=max_j; j+=time_step) {
          
          // normalized center offset in [0,1]
          double cb = (c_max > 1) ? ( (i + 0.5*(win_len-1)) / (double)(c_max - 1) ) : 0.0;
          double cs = (c_min > 1) ? ( (j + 0.5*(win_len-1)) / (double)(c_min - 1) ) : 0.0;
          double dt = std::abs(cb - cs);
          
          double sum_x=0, sum_y=0, sum_xx=0, sum_yy=0, sum_xy=0;
          int N = nrow_eff * win_len;
          
          for (int col_off=0; col_off<win_len; ++col_off) {
            int col_big = i + col_off;
            int col_small = j + col_off;
            int base_big = nrow * col_big;
            int base_small = nrow * col_small;
            
            for (int rr=0; rr<nrow_eff; ++rr) {
              int rB = r_big0 + rr;
              int rS = r_small0 + rr;
              double x = bigPtr[base_big + rB];
              double y = smallPtr[base_small + rS];
              sum_x += x; sum_y += y;
              sum_xx += x*x; sum_yy += y*y;
              sum_xy += x*y;
            }
          }
          
          double mean_x = sum_x / N, mean_y = sum_y / N;
          double var_x = sum_xx - N*mean_x*mean_x;
          double var_y = sum_yy - N*mean_y*mean_y;
          if (var_x <= 0.0 || var_y <= 0.0) continue;
          
          double cov_xy = sum_xy - N*mean_x*mean_y;
          double r = cov_xy / std::sqrt(var_x*var_y);
          if (std::isnan(r)) continue;
          
          double penalty = lambda_t * dt + lambda_f * df;
          double score = r - penalty;
          
          if (!found || score > best_score) {
            best_score = score;
            found = true;
          }
        }
      }
    }
  }
  
  return found ? best_score : NA_REAL;
}
