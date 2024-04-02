#include <Rcpp.h>
using namespace Rcpp;

// Function to calculate row sums for subsets of columns
// [[Rcpp::export]]
NumericMatrix calculateRowSumsInCpp(NumericMatrix mat, IntegerMatrix cols) {
  int nrow = mat.nrow();
  int ncalc = cols.nrow();
  NumericMatrix results(ncalc, nrow); // Result matrix
  
  // Iterate over calculations
  for (int k = 0; k < ncalc; ++k) {
    // Get column indices for current calculation
    IntegerVector cur_cols = cols(k, _);
    int ncol = cur_cols.size();
    
    // Iterate over rows
    for (int i = 0; i < nrow; ++i) {
      double sum = 0;
      // Iterate over selected columns
      for (int j = 0; j < ncol; ++j) {
        sum += mat(i, cur_cols[j] - 1); // Adjust for 0-based indexing
      }
      results(k, i) = sum;
    }
  }
  
  return results;
}