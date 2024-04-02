#include <Rcpp.h>
#include <omp.h>
using namespace Rcpp;

// Function to calculate row sums for subsets of columns and count occurrences of 1 and 2
// for each iteration separately, using OpenMP for parallelization
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
IntegerMatrix countOccurrencesInCppParallel(NumericMatrix mat, IntegerMatrix cols) {
  int nrow = mat.nrow();
  int ncalc = cols.nrow();
  IntegerMatrix counts(ncalc, 2); // Result matrix for counts of occurrences of 1 and 2
  
  // Parallelize loop over calculations using OpenMP
  #pragma omp parallel for
  for (int k = 0; k < ncalc; ++k) {
    // Get column indices for current calculation
    IntegerVector cur_cols = cols(k, _);
    int ncol = cur_cols.size();
    
    // Initialize counts for current iteration
    int count1 = 0;
    int count2 = 0;
    
    // Iterate over rows
    for (int i = 0; i < nrow; ++i) {
      double sum = 0;
      // Iterate over selected columns
      for (int j = 0; j < ncol; ++j) {
        sum += mat(i, cur_cols[j] - 1); // Adjust for 0-based indexing
      }
      // Count occurrences of 1 and 2
      if (sum == 1) {
        count1++;
      } else if (sum == 2) {
        count2++;
      }
    }
    
    // Store counts for current iteration in result matrix
    counts(k, 0) = count1;
    counts(k, 1) = count2;
  }
  
  return counts;
}