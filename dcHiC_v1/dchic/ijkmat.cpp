#include <fstream>
#include <sstream>
#include <string>
#include <Rcpp.h>
#include <numeric>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix ijk2mat(NumericMatrix ijk, int row_size, int col_size) {
  NumericMatrix mat_obj(row_size, col_size);
  for (int u = 0; u < ijk.nrow(); u++) {
    int x = ijk(u,0) - 1;
    int y = ijk(u,1) - 1;
    mat_obj(x,y) = ijk(u,2);
    mat_obj(y,x) = ijk(u,2);
  }
  return(mat_obj);
}
