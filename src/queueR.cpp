#include "queueR.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector complex_multiply(const NumericVector a_real, const NumericVector a_imag,
                               const NumericVector b_real, const NumericVector b_imag) {
  const size_t n = a_real.size();
  NumericVector result_real(n);
  NumericVector result_imag(n);
  
  for (int i = 0; i < n; i++) {
    const std::complex<double> a(a_real[i], a_imag[i]);
    const std::complex<double> b(b_real[i], b_imag[i]);
    const std::complex<double> product = a * b;
    result_real[i] = std::real(product);
    result_imag[i] = std::imag(product);
  }
  
  return Rcpp::cbind(result_real, result_imag);
}
