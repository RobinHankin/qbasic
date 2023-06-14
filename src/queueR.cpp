#include "queueR.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector complex_multiply(const NumericVector a_real, const NumericVector a_imag,
                               const NumericVector q_real, const NumericVector q_imag) {
  const size_t n = a_real.size();
  NumericVector result_real(n);
  NumericVector result_imag(n);
  
  const std::complex<double> one (1,0);

  for (int i = 0; i < n; i++) {
          std::complex<double> out (1,0);
          std::complex<double> outold (0,0);
    const std::complex<double> a (a_real[i], q_imag[i]);
    const std::complex<double> q (q_real[i], q_imag[i]);
          std::complex<double> s (a_real[i], a_imag[i]);  // s == subt

    while(
	  (std::real(outold) != std::real(out)) ||
	  (std::imag(outold) != std::imag(out))
	  ){
      outold = out;
      out = out*(one-s);
      s = s * q;
      i++;
    }
    result_real[i] = std::real(out);
    result_imag[i] = std::imag(out);
  } // i loop closes
  return Rcpp::cbind(result_real, result_imag);
}
