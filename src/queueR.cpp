#include "queueR.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector pochhammer(
			 const NumericVector a_real, const NumericVector a_imag,
			 const NumericVector q_real, const NumericVector q_imag,
			 const NumericVector n) {
  
  if(
     (a_real.size() != a_imag.size()) ||
     (a_real.size() != q_real.size()) ||
     (a_real.size() != a_imag.size()) ||
     (a_real.size() !=      n.size())
     ){
    throw std::invalid_argument("Two vectors are not of the same length!");
  }
  
  const size_t r = a_real.size();
  const int MAX_ITERATIONS = 100000;

  NumericVector result_real(r);
  NumericVector result_imag(r);
  const std::complex<double> one (1,0);
  const std::complex<double> zero (0,0);
  const std::complex<double> nan (R_NaN,R_NaN);
  const std::complex<double> inf (R_PosInf,R_PosInf);

  for (int i = 0; i < r; i++) {
          std::complex<double> out = one;
          std::complex<double> outold (0,0);
    const std::complex<double> a (a_real[i], q_imag[i]);
    const std::complex<double> q (q_real[i], q_imag[i]);
          std::complex<double> s (a_real[i], a_imag[i]);  // s == subt

	  if( (std::real(q)==1) && (std::imag(q) == 0) ){
	    out = zero;
	  } else if(std::abs(q) == 1){
	    out = nan;
	  } else if(std::abs(q) > 1){
	    out = inf;
	  } else {
	    unsigned int f=0;
	    while(
		  ((std::real(outold) != std::real(out)) || (std::imag(outold) != std::imag(out))) &&
		  (f < min((int) n[i], MAX_ITERATIONS))
		  ){
	      outold = out;
	      out *= (one-s);  // the meat
	      s *= q;         // also the meat
	      f++;
	    } // while loop closes
	    if(f>MAX_ITERATIONS){ out = nan; }// not converged
	    result_real[i] = std::real(out);
	    result_imag[i] = std::imag(out);
	  }
  } // i loop closes
  return Rcpp::cbind(result_real, result_imag);
}
