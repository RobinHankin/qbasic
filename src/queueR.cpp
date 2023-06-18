#include "queueR.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector pochhammer(
			 const NumericVector a_real, const NumericVector a_imag,
			 const NumericVector q_real, const NumericVector q_imag,
			 const NumericVector n,
			 const NumericVector maxit) {
  
  if(
     (a_real.size() != a_imag.size()) ||
     (a_real.size() != q_real.size()) ||
     (a_real.size() != a_imag.size()) ||
     (a_real.size() !=      n.size())
     ){
    throw invalid_argument("Two vectors are not of the same length!");
  }

  const int max_iterations = (int) maxit[0];
  if(max_iterations < 0){throw invalid_argument("maxit cannot be negative");}
  
  const size_t r = a_real.size();

  NumericVector result_real(r);
  NumericVector result_imag(r);
  const complex<double> one (1,0);
  const complex<double> zero (0,0);
  const complex<double> nan (R_NaN,R_NaN);
  const complex<double> inf (R_PosInf,R_PosInf);

  for (int i = 0; i < r; i++) {
          complex<double> out = one;
          complex<double> outold (0,0);
    const complex<double> a (a_real[i], q_imag[i]);
    const complex<double> q (q_real[i], q_imag[i]);
          complex<double> s (a_real[i], a_imag[i]);  // s == subt

	  if( (real(q)==1) && (imag(q) == 0) ){
	    out = zero;
	  } else if((abs(q) == 1) && isinf(n[i])){
	    out = nan; 
	  } else if((abs(q)  > 1) && isinf(n[i])){
	    out = inf;
	  } else {
	    unsigned int f=0;
	    while(
		  ((real(outold) != real(out)) || (imag(outold) != imag(out))) &&
		  (f < min((int) n[i], max_iterations))
		  ){
	      outold = out;
	      out *= (one-s);  // the meat
	      s *= q;         // also the meat
	      f++;
	    } // while loop closes
	    if(f > max_iterations){ out = nan; }// not converged
	  } // else clause closes
	  result_real[i] = real(out);
	  result_imag[i] = imag(out);
  } // i loop closes
  return Rcpp::cbind(result_real, result_imag);
}
