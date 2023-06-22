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

  for (int i = 0; i < r; i++) {
          complex<double> out = one;
          complex<double> outold (0,0);
    const complex<double> a (a_real[i], q_imag[i]);
    const complex<double> q (q_real[i], q_imag[i]);
          complex<double> s (a_real[i], a_imag[i]);  // s == subt

	  if( (real(q)==1) && (imag(q) == 0) ){
	    out = zero;
	  } else if((abs(q) == 1) && isinf(n[i])){
	    out = cnan; 
	  } else if((abs(q)  > 1) && isinf(n[i])){
	    out = cinf;
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
	    if(f > max_iterations){ out = cnan; }// not converged
	  } // else clause closes
	  result_real[i] = real(out);
	  result_imag[i] = imag(out);
  } // i loop closes
  return Rcpp::cbind(result_real, result_imag);
}

//[[Rcpp::export]]
NumericVector qexp_C(
		     const NumericVector z_real, const NumericVector z_imag,
		     const NumericVector q_real, const NumericVector q_imag,
		     const NumericVector maxit) {
  
  if(
     (z_real.size() != z_imag.size()) ||
     (z_real.size() != q_real.size()) ||
     (z_real.size() != q_imag.size())
     ){
    throw invalid_argument("Two vectors are not of the same length!");
  }
  
  const int max_iterations = (int) maxit[0];
  if(max_iterations < 0){throw invalid_argument("maxit cannot be negative");}
  
  const size_t r = z_real.size();
  NumericVector result_real(r);
  NumericVector result_imag(r);
  
  for (int i = 0; i < r; i++) {
    complex<double> out    = zero;
    complex<double> outold = one;
    complex<double> term   = one;
    complex<double> qn     = one;
    complex<double> d      = zero;
    const complex<double> z (z_real[i], z_imag[i]);
    const complex<double> q (q_real[i], q_imag[i]);
    unsigned int f=0;
    
    while((real(outold) != real(out) || imag(outold) != imag(out)) && (f < max_iterations)){
      outold = out;// meat (i)
      out += term; // meat (ii)
      d += qn;     // meat (iii)
      qn *= q;     // meat (iv)
      term *= z/d; // meat (v)
      f++;
    } // while loop closes
    if(f > max_iterations){ out = cnan; }// not converged [should not happen here]
    result_real[i] = real(out);
    result_imag[i] = imag(out);
  } // i loop closes
  return Rcpp::cbind(result_real, result_imag);
}

//[[Rcpp::export]]
NumericVector BasicHypergeometric(   // equation 1.2.22 of Gasper
		     const NumericVector z_real, const NumericVector z_imag,
		     const NumericVector q_real, const NumericVector q_imag,
		     const NumericVector a_real, const NumericVector a_imag,
		     const NumericVector b_real, const NumericVector b_imag,
		     const NumericVector maxit) {

  if(
     (z_real.size() != z_imag.size()) ||
     (z_real.size() != q_real.size()) ||
     (z_real.size() != q_imag.size())
     ){
    throw invalid_argument("Two vectors [z,q; Re & Im] are not of the same length!");
  }
   
  if((a_real.size() != a_imag.size())){throw invalid_argument("a_real, a_imag are not of the same length!");}
  if((b_real.size() != b_imag.size())){throw invalid_argument("b_real, b_imag are not of the same length!");}

  const int max_iterations = (int) maxit[0];
  if(max_iterations < 0){throw invalid_argument("maxit cannot be negative");}
  
  const size_t l = z_real.size();
  NumericVector result_real(l);
  NumericVector result_imag(l);

  const size_t r = a_real.size();
  const size_t s = b_real.size();

  const int opsmr = 1 + (int) r - (int) s;
  if(opsmr < 0){throw invalid_argument("r>s: series divergent if z != 0");}

  vector<complex<double>> a(r);
  vector<complex<double>> b(s+1);  // extra "1" for the extra q...

  for(int i=0 ; i<r ; i++){a[i] = complex<double> (a_real[i] , a_imag[i]);}
  for(int i=0 ; i<s ; i++){b[i] = complex<double> (b_real[i] , b_imag[i]);}
  
  for (int i = 0; i < l; i++) {
    complex<double> out    = zero;
    complex<double> outold = one;
    complex<double> term   = one;
    complex<double> qn     = one;
    complex<double> qpc    = one;  // qpc == "q to the power (n choose 2)"
    signed int thesign = 1;
    complex<double> z (z_real[i], z_imag[i]);
    const complex<double> q (q_real[i], q_imag[i]);
    complex<double> squarebracketterm = one;  // not including thesign
    unsigned int f=0;

    b[0] = q;  // ...extra q assigned here!

    while((real(outold) != real(out) || imag(outold) != imag(out)) && (f < max_iterations)){
      outold = out; 
      out += term;
      
      complex<double> numerator   = one;
      complex<double> denominator = one; // sic: (q;q)_n term via extra q
      squarebracketterm *= ((complex <double>) thesign) * qpc;
      
      std::vector<std::complex<double>>::iterator it;
      for (it = a.begin(); it != a.end(); ++it) { numerator   *= one - qn * (*it); }
      for (it = b.begin(); it != b.end(); ++it) { denominator *= one - qn * (*it); }
      
      z *= z;
      out += (numerator/denominator) * pow(squarebracketterm,opsmr) * z * ((complex<double>) thesign);
      thesign *= -1;  
      qn *= q;
      qpc *= qn;
      f++;
    } // while loop closes
    if(f > max_iterations){ out = cnan; }// not converged 
    result_real[i] = real(out);
    result_imag[i] = imag(out);
  } // i loop closes
  return Rcpp::cbind(result_real, result_imag);
}


