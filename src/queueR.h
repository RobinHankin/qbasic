// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-


#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <cmath>


using namespace std;
using namespace Rcpp;


const complex<double> one (1,0);
const complex<double> zero (0,0);
const complex<double> cnan (R_NaN,R_NaN);
const complex<double> cinf (R_PosInf,R_PosInf);
