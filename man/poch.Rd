\name{poch}
\alias{poch}
\alias{pochhammer}
\title{The q-Pochhammer symbol}
\description{\loadmathjax
Pochhammer's \mjseqn{q}-function, written
\mjeqn{\left(a;q\right)}{(a;q)} or \mjeqn{\left(a;q\right)_n}{omitted},
is often
}
\usage{
poch(a, q = a, n = Inf,maxit=1e6)
}
\arguments{
  \item{a,q}{Main input arguments to the Pochhamer function}
  \item{n}{Number of brackets in product}
  \item{maxit}{Maximum number of iterations}
}
\details{
  We define

  \mjdeqn{\left(a;q\right)_n=\prod_{k=0}^{n-1}\left(1-aq^k\right)}{omitted}

  and require \mjeqn{\left(a;q\right)_0=1}{omitted}.  It is a
  \mjseqn{q}-analogue of the Pochhammer symbol
  \mjeqn{\left(x\right)_n=x(x+1)\cdots(x+n-1)}{omitted} because

  \mjdeqn{ \lim_{q\longrightarrow
      1^-}\frac{\left(q^x;q\right)_n}{(1-q)^n}=\left(x\right)_n.}{omitted}

  Provided \mjeqn{\left|q\right| < 1}{|q| < 1}, we may define

    \mjdeqn{\left(a;q\right)_\infty=\left(a;q\right)=\lim_{n\longrightarrow\infty}\left(a;q\right)_n=
      \prod_{k=0}^{\infty}\left(1-aq^k\right)}{omitted}


    The finite product may be defined in terms of the infinite product:

    \mjdeqn{
      \left(a;q\right)_n=\frac{\left(a;q\right)_\infty}{\left(aq^{-n};q\right)_\infty}
    }{
      omitted} 

    whence

    \mjdeqn{
      \left(a;q\right)_{-n}=\prod_{k=1}^{\infty}
      \frac{1}{\left(1-a/q^k\right)}}{omitted}
  }
\value{Returns a complex vector with attributes inherited from the
      appropriate argument}
\author{Robin K. S. Hankin}
\examples{
poch(0.6+0.2i,0.9)
poch(0.6+0.2i,0.9,6)
}


