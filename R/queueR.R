process_args <- function(a,q,n){
  if(any(n<0)){
    stopifnot(all(n==round(n)))
    neg_found <- TRUE
  } else {
    neg_found <- FALSE
  }

  jj <- cbind(a,q,n)
  a <- jj[,1]
  q <- jj[,2]
  n <- Re(jj[,3])
  if(neg_found){
    wanted <- which(n<0)
    n <- abs(n)
    a[wanted] <- a[wanted]/q[wanted]^n[wanted]
  }
  return(cbind(Re(a),Im(a),Re(q),Im(q),n))
}

process_output <- function(n,z){
  if(any(n<0)){
    wanted <- which(n<0)
    z[wanted] <- 1/z[wanted]
  }
  return(z)
}

poch_complex_n <- function(a,q=a,alpha){
  jj <- cbind(Re(a),Im(a),Re(q),Im(q),Inf)
  jj <- pochhammer(a_real=jj[,1],a_imag=jj[,2], q_real=jj[,3],q_imag=jj[,4], n=jj[,5])
  numerator <- jj[,1] + 1i*jj[,2]

  a <- a*q^alpha
  jj <- cbind(Re(a),Im(a),Re(q),Im(q),Inf)
  jj <- pochhammer(a_real=jj[,1],a_imag=jj[,2], q_real=jj[,3],q_imag=jj[,4], n=jj[,5])
  denominator <- jj[,1] + 1i*jj[,2]

  return(numerator/denominator)
}

`poch` <- function(a,q=a,n=Inf){
  if(is.complex(n)){return(poch_complex_n(a,q,n))}
  jj <- process_args(a,q,n)
  jj <- pochhammer(a_real=jj[,1],a_imag=jj[,2], q_real=jj[,3],q_imag=jj[,4], n=jj[,5])
  
  z <- process_output(n, z=jj[,1] + 1i*jj[,2])
  if((!is.complex(q)) && !is.complex(a)){ z <- Re(z) }
  return(z)
}

`pochR` <- function(a,q=a,n=Inf){
  if(any(n<0)){
    wanted <- which(n<0)
    n <- abs(n)
    q[wanted] <- a[wanted]/q[wanted]^n[wanted]
  }

  out <- 1
  outold <- 2
  subt <- a
  i <- 1
  while(any(outold != out) & (i<=n)){
    outold <- out
    out <- out*(1-subt)
    subt <- subt*q
    i <- i+1
  }
  
  if(any(n<0)){ out[wanted] <- 1/out[wanted] }
  return(out)
}

rhs <- function(a,q,x){
  out <- 0
  outold <- 1
  n <- 0
  while(out != outold){
    outold <- out
    out <- out + poch(a,q,n)/poch(q,q,n) * x^n
    n <= n+1
  }
  return(out)
  }
  

disc <- function(a,x,q){

  LHS <- poch(a*x,q)/poch(x,q)
  RHS <- rhs(a,q,x)
  return(c(LHS=LHS,RHS=RHS,diff=RHS-LHS))
}
    
