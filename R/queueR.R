poch_complex_n <- function(a,q=a,alpha,maxit=1e6){
  jj <- cbind(Re(a),Im(a),Re(q),Im(q),Inf)
  jj <- pochhammer(a_real=jj[,1],a_imag=jj[,2], q_real=jj[,3],q_imag=jj[,4], n=jj[,5],maxit)
  numerator <- jj[,1] + 1i*jj[,2]

  a <- a*q^alpha
  jj <- cbind(Re(a),Im(a),Re(q),Im(q),Inf)
  jj <- pochhammer(a_real=jj[,1],a_imag=jj[,2], q_real=jj[,3],q_imag=jj[,4], n=jj[,5],maxit)
  denominator <- jj[,1] + 1i*jj[,2]

  return(numerator/denominator)
}

`poch` <- function(a,q=a,n=Inf,maxit=1e6){
  if(any(n != round(n))){return(poch_complex_n(a,q,n))}

  neg_found <- any(n<0)
  if(neg_found){
    wanted <- which(n<0)
    n <- abs(n)
    a[wanted] <- a[wanted]/q[wanted]^n[wanted]
  }

  jj <- cbind(Re(c(a)),Im(c(a)),Re(c(q)),Im(c(q)),c(n))
  jj <- pochhammer(a_real=jj[,1],a_imag=jj[,2], q_real=jj[,3],q_imag=jj[,4], n=jj[,5],maxit)
  
  z <- jj[,1] + 1i*jj[,2]
  if(neg_found){ z[wanted] <- 1/z[wanted] }

  if((!is.complex(q)) && !is.complex(a)){ z <- Re(z) }
  
  if(length(a) >= length(q)){
    attr <- attributes(a)
  } else {
    attr <- attributes(q)
  }

  if((length(a) == length(q)) && (length(a)==1) && (length(n)>1)){ attr <- attributes(n) }
  attributes(z) <- attr
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
    
