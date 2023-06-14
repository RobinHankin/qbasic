`poch` <- function(a,q=a,n=Inf){  # returns(a;q)_n in Gasper's notation, also known as the q-Pochhammer symbol
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
    
