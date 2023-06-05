`f` <- function(a,q=a,n=Inf){
  out <- 1
  outold <- 2
  subt <- a
  i <- 1
  while(outold != out & (i<=n)){
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
    out <- out + f(a,q,n)/f(q,q,n) * x^n
    n <= n+1
  }
  return(out)
  }
  

disc <- function(a,x,q){

  LHS <- f(a*x,q)/f(x,q)
  RHS <- rhs(a,q,x)
  return(c(LHS=LHS,RHS=RHS,diff=RHS-LHS))
}
    
