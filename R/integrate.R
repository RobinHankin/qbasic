`q_integrate_01` <- function(f, q=1, maxit=1e6, ...){ # 1.11.1
  if(q==1){return(integrate(f,0,1,...))}
  sum <- 0
  sumold <- sum + 1  # ensures while() loop starts
  i <- 0
  while((sum != sumold) && (i<maxit)){
    sumold <- sum
    sum <- sum + f(q^i, ...)*q^i
    i <- i+1
  }
  return(sum*(1-q))
}

`q_integrate_0inf` <- function(f, q=1, maxit=1e6, ...){ # 1.11.4
  if(q==1){return(integrate(f,0,Inf, ...))}
  sum <- f(1) # i=0
  sumold <- sum+1
  i <- 1  # sic, already have f(0) above, don't want to double count
  while((sum != sumold) && (i<maxit)){
    sumold <- sum
    sum <- sum + f(q^i, ...)*q^i + f(q^(-i),...)/q^i # n !=0 in pairs +/-
    i <- i+1
  }
  return(sum*(1-q))
} 

`q_integrate_minusinf_inf` <- function(f, q=1, maxit=1e6,...){ # 1.11.5
  if(q==1){return(integrate(f,-Inf,Inf,...))}
  sum <- 0
  sumold <- sum+1
  i <- 0  # sic; in pairs already [cf q_integrate_0inf() above]
  while((sum != sumold) && (i<maxit)){
    sumold <- sum
    sum <- sum + (f(q^i,...) + f(-q^i,...))*q^i
    i <- i+1
  }
  return(sum*(1-q))
} 

`q_integrate_0a` <- function(f, q=1, a, maxit=1e6, ...){ # 1.11.3
  if(q==1){
    return(integrate(f, lower=0, upper=a, ...))
  } else if(is.infinite(a)){
    return(q_integrate_0inf(f,q=q,...))
  } else {  # default
    return(a*q_integrate_01(f=function(x){f(a*x)},q=q)) 
  }
}

`q_integrate_ab` <- function(f, q=1, a, b, maxit=1e6, ...){  # 1.11.2
  if(q==1){
    return(integrate(f, a, b, ...))
  } else if((a==-Inf) && (b==Inf)){
    return(q_integrate_minusinf_inf(f,q=q,...))
  } else if((a==Inf) && (b==-Inf)){
    return(-q_integrate_minusinf_inf(f,q=q,...))
  } else if(is.infinite(a) || is.infinite(b)){
    return(NaN)
  } else { # default
    return(q_integrate_0a(f, q=q, b,...) - q_integrate_0a(f, q=q, a,...))
  }
}
