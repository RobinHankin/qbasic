library(queueR)

n <- 500
l <- 9

x <- seq(from=-l,to=l,len=n)
y <- seq(from=-l,to=l,len=n)
z <- outer(x,1i*y,"+")
pdf(file="~/f.pdf")
contour(x,y,Im(cos(2*z*(0.5-0.2i))),asp=1,levels=-10:10,drawlabels=F,axes=F)
contour(x,y,Re(cos(2*z*(0.5-0.2i))),asp=1,levels=-10:10,add=T,lty=3,drawlabels=F)
dev.off()

f <- function(x){qexp_incorrect(x,0.9)}
f <- function(x){(x)}
p <- f(z)
par(pty='s')
contour(x,y,Re(p),asp=1,levels=c(seq(from=-10,to=22,len=153),0.01),
       drawlabels=FALSE,axes=TRUE,xlim=c(-l,l),ylim=c(-l,l))

contour(x,y,Im(p),lty=2,asp=1,levels=c(seq(from=-10,to=22,len=153),0.01),
       drawlabels=FALSE,axes=TRUE,xlim=c(-l,l),ylim=c(-l,l),add=T)

#contour(x,x,Re(f(z*(0.5-0.2i))),asp=1,levels=-10:10,add=T,lty=3,drawlabels=F)

