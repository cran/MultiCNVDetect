\name{MultiCNVDetect}
\alias{MultiCNVDetect}
\docType{package}
\title{
Multiple Copy number variation Detection
}
\description{
Detects copy number variation regions with fused lasso and group lasso regression,solved by BCGD algorithm.
}
\details{
\tabular{ll}{
Package: \tab MultiCNVDetect\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2012-10-07\cr
License: \tab GPL(>=2.0)\cr
}
The most important function is \code{las_bycgd} that detects CNV regions.
}

\author{
Hao Lin,Qiang Kou.

Maintainer: Hao Lin <sky_linhao@163.com>
}
\keyword{CNV,BCGD,fused lasso,group lasso}
\examples{
sigma<-0.4;a<-0.5;i<-matrix(1:980,980,1);n<-15
J<-matrix(1:n,ncol=1)
g<-function(i) 0.25*sigma*sin(a*pi*i)
err<-function(i){
  set.seed(i);
  return(rnorm(980,0,sigma))
}
y<-matrix(,n,980);set.seed(123)
a1<-matrix(rnorm(n,0,1),n,300);
a2<-matrix(rnorm(n,0,1),n,380);
y[,1:300]<--0.2
y[,301:480]<-0.8
y[,481:980]<--0.2
y15<-y+t(apply(J,1,err))+apply(J,1,g)
rm(J,i,a,n,sigma,err,g,a1,a2,y)
beta15<-las_bycgd(y15);
plot(y15[1,],col=2,xlab="genome order",ylab="log2 ratio",
     main="n=15(the first sample)",ylim=c(-1.5,2))
lines(beta15[1,],col=4)
}