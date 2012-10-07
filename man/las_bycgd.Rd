\name{las_bycgd}
\alias{las_bycgd}
\title{
Multiple Copy Number Variation Detection
}
\description{
Set up a model for multiple copy number variation detection
with fused lasso and group lasso,the solve it using Block Coordinate Gradient
Descent algorithm. 
}
\usage{
las_bycgd(Y, X = Y, lambda1 = 1, lambda2 =1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
  A n*m matrix,in copy number variation detection problem,n denotes #(n) samples 
  and m denotes #(m) probs.
}
  \item{X}{
  A n*m matrix,the initial estimator of Y,with the same size of Y,often given X=Y.
}
  \item{lambda1}{
A real number larger than 0,the first penalty used in the model,controls the overall
copy number alteration amount of the target chromosome.
}
  \item{lambda2}{
A real number larger than 0,the second penalty used in the model,controls the frequency
of the alterations in the target region.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
%%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Paul Tseng and Sangwoon Yun(2009).A Coordinate Gradient Descent Method for Nonsmooth Separable Minimization.
Math,117,387-423.
}
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
##simulated data:5 samples,980 probs.
sigma<-0.4;a<-0.5;i<-matrix(1:980,980,1);n<-5
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
beta15<-las_bycgd(y15)
plot(y15[1,],col=2,xlab="genome order",ylab="log2 ratio",
     main="n=15(the first sample)",ylim=c(-1.5,2))
lines(beta15[1,],col=4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{BCGD}
\keyword{group lasso}