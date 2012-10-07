\name{Fx}
\alias{Fx}
\title{
Calculate the target function F(x) in BCGD algorithm
}
\description{
F(x)=f(x)+cP(x),the function just calculates F(x) when given x=x0.(f is smooth,P is a 
proper,convex,lower semicontinuous function,c is the tuning parameter,c>0,see more details
in references)
}
\usage{
Fx(Y, X, lambda1, lambda2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
  A n*m matrix.
}
  \item{X}{
  A n*m matrix,the same size of Y.
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
\references{
Paul Tseng and Sangwoon Yun(2009).A Coordinate Gradient Descent Method for Nonsmooth Separable Minimization.
Math,117,387-423.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.