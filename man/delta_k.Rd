\name{delta_k}
\alias{delta_k}

\title{
Calculate delta used in armijo rule.
}
\description{
We choose alpha by armijo rule.It satisfies F(x+alpha*d)<=alpha*sigma*delta+F(x),where alpha is the largest
element of {alpha_init*beta^j},(j=0,1,...).In BCGD algorithm,the
author defines delta by a function.(See more details in Reference)
}
\usage{
delta_k(Y, X, D, lambda1, lambda2)
}
\arguments{
  \item{Y}{
  A n*m matrix.
}
  \item{X}{
  A n*m matrix,should be the same size of Y.
}
  \item{D}{
  A n*m matrix,the direction in each step,can be calculated by function calculate_direction.
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
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
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
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.