\name{P}
\alias{P}
\title{
Compute P(x) in BCGD algorithm
}
\description{
In BCGD algorithm,the target function is F(x)=f(x)+cP(x).Function P calculates P(x) given x=x0.
}
\usage{
P(X, lambda1, lambda2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  n*m matrix.
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