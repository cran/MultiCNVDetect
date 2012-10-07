\name{calculate_direction}
\alias{calculate_direction}
\title{
Calculate direction in each iterative step in BCGD algorithm
}
\description{
In BCGD algorithm,we should calculate the direction d in each iterative step.According to
the multiple copy number variation detection problem,we solve it using chasing method.
}
\usage{
calculate_direction(Y, X, lambda1, lambda2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
  A n*m matrix.
}
  \item{X}{
  A n*m matrix,should be the same size of Y.
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
\value{
It returns a n*m matrix.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{calculate_direction}