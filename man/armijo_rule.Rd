\name{armijo_rule}
\alias{armijo_rule}
\title{
Calculate alpha using armijo_rule in BCGD algorithm 
}
\description{
In BCGD algorithm,armijo_rule is used to get the alpha in each iterative step. 
}
\usage{
armijo_rule(Y,X,D,lambda1,lambda2,alpha)
}

\arguments{
  \item{Y}{
A n*m matrix,in the multiple copy number variation detection problem,n denotes #(n) samples,while m denotes #(m) probs. 
}
  \item{X}{
A n*m matrix,original estimator of Y,give X as the same size of Y,often let X=Y.
}
  \item{D}{
A n*m matrix,original direction in BCGD algorithm.
}
  \item{lambda1}{
A real number larger than 0,the first penalty used in the model,controls the overall
copy number alteration amount of the target chromosome.
}
  \item{lambda2}{
A real number larger than 0,the second penalty used in the model,controls the frequency
of the alterations in the target region.
}
  \item{alpha}{
Original alpha used in the first step.
}
}
\value{
Returns an object of scale.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{armijo_rule}