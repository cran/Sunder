\name{SunderInferenceGauss}
\alias{SunderInferenceGauss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Inference in Sunder package under the Gaussian approximation}
\description{
Computes the mean squared cross validation error using leave-one-out 
cross validaton for each of the three models "e+g", "e" and "g". Prediction is performed using
linear projection based on the normal assumption.
}
\usage{
SunderInferenceGauss(data, e, g, ...)
}

\arguments{
  \item{data}{
An array with dimensions (n,l,a) n: number of geographical locations, l:
  number of loci, a: max number of alleles
}
  \item{e}{
A vector of environmental variables with n entries (n: number of geographical locations)
}
  \item{g}{
A matrix of geographical coordinates with n rows and 2 columns (n: number of geographical locations)
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\value{A vector contianing Mean Square Errors in predicting validation set
  genotypes for the three models}

\author{Casper Eriksen}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Gaussian approximation }
\keyword{Cross-validation }
