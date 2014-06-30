\name{Reformat23}
\alias{Reformat23}
\title{Reformat data}
\description{Reformat data from the 2-way table format (containing allelic
  states) to the 3-way  table format (containing counts of the various alellic states)}
\usage{Reformat23(x)}
\arguments{
  \item{x}{A 2-way table containing allelic states
%%     ~~Describe \code{x} here~~
}
}

\value{A 3-way table with  dimensions (n,l,a) n: number of geographical locations, l:
  number of loci, a: max number of alleles.}
\author{Casper Eriksen}

\keyword{Data format}
