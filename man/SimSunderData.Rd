\name{SimSunderData}
\alias{SimSunderData}

\title{Simulation of data}
\description{Simulate data under the model assumed by the MCMC inference program}
\usage{
SimSunderData(mod, theta, nsite, nloc, hap.pop.size, nalM, nalm, var.par, scale.par)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{A character string equal to "G+E", "G" or "E"}
  \item{theta}{A vector of length 5 containing values for (alpha,beta_G,beta_E,gamma,delta)}
  \item{nsite}{Number of geographical locations}
  \item{nloc}{Number of loci}
  \item{hap.pop.size}{Haploid population size at each combination geographical location x
    locus. In the main inference function of the package
    (\code{SunderInference}), this number can vary across combinations
    of geographical locations and loci. Currently the present simulation
    function accepts only a single number as input.}
  \item{nalM}{Maximum number of alleles over the various loci}
  \item{nalm}{Minimum number of alleles over the various loci}

  \item{var.par}{Variance of the random field model for the environmental variable}
  \item{scale.par}{Scale parameter in the exponential covariance of the
    random field model for the environmental variable}}

\value{A list containing genotypes, geographical distance, environmental
  distances, model and parameters.  }

\author{Filippo Botta}


\keyword{Random field }
\keyword{Covariance function}

\examples{
data(toydata,package='Sunder')
}