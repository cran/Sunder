\name{MLCVGauss} 
\alias{MLCVGauss}
\title{Inference and model selection under the assumption of Gaussian distribution of
  allele counts}
\description{Inference and model selection for analysis of geographical
  genetic variation under the assumption of Gaussian distribution of
  allele counts for bi-allelic loci. Parameter estimation by maximization of the likelihood.}
\usage{MLCVGauss(gen, D_G, D_E,theta.max,theta.min,ntrain,nresamp)}
\arguments{
  \item{gen}{
    A matrix with dimensions (n,l) n: number of geographical locations, l:
    number of loci. }
  \item{D_G}{A matrix of geograpical distances}
  \item{D_E}{A matrix of environmental distances}
  \item{theta.max}{Upper bounds for the vector of parameters in theta. Note that
    in theta, the parameters are assumed to be in this order:
    (alpha,beta_G, beta_E, gamma, delta)}
  \item{theta.min}{Lower bounds for the vector of parameters in theta. Note that
    in theta, the parameters are assumed to be in this order:
    (alpha,beta_G, beta_E, gamma, delta)}
  \item{ntrain}{Number of sites used for training. An integer smaller
    than \code{nrow(gen)}. If \code{ntrain} is equal to the number of
    sampling sites, the function estimates parameters on the whole
    dataset and does not perform cross-validation.}
  \item{nresamp}{Number of resamplings. An integer larger than 1.}
}
\value{A list with either a component named mod.lik (containing likelihoods on
  the   validation set for the various models compared) or a vector of
  estimated parameters (if \code{ntrain} is equal to the number of
    sampling sites).}
\author{Gilles Guillot}
\keyword{Likelihood inference}
\keyword{Model selection}
 
\examples{
\dontrun{
nsite <- 200
nloc <- 1000
hap.pop.size <- 100
theta <- c(runif(n=1,.5,10),
           runif(n=1,.01,10),
           runif(n=1,.01,10),
           runif(n=1,.5,1),
           runif(n=1,.01,.1)
           )
mod <- 'G+E' 
dat <- SimSunderData(mod=mod,
                     theta=theta,
                     nsite=nsite,
                     nloc=nloc,
                     hap.pop.size=hap.pop.size,
                     nalM=2,nalm=2, #bi-allelic loci
                     var.par=1,
                     scale.par=3)
gen <- dat$gen[,,1]
D_G <- dat$D_G
D_E <- dat$D_E

res <- MLCVGauss(gen,D_G,D_E,
                 ntrain=nrow(gen)/2,
                 nresamp=3)

which.max(res$mod.lik)
}
}

  

